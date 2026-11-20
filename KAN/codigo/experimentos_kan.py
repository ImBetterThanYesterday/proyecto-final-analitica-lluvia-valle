"""Experimentos KAN complementarios a las propuestas de Federico, solo con validacion.

A. Combinacion ancho {4,5} x G {3,5} x lambda_D2 {0,10,100}.
B. KAN residual: cabeza = lineal + KAN que arranca en cero.
R. Referencias con el mismo codificador: cabeza lineal y MLP con parametros igualados.

Uso desde la raiz del proyecto:
    python -X utf8 experimentos_kan.py <trabajador> <n_trabajadores>
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

import kan_variantes_core as core
from kan_variantes_core import (
    DATA, K_E, K_S, LAMBDA_W, LR, M, MAX_EPOCHS, N_ENCODER, PATIENCE, R_S_T, ROOT, SEEDS, W_E_T, W_S_T,
    X_T, Y_T, EXTRA_T, KAN, FunctionalKAN, _indices, validation_loss,
)

SALIDA = ROOT / "KAN" / "salidas" / "corridas"
N_IN = N_ENCODER + 2


class CabezaResidual(torch.nn.Module):
    def __init__(self, widths, grid_intervals):
        super().__init__()
        self.lineal = torch.nn.Linear(N_IN, M)
        self.kan = KAN(widths, grid_intervals)
        ultima = self.kan.layers[-1]
        with torch.no_grad():
            ultima.coef.zero_()
            ultima.w_b.zero_()

    def forward(self, x):
        return self.lineal(x) + self.kan(x)

    def second_difference_penalty(self):
        return self.kan.second_difference_penalty()

    def n_parameters(self):
        return int(sum(p.numel() for p in self.parameters()))


class CabezaSimple(torch.nn.Module):
    def __init__(self, red):
        super().__init__()
        self.red = red

    def forward(self, x):
        return self.red(x)

    def second_difference_penalty(self):
        return torch.zeros((), dtype=torch.float64, device=core.DEVICE)

    def n_parameters(self):
        return int(sum(p.numel() for p in self.parameters()))


def red_con_cabeza(tipo, widths, grid_intervals):
    modelo = FunctionalKAN(widths if tipo in ("kan", "residual") else (N_IN, 4, M), grid_intervals)
    if tipo == "residual":
        modelo.head = CabezaResidual(widths, grid_intervals)
    elif tipo == "lineal":
        modelo.head = CabezaSimple(torch.nn.Linear(N_IN, M))
    elif tipo == "mlp":
        oculta = widths[1]
        modelo.head = CabezaSimple(torch.nn.Sequential(torch.nn.Linear(N_IN, oculta), torch.nn.SiLU(), torch.nn.Linear(oculta, M)))
    return modelo.to(core.DEVICE)


def entrenar(tipo, widths, grid_intervals, seed, lambda_d2):
    if tipo == "kan":
        modelo, info = core.fit(widths, grid_intervals, seed, lambda_d2)
        return core.predict_validation(modelo), info
    torch.manual_seed(seed)
    modelo = red_con_cabeza(tipo, widths, grid_intervals)
    opt = torch.optim.Adam(modelo.parameters(), lr=LR)
    tr, va = _indices("train"), _indices("valid")
    x_tr, e_tr, y_tr = X_T[tr], EXTRA_T[tr], Y_T[tr]
    mejor, mejor_ep, estado, ep = np.inf, 0, None, -1
    t0 = time.time()
    for ep in range(MAX_EPOCHS):
        modelo.train()
        opt.zero_grad()
        dif = modelo(x_tr, e_tr) - y_tr
        obj = torch.einsum("ni,ij,nj->n", dif, W_E_T, dif).mean() + LAMBDA_W * modelo.encoder_roughness(R_S_T, W_S_T)
        if lambda_d2:
            obj = obj + lambda_d2 * modelo.head.second_difference_penalty()
        obj.backward()
        opt.step()
        modelo.eval()
        with torch.no_grad():
            actual = float(validation_loss(modelo, va).item())
        if actual < mejor - 1e-10:
            mejor, mejor_ep = actual, ep
            estado = {k: v.detach().clone() for k, v in modelo.state_dict().items()}
        elif ep - mejor_ep > PATIENCE:
            break
    modelo.load_state_dict(estado)
    modelo.eval()
    with torch.no_grad():
        dif = modelo(x_tr, e_tr) - y_tr
        l2_tr = float(torch.einsum("ni,ij,nj->n", dif, W_E_T, dif).mean().item())
        pred = modelo(X_T[va], EXTRA_T[va]).cpu().numpy()
    info = {"epoca_optima": mejor_ep + 1, "epocas_corridas": ep + 1, "L2_entrenamiento": l2_tr, "L2_validacion": mejor,
            "segundos": round(time.time() - t0, 1), "parametros_cabeza": modelo.head.n_parameters(),
            "parametros_totales": modelo.head.n_parameters() + N_ENCODER * K_S}
    return pred, info


def trabajos():
    lista = []
    for ancho in (4, 5):
        for G in (3, 5):
            for lam in (0.0, 10.0, 100.0):
                lista.append(("A", "kan", ancho, G, lam))
    for G in (3, 5):
        for lam in (0.0, 100.0):
            lista.append(("B", "residual", 4, G, lam))
    lista.append(("R", "lineal", 4, 5, 0.0))
    lista.append(("R", "mlp", 37, 5, 0.0))
    lista.append(("R", "mlp", 47, 5, 0.0))
    return [(grupo, tipo, ancho, G, lam, s) for (grupo, tipo, ancho, G, lam) in lista for s in SEEDS]


def nombre(grupo, tipo, ancho, G, lam):
    if tipo == "lineal":
        return "Lineal (mismo codificador)"
    if tipo == "mlp":
        return f"MLP oculta {ancho}"
    base = "KAN residual" if tipo == "residual" else "KAN"
    return f"{base} [10,{ancho},6] G={G} D2={lam:g}"


def main():
    yo, total = int(sys.argv[1]), int(sys.argv[2])
    torch.set_num_threads(max(1, 12 // total // 2))
    SALIDA.mkdir(parents=True, exist_ok=True)
    pares = DATA["valid_pairs"]
    for j, (grupo, tipo, ancho, G, lam, seed) in enumerate(trabajos()):
        if j % total != yo:
            continue
        archivo = SALIDA / f"t{j:03d}.csv"
        if archivo.exists():
            continue
        widths = (N_IN, ancho, M)
        pred, info = entrenar(tipo, widths, G, seed, lam)
        fila = {"grupo": grupo, "modelo": nombre(grupo, tipo, ancho, G, lam), "tipo": tipo, "ancho": ancho, "G": G,
                "lambda_D2": lam, "semilla": seed, **info}
        marco = pd.DataFrame(pred, columns=[f"c{i + 1}" for i in range(K_E)])
        marco.insert(0, "anio", pares["anio"].to_numpy())
        marco.insert(0, "codigo", pares["codigo"].to_numpy())
        marco.insert(0, "fila", DATA["valid_global_rows"] + 1)
        for k, v in reversed(list(fila.items())):
            marco.insert(0, k, v)
        marco.to_csv(archivo, index=False, float_format="%.12g")
        print(j, fila["modelo"], seed, f"L2_val={info['L2_validacion']:.3f}", f"ep={info['epoca_optima']}", f"{info['segundos']}s", flush=True)


if __name__ == "__main__":
    main()

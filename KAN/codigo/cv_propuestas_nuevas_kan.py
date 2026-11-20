"""Propuestas KAN nuevas evaluadas con validacion cruzada por bloques de 5 anios (1981-2000).

Cambios respecto a la KAN de Federico ([10,4,6], G=5, D2=100, lr=0.01), todos con herramientas de clase:
  norm  estandarizar las neuronas ocultas con la media y DE de entrenamiento antes de la capa 2
  D2    penalizacion de segunda diferencia de las B-splines (rugosidad) con otro peso
  lr    tasa de aprendizaje de Adam
  G     intervalos de la rejilla (G=6 queda dentro del +-10 % de parametros)
Para el bloque 1996-2000 (particion original) se guardan los modelos para graficarlos despues.

Uso desde KAN/codigo:
    python -X utf8 cv_propuestas_nuevas_kan.py <trabajador> <n_trabajadores>
"""

from __future__ import annotations

import sys
import time

import numpy as np
import pandas as pd
import torch

import kan_variantes_core as core
from kan_variantes_core import K_E, K_S, LAMBDA_W, M, MAX_EPOCHS, N_ENCODER, PATIENCE, R_S_T, ROOT, SEEDS, W_E_T, W_S_T, X_T, EXTRA_T, KAN
from cv_bloques_kan import ANIO, BLOQUES, FILAS_ELEGIBLES, RedPliegue, perdida

SALIDA = ROOT / "KAN" / "salidas" / "cv_nuevas"
N_IN = N_ENCODER + 2

PROPUESTAS = {
    "N1 norm + D2=100": dict(norm=True, G=5, d2=100.0, lr=1e-2),
    "N2 norm + D2=30": dict(norm=True, G=5, d2=30.0, lr=1e-2),
    "N3 norm + D2=10": dict(norm=True, G=5, d2=10.0, lr=1e-2),
    "N4 norm + D2=100 + lr 0,003": dict(norm=True, G=5, d2=100.0, lr=3e-3),
    "N5 norm + D2=30 + lr 0,003": dict(norm=True, G=5, d2=30.0, lr=3e-3),
    "N6 norm + G=6 + D2=30": dict(norm=True, G=6, d2=30.0, lr=1e-2),
    "N7 D2=100 + lr 0,003 (sin norm)": dict(norm=False, G=5, d2=100.0, lr=3e-3),
}


class RedNorm(RedPliegue):
    def __init__(self, widths, G, idx_ent, norm):
        super().__init__(widths, G, idx_ent)
        self.norm = norm

    def forward(self, c_sat, extra):
        if not self.norm:
            return super().forward(c_sat, extra)
        capa1, capa2 = self.head.layers
        h = capa1(self.head_input(c_sat, extra))
        h_ent = capa1(self.head_input(X_T[self.idx_ent], EXTRA_T[self.idx_ent]))
        h = (h - h_ent.mean(0)) / h_ent.std(0)
        return self.mu + capa2(h) @ self.basis_fpca.T

    def ocultas(self, c_sat, extra):
        capa1 = self.head.layers[0]
        h = capa1(self.head_input(c_sat, extra))
        if self.norm:
            h_ent = capa1(self.head_input(X_T[self.idx_ent], EXTRA_T[self.idx_ent]))
            h = (h - h_ent.mean(0)) / h_ent.std(0)
        return h


def construir(cfg, seed, idx_ent):
    widths = (N_IN, 4, M)
    torch.manual_seed(seed)
    cabeza = KAN(widths, cfg["G"])
    torch.manual_seed(seed)
    modelo = RedNorm(widths, cfg["G"], idx_ent, cfg["norm"])
    modelo.head.load_state_dict(cabeza.state_dict())
    return modelo.to(core.DEVICE)


def entrenar(cfg, seed, idx_ent, idx_fuera):
    modelo = construir(cfg, seed, idx_ent)
    opt = torch.optim.Adam(modelo.parameters(), lr=cfg["lr"])
    mejor, mejor_ep, estado, ep = np.inf, 0, None, -1
    t0 = time.time()
    for ep in range(MAX_EPOCHS):
        modelo.train()
        opt.zero_grad()
        obj = perdida(modelo, idx_ent) + LAMBDA_W * modelo.encoder_roughness(R_S_T, W_S_T) + cfg["d2"] * modelo.head.second_difference_penalty()
        obj.backward()
        opt.step()
        modelo.eval()
        with torch.no_grad():
            actual = float(perdida(modelo, idx_fuera).item())
        if actual < mejor - 1e-10:
            mejor, mejor_ep = actual, ep
            estado = {k: v.detach().clone() for k, v in modelo.state_dict().items()}
        elif ep - mejor_ep > PATIENCE:
            break
    modelo.load_state_dict(estado)
    modelo.eval()
    with torch.no_grad():
        pred = modelo(X_T[idx_fuera], EXTRA_T[idx_fuera]).cpu().numpy()
        h = modelo.ocultas(X_T[idx_ent], EXTRA_T[idx_ent]).cpu().numpy()
    return modelo, pred, {"epoca_optima": mejor_ep + 1, "epocas_corridas": ep + 1, "L2_fuera": mejor, "segundos": round(time.time() - t0, 1),
                          "parametros_totales": modelo.head.n_parameters() + N_ENCODER * K_S,
                          "pct_ocultas_fuera_rejilla": float(100 * (np.abs(h) > 3).mean())}


def trabajos():
    return [(m, b, s) for m in PROPUESTAS for b in BLOQUES for s in SEEDS]


def main():
    yo, total = int(sys.argv[1]), int(sys.argv[2])
    torch.set_num_threads(max(1, 12 // total // 2))
    (SALIDA / "modelos").mkdir(parents=True, exist_ok=True)
    for j, (nombre, bloque, seed) in enumerate(trabajos()):
        if j % total != yo:
            continue
        archivo = SALIDA / f"nv{j:03d}.csv"
        if archivo.exists():
            continue
        a, b = BLOQUES[bloque]
        fuera = (ANIO >= a) & (ANIO <= b)
        idx_ent = torch.as_tensor(np.flatnonzero(~fuera), dtype=torch.long, device=core.DEVICE)
        idx_fuera = torch.as_tensor(np.flatnonzero(fuera), dtype=torch.long, device=core.DEVICE)
        modelo, pred, info = entrenar(PROPUESTAS[nombre], seed, idx_ent, idx_fuera)
        if bloque == "1996-2000":
            torch.save(modelo.state_dict(), SALIDA / "modelos" / f"{nombre.split()[0]}_semilla{seed}.pt")
        marco = pd.DataFrame(pred, columns=[f"c{i + 1}" for i in range(K_E)])
        marco.insert(0, "fila", FILAS_ELEGIBLES[fuera] + 1)
        for clave, valor in reversed(list({"modelo": nombre, "bloque": bloque, "semilla": seed, **info}.items())):
            marco.insert(0, clave, valor)
        marco.to_csv(archivo, index=False, float_format="%.12g")
        print(j, nombre, bloque, seed, f"L2={info['L2_fuera']:.3f}", f"ep={info['epoca_optima']}",
              f"fuera={info['pct_ocultas_fuera_rejilla']:.0f}%", f"{info['segundos']}s", flush=True)


if __name__ == "__main__":
    main()

"""Validacion cruzada por bloques de 5 anios dentro de 1981-2000 para elegir la KAN finalista.

Cada bloque se deja fuera una vez; se entrena con los otros 15 anios y la parada temprana usa el
bloque fuera, igual que el protocolo de validacion del proyecto. La prueba 2001-2003 no se usa.
El bloque 1996-2000 reproduce la particion original.

Uso desde KAN/codigo:
    python -X utf8 cv_bloques_kan.py <trabajador> <n_trabajadores>
"""

from __future__ import annotations

import sys
import time

import numpy as np
import pandas as pd
import torch

import kan_variantes_core as core
from kan_variantes_core import (
    DATA, K_E, K_S, LAMBDA_W, LR, M, MAX_EPOCHS, N_ENCODER, PATIENCE, R_S_T, ROOT, SEEDS, W_E_T, W_S_T,
    X_T, Y_T, EXTRA_T, KAN, FunctionalKAN,
)
from experimentos_kan import CabezaResidual, CabezaSimple

SALIDA = ROOT / "KAN" / "salidas" / "cv_bloques"
N_IN = N_ENCODER + 2
BLOQUES = {"1981-1985": (1981, 1985), "1986-1990": (1986, 1990), "1991-1995": (1991, 1995), "1996-2000": (1996, 2000)}

pares = pd.read_csv(core.INPUT / "pares.csv", dtype={"codigo": str})
FILAS_ELEGIBLES = np.concatenate([np.flatnonzero(pares["particion"].to_numpy() == p) for p in ("entrenamiento", "validación")])
ANIO = pares["anio"].to_numpy()[FILAS_ELEGIBLES]
assert len(FILAS_ELEGIBLES) == X_T.shape[0]

MODELOS = {
    "KAN [10,4,6] G=5 D2=0 (base)": ("kan", 4, 5, 0.0),
    "KAN [10,4,6] G=5 D2=100 (Federico)": ("kan", 4, 5, 100.0),
    "KAN residual [10,4,6] G=3 D2=100": ("residual", 4, 3, 100.0),
    "MLP oculta 34": ("mlp", 34, 5, 0.0),
    "MLP oculta 37": ("mlp", 37, 5, 0.0),
    "Lineal (mismo codificador)": ("lineal", 4, 5, 0.0),
}


class RedPliegue(FunctionalKAN):
    def __init__(self, widths, G, idx_ent):
        super().__init__(widths, G)
        self.register_buffer("idx_ent", idx_ent)

    def head_input(self, c_sat, extra):
        lat = self.latent(X_T[self.idx_ent])
        return torch.cat([(self.latent(c_sat) - lat.mean(0)) / lat.std(0), extra], dim=1)


def construir(tipo, ancho, G, seed, idx_ent):
    widths = (N_IN, ancho if tipo in ("kan", "residual") else 4, M)
    torch.manual_seed(seed)
    cabeza_kan = KAN(widths, G) if tipo == "kan" else None
    torch.manual_seed(seed)
    modelo = RedPliegue(widths, G, idx_ent)
    if tipo == "kan":
        modelo.head.load_state_dict(cabeza_kan.state_dict())
    elif tipo == "residual":
        modelo.head = CabezaResidual(widths, G)
    elif tipo == "lineal":
        modelo.head = CabezaSimple(torch.nn.Linear(N_IN, M))
    else:
        modelo.head = CabezaSimple(torch.nn.Sequential(torch.nn.Linear(N_IN, ancho), torch.nn.SiLU(), torch.nn.Linear(ancho, M)))
    return modelo.to(core.DEVICE)


def perdida(modelo, idx):
    dif = modelo(X_T[idx], EXTRA_T[idx]) - Y_T[idx]
    return torch.einsum("ni,ij,nj->n", dif, W_E_T, dif).mean()


def entrenar(tipo, ancho, G, lam, seed, idx_ent, idx_fuera):
    modelo = construir(tipo, ancho, G, seed, idx_ent)
    opt = torch.optim.Adam(modelo.parameters(), lr=LR)
    mejor, mejor_ep, estado, ep = np.inf, 0, None, -1
    t0 = time.time()
    for ep in range(MAX_EPOCHS):
        modelo.train()
        opt.zero_grad()
        obj = perdida(modelo, idx_ent) + LAMBDA_W * modelo.encoder_roughness(R_S_T, W_S_T)
        if lam:
            obj = obj + lam * modelo.head.second_difference_penalty()
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
    return pred, {"epoca_optima": mejor_ep + 1, "L2_fuera": mejor, "segundos": round(time.time() - t0, 1),
                  "parametros_totales": modelo.head.n_parameters() + N_ENCODER * K_S}


def trabajos():
    return [(m, b, s) for b in BLOQUES for m in MODELOS for s in SEEDS]


def main():
    yo, total = int(sys.argv[1]), int(sys.argv[2])
    torch.set_num_threads(max(1, 12 // total // 2))
    SALIDA.mkdir(parents=True, exist_ok=True)
    lista = sorted(enumerate(trabajos()), key=lambda x: MODELOS[x[1][0]][0] in ("mlp", "lineal"))
    for k, (j, (modelo, bloque, seed)) in enumerate(lista):
        if k % total != yo:
            continue
        archivo = SALIDA / f"cv{j:03d}.csv"
        if archivo.exists():
            continue
        a, b = BLOQUES[bloque]
        fuera = (ANIO >= a) & (ANIO <= b)
        idx_ent = torch.as_tensor(np.flatnonzero(~fuera), dtype=torch.long, device=core.DEVICE)
        idx_fuera = torch.as_tensor(np.flatnonzero(fuera), dtype=torch.long, device=core.DEVICE)
        tipo, ancho, G, lam = MODELOS[modelo]
        pred, info = entrenar(tipo, ancho, G, lam, seed, idx_ent, idx_fuera)
        marco = pd.DataFrame(pred, columns=[f"c{i + 1}" for i in range(K_E)])
        marco.insert(0, "fila", FILAS_ELEGIBLES[fuera] + 1)
        for clave, valor in reversed(list({"modelo": modelo, "bloque": bloque, "semilla": seed, **info}.items())):
            marco.insert(0, clave, valor)
        marco.to_csv(archivo, index=False, float_format="%.12g")
        print(j, modelo, bloque, seed, f"L2={info['L2_fuera']:.3f}", f"ep={info['epoca_optima']}", f"{info['segundos']}s", flush=True)


if __name__ == "__main__":
    main()

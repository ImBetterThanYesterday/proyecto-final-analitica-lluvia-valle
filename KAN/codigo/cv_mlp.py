"""Propuestas de cabeza MLP evaluadas con validacion cruzada por bloques de 5 anios (1981-2000).

Mismo codificador, entradas, decodificador FPCA, semillas, parada temprana y metrica que la KAN de Federico.
Solo cambia la cabeza. Todas las variantes quedan dentro de +-10 % de los 728 parametros de la KAN.
Herramientas de clase: ancho y profundidad, funcion de activacion, penalizacion L2 de los pesos,
dropout, tasa de aprendizaje y mini-lotes.

Uso desde KAN/codigo:
    python -X utf8 cv_mlp.py <ronda> <trabajador> <n_trabajadores>
"""

from __future__ import annotations

import sys
import time

import numpy as np
import pandas as pd
import torch

import kan_variantes_core as core
from kan_variantes_core import K_E, K_S, LAMBDA_W, M, MAX_EPOCHS, N_ENCODER, PATIENCE, R_S_T, ROOT, SEEDS, W_S_T, X_T, EXTRA_T, KAN
from cv_bloques_kan import ANIO, BLOQUES, FILAS_ELEGIBLES, RedPliegue, perdida
from experimentos_kan import CabezaSimple

SALIDA = ROOT / "KAN" / "salidas" / "cv_mlp"
N_IN = N_ENCODER + 2
ACT = {"silu": torch.nn.SiLU, "relu": torch.nn.ReLU, "tanh": torch.nn.Tanh}
BASE_CFG = dict(capas=(37,), act="silu", l2=0.0, drop=0.0, lr=1e-2, lote=0)

RONDAS = {
    "r1": {
        "M00 MLP37 (reproduccion)": {},
        "M01 MLP37 + L2=10": dict(l2=10.0),
        "M02 MLP37 + L2=100": dict(l2=100.0),
        "M03 MLP37 + L2=1000": dict(l2=1000.0),
        "M04 MLP37 + dropout 0,1": dict(drop=0.1),
        "M05 MLP37 + dropout 0,2": dict(drop=0.2),
        "M06 MLP37 tanh": dict(act="tanh"),
        "M07 MLP37 ReLU": dict(act="relu"),
        "M08 MLP37 + lr 0,003": dict(lr=3e-3),
        "M09 MLP 2 capas 18-18": dict(capas=(18, 18)),
        "M10 MLP 2 capas 24-12": dict(capas=(24, 12)),
        "M11 MLP 3 capas 14-14-14": dict(capas=(14, 14, 14)),
    },
    "r2": {
        "M12 MLP37 + L2=30": dict(l2=30.0),
        "M13 MLP37 + L2=50": dict(l2=50.0),
        "M14 MLP37 + L2=200": dict(l2=200.0),
        "M15 MLP37 + L2=300": dict(l2=300.0),
        "M16 MLP37 + dropout 0,3": dict(drop=0.3),
        "M17 MLP37 + L2=100 + dropout 0,1": dict(l2=100.0, drop=0.1),
        "M18 MLP37 + L2=100 + dropout 0,2": dict(l2=100.0, drop=0.2),
        "M19 MLP 3 capas 14-14-14 + L2=100": dict(capas=(14, 14, 14), l2=100.0),
        "M20 MLP 3 capas 14-14-14 + dropout 0,1": dict(capas=(14, 14, 14), drop=0.1),
        "M21 MLP 2 capas 24-12 + L2=100": dict(capas=(24, 12), l2=100.0),
        "M22 MLP37 + L2=100 + lotes de 256": dict(l2=100.0, lote=256),
        "M23 MLP37 + lotes de 256": dict(lote=256),
    },
    "r3": {
        "M24 MLP 3 capas 14-14-14 + L2=50": dict(capas=(14, 14, 14), l2=50.0),
        "M25 MLP 3 capas 14-14-14 + L2=200": dict(capas=(14, 14, 14), l2=200.0),
        "M26 MLP 3 capas 14-14-14 + L2=300": dict(capas=(14, 14, 14), l2=300.0),
        "M27 MLP 3 capas 14-14-14 + L2=100 + dropout 0,1": dict(capas=(14, 14, 14), l2=100.0, drop=0.1),
        "M28 MLP 3 capas 14-14-14 + dropout 0,2": dict(capas=(14, 14, 14), drop=0.2),
        "M29 MLP 3 capas 14-14-14 + L2=100 + lotes de 256": dict(capas=(14, 14, 14), l2=100.0, lote=256),
        "M30 MLP 3 capas 14-14-14 ReLU + L2=100": dict(capas=(14, 14, 14), act="relu", l2=100.0),
        "M31 MLP 4 capas 12-12-12-12 + L2=100": dict(capas=(12, 12, 12, 12), l2=100.0),
        "M32 MLP 3 capas 16-16-8 + L2=100": dict(capas=(16, 16, 8), l2=100.0),
        "M33 MLP 3 capas 20-12-8 + L2=100": dict(capas=(20, 12, 8), l2=100.0),
        "M34 MLP 3 capas 14-14-14 + L2=100 + lr 0,003": dict(capas=(14, 14, 14), l2=100.0, lr=3e-3),
    },
    "r4": {
        "M35 MLP 3 capas 14-14-14 + L2=200 + dropout 0,1": dict(capas=(14, 14, 14), l2=200.0, drop=0.1),
        "M36 MLP 3 capas 14-14-14 + L2=150 + dropout 0,1": dict(capas=(14, 14, 14), l2=150.0, drop=0.1),
        "M37 MLP 3 capas 14-14-14 + L2=150": dict(capas=(14, 14, 14), l2=150.0),
        "M38 MLP 3 capas 14-14-14 + L2=100 + dropout 0,2": dict(capas=(14, 14, 14), l2=100.0, drop=0.2),
        "M39 MLP 3 capas 20-12-8 + L2=200": dict(capas=(20, 12, 8), l2=200.0),
        "M40 MLP 3 capas 20-12-8 + L2=100 + dropout 0,1": dict(capas=(20, 12, 8), l2=100.0, drop=0.1),
        "M41 MLP 3 capas 20-12-8 + L2=200 + dropout 0,1": dict(capas=(20, 12, 8), l2=200.0, drop=0.1),
        "M42 MLP 3 capas 20-12-8 + dropout 0,2": dict(capas=(20, 12, 8), drop=0.2),
        "M43 MLP 3 capas 18-14-10 + L2=100 + dropout 0,1": dict(capas=(18, 14, 10), l2=100.0, drop=0.1),
    },
    "r5": {
        "M44 MLP 3 capas 20-12-8": dict(capas=(20, 12, 8)),
        "M45 MLP 3 capas 18-14-10": dict(capas=(18, 14, 10)),
        "M46 MLP 3 capas 14-14-14 ReLU": dict(capas=(14, 14, 14), act="relu"),
        "M47 MLP 4 capas 12-12-12-12": dict(capas=(12, 12, 12, 12)),
        "M48 MLP 2 capas 18-18 ReLU": dict(capas=(18, 18), act="relu"),
    },
}


def cfg_completa(cambios):
    return {**BASE_CFG, **cambios}


def cabeza_mlp(cfg):
    capas, entrada = [], N_IN
    for ancho in cfg["capas"]:
        capas += [torch.nn.Linear(entrada, ancho), ACT[cfg["act"]]()]
        if cfg["drop"]:
            capas.append(torch.nn.Dropout(cfg["drop"]))
        entrada = ancho
    capas.append(torch.nn.Linear(entrada, M))
    return CabezaSimple(torch.nn.Sequential(*capas))


def construir(cfg, seed, idx_ent):
    widths = (N_IN, 4, M)
    torch.manual_seed(seed)
    torch.manual_seed(seed)
    modelo = RedPliegue(widths, 5, idx_ent)
    modelo.head = cabeza_mlp(cfg)
    return modelo.to(core.DEVICE)


def penalizacion_l2(modelo):
    pesos = [m.weight for m in modelo.head.modules() if isinstance(m, torch.nn.Linear)]
    return torch.cat([w.reshape(-1) for w in pesos]).square().mean()


def entrenar(cfg, seed, idx_ent, idx_fuera):
    modelo = construir(cfg, seed, idx_ent)
    opt = torch.optim.Adam(modelo.parameters(), lr=cfg["lr"])
    gen = torch.Generator().manual_seed(seed)
    mejor, mejor_ep, estado, ep = np.inf, 0, None, -1
    t0 = time.time()
    for ep in range(MAX_EPOCHS):
        modelo.train()
        lotes = torch.randperm(len(idx_ent), generator=gen).split(cfg["lote"]) if cfg["lote"] else [torch.arange(len(idx_ent))]
        for lote in lotes:
            opt.zero_grad()
            obj = perdida(modelo, idx_ent[lote]) + LAMBDA_W * modelo.encoder_roughness(R_S_T, W_S_T)
            if cfg["l2"]:
                obj = obj + cfg["l2"] * penalizacion_l2(modelo)
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
    return modelo, pred, {"epoca_optima": mejor_ep + 1, "epocas_corridas": ep + 1, "L2_fuera": mejor, "segundos": round(time.time() - t0, 1),
                          "parametros_totales": modelo.head.n_parameters() + N_ENCODER * K_S}


def main():
    ronda, yo, total = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    torch.set_num_threads(max(1, 12 // total // 2))
    carpeta = SALIDA / ronda
    (carpeta / "modelos").mkdir(parents=True, exist_ok=True)
    trabajos = [(m, b, s) for m in RONDAS[ronda] for b in BLOQUES for s in SEEDS]
    for j, (nombre, bloque, seed) in enumerate(trabajos):
        if j % total != yo:
            continue
        archivo = carpeta / f"m{j:03d}.csv"
        if archivo.exists():
            continue
        a, b = BLOQUES[bloque]
        fuera = (ANIO >= a) & (ANIO <= b)
        idx_ent = torch.as_tensor(np.flatnonzero(~fuera), dtype=torch.long, device=core.DEVICE)
        idx_fuera = torch.as_tensor(np.flatnonzero(fuera), dtype=torch.long, device=core.DEVICE)
        modelo, pred, info = entrenar(cfg_completa(RONDAS[ronda][nombre]), seed, idx_ent, idx_fuera)
        if bloque == "1996-2000":
            torch.save(modelo.state_dict(), carpeta / "modelos" / f"{nombre.split()[0]}_semilla{seed}.pt")
        marco = pd.DataFrame(pred, columns=[f"c{i + 1}" for i in range(K_E)])
        marco.insert(0, "fila", FILAS_ELEGIBLES[fuera] + 1)
        for clave, valor in reversed(list({"modelo": nombre, "ronda": ronda, "bloque": bloque, "semilla": seed, **info}.items())):
            marco.insert(0, clave, valor)
        marco.to_csv(archivo, index=False, float_format="%.12g")
        print(j, nombre, bloque, seed, f"L2={info['L2_fuera']:.3f}", f"ep={info['epoca_optima']}", f"{info['segundos']}s", flush=True)


if __name__ == "__main__":
    main()

"""Auditoria de los parametros de la KAN oficial con validacion cruzada por bloques de 5 anios (1981-2000).

Cada variante es la configuracion oficial (entrenar_kan_base_oficial.py) con un solo parametro cambiado.
Mismos datos, bloques, semillas, parada temprana y metrica que la validacion cruzada de la ablacion.
La prueba 2001-2003 no se usa. Ademas del error, cada ajuste guarda dos indicadores de interpretabilidad:
la fraccion de entradas de cada capa fuera de la rejilla y cuanto se alejan sus aristas de una recta
(1 - R2 de un ajuste lineal de phi(x) sobre los datos).

Uso desde KAN/codigo:
    python -X utf8 auditoria_parametros_kan.py correr <trabajador> <n_trabajadores> [codigos...]
    python -X utf8 auditoria_parametros_kan.py verificar
"""

from __future__ import annotations

import sys
import time

import numpy as np
import pandas as pd
import torch

import kan_variantes_core as core
from kan_variantes_core import DATA, K_E, K_S, R_S_T, ROOT, SEEDS, W_E_T, W_S_T, X_T, Y_T, EXTRA_T
from cv_bloques_kan import ANIO, BLOQUES, FILAS_ELEGIBLES

BASE = ROOT / "KAN" / "salidas"
SALIDA = BASE / "auditoria_kan"

OFICIAL = dict(neuronas=8, M=6, lambda_w=1e6, ocultas=(4,), G=5, grado=3, rango=(-3.0, 3.0), estandarizar=True,
               d2=0.0, lr=1e-2, lote=0, max_epocas=5000, paciencia=1000)

VARIANTES = {
    "oficial": ("configuración oficial", {}),
    "ancho2": ("ancho oculto", dict(ocultas=(2,))),
    "ancho3": ("ancho oculto", dict(ocultas=(3,))),
    "ancho5": ("ancho oculto", dict(ocultas=(5,))),
    "ancho6": ("ancho oculto", dict(ocultas=(6,))),
    "prof1": ("profundidad", dict(ocultas=())),
    "prof3": ("profundidad", dict(ocultas=(4, 4))),
    "G3": ("intervalos G", dict(G=3)),
    "G4": ("intervalos G", dict(G=4)),
    "G6": ("intervalos G", dict(G=6)),
    "G8": ("intervalos G", dict(G=8)),
    "grado1": ("grado B-spline", dict(grado=1)),
    "grado2": ("grado B-spline", dict(grado=2)),
    "rango2": ("rango de la rejilla", dict(rango=(-2.0, 2.0))),
    "rango4": ("rango de la rejilla", dict(rango=(-4.0, 4.0))),
    "sinestd": ("estandarizar z", dict(estandarizar=False)),
    "d2_10": ("regularización D2", dict(d2=10.0)),
    "d2_30": ("regularización D2", dict(d2=30.0)),
    "d2_100": ("regularización D2", dict(d2=100.0)),
    "lr3e-3": ("tasa de aprendizaje", dict(lr=3e-3)),
    "lr3e-2": ("tasa de aprendizaje", dict(lr=3e-2)),
    "pac300": ("paciencia", dict(paciencia=300)),
    "pac2000": ("paciencia", dict(paciencia=2000)),
    "lote256": ("tamaño de lote", dict(lote=256)),
    "neur6": ("neuronas funcionales", dict(neuronas=6)),
    "neur7": ("neuronas funcionales", dict(neuronas=7)),
    "M4": ("componentes FPCA M", dict(M=4)),
    "M5": ("componentes FPCA M", dict(M=5)),
    "lw0": ("rugosidad lambda_w", dict(lambda_w=0.0)),
    "lw1e5": ("rugosidad lambda_w", dict(lambda_w=1e5)),
    "lw1e7": ("rugosidad lambda_w", dict(lambda_w=1e7)),
    "lw1e8": ("rugosidad lambda_w", dict(lambda_w=1e8)),
}


def cfg_de(codigo):
    return {**OFICIAL, **VARIANTES[codigo][1]}


class CapaKAN(torch.nn.Module):
    def __init__(self, n_in, n_out, G, grado, rango):
        super().__init__()
        self.grado, self.rango = grado, rango
        paso = (rango[1] - rango[0]) / G
        self.register_buffer("grid", torch.arange(-grado, G + grado + 1, dtype=torch.float64) * paso + rango[0])
        self.coef = torch.nn.Parameter(0.1 * torch.randn(n_in, n_out, G + grado))
        self.w_b = torch.nn.Parameter(torch.randn(n_in, n_out) / np.sqrt(n_in))
        self.w_s = torch.nn.Parameter(torch.ones(n_in, n_out))

    def bases(self, x):
        grid, x = self.grid, x.unsqueeze(-1)
        b = ((x >= grid[:-1]) & (x < grid[1:])).to(x.dtype)
        for d in range(1, self.grado + 1):
            b = ((x - grid[:-(d + 1)]) / (grid[d:-1] - grid[:-(d + 1)]) * b[..., :-1]
                 + (grid[d + 1:] - x) / (grid[d + 1:] - grid[1:-d]) * b[..., 1:])
        return b

    def activations(self, x):
        spline = torch.einsum("nig,iog->nio", self.bases(x), self.coef)
        return self.w_b * torch.nn.functional.silu(x).unsqueeze(-1) + self.w_s * spline

    def forward(self, x):
        return self.activations(x).sum(1)


class RedAuditoria(torch.nn.Module):
    def __init__(self, cfg, seed, idx_ent):
        super().__init__()
        anchos = (cfg["neuronas"] + 2, *cfg["ocultas"], cfg["M"])
        torch.manual_seed(seed)
        self.capas = torch.nn.ModuleList(CapaKAN(a, b, cfg["G"], cfg["grado"], cfg["rango"]) for a, b in zip(anchos[:-1], anchos[1:]))
        torch.manual_seed(seed)
        self.encoder = torch.nn.Parameter(0.1 * torch.randn(cfg["neuronas"], K_S))
        self.estandarizar = cfg["estandarizar"]
        self.register_buffer("mu", core._tensor(DATA["mu"]))
        self.register_buffer("psi", core._tensor(DATA["basis_fpca"][:, :cfg["M"]]))
        self.register_buffer("idx_ent", idx_ent)

    def latent(self, c_sat):
        return c_sat @ W_S_T @ self.encoder.T

    def head_input(self, c_sat, extra):
        z = self.latent(c_sat)
        if self.estandarizar:
            lat = self.latent(X_T[self.idx_ent])
            z = (z - lat.mean(0)) / lat.std(0)
        return torch.cat([z, extra], dim=1)

    def forward(self, c_sat, extra):
        v = self.head_input(c_sat, extra)
        for capa in self.capas:
            v = capa(v)
        return self.mu + v @ self.psi.T

    def rugosidad(self):
        e = self.encoder
        return (torch.einsum("ki,ij,kj->k", e, R_S_T, e) / torch.einsum("ki,ij,kj->k", e, W_S_T, e)).sum()

    def penalizacion_d2(self):
        return torch.stack([torch.diff(c.w_s.unsqueeze(-1) * c.coef, n=2, dim=-1).square().mean() for c in self.capas]).mean()

    def n_parametros(self):
        return int(sum(p.numel() for p in self.parameters()))


def perdida(modelo, idx):
    dif = modelo(X_T[idx], EXTRA_T[idx]) - Y_T[idx]
    return torch.einsum("ni,ij,nj->n", dif, W_E_T, dif).mean()


def diagnostico(modelo):
    salida = {}
    with torch.no_grad():
        v = modelo.head_input(X_T, EXTRA_T)
        for i, capa in enumerate(modelo.capas, 1):
            salida[f"fuera_rejilla_capa{i}"] = float(((v < capa.rango[0]) | (v > capa.rango[1])).double().mean())
            phi = capa.activations(v)
            x = v.unsqueeze(-1).expand_as(phi)
            xc, pc = x - x.mean(0), phi - phi.mean(0)
            r2 = (xc * pc).sum(0).square() / (xc.square().sum(0) * pc.square().sum(0)).clamp_min(1e-30)
            salida[f"no_linealidad_capa{i}"] = float((1 - r2).mean())
            v = capa(v)
    return salida


def entrenar(cfg, seed, idx_ent, idx_fuera):
    modelo = RedAuditoria(cfg, seed, idx_ent).to(core.DEVICE)
    opt = torch.optim.Adam(modelo.parameters(), lr=cfg["lr"])
    gen = torch.Generator().manual_seed(seed)
    mejor, mejor_ep, estado, ep = np.inf, 0, None, -1
    t0 = time.time()
    for ep in range(cfg["max_epocas"]):
        modelo.train()
        lotes = torch.randperm(len(idx_ent), generator=gen).split(cfg["lote"]) if cfg["lote"] else [torch.arange(len(idx_ent))]
        for lote in lotes:
            opt.zero_grad()
            obj = perdida(modelo, idx_ent[lote]) + cfg["lambda_w"] * modelo.rugosidad()
            if cfg["d2"]:
                obj = obj + cfg["d2"] * modelo.penalizacion_d2()
            obj.backward()
            opt.step()
        modelo.eval()
        with torch.no_grad():
            actual = float(perdida(modelo, idx_fuera).item())
        if actual < mejor - 1e-10:
            mejor, mejor_ep = actual, ep
            estado = {k: v.detach().clone() for k, v in modelo.state_dict().items()}
        elif ep - mejor_ep > cfg["paciencia"]:
            break
    modelo.load_state_dict(estado)
    modelo.eval()
    with torch.no_grad():
        pred = modelo(X_T[idx_fuera], EXTRA_T[idx_fuera]).cpu().numpy()
    info = {"epoca_optima": mejor_ep + 1, "epocas_corridas": ep + 1, "L2_fuera": mejor, "segundos": round(time.time() - t0, 1),
            "parametros_totales": modelo.n_parametros(), **diagnostico(modelo)}
    return modelo, pred, info


def correr(yo, total, codigos):
    torch.set_num_threads(max(1, 12 // total // 2))
    (SALIDA / "corridas").mkdir(parents=True, exist_ok=True)
    (SALIDA / "modelos").mkdir(parents=True, exist_ok=True)
    trabajos = [(c, b, s) for c in codigos for b in BLOQUES for s in SEEDS]
    for j, (codigo, bloque, seed) in enumerate(trabajos):
        if j % total != yo:
            continue
        archivo = SALIDA / "corridas" / f"{codigo}__{bloque}__{seed}.csv"
        if archivo.exists():
            continue
        a, b = BLOQUES[bloque]
        fuera = (ANIO >= a) & (ANIO <= b)
        idx_ent = torch.as_tensor(np.flatnonzero(~fuera), dtype=torch.long, device=core.DEVICE)
        idx_fuera = torch.as_tensor(np.flatnonzero(fuera), dtype=torch.long, device=core.DEVICE)
        modelo, pred, info = entrenar(cfg_de(codigo), seed, idx_ent, idx_fuera)
        if bloque == "1996-2000":
            torch.save(modelo.state_dict(), SALIDA / "modelos" / f"{codigo}_semilla{seed}.pt")
        marco = pd.DataFrame(pred, columns=[f"c{i + 1}" for i in range(K_E)])
        marco.insert(0, "fila", FILAS_ELEGIBLES[fuera] + 1)
        cab = {"codigo": codigo, "parametro": VARIANTES[codigo][0], "bloque": bloque, "semilla": seed, **info}
        for clave, valor in reversed(list(cab.items())):
            marco.insert(0, clave, valor)
        marco.to_csv(archivo, index=False, float_format="%.12g")
        print(codigo, bloque, seed, f"L2={info['L2_fuera']:.3f}", f"ep={info['epoca_optima']}", f"{info['segundos']}s", flush=True)


def l2_relativo_por_estacion(g):
    pares = pd.read_csv(core.INPUT / "pares.csv", dtype={"codigo": str})
    coef_mm = pd.read_csv(core.INPUT / "coef_estacion_mm.csv").to_numpy(float)
    phi = core._basis_at_nodes(np.arange(1, 366, dtype=float), DATA["tau_e"])
    filas = g["fila"].to_numpy() - 1
    obs = np.maximum(coef_mm[filas] @ phi.T, 0)
    pred = np.maximum(g[[f"c{i + 1}" for i in range(K_E)]].to_numpy(float) @ phi.T, 0) ** 2
    l2 = np.sqrt(((pred - obs) ** 2).sum(1) / (obs ** 2).sum(1))
    return pd.Series(l2).groupby(pares["codigo"].to_numpy()[filas]).mean()


def verificar():
    cols = [f"c{i + 1}" for i in range(K_E)]
    nueva = pd.concat([pd.read_csv(f) for f in sorted((SALIDA / "corridas").glob("oficial__*.csv"))], ignore_index=True)
    vieja = pd.concat([pd.read_csv(f) for f in sorted((BASE / "cv_bloques").glob("cv*.csv"))], ignore_index=True)
    vieja = vieja[vieja.modelo == "KAN [10,4,6] G=5 D2=0 (base)"]
    llave = ["bloque", "semilla", "fila"]
    union = nueva.merge(vieja, on=llave, suffixes=("", "_v"))
    filas = []
    for (b, s), g in union.groupby(["bloque", "semilla"]):
        filas.append({"bloque": b, "semilla": s, "max_dif_pred": float(np.abs(g[cols].to_numpy() - g[[c + "_v" for c in cols]].to_numpy()).max()),
                      "epoca_nueva": int(g.epoca_optima.iloc[0]), "epoca_vieja": int(g.epoca_optima_v.iloc[0]),
                      "L2_rel_nuevo": l2_relativo_por_estacion(g).mean(),
                      "L2_rel_viejo": l2_relativo_por_estacion(g[["fila", *[c + "_v" for c in cols]]].set_axis(["fila", *cols], axis=1)).mean()})
    t = pd.DataFrame(filas)
    cv = t.groupby("semilla").L2_rel_nuevo.mean().mean()
    fed = t.query("bloque == '1996-2000' and semilla == 2026").iloc[0]
    pruebas = pd.DataFrame([
        {"prueba": "12 ajustes (4 bloques x 3 semillas) emparejados con la CV anterior", "valor": len(t), "pasa": len(t) == 12 and len(union) == len(vieja)},
        {"prueba": "max |dif| de predicciones frente a la CV anterior (tolerancia 1e-6, orden de suma de los hilos)", "valor": t.max_dif_pred.max(), "pasa": t.max_dif_pred.max() < 1e-6},
        {"prueba": "misma epoca optima en los 12 ajustes", "valor": int((t.epoca_nueva == t.epoca_vieja).sum()), "pasa": (t.epoca_nueva == t.epoca_vieja).all()},
        {"prueba": "L2 relativo CV (antes 0.4397)", "valor": cv, "pasa": abs(cv - 0.4397) < 5e-5},
        {"prueba": "script oficial de Federico, semilla 2026: L2 0.417875", "valor": fed.L2_rel_nuevo, "pasa": abs(fed.L2_rel_nuevo - 0.417875) < 5e-7},
        {"prueba": "script oficial de Federico, semilla 2026: epoca 107", "valor": fed.epoca_nueva, "pasa": fed.epoca_nueva == 107},
    ])
    t.to_csv(SALIDA / "paso0_equivalencia_por_ajuste.csv", index=False)
    pruebas.to_csv(SALIDA / "paso0_equivalencia.csv", index=False)
    pd.set_option("display.width", 200)
    print(t.round(8).to_string(index=False))
    print(pruebas.to_string(index=False))


if __name__ == "__main__":
    if sys.argv[1] == "verificar":
        verificar()
    else:
        correr(int(sys.argv[2]), int(sys.argv[3]), sys.argv[4:] or list(VARIANTES))

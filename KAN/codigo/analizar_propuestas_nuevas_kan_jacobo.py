"""Compara las propuestas KAN nuevas con la de Federico, el MLP igualado y las referencias (validacion cruzada por bloques)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

import kan_variantes_core as core

BASE = core.ROOT / "KAN" / "salidas"
FED = "KAN [10,4,6] G=5 D2=100 (Federico)"
DIAS = np.arange(1, 366, dtype=float)
REFERENCIAS = [FED, "MLP oculta 37", "KAN [10,4,6] G=5 D2=0 (base)", "Lineal (mismo codificador)"]

pares = pd.read_csv(core.INPUT / "pares.csv", dtype={"codigo": str})
coef_mm = pd.read_csv(core.INPUT / "coef_estacion_mm.csv").to_numpy(float)
Phi = core._basis_at_nodes(DIAS, core.DATA["tau_e"])
cols = [f"c{i + 1}" for i in range(core.K_E)]

viejas = pd.concat([pd.read_csv(f, dtype={"codigo": str}) for f in sorted((BASE / "cv_bloques").glob("cv*.csv"))], ignore_index=True)
nuevas = pd.concat([pd.read_csv(f, dtype={"codigo": str}) for f in sorted((BASE / "cv_nuevas").glob("nv*.csv"))], ignore_index=True)
cv = pd.concat([viejas[viejas.modelo.isin(REFERENCIAS)], nuevas], ignore_index=True)
assert cv.groupby(["modelo", "bloque"])["semilla"].nunique().eq(3).all()
assert cv.groupby("modelo")["bloque"].nunique().eq(4).all()


def l2_por_estacion(g):
    filas = g["fila"].to_numpy() - 1
    obs = np.maximum(coef_mm[filas] @ Phi.T, 0)
    pred = np.maximum(g[cols].to_numpy(float) @ Phi.T, 0) ** 2
    l2 = np.sqrt(((pred - obs) ** 2).sum(1) / (obs ** 2).sum(1))
    return pd.Series(l2).groupby(pares["codigo"].to_numpy()[filas]).mean()


filas, por_est, ens_est = [], {}, {}
for (modelo, bloque, semilla), g in cv.groupby(["modelo", "bloque", "semilla"]):
    est = l2_por_estacion(g)
    por_est[(modelo, bloque, semilla)] = est
    filas.append({"modelo": modelo, "bloque": bloque, "semilla": semilla, "L2": est.mean(), "epoca_optima": g["epoca_optima"].iloc[0],
                  "parametros_totales": g["parametros_totales"].iloc[0],
                  "pct_fuera": g["pct_ocultas_fuera_rejilla"].iloc[0] if "pct_ocultas_fuera_rejilla" in g else np.nan})
ps = pd.DataFrame(filas)
ens = []
for (modelo, bloque), g in cv.groupby(["modelo", "bloque"]):
    e = l2_por_estacion(g.groupby("fila", as_index=False)[cols].mean())
    ens_est[(modelo, bloque)] = e
    ens.append({"modelo": modelo, "bloque": bloque, "L2_ensamble": e.mean()})
ens = pd.DataFrame(ens)

tabla = ps.groupby(["modelo", "bloque"])["L2"].mean().unstack()
res = ps.groupby(["modelo", "semilla"])["L2"].mean().groupby("modelo").agg(L2_CV="mean", DE_semillas="std").reset_index()
res = res.merge(ens.groupby("modelo")["L2_ensamble"].mean().rename("L2_CV_ensamble").reset_index(), on="modelo")
res = res.merge(ps.groupby("modelo").agg(parametros=("parametros_totales", "first"), epoca_media=("epoca_optima", "mean"),
                                         pct_ocultas_fuera_rejilla=("pct_fuera", "mean")).reset_index(), on="modelo")
res = res.merge(tabla.reset_index(), on="modelo")

contraste = []
for modelo in res.modelo:
    if modelo == FED:
        continue
    d = []
    for b in tabla.columns:
        a = pd.concat([por_est[(FED, b, s)] for s in core.SEEDS], axis=1).mean(axis=1)
        m = pd.concat([por_est[(modelo, b, s)] for s in core.SEEDS], axis=1).mean(axis=1)
        d.append((m - a).to_numpy())
    d = np.concatenate(d)
    boot = np.random.default_rng(2026).choice(d, (2000, len(d))).mean(1)
    contraste.append({"modelo": modelo, "dif_vs_Federico": d.mean(), "IC95_inf": np.quantile(boot, .025), "IC95_sup": np.quantile(boot, .975),
                      "pct_estacion_bloque_mejor": 100 * (d < 0).mean(), "p_wilcoxon": wilcoxon(d).pvalue,
                      "bloques_ganados": int((tabla.loc[modelo] < tabla.loc[FED]).sum())})
res = res.merge(pd.DataFrame(contraste), on="modelo", how="left").sort_values("L2_CV")
res.to_csv(BASE / "nuevas_resumen.csv", index=False)
ps.to_csv(BASE / "nuevas_L2_por_bloque_y_semilla.csv", index=False)

bloques = list(tabla.columns)
fig, axes = plt.subplots(1, 2, figsize=(17, 5.8), gridspec_kw={"width_ratios": [1.15, 1]})
ax = axes[0]
orden = res.reset_index(drop=True)
y = np.arange(len(orden))
color = ["#2f5fa8" if m == FED else ("#c0392b" if m.startswith("N") else "#888888") for m in orden.modelo]
ax.errorbar(orden.L2_CV, y, xerr=orden.DE_semillas, fmt="none", ecolor="#999999", capsize=2)
ax.scatter(orden.L2_CV, y, c=color, s=45, zorder=3, label="media de 3 semillas")
ax.scatter(orden.L2_CV_ensamble, y, marker="x", color="#228833", zorder=4, label="promedio de las 3 semillas (ensamble)")
ax.axvline(res.loc[res.modelo == FED, "L2_CV"].item(), color="#2f5fa8", ls=":")
ax.set_yticks(y, [f"{m} ({p})" for m, p in zip(orden.modelo, orden.parametros)], fontsize=8)
ax.invert_yaxis()
ax.set_xlabel("L2 relativo, media de los 4 bloques (menor es mejor)")
ax.set_title("Rojo: propuestas nuevas · Azul: Federico · Gris: referencias", fontsize=10)
ax.legend(fontsize=8, loc="lower right")
ax = axes[1]
for m in orden.modelo:
    destacado = m == FED or m == orden.modelo.iloc[0]
    ax.plot(bloques, [tabla.loc[m, b] for b in bloques], marker="o", lw=2.6 if destacado else 1, alpha=1 if destacado else .55,
            color="#2f5fa8" if m == FED else None, label=m)
ax.set_xlabel("bloque de años dejado fuera")
ax.set_ylabel("L2 relativo")
ax.set_title("Error por bloque (gruesas: Federico y la mejor)", fontsize=10)
ax.legend(fontsize=6.5)
fig.suptitle("Propuestas KAN nuevas frente a la de Federico: validación cruzada por bloques de 5 años (1981–2000), 3 semillas", fontsize=12)
fig.tight_layout()
fig.savefig(BASE / "fig4_propuestas_nuevas_cv.png", dpi=120)

pd.set_option("display.width", 260)
pd.set_option("display.max_columns", 30)
print(res.round(4).to_string(index=False))

"""Compara las propuestas de cabeza MLP de todas las rondas con la KAN de Federico (validacion cruzada por bloques).

Uso desde KAN/codigo:
    python -X utf8 analizar_mlp_jacobo.py [figura]
"""

from __future__ import annotations

import sys

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

import kan_variantes_core as core

BASE = core.ROOT / "KAN" / "salidas"
FED = "KAN [10,4,6] G=5 D2=100 (Federico)"
MLP = "MLP oculta 37"
ELEGIDO = "M25 MLP 3 capas 14-14-14 + L2=200"
DIAS = np.arange(1, 366, dtype=float)

pares = pd.read_csv(core.INPUT / "pares.csv", dtype={"codigo": str})
coef_mm = pd.read_csv(core.INPUT / "coef_estacion_mm.csv").to_numpy(float)
Phi = core._basis_at_nodes(DIAS, core.DATA["tau_e"])
cols = [f"c{i + 1}" for i in range(core.K_E)]

viejas = pd.concat([pd.read_csv(f) for f in sorted((BASE / "cv_bloques").glob("cv*.csv"))], ignore_index=True)
viejas = viejas[viejas.modelo.isin([FED, MLP, "Lineal (mismo codificador)"])].assign(ronda="ref")
nuevas = pd.concat([pd.read_csv(f) for f in sorted((BASE / "cv_mlp").glob("r*/m*.csv"))], ignore_index=True)
nuevas = nuevas[~nuevas.modelo.str.startswith("M00")]
cv = pd.concat([viejas, nuevas], ignore_index=True)
assert cv.groupby(["modelo", "bloque"])["semilla"].nunique().eq(3).all()
BLOQUES = sorted(cv.bloque.unique())


def l2_por_estacion(g):
    filas = g["fila"].to_numpy() - 1
    obs = np.maximum(coef_mm[filas] @ Phi.T, 0)
    pred = np.maximum(g[cols].to_numpy(float) @ Phi.T, 0) ** 2
    l2 = np.sqrt(((pred - obs) ** 2).sum(1) / (obs ** 2).sum(1))
    return pd.Series(l2).groupby(pares["codigo"].to_numpy()[filas]).mean()


por_est, ens_est, filas = {}, {}, []
for (modelo, bloque, semilla), g in cv.groupby(["modelo", "bloque", "semilla"]):
    est = l2_por_estacion(g)
    por_est[(modelo, bloque, semilla)] = est
    filas.append({"modelo": modelo, "ronda": g["ronda"].iloc[0], "bloque": bloque, "semilla": semilla, "L2": est.mean(),
                  "epoca_optima": g["epoca_optima"].iloc[0], "parametros": g["parametros_totales"].iloc[0]})
ps = pd.DataFrame(filas)
for (modelo, bloque), g in cv.groupby(["modelo", "bloque"]):
    ens_est[(modelo, bloque)] = l2_por_estacion(g.groupby("fila", as_index=False)[cols].mean())

tabla = ps.groupby(["modelo", "bloque"])["L2"].mean().unstack()
res = ps.groupby(["modelo", "semilla"])["L2"].mean().groupby("modelo").agg(L2_CV="mean", DE_semillas="std").reset_index()
res["L2_CV_ensamble"] = [np.mean([ens_est[(m, b)].mean() for b in BLOQUES]) for m in res.modelo]
res = res.merge(ps.groupby("modelo").agg(ronda=("ronda", "first"), parametros=("parametros", "first"), epoca_media=("epoca_optima", "mean")).reset_index(), on="modelo")
res = res.merge(tabla.reset_index(), on="modelo")


def contraste(modelo, ref, ensamble=False):
    if ensamble:
        d = np.concatenate([(ens_est[(modelo, b)] - ens_est[(ref, b)]).to_numpy() for b in BLOQUES])
    else:
        d = np.concatenate([(pd.concat([por_est[(modelo, b, s)] for s in core.SEEDS], axis=1).mean(axis=1)
                             - pd.concat([por_est[(ref, b, s)] for s in core.SEEDS], axis=1).mean(axis=1)).to_numpy() for b in BLOQUES])
    boot = np.random.default_rng(2026).choice(d, (2000, len(d))).mean(1)
    return d.mean(), np.quantile(boot, .025), np.quantile(boot, .975), 100 * (d < 0).mean(), wilcoxon(d).pvalue


extra = []
for m in res.modelo:
    fila = {"modelo": m}
    for ref, et in ((FED, "Fed"), (MLP, "MLP37")):
        if m != ref:
            dif, lo, hi, pct, p = contraste(m, ref)
            fila.update({f"dif_vs_{et}": dif, f"IC_inf_{et}": lo, f"IC_sup_{et}": hi, f"pct_mejor_{et}": pct, f"p_{et}": p,
                         f"bloques_ganados_{et}": int((tabla.loc[m] < tabla.loc[ref]).sum())})
    if m != FED:
        dif, lo, hi, pct, p = contraste(m, FED, ensamble=True)
        fila.update(dif_ens_vs_Fed=dif, IC_inf_ens=lo, IC_sup_ens=hi, p_ens=p)
    extra.append(fila)
res = res.merge(pd.DataFrame(extra), on="modelo").sort_values("L2_CV")
res.to_csv(BASE / "mlp_resumen.csv", index=False)
ps.to_csv(BASE / "mlp_L2_por_bloque_y_semilla.csv", index=False)

pd.set_option("display.width", 300)
pd.set_option("display.max_columns", 40)
print(res[["modelo", "parametros", "L2_CV", "DE_semillas", "L2_CV_ensamble", "epoca_media", *BLOQUES, "dif_vs_Fed", "IC_inf_Fed", "IC_sup_Fed",
           "p_Fed", "bloques_ganados_Fed", "dif_vs_MLP37", "p_MLP37", "dif_ens_vs_Fed", "p_ens"]].round(4).to_string(index=False))

if len(sys.argv) > 1:
    finalistas = [m for m in res.modelo if m.split()[0] in ("M43", "M27", "M25", "M39", "M33")]
    pares_f = []
    for a in finalistas:
        for b in finalistas:
            if a < b:
                dif, lo, hi, pct, p = contraste(a, b)
                de, le, he, _, pe = contraste(a, b, ensamble=True)
                pares_f.append({"A": a.split()[0], "B": b.split()[0], "dif_A_menos_B": dif, "IC_inf": lo, "IC_sup": hi, "p": p,
                                "dif_ensamble": de, "IC_inf_ens": le, "IC_sup_ens": he, "p_ens": pe})
    pares_f = pd.DataFrame(pares_f)
    pares_f.to_csv(BASE / "mlp_finalistas_pares.csv", index=False)
    print(pares_f.round(4).to_string(index=False))
    orden = res.reset_index(drop=True)
    fig, axes = plt.subplots(1, 2, figsize=(20, max(6, .32 * len(orden) + 2)), gridspec_kw={"width_ratios": [1.2, 1]})
    ax = axes[0]
    y = np.arange(len(orden))
    color = ["#2f5fa8" if m == FED else ("#888888" if r == "ref" else "#c0392b") for m, r in zip(orden.modelo, orden.ronda)]
    ax.errorbar(orden.L2_CV, y, xerr=orden.DE_semillas, fmt="none", ecolor="#999999", capsize=2)
    ax.scatter(orden.L2_CV, y, c=color, s=40, zorder=3, label="media de 3 semillas")
    ax.scatter(orden.L2_CV_ensamble, y, marker="x", color="#228833", zorder=4, label="promedio de las 3 semillas (ensamble)")
    ax.axvline(res.loc[res.modelo == FED, "L2_CV"].item(), color="#2f5fa8", ls=":")
    ax.set_yticks(y, [f"{m} ({p})" for m, p in zip(orden.modelo, orden.parametros)], fontsize=7.5)
    for etiqueta, m in zip(ax.get_yticklabels(), orden.modelo):
        if m in (FED, ELEGIDO):
            etiqueta.set_fontweight("bold")
    ax.set_xlim(.405, .465)
    ax.invert_yaxis()
    ax.set_xlabel("L2 relativo, media de los 4 bloques (menor es mejor)")
    ax.set_title("Rojo: propuestas MLP · Azul: KAN de Federico · Gris: referencias", fontsize=10)
    ax.legend(fontsize=8, loc="lower right")
    ax = axes[1]
    top = {FED: ("#2f5fa8", 2.6), MLP: ("#888888", 1.4), ELEGIDO: ("#c0392b", 2.6)}
    top.update({m: (c, 1.2) for m, c in zip([m for m in orden.modelo if m.split()[0] in ("M27", "M43", "M33")], ("#e08a3c", "#8e44ad", "#27ae60"))})
    for m, (c, lw) in top.items():
        ax.plot(BLOQUES, [tabla.loc[m, b] for b in BLOQUES], marker="o", lw=lw, color=c, label=m)
    ax.set_xlabel("bloque de años dejado fuera")
    ax.set_ylabel("L2 relativo")
    ax.set_title("Error por bloque: KAN de Federico, MLP37, la elegida (M25, roja) y otras finalistas", fontsize=10)
    ax.legend(fontsize=7.5)
    fig.suptitle("Propuestas de cabeza MLP frente a la KAN de Federico: validación cruzada por bloques de 5 años (1981–2000), 3 semillas", fontsize=12, y=.995)
    fig.subplots_adjust(left=.2, right=.98, top=.94, bottom=.05, wspace=.12)
    fig.savefig(BASE / "fig7_mlp_cv.png", dpi=120)

"""Ablacion de la cabeza con la KAN oficial: KAN, MLP y lineal sin regularizacion en la cabeza.

La KAN oficial (entrenar_kan_base_oficial.py) es la KAN [10,4,6] G=5 sin penalizacion D2; coincide
con la "KAN base" de la validacion cruzada por bloques. Aqui solo entran cabezas sin regularizacion
(ni L2 ni dropout), para que la diferencia mida la arquitectura y no la regularizacion. Las cabezas
regularizadas se reportan aparte como contexto.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

import kan_variantes_core as core
from cv_mlp_jacobo import RONDAS, cfg_completa

BASE = core.ROOT / "KAN" / "salidas"
KAN = "KAN oficial [10,4,6] G=5 (sin regularización)"
LIN = "Lineal (regresión funcional)"
RENOMBRE = {"KAN [10,4,6] G=5 D2=0 (base)": KAN, "Lineal (mismo codificador)": LIN, "MLP oculta 37": "MLP37 (1 capa de 37)",
            "MLP oculta 34": "MLP34 (1 capa de 34)", "KAN [10,4,6] G=5 D2=100 (Federico)": "KAN + D²=100"}
DIAS = np.arange(1, 366, dtype=float)

pares = pd.read_csv(core.INPUT / "pares.csv", dtype={"codigo": str})
coef_mm = pd.read_csv(core.INPUT / "coef_estacion_mm.csv").to_numpy(float)
Phi = core._basis_at_nodes(DIAS, core.DATA["tau_e"])
cols = [f"c{i + 1}" for i in range(core.K_E)]

sin_reg = {n for r in RONDAS.values() for n, c in r.items() if cfg_completa(c)["l2"] == 0 and cfg_completa(c)["drop"] == 0 and not n.startswith("M00")}
viejas = pd.concat([pd.read_csv(f) for f in sorted((BASE / "cv_bloques").glob("cv*.csv"))], ignore_index=True)
viejas = viejas[viejas.modelo.isin(RENOMBRE)].assign(modelo=lambda d: d.modelo.map(RENOMBRE))
nuevas = pd.concat([pd.read_csv(f) for f in sorted((BASE / "cv_mlp").glob("r*/m*.csv"))], ignore_index=True)
nuevas = nuevas[~nuevas.modelo.str.startswith("M00")]
cv = pd.concat([viejas, nuevas], ignore_index=True)
assert cv.groupby(["modelo", "bloque"])["semilla"].nunique().eq(3).all()
BLOQUES = sorted(cv.bloque.unique())


def grupo(m):
    if m == KAN:
        return "KAN oficial"
    if m == LIN:
        return "lineal"
    if m in sin_reg or m.startswith("MLP3"):
        return "MLP sin regularización"
    return "con regularización (aparte)"


def l2_por_estacion(g):
    filas = g["fila"].to_numpy() - 1
    obs = np.maximum(coef_mm[filas] @ Phi.T, 0)
    pred = np.maximum(g[cols].to_numpy(float) @ Phi.T, 0) ** 2
    l2 = np.sqrt(((pred - obs) ** 2).sum(1) / (obs ** 2).sum(1))
    return pd.Series(l2).groupby(pares["codigo"].to_numpy()[filas]).mean()


por_est, ens_est, filas = {}, {}, []
for (m, b, s), g in cv.groupby(["modelo", "bloque", "semilla"]):
    por_est[(m, b, s)] = l2_por_estacion(g)
    filas.append({"modelo": m, "bloque": b, "semilla": s, "L2": por_est[(m, b, s)].mean(),
                  "epoca_optima": g["epoca_optima"].iloc[0], "parametros": g["parametros_totales"].iloc[0]})
for (m, b), g in cv.groupby(["modelo", "bloque"]):
    ens_est[(m, b)] = l2_por_estacion(g.groupby("fila", as_index=False)[cols].mean())
ps = pd.DataFrame(filas)
tabla = ps.groupby(["modelo", "bloque"])["L2"].mean().unstack()
res = ps.groupby(["modelo", "semilla"])["L2"].mean().groupby("modelo").agg(L2_CV="mean", DE_semillas="std").reset_index()
res["L2_CV_ensamble"] = [np.mean([ens_est[(m, b)].mean() for b in BLOQUES]) for m in res.modelo]
res["peor_bloque"] = [tabla.loc[m].max() for m in res.modelo]
res = res.merge(ps.groupby("modelo").agg(parametros=("parametros", "first"), epoca_media=("epoca_optima", "mean")).reset_index(), on="modelo")
res["grupo"] = res.modelo.map(grupo)
res = res.merge(tabla.reset_index(), on="modelo")

contr = []
for m in res.modelo:
    if m == KAN:
        continue
    d = np.concatenate([(pd.concat([por_est[(m, b, s)] for s in core.SEEDS], axis=1).mean(axis=1)
                         - pd.concat([por_est[(KAN, b, s)] for s in core.SEEDS], axis=1).mean(axis=1)).to_numpy() for b in BLOQUES])
    boot = np.random.default_rng(2026).choice(d, (2000, len(d))).mean(1)
    contr.append({"modelo": m, "dif_vs_KAN": d.mean(), "IC_inf": np.quantile(boot, .025), "IC_sup": np.quantile(boot, .975),
                  "pct_estacion_bloque_mejor": 100 * (d < 0).mean(), "p_wilcoxon": wilcoxon(d).pvalue,
                  "bloques_ganados": int((tabla.loc[m] < tabla.loc[KAN]).sum())})
res = res.merge(pd.DataFrame(contr), on="modelo", how="left")
res["dif_pct_parametros_vs_KAN"] = 100 * (res.parametros / res.loc[res.modelo == KAN, "parametros"].item() - 1)
res = res.sort_values(["grupo", "L2_CV"], key=lambda s: s.map({"lineal": 3, "KAN oficial": 1, "MLP sin regularización": 0, "con regularización (aparte)": 4}) if s.name == "grupo" else s)
res.to_csv(BASE / "ablacion_oficial_resumen.csv", index=False)

oficial = res[res.grupo != "con regularización (aparte)"]
mlp_sr = oficial[oficial.grupo == "MLP sin regularización"].sort_values("L2_CV")
mejor_mlp = mlp_sr.iloc[0]
l_lin, l_kan = (res.loc[res.modelo == x, "L2_CV"].item() for x in (LIN, KAN))
costo = {
    "L2_lineal": l_lin, "L2_KAN_oficial": l_kan, "L2_mejor_MLP_sin_reg": mejor_mlp.L2_CV, "mejor_MLP_sin_reg": mejor_mlp.modelo,
    "L2_mediana_MLP_sin_reg": mlp_sr.L2_CV.median(), "n_MLP_sin_reg": len(mlp_sr),
    "MLP_sin_reg_peores_que_KAN": int((mlp_sr.L2_CV > l_kan).sum()),
    "ganancia_no_lineal_KAN": l_lin - l_kan, "ganancia_no_lineal_mejor_MLP": l_lin - mejor_mlp.L2_CV,
    "pct_ganancia_capturada_por_KAN": 100 * (l_lin - l_kan) / (l_lin - mejor_mlp.L2_CV),
    "costo_interpretabilidad_abs": l_kan - mejor_mlp.L2_CV, "costo_interpretabilidad_pct": 100 * (l_kan - mejor_mlp.L2_CV) / l_kan,
}
pd.Series(costo).to_csv(BASE / "ablacion_oficial_costo_interpretabilidad.csv", header=["valor"])

CONTEXTO = ["KAN + D²=100", "M25 MLP 3 capas 14-14-14 + L2=200"]
fig, axes = plt.subplots(1, 2, figsize=(19, 8), gridspec_kw={"width_ratios": [1.1, 1]})
ax = axes[0]
orden = pd.concat([oficial, res[res.modelo.isin(CONTEXTO)]]).reset_index(drop=True)
y = np.arange(len(orden))
col = {"KAN oficial": "#2f5fa8", "MLP sin regularización": "#c0392b", "lineal": "#555555", "con regularización (aparte)": "#aaaaaa"}
ax.errorbar(orden.L2_CV, y, xerr=orden.DE_semillas, fmt="none", ecolor="#aaaaaa", capsize=2)
for g_, c in col.items():
    s = (orden.grupo == g_).to_numpy()
    ax.scatter(orden.L2_CV[s], y[s], color=c, s=50, zorder=3, marker="D" if g_.startswith("con") else "o",
               label=g_ if not g_.startswith("con") else "con regularización (solo contexto)")
ax.axvline(l_kan, color="#2f5fa8", ls=":")
ax.axhline(len(oficial) - .5, color="#888888", lw=.8)
ax.set_yticks(y, [f"{m} ({p})" for m, p in zip(orden.modelo, orden.parametros)], fontsize=9)
for et, m in zip(ax.get_yticklabels(), orden.modelo):
    if m == KAN:
        et.set_fontweight("bold")
ax.invert_yaxis()
ax.set_xlabel("L2 relativo, media de 4 bloques × 3 semillas (menor es mejor); barra = DE entre semillas")
ax.set_title("Ablación oficial: las 16 cabezas sin regularización. Debajo de la línea gris, 2 regularizadas como contexto", fontsize=10)
ax.legend(fontsize=8.5, loc="lower right")
ax = axes[1]
med = costo["L2_mediana_MLP_sin_reg"]
x = np.arange(4)
vals = [l_lin, l_kan, med, mejor_mlp.L2_CV]
ax.bar(x, vals, color=["#555555", "#2f5fa8", "#e8a09a", "#c0392b"], width=.6)
for xi, v in zip(x, vals):
    ax.text(xi, v + .0006, f"{v:.4f}", ha="center", fontsize=11)
ax.annotate("", xy=(1.42, l_kan), xytext=(1.42, l_lin), arrowprops=dict(arrowstyle="<->", color="#2f5fa8", lw=2))
ax.text(1.46, l_lin - .002, f"no linealidad que\ncaptura la KAN:\n{l_lin - l_kan:.4f}", color="#2f5fa8", fontsize=10, va="top")
ax.annotate("", xy=(2.42, med), xytext=(2.42, l_kan), arrowprops=dict(arrowstyle="<->", color="#c0392b", lw=2))
ax.text(2.46, l_kan + .006, f"costo de la\ninterpretabilidad\n(mediana): {l_kan - med:.4f}\n= {100 * (l_kan - med) / l_kan:.1f} % del error", color="#c0392b", fontsize=10, va="top")
ax.set_xticks(x, ["Lineal\n(regresión funcional)", "KAN oficial", "MLP sin reg.\n(mediana de 14)", f"Mejor MLP sin reg.\n({mejor_mlp.modelo.split()[0]}, DE {mejor_mlp.DE_semillas:.3f})"])
ax.set_ylim(.40, .47)
ax.set_ylabel("L2 relativo (validación cruzada)")
ax.set_title("Cuánto de la no linealidad del sesgo captura la KAN y cuánto cuesta su interpretabilidad", fontsize=10)
fig.suptitle("Ablación de la cabeza con la KAN oficial: KAN, MLP y lineal sin regularización, mismos parámetros (±10 %), 4 bloques × 3 semillas", fontsize=12, y=.995)
fig.subplots_adjust(left=.19, right=.98, top=.9, bottom=.12, wspace=.13)
fig.savefig(BASE / "fig9_ablacion_oficial.png", dpi=120)

pd.set_option("display.width", 300)
pd.set_option("display.max_columns", 30)
print(res[["modelo", "grupo", "parametros", "dif_pct_parametros_vs_KAN", "L2_CV", "DE_semillas", "L2_CV_ensamble", "peor_bloque", "epoca_media", *BLOQUES,
           "dif_vs_KAN", "IC_inf", "IC_sup", "pct_estacion_bloque_mejor", "p_wilcoxon", "bloques_ganados"]].round(4).to_string(index=False))
print(pd.Series(costo).to_string())

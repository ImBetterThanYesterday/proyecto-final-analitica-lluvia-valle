"""Resume la validacion cruzada por bloques con el L2 relativo por estacion (misma metrica de evaluar_parametros_kan.R)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

import kan_variantes_core as core

BASE = core.ROOT / "KAN" / "salidas"
A = "KAN [10,4,6] G=5 D2=100 (Federico)"
B = "KAN residual [10,4,6] G=3 D2=100"
DIAS = np.arange(1, 366, dtype=float)

pares = pd.read_csv(core.INPUT / "pares.csv", dtype={"codigo": str})
coef_mm = pd.read_csv(core.INPUT / "coef_estacion_mm.csv").to_numpy(float)
Phi = core._basis_at_nodes(DIAS, core.DATA["tau_e"])
cols = [f"c{i + 1}" for i in range(core.K_E)]

cv = pd.concat([pd.read_csv(f, dtype={"codigo": str}) for f in sorted((BASE / "cv_bloques").glob("cv*.csv"))], ignore_index=True)
assert cv.groupby(["modelo", "bloque"])["semilla"].nunique().eq(3).all()


def l2_por_estacion(g):
    filas = g["fila"].to_numpy() - 1
    obs = np.maximum(coef_mm[filas] @ Phi.T, 0)
    pred = np.maximum(g[cols].to_numpy(float) @ Phi.T, 0) ** 2
    l2 = np.sqrt(((pred - obs) ** 2).sum(1) / (obs ** 2).sum(1))
    return pd.Series(l2).groupby(pares["codigo"].to_numpy()[filas]).mean()


filas, por_est = [], {}
for (modelo, bloque, semilla), g in cv.groupby(["modelo", "bloque", "semilla"]):
    est = l2_por_estacion(g)
    por_est[(modelo, bloque, semilla)] = est
    filas.append({"modelo": modelo, "bloque": bloque, "semilla": semilla, "L2_relativo": est.mean(),
                  "epoca_optima": g["epoca_optima"].iloc[0], "parametros_totales": g["parametros_totales"].iloc[0]})
por_semilla = pd.DataFrame(filas)

ens = []
for (modelo, bloque), g in cv.groupby(["modelo", "bloque"]):
    media = g.groupby("fila", as_index=False)[cols].mean()
    ens.append({"modelo": modelo, "bloque": bloque, "L2_ensamble": l2_por_estacion(media).mean()})
ens = pd.DataFrame(ens)

por_bloque = por_semilla.groupby(["modelo", "bloque"]).agg(L2_media=("L2_relativo", "mean"), L2_DE=("L2_relativo", "std")).reset_index()
por_bloque = por_bloque.merge(ens, on=["modelo", "bloque"])
tabla_bloques = por_bloque.pivot(index="modelo", columns="bloque", values="L2_media")
cv_semilla = por_semilla.groupby(["modelo", "semilla"])["L2_relativo"].mean().reset_index()
resumen = cv_semilla.groupby("modelo").agg(L2_CV_media=("L2_relativo", "mean"), L2_CV_DE=("L2_relativo", "std")).reset_index()
resumen = resumen.merge(ens.groupby("modelo")["L2_ensamble"].mean().rename("L2_CV_ensamble").reset_index(), on="modelo")
resumen = resumen.merge(por_semilla.groupby("modelo")["parametros_totales"].first().reset_index(), on="modelo")
resumen["bloques_ganados_vs_Federico"] = resumen["modelo"].map(
    lambda m: int((tabla_bloques.loc[m] < tabla_bloques.loc[A]).sum()) if m != A else np.nan)
resumen = resumen.merge(tabla_bloques.reset_index(), on="modelo").sort_values("L2_CV_media")

d = []
for bloque in sorted(cv["bloque"].unique()):
    ea = pd.concat([por_est[(A, bloque, s)] for s in core.SEEDS], axis=1).mean(axis=1)
    eb = pd.concat([por_est[(B, bloque, s)] for s in core.SEEDS], axis=1).mean(axis=1)
    d.append(pd.DataFrame({"bloque": bloque, "codigo": ea.index, "dif_B_menos_A": (eb - ea).to_numpy()}))
d = pd.concat(d)
rng = np.random.default_rng(2026)
boot = np.array([rng.choice(d["dif_B_menos_A"].to_numpy(), len(d)).mean() for _ in range(2000)])
contraste = pd.DataFrame([{
    "comparacion": f"{B} menos {A}", "pares_estacion_bloque": len(d), "dif_media": d["dif_B_menos_A"].mean(),
    "IC95_inf": np.quantile(boot, .025), "IC95_sup": np.quantile(boot, .975),
    "pct_B_mejor": 100 * (d["dif_B_menos_A"] < 0).mean(), "p_wilcoxon": wilcoxon(d["dif_B_menos_A"]).pvalue}])

por_semilla.to_csv(BASE / "cv_L2_por_bloque_y_semilla.csv", index=False)
por_bloque.to_csv(BASE / "cv_L2_por_bloque.csv", index=False)
resumen.to_csv(BASE / "cv_resumen.csv", index=False)
contraste.to_csv(BASE / "cv_contraste_finalistas.csv", index=False)

bloques = sorted(cv["bloque"].unique())
fig, ax = plt.subplots(figsize=(9, 4.8))
for m in resumen["modelo"]:
    ax.plot(bloques, [tabla_bloques.loc[m, b] for b in bloques], marker="o", label=m,
            lw=2.5 if m in (A, B) else 1.2, alpha=1 if m in (A, B) else .7)
ax.set_ylabel("L2 relativo (media por estación, 3 semillas)")
ax.set_xlabel("Bloque de años dejado fuera")
ax.set_title("Validación cruzada por bloques de 5 años (1981–2000)")
ax.legend(fontsize=7)
fig.tight_layout()
fig.savefig(BASE / "fig3_cv_bloques.png", dpi=120)

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 20)
print(resumen.round(4).to_string(index=False))
print(contraste.round(4).to_string(index=False))
print(por_bloque.round(4).to_string(index=False))

"""Resume los experimentos KAN con el L2 relativo por estacion en validacion (misma metrica de evaluar_parametros_kan.R)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

import kan_variantes_core as core

BASE = core.ROOT / "KAN" / "salidas"
REFERENCIA = "KAN [10,4,6] G=5 D2=0"
DIAS = np.arange(1, 366, dtype=float)

pares = pd.read_csv(core.INPUT / "pares.csv", dtype={"codigo": str})
coef_mm = pd.read_csv(core.INPUT / "coef_estacion_mm.csv").to_numpy(float)
Phi = core._basis_at_nodes(DIAS, core.DATA["tau_e"])
filas_val = np.flatnonzero(pares["particion"].to_numpy() == "validación")
obs = np.maximum(coef_mm[filas_val] @ Phi.T, 0)
codigos = pares["codigo"].to_numpy()[filas_val]
anios = pares["anio"].to_numpy()[filas_val]
fase = pares["fase"].to_numpy()[filas_val]
cols = [f"c{i + 1}" for i in range(core.K_E)]


def l2_curva(coef_raiz):
    pred = np.maximum(coef_raiz @ Phi.T, 0) ** 2
    return np.sqrt(((pred - obs) ** 2).sum(1) / (obs ** 2).sum(1))


def por_estacion(l2):
    return pd.Series(l2).groupby(codigos).mean()


corridas = pd.concat([pd.read_csv(f, dtype={"codigo": str}) for f in sorted((BASE / "corridas").glob("t*.csv"))], ignore_index=True)
assert (corridas.groupby(["modelo", "semilla"]).size() == len(filas_val)).all()
assert (corridas.groupby(["modelo", "semilla"])["fila"].apply(lambda s: np.array_equal(s.to_numpy() - 1, filas_val))).all()

filas, curvas, ensamble = [], {}, {}
for (modelo, semilla), g in corridas.groupby(["modelo", "semilla"], sort=False):
    l2 = l2_curva(g[cols].to_numpy(float))
    curvas[(modelo, semilla)] = l2
    primera = g.iloc[0]
    filas.append({"modelo": modelo, "grupo": primera["grupo"], "semilla": semilla, "L2_relativo": por_estacion(l2).mean(),
                  "parametros_cabeza": primera["parametros_cabeza"], "parametros_totales": primera["parametros_totales"],
                  "epoca_optima": primera["epoca_optima"], "segundos": primera["segundos"]})
por_semilla = pd.DataFrame(filas)
for modelo, g in corridas.groupby("modelo", sort=False):
    media = g.groupby("fila")[cols].mean().to_numpy(float)
    ensamble[modelo] = l2_curva(media)

resumen = por_semilla.groupby(["grupo", "modelo"], sort=False).agg(
    L2_media=("L2_relativo", "mean"), L2_DE=("L2_relativo", "std"), parametros_cabeza=("parametros_cabeza", "first"),
    parametros_totales=("parametros_totales", "first"), epoca_optima_media=("epoca_optima", "mean")).reset_index()
resumen["L2_ensamble_3_semillas"] = resumen["modelo"].map(lambda m: por_estacion(ensamble[m]).mean())

ref_est = pd.concat([por_estacion(curvas[(REFERENCIA, s)]) for s in core.SEEDS], axis=1).mean(axis=1)
comparacion = []
for modelo in resumen["modelo"]:
    est = pd.concat([por_estacion(curvas[(modelo, s)]) for s in core.SEEDS], axis=1).mean(axis=1)
    d = (est - ref_est).to_numpy()
    rng = np.random.default_rng(2026)
    boot = np.array([rng.choice(d, len(d)).mean() for _ in range(2000)])
    p = wilcoxon(d).pvalue if np.any(d != 0) else np.nan
    comparacion.append({"modelo": modelo, "dif_vs_referencia": d.mean(), "IC95_inf": np.quantile(boot, .025),
                        "IC95_sup": np.quantile(boot, .975), "pct_estaciones_mejor": 100 * (d < 0).mean(), "p_wilcoxon": p})
resumen = resumen.merge(pd.DataFrame(comparacion), on="modelo")

anual = []
for modelo in resumen["modelo"]:
    l2 = np.mean([curvas[(modelo, s)] for s in core.SEEDS], axis=0)
    for a in sorted(set(anios)):
        m = anios == a
        anual.append({"modelo": modelo, "anio": a, "fase": fase[m][0], "L2_relativo": pd.Series(l2[m]).groupby(codigos[m]).mean().mean()})
anual = pd.DataFrame(anual).pivot(index="modelo", columns="anio", values="L2_relativo").reset_index()
fases = {a: fase[anios == a][0] for a in sorted(set(anios))}

BASE.mkdir(parents=True, exist_ok=True)
por_semilla.to_csv(BASE / "L2_relativo_por_semilla.csv", index=False)
resumen.to_csv(BASE / "resumen_experimentos.csv", index=False)
anual.to_csv(BASE / "L2_relativo_por_anio_validacion.csv", index=False)

orden = resumen.sort_values("L2_media")
fig, ax = plt.subplots(figsize=(9, 0.32 * len(orden) + 1.2))
colores = {"A": "#3366aa", "B": "#cc6633", "R": "#777777"}
y = np.arange(len(orden))
ax.errorbar(orden["L2_media"], y, xerr=orden["L2_DE"], fmt="o", color="k", ecolor="#999999", capsize=2, zorder=1)
ax.scatter(orden["L2_media"], y, c=[colores[g] for g in orden["grupo"]], zorder=2)
ax.scatter(orden["L2_ensamble_3_semillas"], y, marker="x", color="#228833", zorder=3, label="Promedio de las 3 semillas")
ax.set_yticks(y, [f"{m} ({p})" for m, p in zip(orden["modelo"], orden["parametros_totales"])], fontsize=8)
ax.axvline(resumen.loc[resumen.modelo == REFERENCIA, "L2_media"].item(), ls=":", color="#3366aa")
ax.invert_yaxis()
ax.set_xlabel("L2 relativo en validación (media por estación; media ± DE de 3 semillas)")
ax.set_title("Azul: KAN (A)   Naranja: KAN residual (B)   Gris: referencias", fontsize=10)
ax.legend(loc="lower right", fontsize=8)
fig.tight_layout()
fig.savefig(BASE / "fig1_L2_validacion.png", dpi=120)

top = ["MLP oculta 37", "KAN residual [10,4,6] G=3 D2=100", "KAN [10,4,6] G=5 D2=100", REFERENCIA, "Lineal (mismo codificador)"]
top = list(dict.fromkeys(top))
fig, ax = plt.subplots(figsize=(9, 4.5))
for m in top:
    fila = anual[anual.modelo == m].iloc[0]
    ax.plot(list(fases), [fila[a] for a in fases], marker="o", label=m)
ax.set_xticks(list(fases), [f"{a}\n{f}" for a, f in fases.items()])
ax.set_ylabel("L2 relativo (media por estación)")
ax.set_title("Validación por año y fase ENSO")
ax.legend(fontsize=7)
fig.tight_layout()
fig.savefig(BASE / "fig2_L2_por_anio.png", dpi=120)

pd.set_option("display.width", 220)
pd.set_option("display.max_columns", 20)
print(orden[["grupo", "modelo", "parametros_totales", "L2_media", "L2_DE", "L2_ensamble_3_semillas", "dif_vs_referencia",
             "IC95_inf", "IC95_sup", "pct_estaciones_mejor", "p_wilcoxon", "epoca_optima_media"]].round(4).to_string(index=False))
print(anual.round(3).to_string(index=False))
print(fases)

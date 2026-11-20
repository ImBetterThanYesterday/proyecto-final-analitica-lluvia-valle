"""Analisis de la auditoria de parametros de la KAN oficial (validacion cruzada por bloques, 1981-2000).

Para cada parametro compara cada valor probado con el oficial: L2 relativo CV, DE entre semillas,
diferencia pareada por estacion-bloque con IC bootstrap del 95 % y Wilcoxon, bloques ganados,
parametros, epocas e indicadores de interpretabilidad. Regla fijada antes de correr: un valor
alternativo solo es candidato a cambio si el IC de su diferencia frente al oficial excluye el cero
a su favor; si empata, se mantiene el oficial o el mas simple.

Uso desde KAN/codigo:
    python -X utf8 analizar_auditoria_kan.py
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

import kan_variantes_core as core
from auditoria_parametros_kan import OFICIAL, SALIDA, VARIANTES, cfg_de

DIAS = np.arange(1, 366, dtype=float)
COLS = [f"c{i + 1}" for i in range(core.K_E)]
CLAVE = {"ancho oculto": "ocultas", "profundidad": "ocultas", "intervalos G": "G", "grado B-spline": "grado",
         "rango de la rejilla": "rango", "estandarizar z": "estandarizar", "regularización D2": "d2",
         "tasa de aprendizaje": "lr", "paciencia": "paciencia", "tamaño de lote": "lote",
         "neuronas funcionales": "neuronas", "componentes FPCA M": "M", "rugosidad lambda_w": "lambda_w"}
ORIGEN = {"neuronas funcionales": "Rango del enunciado (6 a 8)", "componentes FPCA M": "Rango del enunciado (4 a 6)"}
ORDEN = list(CLAVE)


def texto(param, cfg):
    v = cfg[CLAVE[param]]
    if param == "ancho oculto":
        return str(v[0])
    if param == "profundidad":
        return "[" + ",".join(map(str, (cfg["neuronas"] + 2, *v, cfg["M"]))) + "]"
    if param == "rango de la rejilla":
        return f"[{v[0]:g},{v[1]:g}]"
    if param == "estandarizar z":
        return "sí" if v else "no"
    if param == "tamaño de lote":
        return "completo" if v == 0 else str(v)
    return f"{v:g}"


def orden_valor(clave, v):
    if clave == "ocultas":
        return 100 * len(v) + (v[0] if v else 0)
    if clave == "rango":
        return v[1]
    return float(v)


pares = pd.read_csv(core.INPUT / "pares.csv", dtype={"codigo": str})
coef_mm = pd.read_csv(core.INPUT / "coef_estacion_mm.csv").to_numpy(float)
Phi = core._basis_at_nodes(DIAS, core.DATA["tau_e"])


def l2_por_estacion(g):
    filas = g["fila"].to_numpy() - 1
    obs = np.maximum(coef_mm[filas] @ Phi.T, 0)
    pred = np.maximum(g[COLS].to_numpy(float) @ Phi.T, 0) ** 2
    l2 = np.sqrt(((pred - obs) ** 2).sum(1) / (obs ** 2).sum(1))
    return pd.Series(l2).groupby(pares["codigo"].to_numpy()[filas]).mean()


cv = pd.concat([pd.read_csv(f) for f in sorted((SALIDA / "corridas").glob("*.csv"))], ignore_index=True)
completos = cv.groupby("codigo").apply(lambda d: d[["bloque", "semilla"]].drop_duplicates().shape[0]).eq(12)
cv = cv[cv.codigo.isin(completos[completos].index)]
BLOQUES = sorted(cv.bloque.unique())
diag = [c for c in cv.columns if c.startswith(("fuera_rejilla", "no_linealidad"))]

por_est, filas = {}, []
for (c, b, s), g in cv.groupby(["codigo", "bloque", "semilla"]):
    por_est[(c, b, s)] = l2_por_estacion(g)
    filas.append({"codigo": c, "bloque": b, "semilla": s, "L2": por_est[(c, b, s)].mean(), **g.iloc[0][["epoca_optima", "epocas_corridas", "parametros_totales", *diag]].to_dict()})
ps = pd.DataFrame(filas)
ps.to_csv(SALIDA / "auditoria_L2_por_bloque_y_semilla.csv", index=False)
tabla = ps.groupby(["codigo", "bloque"])["L2"].mean().unstack()

res = ps.groupby(["codigo", "semilla"])["L2"].mean().groupby("codigo").agg(L2_CV="mean", DE_semillas="std").reset_index()
extra = ps.groupby("codigo").agg(parametros=("parametros_totales", "first"), epoca_media=("epoca_optima", "mean"),
                                 pct_tope_epocas=("epocas_corridas", lambda e: 100 * (e >= OFICIAL["max_epocas"]).mean()),
                                 **{d: (d, "mean") for d in diag}).reset_index()
res = res.merge(extra, on="codigo").merge(tabla.reset_index(), on="codigo")


def media_semillas(c, b):
    return pd.concat([por_est[(c, b, s)] for s in core.SEEDS], axis=1).mean(axis=1)


contr = []
for c in res.codigo:
    if c == "oficial":
        continue
    d = np.concatenate([(media_semillas(c, b) - media_semillas("oficial", b)).to_numpy() for b in BLOQUES])
    boot = np.random.default_rng(2026).choice(d, (2000, len(d))).mean(1)
    contr.append({"codigo": c, "dif_vs_oficial": d.mean(), "IC_inf": np.quantile(boot, .025), "IC_sup": np.quantile(boot, .975),
                  "p_wilcoxon": wilcoxon(d).pvalue, "bloques_ganados": int((tabla.loc[c] < tabla.loc["oficial"]).sum())})
res = res.merge(pd.DataFrame(contr), on="codigo", how="left")
res["parametro"] = res.codigo.map(lambda c: VARIANTES[c][0])
res["valor"] = [texto(p, cfg_de(c)) if c != "oficial" else "oficial" for c, p in zip(res.codigo, res.parametro)]
res["veredicto"] = np.select([res.codigo == "oficial", res.IC_sup < 0, res.IC_inf > 0],
                             ["referencia", "mejor que el oficial (IC < 0)", "peor que el oficial (IC > 0)"], "empate (IC incluye 0)")
res["orden"] = res.parametro.map({p: i for i, p in enumerate(["configuración oficial", *ORDEN])})
res = res.sort_values(["orden", "L2_CV"]).drop(columns="orden")
res.to_csv(SALIDA / "auditoria_resumen.csv", index=False)

of = res[res.codigo == "oficial"].iloc[0]
maestra = []
for p in ORDEN:
    sub = res[res.parametro == p]
    if sub.empty:
        continue
    mejores = sub[sub.IC_sup < 0]
    peores = sub[sub.IC_inf > 0]
    maestra.append({
        "parametro": p, "valor_oficial": texto(p, OFICIAL), "origen": ORIGEN.get(p, "Decisión del grupo"),
        "valores_probados": ", ".join(f"{v} ({l:.4f})" for v, l in zip(sub.valor, sub.L2_CV)),
        "L2_CV_oficial": of.L2_CV,
        "mejor_alternativa": sub.iloc[0].valor, "dif_mejor_alternativa": sub.iloc[0].dif_vs_oficial,
        "IC_mejor_alternativa": f"[{sub.iloc[0].IC_inf:+.4f}; {sub.iloc[0].IC_sup:+.4f}]",
        "alternativas_significativamente_mejores": ", ".join(mejores.valor) or "ninguna",
        "alternativas_significativamente_peores": ", ".join(peores.valor) or "ninguna",
        "veredicto": "revisar cambio" if len(mejores) else "se mantiene el valor oficial",
    })
maestra = pd.DataFrame(maestra)
maestra.to_csv(SALIDA / "auditoria_tabla_maestra.csv", index=False)

params = [p for p in ORDEN if (res.parametro == p).any()]
nc = 4
nf = int(np.ceil(len(params) / nc))
fig, axes = plt.subplots(nf, nc, figsize=(18, 3.6 * nf), squeeze=False)
color = {"referencia": "#2f5fa8", "mejor que el oficial (IC < 0)": "#1e8449", "peor que el oficial (IC > 0)": "#c0392b", "empate (IC incluye 0)": "#999999"}
for ax, p in zip(axes.ravel(), params):
    sub = pd.concat([res[res.codigo == "oficial"].assign(valor=texto(p, OFICIAL) + "\n(oficial)"), res[res.parametro == p]])
    sub = sub.assign(_o=[orden_valor(CLAVE[p], cfg_de(c)[CLAVE[p]]) for c in sub.codigo]).sort_values("_o")
    x = np.arange(len(sub))
    ax.bar(x, sub.L2_CV, yerr=sub.DE_semillas, color=[color[v] for v in sub.veredicto], capsize=3, width=.65)
    ax.axhline(of.L2_CV, color="#2f5fa8", ls=":", lw=1)
    ax.set_xticks(x, sub.valor, fontsize=8)
    lo, hi = sub.L2_CV.min(), sub.L2_CV.max()
    ax.set_ylim(lo - .01, hi + .012)
    for xi, (l, dv) in enumerate(zip(sub.L2_CV, sub.dif_vs_oficial)):
        ax.text(xi, l + .002, f"{l:.4f}" + ("" if pd.isna(dv) else f"\n{dv:+.4f}"), ha="center", fontsize=7)
    ax.set_title(p, fontsize=10)
for ax in axes.ravel()[len(params):]:
    ax.axis("off")
handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in color.values()]
fig.legend(handles, ["oficial", "mejor (IC 95 % < 0)", "peor (IC 95 % > 0)", "empate (IC incluye 0)"], loc="lower right", ncol=4, fontsize=9)
fig.suptitle("Auditoría de parámetros de la KAN oficial: L2 relativo CV (4 bloques × 3 semillas), un parámetro a la vez; "
             "debajo, diferencia frente al oficial", fontsize=12)
fig.tight_layout(rect=(0, .03, 1, .97))
fig.savefig(SALIDA / "fig_auditoria_parametros.png", dpi=120)

pd.set_option("display.width", 320)
pd.set_option("display.max_columns", 40)
print(res[["parametro", "valor", "parametros", "L2_CV", "DE_semillas", *BLOQUES, "dif_vs_oficial", "IC_inf", "IC_sup", "p_wilcoxon",
           "bloques_ganados", "epoca_media", "pct_tope_epocas", *diag, "veredicto"]].round(4).to_string(index=False))
print(maestra.drop(columns="valores_probados").to_string(index=False))

"""Figuras y tablas para justificar las propuestas KAN nuevas.

1. Promedio de modelos (ensamble) de Federico con cada propuesta, con contraste pareado por estacion-bloque.
2. Uso de la rejilla [-3, 3] y rugosidad de las B-splines: Federico frente a N2 y N6 (bloque 1996-2000, semilla 2026).
3. B-splines al estilo de clase de la capa 2: Federico frente a N6.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import torch
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

import kan_variantes_core as core
from kan_variantes_core import EXTRA_T, X_T
import analizar_propuestas_nuevas_kan as an
from bsplines_kan import entradas_de_capas, entrenar_o_cargar
from bspline_estilo_clase_kan import panel
from cv_bloques_kan import ANIO, BLOQUES
from cv_propuestas_nuevas_kan import PROPUESTAS, construir
from interpretabilidad_kan import OCULTAS, SALIDAS

BASE = an.BASE
FED = an.FED
N2, N6, N7 = "N2 norm + D2=30", "N6 norm + G=6 + D2=30", "N7 D2=100 + lr 0,003 (sin norm)"
COMBOS = {
    "Federico (3 semillas)": [FED],
    "Federico + N7 (6 modelos, control)": [FED, N7],
    "Federico + MLP37 (6 modelos, control)": [FED, "MLP oculta 37"],
    "Federico + N2 (6 modelos)": [FED, N2],
    "Federico + N6 (6 modelos)": [FED, N6],
}


def ensamble_por_estacion(modelos):
    out = {}
    for b in BLOQUES:
        g = an.cv[(an.cv.bloque == b) & an.cv.modelo.isin(modelos)]
        out[b] = an.l2_por_estacion(g.groupby("fila", as_index=False)[an.cols].mean())
    return out


def ensambles():
    ref = ensamble_por_estacion([FED])
    filas = []
    for nombre, modelos in COMBOS.items():
        e = ensamble_por_estacion(modelos)
        d = np.concatenate([(e[b] - ref[b]).to_numpy() for b in BLOQUES])
        boot = np.random.default_rng(2026).choice(d, (2000, len(d))).mean(1)
        fila = {"combinacion": nombre, "parametros": sum(int(an.res.loc[an.res.modelo == m, "parametros"].item()) for m in modelos),
                "L2_CV": np.mean([e[b].mean() for b in BLOQUES]), **{b: e[b].mean() for b in BLOQUES}}
        if len(modelos) > 1:
            fila.update(dif_vs_Federico=d.mean(), IC95_inf=np.quantile(boot, .025), IC95_sup=np.quantile(boot, .975),
                        pct_estacion_bloque_mejor=100 * (d < 0).mean(), p_wilcoxon=wilcoxon(d).pvalue,
                        bloques_ganados=int(sum(e[b].mean() < ref[b].mean() for b in BLOQUES)))
        filas.append(fila)
    tabla = pd.DataFrame(filas)
    tabla.to_csv(BASE / "nuevas_ensambles.csv", index=False)
    fig, axes = plt.subplots(1, 2, figsize=(15, 5), gridspec_kw={"width_ratios": [1, 1.2]})
    colores = ["#2f5fa8", "#999999", "#bbbbbb", "#e08a3c", "#c0392b"]
    ax = axes[0]
    ax.barh(tabla.combinacion, tabla.L2_CV, color=colores)
    for y, v in enumerate(tabla.L2_CV):
        ax.text(v + .0003, y, f"{v:.4f}", va="center", fontsize=9)
    ax.set_xlim(tabla.L2_CV.min() - .006, tabla.L2_CV.max() + .004)
    ax.invert_yaxis()
    ax.set_xlabel("L2 relativo, media de los 4 bloques (menor es mejor)")
    ax.set_title("Promedio de las predicciones de varios modelos", fontsize=10)
    ax = axes[1]
    for (_, f), c in zip(tabla.iterrows(), colores):
        ax.plot(list(BLOQUES), [f[b] for b in BLOQUES], marker="o", color=c, lw=2.4 if c in ("#2f5fa8", "#c0392b") else 1.2, label=f.combinacion)
    ax.set_xlabel("bloque de años dejado fuera")
    ax.set_ylabel("L2 relativo")
    ax.set_title("Por bloque: N6 complementa a Federico en los 4 bloques", fontsize=10)
    ax.legend(fontsize=8)
    fig.suptitle("Combinar la KAN de Federico con una KAN normalizada: la ganancia viene de la diversidad, no de promediar más modelos", fontsize=12)
    fig.tight_layout()
    fig.savefig(BASE / "fig5_ensambles_cv.png", dpi=120)
    plt.close(fig)
    return tabla


def cargar_nueva(nombre, seed=2026):
    a, b = BLOQUES["1996-2000"]
    idx_ent = torch.as_tensor(np.flatnonzero(~((ANIO >= a) & (ANIO <= b))), dtype=torch.long)
    modelo = construir(PROPUESTAS[nombre], seed, idx_ent)
    modelo.load_state_dict(torch.load(BASE / "cv_nuevas" / "modelos" / f"{nombre.split()[0]}_semilla{seed}.pt"))
    modelo.eval()
    with torch.no_grad():
        x1 = modelo.head_input(X_T[idx_ent], EXTRA_T[idx_ent])
        x2 = modelo.ocultas(X_T[idx_ent], EXTRA_T[idx_ent])
    return modelo, x1.numpy(), x2.numpy()


def rugosidad(capa):
    with torch.no_grad():
        return float(torch.diff(capa.w_s.unsqueeze(-1) * capa.coef, n=2, dim=-1).square().mean())


def importancia(capa, datos):
    with torch.no_grad():
        return capa.activations(torch.as_tensor(datos)).numpy().std(0)


def uso_rejilla():
    fed = entrenar_o_cargar("federico")
    x1f, x2f = [d.numpy() for d in entradas_de_capas(fed)]
    modelos = {"Federico (D2=100, sin norm)": (fed, x1f, x2f)}
    for n in (N2, N6):
        modelos[n] = cargar_nueva(n)
    filas = []
    for nombre, (m, x1, x2) in modelos.items():
        c1, c2 = m.head.layers
        imp2 = importancia(c2, x2)
        fila = {"modelo": nombre, "pct_fuera_capa2": 100 * (np.abs(x2) > 3).mean(),
                **{f"pct_fuera_{h}": 100 * (np.abs(x2[:, k]) > 3).mean() for k, h in enumerate(OCULTAS)},
                "rango_h_min": x2.min(), "rango_h_max": x2.max(),
                "rugosidad_capa1": rugosidad(c1), "rugosidad_capa2": rugosidad(c2),
                "pct_importancia_max_capa2": 100 * imp2.max() / imp2.sum()}
        filas.append(fila)
    tabla = pd.DataFrame(filas)
    tabla.to_csv(BASE / "nuevas_uso_rejilla.csv", index=False)
    return modelos, tabla


def figura_capa2(modelos):
    elegidos = ["Federico (D2=100, sin norm)", N6]
    fig, axes = plt.subplots(2, 4, figsize=(19, 9))
    for fila, nombre in zip(axes, elegidos):
        m, _, x2 = modelos[nombre]
        capa2 = m.head.layers[1]
        imp = importancia(capa2, x2)
        orden = np.dstack(np.unravel_index(np.argsort(-imp, axis=None), imp.shape))[0][:4]
        for ax, (i, o) in zip(fila, orden):
            fuera = 100 * (np.abs(x2[:, i]) > 3).mean()
            panel(ax, capa2, i, o, x2[:, i], f"{nombre.split(' (')[0].split(' norm')[0]}: {OCULTAS[i]} → {SALIDAS[o]}\n"
                                               f"{100 * imp[i, o] / imp.sum():.0f} % de la capa; {fuera:.0f} % de datos fuera de [−3, 3]",
                  leyenda=ax is axes[0, 0])
            ax.hist(x2[:, i], bins=40, weights=np.full(len(x2), .25 * (ax.get_ylim()[1] - ax.get_ylim()[0]) / len(x2) * 8),
                    bottom=ax.get_ylim()[0], color="#9ecae1", alpha=.5, zorder=0)
    axes[0, 0].legend(fontsize=7, loc="best")
    fig.suptitle("Capa 2 (neuronas ocultas → puntajes FPCA), las 4 conexiones más importantes\n"
                 "Arriba: Federico, el 60 % de los datos de la conexión principal cae fuera de la rejilla y ahí solo actúa SiLU. "
                 "Abajo: N6, al estandarizar las neuronas casi todos los datos caen dentro y la importancia se reparte. Barras azules: dónde caen los datos", fontsize=11)
    fig.tight_layout()
    fig.savefig(BASE / "fig6_capa2_federico_vs_N6.png", dpi=115)
    plt.close(fig)


def main():
    torch.set_num_threads(4)
    pd.set_option("display.width", 260)
    pd.set_option("display.max_columns", 30)
    print(ensambles().round(4).to_string(index=False))
    modelos, tabla = uso_rejilla()
    print(tabla.round(4).to_string(index=False))
    figura_capa2(modelos)


if __name__ == "__main__":
    main()

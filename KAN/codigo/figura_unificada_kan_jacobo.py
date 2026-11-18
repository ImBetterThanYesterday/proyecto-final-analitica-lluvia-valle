"""Figura unificada de las funciones aprendidas en la capa 1: KAN base frente a KAN D2=100 (modelos guardados, semilla 2026)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import torch
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from bsplines_kan_jacobo import ENTRADAS, OCULTAS, SALIDA, capas_kan, entradas_de_capas, entrenar_o_cargar, funcion_arista

BLOQUES = [("base", "KAN base [10,4,6] G=5, sin D²", "#777777"),
           ("federico", "KAN [10,4,6] G=5, D²=100 (Federico)", "#2f5fa8")]


def main():
    rug = pd.read_csv(SALIDA / "rugosidad_por_capa.csv")
    aristas = pd.read_csv(SALIDA / "aristas_importancia_no_linealidad.csv")
    n_in, n_h = len(ENTRADAS), len(OCULTAS)
    fig = plt.figure(figsize=(17, 11.5))
    sub = fig.subfigures(2, 1, hspace=0.04)
    for b, (clave, titulo, color) in enumerate(BLOQUES):
        modelo = entrenar_o_cargar(clave)
        capa = capas_kan(modelo)[0]
        datos = entradas_de_capas(modelo)[0]
        r = rug[(rug.modelo == clave) & (rug.capa == 1)].iloc[0]
        nl = aristas[(aristas.modelo == clave) & (aristas.capa == 1)]
        sub[b].suptitle(f"{titulo}   |   rugosidad de las splines: {r.rugosidad_media_dif2:.3f}   |   no linealidad media: {r.no_linealidad_media:.2f}",
                        fontsize=12, color=color, fontweight="bold")
        axes = sub[b].subplots(n_h, n_in)
        for o in range(n_h):
            for i in range(n_in):
                d = datos[:, i].numpy()
                lo, hi = np.quantile(d, [.01, .99])
                x = np.linspace(lo, hi, 300)
                total, _, _ = funcion_arista(capa, i, o, x)
                ax = axes[o, i]
                ax.plot(x, total, color=color, lw=2)
                ax.axhline(0, color="#bbbbbb", lw=.6)
                ax.set_xticks([]); ax.set_yticks([])
                for s in ax.spines.values():
                    s.set_color("#cccccc")
                v = nl[(nl.entrada == ENTRADAS[i]) & (nl.salida == OCULTAS[o])].no_linealidad_1_menos_R2.iloc[0]
                ax.text(.97, .95, f"{v:.2f}", transform=ax.transAxes, ha="right", va="top", fontsize=7, color="#666666",
                        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=.85))
                if o == 0:
                    ax.set_title(ENTRADAS[i], fontsize=10)
                if i == 0:
                    ax.set_ylabel(OCULTAS[o], fontsize=10, rotation=0, labelpad=14)
    fig.suptitle("Funciones de activación aprendidas por la KAN (capa 1: 10 entradas → 4 neuronas ocultas), semilla 2026", fontsize=14, y=1.045)
    fig.text(.5, -.012,
             "Cada recuadro es una conexión: su curva φ(x) = w_b·SiLU(x) + w_s·Σ c_g B_g(x), dibujada en el rango p1–p99 de su entrada en entrenamiento. "
             "Columnas: entradas (z1–z8 del codificador funcional, altitud estandarizada y ONI). Filas: neuronas ocultas.\n"
             "Número gris = no linealidad (1 − R² de una recta; 0 = lineal). Las neuronas de los dos modelos no se corresponden entre sí: "
             "se compara la forma de las curvas, no conexión por conexión. D²=100 deja curvas suaves y monótonas o con un solo máximo; sin D² oscilan.",
             ha="center", va="top", fontsize=9, color="#444444")
    fig.savefig(SALIDA / "B2_activaciones_capa1_base_vs_D2.png", dpi=130, bbox_inches="tight")
    print(SALIDA / "B2_activaciones_capa1_base_vs_D2.png")


if __name__ == "__main__":
    torch.set_num_threads(4)
    main()

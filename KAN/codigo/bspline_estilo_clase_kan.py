"""B-splines de las conexiones KAN al estilo de clase: puntos de control, nodos y un polinomio por intervalo.

Modelo: el de interpretabilidad_kan.py segun el argumento federico u oficial, semilla 2026.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import torch
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from bsplines_kan import capas_kan, entradas_de_capas, entrenar_o_cargar
from interpretabilidad_kan import ENTRADAS, ETIQUETA, OCULTAS, SALIDAS, SALIDA, VERSION

COLORES = plt.get_cmap("tab10").colors


def partes(capa, i, o, x):
    with torch.no_grad():
        xx = torch.as_tensor(np.asarray(x, dtype=float)).reshape(-1, 1).repeat(1, capa.coef.shape[0])
        B = capa.bases(xx)[:, i, :].numpy()
        c = (capa.w_s[i, o] * capa.coef[i, o]).numpy()
        silu = float(capa.w_b[i, o]) * torch.nn.functional.silu(xx[:, i]).numpy()
    return B @ c, silu, c


def panel(ax, capa, i, o, datos, titulo, leyenda=False):
    t = capa.grid.numpy()
    orden = capa.order
    lo_r, hi_r = capa.range
    _, _, c = partes(capa, i, o, [0.0])
    greville = np.array([t[j + 1:j + orden + 1].mean() for j in range(len(c))])
    nodos = t[(t >= lo_r - 1e-9) & (t <= hi_r + 1e-9)]
    p1, p99 = np.quantile(datos, [.01, .99])
    ax.axvspan(max(p1, lo_r), min(p99, hi_r), color="#f3f3f3", zorder=0, label="rango del 98 % de los datos" if leyenda else None)
    ax.plot(greville, c, ":", color="#d62728", lw=1.2, zorder=2)
    ax.plot(greville, c, "o", color="#d62728", ms=5, zorder=3, label="puntos de control (w_s·c_g)" if leyenda else None)
    for k in range(len(nodos) - 1):
        x = np.linspace(nodos[k], nodos[k + 1], 80)
        s, _, _ = partes(capa, i, o, x)
        ax.plot(x, s, color=COLORES[k % 10], lw=2.4, zorder=4, label=f"polinomio {k}" if leyenda else None)
    s_n, _, _ = partes(capa, i, o, nodos)
    ax.plot(nodos, s_n, "o", color="#2ca02c", ms=8, zorder=5, label="nodos" if leyenda else None)
    x = np.linspace(min(lo_r, p1), max(hi_r, p99), 300)
    s, silu, _ = partes(capa, i, o, x)
    ax.plot(x, s + silu, "--", color="#555555", lw=1.2, zorder=1, label="φ(x) total = spline + w_b·SiLU" if leyenda else None)
    ax.grid(alpha=.35)
    ax.set_title(titulo, fontsize=10)
    ax.set_xlabel("x (entrada de la conexión)", fontsize=8)
    ax.tick_params(labelsize=7)


def main(clave="federico", archivo="G8_bsplines_estilo_clase.png", nombre="KAN [10,4,6] G=5 D²=100"):
    modelo = entrenar_o_cargar(clave)
    capa1, capa2 = capas_kan(modelo)
    datos1, datos2 = [d.numpy() for d in entradas_de_capas(modelo)]
    filas = []
    for k, (capa, datos, ni, no) in enumerate(((capa1, datos1, ENTRADAS, OCULTAS), (capa2, datos2, OCULTAS, SALIDAS))):
        with torch.no_grad():
            a = capa.activations(torch.as_tensor(datos)).numpy().std(0)
        for i in range(a.shape[0]):
            for o in range(a.shape[1]):
                filas.append({"capa": k + 1, "desde": ni[i], "hacia": no[o], "importancia_DE": a[i, o]})
    imp = pd.DataFrame(filas)
    imp["pct_de_su_capa"] = 100 * imp.importancia_DE / imp.groupby("capa").importancia_DE.transform("sum")
    top1 = imp[imp.capa == 1].sort_values("importancia_DE", ascending=False).head(8)
    top2 = imp[imp.capa == 2].sort_values("importancia_DE", ascending=False).head(4)
    fig, axes = plt.subplots(3, 4, figsize=(19, 13))
    for ax, (_, a) in zip(axes.flat[:8], top1.iterrows()):
        i, o = ENTRADAS.index(a.desde), OCULTAS.index(a.hacia)
        panel(ax, capa1, i, o, datos1[:, i], f"Capa 1: {a.desde} → {a.hacia}  ({a.pct_de_su_capa:.0f} % de la capa)", leyenda=ax is axes.flat[0])
    for ax, (_, a) in zip(axes.flat[8:], top2.iterrows()):
        i, o = OCULTAS.index(a.desde), SALIDAS.index(a.hacia)
        fuera = 100 * (np.abs(datos2[:, i]) > 3).mean()
        panel(ax, capa2, i, o, datos2[:, i], f"Capa 2: {a.desde} → {a.hacia}  ({a.pct_de_su_capa:.0f} %; {fuera:.0f} % de datos fuera)")
    axes.flat[0].legend(fontsize=7, loc="best")
    fig.suptitle(f"B-splines de la {nombre} al estilo de clase: las 8 conexiones más importantes de la capa 1 y las 4 de la capa 2\n"
                 "Cada color es un polinomio cúbico entre dos nodos (G=5 → 5 tramos en [−3, 3]); los puntos rojos son los 8 coeficientes de la conexión. "
                 "La línea gris discontinua es la función completa (spline + SiLU)", fontsize=12)
    fig.tight_layout()
    fig.savefig(SALIDA / archivo, dpi=115)
    plt.close(fig)
    print(SALIDA / archivo)


if __name__ == "__main__":
    torch.set_num_threads(4)
    if VERSION == "oficial":
        main("base", "G8_bsplines_estilo_clase.png", ETIQUETA)
    else:
        main()
        main("base", "G8b_bsplines_estilo_clase_base_sin_D2.png", "KAN base [10,4,6] G=5 sin D²")

"""Grafica las bases B-spline y las funciones aprendidas en cada conexion de las mejores KAN (semilla 2026)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import torch
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import kan_variantes_core as core
from kan_variantes_core import DATA, EXTRA_T, X_T, KANLayer, _indices
from experimentos_kan_jacobo import entrenar, red_con_cabeza

SALIDA = core.ROOT / "KAN" / "salidas" / "bsplines"
SALIDA.mkdir(parents=True, exist_ok=True)
SEMILLA = 2026
ENTRADAS = [f"z{k + 1}" for k in range(core.N_ENCODER)] + ["altitud", "ONI"]
OCULTAS = ["h1", "h2", "h3", "h4"]
SALIDAS = [f"FPCA {m + 1}" for m in range(core.M)]
MODELOS = {
    "federico": ("KAN [10,4,6] G=5 D2=100 (Federico)", "kan", 5, 100.0),
    "residual": ("KAN residual [10,4,6] G=3 D2=100", "residual", 3, 100.0),
    "base": ("KAN [10,4,6] G=5 D2=0 (base)", "kan", 5, 0.0),
}


def entrenar_o_cargar(clave):
    nombre, tipo, G, lam = MODELOS[clave]
    archivo = SALIDA / f"modelo_{clave}.pt"
    modelo = red_con_cabeza(tipo, (core.N_ENCODER + 2, 4, core.M), G) if tipo == "residual" else core.FunctionalKAN((core.N_ENCODER + 2, 4, core.M), G)
    if archivo.exists():
        modelo.load_state_dict(torch.load(archivo))
    else:
        if tipo == "kan":
            modelo, info = core.fit((core.N_ENCODER + 2, 4, core.M), G, SEMILLA, lam)
        else:
            torch.manual_seed(SEMILLA)
            modelo = red_con_cabeza(tipo, (core.N_ENCODER + 2, 4, core.M), G)
            modelo = entrenar_modelo(modelo, lam)
        torch.save(modelo.state_dict(), archivo)
    modelo.eval()
    return modelo


def entrenar_modelo(modelo, lam):
    opt = torch.optim.Adam(modelo.parameters(), lr=core.LR)
    tr, va = _indices("train"), _indices("valid")
    mejor, mejor_ep, estado = np.inf, 0, None
    for ep in range(core.MAX_EPOCHS):
        modelo.train()
        opt.zero_grad()
        dif = modelo(X_T[tr], EXTRA_T[tr]) - core.Y_T[tr]
        obj = torch.einsum("ni,ij,nj->n", dif, core.W_E_T, dif).mean() + core.LAMBDA_W * modelo.encoder_roughness(core.R_S_T, core.W_S_T)
        if lam:
            obj = obj + lam * modelo.head.second_difference_penalty()
        obj.backward()
        opt.step()
        modelo.eval()
        with torch.no_grad():
            actual = float(core.validation_loss(modelo, va).item())
        if actual < mejor - 1e-10:
            mejor, mejor_ep = actual, ep
            estado = {k: v.detach().clone() for k, v in modelo.state_dict().items()}
        elif ep - mejor_ep > core.PATIENCE:
            break
    modelo.load_state_dict(estado)
    return modelo


def capas_kan(modelo):
    cabeza = modelo.head
    return list(cabeza.kan.layers) if hasattr(cabeza, "kan") else list(cabeza.layers)


def entradas_de_capas(modelo):
    tr = _indices("train")
    with torch.no_grad():
        x = modelo.head_input(X_T[tr], EXTRA_T[tr])
        lista = []
        for capa in capas_kan(modelo):
            lista.append(x)
            x = capa(x)
    return lista


def funcion_arista(capa, i, o, x):
    with torch.no_grad():
        xx = torch.as_tensor(x, dtype=torch.float64).reshape(-1, 1).repeat(1, capa.coef.shape[0])
        b = capa.bases(xx)
        spline = (b[:, i, :] @ capa.coef[i, o, :]).numpy() * float(capa.w_s[i, o])
        silu = float(capa.w_b[i, o]) * torch.nn.functional.silu(xx[:, i]).numpy()
    return spline + silu, spline, silu


def figura_bases():
    fig, axes = plt.subplots(1, 2, figsize=(12, 3.8), sharey=True)
    x = np.linspace(-3, 3, 600)
    for ax, G in zip(axes, (5, 3)):
        capa = KANLayer(1, 1, G)
        with torch.no_grad():
            B = capa.bases(torch.as_tensor(x).reshape(-1, 1))[:, 0, :].numpy()
        for g in range(B.shape[1]):
            ax.plot(x, B[:, g], lw=1.8, label=f"B{g + 1}")
        nodos = capa.grid.numpy()
        for n in nodos[(nodos >= -3) & (nodos <= 3)]:
            ax.axvline(n, color="#bbbbbb", lw=.8, ls=":")
        ax.set_title(f"G = {G}: {G} intervalos en [−3, 3] → {G + 3} bases cúbicas por conexión")
        ax.set_xlabel("valor de la entrada x")
        ax.legend(fontsize=7, ncol=2)
    axes[0].set_ylabel("B_g(x)")
    fig.suptitle("Bases B-spline cúbicas de cada conexión KAN (líneas punteadas = nodos de la rejilla)")
    fig.tight_layout()
    fig.savefig(SALIDA / "B0_bases_bspline_G5_G3.png", dpi=120)
    plt.close(fig)


def figura_capa(modelo, clave, k, nombres_in, nombres_out, datos, color):
    capa = capas_kan(modelo)[k]
    n_in, n_out = capa.coef.shape[:2]
    fig, axes = plt.subplots(n_out, n_in, figsize=(1.65 * n_in + 1, 1.45 * n_out + 1), squeeze=False)
    filas = []
    for o in range(n_out):
        for i in range(n_in):
            d = datos[:, i].numpy()
            lo, hi = np.quantile(d, [.01, .99])
            x = np.linspace(min(lo, -3), max(hi, 3), 300)
            total, spline, silu = funcion_arista(capa, i, o, x)
            dentro = (x >= lo) & (x <= hi)
            ax = axes[o, i]
            ax.axvspan(lo, hi, color="#eeeeee", zorder=0)
            ax.plot(x, total, color=color, lw=1.6)
            ax.plot(x[~dentro], total[~dentro], color="white", lw=2.2, alpha=.6)
            ax.axhline(0, color="#999999", lw=.5)
            ax.set_xticks([]); ax.set_yticks([])
            if o == 0:
                ax.set_title(nombres_in[i], fontsize=8)
            if i == 0:
                ax.set_ylabel(nombres_out[o], fontsize=8)
            td, _, _ = funcion_arista(capa, i, o, d)
            r2 = np.corrcoef(d, td)[0, 1] ** 2 if np.std(td) > 0 else 1.0
            filas.append({"modelo": clave, "capa": k + 1, "entrada": nombres_in[i], "salida": nombres_out[o],
                          "DE_salida": float(np.std(td)), "no_linealidad_1_menos_R2": 1 - r2})
    titulo = f"{MODELOS[clave][0]}: capa {k + 1} ({n_in} → {n_out}), φ(x) de cada conexión"
    fig.suptitle(titulo + "\nfondo gris = rango p1–p99 de los datos de entrenamiento", fontsize=10)
    fig.tight_layout()
    fig.savefig(SALIDA / f"{clave}_capa{k + 1}.png", dpi=110)
    plt.close(fig)
    return filas


def figura_top(modelos, aristas):
    fig, axes = plt.subplots(len(modelos), 4, figsize=(14, 3.1 * len(modelos)), squeeze=False)
    for r, (clave, modelo) in enumerate(modelos.items()):
        datos = entradas_de_capas(modelo)
        top = aristas[aristas.modelo == clave].sort_values("DE_salida", ascending=False).head(4)
        nombres = {1: (ENTRADAS, OCULTAS), 2: (OCULTAS, SALIDAS)}
        for c, (_, a) in enumerate(top.iterrows()):
            k = int(a.capa) - 1
            capa = capas_kan(modelo)[k]
            i = nombres[k + 1][0].index(a.entrada)
            o = nombres[k + 1][1].index(a.salida)
            d = datos[k][:, i].numpy()
            lo, hi = np.quantile(d, [.01, .99])
            x = np.linspace(lo, hi, 300)
            total, spline, silu = funcion_arista(capa, i, o, x)
            ax = axes[r, c]
            ax.plot(x, total, color="k", lw=2, label="φ(x) total")
            ax.plot(x, spline, color="#3366aa", lw=1.2, ls="--", label="w_s·Σ c_g B_g(x)")
            ax.plot(x, silu, color="#cc6633", lw=1.2, ls=":", label="w_b·SiLU(x)")
            ax.hist(d, bins=40, weights=np.full(len(d), (total.max() - total.min() + 1e-9) * .35 / len(d) * 40), bottom=total.min(),
                    color="#dddddd", zorder=0)
            ax.set_title(f"capa {k + 1}: {a.entrada} → {a.salida}\nDE {a.DE_salida:.2f}, no linealidad {a.no_linealidad_1_menos_R2:.2f}", fontsize=8)
            if c == 0:
                ax.set_ylabel(MODELOS[clave][0], fontsize=8)
            if r == 0 and c == 0:
                ax.legend(fontsize=7)
    fig.suptitle("Las 4 conexiones más importantes de cada KAN: la función total y sus dos partes (spline y SiLU); en gris, la distribución de la entrada")
    fig.tight_layout()
    fig.savefig(SALIDA / "B1_conexiones_principales.png", dpi=115)
    plt.close(fig)


def main():
    torch.set_num_threads(6)
    figura_bases()
    modelos = {k: entrenar_o_cargar(k) for k in MODELOS}
    filas = []
    colores = {"federico": "#3366aa", "residual": "#cc6633", "base": "#777777"}
    for clave, modelo in modelos.items():
        datos = entradas_de_capas(modelo)
        filas += figura_capa(modelo, clave, 0, ENTRADAS, OCULTAS, datos[0], colores[clave])
        filas += figura_capa(modelo, clave, 1, OCULTAS, SALIDAS, datos[1], colores[clave])
    aristas = pd.DataFrame(filas)
    aristas.to_csv(SALIDA / "aristas_importancia_no_linealidad.csv", index=False)
    figura_top(modelos, aristas)
    rug = []
    for clave, modelo in modelos.items():
        for k, capa in enumerate(capas_kan(modelo)):
            c = (capa.w_s.unsqueeze(-1) * capa.coef).detach()
            rug.append({"modelo": clave, "capa": k + 1, "rugosidad_media_dif2": float(torch.diff(c, n=2, dim=-1).square().mean()),
                        "no_linealidad_media": aristas[(aristas.modelo == clave) & (aristas.capa == k + 1)].no_linealidad_1_menos_R2.mean(),
                        "pct_entradas_fuera_de_[-3,3]": 100 * float(((entradas_de_capas(modelo)[k].abs() > 3).double().mean()))})
    rug = pd.DataFrame(rug)
    rug.to_csv(SALIDA / "rugosidad_por_capa.csv", index=False)
    print(rug.round(4).to_string(index=False))
    print(aristas.sort_values("DE_salida", ascending=False).groupby("modelo").head(5).round(3).to_string(index=False))


if __name__ == "__main__":
    main()

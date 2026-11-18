"""Graficas de interpretabilidad de la KAN pedidas en la asesoria (semilla 2026).

Uso: python -X utf8 interpretabilidad_kan_jacobo.py [federico|oficial]
federico = KAN [10,4,6] G=5 D2=100 (carpeta interpretabilidad); oficial = KAN oficial [10,4,6] G=5 sin D2,
identica a entrenar_kan_base_oficial.py (carpeta interpretabilidad_oficial).

G1 prediccion con la historia de la estacion; G2 casos de lluvia alta, baja y variacional con sus B-splines y
nodos de suma; G3 B-splines con todos los casos coloreados por lluvia; G4 diagrama de la red; G5 frecuencia de
activacion KAN frente a MLP; G6 perdida y gradientes por capa; G7 desfase de los picos.
Solo entrenamiento y validacion (1981-2000); la prueba 2001-2003 no se usa.
"""

from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd
import torch
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from matplotlib.lines import Line2D

import kan_variantes_core as core
from kan_variantes_core import DATA, EXTRA_T, X_T, Y_T, KAN, FunctionalKAN, _indices
from experimentos_kan_jacobo import red_con_cabeza
from bsplines_kan_jacobo import SALIDA as DIR_MODELOS, capas_kan, entrenar_o_cargar

CONFIG = {"federico": ("interpretabilidad", "federico", "KAN [10,4,6] G=5 D²=100", 100.0),
          "oficial": ("interpretabilidad_oficial", "base", "KAN oficial [10,4,6] G=5 sin D²", 0.0)}
VERSION = next((a for a in sys.argv[1:] if a in CONFIG), "federico")
CARPETA, CLAVE_MODELO, ETIQUETA, LAM_KAN = CONFIG[VERSION]
SALIDA = core.ROOT / "KAN" / "salidas" / CARPETA
SALIDA.mkdir(parents=True, exist_ok=True)
SEMILLA = 2026
DIAS = np.arange(1, 366, dtype=float)
ENTRADAS = [f"z{k + 1}" for k in range(core.N_ENCODER)] + ["altitud", "ONI"]
OCULTAS = ["h1", "h2", "h3", "h4"]
SALIDAS = [f"FPCA{m + 1}" for m in range(core.M)]
MESES = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]
NOM_MESES = ["ene", "feb", "mar", "abr", "may", "jun", "jul", "ago", "sep", "oct", "nov", "dic"]
C_OBS, C_SAT, C_KAN = "#111111", "#d08c2b", "#2f5fa8"
C_NEUR = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a"]

pares = pd.read_csv(core.INPUT / "pares.csv", dtype={"codigo": str})
FILAS = np.concatenate([np.flatnonzero(pares["particion"].to_numpy() == p) for p in ("entrenamiento", "validación")])
P = pares.iloc[FILAS].reset_index(drop=True)
P["indice"] = np.arange(len(P))
Phi_e = core._basis_at_nodes(DIAS, DATA["tau_e"])
Phi_s = core._basis_at_nodes(DIAS, DATA["tau_s"])
coef_mm = pd.read_csv(core.INPUT / "coef_estacion_mm.csv").to_numpy(float)
coef_sat_mm = pd.read_csv(core.INPUT / "coef_satelite_mm.csv").to_numpy(float)
OBS_TODAS = np.maximum(coef_mm @ Phi_e.T, 0)
SAT_TODAS = np.maximum(coef_sat_mm @ Phi_s.T, 0)
OBS = OBS_TODAS[FILAS]
SAT = SAT_TODAS[FILAS]
MU_T = Phi_e @ DATA["mu"]
PSI_T = Phi_e @ DATA["basis_fpca"]
P["anual_obs"] = OBS.mean(1) * 365
clim = P[P.particion == "entrenamiento"].groupby("codigo")["anual_obs"].mean()
P["anom_pct"] = 100 * (P["anual_obs"] / P["codigo"].map(clim) - 1)
P["categoria"] = pd.qcut(P["anom_pct"], 3, labels=["baja", "normal", "alta"])


def ejes_meses(ax, etiquetas=True):
    ax.set_xticks(MESES, NOM_MESES if etiquetas else [""] * 12, fontsize=8)
    ax.set_xlim(1, 365)


def descomponer(modelo, idx):
    capa1, capa2 = capas_kan(modelo)
    with torch.no_grad():
        hin = modelo.head_input(X_T[idx], EXTRA_T[idx])
        a1 = capa1.activations(hin)
        h = a1.sum(1)
        a2 = capa2.activations(h)
        sc = a2.sum(1)
        assert torch.allclose(sc, modelo.head(hin), atol=1e-9)
    raiz = MU_T + sc.numpy() @ PSI_T.T
    return {"hin": hin.numpy(), "a1": a1.numpy(), "h": h.numpy(), "a2": a2.numpy(), "scores": sc.numpy(),
            "raiz": raiz, "mm": np.maximum(raiz, 0) ** 2,
            "aporte_neurona": np.einsum("nom,tm->not", a2.numpy(), PSI_T)}


def phi(capa, i, o, x):
    with torch.no_grad():
        xx = torch.as_tensor(np.asarray(x, dtype=float)).reshape(-1, 1).repeat(1, capa.coef.shape[0])
        return capa.activations(xx)[:, i, o].numpy()


def l2(pred, obs):
    return np.sqrt(((pred - obs) ** 2).sum(-1) / (obs ** 2).sum(-1))


def elegir_casos(D):
    v = P[P.particion == "validación"].copy()
    v["L2"] = l2(D["mm"][v.indice], OBS[v.indice])
    buenos = v[v.L2 <= v.L2.quantile(.6)]
    alta = buenos.sort_values("anom_pct").iloc[-1]
    baja = buenos.sort_values("anom_pct").iloc[0]
    puntajes = []
    for _, r in buenos.iterrows():
        hist = P[(P.codigo == r.codigo) & (P.particion == "entrenamiento")].indice
        a = OBS[r.indice] - OBS[hist].mean(0)
        k = np.convolve(a, np.ones(45) / 45, mode="same")
        puntajes.append(min(k.max(), -k.min()) / OBS[hist].mean())
    variacional = buenos.iloc[int(np.argmax(puntajes))]
    return {"Lluvia alta": alta, "Lluvia baja": baja, "Variacional (sube y baja frente a su clima)": variacional}


def panel_historia(ax, r, D, titulo=True):
    hist = P[(P.codigo == r.codigo) & (P.anio < r.anio)]
    for _, h in hist.iterrows():
        ax.plot(DIAS, OBS[h.indice], color="#bbbbbb", lw=.6, alpha=.6)
    ent = hist[hist.particion == "entrenamiento"].indice
    ax.plot(DIAS, OBS[ent].mean(0), color="#666666", lw=1.6, ls="--", label=f"Clima de la estación 1981–1995 (media de {len(ent)} años)")
    ax.plot(DIAS, SAT[r.indice], color=C_SAT, lw=1.4, label="Satélite (entrada)")
    ax.plot(DIAS, OBS[r.indice], color=C_OBS, lw=2.2, label=f"Observada {r.anio}")
    ax.plot(DIAS, D["mm"][r.indice], color=C_KAN, lw=2.2, label=f"KAN {r.anio} (L² {l2(D['mm'][r.indice], OBS[r.indice]):.2f})")
    ax.set_ylabel("mm/día")
    ejes_meses(ax)
    if titulo:
        ax.set_title(f"Estación {r.codigo} ({r.vertiente}, {r.altitud:.0f} m) — {r.anio}, ONI {r.ONI_anual:+.2f} ({r.fase}); "
                     f"{r.anual_obs:.0f} mm/año, {r.anom_pct:+.0f} % frente a su clima", fontsize=9)
    ax.legend(fontsize=7, loc="upper left")


def g1_historia(D, casos):
    fig, axes = plt.subplots(len(casos), 1, figsize=(13, 3.6 * len(casos)))
    for ax, (nombre, r) in zip(axes, casos.items()):
        panel_historia(ax, r, D)
        ax.text(.995, .97, nombre, transform=ax.transAxes, ha="right", va="top", fontsize=10, fontweight="bold")
    fig.suptitle("G1. Predicción del año T (validación) sobre todo el pasado de la estación: años anteriores en gris, clima de entrenamiento, satélite, observada y KAN",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(SALIDA / "G1_prediccion_con_historia.png", dpi=120)
    plt.close(fig)


def g2_caso(modelo, D, nombre, r, archivo):
    capa1, capa2 = capas_kan(modelo)
    i = int(r.indice)
    fig = plt.figure(figsize=(20, 15))
    gs = fig.add_gridspec(4, 12, height_ratios=[1.25, 1.0, 1.0, 1.05], hspace=.55, wspace=.45)
    ax = fig.add_subplot(gs[0, :6])
    panel_historia(ax, r, D, titulo=False)
    ax.set_title("Predicción del año sobre la historia de la estación", fontsize=10)
    ax = fig.add_subplot(gs[0, 6:])
    ax.plot(DIAS, MU_T, color="#888888", lw=1.5, ls="--", label="μ(t): curva media (raíz)")
    for o in range(4):
        ax.plot(DIAS, D["aporte_neurona"][i, o], color=C_NEUR[o], lw=1.8, label=f"aporte de {OCULTAS[o]}")
    ax.plot(DIAS, D["raiz"][i], color=C_KAN, lw=2.4, label="predicción KAN = μ + Σ aportes")
    ax.plot(DIAS, np.sqrt(OBS[i]), color=C_OBS, lw=1.4, alpha=.7, label="observada (raíz)")
    ax.axhline(0, color="#bbbbbb", lw=.6)
    ax.set_ylabel("raíz(mm/día)")
    ejes_meses(ax)
    ax.set_title("Nodo de salida: en qué época del año cada neurona empuja hacia más (+) o menos (−) lluvia", fontsize=10)
    ax.legend(fontsize=7, ncol=2)
    sub = gs[1:3, :].subgridspec(4, 10, hspace=.35, wspace=.25)
    for o in range(4):
        for k in range(10):
            ax = fig.add_subplot(sub[o, k])
            d = D["hin"][P.particion == "entrenamiento"][:, k]
            xv, yv = D["hin"][i, k], D["a1"][i, k, o]
            lo, hi = np.quantile(d, [.01, .99])
            x = np.linspace(min(lo, xv), max(hi, xv), 200)
            ax.axvspan(lo, hi, color="#f0f0f0", zorder=0)
            ax.plot(x, phi(capa1, k, o, x), color="#9db4d8", lw=1.6)
            ax.plot([xv], [yv], "o", color="#c0392b" if yv >= 0 else "#2471a3", ms=6)
            ax.axhline(0, color="#cccccc", lw=.6)
            ax.set_xticks([]); ax.set_yticks([])
            ax.text(.03, .95, f"{yv:+.2f}", transform=ax.transAxes, fontsize=7, va="top",
                    color="#c0392b" if yv >= 0 else "#2471a3")
            if o == 0:
                ax.set_title(ENTRADAS[k], fontsize=9)
            if k == 0:
                ax.set_ylabel(f"→ {OCULTAS[o]}", fontsize=9, color=C_NEUR[o], fontweight="bold")
    sub2 = gs[3, :].subgridspec(1, 5, wspace=.35, width_ratios=[1, 1, 1, 1, 2.2])
    for o in range(4):
        ax = fig.add_subplot(sub2[0, o])
        v = D["a1"][i, :, o]
        ax.barh(range(10), v, color=["#c0392b" if x >= 0 else "#2471a3" for x in v])
        ax.set_yticks(range(10), ENTRADAS, fontsize=7)
        ax.invert_yaxis()
        ax.axvline(0, color="#999999", lw=.6)
        ax.set_title(f"Nodo de suma {OCULTAS[o]} = {D['h'][i, o]:+.2f}", fontsize=9, color=C_NEUR[o])
    ax = fig.add_subplot(sub2[0, 4])
    M = D["a2"][i]
    lim = np.abs(M).max()
    im = ax.imshow(M, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    for o in range(4):
        for m in range(core.M):
            ax.text(m, o, f"{M[o, m]:+.1f}", ha="center", va="center", fontsize=7)
    ax.set_xticks(range(core.M), [f"{s}\n={D['scores'][i, m]:+.1f}" for m, s in enumerate(SALIDAS)], fontsize=7)
    ax.set_yticks(range(4), OCULTAS, fontsize=8)
    ax.set_title("Capa 2: φ(h) de cada neurona hacia cada puntaje FPCA (suma por columna = puntaje)", fontsize=9)
    fig.colorbar(im, ax=ax, fraction=.04)
    fig.suptitle(f"G2. {nombre}: estación {r.codigo} ({r.vertiente}, {r.altitud:.0f} m), año {r.anio} ({r.fase}, ONI {r.ONI_anual:+.2f}), "
                 f"{r.anual_obs:.0f} mm/año ({r.anom_pct:+.0f} % frente a su clima)\n"
                 "Filas centrales: B-spline φ(x) de cada conexión de la capa 1 y, en rojo/azul, dónde cae este caso (valor = aporte al nodo). "
                 "Abajo: suma en cada nodo y capa 2.", fontsize=11)
    fig.savefig(SALIDA / archivo, dpi=105, bbox_inches="tight")
    plt.close(fig)


def g3_generalizado(modelo, D):
    capa1, capa2 = capas_kan(modelo)
    norma = colors.TwoSlopeNorm(vmin=-40, vcenter=0, vmax=40)
    cmap = plt.get_cmap("BrBG")
    c = cmap(norma(np.clip(P["anom_pct"].to_numpy(), -40, 40)))
    orden = np.random.default_rng(0).permutation(len(P))
    fig, axes = plt.subplots(4, 10, figsize=(20, 8.5))
    for o in range(4):
        for k in range(10):
            ax = axes[o, k]
            x = np.linspace(*np.quantile(D["hin"][:, k], [.005, .995]), 200)
            ax.plot(x, phi(capa1, k, o, x), color="#bbbbbb", lw=4, zorder=1)
            ax.scatter(D["hin"][orden, k], D["a1"][orden, k, o], c=c[orden], s=7, alpha=.85, linewidths=0, zorder=2)
            ax.set_xticks([]); ax.set_yticks([])
            if o == 0:
                ax.set_title(ENTRADAS[k], fontsize=9)
            if k == 0:
                ax.set_ylabel(f"→ {OCULTAS[o]}", fontsize=9)
    sm = cm.ScalarMappable(norm=norma, cmap=cmap)
    fig.colorbar(sm, ax=axes, fraction=.015, pad=.01, label="lluvia del año frente al clima de la estación (%)")
    fig.suptitle("G3a. B-splines de la capa 1 con los 2 929 casos de 1981–2000: cada punto es un estación-año, "
                 "café = año seco, verde = año lluvioso. Muestra en qué tramo de cada spline caen los años secos y los lluviosos", fontsize=11)
    fig.savefig(SALIDA / "G3a_bsplines_capa1_generalizado.png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    fig, axes = plt.subplots(core.M, 4, figsize=(11, 13))
    for m in range(core.M):
        for o in range(4):
            ax = axes[m, o]
            x = np.linspace(*np.quantile(D["h"][:, o], [.005, .995]), 200)
            ax.plot(x, phi(capa2, o, m, x), color="#bbbbbb", lw=4, zorder=1)
            ax.scatter(D["h"][orden, o], D["a2"][orden, o, m], c=c[orden], s=7, alpha=.85, linewidths=0, zorder=2)
            ax.axvspan(-3, 3, color="#f2f2f2", zorder=0)
            ax.set_xticks([]); ax.set_yticks([])
            if m == 0:
                ax.set_title(f"desde {OCULTAS[o]}", fontsize=9)
            if o == 0:
                ax.set_ylabel(f"→ {SALIDAS[m]}", fontsize=9)
    fig.colorbar(sm, ax=axes, fraction=.03, pad=.02, label="lluvia del año frente al clima (%)")
    fig.suptitle("G3b. B-splines de la capa 2 (neuronas → puntajes FPCA) con todos los casos.\n"
                 "Fondo gris = rejilla [−3, 3]: fuera de ella la spline vale 0 y solo actúa SiLU", fontsize=11)
    fig.savefig(SALIDA / "G3b_bsplines_capa2_generalizado.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


def g4_diagrama(modelo, D):
    capa1, capa2 = capas_kan(modelo)
    ent = (P.particion == "entrenamiento").to_numpy()
    imp1 = D["a1"][ent].std(0)
    imp2 = D["a2"][ent].std(0)
    fig = plt.figure(figsize=(22, 13))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    y_in, y_h, y_out = .06, .52, .94
    x_in = np.linspace(.05, .95, 10)
    x_h = np.linspace(.2, .8, 4)
    x_out = np.linspace(.15, .85, core.M)
    caja = .016
    def dibujar_capa(xs, ys, xt, yt, imp, capa, datos, y_cajas, n_in, n_out, nombres_in):
        rel = imp / imp.max()
        orden = [(i, o) for o in range(n_out) for i in range(n_in)]
        xc = np.linspace(.02, .98, len(orden))
        for k, (i, o) in enumerate(orden):
            a = .12 + .88 * rel[i, o]
            ax.plot([xs[i], xc[k]], [ys + .02, y_cajas - caja * 1.6], color=C_NEUR[o % 4] if n_out == 4 else "#555555", lw=.4 + 3 * rel[i, o], alpha=a)
            ax.plot([xc[k], xt[o]], [y_cajas + caja * 1.6, yt - .02], color=C_NEUR[o % 4] if n_out == 4 else C_NEUR[i % 4], lw=.4 + 3 * rel[i, o], alpha=a)
            sub = fig.add_axes([xc[k] - caja * .7, y_cajas - caja * 1.5, caja * 1.4, caja * 3])
            x = np.linspace(*np.quantile(datos[:, i], [.01, .99]), 100)
            sub.plot(x, phi(capa, i, o, x), color="k", lw=1, alpha=a)
            sub.set_xticks([]); sub.set_yticks([])
            for s in sub.spines.values():
                s.set_alpha(a)
                s.set_linewidth(.5)
    dibujar_capa(x_in, y_in, x_h, y_h, imp1, capa1, D["hin"][ent], (y_in + y_h) / 2 + .02, 10, 4, ENTRADAS)
    dibujar_capa(x_h, y_h, x_out, y_out, imp2, capa2, D["h"][ent], (y_h + y_out) / 2 + .02, 4, core.M, OCULTAS)
    for xs, y, nombres, col in ((x_in, y_in, ENTRADAS, ["#444444"] * 10), (x_h, y_h, OCULTAS, C_NEUR), (x_out, y_out, SALIDAS, ["#444444"] * core.M)):
        for x, n, cc in zip(xs, nombres, col):
            ax.add_patch(plt.Circle((x, y), .014, color=cc, zorder=5))
            ax.text(x, y - .035 if y < .5 else y + .03, n, ha="center", fontsize=11, fontweight="bold")
    ax.text(.5, .985, f"G4. Diagrama de la {ETIQUETA}: cada caja es la B-spline φ(x) de una conexión; "
            "grosor y opacidad = importancia (DE de su aporte en entrenamiento)", ha="center", fontsize=13)
    ax.text(.01, y_h, "nodos de suma\n(neuronas ocultas)", fontsize=10, va="center")
    fig.savefig(SALIDA / "G4_diagrama_red_kan.png", dpi=110)
    plt.close(fig)
    t = []
    for i in range(10):
        for o in range(4):
            t.append({"capa": 1, "desde": ENTRADAS[i], "hacia": OCULTAS[o], "importancia_DE": imp1[i, o]})
    for o in range(4):
        for m in range(core.M):
            t.append({"capa": 2, "desde": OCULTAS[o], "hacia": SALIDAS[m], "importancia_DE": imp2[o, m]})
    t = pd.DataFrame(t)
    t["pct_de_su_capa"] = 100 * t.importancia_DE / t.groupby("capa").importancia_DE.transform("sum")
    t.sort_values(["capa", "importancia_DE"], ascending=[True, False]).to_csv(SALIDA / "G4_importancia_conexiones.csv", index=False)
    return t


def entrenar_registrando(tipo, ancho=37):
    widths = (core.N_ENCODER + 2, 4, core.M)
    torch.manual_seed(SEMILLA)
    if tipo == "kan":
        cab = KAN(widths, 5)
        torch.manual_seed(SEMILLA)
        modelo = FunctionalKAN(widths, 5)
        modelo.head.load_state_dict(cab.state_dict())
        lam = LAM_KAN
    else:
        modelo = red_con_cabeza("mlp", (core.N_ENCODER + 2, ancho, core.M), 5)
        lam = 0.0
    grupos = {"codificador": [modelo.encoder]}
    if tipo == "kan":
        for k, capa in enumerate(modelo.head.layers):
            grupos[f"capa {k + 1}: splines (coef)"] = [capa.coef]
            grupos[f"capa {k + 1}: w_b y w_s"] = [capa.w_b, capa.w_s]
    else:
        red = modelo.head.red
        grupos["capa oculta"] = list(red[0].parameters())
        grupos["capa de salida"] = list(red[2].parameters())
    opt = torch.optim.Adam(modelo.parameters(), lr=core.LR)
    tr, va = _indices("train"), _indices("valid")
    reg, mejor, mejor_ep, estado = [], np.inf, 0, None
    for ep in range(core.MAX_EPOCHS):
        modelo.train()
        opt.zero_grad()
        dif = modelo(X_T[tr], EXTRA_T[tr]) - Y_T[tr]
        l_tr = torch.einsum("ni,ij,nj->n", dif, core.W_E_T, dif).mean()
        obj = l_tr + core.LAMBDA_W * modelo.encoder_roughness(core.R_S_T, core.W_S_T)
        if lam:
            obj = obj + lam * modelo.head.second_difference_penalty()
        obj.backward()
        fila = {"epoca": ep + 1, "L2_entrenamiento": float(l_tr.detach())}
        for g, ps in grupos.items():
            fila[g] = float(torch.sqrt(sum((p.grad ** 2).sum() for p in ps)))
        opt.step()
        modelo.eval()
        with torch.no_grad():
            fila["L2_validacion"] = float(core.validation_loss(modelo, va))
        reg.append(fila)
        if fila["L2_validacion"] < mejor - 1e-10:
            mejor, mejor_ep = fila["L2_validacion"], ep
            estado = {k: v.detach().clone() for k, v in modelo.state_dict().items()}
        elif ep - mejor_ep > core.PATIENCE:
            break
    modelo.load_state_dict(estado)
    modelo.eval()
    return modelo, pd.DataFrame(reg), mejor_ep + 1


def g5_frecuencia(modelo_kan, D, mlp):
    cats = ["baja", "normal", "alta"]
    col = {"baja": "#a6611a", "normal": "#999999", "alta": "#018571"}
    with torch.no_grad():
        hin_m = mlp.head_input(X_T, EXTRA_T)
        act_m = mlp.head.red[1](mlp.head.red[0](hin_m)).numpy()
    fig = plt.figure(figsize=(18, 10))
    gs = fig.add_gridspec(2, 5, height_ratios=[1, 1.1], hspace=.6, wspace=.35)
    for o in range(4):
        ax = fig.add_subplot(gs[0, o])
        for cname in cats:
            m = (P.categoria == cname).to_numpy()
            ax.hist(D["h"][m, o], bins=40, histtype="step", lw=1.6, color=col[cname], density=True, label=f"años de lluvia {cname}")
        ax.axvspan(-3, 3, color="#f2f2f2", zorder=0)
        fuera = 100 * (np.abs(D["h"][:, o]) > 3).mean()
        ax.set_title(f"KAN: activación de {OCULTAS[o]}\n{fuera:.0f} % fuera de la rejilla\nde la capa 2 (fondo gris = dentro)", fontsize=9, color=C_NEUR[o])
        if o == 0:
            ax.legend(fontsize=7)
    ax = fig.add_subplot(gs[0, 4])
    dec = pd.qcut(P["anom_pct"], 10, labels=False).to_numpy()
    H = np.array([[D["h"][dec == d, o].mean() for d in range(10)] for o in range(4)])
    Hz = (H - H.mean(1, keepdims=True)) / H.std(1, keepdims=True)
    im = ax.imshow(Hz, cmap="RdBu_r", aspect="auto", vmin=-2, vmax=2)
    ax.set_yticks(range(4), OCULTAS)
    ax.set_xticks(range(10), [f"D{d + 1}" for d in range(10)], fontsize=7)
    ax.set_xlabel("decil de lluvia del año (D1 = más seco frente a su clima)")
    ax.set_title("KAN: activación media por decil\n(estandarizada por neurona)", fontsize=9)
    fig.colorbar(im, ax=ax, fraction=.05)
    ax = fig.add_subplot(gs[1, :3])
    Am = np.array([[act_m[dec == d, j].mean() for d in range(10)] for j in range(act_m.shape[1])])
    sd = act_m.std(0)
    orden = np.argsort(-np.corrcoef(np.vstack([Am, np.arange(10)]))[-1, :-1])
    Az = (Am - Am.mean(1, keepdims=True)) / (Am.std(1, keepdims=True) + 1e-12)
    im = ax.imshow(Az[orden], cmap="RdBu_r", aspect="auto", vmin=-2, vmax=2)
    ax.set_yticks(range(len(orden)), [f"n{j + 1}" for j in orden], fontsize=5)
    ax.set_xticks(range(10), [f"D{d + 1}" for d in range(10)], fontsize=8)
    ax.set_xlabel("decil de lluvia del año")
    ax.set_title(f"MLP de {act_m.shape[1]} neuronas (contraste): activación media por decil de lluvia\n(estandarizada; ordenadas de 'se activa en años secos' a 'en años lluviosos')", fontsize=9)
    fig.colorbar(im, ax=ax, fraction=.02)
    ax = fig.add_subplot(gs[1, 3:])
    muertas = (sd < .02 * sd.max()).sum()
    ax.bar(range(len(sd)), np.sort(sd)[::-1], color="#888888")
    ax.set_title(f"MLP: DE de la activación de cada neurona\nneuronas casi muertas: {muertas} de {len(sd)}", fontsize=9)
    ax.set_xlabel("neuronas ordenadas")
    fig.suptitle("G5. Diagrama de frecuencia de activación según la lluvia del año (1981–2000, terciles y deciles de anomalía frente al clima de cada estación)",
                 fontsize=12)
    fig.savefig(SALIDA / "G5_frecuencia_activacion.png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    return {"mlp_neuronas_casi_muertas": int(muertas), "kan_pct_fuera_rejilla_capa2": [float(100 * (np.abs(D["h"][:, o]) > 3).mean()) for o in range(4)],
            "kan_corr_activacion_vs_anomalia": [float(np.corrcoef(D["h"][:, o], P["anom_pct"])[0, 1]) for o in range(4)]}


def g6_gradientes(reg_kan, ep_kan, reg_mlp, ep_mlp):
    fig, axes = plt.subplots(2, 2, figsize=(15, 9))
    for c, (reg, ep, nombre) in enumerate(((reg_kan, ep_kan, ETIQUETA), (reg_mlp, ep_mlp, "MLP de 37 neuronas"))):
        ax = axes[0, c]
        ax.plot(reg.epoca, reg.L2_entrenamiento, label="entrenamiento", color="#555555")
        ax.plot(reg.epoca, reg.L2_validacion, label="validación", color=C_KAN)
        ax.axvline(ep, color="#c0392b", ls=":", label=f"época elegida ({ep})")
        ax.set_yscale("log"); ax.set_title(f"{nombre}: pérdida L² por época", fontsize=10); ax.legend(fontsize=8)
        ax = axes[1, c]
        for g in [k for k in reg.columns if k not in ("epoca", "L2_entrenamiento", "L2_validacion")]:
            ax.plot(reg.epoca, reg[g], lw=1, label=g)
        ax.axvline(ep, color="#c0392b", ls=":")
        ax.set_yscale("log"); ax.set_xlabel("época"); ax.set_ylabel("norma del gradiente")
        ax.set_title(f"{nombre}: norma del gradiente por capa\n(explosión = crece sin control; desvanecimiento = cae a ~0)", fontsize=9)
        ax.legend(fontsize=7)
    fig.suptitle("G6. Entrenamiento: pérdida y gradientes por capa", fontsize=12)
    fig.tight_layout()
    fig.savefig(SALIDA / "G6_perdida_y_gradientes.png", dpi=110)
    plt.close(fig)


def desfase(pred, obs, max_lag=60):
    a = pred - pred.mean()
    b = obs - obs.mean()
    lags = np.arange(-max_lag, max_lag + 1)
    cc = [np.corrcoef(np.roll(a, -L), b)[0, 1] for L in lags]
    return lags[int(np.argmax(cc))]


def fase_armonica(y, k):
    w = 2 * np.pi * k * DIAS / 365
    A = np.column_stack([np.ones_like(DIAS), np.cos(w), np.sin(w)])
    c = np.linalg.lstsq(A, y, rcond=None)[0]
    return np.arctan2(c[2], c[1])


def g7_desfase(D):
    res = []
    for _, r in P.iterrows():
        i = r.indice
        f = {"codigo": r.codigo, "anio": r.anio, "particion": r.particion,
             "desfase_KAN_dias": desfase(D["mm"][i], OBS[i]), "desfase_satelite_dias": desfase(SAT[i], OBS[i])}
        for k, per in ((1, 365), (2, 182.5)):
            d = (fase_armonica(D["mm"][i], k) - fase_armonica(OBS[i], k) + np.pi) % (2 * np.pi) - np.pi
            f[f"desfase_armonico_{k}_dias"] = d / (2 * np.pi) * per
        res.append(f)
    res = pd.DataFrame(res)
    ent, val = res.particion == "entrenamiento", res.particion == "validación"
    lag = int(np.round(res.loc[ent, "desfase_KAN_dias"].median()))
    vi = P[P.particion == "validación"].indice.to_numpy()
    corr = np.roll(D["mm"][vi], -lag, axis=1)
    antes, despues = l2(D["mm"][vi], OBS[vi]), l2(corr, OBS[vi])
    cod = P.loc[vi, "codigo"].to_numpy()
    resumen = {"desfase_mediano_KAN_entrenamiento_dias": lag,
               "desfase_mediano_KAN_validacion_dias": float(res.loc[val, "desfase_KAN_dias"].median()),
               "desfase_mediano_satelite_validacion_dias": float(res.loc[val, "desfase_satelite_dias"].median()),
               "desfase_armonico_anual_mediano_validacion_dias": float(res.loc[val, "desfase_armonico_1_dias"].median()),
               "desfase_armonico_semestral_mediano_validacion_dias": float(res.loc[val, "desfase_armonico_2_dias"].median()),
               "L2_validacion_sin_correccion": float(pd.Series(antes).groupby(cod).mean().mean()),
               "L2_validacion_con_correccion": float(pd.Series(despues).groupby(cod).mean().mean())}
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))
    ax = axes[0]
    b = np.arange(-60.5, 61.5, 3)
    ax.hist(res.loc[val, "desfase_satelite_dias"], bins=b, alpha=.6, color=C_SAT, label="satélite frente a observada")
    ax.hist(res.loc[val, "desfase_KAN_dias"], bins=b, alpha=.6, color=C_KAN, label="KAN frente a observada")
    ax.axvline(0, color="k", lw=.8)
    ax.set_xlabel("desfase (días; + = la curva va atrasada)")
    ax.set_title(f"Desfase por correlación cruzada (validación)\nmediana KAN {resumen['desfase_mediano_KAN_validacion_dias']:+.0f} d, "
                 f"satélite {resumen['desfase_mediano_satelite_validacion_dias']:+.0f} d", fontsize=9)
    ax.legend(fontsize=8)
    ax = axes[1]
    ax.hist(res.loc[val, "desfase_armonico_1_dias"], bins=np.arange(-60, 61, 4), alpha=.6, color="#5e3c99", label="armónico anual (1 ciclo/año)")
    ax.hist(res.loc[val, "desfase_armonico_2_dias"], bins=np.arange(-60, 61, 4), alpha=.6, color="#e66101", label="armónico semestral (2 ciclos/año)")
    ax.axvline(0, color="k", lw=.8)
    ax.set_xlabel("desfase de fase, KAN − observada (días)")
    ax.set_title(f"Ángulo de desfase de los armónicos (validación)\nmedianas: anual {resumen['desfase_armonico_anual_mediano_validacion_dias']:+.1f} d, "
                 f"semestral {resumen['desfase_armonico_semestral_mediano_validacion_dias']:+.1f} d", fontsize=9)
    ax.legend(fontsize=8)
    ax = axes[2]
    ax.plot(DIAS, OBS[vi].mean(0), color=C_OBS, lw=2, label="observada (media de validación)")
    ax.plot(DIAS, D["mm"][vi].mean(0), color=C_KAN, lw=2, label="KAN")
    ax.plot(DIAS, SAT[vi].mean(0), color=C_SAT, lw=1.5, label="satélite")
    if lag:
        ax.plot(DIAS, corr.mean(0), color=C_KAN, lw=1.2, ls="--", label=f"KAN corregida {lag:+d} d")
    ejes_meses(ax)
    ax.set_ylabel("mm/día")
    ax.set_title(f"Curva media: corrección con el desfase de entrenamiento ({lag:+d} d)\nL² {resumen['L2_validacion_sin_correccion']:.3f} → "
                 f"{resumen['L2_validacion_con_correccion']:.3f}", fontsize=9)
    ax.legend(fontsize=7)
    fig.suptitle("G7. Análisis de desfase de la predicción KAN", fontsize=12)
    fig.tight_layout()
    fig.savefig(SALIDA / "G7_desfase.png", dpi=110)
    plt.close(fig)
    res.to_csv(SALIDA / "G7_desfase_por_curva.csv", index=False)
    return resumen


def main():
    torch.set_num_threads(6)
    modelo = entrenar_o_cargar(CLAVE_MODELO)
    idx = torch.arange(X_T.shape[0])
    D = descomponer(modelo, idx)
    casos = elegir_casos(D)
    g1_historia(D, casos)
    for k, (nombre, r) in enumerate(casos.items()):
        g2_caso(modelo, D, nombre, r, f"G2_{'abc'[k]}_{['lluvia_alta', 'lluvia_baja', 'variacional'][k]}.png")
    g3_generalizado(modelo, D)
    imp = g4_diagrama(modelo, D)
    archivo_mlp = DIR_MODELOS / "modelo_mlp37.pt"
    _, reg_kan, ep_kan = entrenar_registrando("kan")
    mlp, reg_mlp, ep_mlp = entrenar_registrando("mlp")
    torch.save(mlp.state_dict(), archivo_mlp)
    reg_kan.to_csv(SALIDA / "G6_registro_kan.csv", index=False)
    reg_mlp.to_csv(SALIDA / "G6_registro_mlp.csv", index=False)
    g6_gradientes(reg_kan, ep_kan, reg_mlp, ep_mlp)
    act = g5_frecuencia(modelo, D, mlp)
    des = g7_desfase(D)
    resumen = {"casos": {k: {"codigo": r.codigo, "anio": int(r.anio), "anual_mm": float(r.anual_obs), "anomalia_pct": float(r.anom_pct)} for k, r in casos.items()},
               "activacion": act, "desfase": des,
               "top_conexiones": imp.sort_values("importancia_DE", ascending=False).groupby("capa").head(5).to_dict("records"),
               "gradientes_kan_inicio_fin": {c: [float(reg_kan[c].iloc[:10].mean()), float(reg_kan[c].iloc[ep_kan - 10:ep_kan].mean())] for c in reg_kan.columns if c not in ("epoca", "L2_entrenamiento", "L2_validacion")},
               "epocas": {"kan": ep_kan, "mlp": ep_mlp}}
    (SALIDA / "resumen_interpretabilidad.json").write_text(json.dumps(resumen, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    print(json.dumps(resumen, ensure_ascii=False, indent=1, default=str))


if __name__ == "__main__":
    main()

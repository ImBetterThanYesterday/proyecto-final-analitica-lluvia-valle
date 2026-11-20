"""Diagnostico del MLP elegido (M25: 3 capas 14-14-14 + L2=200) frente al MLP37 sin regularizar.

Bloque 1996-2000 (particion original), semilla 2026. Se reentrena registrando por epoca la perdida
de entrenamiento y validacion y la norma del gradiente de cada capa; despues se miran los pesos
y la activacion de las neuronas ocultas en entrenamiento.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import torch
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import kan_variantes_core as core
from kan_variantes_core import EXTRA_T, LAMBDA_W, MAX_EPOCHS, PATIENCE, R_S_T, W_S_T, X_T
from cv_bloques_kan import ANIO, BLOQUES, perdida
from cv_mlp import RONDAS, cfg_completa, construir, penalizacion_l2

BASE = core.ROOT / "KAN" / "salidas"
SEMILLA = 2026
MODELOS = {"MLP37 sin regularizar": cfg_completa({}), "M25 3 capas 14-14-14 + L2=200": cfg_completa(RONDAS["r3"]["M25 MLP 3 capas 14-14-14 + L2=200"])}
COLOR = {"MLP37 sin regularizar": "#888888", "M25 3 capas 14-14-14 + L2=200": "#c0392b"}


def indices():
    a, b = BLOQUES["1996-2000"]
    fuera = (ANIO >= a) & (ANIO <= b)
    return torch.as_tensor(np.flatnonzero(~fuera), dtype=torch.long), torch.as_tensor(np.flatnonzero(fuera), dtype=torch.long)


def entrenar_registrando(cfg, idx_ent, idx_fuera):
    modelo = construir(cfg, SEMILLA, idx_ent)
    lineales = [m for m in modelo.head.modules() if isinstance(m, torch.nn.Linear)]
    opt = torch.optim.Adam(modelo.parameters(), lr=cfg["lr"])
    reg, mejor, mejor_ep, estado = [], np.inf, 0, None
    for ep in range(MAX_EPOCHS):
        modelo.train()
        opt.zero_grad()
        l_tr = perdida(modelo, idx_ent)
        obj = l_tr + LAMBDA_W * modelo.encoder_roughness(R_S_T, W_S_T)
        if cfg["l2"]:
            obj = obj + cfg["l2"] * penalizacion_l2(modelo)
        obj.backward()
        fila = {"epoca": ep + 1, "entrenamiento": float(l_tr.detach()), "codificador": float(modelo.encoder.grad.norm())}
        fila.update({f"capa {k + 1}": float(c.weight.grad.norm()) for k, c in enumerate(lineales)})
        opt.step()
        modelo.eval()
        with torch.no_grad():
            fila["validacion"] = float(perdida(modelo, idx_fuera))
        reg.append(fila)
        if fila["validacion"] < mejor - 1e-10:
            mejor, mejor_ep = fila["validacion"], ep
            estado = {k: v.detach().clone() for k, v in modelo.state_dict().items()}
        elif ep - mejor_ep > PATIENCE:
            break
    modelo.load_state_dict(estado)
    modelo.eval()
    return modelo, pd.DataFrame(reg), mejor_ep + 1


def activaciones(modelo, idx_ent):
    salida, x = [], modelo.head_input(X_T[idx_ent], EXTRA_T[idx_ent])
    with torch.no_grad():
        for capa in modelo.head.red:
            x = capa(x)
            if isinstance(capa, (torch.nn.SiLU, torch.nn.ReLU, torch.nn.Tanh)):
                salida.append(x.numpy())
    return salida


def main():
    torch.set_num_threads(4)
    idx_ent, idx_fuera = indices()
    res = {}
    for nombre, cfg in MODELOS.items():
        modelo, reg, ep = entrenar_registrando(cfg, idx_ent, idx_fuera)
        pesos = np.concatenate([m.weight.detach().numpy().ravel() for m in modelo.head.modules() if isinstance(m, torch.nn.Linear)])
        res[nombre] = (modelo, reg, ep, pesos, activaciones(modelo, idx_ent))
        print(nombre, "época", ep, "validación", round(reg.validacion.min(), 3))

    fig, axes = plt.subplots(2, 3, figsize=(19, 10))
    ax = axes[0, 0]
    for nombre, (_, reg, ep, _, _) in res.items():
        ax.plot(reg.epoca, reg.entrenamiento, color=COLOR[nombre], lw=1, alpha=.5)
        ax.plot(reg.epoca, reg.validacion, color=COLOR[nombre], lw=2, label=f"{nombre} (mejor época {ep})")
        ax.axvline(ep, color=COLOR[nombre], ls=":")
    ax.set_yscale("log")
    ax.set_ylim(top=300)
    ax.set_xlabel("época")
    ax.set_ylabel("pérdida (delgada: entrenamiento, gruesa: validación)")
    ax.set_title("Pérdida: el MLP37 sobreajusta en pocas épocas;\nM25 aprende más despacio y llega más bajo en validación", fontsize=10)
    ax.legend(fontsize=8)
    for ax, nombre in zip(axes[0, 1:], res):
        _, reg, ep, _, _ = res[nombre]
        for c in [c for c in reg.columns if c.startswith("capa") or c == "codificador"]:
            ax.plot(reg.epoca, reg[c], lw=.9, label=c)
        ax.axvline(ep, color="#c0392b", ls=":")
        ax.set_yscale("log")
        ax.set_xlabel("época")
        ax.set_ylabel("norma del gradiente")
        ax.set_title(f"Gradientes por capa: {nombre}", fontsize=10)
        ax.legend(fontsize=7)
    ax = axes[1, 0]
    for nombre, (_, _, _, pesos, _) in res.items():
        ax.hist(pesos, bins=80, range=(-2.5, 2.5), color=COLOR[nombre], alpha=.55, density=True,
                label=f"{nombre}: media |w| = {np.abs(pesos).mean():.2f}")
    ax.set_xlabel("valor del peso")
    ax.set_title("Pesos de la cabeza: la penalización L2 los encoge hacia 0", fontsize=10)
    ax.legend(fontsize=8)
    ax = axes[1, 1]
    acts = res["M25 3 capas 14-14-14 + L2=200"][4]
    ax.boxplot([a.std(0) for a in acts] + [res["MLP37 sin regularizar"][4][0].std(0)],
               tick_labels=["M25 capa 1", "M25 capa 2", "M25 capa 3", "MLP37 capa 1"])
    ax.set_ylabel("DE de la activación de cada neurona")
    ax.set_title("Variabilidad de cada neurona oculta (entrenamiento):\nuna neurona con DE ≈ 0 estaría muerta", fontsize=10)
    ax = axes[1, 2]
    filas = []
    for nombre, (_, _, _, _, acts) in res.items():
        for k, a in enumerate(acts):
            filas.append({"modelo": nombre, "capa": k + 1, "neuronas": a.shape[1], "muertas_DE<0,01": int((a.std(0) < .01).sum()),
                          "pct_activas_x>0": 100 * (a > 0).mean()})
    tabla = pd.DataFrame(filas)
    tabla.to_csv(BASE / "mlp_diagnostico_neuronas.csv", index=False)
    ax.axis("off")
    t = ax.table(cellText=tabla.assign(modelo=tabla.modelo.str.split().str[0]).round(1).values, colLabels=["modelo", "capa", "neuronas", "muertas (DE<0,01)", "% activación > 0"], loc="center")
    t.auto_set_font_size(False)
    t.set_fontsize(11)
    t.scale(1, 2.2)
    ax.set_title("Neuronas muertas y porcentaje de activación positiva", fontsize=10)
    fig.suptitle("Diagnóstico del MLP elegido (M25) frente al MLP37 sin regularizar: bloque 1996–2000, semilla 2026", fontsize=12, y=.995)
    fig.tight_layout(rect=(0, 0, 1, .98))
    fig.savefig(BASE / "fig8_mlp_diagnostico.png", dpi=115)
    print(tabla.to_string(index=False))


if __name__ == "__main__":
    main()

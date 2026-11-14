"""
LLUVIA EN EL VALLE DEL CAUCA 1981-2003 — RED NEURONAL FUNCIONAL Y ABLACIÓN DE LA CABEZA
Proyecto final de Analítica de Datos (UAO, 2026-2S)

Reconstruye la curva anual de una estación a partir de la curva anual del satélite en su píxel.

    codificador  z_k = ∫ w_k(s) X(s) ds = a_kᵀ W_s c,   k = 1..N_NEURONAS        (rúbrica 4.2)
                 + altitud de la estación (estandarizada con entrenamiento) + ONI del año
    cabeza       lineal | KAN | MLP  ->  M puntajes del FPCA                        (rúbrica 4.3)
    decodificador ŷ(t) = μ(t) + Σ_m s_m ψ_m(t), μ y ψ del FPCA de entrenamiento, congelados
    pérdida      (ĉ − c)ᵀ W_e (ĉ − c) + λ_w · Σ_k (a_kᵀ R a_k) / (a_kᵀ W_s a_k)

Todo en escala √P (como la regresión de R); la curva final en mm es max(ŷ, 0)².

Entrada única: salidas_R_funcional/para_python/ (sale de lluvia_valle_funcional.R).
Salidas:       salidas_python/tablas y salidas_python/graficos, con prefijos P1_ a P6_ por sección.
Uso, desde la raíz del repositorio:  python lluvia_valle_red.py
Librerías del apartado (enunciado): PyTorch. La KAN se programa aquí con la forma del enunciado,
φ(x) = w_b·SiLU(x) + w_s·Σ_i c_i B_i(x), para poder graficar y verificar cada arista.

Secciones:
  1. Insumos de R y verificación de que Python reproduce lo mismo
  2. La red funcional (codificador, decodificador y las tres cabezas)
  3. Cabeza lineal = regresión funcional: demostración y verificación numérica
  4. KAN verificada con una función sintética y MLP con los mismos parámetros
  5. Ablación: tres cabezas, mismo presupuesto, tres semillas
  6. Pruebas y cierre
"""

from __future__ import annotations

import json
import platform
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from scipy import stats

torch.set_default_dtype(torch.float64)
torch.set_num_threads(4)

RAIZ = Path(__file__).resolve().parent
ENTRADA = RAIZ / "salidas_R_funcional" / "para_python"
SALIDA = RAIZ / "salidas_python"
TAB, GRA = SALIDA / "tablas", SALIDA / "graficos"
for d in (TAB, GRA):
    d.mkdir(parents=True, exist_ok=True)

DIAS = np.arange(1, 366, dtype=float)
SEMILLAS = (2026, 2027, 2028)          # rúbrica 4.3: tres semillas por configuración
N_NEURONAS = 8                          # enunciado: entre 6 y 8 (sección 3 justifica 8)
# Optimizador: con lr = 3e-3 y paciencia 300 la cabeza lineal se detenía antes de converger (pérdida de
# validación 142, peor que mínimos cuadrados sin penalizar, 137). Con lr = 1e-2 y paciencia 1 000 llega al
# mínimo de entrenamiento de mínimos cuadrados (110), y con λ_w = 1e8 a la validación de la regresión de R.
EPOCAS_MAX, PACIENCIA, LR = 5000, 1000, 1e-2
# λ_w: la rugosidad relativa ∫w''²/∫w² de una w_k cualquiera es del orden de 1e-6 (día^-4), así que la
# rejilla va de 1e6 (penalización leve) a 1e9 (w_k casi lineales). Misma rejilla para las tres cabezas.
LAMBDAS_W = (0.0, 1e6, 1e7, 1e8, 1e9)
N_BOOT = 2000
PRUEBAS = {}


def titulo(texto):
    print("\n" + "=" * 78 + "\n" + texto + "\n" + "=" * 78, flush=True)


def tabla(df, nombre, indice=False):
    df.to_csv(TAB / f"{nombre}.csv", index=indice)
    return df


def prueba(nombre, ok, detalle=""):
    PRUEBAS[nombre] = (bool(ok), detalle)
    print(("  OK  " if ok else "  FALLA ") + nombre + (f"  ({detalle})" if detalle else ""), flush=True)


# =============================================================================================
# Base B-spline (Cox–de Boor), misma definición que en R
# =============================================================================================
def base_en_nodos(t, tau, orden=4):
    t = np.asarray(t, float)
    B = np.zeros((len(t), len(tau) - 1))
    for j in range(len(tau) - 1):
        B[:, j] = (t >= tau[j]) & (t < tau[j + 1])
    B[t == tau[-1], np.max(np.flatnonzero(tau[:-1] < tau[1:]))] = 1
    for m in range(2, orden + 1):
        Bn = np.zeros((len(t), len(tau) - m))
        for j in range(len(tau) - m):
            d1, d2 = tau[j + m - 1] - tau[j], tau[j + m] - tau[j + 1]
            Bn[:, j] = ((t - tau[j]) / d1 * B[:, j] if d1 > 0 else 0) + ((tau[j + m] - t) / d2 * B[:, j + 1] if d2 > 0 else 0)
        B = Bn
    return B


def cuadratura(tau, n_gl=4):
    """Gauss–Legendre en cada tramo entre nodos: exacta para productos de cúbicas (grado 6 ≤ 2·4 − 1)."""
    u = np.unique(tau)
    x, w = np.polynomial.legendre.leggauss(n_gl)
    a, b = u[:-1, None], u[1:, None]
    return ((b - a) / 2 * x + (a + b) / 2).ravel(), ((b - a) / 2 * w).ravel()


def gram(tau):
    tq, wq = cuadratura(tau)
    P = base_en_nodos(tq, tau)
    return P.T @ (wq[:, None] * P)


# =============================================================================================
# 1. INSUMOS DE R Y VERIFICACIÓN
# =============================================================================================
titulo("1. Insumos de R y verificación de que Python reproduce lo mismo")
spec = json.loads((ENTRADA / "especificacion.json").read_text(encoding="utf-8"))
K_E, K_S, M, LAMBDA_R = spec["K_estaciones"], spec["K_satelite"], spec["M"], spec["lambda"]
tau_e, tau_s = np.array(spec["nodos_estaciones"]), np.array(spec["nodos_satelite"])
pares = pd.read_csv(ENTRADA / "pares.csv", dtype={"codigo": str})
leer = lambda n, **kw: pd.read_csv(ENTRADA / n, **kw).to_numpy(float)
C_e_raiz, C_s_raiz = leer("coef_estacion_raiz.csv"), leer("coef_satelite_raiz.csv")
C_e_mm, C_s_mm = leer("coef_estacion_mm.csv"), leer("coef_satelite_mm.csv")
W_e_R, W_s_R, R_s = (leer(n, header=None) for n in ("W_e.csv", "W_s.csv", "R_s.csv"))
fp = pd.read_csv(ENTRADA / "fpca_estaciones.csv")
mu_R, psi_R = fp["mu"].to_numpy(), fp.drop(columns="mu").to_numpy()
B_R = leer("coeficientes_regresion.csv", header=None)

Phi_e, Phi_s = base_en_nodos(DIAS, tau_e), base_en_nodos(DIAS, tau_s)
W_e, W_s = gram(tau_e), gram(tau_s)
ent = (pares.particion == "entrenamiento").to_numpy()
val = (pares.particion == "validación").to_numpy()
pru = (pares.particion == "prueba").to_numpy()
print(pares.groupby("particion").agg(pares=("anio", "size"), estaciones=("codigo", "nunique")))
prueba("Pares: 3 299 = 2 257 + 672 + 370 y todos usables", len(pares) == 3299 and (ent.sum(), val.sum(), pru.sum()) == (2257, 672, 370)
       and pares.usar.all())
prueba("Partición temporal sin cruces (1981-1995 / 1996-2000 / 2001-2003)",
       pares.anio[ent].max() <= 1995 and pares.anio[val].between(1996, 2000).all() and pares.anio[pru].min() >= 2001)
prueba("W calculada en Python por cuadratura = W de R", max(np.abs(W_e - W_e_R).max(), np.abs(W_s - W_s_R).max()) < 1e-10,
       f"dif. máx. {max(np.abs(W_e - W_e_R).max(), np.abs(W_s - W_s_R).max()):.1e}")


def fpca(C, W):
    """FPCA con producto interno de L²: autovectores de W^½ Σ W^½; ψ = W^-½ v (ψᵀWψ = I)."""
    mu = C.mean(0)
    S = np.cov(C - mu, rowvar=False)                   # divisor N − 1, como R
    ev, U = np.linalg.eigh(W)
    W_h, W_mh = U @ np.diag(ev ** .5) @ U.T, U @ np.diag(ev ** -.5) @ U.T
    lam, V = np.linalg.eigh(W_h @ S @ W_h)
    o = np.argsort(lam)[::-1]
    return mu, W_mh @ V[:, o], lam[o]


mu_py, psi_py, lam_py = fpca(C_e_raiz[ent], W_e)
signo = np.sign((psi_py[:, :M] * psi_R).sum(0))          # el signo de un autovector es arbitrario
dif_psi = np.abs(psi_py[:, :M] * signo - psi_R).max()
var_acum = np.cumsum(lam_py) / lam_py.sum()
tabla(pd.DataFrame({"componente": np.arange(1, K_E + 1), "autovalor": lam_py, "pct": 100 * lam_py / lam_py.sum(),
                    "acumulado_pct": 100 * var_acum}), "P1_fpca_varianza")
prueba("FPCA de Python (solo entrenamiento) = FPCA de R", np.abs(mu_py - mu_R).max() < 1e-10 and dif_psi < 1e-8,
       f"μ {np.abs(mu_py - mu_R).max():.1e}, ψ {dif_psi:.1e} (salvo signo)")
print(f"M = {M}: {100 * var_acum[M - 1]:.1f} % de la varianza (M y λ elegidos en validación por R)")
mu, D = mu_R, psi_R                                       # decodificador: el de R, congelado

alt = ((pares.altitud - spec["altitud_media"]) / spec["altitud_de"]).to_numpy()
oni = pares.ONI_anual.to_numpy()
prueba("Altitud estandarizada con media y DE de entrenamiento",
       abs(pares.altitud[ent].mean() - spec["altitud_media"]) < 1e-4 and abs(pares.altitud[ent].std() - spec["altitud_de"]) < 1e-4)

# Curvas observadas y del satélite sin corregir en mm/día (la referencia de toda comparación)
Y_mm = np.clip(C_e_mm @ Phi_e.T, 0, None)
X_mm = np.clip(C_s_mm @ Phi_s.T, 0, None)


def a_mm(C_raiz):
    y = C_raiz @ Phi_e.T
    return np.clip(y, 0, None) ** 2, y < 0


def pico(Y, v):
    return np.argmax(Y[:, v[0] - 1:v[1]], axis=1) + v[0]


def metricas(P, i, modelo, semilla=0):
    """Las cuatro métricas de la rúbrica 4.4 por curva, con las mismas fórmulas que la sección 5.5 de R."""
    Yo = Y_mm[i]
    r = np.array([np.corrcoef(a, b)[0, 1] for a, b in zip(P, Yo)])
    kge = 1 - np.sqrt((r - 1) ** 2 + (P.std(1, ddof=1) / Yo.std(1, ddof=1) - 1) ** 2 + (P.mean(1) / Yo.mean(1) - 1) ** 2)
    return pd.DataFrame({"modelo": modelo, "semilla": semilla, "codigo": pares.codigo[i].to_numpy(), "anio": pares.anio[i].to_numpy(),
                         "fase": pares.fase[i].to_numpy(), "particion": pares.particion[i].to_numpy(),
                         "l2_rel": np.sqrt(((P - Yo) ** 2).sum(1) / (Yo ** 2).sum(1)), "rmse": np.sqrt(((P - Yo) ** 2).mean(1)),
                         "kge": kge, "sesgo_pct": 100 * (P.mean(1) - Yo.mean(1)) / Yo.mean(1),
                         "err_pico_1": np.abs(pico(P, (60, 200)) - pico(Yo, (60, 200))),
                         "err_pico_2": np.abs(pico(P, (230, 350)) - pico(Yo, (230, 350)))})


METRICAS = ["l2_rel", "rmse", "kge", "sesgo_pct", "err_pico_1", "err_pico_2"]

# La regresión de R, evaluada en Python con sus coeficientes: debe dar las mismas métricas de prueba
Zx = C_s_raiz @ W_s                                       # ∫ φ_j(s) X(s) ds, j = 1..K_s
diseno_R = np.column_stack([np.ones(len(pares)), Zx, alt, oni])
C_R = mu + (diseno_R @ B_R) @ D.T                         # coeficientes √P reconstruidos por la regresión de R
i_pru = np.flatnonzero(pru)
m_R = metricas(a_mm(C_R[i_pru])[0], i_pru, "Regresión funcional (R)")
m_sat = metricas(X_mm[i_pru], i_pru, "Satélite sin corregir")
ref_R = pd.read_csv(ENTRADA.parent / "tablas" / "F5_metricas_prueba.csv")
punto = lambda s: float(str(s).split(" ")[0])
comp_R = pd.DataFrame({
    "metrica": ["L2_relativo", "KGE", "sesgo_anual_pct"],
    "R": [punto(ref_R.loc[ref_R.modelo == "Regresión funcional", c].iloc[0]) for c in ("L2_relativo", "KGE", "sesgo_anual_pct")],
    "Python": [m_R.groupby("codigo")[c].mean().mean() for c in ("l2_rel", "kge", "sesgo_pct")]})
print(tabla(comp_R, "P1_regresion_R_reproducida"))
comp_R["decimales_R"] = [2, 2, 1]
comp_R["coincide"] = [abs(r - round(p, d)) < 1e-9 for r, p, d in zip(comp_R.R, comp_R.Python, comp_R.decimales_R)]
prueba("Python reproduce las métricas de prueba de la regresión de R", comp_R.coincide.all(), "a los decimales que publica R")

# QUÉ DICEN LOS DATOS (sección 1)
#   Python lee exactamente lo que dejó R: la W por cuadratura coincide (< 1e-10), el FPCA
#   recalculado solo con entrenamiento coincide con el de R salvo el signo de cada ψ, y con
#   los coeficientes de R la regresión da las mismas métricas de prueba. Desde aquí, todo
#   lo que cambie en la red es efecto de la cabeza, no de los datos.

# =============================================================================================
# 2. LA RED FUNCIONAL
# =============================================================================================
titulo("2. La red funcional: codificador, decodificador y cabezas")
T = lambda a: torch.as_tensor(np.array(a, dtype=float))
W_s_t, W_e_t, R_s_t = T(W_s), T(W_e), T(R_s)
X_t, E_t, Y_t = T(C_s_raiz), T(np.column_stack([alt, oni])), T(C_e_raiz)
X_ent = X_t[torch.as_tensor(ent)]


class KANLayer(torch.nn.Module):
    """Una capa KAN: cada arista i→j aplica φ_ij(x) = w_b·SiLU(x) + w_s·Σ_g c_g B_g(x) (B-splines cúbicas
    uniformes en `rango`, G intervalos), y la salida j suma sus aristas: y_j = Σ_i φ_ij(x_i)."""

    def __init__(self, n_in, n_out, G=5, rango=(-3.0, 3.0), orden=3):
        super().__init__()
        self.G, self.k, self.rango = G, orden, rango
        h = (rango[1] - rango[0]) / G
        self.register_buffer("grid", torch.arange(-orden, G + orden + 1, dtype=torch.float64) * h + rango[0])
        self.coef = torch.nn.Parameter(0.1 * torch.randn(n_in, n_out, G + orden))
        self.w_b = torch.nn.Parameter(torch.randn(n_in, n_out) / np.sqrt(n_in))
        self.w_s = torch.nn.Parameter(torch.ones(n_in, n_out))

    def bases(self, x):
        g, x = self.grid, x.unsqueeze(-1)
        B = ((x >= g[:-1]) & (x < g[1:])).to(x.dtype)
        for p in range(1, self.k + 1):
            B = (x - g[:-(p + 1)]) / (g[p:-1] - g[:-(p + 1)]) * B[..., :-1] + (g[p + 1:] - x) / (g[p + 1:] - g[1:-p]) * B[..., 1:]
        return B

    def activaciones(self, x):                      # (n, n_in, n_out): φ_ij(x_i)
        spl = torch.einsum("nig,iog->nio", self.bases(x), self.coef)
        return self.w_b * torch.nn.functional.silu(x).unsqueeze(-1) + self.w_s * spl

    def forward(self, x):
        return self.activaciones(x).sum(1)


class KAN(torch.nn.Module):
    def __init__(self, ancho, G=5, rangos=None):
        super().__init__()
        rangos = rangos or [(-3.0, 3.0)] * (len(ancho) - 1)
        self.capas = torch.nn.ModuleList(KANLayer(a, b, G, r) for a, b, r in zip(ancho[:-1], ancho[1:], rangos))

    def forward(self, x):
        for c in self.capas:
            x = c(x)
        return x


def mlp(n_in, oculta, n_out):
    return torch.nn.Sequential(torch.nn.Linear(n_in, oculta), torch.nn.SiLU(), torch.nn.Linear(oculta, n_out))


def n_param(m):
    return int(sum(p.numel() for p in m.parameters()))


class RedFuncional(torch.nn.Module):
    """Codificador funcional + cabeza + decodificador FPCA congelado.
    z_k = a_kᵀ W_s c. Antes de la cabeza, z se estandariza con la media y la DE de z en las curvas de
    ENTRENAMIENTO (calculadas en cada paso con los pesos actuales): así la KAN recibe valores dentro de su
    rejilla y ningún dato de validación o prueba entra a la escala. Es afín, así que no rompe la
    equivalencia de la cabeza lineal con la regresión funcional."""

    def __init__(self, cabeza, n_neuronas):
        super().__init__()
        self.A = torch.nn.Parameter(0.1 * torch.randn(n_neuronas, K_S))
        self.cabeza = cabeza
        self.register_buffer("mu", T(mu))
        self.register_buffer("D", T(D))

    def z(self, c_sat):
        return c_sat @ W_s_t @ self.A.T

    def escala_z(self):
        z = self.z(X_ent)
        return z.mean(0), z.std(0)

    def entrada_cabeza(self, c_sat, extra):
        m, s = self.escala_z()
        return torch.cat([(self.z(c_sat) - m) / s, extra], 1)

    def forward(self, c_sat, extra):
        s = self.cabeza(self.entrada_cabeza(c_sat, extra))
        return self.mu + s @ self.D.T, s

    def rugosidad_relativa(self):
        """Σ_k ∫ w_k''² / ∫ w_k²: no cambia si se multiplica w_k por una constante (la estandarización de z
        tampoco), así la penalización no se puede "ganar" encogiendo las w_k."""
        return (torch.einsum("ki,ij,kj->k", self.A, R_s_t, self.A) / torch.einsum("ki,ij,kj->k", self.A, W_s_t, self.A)).sum()


def construir(tipo, semilla, cfg):
    torch.manual_seed(semilla)
    n_in = N_NEURONAS + 2
    cab = {"lineal": lambda: torch.nn.Linear(n_in, M), "kan": lambda: KAN([n_in, cfg["kan_oculta"], M], G=cfg["G"]),
           "mlp": lambda: mlp(n_in, cfg["mlp_oculta"], M)}[tipo]()
    torch.manual_seed(semilla)
    return RedFuncional(cab, N_NEURONAS)


def perdida_L2(red, idx):
    c_hat, _ = red(X_t[idx], E_t[idx])
    d = c_hat - Y_t[idx]
    return torch.einsum("ni,ij,nj->n", d, W_e_t, d).mean()


I_ENT, I_VAL = torch.as_tensor(np.flatnonzero(ent)), torch.as_tensor(np.flatnonzero(val))


def entrenar(tipo, semilla, lam_w, cfg):
    """Mismo presupuesto para todas las cabezas: Adam (lr fijo), lote completo de entrenamiento, a lo sumo
    EPOCAS_MAX épocas y parada temprana con la pérdida L² de VALIDACIÓN (paciencia PACIENCIA)."""
    red = construir(tipo, semilla, cfg)
    opt = torch.optim.Adam(red.parameters(), lr=LR)
    mejor, mejor_ep, estado = np.inf, 0, None
    t0 = time.time()
    for ep in range(EPOCAS_MAX):
        red.train()
        opt.zero_grad()
        l2 = perdida_L2(red, I_ENT)
        (l2 + lam_w * red.rugosidad_relativa()).backward()
        opt.step()
        red.eval()
        with torch.no_grad():
            lv = perdida_L2(red, I_VAL).item()
        if lv < mejor - 1e-10:
            mejor, mejor_ep, estado = lv, ep, {k: v.detach().clone() for k, v in red.state_dict().items()}
        elif ep - mejor_ep > PACIENCIA:
            break
    red.load_state_dict(estado)
    red.eval()
    with torch.no_grad():
        le = perdida_L2(red, I_ENT).item()
    return red, {"cabeza": tipo, "semilla": semilla, "lambda_w": lam_w, "epoca_optima": mejor_ep + 1, "epocas_corridas": ep + 1,
                 "L2_entrenamiento": le, "L2_validacion": mejor, "segundos": round(time.time() - t0, 1)}


def predecir(red, idx):
    with torch.no_grad():
        c, _ = red(X_t[idx], E_t[idx])
    return c.numpy()


# =============================================================================================
# 3. CABEZA LINEAL = REGRESIÓN FUNCIONAL
# =============================================================================================
titulo("3. Cabeza lineal = regresión funcional: demostración y verificación")
# DEMOSTRACIÓN. Con cabeza lineal, s_m = b_m + Σ_k B_mk (z_k − m_k)/σ_k + γ_m·alt + δ_m·ONI, y
# z_k = ∫ w_k(s) X(s) ds. Sustituyendo en ŷ(t) = μ(t) + Σ_m s_m ψ_m(t):
#   ŷ(t) = α(t) + ∫ β(s,t) X(s) ds + γ(t)·alt + δ(t)·ONI, con
#   β(s,t) = Σ_k Σ_m w_k(s) (B_mk/σ_k) ψ_m(t),   α(t) = μ(t) + Σ_m (b_m − Σ_k B_mk m_k/σ_k) ψ_m(t),
#   γ(t) = Σ_m γ_m ψ_m(t),   δ(t) = Σ_m δ_m ψ_m(t).
# Es una regresión funcional función a función con β(s,t) de rango ≤ min(N_NEURONAS, M) = M.
# Como M = 6 ≤ 8 neuronas, la regresión de R (β_m(s) libre en la base del satélite) cabe exactamente
# en la red: basta poner w_m(s) = β_m(s). Eso da una segunda verificación, más fuerte: la red con los
# coeficientes de R tiene que dar las predicciones de R.


def beta_de_red(red):
    """β(s,t), α(t), γ(t), δ(t) de una red con cabeza lineal, evaluadas en la cuadratura del satélite."""
    lin = red.cabeza
    Bw, b0 = lin.weight.detach().numpy(), lin.bias.detach().numpy()
    with torch.no_grad():
        m_z, s_z = (v.numpy() for v in red.escala_z())
    B_z, B_e = Bw[:, :N_NEURONAS], Bw[:, N_NEURONAS:]
    psi_t = D.T @ Phi_e.T                                   # (M, 365)
    tq, wq = cuadratura(tau_s)
    w_q = red.A.detach().numpy() @ base_en_nodos(tq, tau_s).T  # w_k(s) en la cuadratura
    beta = np.einsum("ks,mk,mt->st", w_q, B_z / s_z, psi_t)
    alfa = mu @ Phi_e.T + (b0 - B_z @ (m_z / s_z)) @ psi_t
    return beta, alfa, B_e.T @ psi_t, tq, wq


def prediccion_por_beta(red, idx):
    beta, alfa, ge, tq, wq = beta_de_red(red)
    X_q = C_s_raiz[idx] @ base_en_nodos(tq, tau_s).T        # X(s) en la cuadratura
    return alfa + (X_q * wq) @ beta + np.column_stack([alt, oni])[idx] @ ge


# (a) La red con los coeficientes de la regresión de R
def red_con_coeficientes_de_R():
    red = construir("lineal", SEMILLAS[0], {})
    with torch.no_grad():
        A = np.zeros((N_NEURONAS, K_S))
        A[:M] = B_R[1:1 + K_S].T                            # w_m(s) = β_m(s) de R
        A[M:] = np.random.default_rng(0).normal(size=(N_NEURONAS - M, K_S))   # neuronas sobrantes: peso 0 en la cabeza
        red.A.copy_(T(A))
        m_z, s_z = (v.numpy() for v in red.escala_z())
        Wl = np.zeros((M, N_NEURONAS + 2))
        Wl[:, :M] = np.diag(s_z[:M])                         # deshace la estandarización: s_m = z_m + ...
        Wl[:, N_NEURONAS:] = B_R[1 + K_S:].T                 # altitud y ONI
        red.cabeza.weight.copy_(T(Wl))
        red.cabeza.bias.copy_(T(B_R[0] + m_z[:M]))
    return red


red_R = red_con_coeficientes_de_R()
todos = np.arange(len(pares))
dif_puente = np.abs(predecir(red_R, todos) - C_R).max()
prueba("La red lineal con los coeficientes de R da las predicciones de R", dif_puente < 1e-8, f"dif. máx. {dif_puente:.1e} en coef. √P")

with torch.no_grad():
    L2_val_R = perdida_L2(red_R, I_VAL).item()
print(f"Pérdida L² de validación de la regresión de R: {L2_val_R:.2f} (referencia para la cabeza lineal entrenada)")

# QUÉ DICEN LOS DATOS (sección 3)
#   La red lineal cargada con los coeficientes de la regresión de R da sus mismas predicciones
#   (3,6e-15): el codificador y el decodificador de Python son la misma matemática que R. Las
#   tres redes lineales entrenadas cumplen la identidad con β(s,t) a 1,5e-14 (criterio 1e-4).
#   Entrenada, la cabeza lineal llega a una pérdida L² de validación de 124,9 a 126,5 según la
#   semilla, frente a 126,9 de la regresión de R: llega al mismo lugar por otro camino.
#   Con 6 neuronas la validación es igual (126,8 frente a 126,3 con 8), como predice el rango
#   ≤ M; se dejan 8 (el máximo del enunciado) para no limitar a las cabezas no lineales.

# =============================================================================================
# 4. KAN VERIFICADA Y MLP CON LOS MISMOS PARÁMETROS
# =============================================================================================
titulo("4. KAN verificada con una función sintética y MLP con los mismos parámetros")


def kan_sintetica(semilla=0):
    """Antes de los datos reales: la KAN [2,5,1] debe aprender f(x,y) = exp(sin(πx) + y²) en [−1,1]² con
    RMSE de prueba ≤ 1e-2 (rúbrica 4.3). Es la función de ejemplo de la documentación de pykan."""
    torch.manual_seed(semilla)
    rng = np.random.default_rng(semilla)
    f = lambda X: np.exp(np.sin(np.pi * X[:, 0]) + X[:, 1] ** 2)
    Xa, Xb = rng.uniform(-1, 1, (3000, 2)), rng.uniform(-1, 1, (1000, 2))
    xa, ya, xb, yb = T(Xa), T(f(Xa))[:, None], T(Xb), T(f(Xb))[:, None]
    kan = KAN([2, 5, 1], G=10, rangos=[(-1.0, 1.0), (-4.0, 4.0)])
    opt = torch.optim.Adam(kan.parameters(), lr=1e-2)
    for _ in range(3000):
        opt.zero_grad()
        ((kan(xa) - ya) ** 2).mean().backward()
        opt.step()
    lb = torch.optim.LBFGS(kan.parameters(), lr=1, max_iter=500, history_size=50, line_search_fn="strong_wolfe")

    def cierre():
        lb.zero_grad()
        l = ((kan(xa) - ya) ** 2).mean()
        l.backward()
        return l

    for _ in range(4):
        lb.step(cierre)
    with torch.no_grad():
        return float(torch.sqrt(((kan(xb) - yb) ** 2).mean()))


rmse_sint = kan_sintetica()
tabla(pd.DataFrame([{"funcion": "exp(sin(pi x) + y^2) en [-1,1]^2", "arquitectura": "KAN [2,5,1], G = 10", "RMSE_prueba": rmse_sint,
                     "umbral": 1e-2}]), "P4_kan_sintetica")
prueba("KAN sintética: RMSE de prueba ≤ 1e-2", rmse_sint <= 1e-2, f"RMSE {rmse_sint:.1e}")

n_in = N_NEURONAS + 2
CFG = {"kan_oculta": 4, "G": 5}
n_kan = n_param(KAN([n_in, CFG["kan_oculta"], M], G=CFG["G"]))
CFG["mlp_oculta"] = min(range(2, 200), key=lambda h: abs(n_param(mlp(n_in, h, M)) - n_kan))
n_mlp = n_param(mlp(n_in, CFG["mlp_oculta"], M))
n_lin, n_cod = n_param(torch.nn.Linear(n_in, M)), N_NEURONAS * K_S
param = tabla(pd.DataFrame([
    {"cabeza": "lineal", "arquitectura": f"Linear({n_in}→{M})", "param_cabeza": n_lin},
    {"cabeza": "kan", "arquitectura": f"KAN[{n_in},{CFG['kan_oculta']},{M}], G = {CFG['G']}, cúbica", "param_cabeza": n_kan},
    {"cabeza": "mlp", "arquitectura": f"MLP[{n_in},{CFG['mlp_oculta']},{M}], SiLU", "param_cabeza": n_mlp}]).assign(
    param_codificador=n_cod, param_total=lambda d: d.param_cabeza + n_cod), "P4_parametros_cabezas")
print(param)
dif_par = abs(n_kan - n_mlp) / n_kan
prueba("KAN y MLP con parámetros igualados (sum(p.numel()), ≤ 10 %)", dif_par <= .10, f"{n_kan} vs {n_mlp}: {100 * dif_par:.1f} %")

# QUÉ DICEN LOS DATOS (sección 4)
#   La KAN programada aquí aprende la función de prueba con RMSE 2e-3 (criterio 1e-2), así que
#   la capa funciona antes de tocar datos reales. KAN [10,4,6] con G = 5 tiene 640 parámetros;
#   el MLP más cercano es [10,37,6] con 635 (0,8 % de diferencia). El codificador (88) es igual.

# =============================================================================================
# 5. ABLACIÓN: TRES CABEZAS, MISMO PRESUPUESTO, TRES SEMILLAS
# =============================================================================================
titulo("5. Ablación de la cabeza")
# Presupuesto idéntico: cada cabeza prueba los mismos 4 valores de λ_w con la primera semilla y se queda
# con el de menor L² de validación (4 pruebas de hiperparámetros por cabeza); luego 3 semillas.
seleccion, final, redes = [], [], {}
for tipo in ("lineal", "kan", "mlp"):
    for lam_w in LAMBDAS_W:
        _, info = entrenar(tipo, SEMILLAS[0], lam_w, CFG)
        seleccion.append(info)
        print(tipo, {k: info[k] for k in ("lambda_w", "L2_validacion", "epoca_optima", "segundos")}, flush=True)
    lam_t = min((s for s in seleccion if s["cabeza"] == tipo), key=lambda s: s["L2_validacion"])["lambda_w"]
    for s in SEMILLAS:
        red, info = entrenar(tipo, s, lam_t, CFG)
        redes[(tipo, s)] = red
        final.append({**info, "param_total": n_param(red)})
        print(tipo, s, {k: info[k] for k in ("L2_validacion", "epoca_optima", "segundos")}, flush=True)
tabla(pd.DataFrame(seleccion), "P5_seleccion_lambda_por_cabeza")
LAM_LINEAL = min((s for s in seleccion if s["cabeza"] == "lineal"), key=lambda s: s["L2_validacion"])["lambda_w"]

# La identidad de la sección 3 para las redes lineales ENTRENADAS, en las curvas de prueba (≤ 1e-4)
eq = []
for s in SEMILLAS:
    red = redes[("lineal", s)]
    eq.append({"verificacion": f"red lineal entrenada (semilla {s}) vs ∫β(s,t)X(s)ds por cuadratura, prueba",
               "dif_max": np.abs(predecir(red, i_pru) @ Phi_e.T - prediccion_por_beta(red, i_pru)).max(), "umbral": 1e-4})
eq = pd.DataFrame([{"verificacion": "red lineal con los coeficientes de R vs predicciones de R", "dif_max": dif_puente, "umbral": 1e-8}] + eq)
print(tabla(eq, "P3_equivalencia_lineal").to_string())
prueba("Cabeza lineal entrenada = regresión funcional con β(s,t) (≤ 1e-4, 3 semillas)", (eq.dif_max.iloc[1:] <= 1e-4).all(),
       f"dif. máx. {eq.dif_max.iloc[1:].max():.1e} √(mm/día)")

# ¿6 u 8 neuronas? Con la cabeza lineal el rango de β(s,t) es ≤ M = 6, así que 6 neuronas bastan; las 2
# extra solo sirven a las cabezas no lineales. Se comprueba que a la lineal no le cambian nada.
N_NEURONAS = 6
neur = [{"neuronas": 6, **entrenar("lineal", SEMILLAS[0], LAM_LINEAL, CFG)[1]}]
N_NEURONAS = 8
neur.append({"neuronas": 8, **[f for f in final if f["cabeza"] == "lineal" and f["semilla"] == SEMILLAS[0]][0]})
neur = tabla(pd.DataFrame(neur), "P3_neuronas_6_vs_8_lineal")
print(neur[["neuronas", "lambda_w", "L2_entrenamiento", "L2_validacion"]], f"\nRegresión de R: L² de validación {L2_val_R:.2f}")
final = tabla(pd.DataFrame(final), "P5_entrenamiento_por_semilla")
tabla(pd.DataFrame([{"optimizador": "Adam", "lr": LR, "lote": "completo (entrenamiento)", "epocas_max": EPOCAS_MAX,
                     "parada_temprana": f"L² de validación, paciencia {PACIENCIA}", "pruebas_hiperparametros_por_cabeza": len(LAMBDAS_W),
                     "rejilla_lambda_w": str(LAMBDAS_W), "semillas": str(SEMILLAS)}]), "P5_presupuesto")

# Métricas en validación y prueba, contra el satélite sin corregir y la regresión de R
filas, recorte = [], []
for part, mascara in (("validación", val), ("prueba", pru)):
    idx = np.flatnonzero(mascara)
    filas += [metricas(X_mm[idx], idx, "satélite sin corregir"), metricas(a_mm(C_R[idx])[0], idx, "regresión (R)")]
    for (tipo, s), red in redes.items():
        P, neg = a_mm(predecir(red, idx))
        filas.append(metricas(P, idx, tipo, s))
        recorte.append({"particion": part, "cabeza": tipo, "semilla": s, "pct_valores_recortados": 100 * neg.mean(),
                        "pct_curvas_con_recorte": 100 * neg.any(1).mean()})
met = pd.concat(filas, ignore_index=True)
met.to_csv(TAB / "P5_metricas_por_curva.csv", index=False)
tabla(pd.DataFrame(recorte).groupby(["particion", "cabeza"])[["pct_valores_recortados", "pct_curvas_con_recorte"]].mean().round(3).reset_index(),
      "P5_recorte_no_negatividad")

# Media ± DE entre semillas (cada semilla: media de las medias por estación, como en R)
por_sem = met.groupby(["particion", "modelo", "semilla", "codigo"])[METRICAS].mean().groupby(["particion", "modelo", "semilla"]).mean()
res = por_sem.groupby(["particion", "modelo"]).agg(["mean", "std"]).round(3)
res.columns = [f"{a}_{b}" for a, b in res.columns]
print(tabla(res.reset_index(), "P5_ablacion_media_DE_semillas").to_string())

# Comparaciones pareadas por estación en PRUEBA (métricas promediadas sobre las 3 semillas):
# IC 95 % bootstrap de la diferencia (2 000 remuestreos de estaciones) y Wilcoxon pareado.
rng = np.random.default_rng(20261114)
pe = met[met.particion == "prueba"].groupby(["modelo", "codigo"])[METRICAS].mean()
comp = []
for a, b in (("kan", "lineal"), ("mlp", "lineal"), ("kan", "mlp"), ("lineal", "satélite sin corregir"),
             ("kan", "satélite sin corregir"), ("mlp", "satélite sin corregir"), ("lineal", "regresión (R)")):
    for mt in METRICAS:
        va, vb = pe.loc[a, mt], pe.loc[b, mt]
        if mt == "sesgo_pct":
            va, vb = va.abs(), vb.abs()
        d = (va - vb).to_numpy()
        boot = d[rng.integers(0, len(d), (N_BOOT, len(d)))].mean(1)
        q = np.quantile(boot, [.025, .975])
        mejor_es_mayor = mt == "kge"
        comp.append({"comparacion": f"{a} − {b}", "metrica": "|sesgo_pct|" if mt == "sesgo_pct" else mt, "diferencia": d.mean(),
                     "ic_inf": q[0], "ic_sup": q[1], "ic_excluye_0": q[0] > 0 or q[1] < 0,
                     "gana": (a if (d.mean() > 0) == mejor_es_mayor else b) if (q[0] > 0 or q[1] < 0) else "ninguna (IC con 0)",
                     "p_wilcoxon": stats.wilcoxon(va, vb).pvalue if np.any(d != 0) else 1.0})
comp = tabla(pd.DataFrame(comp).round(4), "P5_comparaciones_pareadas_prueba")
print(comp.to_string())

# QUÉ DICEN LOS DATOS (sección 5)
#   Elección de λ_w (validación, semilla 2026): lineal 1e8, KAN 1e6, MLP 1e7.
#   Validación (pérdida L² en √P): lineal 124,9-126,5; KAN 106,3-124,3; MLP 100,0-108,2. Las
#   no lineales se detienen hacia la época 100 (la lineal hacia la 2 000-3 800) y quedan con una
#   brecha entrenamiento-validación de 15 puntos: aprenden rápido y sobreajustan rápido. La KAN
#   es menos estable entre semillas (una de las tres queda en 124).
#   Prueba (2001-2003, 370 curvas, 137 estaciones; media ± DE de 3 semillas):
#     L² relativo  satélite 0,513 | lineal 0,440 ± 0,002 | KAN 0,434 ± 0,002 | MLP 0,420 ± 0,010
#     KGE          satélite 0,382 | lineal 0,454 ± 0,005 | KAN 0,471 ± 0,009 | MLP 0,485 ± 0,020
#     sesgo anual  satélite +26 % | lineal +9,8 %        | KAN +9,5 %        | MLP +10,2 %
#   Las tres cabezas mejoran al satélite sin corregir en L², KGE, sesgo y pico 1 (IC 95 % de la
#   diferencia lejos del 0). Entre cabezas (IC 95 % bootstrap pareado por estación):
#     - MLP supera a la lineal en L² relativo (−0,020 [−0,040, −0,000]) y en KGE (+0,030
#       [0,005, 0,056]), por poco; la lineal ubica mejor el pico 1 (0,7 días).
#     - KAN no se distingue de la lineal en L² ni en KGE; reduce el |sesgo| 4 puntos y ubica
#       peor el pico 2 (2 días).
#     - KAN frente a MLP: no se distinguen en L², KGE, sesgo ni pico 1; el MLP ubica mejor el
#       pico 2 (1,8 días).
#   Lectura: hay algo de no linealidad en el sesgo (el MLP gana), pero es pequeña y se ve sobre
#   todo en validación, que tiene un Niño fuerte (1997) y dos Niñas (1999, 2000); la prueba
#   solo tiene años neutros y un Niño moderado (2002). El
#   costo de la interpretabilidad (KAN frente a MLP) no es significativo en prueba. Ninguna
#   predicción queda negativa en la escala √P (0 % de valores recortados).

# =============================================================================================
# 6. PRUEBAS Y CIERRE
# =============================================================================================
titulo("6. Pruebas y cierre")
prueba("Las tres cabezas tienen el mismo codificador y el mismo decodificador",
       len({(r.A.shape, tuple(r.D.shape)) for r in redes.values()}) == 1 and all(torch.equal(r.D, redes[("lineal", SEMILLAS[0])].D) for r in redes.values()))
prueba("Tres semillas por cabeza", final.groupby("cabeza").semilla.nunique().eq(3).all())
prueba("Reconstrucción en mm no negativa", all((a_mm(predecir(r, i_pru))[0] >= 0).all() for r in redes.values()))
prueba("Bootstrap con al menos 1 000 remuestreos", N_BOOT >= 1000)
pd.DataFrame([{"prueba": k, "pasa": v[0], "detalle": v[1]} for k, v in PRUEBAS.items()]).to_csv(TAB / "P6_pruebas.csv", index=False)
(SALIDA / "versiones.txt").write_text(f"python {sys.version.split()[0]} | {platform.platform()}\nnumpy {np.__version__} | pandas {pd.__version__} | "
                                      f"torch {torch.__version__} | scipy {__import__('scipy').__version__}\n", encoding="utf-8")
n_ok = sum(v[0] for v in PRUEBAS.values())
print(f"\nPruebas: {n_ok} de {len(PRUEBAS)} pasan")
assert n_ok == len(PRUEBAS), "Alguna prueba falla"

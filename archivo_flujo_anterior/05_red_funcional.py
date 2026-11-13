"""Fase 4: red neuronal funcional, ablación de la cabeza, evaluación e interpretación.

Reconstruye la curva anual de una estación a partir de la curva satelital de su píxel.

    codificador funcional   z_k = (1/|T|) ∫ w_k(t) X(t) dt = a_kᵀ W_sat c / |T|,  k = 1..8
                            + altitud de la estación (estandarizada con entrenamiento) + ONI anual
    cabeza (ablación)       lineal | KAN | MLP      ->  M puntajes
    decodificador FPCA      ŷ(t) = μ(t) + Σ_m s_m ψ_m(t)   (μ, ψ solo con entrenamiento, congelados)
    pérdida                 (ĉ − c)ᵀ W_est (ĉ − c) + λ Σ_k a_kᵀ R a_k,  R = ∫ φ''(t) φ''(t)ᵀ dt

Todo se modela en escala √P (coeficientes de `coeficientes_transformados.npz`) y la curva final en mm
es max(ŷ, 0)²; se reporta el porcentaje de puntos recortados.

Partición (fe de erratas): entrenamiento 1981–1995, validación 1996–2000, prueba 2001–2003.
La validación solo se usa para elegir la época (parada temprana) y λ; la prueba solo se evalúa al final.
Si existe `estaciones_reservadas.csv` (columna station_id), esas estaciones se excluyen de
entrenamiento y validación y se evalúan aparte.

Librerías del apartado de modelos: PyTorch (autorizada). La KAN se programa aquí con la forma del
enunciado: φ(x) = w_b·SiLU(x) + w_s·Σ_i c_i B_i(x) en cada arista.

Uso, desde la raíz del repo:  python 05_red_funcional.py
Salidas: fda_fase4/tablas, fda_fase4/figuras, fda_fase4/predicciones.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from scipy import stats
from scipy.optimize import curve_fit

torch.set_default_dtype(torch.float64)
torch.set_num_threads(4)

ROOT = Path(__file__).resolve().parent
F1, F2 = ROOT / "fda_fase1", ROOT / "fda_fase2"
OUT = ROOT / "fda_fase4"
TAB, FIG, PRED = OUT / "tablas", OUT / "figuras", OUT / "predicciones"
for d in (TAB, FIG, PRED):
    d.mkdir(parents=True, exist_ok=True)

ENTREN_FIN, VALID_FIN = 1995, 2000
SEMILLAS = (2026, 2027, 2028)
N_NEURONAS = 8
EPOCAS, PACIENCIA, LR = 6000, 500, 3e-3
# Penalización de rugosidad sobre R normalizada (R / promedio de su diagonal), así λ·aᵀR̃a es comparable a ‖a‖²;
# λ se elige con validación usando la cabeza lineal.
LAMBDAS = (1.0, 10.0, 100.0, 1e3, 1e4, 1e5)
N_BOOT = 2000
T_LEN = 365.0
DIAS = np.arange(1, 366, dtype=float)
VENTANAS_PICO = ((60, 200), (230, 365))   # mismas ventanas que el EDA (fechas candidatas)


# ---------------------------------------------------------------------------------------------
# Base B-spline (Cox–de Boor, igual que la fase 1), Gram y rugosidad
# ---------------------------------------------------------------------------------------------
def nodos(rango, K, orden=4):
    a, b = rango
    return np.r_[[a] * orden, np.linspace(a, b, K - orden + 2)[1:-1], [b] * orden]


def base_bspline(t, rango, K, orden=4):
    t = np.asarray(t, float)
    tau = nodos(rango, K, orden)
    B = np.zeros((len(t), len(tau) - 1))
    for j in range(len(tau) - 1):
        B[:, j] = (t >= tau[j]) & (t < tau[j + 1])
    B[t == rango[1], np.max(np.flatnonzero(tau[:-1] < tau[1:]))] = 1
    for m in range(2, orden + 1):
        Bn = np.zeros((len(t), len(tau) - m))
        for j in range(len(tau) - m):
            d1, d2 = tau[j + m - 1] - tau[j], tau[j + m] - tau[j + 1]
            Bn[:, j] = ((t - tau[j]) / d1 * B[:, j] if d1 > 0 else 0) + \
                       ((tau[j + m] - t) / d2 * B[:, j + 1] if d2 > 0 else 0)
        B = Bn
    return B


def cuadratura(K, rango=(1, 365), n_gl=4):
    tau = np.unique(nodos(rango, K))
    x, w = np.polynomial.legendre.leggauss(n_gl)
    a, b = tau[:-1, None], tau[1:, None]
    return ((b - a) / 2 * x + (a + b) / 2).ravel(), ((b - a) / 2 * w).ravel()


def gram(K):
    tq, wq = cuadratura(K)
    P = base_bspline(tq, (1, 365), K)
    return P.T @ (wq[:, None] * P)


def rugosidad(K, h=1e-3):
    """R = ∫ φ''φ''ᵀ: dentro de cada intervalo la base es cúbica y la diferencia central es exacta."""
    tq, wq = cuadratura(K)
    d2 = (base_bspline(tq + h, (1, 365), K) - 2 * base_bspline(tq, (1, 365), K)
          + base_bspline(tq - h, (1, 365), K)) / h ** 2
    return d2.T @ (wq[:, None] * d2)


# ---------------------------------------------------------------------------------------------
# Datos: pares estación–píxel, coeficientes √P y covariables
# ---------------------------------------------------------------------------------------------
def cargar_datos():
    spec = json.loads((F1 / "base_especificacion.json").read_text(encoding="utf-8"))
    K_est, K_sat = spec["K_estaciones"], spec["K_satelite"]
    ce = np.load(F1 / "coeficientes_estaciones.npz")
    cs = np.load(F1 / "coeficientes_satelite.npz")
    tr = np.load(F1 / "coeficientes_transformados.npz")
    est = pd.DataFrame({"STATION_ID": ce["station_id"].astype(str), "YEAR": ce["anio"], "i_est": np.arange(len(ce["anio"]))})
    sat = pd.DataFrame({"YEAR": cs["anio"], "fila": cs["fila"], "col": cs["col"], "i_sat": np.arange(len(cs["anio"]))})

    pares = pd.read_csv(F2 / "tablas_R" / "pares_estacion_pixel.csv", dtype={"STATION_ID": str})
    pares = pares[(pares.estado == "pareada") & pares.altitud_m.notna()].copy()
    pares = pares.merge(est, on=["STATION_ID", "YEAR"], how="inner").merge(sat, on=["YEAR", "fila", "col"], how="inner")
    oni = pd.read_csv(ROOT / "eda_fase0" / "oni_anual_fase.csv").rename(columns={"anio": "YEAR"})
    pares = pares.merge(oni, on="YEAR", how="left")
    assert pares.ONI_anual.notna().all(), "Falta ONI para algún año"
    pares["particion"] = np.select([pares.YEAR <= ENTREN_FIN, pares.YEAR <= VALID_FIN], ["entrenamiento", "validacion"], "prueba")

    reservadas_path = ROOT / "estaciones_reservadas.csv"
    reservadas = set()
    if reservadas_path.exists():
        reservadas = set(pd.read_csv(reservadas_path, dtype=str).iloc[:, 0])
        pares.loc[pares.STATION_ID.isin(reservadas), "particion"] = "reservada"
    pares = pares.sort_values(["particion", "STATION_ID", "YEAR"]).reset_index(drop=True)

    datos = {
        "pares": pares, "K_est": K_est, "K_sat": K_sat, "reservadas": sorted(reservadas),
        "c_est_raiz": tr["raiz_est"][pares.i_est.to_numpy()].astype(float),
        "c_sat_raiz": tr["raiz_sat"][pares.i_sat.to_numpy()].astype(float),
        "c_est_mm": ce["coef"][pares.i_est.to_numpy()].astype(float),
        "c_sat_mm": cs["coef"][pares.i_sat.to_numpy()].astype(float),
        "W_est": gram(K_est), "W_sat": gram(K_sat), "R_sat": rugosidad(K_sat) / np.diag(rugosidad(K_sat)).mean(),
        "Phi_est": base_bspline(DIAS, (1, 365), K_est), "Phi_sat": base_bspline(DIAS, (1, 365), K_sat),
    }
    W_guardada = np.load(F1 / "matrices_gram.npz")
    assert np.abs(datos["W_est"] - W_guardada["W_est"]).max() < 1e-10
    assert np.abs(datos["W_sat"] - W_guardada["W_sat"]).max() < 1e-10
    return datos


# ---------------------------------------------------------------------------------------------
# FPCA funcional (solo entrenamiento)
# ---------------------------------------------------------------------------------------------
def fpca(C, W):
    """Devuelve μ (coef.), D (coef. de ψ con DᵀWD = I) y autovalores, con C = coeficientes de entrenamiento."""
    mu = C.mean(0)
    S = np.cov(C - mu, rowvar=False, bias=True)
    ev, U = np.linalg.eigh(W)
    W_mhalf, W_half = U @ np.diag(ev ** -.5) @ U.T, U @ np.diag(ev ** .5) @ U.T
    lam, V = np.linalg.eigh(W_half @ S @ W_half)
    orden = np.argsort(lam)[::-1]
    lam, V = lam[orden], V[:, orden]
    D = W_mhalf @ V
    for m in range(D.shape[1]):                 # signo: que ψ_m tenga integral positiva (solo convención gráfica)
        if (np.ones(len(mu)) @ W @ D[:, m]) < 0:
            D[:, m] *= -1
    return mu, D, np.clip(lam, 0, None)


# ---------------------------------------------------------------------------------------------
# Capas: KAN programada con la forma del enunciado
# ---------------------------------------------------------------------------------------------
class KANLayer(torch.nn.Module):
    """φ_ij(x) = w_b·SiLU(x) + w_s·Σ_g c_g B_g(x); salida_j = Σ_i φ_ij(x_i). B-splines cúbicas uniformes."""

    def __init__(self, n_in, n_out, G=5, rango=(-3.0, 3.0), orden=3, semilla_ruido=0.1):
        super().__init__()
        self.n_in, self.n_out, self.G, self.k = n_in, n_out, G, orden
        h = (rango[1] - rango[0]) / G
        self.register_buffer("grid", torch.arange(-orden, G + orden + 1, dtype=torch.float64) * h + rango[0])
        self.coef = torch.nn.Parameter(semilla_ruido * torch.randn(n_in, n_out, G + orden))
        self.w_b = torch.nn.Parameter(torch.randn(n_in, n_out) / np.sqrt(n_in))
        self.w_s = torch.nn.Parameter(torch.ones(n_in, n_out))
        self.rango = rango

    def bases(self, x):                      # x: (n, n_in) -> (n, n_in, G+k)
        g = self.grid
        x = x.unsqueeze(-1)
        B = ((x >= g[:-1]) & (x < g[1:])).to(x.dtype)
        for p in range(1, self.k + 1):
            B = (x - g[:-(p + 1)]) / (g[p:-1] - g[:-(p + 1)]) * B[..., :-1] + \
                (g[p + 1:] - x) / (g[p + 1:] - g[1:-p]) * B[..., 1:]
        return B

    def activaciones(self, x):               # (n, n_in, n_out)
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


# ---------------------------------------------------------------------------------------------
# Red funcional completa
# ---------------------------------------------------------------------------------------------
class RedFuncional(torch.nn.Module):
    def __init__(self, cabeza, K_sat, W_sat, mu, D_M, n_neuronas=N_NEURONAS):
        super().__init__()
        self.A = torch.nn.Parameter(0.1 * torch.randn(n_neuronas, K_sat))          # coeficientes de w_k(t)
        self.register_buffer("W_sat", torch.as_tensor(W_sat))
        # Normalización de z con media/varianza del conjunto de entrenamiento (momentum=1: el último paso
        # completo sobre entrenamiento fija los estadísticos; ver refrescar_norm)
        self.norm = torch.nn.BatchNorm1d(n_neuronas, affine=False, momentum=1.0)
        self.cabeza = cabeza
        self.register_buffer("mu", torch.as_tensor(mu))
        self.register_buffer("D_M", torch.as_tensor(D_M))

    def codificar(self, c_sat):
        return c_sat @ self.W_sat @ self.A.T / T_LEN

    def entrada_cabeza(self, c_sat, extra):
        return torch.cat([self.norm(self.codificar(c_sat)), extra], 1)

    def forward(self, c_sat, extra):
        s = self.cabeza(self.entrada_cabeza(c_sat, extra))
        return self.mu + s @ self.D_M.T, s

    @torch.no_grad()
    def refrescar_norm(self, c_sat_entren):
        """Fija media/varianza de z con TODO el entrenamiento y los pesos actuales (sin tocar validación/prueba)."""
        modo = self.training
        self.norm.train()
        self.norm(self.codificar(c_sat_entren))
        self.train(modo)


def construir_cabeza(tipo, n_in, M, semilla, oculta_mlp=None, kan_ancho=None, G=5):
    torch.manual_seed(semilla)
    if tipo == "lineal":
        return torch.nn.Linear(n_in, M)
    if tipo == "kan":
        return KAN([n_in, *kan_ancho, M], G=G)
    return mlp(n_in, oculta_mlp, M)


def entrenar(tipo, datos, tensores, M, mu, D_M, lam, semilla, cfg):
    torch.manual_seed(semilla)
    np.random.seed(semilla)
    n_in = N_NEURONAS + 2
    cabeza = construir_cabeza(tipo, n_in, M, semilla, cfg.get("oculta_mlp"), cfg.get("kan_ancho"), cfg.get("G", 5))
    torch.manual_seed(semilla)
    red = RedFuncional(cabeza, datos["K_sat"], datos["W_sat"], mu, D_M)
    W_e = torch.as_tensor(datos["W_est"])
    R = torch.as_tensor(datos["R_sat"])
    opt = torch.optim.Adam(red.parameters(), lr=LR)
    (Xtr, Etr, Ytr), (Xva, Eva, Yva) = tensores["entrenamiento"], tensores["validacion"]

    def perdida_L2(red_, X, E, Y):
        c_hat, _ = red_(X, E)
        d = c_hat - Y
        return torch.einsum("ni,ij,nj->n", d, W_e, d).mean()

    mejor, mejor_ep, estado, hist = np.inf, 0, None, []
    for ep in range(EPOCAS):
        red.train()
        opt.zero_grad()
        l2 = perdida_L2(red, Xtr, Etr, Ytr)
        rug = lam * torch.einsum("ki,ij,kj->", red.A, R, red.A)
        (l2 + rug).backward()
        opt.step()
        red.refrescar_norm(Xtr)
        red.eval()
        with torch.no_grad():
            lva = perdida_L2(red, Xva, Eva, Yva).item()
        hist.append((ep, l2.item(), lva))
        if lva < mejor - 1e-9:
            mejor, mejor_ep = lva, ep
            estado = {k: v.detach().clone() for k, v in red.state_dict().items()}
        elif ep - mejor_ep > PACIENCIA:
            break
    red.load_state_dict(estado)
    red.eval()
    return red, {"epoca_optima": mejor_ep, "L2_validacion": mejor, "epocas_corridas": ep + 1}, hist


# ---------------------------------------------------------------------------------------------
# Métricas
# ---------------------------------------------------------------------------------------------
def a_mm(c_raiz, Phi):
    y = c_raiz @ Phi.T
    return np.clip(y, 0, None) ** 2, (y < 0)


def picos(y):
    out = []
    for a, b in VENTANAS_PICO:
        m = (DIAS >= a) & (DIAS <= b)
        out.append(DIAS[m][np.argmax(y[:, m], axis=1)])
    return np.column_stack(out)


def metricas(y_hat, y):
    """Por curva: L² relativo, KGE, sesgo del acumulado anual (%) y error medio en la fecha de los dos picos (días)."""
    l2 = np.sqrt(((y_hat - y) ** 2).sum(1)) / np.sqrt((y ** 2).sum(1))
    r = np.array([np.corrcoef(a, b)[0, 1] if a.std() > 0 else 0.0 for a, b in zip(y_hat, y)])
    alfa = y_hat.std(1) / y.std(1)
    beta = y_hat.mean(1) / y.mean(1)
    kge = 1 - np.sqrt((r - 1) ** 2 + (alfa - 1) ** 2 + (beta - 1) ** 2)
    sesgo = 100 * (y_hat.sum(1) - y.sum(1)) / y.sum(1)
    pico = np.abs(picos(y_hat) - picos(y)).mean(1)
    return pd.DataFrame({"L2_rel": l2, "KGE": kge, "sesgo_pct": sesgo, "error_pico_dias": pico})


def media_por_estacion(df, cols):
    return df.groupby("STATION_ID")[cols].mean()


def bootstrap_estaciones(valores_por_estacion, B=N_BOOT, semilla=20261111):
    """IC 95 % de la media remuestreando estaciones (conglomerados) con reemplazo."""
    v = np.asarray(valores_por_estacion, float)
    rng = np.random.default_rng(semilla)
    idx = rng.integers(0, len(v), size=(B, len(v)))
    medias = v[idx].mean(1)
    return float(v.mean()), float(np.quantile(medias, .025)), float(np.quantile(medias, .975))


# ---------------------------------------------------------------------------------------------
# Verificaciones exigidas
# ---------------------------------------------------------------------------------------------
def verificar_kan_sintetica(semilla=0):
    """KAN [2,5,1] sobre f(x,y) = exp(sin(πx) + y²); exige RMSE de prueba ≤ 1e-2."""
    torch.manual_seed(semilla)
    rng = np.random.default_rng(semilla)
    f = lambda X: np.exp(np.sin(np.pi * X[:, 0]) + X[:, 1] ** 2)
    Xtr, Xte = rng.uniform(-1, 1, (3000, 2)), rng.uniform(-1, 1, (1000, 2))
    xt, yt = torch.as_tensor(Xtr), torch.as_tensor(f(Xtr))[:, None]
    xe, ye = torch.as_tensor(Xte), torch.as_tensor(f(Xte))[:, None]
    kan = KAN([2, 5, 1], G=10, rangos=[(-1.0, 1.0), (-4.0, 4.0)])
    opt = torch.optim.Adam(kan.parameters(), lr=1e-2)
    for ep in range(6000):
        opt.zero_grad()
        loss = ((kan(xt) - yt) ** 2).mean()
        loss.backward()
        opt.step()
    lbfgs = torch.optim.LBFGS(kan.parameters(), lr=1, max_iter=500, history_size=50, line_search_fn="strong_wolfe")

    def cierre():
        lbfgs.zero_grad()
        l = ((kan(xt) - yt) ** 2).mean()
        l.backward()
        return l

    for _ in range(4):
        lbfgs.step(cierre)
    with torch.no_grad():
        rmse = float(torch.sqrt(((kan(xe) - ye) ** 2).mean()))
    return {"funcion": "exp(sin(pi x) + y^2) en [-1,1]^2", "arquitectura": "[2,5,1], G=10",
            "RMSE_prueba": rmse, "umbral": 1e-2, "cumple": rmse <= 1e-2}


def verificar_equivalencia_lineal(red, datos, X, E, D_M, mu):
    """La red con cabeza lineal = regresión funcional ŷ(t) = α(t) + ∫ β(s,t) X(s) ds + γ(t)ᵀ e,
    β(s,t) = Σ_k Σ_m w_k(s) B̃_mk ψ_m(t), con B̃ = B·diag(1/(|T|·σ_k)) por la normalización de z.
    Se evalúa con cuadratura de Gauss–Legendre (exacta para el producto de B-splines)."""
    with torch.no_grad():
        c_red, _ = red(X, E)
    y_red = c_red.numpy() @ datos["Phi_est"].T
    lin = red.cabeza
    Bw = lin.weight.detach().numpy()            # (M, 8 + 2)
    b0 = lin.bias.detach().numpy()
    m_z = red.norm.running_mean.numpy()
    s_z = np.sqrt(red.norm.running_var.numpy() + red.norm.eps)
    B_z, B_e = Bw[:, :N_NEURONAS], Bw[:, N_NEURONAS:]
    B_til = B_z / (T_LEN * s_z)                 # (M, 8)
    psi = D_M.T @ datos["Phi_est"].T            # (M, 365): ψ_m(t)
    tq, wq = cuadratura(datos["K_sat"])
    Phi_q = base_bspline(tq, (1, 365), datos["K_sat"])
    w_s = red.A.detach().numpy() @ Phi_q.T      # (8, nq): w_k(s)
    beta = np.einsum("ks,mk,mt->st", w_s, B_til, psi)      # β(s,t)
    X_s = X.numpy() @ Phi_q.T                                # X(s) en los nodos de cuadratura
    integral = (X_s * wq) @ beta                             # ∫ β(s,t) X(s) ds
    alfa = mu @ datos["Phi_est"].T + (b0 - B_z @ (m_z / s_z)) @ psi
    gamma = B_e.T @ psi                                      # (2, 365)
    y_rf = alfa + integral + E.numpy() @ gamma
    return float(np.abs(y_red - y_rf).max()), beta, tq


def ajuste_simbolico(x, y):
    candidatos = {
        "a*x + b": (lambda x, a, b: a * x + b, (1, 0)),
        "a*x^2 + b*x + c": (lambda x, a, b, c: a * x ** 2 + b * x + c, (0, 1, 0)),
        "a*x^3 + b*x^2 + c*x + d": (lambda x, a, b, c, d: a * x ** 3 + b * x ** 2 + c * x + d, (0, 0, 1, 0)),
        "a*tanh(b*x + c) + d": (lambda x, a, b, c, d: a * np.tanh(b * x + c) + d, (1, 1, 0, 0)),
        "a*sin(b*x + c) + d": (lambda x, a, b, c, d: a * np.sin(b * x + c) + d, (1, 1, 0, 0)),
        "a*exp(b*x) + c": (lambda x, a, b, c: a * np.exp(b * x) + c, (1, .1, 0)),
    }
    filas = []
    for nombre, (fn, p0) in candidatos.items():
        try:
            p, _ = curve_fit(fn, x, y, p0=p0, maxfev=20000)
            r2 = 1 - ((y - fn(x, *p)) ** 2).sum() / ((y - y.mean()) ** 2).sum()
            filas.append({"forma": nombre, "parametros": np.round(p, 5).tolist(), "R2": float(r2), "n_param": len(p)})
        except (RuntimeError, ValueError):
            continue
    return pd.DataFrame(filas).sort_values(["R2"], ascending=False)


# ---------------------------------------------------------------------------------------------
# Programa principal
# ---------------------------------------------------------------------------------------------
def main():
    t0 = time.time()
    resumen = {}
    datos = cargar_datos()
    P = datos["pares"]
    tr = (P.particion == "entrenamiento").to_numpy()
    resumen["pares"] = P.groupby("particion").agg(pares=("YEAR", "size"), estaciones=("STATION_ID", "nunique")).to_dict("index")
    resumen["estaciones_reservadas"] = datos["reservadas"] or "no entregadas por el docente"
    # Sin fuga: las estaciones/años de validación y prueba no entran a ningún ajuste
    assert P.loc[tr, "YEAR"].max() <= ENTREN_FIN

    # Covariables: altitud estandarizada con entrenamiento; ONI anual tal cual
    alt_m, alt_s = P.loc[tr, "altitud_m"].mean(), P.loc[tr, "altitud_m"].std()
    E_all = np.column_stack([(P.altitud_m - alt_m) / alt_s, P.ONI_anual])
    resumen["estandarizacion_altitud"] = {"media_m": float(alt_m), "de_m": float(alt_s), "fuente": "solo entrenamiento"}

    # FPCA de entrenamiento y elección de M (4–6)
    mu, D, lam = fpca(datos["c_est_raiz"][tr], datos["W_est"])
    var_acum = np.cumsum(lam) / lam.sum()
    M = next((m for m in (4, 5, 6) if var_acum[m - 1] >= .95), 6)
    D_M = D[:, :M]
    pd.DataFrame({"componente": np.arange(1, len(lam) + 1), "autovalor": lam, "varianza_explicada": lam / lam.sum(),
                  "varianza_acumulada": var_acum}).to_csv(TAB / "fpca_varianza_explicada.csv", index=False)
    assert np.allclose(D.T @ datos["W_est"] @ D, np.eye(D.shape[1]), atol=1e-10), "ψ no ortonormales en L²"
    resumen["FPCA"] = {"M": M, "varianza_acumulada_M": float(var_acum[M - 1]),
                       "criterio": "menor M en [4,6] con >= 95 % de varianza (entrenamiento)"}

    tens = {}
    for part in ("entrenamiento", "validacion", "prueba", "reservada"):
        m = (P.particion == part).to_numpy()
        if m.any():
            tens[part] = (torch.as_tensor(datos["c_sat_raiz"][m]), torch.as_tensor(E_all[m]),
                          torch.as_tensor(datos["c_est_raiz"][m]))

    # Parámetros igualados KAN / MLP (±10 %)
    n_in = N_NEURONAS + 2
    cfg = {"kan_ancho": [4], "G": 5}
    n_kan = n_param(KAN([n_in, 4, M], G=5))
    oculta = min(range(2, 200), key=lambda h: abs(n_param(mlp(n_in, h, M)) - n_kan))
    cfg["oculta_mlp"] = oculta
    n_mlp = n_param(mlp(n_in, oculta, M))
    n_lin = n_param(torch.nn.Linear(n_in, M))
    n_enc = N_NEURONAS * datos["K_sat"]
    tabla_param = pd.DataFrame([
        {"cabeza": "lineal", "arquitectura": f"Linear({n_in}->{M})", "param_cabeza": n_lin, "param_total": n_lin + n_enc},
        {"cabeza": "kan", "arquitectura": f"KAN[{n_in},4,{M}], G=5, cúbica", "param_cabeza": n_kan, "param_total": n_kan + n_enc},
        {"cabeza": "mlp", "arquitectura": f"MLP[{n_in},{oculta},{M}], SiLU", "param_cabeza": n_mlp, "param_total": n_mlp + n_enc},
    ])
    dif_rel = abs(n_kan - n_mlp) / n_kan
    tabla_param["dif_rel_kan_mlp"] = dif_rel
    tabla_param.to_csv(TAB / "parametros_cabezas.csv", index=False)
    assert dif_rel <= .10, f"KAN y MLP difieren {dif_rel:.1%} en parámetros"

    # Verificación sintética de la KAN antes de usar datos reales
    sint = verificar_kan_sintetica()
    resumen["KAN_sintetica"] = sint
    assert sint["cumple"], f"KAN sintética no alcanza RMSE ≤ 1e-2: {sint['RMSE_prueba']}"
    print("KAN sintética:", sint, flush=True)

    # Elección de λ con validación (cabeza lineal, primera semilla)
    sel = []
    for lam_r in LAMBDAS:
        _, info, _ = entrenar("lineal", datos, tens, M, mu, D_M, lam_r, SEMILLAS[0], cfg)
        sel.append({"lambda": lam_r, **info})
    sel = pd.DataFrame(sel)
    sel.to_csv(TAB / "seleccion_lambda_validacion.csv", index=False)
    LAM = float(sel.loc[sel.L2_validacion.idxmin(), "lambda"])
    resumen["lambda_rugosidad"] = LAM
    print("λ elegido con validación:", LAM, flush=True)

    # Ablación: 3 cabezas × 3 semillas
    modelos, info_ent, historiales = {}, [], {}
    for tipo in ("lineal", "kan", "mlp"):
        for s in SEMILLAS:
            red, info, hist = entrenar(tipo, datos, tens, M, mu, D_M, LAM, s, cfg)
            modelos[(tipo, s)] = red
            historiales[(tipo, s)] = hist
            info_ent.append({"cabeza": tipo, "semilla": s, "param_total": n_param(red), **info})
            print(tipo, s, info, flush=True)
    pd.DataFrame(info_ent).to_csv(TAB / "entrenamiento_por_semilla.csv", index=False)

    # Equivalencia de la cabeza lineal con la regresión funcional
    X_te, E_te, _ = tens["prueba"]
    eq = []
    for s in SEMILLAS:
        err, beta, tq = verificar_equivalencia_lineal(modelos[("lineal", s)], datos, X_te, E_te, D_M, mu)
        eq.append({"semilla": s, "max_abs_diferencia_sqrt_mm": err, "umbral": 1e-4, "cumple": err <= 1e-4})
        if s == SEMILLAS[0]:
            beta_0, tq_0 = beta, tq
    pd.DataFrame(eq).to_csv(TAB / "equivalencia_lineal_regresion_funcional.csv", index=False)
    assert all(e["cumple"] for e in eq), "La cabeza lineal no reproduce la regresión funcional"
    resumen["equivalencia_lineal"] = eq

    # Evaluación: validación y prueba; baseline = curva satelital sin corregir
    filas, pred_largas, recortes = [], [], []
    for part in ("validacion", "prueba", "reservada"):
        if part not in tens:
            continue
        m = (P.particion == part).to_numpy()
        base = P.loc[m, ["STATION_ID", "YEAR", "fase"]].reset_index(drop=True)
        y_obs = np.clip(datos["c_est_mm"][m] @ datos["Phi_est"].T, 0, None)
        y_sat = np.clip(datos["c_sat_mm"][m] @ datos["Phi_sat"].T, 0, None)
        filas.append(base.assign(particion=part, modelo="satelite_sin_corregir", semilla=0, **metricas(y_sat, y_obs)))
        X, E, _ = tens[part]
        for (tipo, s), red in modelos.items():
            with torch.no_grad():
                c_hat, _ = red(X, E)
            y_hat, neg = a_mm(c_hat.numpy(), datos["Phi_est"])
            recortes.append({"particion": part, "cabeza": tipo, "semilla": s, "pct_puntos_recortados": 100 * neg.mean(),
                             "pct_curvas_con_recorte": 100 * neg.any(1).mean()})
            filas.append(base.assign(particion=part, modelo=tipo, semilla=s, **metricas(y_hat, y_obs)))
            if part != "validacion":
                pred_largas.append((part, tipo, s, base, y_hat))
    met = pd.concat(filas, ignore_index=True)
    met.to_csv(TAB / "metricas_por_curva.csv", index=False)
    pd.DataFrame(recortes).to_csv(TAB / "recorte_no_negatividad.csv", index=False)

    cols = ["L2_rel", "KGE", "sesgo_pct", "error_pico_dias"]
    # media ± DE entre semillas (media de curvas por semilla)
    por_semilla = met.groupby(["particion", "modelo", "semilla"])[cols].mean().reset_index()
    tabla_sem = por_semilla.groupby(["particion", "modelo"])[cols].agg(["mean", "std"])
    tabla_sem.columns = [f"{a}_{b}" for a, b in tabla_sem.columns]
    tabla_sem.reset_index().to_csv(TAB / "ablacion_media_de_semillas.csv", index=False)

    # IC bootstrap por estación (métricas promediadas sobre semillas) y comparaciones pareadas
    met_s = met.groupby(["particion", "modelo", "STATION_ID", "YEAR", "fase"])[cols].mean().reset_index()
    ic, comp = [], []
    for part in met_s.particion.unique():
        for fase in ["todas", "Niña", "Neutro", "Niño"]:
            d = met_s[met_s.particion == part]
            if fase != "todas":
                d = d[d.fase == fase]
            if d.empty:
                continue
            por_est = {mod: media_por_estacion(g, cols) for mod, g in d.groupby("modelo")}
            for mod, g in por_est.items():
                for c in cols:
                    media, lo, hi = bootstrap_estaciones(g[c])
                    ic.append({"particion": part, "fase_ENSO": fase, "modelo": mod, "metrica": c, "media": media,
                               "IC95_inf": lo, "IC95_sup": hi, "estaciones": len(g), "curvas": int((d.modelo == mod).sum())})
            for a, b in [("kan", "lineal"), ("kan", "mlp"), ("mlp", "lineal"), ("kan", "satelite_sin_corregir"),
                         ("lineal", "satelite_sin_corregir"), ("mlp", "satelite_sin_corregir")]:
                ga, gb = por_est[a], por_est[b]
                comunes = ga.index.intersection(gb.index)
                for c in cols:
                    # error absoluto: para sesgo se compara |sesgo|; para KGE, mayor es mejor
                    va, vb = ga.loc[comunes, c], gb.loc[comunes, c]
                    if c == "sesgo_pct":
                        va, vb = va.abs(), vb.abs()
                    diff = (va - vb).to_numpy()
                    media, lo, hi = bootstrap_estaciones(diff)
                    try:
                        w = stats.wilcoxon(va, vb)
                        p_w = float(w.pvalue)
                    except ValueError:
                        p_w = np.nan
                    mejor = "KGE" == c
                    gana = (lo > 0) if mejor else (hi < 0)
                    pierde = (hi < 0) if mejor else (lo > 0)
                    comp.append({"particion": part, "fase_ENSO": fase, "A": a, "B": b,
                                 "metrica": c if c != "sesgo_pct" else "|sesgo_pct|", "media_A_menos_B": media,
                                 "IC95_inf": lo, "IC95_sup": hi, "wilcoxon_p": p_w, "estaciones": len(comunes),
                                 "conclusion": "A mejor" if gana else ("B mejor" if pierde else "sin diferencia (IC incluye 0)")})
    pd.DataFrame(ic).to_csv(TAB / "metricas_IC_bootstrap.csv", index=False)
    pd.DataFrame(comp).to_csv(TAB / "comparaciones_pareadas_bootstrap_wilcoxon.csv", index=False)

    # Predicciones en el formato exigido (station_id, year, day, value_mm): cabeza elegida por validación
    val = por_semilla[por_semilla.particion == "validacion"].groupby("modelo").L2_rel.mean()
    elegida = val.drop("satelite_sin_corregir").idxmin()
    resumen["cabeza_elegida_por_validacion"] = elegida
    for part in ("prueba", "reservada"):
        sel_p = [(b, y) for (p, t, s, b, y) in pred_largas if p == part and t == elegida]
        if not sel_p:
            continue
        base = sel_p[0][0]
        y_ens = np.mean([y for _, y in sel_p], axis=0)        # promedio de las 3 semillas
        largo = pd.DataFrame({"station_id": np.repeat(base.STATION_ID.to_numpy(), 365),
                              "year": np.repeat(base.YEAR.to_numpy(), 365),
                              "day": np.tile(np.arange(1, 366), len(base)),
                              "value_mm": y_ens.ravel()})
        largo.to_csv(PRED / f"predicciones_{part}_{elegida}.csv", index=False, float_format="%.4f")

    # -----------------------------------------------------------------------------------------
    # Interpretación
    # -----------------------------------------------------------------------------------------
    s0 = SEMILLAS[0]
    tq_d = DIAS
    Phi_sat = datos["Phi_sat"]
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.6), sharey=False)
    for ax, tipo in zip(axs, ("lineal", "kan")):
        Wk = modelos[(tipo, s0)].A.detach().numpy() @ Phi_sat.T
        for k in range(N_NEURONAS):
            ax.plot(tq_d, Wk[k], lw=1, label=f"w{k + 1}")
        ax.axhline(0, c="k", lw=.5)
        for a, b in VENTANAS_PICO:
            ax.axvspan(a, b, color="0.9", zorder=0)
        ax.set(title=f"Funciones de peso w_k(t) — cabeza {tipo}", xlabel="día del año")
    axs[1].legend(ncol=4, fontsize=7)
    plt.tight_layout()
    fig.savefig(FIG / "pesos_wk.png", dpi=130)
    plt.close(fig)
    pd.DataFrame(modelos[("lineal", s0)].A.detach().numpy() @ Phi_sat.T, index=[f"w{k+1}" for k in range(N_NEURONAS)]) \
        .to_csv(TAB / "pesos_wk_lineal.csv")
    pd.DataFrame(modelos[("kan", s0)].A.detach().numpy() @ Phi_sat.T, index=[f"w{k+1}" for k in range(N_NEURONAS)]) \
        .to_csv(TAB / "pesos_wk_kan.csv")

    # Activaciones de la primera capa KAN en el rango de los datos de entrenamiento
    red_k = modelos[("kan", s0)]
    with torch.no_grad():
        H_tr = red_k.entrada_cabeza(*tens["entrenamiento"][:2])
    capa = red_k.cabeza.capas[0]
    nombres = [f"z{k + 1}" for k in range(N_NEURONAS)] + ["altitud", "ONI"]
    lo, hi = H_tr.min(0).values, H_tr.max(0).values
    rejilla = torch.stack([torch.linspace(float(a), float(b), 200) for a, b in zip(lo, hi)], 1)
    with torch.no_grad():
        act = capa.activaciones(rejilla).numpy()          # (200, n_in, n_out)
    n_out = act.shape[2]
    fig, axs = plt.subplots(n_out, n_in, figsize=(1.6 * n_in, 1.4 * n_out), sharex="col")
    for j in range(n_out):
        for i in range(n_in):
            ax = axs[j, i]
            ax.plot(rejilla[:, i], act[:, i, j], lw=1)
            ax.tick_params(labelsize=5)
            if j == 0:
                ax.set_title(nombres[i], fontsize=7)
    fig.suptitle("Activaciones φ_ij de la primera capa KAN (rango observado en entrenamiento)", fontsize=9)
    plt.tight_layout()
    fig.savefig(FIG / "activaciones_kan_capa1.png", dpi=130)
    plt.close(fig)

    # Contraste con la cabeza lineal, invariante al codificador. Las z_k de redes entrenadas por separado no son
    # las mismas funciones (y con mucha rugosidad las w_k son casi colineales), así que no se comparan pesos
    # neurona a neurona. Se compara el efecto funcional: para cualquier cabeza,
    #   β_ef(s,t) = Σ_k w_k(s) · E_entren[∂ŷ(t)/∂z_k] / |T|      (derivada funcional media de ŷ(t) respecto de X(s))
    # y para altitud/ONI la curva E[∂ŷ(t)/∂e]. En la cabeza lineal la derivada es constante y β_ef es exactamente
    # su β(s,t); en la KAN/MLP es la mejor aproximación lineal media. La no linealidad se mide con la variación
    # del jacobiano entre curvas de entrenamiento: CV = DE/|media| de ∂ŷ/∂entrada.
    amplitud = act.max(0) - act.min(0)                     # (n_in, n_out) de la primera capa
    psi_t = (D_M.T @ datos["Phi_est"].T)                   # (M, 365)
    tq_s, _ = cuadratura(datos["K_sat"])
    Phi_q = base_bspline(tq_s, (1, 365), datos["K_sat"])

    def efecto_funcional(red):
        X, E = tens["entrenamiento"][0], tens["entrenamiento"][1]
        z = red.codificar(X).detach().requires_grad_(True)
        e = E.clone().requires_grad_(True)
        s = red.cabeza(torch.cat([red.norm(z), e], 1))
        J = []
        for m in range(M):
            gz, ge = torch.autograd.grad(s[:, m].sum(), (z, e), retain_graph=True)
            J.append(torch.cat([gz, ge], 1).numpy())
        J = np.stack(J, 1)                                  # (n, M, 8+2): ∂s_m/∂entrada (z en escala cruda)
        Jm = J.mean(0)
        w_s = red.A.detach().numpy() @ Phi_q.T              # (8, nq)
        beta_ef = np.einsum("ks,mk,mt->st", w_s, Jm[:, :N_NEURONAS], psi_t) / T_LEN
        g_e = Jm[:, N_NEURONAS:].T @ psi_t                  # (2, 365): efecto de altitud y ONI sobre ŷ(t)
        # no linealidad por entrada: variación del efecto sobre ŷ entre curvas (norma L² en t)
        efecto_por_curva = np.einsum("nmi,mt->nit", J, psi_t)
        cv = np.linalg.norm(efecto_por_curva.std(0), axis=1) / np.maximum(np.linalg.norm(efecto_por_curva.mean(0), axis=1), 1e-12)
        return beta_ef, g_e, cv

    efectos = {tipo: efecto_funcional(modelos[(tipo, s0)]) for tipo in ("lineal", "kan", "mlp")}
    with torch.no_grad():
        _, beta_lin_ref, _ = verificar_equivalencia_lineal(modelos[("lineal", s0)], datos, *tens["prueba"][:2], D_M, mu)
    assert np.abs(efectos["lineal"][0] - beta_lin_ref).max() < 1e-10, "β_ef de la cabeza lineal debe ser su β(s,t)"
    b_lin = efectos["lineal"][0]
    filas_c = []
    for tipo in ("kan", "mlp"):
        b, g, cv = efectos[tipo]
        filas_c.append({"cabeza": tipo,
                        "corr_beta_con_lineal": float(np.corrcoef(b.ravel(), b_lin.ravel())[0, 1]),
                        "dif_rel_beta_L2": float(np.linalg.norm(b - b_lin) / np.linalg.norm(b_lin)),
                        "corr_efecto_altitud_con_lineal": float(np.corrcoef(g[0], efectos["lineal"][1][0])[0, 1]),
                        "corr_efecto_ONI_con_lineal": float(np.corrcoef(g[1], efectos["lineal"][1][1])[0, 1]),
                        "CV_jacobiano_z_medio": float(cv[:N_NEURONAS].mean()),
                        "CV_jacobiano_altitud": float(cv[N_NEURONAS]), "CV_jacobiano_ONI": float(cv[N_NEURONAS + 1])})
    filas_c.append({"cabeza": "lineal", "corr_beta_con_lineal": 1.0, "dif_rel_beta_L2": 0.0,
                    "corr_efecto_altitud_con_lineal": 1.0, "corr_efecto_ONI_con_lineal": 1.0,
                    "CV_jacobiano_z_medio": float(efectos["lineal"][2][:N_NEURONAS].mean()),
                    "CV_jacobiano_altitud": float(efectos["lineal"][2][N_NEURONAS]),
                    "CV_jacobiano_ONI": float(efectos["lineal"][2][N_NEURONAS + 1])})
    contraste = pd.DataFrame(filas_c)
    contraste.to_csv(TAB / "contraste_kan_vs_lineal.csv", index=False)
    resumen["contraste_kan_vs_lineal"] = contraste.to_dict("records")
    pd.DataFrame(modelos[("lineal", s0)].cabeza.weight.detach().numpy().T, index=nombres,
                 columns=[f"s{m + 1}" for m in range(M)]).to_csv(TAB / "cabeza_lineal_coeficientes.csv")
    pd.DataFrame({"dia": DIAS, **{f"efecto_{v}_{tipo}": efectos[tipo][1][i] for tipo in efectos
                                  for i, v in enumerate(("altitud", "ONI"))}}).to_csv(TAB / "efecto_altitud_ONI_por_cabeza.csv", index=False)

    fig, axs = plt.subplots(1, 4, figsize=(15, 3.6))
    vmax = max(np.abs(efectos[t][0]).max() for t in ("lineal", "kan"))
    for ax, tipo in zip(axs[:2], ("lineal", "kan")):
        im = ax.pcolormesh(DIAS, tq_s, efectos[tipo][0], cmap="RdBu_r", shading="auto", vmin=-vmax, vmax=vmax)
        ax.set(title=f"β efectivo(s,t), cabeza {tipo}", xlabel="t (estación)", ylabel="s (satélite)")
        ax.title.set_fontsize(9)
    plt.colorbar(im, ax=axs[1])
    for ax, i, v in ((axs[2], 0, "altitud (+1 DE)"), (axs[3], 1, "ONI (+1 °C)")):
        for tipo in ("lineal", "kan", "mlp"):
            ax.plot(DIAS, efectos[tipo][1][i], label=tipo)
        ax.axhline(0, c="k", lw=.5)
        ax.set(title=f"Efecto medio de {v} sobre √ŷ(t)", xlabel="día del año")
        ax.title.set_fontsize(9)
    axs[3].legend(fontsize=7)
    plt.tight_layout()
    fig.savefig(FIG / "contraste_kan_lineal.png", dpi=130)
    plt.close(fig)

    # Ajuste simbólico de la activación de mayor amplitud (y de la más no lineal)
    i_max, j_max = np.unravel_index(np.argmax(amplitud), amplitud.shape)
    no_lin = np.array([[1 - np.corrcoef(rejilla[:, i], act[:, i, j])[0, 1] ** 2 for j in range(n_out)] for i in range(n_in)])
    no_lin[amplitud < .25 * amplitud.max()] = -1
    i_nl, j_nl = np.unravel_index(np.argmax(no_lin), no_lin.shape)
    simb = []
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.2))
    for ax, (i, j, etiqueta) in zip(axs, [(i_max, j_max, "mayor amplitud"), (i_nl, j_nl, "más no lineal")]):
        x, y = rejilla[:, i].numpy(), act[:, i, j]
        aj = ajuste_simbolico(x, y)
        mejor = aj.iloc[0]
        simb.append({"activacion": f"φ({nombres[i]} -> h{j + 1})", "criterio": etiqueta, **mejor.to_dict()})
        aj.assign(activacion=f"φ({nombres[i]} -> h{j + 1})").to_csv(TAB / f"ajuste_simbolico_{nombres[i]}_h{j + 1}.csv", index=False)
        ax.plot(x, y, lw=2, label="activación KAN")
        ax.set(title=f"φ({nombres[i]}→h{j + 1}), {etiqueta}\nmejor: {mejor.forma}, R²={mejor.R2:.4f}", xlabel=nombres[i])
        ax.title.set_fontsize(8)
    plt.tight_layout()
    fig.savefig(FIG / "ajuste_simbolico.png", dpi=130)
    plt.close(fig)
    pd.DataFrame(simb).to_csv(TAB / "ajuste_simbolico_resumen.csv", index=False)
    resumen["ajuste_simbolico"] = simb

    # Superficie β(s,t) de la regresión funcional equivalente (cabeza lineal)
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    im = ax.pcolormesh(DIAS, tq_0, beta_0, cmap="RdBu_r", shading="auto",
                       vmin=-np.abs(beta_0).max(), vmax=np.abs(beta_0).max())
    ax.set(xlabel="t: día de la curva de estación", ylabel="s: día de la curva satelital",
           title="β(s,t) de la cabeza lineal (regresión funcional equivalente)")
    ax.title.set_fontsize(8)
    plt.colorbar(im, ax=ax)
    plt.tight_layout()
    fig.savefig(FIG / "beta_regresion_funcional.png", dpi=130)
    plt.close(fig)

    # FPCA: media y componentes
    fig, ax = plt.subplots(figsize=(7, 3.2))
    psi = D_M.T @ datos["Phi_est"].T
    for m in range(M):
        ax.plot(DIAS, psi[m], label=f"ψ{m + 1} ({lam[m] / lam.sum():.0%})")
    ax.plot(DIAS, mu @ datos["Phi_est"].T / 10, "k--", label="μ(t)/10")
    ax.legend(fontsize=7, ncol=3)
    ax.set(title="FPCA de entrenamiento (escala √P)", xlabel="día del año")
    plt.tight_layout()
    fig.savefig(FIG / "fpca_componentes.png", dpi=130)
    plt.close(fig)

    # Peores casos de prueba (cabeza elegida)
    mp = met_s[(met_s.particion == "prueba") & (met_s.modelo == elegida)].nlargest(4, "L2_rel")
    m_te = (P.particion == "prueba").to_numpy()
    base_te = P.loc[m_te, ["STATION_ID", "YEAR"]].reset_index(drop=True)
    y_obs_te = np.clip(datos["c_est_mm"][m_te] @ datos["Phi_est"].T, 0, None)
    y_sat_te = np.clip(datos["c_sat_mm"][m_te] @ datos["Phi_sat"].T, 0, None)
    y_el = np.mean([y for (p, t, s, b, y) in pred_largas if p == "prueba" and t == elegida], axis=0)
    fig, axs = plt.subplots(1, len(mp), figsize=(3.2 * len(mp), 3), sharey=False)
    for ax, (_, r) in zip(np.atleast_1d(axs), mp.iterrows()):
        i = base_te.index[(base_te.STATION_ID == r.STATION_ID) & (base_te.YEAR == r.YEAR)][0]
        ax.plot(DIAS, y_obs_te[i], "k", lw=1.5, label="estación")
        ax.plot(DIAS, y_sat_te[i], c="tab:orange", lw=1, label="satélite")
        ax.plot(DIAS, y_el[i], c="tab:blue", lw=1, label=elegida)
        ax.set(title=f"{r.STATION_ID} {r.YEAR} (L²rel={r.L2_rel:.2f})", xlabel="día")
        ax.title.set_fontsize(8)
    np.atleast_1d(axs)[0].legend(fontsize=7)
    plt.tight_layout()
    fig.savefig(FIG / "peores_casos_prueba.png", dpi=130)
    plt.close(fig)
    resumen["peores_casos_prueba"] = mp[["STATION_ID", "YEAR", "fase", "L2_rel", "KGE", "sesgo_pct"]].to_dict("records")

    # Curvas de aprendizaje
    fig, ax = plt.subplots(figsize=(6, 3))
    for (tipo, s), h in historiales.items():
        if s == s0:
            h = np.array(h)
            ax.plot(h[:, 0], h[:, 2], label=f"{tipo} (validación)")
    ax.set(yscale="log", xlabel="época", ylabel="L² funcional", title="Pérdida de validación (semilla 1)")
    ax.legend(fontsize=7)
    plt.tight_layout()
    fig.savefig(FIG / "curvas_aprendizaje.png", dpi=130)
    plt.close(fig)

    resumen["versiones"] = {"python": __import__("sys").version.split()[0], "numpy": np.__version__, "pandas": pd.__version__,
                            "torch": torch.__version__, "scipy": __import__("scipy").__version__}
    resumen["semillas"] = list(SEMILLAS)
    resumen["segundos"] = round(time.time() - t0, 1)
    (OUT / "resumen_fase4.json").write_text(json.dumps(resumen, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    print(json.dumps({k: resumen[k] for k in ("pares", "FPCA", "lambda_rugosidad", "cabeza_elegida_por_validacion", "segundos")},
                     indent=2, ensure_ascii=False, default=str))


if __name__ == "__main__":
    main()

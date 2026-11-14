"""Red funcional con cabeza KAN sobre las salidas verificadas de Jacobo.

El script no lee resultados de notebooks KAN anteriores. Parte de:
  - eda_fase0/pares_estacion_anio_validos_preliminar.csv
  - fda_fase1/coeficientes_*.npz y matrices_gram.npz
  - Datos/CVC_Wheather_Stations_Long_Lat.* y CatalogoEstacionesElevacion.xlsx
  - datos_web/ONI_NOAA_CPC_1981_2013.csv

Entrenamiento: 1981-1995. Validación: 1996-2000. La prueba 2001-2003
permanece cerrada y no se evalúa en este script.
"""

from __future__ import annotations

import copy
import json
import math
import random
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from scipy.stats import wilcoxon
from torch.utils.data import DataLoader, TensorDataset


KAN_ROOT = Path(__file__).resolve().parent
ROOT = KAN_ROOT.parent
E0 = ROOT / "eda_fase0"
F1 = ROOT / "fda_fase1"
DATA = ROOT / "Datos"
WEB = ROOT / "datos_web"
OUT = KAN_ROOT / "resultados"
OUT.mkdir(exist_ok=True)

DEVICE = torch.device("cpu")
SEEDS = [2026, 2027, 2028]
N_NEURONAS_FUNCIONALES = 8
N_COMPONENTES_FPCA = 6
KAN_BASES = 8
MLP_HIDDEN = 31
MAX_EPOCHS = 220
PATIENCE = 30
BATCH_SIZE = 64
LEARNING_RATE = 2e-3
LAMBDA_RUGOSIDAD = 1e-4


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.use_deterministic_algorithms(True)


def normalizar_codigo(serie: pd.Series) -> pd.Series:
    return pd.to_numeric(serie, errors="coerce").round().astype("Int64").astype("string")


def particion(anio: pd.Series) -> pd.Series:
    return pd.Series(
        np.select(
            [anio.between(1981, 1995), anio.between(1996, 2000)],
            ["entrenamiento", "validacion"], default="prueba",
        ),
        index=anio.index,
    )


def pixel_mas_cercano(lat: np.ndarray, lon: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    fila = np.rint(((5.0 - 0.025) - lat) / 0.05).astype(int)
    col = np.rint((lon - (-77.55 + 0.025)) / 0.05).astype(int)
    return np.clip(fila, 0, 38), np.clip(col, 0, 36)


def base_bspline(t: np.ndarray, nodos: np.ndarray, orden: int) -> np.ndarray:
    """Cox-de Boor, igual a la construcción paso a paso de Jacobo."""
    t = np.asarray(t, float)
    B = np.zeros((len(t), len(nodos) - 1))
    for j in range(len(nodos) - 1):
        B[:, j] = (t >= nodos[j]) & (t < nodos[j + 1])
    B[t == nodos[-1], np.max(np.flatnonzero(nodos[:-1] < nodos[1:]))] = 1
    for k in range(2, orden + 1):
        Bn = np.zeros((len(t), len(nodos) - k))
        for j in range(len(nodos) - k):
            d1 = nodos[j + k - 1] - nodos[j]
            d2 = nodos[j + k] - nodos[j + 1]
            if d1 > 0:
                Bn[:, j] += (t - nodos[j]) / d1 * B[:, j]
            if d2 > 0:
                Bn[:, j] += (nodos[j + k] - t) / d2 * B[:, j + 1]
        B = Bn
    return B


def cargar_datos() -> dict:
    pares = pd.read_csv(
        E0 / "pares_estacion_anio_validos_preliminar.csv",
        dtype={"STATION_ID": "string"},
    ).rename(columns={"STATION_ID": "station_id", "YEAR": "year"})
    pares["year"] = pares["year"].astype(int)
    pares["particion"] = particion(pares["year"])

    coords = pd.read_csv(DATA / "CVC_Wheather_Stations_Long_Lat.csv", encoding="latin1")
    coords["station_id"] = normalizar_codigo(coords["CODIGO"])
    coords = coords.rename(columns={"STATION_NAME": "station_name", "Lon": "longitude", "lat": "latitude"})
    coords = coords[["station_id", "station_name", "longitude", "latitude"]].drop_duplicates("station_id")

    # El DTA se usa únicamente como verificación del CSV suministrado.
    coords_dta = pd.read_stata(DATA / "CVC_Wheather_Stations_Long_Lat.dta")
    coords_dta["station_id"] = normalizar_codigo(coords_dta["station_id"])
    assert set(coords.station_id) == set(coords_dta.station_id)

    elev = pd.read_excel(DATA / "CatalogoEstacionesElevacion.xlsx")
    elev["station_id"] = normalizar_codigo(elev["codigo"])
    elev = elev.rename(columns={"alt": "elevation_m"})[["station_id", "elevation_m"]]
    elev = elev.drop_duplicates("station_id")

    base = pares.merge(coords, on="station_id", how="left", validate="many_to_one")
    base = base.merge(elev, on="station_id", how="left", validate="many_to_one")
    tiene_coordenadas = base.latitude.notna() & base.longitude.notna()
    tiene_altitud = base.elevation_m.notna()
    dentro = base.latitude.between(3.05, 5.0) & base.longitude.between(-77.55, -75.70)
    trazabilidad = []
    for etapa, mascara in [
        ("pares_validos", np.ones(len(base), dtype=bool)),
        ("coordenadas_disponibles", tiene_coordenadas),
        ("coordenadas_y_altitud", tiene_coordenadas & tiene_altitud),
        ("coordenadas_altitud_y_dominio", tiene_coordenadas & tiene_altitud & dentro),
    ]:
        subset = base.loc[mascara]
        trazabilidad.append({
            "etapa": etapa,
            "pares_estacion_anio": len(subset),
            "estaciones": subset.station_id.nunique(),
        })
    base["motivo_exclusion"] = np.select(
        [
            ~tiene_coordenadas,
            tiene_coordenadas & ~tiene_altitud,
            tiene_coordenadas & tiene_altitud & ~dentro,
        ],
        ["sin_coordenadas", "sin_altitud", "fuera_del_dominio"],
        default="incluida",
    )
    estaciones_excluidas = (
        base.loc[base.motivo_exclusion.ne("incluida"), ["station_id", "motivo_exclusion"]]
        .drop_duplicates()
        .sort_values("station_id")
        .reset_index(drop=True)
    )
    base = base.loc[base.motivo_exclusion.eq("incluida")].copy()
    base["pixel_row"], base["pixel_col"] = pixel_mas_cercano(
        base.latitude.to_numpy(), base.longitude.to_numpy()
    )
    base["pixel_latitude"] = 4.975 - 0.05 * base.pixel_row
    base["pixel_longitude"] = -77.525 + 0.05 * base.pixel_col
    lat1 = np.radians(base.latitude.to_numpy())
    lat2 = np.radians(base.pixel_latitude.to_numpy())
    dlat = lat2 - lat1
    dlon = np.radians(base.pixel_longitude.to_numpy() - base.longitude.to_numpy())
    haversine = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    base["distancia_pixel_km"] = 6371.0 * 2 * np.arcsin(np.sqrt(haversine))

    est = np.load(F1 / "coeficientes_estaciones.npz")
    sat = np.load(F1 / "coeficientes_satelite.npz")
    trans = np.load(F1 / "coeficientes_transformados.npz")
    gram = np.load(F1 / "matrices_gram.npz")
    spec = json.loads((F1 / "base_especificacion.json").read_text(encoding="utf-8"))

    idx_est = pd.DataFrame({
        "station_id": pd.Series(est["station_id"]).astype("int64").astype("string"),
        "year": est["anio"].astype(int),
        "idx_est": np.arange(len(est["anio"])),
    })
    idx_sat = pd.DataFrame({
        "year": sat["anio"].astype(int), "pixel_row": sat["fila"].astype(int),
        "pixel_col": sat["col"].astype(int), "idx_sat": np.arange(len(sat["anio"])),
    })
    base = base.merge(idx_est, on=["station_id", "year"], validate="one_to_one")
    base = base.merge(idx_sat, on=["year", "pixel_row", "pixel_col"], validate="many_to_one")
    trazabilidad.append({
        "etapa": "curvas_estacion_y_pixel_emparejadas",
        "pares_estacion_anio": len(base),
        "estaciones": base.station_id.nunique(),
    })

    oni = pd.read_csv(WEB / "ONI_NOAA_CPC_1981_2013.csv")
    columnas_oni = [c for c in oni.columns if c != "anio"]
    oni["oni_anual"] = oni[columnas_oni].mean(axis=1)
    oni["fase_enso"] = np.select(
        [oni.oni_anual.ge(.5), oni.oni_anual.le(-.5)], ["Nino", "Nina"], default="Neutral"
    )
    base = base.merge(oni[["anio", "oni_anual", "fase_enso"]], left_on="year", right_on="anio", validate="many_to_one")

    iest, isat = base.idx_est.to_numpy(), base.idx_sat.to_numpy()
    t = np.arange(1, 366, dtype=float)
    phi_est = base_bspline(t, np.asarray(spec["nodos_estaciones"]), int(spec["orden"]))
    phi_sat = base_bspline(t, np.asarray(spec["nodos_satelite"]), int(spec["orden"]))

    resultado = {
        "ids": base.reset_index(drop=True),
        "trazabilidad_emparejamiento": pd.DataFrame(trazabilidad),
        "estaciones_excluidas": estaciones_excluidas,
        "mapa_estacion_pixel": (
            base[[
                "station_id", "station_name", "latitude", "longitude", "elevation_m",
                "pixel_row", "pixel_col", "pixel_latitude", "pixel_longitude",
                "distancia_pixel_km",
            ]]
            .drop_duplicates("station_id")
            .sort_values("station_id")
            .reset_index(drop=True)
        ),
        "x_sqrt": trans["raiz_sat"][isat].astype("float32"),
        "y_sqrt": trans["raiz_est"][iest].astype("float32"),
        "x_mm_coef": sat["coef"][isat].astype("float32"),
        "y_mm_coef": est["coef"][iest].astype("float32"),
        "covariables": base[["elevation_m", "oni_anual"]].to_numpy("float32"),
        "W_sat": gram["W_sat"].astype("float32"),
        "W_est": gram["W_est"].astype("float32"),
        "phi_sat": phi_sat.astype("float32"),
        "phi_est": phi_est.astype("float32"),
        "spec": spec,
    }
    assert len(base) == 2051
    assert base.loc[base.particion.eq("entrenamiento")].shape[0] == 1349
    assert base.loc[base.particion.eq("validacion")].shape[0] == 445
    return resultado


def ajustar_fpca(y_train: np.ndarray, W: np.ndarray, m: int = 6) -> dict:
    media = y_train.mean(axis=0)
    L = np.linalg.cholesky(W)
    z = (y_train - media) @ L
    _, s, vt = np.linalg.svd(z, full_matrices=False)
    var = s**2 / (len(z) - 1)
    proporcion = var / var.sum()
    direcciones = vt[:m] @ np.linalg.inv(L)
    return {
        "media": media.astype("float32"),
        "direcciones": direcciones.astype("float32"),
        "proporcion": proporcion.astype(float),
    }


class CodificadorFuncional(nn.Module):
    W: torch.Tensor
    phi: torch.Tensor

    def __init__(self, W_sat: np.ndarray, phi_sat: np.ndarray, n_neuronas: int = 8):
        super().__init__()
        self.register_buffer("W", torch.tensor(W_sat, dtype=torch.float32))
        self.register_buffer("phi", torch.tensor(phi_sat, dtype=torch.float32))
        self.A = nn.Parameter(torch.randn(n_neuronas, W_sat.shape[0]) * 0.08)

    def forward(self, coef: torch.Tensor) -> torch.Tensor:
        return coef @ self.W @ self.A.T

    def rugosidad(self) -> torch.Tensor:
        pesos = self.phi @ self.A.T
        segunda = pesos[2:] - 2 * pesos[1:-1] + pesos[:-2]
        return segunda.square().mean()


class CapaKANBSpline(nn.Module):
    """Capa KAN: una función B-spline aprendible en cada arista."""
    def __init__(self, n_in: int, n_out: int, n_bases: int = 8, grado: int = 3, rango=(-3.0, 3.0)):
        super().__init__()
        self.n_in, self.n_out, self.n_bases, self.grado = n_in, n_out, n_bases, grado
        internos = n_bases - grado - 1
        puntos = torch.linspace(rango[0], rango[1], internos + 2)[1:-1] if internos > 0 else torch.empty(0)
        nodos = torch.cat([torch.full((grado + 1,), rango[0]), puntos, torch.full((grado + 1,), rango[1])])
        self.register_buffer("nodos", nodos)
        self.base_weight = nn.Parameter(torch.empty(n_out, n_in))
        self.spline_weight = nn.Parameter(torch.randn(n_out, n_in, n_bases) * 0.05)
        nn.init.kaiming_uniform_(self.base_weight, a=math.sqrt(5))

    @torch.no_grad()
    def ajustar_rango(self, caracteristicas: torch.Tensor, margen: float = 0.05) -> None:
        minimo = float(caracteristicas.min())
        maximo = float(caracteristicas.max())
        amplitud = max(maximo - minimo, 1e-3)
        minimo -= margen * amplitud
        maximo += margen * amplitud
        internos = self.n_bases - self.grado - 1
        puntos = (
            torch.linspace(minimo, maximo, internos + 2, device=self.nodos.device)[1:-1]
            if internos > 0 else torch.empty(0, device=self.nodos.device)
        )
        self.nodos = torch.cat([
            torch.full((self.grado + 1,), minimo, device=self.nodos.device),
            puntos,
            torch.full((self.grado + 1,), maximo, device=self.nodos.device),
        ])

    def basis(self, x: torch.Tensor) -> torch.Tensor:
        x = x.clamp(float(self.nodos[0]), float(self.nodos[-1]) - 1e-6)
        B = ((x[..., None] >= self.nodos[:-1]) & (x[..., None] < self.nodos[1:])).float()
        for p in range(1, self.grado + 1):
            n_cols = len(self.nodos) - p - 1
            izquierda = torch.zeros(*x.shape, n_cols, device=x.device)
            derecha = torch.zeros_like(izquierda)
            d1 = self.nodos[p:p+n_cols] - self.nodos[:n_cols]
            d2 = self.nodos[p+1:p+1+n_cols] - self.nodos[1:1+n_cols]
            m1, m2 = d1 > 0, d2 > 0
            izquierda[..., m1] = ((x[..., None] - self.nodos[:n_cols])[..., m1] / d1[m1]) * B[..., :n_cols][..., m1]
            derecha[..., m2] = ((self.nodos[p+1:p+1+n_cols] - x[..., None])[..., m2] / d2[m2]) * B[..., 1:n_cols+1][..., m2]
            B = izquierda + derecha
        return B

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        base = F.silu(x) @ self.base_weight.T
        spline = torch.einsum("big,oig->bo", self.basis(x), self.spline_weight)
        return base + spline

    @torch.no_grad()
    def edge_function(self, salida: int, entrada: int, rejilla: torch.Tensor) -> torch.Tensor:
        b = self.basis(rejilla[:, None])[:, 0, :]
        return self.base_weight[salida, entrada] * F.silu(rejilla) + b @ self.spline_weight[salida, entrada]


class RedFuncional(nn.Module):
    encoder: CodificadorFuncional
    head: nn.Module
    media_caracteristicas: torch.Tensor
    escala_caracteristicas: torch.Tensor
    media_y: torch.Tensor
    direcciones: torch.Tensor

    def __init__(self, tipo: str, W_sat, phi_sat, media_y, direcciones, mlp_hidden=31):
        super().__init__()
        self.tipo = tipo
        self.encoder = CodificadorFuncional(W_sat, phi_sat, N_NEURONAS_FUNCIONALES)
        n_in = N_NEURONAS_FUNCIONALES + 2
        self.register_buffer("media_caracteristicas", torch.zeros(n_in))
        self.register_buffer("escala_caracteristicas", torch.ones(n_in))
        m = direcciones.shape[0]
        if tipo == "lineal":
            self.head = nn.Linear(n_in, m)
        elif tipo == "mlp":
            self.head = nn.Sequential(nn.Linear(n_in, mlp_hidden), nn.GELU(), nn.Linear(mlp_hidden, m))
        elif tipo == "kan":
            self.head = CapaKANBSpline(n_in, m, KAN_BASES)
        else:
            raise ValueError(tipo)
        self.register_buffer("media_y", torch.tensor(media_y, dtype=torch.float32))
        self.register_buffer("direcciones", torch.tensor(direcciones, dtype=torch.float32))

    def forward(self, x: torch.Tensor, cov: torch.Tensor, devolver_scores=False):
        z = self.encoder(x)
        caracteristicas = (
            torch.cat([z, cov], dim=1) - self.media_caracteristicas
        ) / self.escala_caracteristicas
        scores = self.head(caracteristicas)
        coef = self.media_y + scores @ self.direcciones
        return (coef, scores, caracteristicas) if devolver_scores else coef


def perdida_funcional(pred, real, W, modelo):
    d = pred - real
    l2 = torch.einsum("bi,ij,bj->b", d, W, d).mean()
    return l2 + LAMBDA_RUGOSIDAD * modelo.encoder.rugosidad(), l2


def preparar_tensores(datos: dict) -> dict:
    ids = datos["ids"]
    train = ids.particion.eq("entrenamiento").to_numpy()
    val = ids.particion.eq("validacion").to_numpy()
    test = ids.particion.eq("prueba").to_numpy()
    media_x = datos["x_sqrt"][train].mean(axis=0)
    xc = datos["x_sqrt"][train] - media_x
    escala_x = float(np.sqrt(np.mean(np.einsum("bi,ij,bj->b", xc, datos["W_sat"], xc))))
    x = (datos["x_sqrt"] - media_x) / escala_x
    media_cov = datos["covariables"][train].mean(axis=0)
    sd_cov = datos["covariables"][train].std(axis=0)
    cov = (datos["covariables"] - media_cov) / sd_cov
    return {
        "x": x.astype("float32"), "cov": cov.astype("float32"), "y": datos["y_sqrt"],
        "train": train, "val": val, "test": test,
        "media_x": media_x, "escala_x": escala_x,
        "media_cov": media_cov, "sd_cov": sd_cov,
    }


def loader(tensores, mascara, shuffle, seed):
    ds = TensorDataset(
        torch.tensor(tensores["x"][mascara]),
        torch.tensor(tensores["cov"][mascara]),
        torch.tensor(tensores["y"][mascara]),
    )
    gen = torch.Generator().manual_seed(seed)
    return DataLoader(ds, batch_size=BATCH_SIZE, shuffle=shuffle, generator=gen)


def entrenar(tipo, seed, datos, tensores, fpca):
    set_seed(seed)
    modelo = RedFuncional(tipo, datos["W_sat"], datos["phi_sat"], fpca["media"], fpca["direcciones"]).to(DEVICE)
    x_train = torch.tensor(tensores["x"][tensores["train"]], device=DEVICE)
    cov_train = torch.tensor(tensores["cov"][tensores["train"]], device=DEVICE)

    @torch.no_grad()
    def actualizar_normalizacion():
        caracteristicas_train = torch.cat([modelo.encoder(x_train), cov_train], dim=1)
        media = caracteristicas_train.mean(dim=0)
        escala = caracteristicas_train.std(dim=0, unbiased=False).clamp_min(1e-6)
        modelo.media_caracteristicas = media
        modelo.escala_caracteristicas = escala

    actualizar_normalizacion()
    if tipo == "kan":
        if not isinstance(modelo.head, CapaKANBSpline):
            raise TypeError("La cabeza KAN debe ser una CapaKANBSpline.")
        with torch.no_grad():
            caracteristicas_train = torch.cat([modelo.encoder(x_train), cov_train], dim=1)
            caracteristicas_estandarizadas = (
                caracteristicas_train - modelo.media_caracteristicas
            ) / modelo.escala_caracteristicas
            modelo.head.ajustar_rango(caracteristicas_estandarizadas)
    W = torch.tensor(datos["W_est"], dtype=torch.float32, device=DEVICE)
    opt = torch.optim.Adam(modelo.parameters(), lr=LEARNING_RATE)
    tr = loader(tensores, tensores["train"], True, seed)
    va = loader(tensores, tensores["val"], False, seed)
    mejor, mejor_epoca, espera = float("inf"), 0, 0
    mejor_estado = None
    historia = []
    for epoca in range(1, MAX_EPOCHS + 1):
        actualizar_normalizacion()
        modelo.train(); suma, n = 0.0, 0
        for xb, cb, yb in tr:
            xb, cb, yb = xb.to(DEVICE), cb.to(DEVICE), yb.to(DEVICE)
            opt.zero_grad()
            pred = modelo(xb, cb)
            loss, l2 = perdida_funcional(pred, yb, W, modelo)
            loss.backward(); opt.step()
            suma += float(l2.detach()) * len(xb); n += len(xb)
        actualizar_normalizacion()
        modelo.eval(); vsuma, vn = 0.0, 0
        with torch.no_grad():
            for xb, cb, yb in va:
                xb, cb, yb = xb.to(DEVICE), cb.to(DEVICE), yb.to(DEVICE)
                _, l2 = perdida_funcional(modelo(xb, cb), yb, W, modelo)
                vsuma += float(l2) * len(xb); vn += len(xb)
        tr_l2, va_l2 = suma/n, vsuma/vn
        historia.append({"tipo": tipo, "seed": seed, "epoca": epoca, "train_l2": tr_l2, "val_l2": va_l2})
        if va_l2 < mejor - 1e-6:
            mejor, mejor_epoca, espera = va_l2, epoca, 0
            mejor_estado = copy.deepcopy(modelo.state_dict())
        else:
            espera += 1
            if espera >= PATIENCE:
                break
    if mejor_estado is None:
        raise RuntimeError(f"El entrenamiento de {tipo} no guardó un estado de validación.")
    modelo.load_state_dict(mejor_estado)
    return modelo, pd.DataFrame(historia), mejor, mejor_epoca


@torch.no_grad()
def predecir_coef(modelo, tensores, mascara):
    modelo.eval()
    x = torch.tensor(tensores["x"][mascara], device=DEVICE)
    c = torch.tensor(tensores["cov"][mascara], device=DEVICE)
    return modelo(x, c).cpu().numpy()


def metricas_curva(real, pred):
    eps = 1e-8
    rel = np.sqrt(((pred-real)**2).sum(1) / (real**2).sum(1).clip(eps))
    sesgo = 100 * (pred.sum(1)-real.sum(1)) / real.sum(1).clip(eps)
    kge = []
    picos = []
    for y, p in zip(real, pred):
        r = np.corrcoef(y, p)[0, 1] if np.std(y) > 0 and np.std(p) > 0 else 0.0
        alpha = np.std(p) / max(np.std(y), eps)
        beta = np.mean(p) / max(np.mean(y), eps)
        kge.append(1 - math.sqrt((r-1)**2 + (alpha-1)**2 + (beta-1)**2))
        p1y, p1p = int(np.argmax(y[:182])), int(np.argmax(p[:182]))
        p2y, p2p = 182 + int(np.argmax(y[182:])), 182 + int(np.argmax(p[182:]))
        picos.append((abs(p1p-p1y) + abs(p2p-p2y))/2)
    return pd.DataFrame({
        "l2_relativo": rel, "kge": kge, "sesgo_anual_pct": sesgo,
        "error_fecha_picos_dias": picos,
    })


def verificar_kan_sintetica() -> float:
    set_seed(2026)
    x = torch.rand(2400, 2) * 4 - 2
    y = 0.5*x[:, 0] + x[:, 1]**2
    tr, va = torch.arange(1900), torch.arange(1900, 2400)
    capa = CapaKANBSpline(2, 1, n_bases=12, rango=(-2, 2)).to(DEVICE)
    opt = torch.optim.Adam(capa.parameters(), lr=1e-2)
    rmse = float("inf")
    for epoca in range(2000):
        opt.zero_grad(); loss = (capa(x[tr]).squeeze()-y[tr]).square().mean(); loss.backward(); opt.step()
        if (epoca + 1) % 10 == 0:
            with torch.no_grad():
                rmse = torch.sqrt((capa(x[va]).squeeze()-y[va]).square().mean()).item()
            if rmse <= 5e-3:
                break
    if rmse > 1e-2:
        raise RuntimeError(f"KAN sintética no alcanzó RMSE <= 1e-2: {rmse:.5f}")
    return rmse


def inferencia_pareada_por_estacion(tabla, metrica, modelo_a, modelo_b, n_boot=2000, seed=2026):
    columnas = [f"{metrica}_{modelo_a}", f"{metrica}_{modelo_b}"]
    por_estacion = tabla.groupby("station_id")[columnas].mean().dropna()
    diferencias = (por_estacion.iloc[:, 0] - por_estacion.iloc[:, 1]).to_numpy(dtype=float)
    if not len(diferencias):
        raise ValueError(f"No hay estaciones emparejadas para {modelo_a} vs {modelo_b}")
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, len(diferencias), size=(n_boot, len(diferencias)))
    bootstrap = diferencias[indices].mean(axis=1)
    if np.allclose(diferencias, 0):
        estadistico, p_value = 0.0, 1.0
    else:
        estadistico, p_value = np.asarray(wilcoxon(diferencias), dtype=float).tolist()
    return {
        "diferencia_media_por_estacion": float(np.mean(diferencias)),
        "ic95": np.quantile(bootstrap, [0.025, 0.975]).astype(float).tolist(),
        "n_estaciones": int(len(diferencias)),
        "bootstrap_replicas": int(n_boot),
        "wilcoxon_estadistico": float(estadistico),
        "wilcoxon_p_value": float(p_value),
    }


def ejecutar_experimento() -> dict:
    datos = cargar_datos()
    tensores = preparar_tensores(datos)
    fpca = ajustar_fpca(tensores["y"][tensores["train"]], datos["W_est"], N_COMPONENTES_FPCA)
    datos["trazabilidad_emparejamiento"].to_csv(
        OUT / "trazabilidad_emparejamiento.csv", index=False
    )
    datos["estaciones_excluidas"].to_csv(
        OUT / "estaciones_excluidas_emparejamiento.csv", index=False
    )
    datos["mapa_estacion_pixel"].to_csv(
        OUT / "mapa_estacion_pixel.csv", index=False
    )
    estaciones_por_pixel = (
        datos["mapa_estacion_pixel"]
        .groupby(["pixel_row", "pixel_col"], as_index=False)
        .agg(estaciones=("station_id", "nunique"))
    )
    estaciones_por_pixel.to_csv(OUT / "estaciones_por_pixel.csv", index=False)
    rmse_sintetica = verificar_kan_sintetica()

    modelos, historias, filas = {}, [], []
    pred_val = {"lineal": [], "mlp": [], "kan": []}
    for tipo in ["lineal", "mlp", "kan"]:
        for seed in SEEDS:
            modelo, hist, mejor, epoca = entrenar(tipo, seed, datos, tensores, fpca)
            modelos[(tipo, seed)] = modelo
            historias.append(hist)
            pred = predecir_coef(modelo, tensores, tensores["val"])
            pred_val[tipo].append(pred)
            torch.save(modelo.state_dict(), OUT / f"modelo_{tipo}_seed{seed}.pt")
            filas.append({
                "modelo": tipo, "seed": seed, "mejor_epoca": epoca,
                "val_l2_sqrt": mejor,
                "parametros_totales": sum(p.numel() for p in modelo.parameters() if p.requires_grad),
                "parametros_cabeza": sum(p.numel() for p in modelo.head.parameters() if p.requires_grad),
            })
    historia = pd.concat(historias, ignore_index=True)
    ejecuciones = pd.DataFrame(filas)

    mask_val = tensores["val"]
    y_real = np.maximum(datos["y_mm_coef"][mask_val] @ datos["phi_est"].T, 0)
    sat_base = np.maximum(datos["x_mm_coef"][mask_val] @ datos["phi_sat"].T, 0)
    ids_val = datos["ids"].loc[mask_val].reset_index(drop=True)

    metricas_semilla = []
    pred_mm = {}
    for tipo, lista in pred_val.items():
        pred_mm[tipo] = []
        for seed, coef in zip(SEEDS, lista):
            raiz = coef @ datos["phi_est"].T
            mm = np.maximum(raiz, 0)**2
            pred_mm[tipo].append(mm)
            m = metricas_curva(y_real, mm)
            m["modelo"], m["seed"] = tipo, seed
            metricas_semilla.append(m)
    metricas_semilla = pd.concat(metricas_semilla, ignore_index=True)
    met_sat = metricas_curva(y_real, sat_base)
    met_sat["modelo"], met_sat["seed"] = "satelite_sin_corregir", 0
    metricas_semilla = pd.concat([metricas_semilla, met_sat], ignore_index=True)

    resumen = (metricas_semilla.groupby(["modelo", "seed"])
               .agg(l2_relativo=("l2_relativo", "mean"), kge=("kge", "mean"),
                    sesgo_anual_pct=("sesgo_anual_pct", "mean"),
                    error_fecha_picos_dias=("error_fecha_picos_dias", "mean"))
               .reset_index())

    comparacion = ids_val[["station_id", "year", "fase_enso"]].copy()
    pred_ensamble = {
        tipo: np.mean(predicciones, axis=0)
        for tipo, predicciones in pred_mm.items()
    }
    pred_ensamble["satelite"] = sat_base
    metricas_ensamble = {}
    filas_ensamble = []
    for modelo, pred in pred_ensamble.items():
        metricas_modelo = metricas_curva(y_real, pred)
        metricas_ensamble[modelo] = metricas_modelo
        for metrica in metricas_modelo.columns:
            comparacion[f"{metrica}_{modelo}"] = metricas_modelo[metrica].to_numpy()
        filas_ensamble.append({
            "modelo": modelo,
            **{metrica: float(metricas_modelo[metrica].mean())
               for metrica in metricas_modelo.columns},
        })
    ensamble_df = pd.DataFrame(filas_ensamble)

    pares_modelos = [("kan", "mlp"), ("kan", "lineal"), ("kan", "satelite")]
    bootstrap = {}
    filas_inferencia = []
    for modelo_a, modelo_b in pares_modelos:
        clave = f"{modelo_a}_menos_{modelo_b}"
        bootstrap[clave] = {}
        for metrica in metricas_ensamble["kan"].columns:
            resultado = inferencia_pareada_por_estacion(
                comparacion, metrica, modelo_a, modelo_b
            )
            bootstrap[clave][metrica] = resultado
            filas_inferencia.append({
                "comparacion": clave,
                "metrica": metrica,
                **resultado,
                "ic95_inferior": resultado["ic95"][0],
                "ic95_superior": resultado["ic95"][1],
            })
    inferencia_df = pd.DataFrame(filas_inferencia)

    por_enso = []
    for fase, grupo in comparacion.groupby("fase_enso"):
        for modelo in ["lineal", "mlp", "kan", "satelite"]:
            por_enso.append({
                "fase_enso": fase,
                "modelo": modelo,
                "n_curvas": len(grupo),
                "n_estaciones": int(grupo.station_id.nunique()),
                **{
                    metrica: float(grupo[f"{metrica}_{modelo}"].mean())
                    for metrica in metricas_ensamble[modelo].columns
                },
            })
    por_enso = pd.DataFrame(por_enso)

    # Verificación independiente de la forma funcional lineal, incluida la escala afín fija.
    lineal = modelos[("lineal", SEEDS[0])]
    if not isinstance(lineal.head, nn.Linear):
        raise TypeError("La cabeza del modelo lineal debe ser nn.Linear.")
    x_check = torch.tensor(tensores["x"][mask_val][:64], device=DEVICE)
    cov_check = torch.tensor(tensores["cov"][mask_val][:64], device=DEVICE)
    with torch.no_grad():
        pred_red = lineal(x_check, cov_check)
        n_funcionales = N_NEURONAS_FUNCIONALES
        peso_z = lineal.head.weight[:, :n_funcionales]
        peso_cov = lineal.head.weight[:, n_funcionales:]
        escala_z = lineal.escala_caracteristicas[:n_funcionales]
        escala_cov = lineal.escala_caracteristicas[n_funcionales:]
        peso_z_efectivo = peso_z / escala_z[None, :]
        peso_cov_efectivo = peso_cov / escala_cov[None, :]
        intercepto_efectivo = (
            lineal.head.bias
            - (lineal.media_caracteristicas[:n_funcionales] / escala_z) @ peso_z.T
            - (lineal.media_caracteristicas[n_funcionales:] / escala_cov) @ peso_cov.T
        )
        scores_directos = (
            x_check
            @ lineal.encoder.W
            @ (peso_z_efectivo @ lineal.encoder.A).T
            + cov_check @ peso_cov_efectivo.T
            + intercepto_efectivo
        )
        pred_directa = lineal.media_y + scores_directos @ lineal.direcciones
    max_diff_lineal = float((pred_red-pred_directa).abs().max())

    # Porcentaje de raíces negativas recortadas a cero.
    raiz_kan = np.mean(pred_val["kan"], axis=0) @ datos["phi_est"].T
    porcentaje_recortado = float(100 * (raiz_kan < 0).mean())

    mejor_seed = int(
        ejecuciones.loc[ejecuciones.modelo.eq("kan")]
        .sort_values("val_l2_sqrt")
        .iloc[0].seed
    )
    mejor_kan = modelos[("kan", mejor_seed)]
    x_train = torch.tensor(tensores["x"][tensores["train"]], device=DEVICE)
    cov_train = torch.tensor(tensores["cov"][tensores["train"]], device=DEVICE)
    x_val = torch.tensor(tensores["x"][tensores["val"]], device=DEVICE)
    cov_val = torch.tensor(tensores["cov"][tensores["val"]], device=DEVICE)
    with torch.no_grad():
        feat_train = mejor_kan(
            x_train, cov_train, devolver_scores=True
        )[2].cpu().numpy()
        feat_val = mejor_kan(
            x_val, cov_val, devolver_scores=True
        )[2].cpu().numpy()
    limite_inferior = float(mejor_kan.head.nodos[0])
    limite_superior = float(mejor_kan.head.nodos[-1])
    porcentaje_fuera_rango = {
        "entrenamiento": float(
            100 * np.mean((feat_train < limite_inferior) | (feat_train > limite_superior))
        ),
        "validacion": float(
            100 * np.mean((feat_val < limite_inferior) | (feat_val > limite_superior))
        ),
    }

    historia.to_csv(OUT / "historia_entrenamiento.csv", index=False)
    ejecuciones.to_csv(OUT / "ejecuciones.csv", index=False)
    resumen.to_csv(OUT / "metricas_validacion_por_semilla.csv", index=False)
    ensamble_df.to_csv(OUT / "metricas_validacion_ensamble.csv", index=False)
    inferencia_df.to_csv(OUT / "comparaciones_inferenciales.csv", index=False)
    comparacion.to_csv(OUT / "metricas_validacion_por_curva.csv", index=False)
    por_enso.to_csv(OUT / "metricas_validacion_por_enso.csv", index=False)
    predicciones_kan = pd.DataFrame(
        pred_ensamble["kan"],
        columns=[f"dia_{i+1}" for i in range(pred_ensamble["kan"].shape[1])],
    )
    pd.concat([ids_val, predicciones_kan], axis=1).to_csv(
        OUT / "predicciones_validacion_kan.csv.gz",
        index=False,
        compression="gzip",
    )

    # Figuras de diagnóstico e interpretación.
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for seed in SEEDS:
        h = historia[(historia.tipo == "kan") & (historia.seed == seed)]
        ax.plot(h.epoca, h.val_l2, label=f"validación seed {seed}")
    ax.set(title="KAN funcional: pérdida de validación", xlabel="Época", ylabel="L² en escala raíz")
    ax.legend(); ax.grid(alpha=.25); fig.tight_layout()
    fig.savefig(OUT / "perdida_validacion_kan.png", dpi=160); plt.close(fig)

    ejemplo_idx = np.argsort(
        metricas_ensamble["kan"]["l2_relativo"].to_numpy()
    )[-3:][::-1]
    fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
    for ax, i in zip(axes, ejemplo_idx):
        ax.plot(y_real[i], "k", lw=2, label="estación")
        ax.plot(sat_base[i], color="#e76f51", alpha=.8, label="satélite sin corregir")
        ax.plot(pred_ensamble["kan"][i], color="#2a9d8f", label="KAN")
        ax.set_ylabel("mm/día")
        ax.set_title(
            f"{ids_val.loc[i, 'station_id']} — {ids_val.loc[i, 'year']} "
            f"(L2 rel. {metricas_ensamble['kan'].loc[i, 'l2_relativo']:.3f})"
        )
    axes[0].legend(ncol=3); axes[-1].set_xlabel("Día del año")
    fig.tight_layout(); fig.savefig(OUT / "ejemplos_validacion.png", dpi=160); plt.close(fig)
    comparacion.iloc[ejemplo_idx].to_csv(OUT / "peores_casos_validacion.csv", index=False)

    with torch.no_grad():
        w = (mejor_kan.encoder.phi @ mejor_kan.encoder.A.T).cpu().numpy()
    fig, ax = plt.subplots(figsize=(11, 5))
    ax.plot(np.arange(1, 366), w)
    ax.set(title="Funciones de peso del codificador funcional KAN", xlabel="Día del año", ylabel="w_k(t)")
    ax.grid(alpha=.25); fig.tight_layout(); fig.savefig(OUT / "pesos_funcionales.png", dpi=160); plt.close(fig)

    if not isinstance(mejor_kan.head, CapaKANBSpline):
        raise TypeError("La cabeza KAN debe ser una CapaKANBSpline.")
    capa = mejor_kan.head
    nombres = [f"z{k+1}" for k in range(N_NEURONAS_FUNCIONALES)] + ["altitud", "ONI"]
    rangos_caracteristicas = []
    for i, nombre in enumerate(nombres):
        inferior = float(feat_train[:, i].min())
        superior = float(feat_train[:, i].max())
        rangos_caracteristicas.append({
            "caracteristica": nombre,
            "minimo_entrenamiento": inferior,
            "maximo_entrenamiento": superior,
        })

    fig, axes = plt.subplots(capa.n_in, 1, figsize=(12, 2.4 * capa.n_in), sharex=False)
    activaciones = []
    with torch.no_grad():
        for i, ax in enumerate(axes):
            rango = torch.linspace(
                rangos_caracteristicas[i]["minimo_entrenamiento"],
                rangos_caracteristicas[i]["maximo_entrenamiento"],
                240,
            )
            for o in range(capa.n_out):
                valores = capa.edge_function(o, i, rango).cpu().numpy()
                activaciones.append({
                    "amplitud": float(np.ptp(valores)),
                    "salida": o,
                    "entrada": i,
                    "x": rango.cpu().numpy(),
                    "y": valores,
                })
                ax.plot(rango.cpu().numpy(), valores, label=f"score {o+1}")
            ax.set_title(f"Entrada {nombres[i]}: activaciones hacia los scores FPCA")
            ax.set_ylabel("φ(x)")
            ax.grid(alpha=.25)
            ax.legend(ncol=3, fontsize=8)
    axes[-1].set_xlabel("Característica estandarizada")
    fig.tight_layout(); fig.savefig(OUT / "activaciones_kan.png", dpi=160); plt.close(fig)

    activacion_top = max(activaciones, key=lambda fila: fila["amplitud"])
    x_activacion = activacion_top["x"]
    y_activacion = activacion_top["y"]
    grados = range(1, min(8, len(x_activacion) // 2))
    mejor_ajuste = None
    for grado in grados:
        polinomio = np.polynomial.Polynomial.fit(
            x_activacion[::2], y_activacion[::2], deg=grado
        ).convert()
        pred_holdout = polinomio(x_activacion[1::2])
        observado = y_activacion[1::2]
        ss_tot = np.sum((observado - observado.mean()) ** 2)
        r2 = 1.0 - np.sum((observado - pred_holdout) ** 2) / max(ss_tot, 1e-12)
        if mejor_ajuste is None or r2 > mejor_ajuste["r2_holdout"]:
            mejor_ajuste = {
                "grado": grado,
                "r2_holdout": float(r2),
                "coeficientes_polinomio": polinomio.coef.tolist(),
            }
    mejor_ajuste.update({
        "entrada": nombres[activacion_top["entrada"]],
        "score_fpca": int(activacion_top["salida"] + 1),
        "rango_entrenamiento": [
            rangos_caracteristicas[activacion_top["entrada"]]["minimo_entrenamiento"],
            rangos_caracteristicas[activacion_top["entrada"]]["maximo_entrenamiento"],
        ],
        "nota": "Ajuste polinómico de la activación; R2 calculado en puntos alternos retenidos.",
    })
    (OUT / "ajuste_simbolico_activacion.json").write_text(
        json.dumps(mejor_ajuste, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    pd.DataFrame(rangos_caracteristicas).to_csv(
        OUT / "rangos_caracteristicas_kan.csv", index=False
    )

    cabeza_lineal = modelos[("lineal", SEEDS[0])]
    if not isinstance(cabeza_lineal.head, nn.Linear):
        raise TypeError("La cabeza del modelo lineal debe ser nn.Linear.")
    pesos_efectivos = (
        cabeza_lineal.head.weight.detach().cpu().numpy()
        / cabeza_lineal.escala_caracteristicas.detach().cpu().numpy()[None, :]
    )
    coef_lineales = pd.DataFrame([
        {
            "caracteristica": nombre,
            "score_fpca": salida + 1,
            "coeficiente_efectivo": float(pesos_efectivos[salida, entrada]),
        }
        for entrada, nombre in enumerate(nombres)
        for salida in range(pesos_efectivos.shape[0])
    ])
    coef_lineales.to_csv(OUT / "coeficientes_cabeza_lineal.csv", index=False)
    fig, ax = plt.subplots(figsize=(11, 5))
    matriz_coef = pesos_efectivos.T
    imagen = ax.imshow(matriz_coef, aspect="auto", cmap="coolwarm")
    ax.set_xticks(np.arange(matriz_coef.shape[1]), labels=[f"score {i+1}" for i in range(matriz_coef.shape[1])])
    ax.set_yticks(np.arange(len(nombres)), labels=nombres)
    ax.set_title("Coeficientes efectivos de la cabeza lineal (contraste con KAN)")
    fig.colorbar(imagen, ax=ax, label="Coeficiente")
    fig.tight_layout(); fig.savefig(OUT / "coeficientes_cabeza_lineal.png", dpi=160); plt.close(fig)

    params = ejecuciones.groupby("modelo").first()[["parametros_totales", "parametros_cabeza"]].to_dict("index")
    diferencia_parametros = 100 * abs(
        params["kan"]["parametros_totales"] - params["mlp"]["parametros_totales"]
    ) / params["kan"]["parametros_totales"]
    diferencia_parametros_cabeza = 100 * abs(
        params["kan"]["parametros_cabeza"] - params["mlp"]["parametros_cabeza"]
    ) / params["kan"]["parametros_cabeza"]
    assert diferencia_parametros <= 10
    assert diferencia_parametros_cabeza <= 10
    assert max_diff_lineal <= 1e-4

    salida = {
        "fuentes": ["eda_fase0", "fda_fase1", "Datos", "datos_web/ONI_NOAA_CPC_1981_2013.csv"],
        "pares": {
            "total_emparejado": int(len(datos["ids"])),
            "entrenamiento": int(tensores["train"].sum()),
            "validacion": int(tensores["val"].sum()),
            "prueba_cerrada": int(tensores["test"].sum()),
            "estaciones": int(datos["ids"].station_id.nunique()),
            "estaciones_reservadas_docente": "PENDIENTE: no se encontró una lista separada",
        },
        "emparejamiento": {
            "pendiente_confirmacion_docente": True,
            "distancia_pixel_km_mediana": float(datos["mapa_estacion_pixel"].distancia_pixel_km.median()),
            "distancia_pixel_km_maxima": float(datos["mapa_estacion_pixel"].distancia_pixel_km.max()),
            "estaciones_que_comparten_pixel": int((estaciones_por_pixel.estaciones > 1).sum()),
            "trazabilidad": datos["trazabilidad_emparejamiento"].to_dict("records"),
        },
        "arquitectura": {
            "neuronas_funcionales": N_NEURONAS_FUNCIONALES,
            "componentes_fpca": N_COMPONENTES_FPCA,
            "varianza_fpca_acumulada": float(fpca["proporcion"][:N_COMPONENTES_FPCA].sum()),
            "bases_kan_por_arista": KAN_BASES,
            "lambda_rugosidad": LAMBDA_RUGOSIDAD,
            "parametros": params,
            "diferencia_parametros_kan_mlp_pct": float(diferencia_parametros),
            "diferencia_parametros_cabeza_kan_mlp_pct": float(diferencia_parametros_cabeza),
        },
        "verificaciones": {
            "rmse_kan_sintetica": rmse_sintetica,
            "max_diferencia_equivalencia_lineal": max_diff_lineal,
            "rango_kan_fijado_solo_con_entrenamiento": [limite_inferior, limite_superior],
            "porcentaje_caracteristicas_fuera_del_rango_kan": porcentaje_fuera_rango,
            "porcentaje_valores_raiz_recortados": porcentaje_recortado,
            "inferencias_pareadas_por_estacion": bootstrap,
            "prueba_utilizada": False,
            "fase_enso": {
                "fuente": "datos_web/ONI_NOAA_CPC_1981_2013.csv",
                "regla": "promedio anual de los 12 valores; Nino >= 0.5, Nina <= -0.5, Neutral en otro caso",
            },
        },
        "ajuste_simbolico_activacion": mejor_ajuste,
        "mejor_seed_kan_validacion": mejor_seed,
    }
    (OUT / "resumen_experimento.json").write_text(json.dumps(salida, ensure_ascii=False, indent=2), encoding="utf-8")
    return {"resumen": salida, "ejecuciones": ejecuciones, "metricas": resumen, "enso": por_enso}


if __name__ == "__main__":
    resultado = ejecutar_experimento()
    print(json.dumps(resultado["resumen"], ensure_ascii=True, indent=2))
    print("\nMétricas de validación por semilla:\n", resultado["metricas"].to_string(index=False))
    print("\nMétricas por fase ENSO:\n", resultado["enso"].to_string(index=False))

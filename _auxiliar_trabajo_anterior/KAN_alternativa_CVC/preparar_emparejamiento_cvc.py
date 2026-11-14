"""Prepara metadatos geograficos CVC y empareja estaciones con la grilla satelital.

La fuente externa aporta exclusivamente codigo, nombre y localizacion. No se
incorporan registros de precipitacion ajenos a los archivos del proyecto.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


KAN_ROOT = Path(__file__).resolve().parent
ROOT = KAN_ROOT.parent
DATA_DIR = ROOT / "Datos"
SOURCE = DATA_DIR / "fuentes_externas" / "catalogo_estaciones_cvc_datos_gov_co.csv"
PAIRS = ROOT / "fda_fase1" / "pares_estacion_anio_incluidos.csv"
CATALOG_OUT = DATA_DIR / "fuentes_externas" / "catalogo_estaciones_cvc_preparado.csv"
PAIR_OUT = ROOT / "fda_fase1" / "pares_estacion_pixel_cvc.csv"
TRACE_OUT = DATA_DIR / "fuentes_externas" / "trazabilidad_catalogo_cvc.json"

SOURCE_URL = "https://www.datos.gov.co/resource/kfsz-ija2.json"
DATASET_PAGE = (
    "https://www.datos.gov.co/Ambiente-y-Desarrollo-Sostenible/"
    "Cat-logo-de-estaciones-hidroclimatol-gicas-de-la-C/kfsz-ija2/about_data"
)

# MAGNA-SIRGAS / Colombia West zone (EPSG:3115), segun cartografia CVC.
A = 6_378_137.0
INV_F = 298.257222101
LAT0_DEG = 4 + 35 / 60 + 46.3215 / 3600
LON0_DEG = -(77 + 4 / 60 + 39.0285 / 3600)
FALSE_EASTING = 1_000_000.0
FALSE_NORTHING = 1_000_000.0

# Georreferenciacion verificada de las imagenes del proyecto.
GRID_XMIN = -77.55
GRID_YMAX = 5.00
GRID_RES = 0.05
GRID_ROWS = 39
GRID_COLS = 37


def inverse_epsg3115(easting: float, northing: float) -> tuple[float, float]:
    """Convierte MAGNA Colombia Oeste a latitud/longitud (GRS80)."""
    f = 1.0 / INV_F
    e2 = 2 * f - f * f
    ep2 = e2 / (1 - e2)
    lat0 = math.radians(LAT0_DEG)
    lon0 = math.radians(LON0_DEG)

    a0 = 1 - e2 / 4 - 3 * e2**2 / 64 - 5 * e2**3 / 256
    a2 = 3 * e2 / 8 + 3 * e2**2 / 32 + 45 * e2**3 / 1024
    a4 = 15 * e2**2 / 256 + 45 * e2**3 / 1024
    a6 = 35 * e2**3 / 3072
    m0 = A * (
        a0 * lat0
        - a2 * math.sin(2 * lat0)
        + a4 * math.sin(4 * lat0)
        - a6 * math.sin(6 * lat0)
    )
    m = m0 + (northing - FALSE_NORTHING)
    mu = m / (A * a0)
    e1 = (1 - math.sqrt(1 - e2)) / (1 + math.sqrt(1 - e2))
    phi1 = (
        mu
        + (3 * e1 / 2 - 27 * e1**3 / 32) * math.sin(2 * mu)
        + (21 * e1**2 / 16 - 55 * e1**4 / 32) * math.sin(4 * mu)
        + 151 * e1**3 / 96 * math.sin(6 * mu)
        + 1097 * e1**4 / 512 * math.sin(8 * mu)
    )
    sin_phi = math.sin(phi1)
    cos_phi = math.cos(phi1)
    tan_phi = math.tan(phi1)
    n1 = A / math.sqrt(1 - e2 * sin_phi**2)
    r1 = A * (1 - e2) / (1 - e2 * sin_phi**2) ** 1.5
    t1 = tan_phi**2
    c1 = ep2 * cos_phi**2
    d = (easting - FALSE_EASTING) / n1

    lat = phi1 - (n1 * tan_phi / r1) * (
        d**2 / 2
        - (5 + 3 * t1 + 10 * c1 - 4 * c1**2 - 9 * ep2) * d**4 / 24
        + (61 + 90 * t1 + 298 * c1 + 45 * t1**2 - 252 * ep2 - 3 * c1**2)
        * d**6
        / 720
    )
    lon = lon0 + (
        d
        - (1 + 2 * t1 + c1) * d**3 / 6
        + (5 - 2 * c1 + 28 * t1 - 3 * c1**2 + 8 * ep2 + 24 * t1**2)
        * d**5
        / 120
    ) / cos_phi
    return math.degrees(lat), math.degrees(lon)


def haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    radius = 6371.0088
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp = p2 - p1
    dl = math.radians(lon2 - lon1)
    h = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * radius * math.asin(math.sqrt(h))


def nearest_pixel(lat: float, lon: float) -> tuple[int, int, float, float, float]:
    col = int(round((lon - (GRID_XMIN + GRID_RES / 2)) / GRID_RES))
    row = int(round(((GRID_YMAX - GRID_RES / 2) - lat) / GRID_RES))
    col = int(np.clip(col, 0, GRID_COLS - 1))
    row = int(np.clip(row, 0, GRID_ROWS - 1))
    center_lon = GRID_XMIN + (col + 0.5) * GRID_RES
    center_lat = GRID_YMAX - (row + 0.5) * GRID_RES
    distance = haversine_km(lat, lon, center_lat, center_lon)
    return row, col, center_lat, center_lon, distance


def first_non_null(series: pd.Series):
    values = series.dropna()
    return values.iloc[0] if len(values) else np.nan


def main() -> None:
    raw = pd.read_csv(SOURCE, dtype={"codigo_estacion": "string"})
    raw["codigo_estacion"] = raw["codigo_estacion"].str.strip()
    raw["este_x"] = pd.to_numeric(raw["este_x"], errors="coerce")
    raw["norte_y"] = pd.to_numeric(raw["norte_y"], errors="coerce")

    coordinate_variants = (
        raw.dropna(subset=["codigo_estacion", "este_x", "norte_y"])
        .groupby("codigo_estacion")[["este_x", "norte_y"]]
        .apply(lambda x: len(x.drop_duplicates()))
    )
    conflicting_codes = coordinate_variants[coordinate_variants > 1].index.tolist()

    grouped = (
        raw.groupby("codigo_estacion", as_index=False)
        .agg(
            nombre_estacion=("nombre_estacion", first_non_null),
            departamento=("departamento", first_non_null),
            municipio=("municipio", first_non_null),
            cuenca=("cuenca", first_non_null),
            entidad=("entidad", first_non_null),
            tipo_estacion=("tipo_estacion", first_non_null),
            estado=("estado", first_non_null),
            este_x=("este_x", first_non_null),
            norte_y=("norte_y", first_non_null),
            filas_catalogo=("codigo_estacion", "size"),
        )
    )
    grouped["coordenada_conflictiva"] = grouped["codigo_estacion"].isin(conflicting_codes)

    lat_lon = [
        inverse_epsg3115(e, n) if np.isfinite(e) and np.isfinite(n) else (np.nan, np.nan)
        for e, n in zip(grouped["este_x"], grouped["norte_y"])
    ]
    grouped["latitud"] = [x[0] for x in lat_lon]
    grouped["longitud"] = [x[1] for x in lat_lon]
    grouped["crs_origen"] = "MAGNA-SIRGAS / Colombia West zone (EPSG:3115)"
    grouped["fuente_url"] = DATASET_PAGE
    grouped.to_csv(CATALOG_OUT, index=False)

    pairs = pd.read_csv(PAIRS, dtype={"station_id": "string"})
    matched = pairs.merge(
        grouped.rename(columns={"codigo_estacion": "station_id"}),
        on="station_id",
        how="left",
        validate="many_to_one",
        indicator=True,
    )
    matched["catalogo_cvc_encontrado"] = matched["_merge"].eq("both")
    matched = matched.drop(columns="_merge")

    pixels = []
    for lat, lon in zip(matched["latitud"], matched["longitud"]):
        if np.isfinite(lat) and np.isfinite(lon):
            pixels.append(nearest_pixel(float(lat), float(lon)))
        else:
            pixels.append((np.nan, np.nan, np.nan, np.nan, np.nan))
    matched["pixel_fila"] = [x[0] for x in pixels]
    matched["pixel_col"] = [x[1] for x in pixels]
    matched["pixel_centro_lat"] = [x[2] for x in pixels]
    matched["pixel_centro_lon"] = [x[3] for x in pixels]
    matched["distancia_centro_pixel_km"] = [x[4] for x in pixels]
    matched["dentro_extension_satelital"] = (
        matched["latitud"].between(GRID_YMAX - GRID_ROWS * GRID_RES, GRID_YMAX)
        & matched["longitud"].between(GRID_XMIN, GRID_XMIN + GRID_COLS * GRID_RES)
    )
    matched.to_csv(PAIR_OUT, index=False)

    station_match = matched.drop_duplicates("station_id")
    usable_pairs = matched[
        matched["catalogo_cvc_encontrado"]
        & ~matched["coordenada_conflictiva"].fillna(True)
        & matched["dentro_extension_satelital"]
    ]
    trace = {
        "uso_fuente_externa": "solo metadatos de localizacion; no precipitacion externa",
        "fuente_api": SOURCE_URL,
        "pagina_catalogo": DATASET_PAGE,
        "filas_catalogo_descargadas": int(len(raw)),
        "codigos_unicos_catalogo": int(grouped["codigo_estacion"].nunique()),
        "codigos_con_coordenadas_conflictivas": int(len(conflicting_codes)),
        "estaciones_proyecto": int(station_match["station_id"].nunique()),
        "estaciones_encontradas_cvc": int(station_match["catalogo_cvc_encontrado"].sum()),
        "estaciones_sin_catalogo_cvc": int((~station_match["catalogo_cvc_encontrado"]).sum()),
        "estaciones_encontradas_fuera_grilla": int(
            (station_match["catalogo_cvc_encontrado"] & ~station_match["dentro_extension_satelital"]).sum()
        ),
        "pares_estacion_anio_utilizables": int(len(usable_pairs)),
        "pares_por_particion": {
            str(k): int(v) for k, v in usable_pairs["particion"].value_counts().items()
        },
        "distancia_pixel_km": {
            "mediana": float(usable_pairs["distancia_centro_pixel_km"].median()),
            "maximo": float(usable_pairs["distancia_centro_pixel_km"].max()),
        },
    }
    TRACE_OUT.write_text(json.dumps(trace, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(trace, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

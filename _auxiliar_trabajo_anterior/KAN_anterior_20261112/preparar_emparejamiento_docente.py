"""Integra los tres archivos de localizacion entregados por el docente.

El CSV es la fuente operativa de latitud/longitud, el archivo Stata verifica
que el contenido geográfico sea equivalente y el Excel aporta elevación.
No se estiman coordenadas ni elevaciones faltantes.
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
COORD_CSV = DATA_DIR / "CVC_Wheather_Stations_Long_Lat.csv"
COORD_DTA = DATA_DIR / "CVC_Wheather_Stations_Long_Lat.dta"
ELEV_XLSX = DATA_DIR / "CatalogoEstacionesElevacion.xlsx"
F1 = ROOT / "fda_fase1"
E0 = ROOT / "eda_fase0"
F2 = KAN_ROOT / "resultados"

GRID_XMIN = -77.55
GRID_YMAX = 5.00
GRID_RES = 0.05
GRID_ROWS = 39
GRID_COLS = 37


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
    return row, col, center_lat, center_lon, haversine_km(lat, lon, center_lat, center_lon)


def normalize_code(series: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce").round().astype("Int64")
    return numeric.astype("string")


def main() -> None:
    F2.mkdir(exist_ok=True)
    coords_csv = pd.read_csv(COORD_CSV, encoding="latin1")
    coords_csv["station_id"] = normalize_code(coords_csv["CODIGO"])
    coords_csv = coords_csv.rename(
        columns={
            "STATION_NAME": "station_name",
            "CATEGORY": "category",
            "Municipality": "municipality",
            "Lon": "longitude",
            "lat": "latitude",
        }
    )

    coords_dta = pd.read_stata(COORD_DTA)
    coords_dta["station_id"] = normalize_code(coords_dta["station_id"])
    if set(coords_csv["station_id"]) != set(coords_dta["station_id"]):
        raise RuntimeError("CSV y DTA no contienen los mismos códigos")
    check = coords_csv[["station_id", "longitude", "latitude"]].merge(
        coords_dta[["station_id", "lon", "lat"]], on="station_id", validate="one_to_one"
    )
    max_lon_diff = float((check["longitude"] - check["lon"]).abs().max())
    max_lat_diff = float((check["latitude"] - check["lat"]).abs().max())
    if max_lon_diff > 5e-6 or max_lat_diff > 5e-6:
        raise RuntimeError("CSV y DTA difieren más que el redondeo esperado")

    elevation = pd.read_excel(ELEV_XLSX, sheet_name="Hoja1")
    elevation["station_id"] = normalize_code(elevation["codigo"])
    elevation = elevation.rename(
        columns={
            "cuenca": "basin_elevation_file",
            "nombre": "station_name_elevation_file",
            "municipio": "municipality_elevation_file",
            "alt": "elevation_m",
        }
    )
    if elevation["station_id"].duplicated().any():
        raise RuntimeError("El catálogo de elevación tiene códigos duplicados")
    if elevation["elevation_m"].isna().any():
        raise RuntimeError("El catálogo de elevación tiene elevaciones faltantes")

    metadata = coords_csv[
        [
            "station_id",
            "station_name",
            "category",
            "municipality",
            "CORREGIMIENTO",
            "BASIN",
            "STATUS",
            "FIRST_DATE",
            "LAST_DATE",
            "longitude",
            "latitude",
        ]
    ].merge(
        elevation[
            [
                "station_id",
                "elevation_m",
                "basin_elevation_file",
                "station_name_elevation_file",
                "municipality_elevation_file",
            ]
        ],
        on="station_id",
        how="outer",
        validate="one_to_one",
        indicator=True,
    )
    metadata["has_coordinates"] = metadata["longitude"].notna() & metadata["latitude"].notna()
    metadata["has_elevation"] = metadata["elevation_m"].notna()
    metadata["source_coordinates"] = COORD_CSV.name
    metadata["source_coordinates_check"] = COORD_DTA.name
    metadata["source_elevation"] = ELEV_XLSX.name
    metadata.to_csv(F2 / "catalogo_localizacion_docente_preparado.csv", index=False)

    # Pares finales producidos por el EDA/preprocesamiento de Jacobo.
    pairs = pd.read_csv(
        E0 / "pares_estacion_anio_validos_preliminar.csv",
        dtype={"STATION_ID": "string"},
    ).rename(columns={"STATION_ID": "station_id", "YEAR": "year"})
    pairs["year"] = pd.to_numeric(pairs["year"], errors="raise").astype(int)
    pairs["particion"] = np.select(
        [pairs["year"].between(1981, 1995), pairs["year"].between(1996, 2000)],
        ["entrenamiento", "validacion"],
        default="prueba",
    )
    matched = pairs.merge(metadata, on="station_id", how="left", validate="many_to_one")
    matched["inside_satellite_extent"] = (
        matched["latitude"].between(GRID_YMAX - GRID_ROWS * GRID_RES, GRID_YMAX)
        & matched["longitude"].between(GRID_XMIN, GRID_XMIN + GRID_COLS * GRID_RES)
    )
    matched["eligible_kan_all_covariates"] = (
        matched["has_coordinates"].fillna(False)
        & matched["has_elevation"].fillna(False)
        & matched["inside_satellite_extent"]
    )

    pixels = []
    for lat, lon in zip(matched["latitude"], matched["longitude"]):
        if np.isfinite(lat) and np.isfinite(lon):
            pixels.append(nearest_pixel(float(lat), float(lon)))
        else:
            pixels.append((np.nan, np.nan, np.nan, np.nan, np.nan))
    matched["pixel_row"] = [x[0] for x in pixels]
    matched["pixel_col"] = [x[1] for x in pixels]
    matched["pixel_center_lat"] = [x[2] for x in pixels]
    matched["pixel_center_lon"] = [x[3] for x in pixels]
    matched["distance_to_pixel_center_km"] = [x[4] for x in pixels]
    matched.to_csv(F2 / "pares_estacion_pixel_docente.csv", index=False)

    usable = matched[matched["eligible_kan_all_covariates"]].copy()
    usable["pixel_row"] = usable["pixel_row"].astype(int)
    usable["pixel_col"] = usable["pixel_col"].astype(int)
    ground = pd.read_csv(
        F1 / "curvas_estaciones_mm.csv.gz", dtype={"STATION_ID": "string"}
    ).rename(columns={"STATION_ID": "station_id", "YEAR": "year"})
    satellite = pd.read_csv(F1 / "curvas_satelite_mm.csv.gz").rename(
        columns={"anio": "year", "fila": "pixel_row", "col": "pixel_col"}
    )
    ids = [
        "station_id",
        "year",
        "particion",
        "station_name",
        "latitude",
        "longitude",
        "elevation_m",
        "pixel_row",
        "pixel_col",
        "pixel_center_lat",
        "pixel_center_lon",
        "distance_to_pixel_center_km",
    ]
    joined = usable[ids].merge(
        ground, on=["station_id", "year"], validate="one_to_one"
    ).merge(
        satellite,
        on=["year", "pixel_row", "pixel_col"],
        validate="many_to_one",
    )
    if len(joined) != len(usable):
        raise RuntimeError("Se perdieron filas al construir los pares")
    ground_cols = [str(i) for i in range(1, 366)]
    satellite_cols = [c for c in satellite.columns if c.startswith("t")]
    joined[ids].to_csv(F2 / "identificadores_pares_kan_docente.csv", index=False)
    joined[ids + satellite_cols].to_csv(
        F2 / "entradas_satelitales_72_docente.csv.gz", index=False, compression="gzip"
    )
    joined[ids + ground_cols].to_csv(
        F2 / "salidas_terrestres_365_docente.csv.gz", index=False, compression="gzip"
    )

    station_level = matched.drop_duplicates("station_id")
    summary = {
        "uso": "archivos suministrados por el docente; sin estimar faltantes",
        "csv_coordinate_rows": int(len(coords_csv)),
        "dta_coordinate_rows": int(len(coords_dta)),
        "same_coordinate_code_set": True,
        "max_csv_dta_longitude_difference_degrees": max_lon_diff,
        "max_csv_dta_latitude_difference_degrees": max_lat_diff,
        "elevation_rows": int(len(elevation)),
        "project_stations_after_preprocessing": int(station_level["station_id"].nunique()),
        "project_stations_with_coordinates": int(station_level["has_coordinates"].fillna(False).sum()),
        "project_stations_with_elevation": int(station_level["has_elevation"].fillna(False).sum()),
        "project_stations_with_all_covariates_inside_grid": int(
            station_level["eligible_kan_all_covariates"].sum()
        ),
        "usable_station_year_pairs": int(len(joined)),
        "usable_pairs_by_partition": {
            str(k): int(v) for k, v in joined["particion"].value_counts().items()
        },
        "unique_satellite_pixels": int(
            joined[["pixel_row", "pixel_col"]].drop_duplicates().shape[0]
        ),
        "median_pixel_distance_km": float(joined["distance_to_pixel_center_km"].median()),
        "max_pixel_distance_km": float(joined["distance_to_pixel_center_km"].max()),
        "elevation_range_m": [float(joined["elevation_m"].min()), float(joined["elevation_m"].max())],
    }
    (F2 / "resumen_pares_kan_docente.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

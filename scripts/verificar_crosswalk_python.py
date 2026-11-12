"""Verificación independiente (Python + rasterio) del cruce estación–píxel hecho en R.

Comprueba, sin reutilizar el código R:
1. que la fila/columna de cada estación coincide con la que da la geotransformada del TIFF
   (rasterio.index) y que esa celda es de tierra en las curvas satelitales;
2. la concordancia de nombres entre el CSV CVC y el catálogo de elevación;
3. qué estaciones sin coordenada CVC tienen al menos altitud en el catálogo.

Falla (AssertionError) si alguna fila/columna difiere. No modifica los datos fuente.
Uso, desde la raíz del repo:  python scripts/verificar_crosswalk_python.py
"""

from __future__ import annotations

import json
import re
import unicodedata
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "fda_fase2" / "tablas_R"


def normalizar(nombre) -> str:
    s = unicodedata.normalize("NFKD", str(nombre)).encode("ascii", "ignore").decode().upper()
    return " ".join(re.sub(r"[^A-Z0-9 ]", " ", s).split())


def main() -> dict:
    cw = pd.read_csv(TAB / "crosswalk_estacion_pixel.csv", dtype={"STATION_ID": str})
    cvc = pd.read_csv(ROOT / "CVC_Wheather_Stations_Long_Lat.csv", encoding="latin-1", dtype={"CODIGO": str})
    cat = pd.read_excel(ROOT / "CatalogoEstacionesElevacion.xlsx", sheet_name=0, dtype={"codigo": str})
    # 0. El CSV y el DTA de la CVC deben traer los mismos códigos y coordenadas (en R falta `haven`)
    dta = pd.read_stata(ROOT / "CVC_Wheather_Stations_Long_Lat.dta")
    dta["station_id"] = dta["station_id"].astype("int64").astype(str)
    assert set(dta.station_id) == set(cvc.CODIGO) and len(dta) == len(cvc), "CSV y DTA CVC con códigos distintos"
    m = cvc.merge(dta, left_on="CODIGO", right_on="station_id")
    err_lon, err_lat = float((m.Lon - m.lon).abs().max()), float((m.lat_x - m.lat_y).abs().max())
    assert err_lon < 1e-4 and err_lat < 1e-4, f"Coordenadas CSV/DTA distintas: {err_lon}, {err_lat}"
    sat = pd.read_csv(ROOT / "fda_fase1" / "curvas_satelite_mm.csv.gz", usecols=["anio", "fila", "col"])
    tierra = set(map(tuple, sat[["fila", "col"]].drop_duplicates().to_numpy()))

    # 1. fila/col con la geotransformada del raster
    with rasterio.open(ROOT / "ImagenesSatelitePentadal" / "ppt198101.tif") as r:
        assert r.crs.to_epsg() == 4326, r.crs
        geo = cw[cw.en_malla.astype(bool)].copy()
        rc = [r.index(lon, lat) for lon, lat in zip(geo.lon, geo.lat)]
        geo["fila_rio"] = [a for a, _ in rc]
        geo["col_rio"] = [b for _, b in rc]
        alto, ancho = r.height, r.width
    dif = geo[(geo.fila_rio != geo.fila) | (geo.col_rio != geo.col)]
    assert dif.empty, f"{len(dif)} estaciones con fila/col distinta:\n{dif[['STATION_ID','fila','col','fila_rio','col_rio']]}"
    geo["tierra_py"] = [(f, c) in tierra for f, c in zip(geo.fila_rio, geo.col_rio)]
    assert (geo.tierra_py == geo.en_tierra.astype(bool)).all(), "La máscara de tierra no coincide con la del cruce R"

    # 2. nombres CVC frente a catálogo (solo estaciones con ambos)
    ids = set(cw.STATION_ID)
    ambos = (cvc[cvc.CODIGO.isin(ids)].merge(cat, left_on="CODIGO", right_on="codigo"))
    ambos["mismo_nombre"] = ambos.STATION_NAME.map(normalizar) == ambos.nombre.map(normalizar)
    ambos[["CODIGO", "STATION_NAME", "nombre", "Municipality", "municipio", "alt", "mismo_nombre"]] \
        .to_csv(TAB / "concordancia_nombres_cvc_catalogo.csv", index=False)

    # 3. cobertura de metadatos de las 248 estaciones
    con_cvc = ids & set(cvc.CODIGO)
    con_cat = ids & set(cat.codigo)
    resumen = {
        "cvc_csv_vs_dta": {"codigos": len(cvc), "error_max_lon": err_lon, "error_max_lat": err_lat},
        "estaciones_en_curvas": len(ids),
        "con_coordenada_CVC": len(con_cvc),
        "con_altitud_catalogo": len(con_cat),
        "con_coordenada_o_altitud": len(con_cvc | con_cat),
        "solo_altitud_sin_coordenada": len(con_cat - con_cvc),
        "sin_ningun_metadato": len(ids - con_cvc - con_cat),
        "en_malla": int(len(geo)),
        "en_malla_y_tierra": int(geo.tierra_py.sum()),
        "fila_col_iguales_a_rasterio": int(len(geo) - len(dif)),
        "raster": {"filas": alto, "columnas": ancho, "epsg": 4326},
        "cvc_y_catalogo_con_mismo_nombre": int(ambos.mismo_nombre.sum()),
        "cvc_y_catalogo_total": int(len(ambos)),
        "ids_sin_ningun_metadato": sorted(ids - con_cvc - con_cat),
        "ids_solo_altitud": sorted(con_cat - con_cvc),
    }
    (TAB / "verificacion_crosswalk_python.json").write_text(json.dumps(resumen, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({k: v for k, v in resumen.items() if not k.startswith("ids_")}, indent=2, ensure_ascii=False))
    return resumen


if __name__ == "__main__":
    main()

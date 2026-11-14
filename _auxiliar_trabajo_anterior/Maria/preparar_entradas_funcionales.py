from __future__ import annotations

import calendar
import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image


ROOT = Path(__file__).resolve().parent
if not (ROOT / "P_1980_1990.csv").exists():
    ROOT = ROOT.parent
if not (ROOT / "P_1980_1990.csv").exists():
    raise FileNotFoundError("No se encontró la raíz del proyecto junto al script ni en su carpeta padre")
SAT_DIR = ROOT / "ImagenesSatelitePentadal" / "ImagenesSatelitePentadal"
OUT1 = ROOT / "fda_fase1"
OUT2 = ROOT / "fda_fase2"
OUT1.mkdir(exist_ok=True)
OUT2.mkdir(exist_ok=True)

TRAIN_YEARS = range(1981, 1996)
OFFICIAL_YEARS = range(1981, 2004)


def partition(year: pd.Series) -> np.ndarray:
    return np.select(
        [year.between(1981, 1995), year.between(1996, 2000), year.between(2001, 2003)],
        ["ENTRENAMIENTO", "VALIDACION", "PRUEBA"],
        default="FUERA_PERIODO",
    )


def normalize_station_id(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").astype("Int64").astype("string")


def max_inner_gap(group: pd.DataFrame) -> int:
    dates = group["fecha"].sort_values().drop_duplicates()
    if len(dates) < 2:
        return 0
    return int(dates.diff().dt.days.sub(1).clip(lower=0).max())


def prepare_ground() -> dict:
    pieces = []
    for name in ["P_1980_1990.csv", "P_1991_2003.csv"]:
        frame = pd.read_csv(
            ROOT / name,
            usecols=["STATION_ID", "YEAR", "MONTH", "DAY", "VALUE"],
            low_memory=False,
        )
        frame["archivo_origen"] = name
        pieces.append(frame)
    raw = pd.concat(pieces, ignore_index=True)
    raw.columns = [str(c).lower() for c in raw.columns]
    raw = raw.rename(
        columns={
            "station_id": "station_id",
            "year": "year",
            "month": "month",
            "day": "day",
            "value": "precipitacion_original_mm",
        }
    )
    raw["station_id"] = normalize_station_id(raw["station_id"])
    for col in ["year", "month", "day", "precipitacion_original_mm"]:
        raw[col] = pd.to_numeric(raw[col], errors="coerce")
    raw = raw.loc[raw.year.between(1981, 2003)].copy()
    raw["year"] = raw.year.astype("Int64")
    raw["fecha"] = pd.to_datetime(
        dict(year=raw.year, month=raw.month, day=raw.day), errors="coerce"
    )
    raw["particion"] = partition(raw.year.astype(int))
    raw["precipitacion_limpia_mm"] = raw.precipitacion_original_mm
    known_error = (
        raw.station_id.eq("2611600101")
        & raw.fecha.eq(pd.Timestamp("1984-01-10"))
        & raw.precipitacion_original_mm.eq(859375)
    )
    if int(known_error.sum()) != 1:
        raise RuntimeError(f"Se esperaban 1 registro 859375; encontrados {known_error.sum()}")
    raw.loc[known_error, "precipitacion_limpia_mm"] = np.nan

    pair_catalog = raw[["station_id", "year", "particion"]].dropna().drop_duplicates()
    daily = raw.loc[
        raw.station_id.notna()
        & raw.fecha.notna()
        & raw.precipitacion_limpia_mm.notna(),
        ["station_id", "year", "particion", "fecha", "precipitacion_limpia_mm"],
    ].copy()
    duplicates = daily.duplicated(["station_id", "year", "fecha"], keep=False)
    if duplicates.any():
        conflicts = (
            daily.loc[duplicates]
            .groupby(["station_id", "year", "fecha"])
            .precipitacion_limpia_mm.nunique()
        )
        if conflicts.gt(1).any():
            raise RuntimeError("Hay duplicados diarios conflictivos; no se consolidan silenciosamente")
        daily = daily.drop_duplicates(["station_id", "year", "fecha"])

    grouped = daily.groupby(["station_id", "year"], sort=False)
    diag = grouped.agg(
        dias_validos=("fecha", "nunique"),
        primera_fecha=("fecha", "min"),
        ultima_fecha=("fecha", "max"),
    ).reset_index()
    gaps = grouped.apply(max_inner_gap, include_groups=False).rename("hueco_interior_max").reset_index()
    diag = diag.merge(gaps, on=["station_id", "year"], how="left")
    diag = pair_catalog.merge(diag, on=["station_id", "year"], how="left")
    diag["dias_esperados"] = diag.year.map(lambda y: 366 if calendar.isleap(int(y)) else 365)
    diag["dias_validos"] = diag.dias_validos.fillna(0).astype(int)
    diag["cobertura_anual"] = diag.dias_validos / diag.dias_esperados
    start = pd.to_datetime(diag.year.astype(str) + "-01-01")
    end = pd.to_datetime(diag.year.astype(str) + "-12-31")
    diag["faltantes_inicio"] = (diag.primera_fecha - start).dt.days.fillna(diag.dias_esperados).astype(int)
    diag["faltantes_final"] = (end - diag.ultima_fecha).dt.days.fillna(diag.dias_esperados).astype(int)
    diag["hueco_interior_max"] = diag.hueco_interior_max.fillna(diag.dias_esperados).astype(int)
    diag["cumple_continuidad"] = (
        diag.cobertura_anual.ge(0.90)
        & diag.hueco_interior_max.le(30)
        & diag.faltantes_inicio.le(15)
        & diag.faltantes_final.le(15)
    )
    valid_train = diag.loc[
        diag.cumple_continuidad & diag.particion.eq("ENTRENAMIENTO")
    ]
    years_per_station = valid_train.groupby("station_id").year.nunique()
    selected_stations = years_per_station.loc[years_per_station.ge(5)].index
    diag["cumple_minimo_5_anios_train"] = diag.station_id.isin(selected_stations)
    diag["par_elegible_actual"] = diag.cumple_continuidad & diag.cumple_minimo_5_anios_train
    selected_pairs = diag.loc[
        diag.par_elegible_actual, ["station_id", "year", "particion"]
    ].copy()

    selected_daily = daily.merge(
        selected_pairs[["station_id", "year"]], on=["station_id", "year"], how="inner"
    )
    is_feb29 = selected_daily.fecha.dt.month.eq(2) & selected_daily.fecha.dt.day.eq(29)
    feb29_count = int(is_feb29.sum())
    functional_daily = selected_daily.loc[~is_feb29].copy()
    doy = functional_daily.fecha.dt.dayofyear.astype(int)
    after_feb = functional_daily.fecha.dt.month.gt(2)
    leap = functional_daily.year.map(lambda y: calendar.isleap(int(y)))
    functional_daily["dia_funcional_365"] = doy - (leap & after_feb).astype(int)
    if not functional_daily.dia_funcional_365.between(1, 365).all():
        raise RuntimeError("Índice funcional fuera de 1..365")

    matrix = functional_daily.pivot(
        index=["station_id", "year"],
        columns="dia_funcional_365",
        values="precipitacion_limpia_mm",
    ).reindex(columns=range(1, 366))
    selected_index = pd.MultiIndex.from_frame(selected_pairs[["station_id", "year"]])
    matrix = matrix.reindex(selected_index)
    matrix = matrix.reset_index().rename(columns={"station_id": "STATION_ID", "year": "YEAR"})
    matrix.columns = [str(c) for c in matrix.columns]
    if len(matrix) != len(selected_pairs):
        raise RuntimeError("La matriz terrestre perdió pares elegibles")

    ground_path = OUT1 / "curvas_estaciones_mm.csv.gz"
    matrix.to_csv(ground_path, index=False, compression={"method": "gzip", "compresslevel": 6})
    trace = diag.copy()
    trace.to_csv(OUT1 / "trazabilidad_estaciones.csv", index=False)
    selected_pairs.to_csv(OUT1 / "pares_estacion_anio_incluidos.csv", index=False)

    return {
        "archivo": str(ground_path.relative_to(ROOT)),
        "filas": int(len(matrix)),
        "estaciones": int(matrix.STATION_ID.nunique()),
        "columnas_temporales": 365,
        "valores_faltantes": int(matrix.iloc[:, 2:].isna().sum().sum()),
        "feb29_excluidos_representacion_funcional": feb29_count,
        "pares_por_particion": selected_pairs.particion.value_counts().to_dict(),
        "reglas": [
            "periodo 1981-2003",
            "859375 mm -> NA en columna limpia",
            "cobertura diaria >= 90%",
            "hueco interior <= 30 dias",
            "cada borde <= 15 dias",
            "estacion con >=5 anios validos en train",
            "29-feb excluido solo de la malla funcional de 365 dias",
            "totales mensuales y ceros falsos aun no resueltos: no se alteran",
        ],
    }


def pentad_support_nonleap() -> pd.DataFrame:
    rows = []
    offsets = np.cumsum([0] + [calendar.monthrange(2001, m)[1] for m in range(1, 13)])[:-1]
    for p in range(1, 73):
        month = (p - 1) // 6 + 1
        slot = (p - 1) % 6 + 1
        local_start = 1 + 5 * (slot - 1)
        local_end = 5 * slot if slot < 6 else calendar.monthrange(2001, month)[1]
        midpoint = offsets[month - 1] + (local_start + local_end) / 2
        rows.append(
            {
                "pentada": p,
                "mes": month,
                "posicion_mes": slot,
                "dia_inicio_mes": local_start,
                "dia_fin_mes_no_bisiesto": local_end,
                "dia_anual_representativo": midpoint,
                "nombre_columna": f"t{midpoint:g}",
            }
        )
    result = pd.DataFrame(rows)
    if not result.dia_anual_representativo.is_monotonic_increasing:
        raise RuntimeError("Los soportes pentadales no son crecientes")
    return result


def prepare_satellite() -> dict:
    pattern = re.compile(r"^ppt(?P<year>\d{4})(?P<pentad>\d{2})\.tif$", re.I)
    records = []
    auxiliary = 0
    for path in SAT_DIR.glob("*.tif"):
        if path.name.startswith("._"):
            auxiliary += 1
            continue
        match = pattern.match(path.name)
        if match:
            year = int(match.group("year"))
            pentad = int(match.group("pentad"))
            if year in OFFICIAL_YEARS:
                records.append((year, pentad, path))
    inventory = pd.DataFrame(records, columns=["anio", "pentada", "ruta"]).sort_values(["anio", "pentada"])
    counts = inventory.groupby("anio").pentada.nunique()
    if not counts.reindex(list(OFFICIAL_YEARS)).eq(72).all():
        raise RuntimeError(f"Cobertura satelital oficial incompleta:\n{counts}")

    support = pentad_support_nonleap()
    time_columns = support.nombre_columna.tolist()
    frames = []
    mask_frames = []
    nodata_total = 0
    negative_total = 0
    cap200_total = 0
    shape = None
    for year in OFFICIAL_YEARS:
        subset = inventory.loc[inventory.anio.eq(year)].sort_values("pentada")
        layers = []
        for row in subset.itertuples(index=False):
            with Image.open(row.ruta) as image:
                array = np.asarray(image).astype(np.float32)
                nodata_value = float(image.tag_v2.get(42113, -9999))
            if shape is None:
                shape = array.shape
            if array.shape != shape:
                raise RuntimeError("Las imágenes no comparten la misma forma")
            nodata_mask = array == nodata_value
            nodata_total += int(nodata_mask.sum())
            array[nodata_mask] = np.nan
            negative_total += int(np.nansum(array < 0))
            cap200_total += int(np.nansum(array == 200))
            layers.append(array)
        cube = np.stack(layers, axis=0)  # 72 x rows x cols
        n_rows, n_cols = shape
        values = cube.transpose(1, 2, 0).reshape(n_rows * n_cols, 72)
        ids = pd.DataFrame(
            {
                "anio": year,
                "fila": np.repeat(np.arange(n_rows, dtype=int), n_cols),
                "col": np.tile(np.arange(n_cols, dtype=int), n_rows),
            }
        )
        frames.append(pd.concat([ids, pd.DataFrame(values, columns=time_columns)], axis=1))
        valid_pixel = np.isfinite(values).any(axis=1)
        mask_frames.append(ids.loc[valid_pixel].copy())

    satellite = pd.concat(frames, ignore_index=True)
    land_mask = pd.concat(mask_frames, ignore_index=True)
    sat_path = OUT1 / "curvas_satelite_mm.csv.gz"
    mask_path = OUT2 / "coef_satelite_tierra.csv.gz"
    satellite.to_csv(sat_path, index=False, compression={"method": "gzip", "compresslevel": 6})
    land_mask.to_csv(mask_path, index=False, compression={"method": "gzip", "compresslevel": 6})
    support.to_csv(OUT1 / "soporte_temporal_satelite.csv", index=False)
    inventory.drop(columns="ruta").to_csv(OUT1 / "inventario_satelite_oficial.csv", index=False)

    return {
        "archivo": str(sat_path.relative_to(ROOT)),
        "archivo_mascara": str(mask_path.relative_to(ROOT)),
        "filas": int(len(satellite)),
        "pixeles_por_anio": int(np.prod(shape)),
        "forma_raster": list(shape),
        "columnas_temporales": 72,
        "anios": [min(OFFICIAL_YEARS), max(OFFICIAL_YEARS)],
        "nodata": nodata_total,
        "negativos": negative_total,
        "valores_igual_200": cap200_total,
        "archivos_auxiliares_ignorados": auxiliary,
        "reglas": [
            "solo TIFF con nombre pptYYYYPP.tif",
            "archivos ._ ignorados",
            "periodo 1981-2003",
            "72 capas por anio",
            "escala original conservada; unidad pendiente de confirmar",
            "valor 200 conservado y reportado; no se corrige",
            "soporte temporal provisional: seis periodos por mes, punto medio no bisiesto",
            "fila y columna con indice base cero",
        ],
    }


def main() -> None:
    ground = prepare_ground()
    satellite = prepare_satellite()
    manifest = {
        "proposito": "Entradas reproducibles para construccion y EDA funcional en R",
        "estado": "PROVISIONAL: faltan decisiones sobre totales mensuales, ceros falsos y unidad/calendario satelital",
        "estaciones": ground,
        "satelite": satellite,
        "oni": {
            "archivo_esperado": "eda_fase0/oni_anual_fase.csv",
            "estado": "NO CREADO: requiere fuente oficial y regla de clasificacion documentada",
        },
    }
    (OUT1 / "manifiesto_preparacion.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    summary = pd.DataFrame(
        [
            {
                "fuente": "estaciones",
                "filas": ground["filas"],
                "sitios": ground["estaciones"],
                "tiempos": ground["columnas_temporales"],
                "estado_unidad": "mm/dia",
            },
            {
                "fuente": "satelite",
                "filas": satellite["filas"],
                "sitios": satellite["pixeles_por_anio"],
                "tiempos": satellite["columnas_temporales"],
                "estado_unidad": "escala_original_pendiente",
            },
        ]
    )
    summary.to_csv(OUT1 / "resumen_matrices.csv", index=False)
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

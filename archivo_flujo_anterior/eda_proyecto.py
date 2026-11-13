"""EDA inicial para los datos del proyecto de precipitación.

Ejecutar desde la carpeta del proyecto:
    python eda_proyecto.py

Genera tablas CSV y un resumen Markdown en ./eda_resultados/.
No modifica los archivos originales.
"""

from collections import Counter
from pathlib import Path
import zipfile

import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "eda_resultados"
OUT.mkdir(exist_ok=True)

FILES = [
    ("estaciones_1980_1990", ROOT / "P_1980_1990.csv", "csv"),
    ("estaciones_1991_2003", ROOT / "P_1991_2003.csv", "csv"),
    ("estaciones_1930_1979", ROOT / "P_1930_1979.xlsx", "xlsx"),
]


def inspect_table(name, path, kind):
    if kind == "csv":
        data = pd.read_csv(path, dtype={"STATION_ID": "string"})
    else:
        data = pd.read_excel(path, sheet_name="Hoja1", dtype={"STATION_ID": "string"})

    required = {"STATION_ID", "YEAR", "MONTH", "DAY", "VALUE"}
    missing = required - set(data.columns)
    if missing:
        raise ValueError(f"{path.name}: faltan columnas {sorted(missing)}")

    for col in ["YEAR", "MONTH", "DAY", "VALUE"]:
        data[col] = pd.to_numeric(data[col], errors="coerce")
    valid = data.dropna(subset=["STATION_ID", "YEAR", "MONTH", "DAY", "VALUE"]).copy()
    valid["YEAR"] = valid["YEAR"].astype(int)
    valid["STATION_ID"] = valid["STATION_ID"].astype(str)
    key = ["STATION_ID", "YEAR", "MONTH", "DAY"]
    duplicate_rows = int(valid.duplicated(key, keep=False).sum())

    valid["DATE"] = pd.to_datetime(
        dict(year=valid.YEAR, month=valid.MONTH, day=valid.DAY), errors="coerce"
    )
    bad_dates = int(valid["DATE"].isna().sum())
    valid["MONTH"] = valid["MONTH"].astype(int)

    summary = {
        "source": name,
        "file": path.name,
        "rows_raw": len(data),
        "rows_valid": len(valid),
        "stations": valid.STATION_ID.nunique(),
        "year_min": int(valid.YEAR.min()),
        "year_max": int(valid.YEAR.max()),
        "value_min": float(valid.VALUE.min()),
        "value_median": float(valid.VALUE.median()),
        "value_mean": float(valid.VALUE.mean()),
        "value_p99": float(valid.VALUE.quantile(.99)),
        "value_p999": float(valid.VALUE.quantile(.999)),
        "value_max": float(valid.VALUE.max()),
        "zero_count": int((valid.VALUE == 0).sum()),
        "zero_percent": float((valid.VALUE == 0).mean() * 100),
        "negative_count": int((valid.VALUE < 0).sum()),
        "over_500_count": int((valid.VALUE > 500).sum()),
        "duplicate_key_rows": duplicate_rows,
        "invalid_calendar_dates": bad_dates,
        "missing_required_values": int(len(data) - len(valid)),
    }

    # Cada fila de la cobertura representa una estación en un año.
    coverage = valid.groupby(["STATION_ID", "YEAR"], as_index=False).agg(
        records=("VALUE", "size"),
        unique_dates=("DATE", "nunique"),
        zero_days=("VALUE", lambda x: int((x == 0).sum())),
        annual_total_observed=("VALUE", "sum"),
    )
    coverage["expected_days"] = coverage.YEAR.map(
        lambda year: 366 if pd.Timestamp(year=int(year), month=12, day=31).dayofyear == 366 else 365
    )
    coverage["coverage_percent"] = coverage.unique_dates / coverage.expected_days * 100
    coverage["missing_days_vs_calendar"] = coverage.expected_days - coverage.unique_dates
    coverage.insert(0, "source", name)
    coverage.to_csv(OUT / f"cobertura_{name}.csv", index=False)

    # Valores muy altos se marcan como alertas exploratorias, no como errores confirmados.
    alerts = valid.loc[valid.VALUE > 500, ["STATION_ID", "YEAR", "MONTH", "DAY", "VALUE"]].copy()
    alerts.insert(0, "source", name)
    alerts.to_csv(OUT / f"alertas_mayores_500_{name}.csv", index=False)

    return summary, coverage


def inspect_satellite(path):
    with zipfile.ZipFile(path) as archive:
        files = [n for n in archive.namelist() if n.lower().endswith(".tif") and "/._" not in n]
    per_year = Counter(Path(n).stem[3:7] for n in files)
    rows = [{"year": int(year), "tiff_count": count, "expected_72": count == 72}
            for year, count in sorted(per_year.items())]
    pd.DataFrame(rows).to_csv(OUT / "conteo_tiff_por_anio.csv", index=False)
    return {
        "source": "satelite_pentadal",
        "file": path.name,
        "data_tiff_count": len(files),
        "year_min": min(map(int, per_year)),
        "year_max": max(map(int, per_year)),
        "years_with_72_tiffs": sum(count == 72 for count in per_year.values()),
        "years_with_other_tiff_count": sum(count != 72 for count in per_year.values()),
        "tiff_count_2013": per_year.get("2013", 0),
        "note": "El script cuenta archivos; no interpreta valores raster ni sus unidades.",
    }


def main():
    summaries = []
    all_coverage = []
    for name, path, kind in FILES:
        if not path.exists():
            print(f"AVISO: no se encontró {path.name}; se omite")
            continue
        print(f"Leyendo {path.name}...")
        summary, coverage = inspect_table(name, path, kind)
        summaries.append(summary)
        all_coverage.append(coverage)

    satellite = ROOT / "ImagenesSatelitePentadal.zip"
    if satellite.exists():
        print(f"Inventariando {satellite.name}...")
        summaries.append(inspect_satellite(satellite))

    pd.DataFrame(summaries).to_csv(OUT / "resumen_fuentes.csv", index=False)
    if all_coverage:
        coverage = pd.concat(all_coverage, ignore_index=True)
        coverage.to_csv(OUT / "cobertura_estaciones_todos_los_archivos.csv", index=False)

    summary_df = pd.DataFrame(summaries)
    table_lines = []
    if not summary_df.empty:
        columns = list(summary_df.columns)
        table_lines = [
            "| " + " | ".join(columns) + " |",
            "| " + " | ".join("---" for _ in columns) + " |",
        ]
        for row in summary_df.itertuples(index=False, name=None):
            table_lines.append("| " + " | ".join(str(value) for value in row) + " |")
    lines = [
        "# EDA inicial del proyecto de precipitación", "",
        "Este inventario describe estructura y calidad básica. No elimina ni corrige datos.",
        "Los valores atípicos son alertas que requieren revisión con metadatos.", "",
        "## Resumen de fuentes", "",
        *table_lines, "",
        "## Archivos generados", "",
        "- `resumen_fuentes.csv`: métricas generales por archivo.",
        "- `cobertura_estaciones_todos_los_archivos.csv`: cobertura por estación y año.",
        "- `conteo_tiff_por_anio.csv`: cantidad de imágenes satelitales por año.",
        "- `alertas_mayores_500_*.csv`: observaciones superiores a 500, para inspección.", "",
        "Nota: el número de días esperados es calendario; no significa que una estación deba tener todos esos registros."
    ]
    (OUT / "informe_eda.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"\nEDA terminado. Resultados en: {OUT}")


if __name__ == "__main__":
    main()

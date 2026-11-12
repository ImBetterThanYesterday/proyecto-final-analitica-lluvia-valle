"""Compare generated R constructor outputs with locally stored Python references.

The satellite NPZ is float32, so its comparison tolerance reflects serialization
precision; a freshly generated `coef_satelite_python.csv.gz` is preferred.
"""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
F1 = ROOT / "fda_fase1"


def table(path: Path) -> tuple[list[str], list[list[str]]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as f:
        reader = csv.reader(f)
        return next(reader), list(reader)


def key_text(values: list[str]) -> tuple[str, ...]:
    return tuple(str(int(float(v))) if v not in {"", "NA"} else v for v in values)


def compare_rows(py_rows, r_rows, key_count, value_tolerance, label):
    assert len(py_rows) == len(r_rows), f"{label}: row count mismatch {len(py_rows)} != {len(r_rows)}"
    py_keys = [key_text(row[:key_count]) for row in py_rows]
    r_keys = [key_text(row[:key_count]) for row in r_rows]
    assert py_keys == r_keys, f"{label}: keys/order differ"
    py_vals = np.asarray([[float(v) for v in row[key_count:]] for row in py_rows])
    r_vals = np.asarray([[float(v) for v in row[key_count:]] for row in r_rows])
    error = float(np.max(np.abs(py_vals - r_vals)))
    assert error < value_tolerance, f"{label}: max absolute error {error} >= {value_tolerance}"
    return error


def main():
    _, py_est = table(F1 / "coef_estaciones.csv")
    _, r_est = table(F1 / "coef_estaciones_R.csv")
    err_est = compare_rows(py_est, r_est, 2, 1e-6, "station coefficients")

    sat_csv = F1 / "coef_satelite_python.csv.gz"
    if sat_csv.exists():
        _, py_sat = table(sat_csv)
        sat_tol = 1e-7
        sat_source = "fresh Python CSV"
    else:
        ref = np.load(F1 / "coeficientes_satelite.npz")
        keys = np.column_stack([ref["anio"], ref["fila"], ref["col"]]).astype(str)
        vals = ref["coef"].astype(float)
        py_sat = [list(k) + list(map(str, v)) for k, v in zip(keys, vals)]
        sat_tol = 2e-5
        sat_source = "committed Python float32 NPZ (notebook rerun pending)"
    _, r_sat = table(F1 / "coef_satelite_R.csv.gz")
    err_sat = compare_rows(py_sat, r_sat, 3, sat_tol, "satellite coefficients")

    spec = (F1 / "base_especificacion.json").read_text(encoding="utf-8")
    import json
    expected_k = json.loads(spec)
    summary_head, summary_rows = table(F1 / "resumen_seleccion_K.csv")
    summary = [dict(zip(summary_head, r)) for r in summary_rows]
    by_source = [r for r in summary if r.get("escala", "").startswith("mm/día")]
    station_row = next(r for r in by_source if r["fuente"].lower().startswith("estaciones"))
    satellite_row = next(r for r in by_source if r["fuente"].lower().startswith("sat"))
    assert int(expected_k["K_estaciones"]) == int(station_row["K por VCG"])
    assert int(expected_k["K_satelite"]) == int(satellite_row["K por VCG"])

    W_npz = np.load(F1 / "matrices_gram.npz")
    W_est_r = np.loadtxt(F1 / "gram_estaciones_R.csv", delimiter=",")
    W_sat_r = np.loadtxt(F1 / "gram_satelite_R.csv", delimiter=",")
    err_w_est = float(np.max(np.abs(W_est_r - W_npz["W_est"])))
    err_w_sat = float(np.max(np.abs(W_sat_r - W_npz["W_sat"])))
    assert err_w_est < 1e-8, f"W stations mismatch: {err_w_est}"
    assert err_w_sat < 1e-8, f"W satellite mismatch: {err_w_sat}"
    print(f"OK keys/K | station coeff maxerr={err_est:.3g} | satellite maxerr={err_sat:.3g} ({sat_source})")
    print(f"OK W | station maxerr={err_w_est:.3g} | satellite maxerr={err_w_sat:.3g}")


if __name__ == "__main__":
    main()

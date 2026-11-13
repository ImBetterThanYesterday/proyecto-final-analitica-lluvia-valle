"""Compare generated R constructor outputs with locally stored Python references.

Satellite: same-input comparison (<1e-8) plus an exact bound for the 5-decimal rounding
of the exported curves (notebook fits the float32 TIFF values).
"""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
F1 = ROOT / "fda_fase1"


def table(path: Path) -> tuple[list[str], list[list[str]]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as f:
        reader = csv.reader(f)
        return next(reader), list(reader)


def base_bspline(t, rango, K, orden=4):
    """Cox–de Boor (copia independiente de la del notebook 03, solo para esta prueba)."""
    a, b = rango
    tau = np.r_[[a] * orden, np.linspace(a, b, K - orden + 2)[1:-1], [b] * orden]
    t = np.asarray(t, float)
    B = np.zeros((len(t), len(tau) - 1))
    for j in range(len(tau) - 1):
        B[:, j] = (t >= tau[j]) & (t < tau[j + 1])
    B[t == b, np.max(np.flatnonzero(tau[:-1] < tau[1:]))] = 1
    for m in range(2, orden + 1):
        Bn = np.zeros((len(t), len(tau) - m))
        for j in range(len(tau) - m):
            d1, d2 = tau[j + m - 1] - tau[j], tau[j + m] - tau[j + 1]
            Bn[:, j] = ((t - tau[j]) / d1 * B[:, j] if d1 > 0 else 0) + \
                       ((tau[j + m] - t) / d2 * B[:, j + 1] if d2 > 0 else 0)
        B = Bn
    return B


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

    # Satélite. El notebook ajusta con los valores float32 del TIFF; R ajusta con
    # `curvas_satelite_mm.csv.gz`, redondeado a 5 decimales. Por eso se hacen dos pruebas estrictas:
    #  (a) misma entrada: Python (Cox–de Boor + y2cMap propios) sobre el CSV vs R  -> < 1e-8
    #  (b) notebook vs R: |Δc| ≤ max_k Σ_j |y2cMap_kj| · 5e-6 (cota exacta del redondeo de la entrada)
    sat_csv = F1 / "coef_satelite_python.csv.gz"
    assert sat_csv.exists(), "Falta coef_satelite_python.csv.gz: ejecute el notebook 03"
    _, py_sat = table(sat_csv)
    _, r_sat = table(F1 / "coef_satelite_R.csv.gz")
    spec_sat = json.loads((F1 / "base_especificacion.json").read_text(encoding="utf-8"))
    head_c, curvas = table(F1 / "curvas_satelite_mm.csv.gz")
    t_sat = np.array([float(c[1:]) for c in head_c[3:]])
    Y_sat = np.asarray([[float(v) for v in row[3:]] for row in curvas])
    assert np.isfinite(Y_sat).all(), "Curvas satelitales con NA: la prueba de misma entrada supone curvas completas"
    Phi = base_bspline(t_sat, (spec_sat["rango"][0], spec_sat["rango"][1]), int(spec_sat["K_satelite"]))
    mapa = np.linalg.solve(Phi.T @ Phi, Phi.T)
    misma_entrada = [list(row[:3]) + [repr(float(v)) for v in c] for row, c in zip(curvas, Y_sat @ mapa.T)]
    err_sat_misma = compare_rows(misma_entrada, r_sat, 3, 1e-8, "satellite coefficients, same input")
    cota = float(np.abs(mapa).sum(1).max() * 5e-6)
    err_sat = compare_rows(py_sat, r_sat, 3, cota, "satellite coefficients, notebook vs R (rounding bound)")
    sat_source = f"same-input maxerr={err_sat_misma:.3g} (<1e-8); notebook float32 input vs CSV input bound={cota:.3g}"

    spec = (F1 / "base_especificacion.json").read_text(encoding="utf-8")
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

"""Reconcile baseline/Jacobo station-year sets using only local Git/data files."""

from __future__ import annotations

import csv
import gzip
import io
import subprocess
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MATRIX_REL = "fda_fase1/curvas_estaciones_mm.csv.gz"


def _keys_from_csv_gzip(data: bytes) -> set[tuple[str, str]]:
    with gzip.GzipFile(fileobj=io.BytesIO(data)) as stream:
        rows = csv.reader(io.TextIOWrapper(stream, encoding="utf-8"))
        next(rows)
        return {(row[0], str(int(row[1]))) for row in rows}


def cargar_claves(root: Path, baseline_ref: str = "main") -> tuple[set[tuple[str, str]], set[tuple[str, str]]]:
    """Load matrix keys from the local baseline Git object and current Jacobo checkout."""
    baseline_blob = subprocess.check_output(
        ["git", "show", f"{baseline_ref}:{MATRIX_REL}"], cwd=root
    )
    baseline = _keys_from_csv_gzip(baseline_blob)
    with (root / MATRIX_REL).open("rb") as f:
        jacobo = _keys_from_csv_gzip(f.read())
    return baseline, jacobo


def _coverage_index(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        return {(r["STATION_ID"], str(int(r["YEAR"]))): r for r in csv.DictReader(f)}


def _flags_by_pair(path: Path) -> dict[tuple[str, str], set[str]]:
    out: dict[tuple[str, str], set[str]] = defaultdict(set)
    with path.open(encoding="utf-8-sig", newline="") as f:
        for r in csv.DictReader(f):
            year = r["DATE"][:4]
            out[(r["STATION_ID"], year)].add(r["flag"])
    return out


def _selected_stations(path: Path) -> set[str]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        return {r["STATION_ID"] for r in csv.DictReader(f)}


def reconcile(root: Path = ROOT, baseline_ref: str = "main") -> tuple[list[dict[str, str]], dict[str, int]]:
    baseline, jacobo = cargar_claves(root, baseline_ref)
    coverage = _coverage_index(root / "eda_fase0/cobertura_estacion_anio_1981_2003_tras_depuracion.csv")
    flags = _flags_by_pair(root / "eda_fase0/dias_marcados.csv")
    selected_stations = _selected_stations(root / "eda_fase0/pares_estacion_anio_validos_preliminar.csv")
    baseline_stations = {s for s, _ in baseline}
    records: list[dict[str, str]] = []
    reasons: Counter[str] = Counter()

    for station, year in sorted(baseline - jacobo):
        row = coverage.get((station, year))
        pair_flags = sorted(flags.get((station, year), set()))
        if row is not None and row.get("ok", "").lower() == "false":
            reason = "P5_no_cumple_cobertura_hueco_o_meses"
        elif row is not None and station not in selected_stations:
            reason = "P6_estacion_no_alcanza_minimo_de_anios_de_entrenamiento"
        elif row is not None:
            reason = "revisar_regla_de_seleccion_porque_P5_pasa"
        elif pair_flags:
            reason = "sin_curva_tras_banderas:" + "+".join(pair_flags)
        else:
            reason = "no_aparece_en_cobertura_ni_banderas_revisar_P0_a_P4"
        reasons[reason] += 1
        records.append({"STATION_ID": station, "YEAR": year, "pertenece": "solo_main",
                        "diagnostico": reason, "flags_locales": ";".join(pair_flags),
                        "cobertura": row.get("cob", "") if row else "",
                        "hueco_max": row.get("hueco", "") if row else "",
                        "meses20": row.get("meses20", "") if row else ""})

    for station, year in sorted(jacobo - baseline):
        row = coverage.get((station, year))
        reason = ("estacion_nueva_en_Jacobo_reconciliar_con_reglas_baseline"
                  if station not in baseline_stations else "nuevo_anio_de_estacion_existente_reconciliar_con_reglas_baseline")
        reasons[reason] += 1
        records.append({"STATION_ID": station, "YEAR": year, "pertenece": "solo_jacobo",
                        "diagnostico": reason, "flags_locales": ";".join(sorted(flags.get((station, year), set()))),
                        "cobertura": row.get("cob", "") if row else "",
                        "hueco_max": row.get("hueco", "") if row else "",
                        "meses20": row.get("meses20", "") if row else ""})

    summary = {
        "baseline_curves": len(baseline),
        "jacobo_curves": len(jacobo),
        "common": len(baseline & jacobo),
        "baseline_only": len(baseline - jacobo),
        "jacobo_only": len(jacobo - baseline),
        "baseline_stations": len(baseline_stations),
        "jacobo_stations": len({s for s, _ in jacobo}),
        "jacobo_station_ids_not_in_baseline": len({s for s, _ in jacobo} - baseline_stations),
        "baseline_station_ids_not_in_jacobo": len(baseline_stations - {s for s, _ in jacobo}),
    }
    return records, {**summary, **{f"reason::{k}": v for k, v in sorted(reasons.items())}}


def main() -> None:
    records, summary = reconcile()
    out = ROOT / "docs/diff_claves_estacion_main_jacobo.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    fields = ["STATION_ID", "YEAR", "pertenece", "diagnostico", "flags_locales", "cobertura", "hueco_max", "meses20"]
    with out.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        writer.writeheader(); writer.writerows(records)
    for key, value in summary.items():
        print(f"{key}: {value}")
    print("Detalle:", out.relative_to(ROOT))


if __name__ == "__main__":
    main()

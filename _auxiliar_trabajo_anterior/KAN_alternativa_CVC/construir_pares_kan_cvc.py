"""Construye matrices emparejadas para el piloto KAN con metadatos CVC."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


KAN_ROOT = Path(__file__).resolve().parent
ROOT = KAN_ROOT.parent
F1 = ROOT / "fda_fase1"
F2 = ROOT / "fda_fase2"


def main() -> None:
    F2.mkdir(exist_ok=True)
    mapping = pd.read_csv(
        F1 / "pares_estacion_pixel_cvc.csv", dtype={"station_id": "string"}
    )
    usable = mapping[
        mapping["catalogo_cvc_encontrado"]
        & ~mapping["coordenada_conflictiva"].fillna(True)
        & mapping["dentro_extension_satelital"]
    ].copy()
    usable["pixel_fila"] = usable["pixel_fila"].astype(int)
    usable["pixel_col"] = usable["pixel_col"].astype(int)

    ground = pd.read_csv(
        F1 / "curvas_estaciones_mm.csv.gz", dtype={"STATION_ID": "string"}
    ).rename(columns={"STATION_ID": "station_id", "YEAR": "year"})
    satellite = pd.read_csv(F1 / "curvas_satelite_mm.csv.gz").rename(
        columns={"anio": "year", "fila": "pixel_fila", "col": "pixel_col"}
    )

    ids = [
        "station_id",
        "year",
        "particion",
        "nombre_estacion",
        "latitud",
        "longitud",
        "pixel_fila",
        "pixel_col",
        "pixel_centro_lat",
        "pixel_centro_lon",
        "distancia_centro_pixel_km",
    ]
    base = usable[ids].merge(
        ground,
        on=["station_id", "year"],
        how="inner",
        validate="one_to_one",
    )
    joined = base.merge(
        satellite,
        on=["year", "pixel_fila", "pixel_col"],
        how="inner",
        validate="many_to_one",
        suffixes=("_tierra", "_satelite"),
    )
    if len(joined) != len(usable):
        raise RuntimeError(
            f"Se esperaban {len(usable)} pares y se construyeron {len(joined)}"
        )

    ground_cols = [str(i) for i in range(1, 366)]
    satellite_cols = [col for col in satellite.columns if col.startswith("t")]
    key_cols = ids

    identifiers = joined[key_cols].copy()
    satellite_out = joined[key_cols + satellite_cols].copy()
    ground_out = joined[key_cols + ground_cols].copy()

    identifiers.to_csv(F2 / "identificadores_pares_kan_cvc.csv", index=False)
    satellite_out.to_csv(
        F2 / "entradas_satelitales_72_cvc.csv.gz", index=False, compression="gzip"
    )
    ground_out.to_csv(
        F2 / "salidas_terrestres_365_cvc.csv.gz", index=False, compression="gzip"
    )

    summary = {
        "pares_totales": int(len(joined)),
        "estaciones": int(joined["station_id"].nunique()),
        "pixeles_unicos": int(
            joined[["pixel_fila", "pixel_col"]].drop_duplicates().shape[0]
        ),
        "pares_por_particion": {
            str(k): int(v) for k, v in joined["particion"].value_counts().items()
        },
        "estaciones_por_particion": {
            str(k): int(v)
            for k, v in joined.groupby("particion")["station_id"].nunique().items()
        },
        "variables_satelitales": int(len(satellite_cols)),
        "variables_terrestres": int(len(ground_cols)),
        "regla_espacial": "pixel cuyo centro es el mas cercano a la estacion",
        "covariables_piloto": "ninguna; no se usan altitud ni ONI",
    }
    (F2 / "resumen_pares_kan_cvc.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

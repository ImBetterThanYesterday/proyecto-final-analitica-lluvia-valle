from __future__ import annotations

from pathlib import Path
import calendar

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
FDA = ROOT / "fda_fase1"

# Las coordenadas se usan solo para localizar el pixel. No se incorporan
# precipitaciones de fuentes externas ni se reemplazan observaciones terrestres.
STATIONS = {
    "2612700202": {"lat": 3.4971436609, "lon": -76.3135887},
    "2612800104": {"lat": 3.622057350, "lon": -76.10813753},
    "2622330105": {"lat": 3.427707383, "lon": -76.56881179},
    "2623500101": {"lat": 3.929589828, "lon": -76.41868528},
}

XMIN = -77.55
YMAX = 5.00
RES = 0.05


def station_pixel(lat: float, lon: float) -> tuple[int, int]:
    row = int(np.floor((YMAX - lat) / RES))
    col = int(np.floor((lon - XMIN) / RES))
    return row, col


def month_day_columns(month: int) -> list[str]:
    start = sum(calendar.monthrange(2001, m)[1] for m in range(1, month)) + 1
    end = start + calendar.monthrange(2001, month)[1] - 1
    return [str(day) for day in range(start, end + 1)]


def six_ground_blocks(values: np.ndarray) -> np.ndarray:
    # Seis bloques compatibles con las seis imagenes mensuales:
    # 1-5, 6-10, 11-15, 16-20, 21-25 y 26-fin de mes.
    cuts = [(0, 5), (5, 10), (10, 15), (15, 20), (20, 25), (25, len(values))]
    return np.asarray([np.nansum(values[a:b]) for a, b in cuts], dtype=float)


def shares(values: np.ndarray) -> np.ndarray:
    total = float(np.nansum(values))
    if total <= 0:
        return np.full(6, np.nan)
    return values / total


def main() -> None:
    audit = pd.read_csv(FDA / "auditoria_posibles_totales_mensuales.csv", dtype={"station_id": str})
    ground = pd.read_csv(FDA / "curvas_estaciones_mm.csv.gz", dtype={"STATION_ID": str})
    satellite = pd.read_csv(FDA / "curvas_satelite_mm.csv.gz")
    sat_cols = satellite.columns[3:].tolist()

    rows: list[dict[str, object]] = []
    detail_rows: list[dict[str, object]] = []

    for candidate in audit.itertuples(index=False):
        sid = str(candidate.station_id)
        year = int(candidate.year)
        month = int(candidate.month)
        lat = STATIONS[sid]["lat"]
        lon = STATIONS[sid]["lon"]
        row, col = station_pixel(lat, lon)

        month_sat_cols = sat_cols[(month - 1) * 6 : month * 6]
        center = satellite.loc[
            (satellite["anio"] == year)
            & (satellite["fila"] == row)
            & (satellite["col"] == col),
            month_sat_cols,
        ].iloc[0].to_numpy(dtype=float)

        neighborhood = satellite.loc[
            (satellite["anio"] == year)
            & satellite["fila"].between(row - 1, row + 1)
            & satellite["col"].between(col - 1, col + 1),
            month_sat_cols,
        ]
        neigh_median = neighborhood.median(axis=0).to_numpy(dtype=float)

        g = ground.loc[(ground["STATION_ID"] == sid) & (ground["YEAR"] == year)]
        if g.empty:
            raise RuntimeError(f"No se encontro la curva terrestre {sid}-{year}")
        ground_days = g[month_day_columns(month)].iloc[0].to_numpy(dtype=float)
        ground_blocks = six_ground_blocks(ground_days)

        g_share = shares(ground_blocks)
        c_share = shares(center)
        n_share = shares(neigh_median)
        event_block = min((int(candidate.dominant_day) - 1) // 5, 5)

        historical_concentration = []
        historical_event_block_values = []
        for hist_year in range(1981, 2004):
            hist = satellite.loc[
                (satellite["anio"] == hist_year)
                & (satellite["fila"] == row)
                & (satellite["col"] == col),
                month_sat_cols,
            ]
            if not hist.empty:
                h = hist.iloc[0].to_numpy(dtype=float)
                historical_concentration.append(float(np.nanmax(shares(h))))
                historical_event_block_values.append(float(h[event_block]))
        center_dom = float(np.nanmax(c_share))
        concentration_percentile = float(
            100.0 * np.mean(np.asarray(historical_concentration) <= center_dom)
        )
        event_value_percentile = float(
            100.0
            * np.mean(np.asarray(historical_event_block_values) <= float(center[event_block]))
        )

        sat_peak_block = int(np.nanargmax(n_share))
        event_alignment = abs(sat_peak_block - event_block) <= 1
        spread_blocks = int(np.sum(n_share >= 0.15))
        rows.append(
            {
                "station_id": sid,
                "year": year,
                "month": month,
                "partition": candidate.particion,
                "lat": lat,
                "lon": lon,
                "pixel_row": row,
                "pixel_col": col,
                "ground_dominant_share": float(np.nanmax(g_share)),
                "sat_center_dominant_share": center_dom,
                "sat_neighborhood_dominant_share": float(np.nanmax(n_share)),
                "sat_neighborhood_blocks_ge_15pct": spread_blocks,
                "ground_event_block_1_6": event_block + 1,
                "sat_peak_block_1_6": sat_peak_block + 1,
                "peak_aligned_within_one_block": event_alignment,
                "sat_concentration_percentile_same_month_1981_2003": concentration_percentile,
                "sat_event_block_value_original_scale": float(center[event_block]),
                "sat_event_block_value_percentile_same_month_1981_2003": event_value_percentile,
            }
        )

        for block in range(6):
            detail_rows.append(
                {
                    "station_id": sid,
                    "year": year,
                    "month": month,
                    "block_1_6": block + 1,
                    "ground_mm": ground_blocks[block],
                    "ground_share": g_share[block],
                    "sat_center_original_scale": center[block],
                    "sat_center_share": c_share[block],
                    "sat_neighborhood_median_original_scale": neigh_median[block],
                    "sat_neighborhood_share": n_share[block],
                }
            )

    summary = pd.DataFrame(rows)
    detail = pd.DataFrame(detail_rows)
    summary.to_csv(FDA / "contraste_totales_mensuales_satelite.csv", index=False)
    detail.to_csv(FDA / "detalle_totales_mensuales_satelite.csv", index=False)
    print(summary.to_string(index=False))
    print("\nDetalle de proporciones por bloque:")
    print(
        detail[
            [
                "station_id",
                "year",
                "month",
                "block_1_6",
                "ground_share",
                "sat_neighborhood_share",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()

"""Piloto reproducible satelite -> estacion con cabeza KAN B-spline.

Este archivo ejecuta los calculos que se documentan en el notebook
04_emparejamiento_docente_y_piloto_kan.ipynb. La prueba es metodologica: usa
validacion para comparar modelos y mantiene prueba cerrada.
"""

from __future__ import annotations

import copy
import json
import math
import random
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset


KAN_ROOT = Path(__file__).resolve().parent
ROOT = KAN_ROOT.parent
F2 = KAN_ROOT / "resultados"
OUT = F2 / "piloto_kan_docente"
SEED = 2026
N_SAT_BASIS = 11
N_GROUND_BASIS = 8
DEGREE = 3
N_FUNCTIONAL = 8
KAN_EDGE_BASIS = 10
MAX_EPOCHS = 180
PATIENCE = 25
BATCH_SIZE = 64
LEARNING_RATE = 2e-3


def set_seed(seed: int = SEED) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


def open_uniform_knots(
    n_basis: int, degree: int = 3, low: float = 0.0, high: float = 1.0
) -> np.ndarray:
    n_internal = n_basis - degree - 1
    internal = (
        np.linspace(low, high, n_internal + 2)[1:-1]
        if n_internal > 0
        else np.empty(0)
    )
    return np.r_[np.repeat(low, degree + 1), internal, np.repeat(high, degree + 1)]


def bspline_basis_numpy(
    x: np.ndarray, knots: np.ndarray, degree: int = 3
) -> np.ndarray:
    x = np.asarray(x, dtype=float).reshape(-1)
    upper = np.nextafter(float(knots[-1]), float(knots[0]))
    xe = np.clip(x, float(knots[0]), upper)
    basis = ((xe[:, None] >= knots[:-1]) & (xe[:, None] < knots[1:])).astype(float)
    for p in range(1, degree + 1):
        n_cols = len(knots) - p - 1
        left_den = knots[p : p + n_cols] - knots[:n_cols]
        right_den = knots[p + 1 : p + 1 + n_cols] - knots[1 : 1 + n_cols]
        left = np.zeros((len(xe), n_cols))
        right = np.zeros((len(xe), n_cols))
        valid_left = np.abs(left_den) > 0
        valid_right = np.abs(right_den) > 0
        left[:, valid_left] = (
            (xe[:, None] - knots[:n_cols])[:, valid_left]
            / left_den[valid_left]
            * basis[:, :n_cols][:, valid_left]
        )
        right[:, valid_right] = (
            (knots[p + 1 : p + 1 + n_cols] - xe[:, None])[:, valid_right]
            / right_den[valid_right]
            * basis[:, 1 : n_cols + 1][:, valid_right]
        )
        basis = left + right
    return basis


def curves_to_coefficients(curves: np.ndarray, design: np.ndarray) -> np.ndarray:
    curves = np.asarray(curves, dtype=float)
    if np.isfinite(curves).all():
        lhs = design.T @ design + 1e-8 * np.eye(design.shape[1])
        operator = np.linalg.solve(lhs, design.T)
        return curves @ operator.T
    coefficients = []
    for curve in curves:
        observed = np.isfinite(curve)
        if observed.sum() < design.shape[1]:
            raise RuntimeError("Una curva no tiene observaciones suficientes para la base")
        b = design[observed]
        lhs = b.T @ b + 1e-8 * np.eye(design.shape[1])
        coefficients.append(np.linalg.solve(lhs, b.T @ curve[observed]))
    return np.asarray(coefficients)


def gram_matrix(n_basis: int, degree: int) -> tuple[np.ndarray, np.ndarray]:
    knots = open_uniform_knots(n_basis, degree)
    grid = np.linspace(0, 1, 2001)
    phi = bspline_basis_numpy(grid, knots, degree)
    dx = 1.0 / (len(grid) - 1)
    weights = np.full(len(grid), dx)
    weights[[0, -1]] *= 0.5
    gram = phi.T @ (weights[:, None] * phi)
    return knots, gram


def fit_weighted_fpca(coeffs: np.ndarray, gram: np.ndarray) -> dict[str, np.ndarray]:
    mean = coeffs.mean(axis=0)
    chol = np.linalg.cholesky(gram)
    z = (coeffs - mean) @ chol
    _, singular, vh = np.linalg.svd(z, full_matrices=False)
    basis = vh @ np.linalg.inv(chol)
    signs = np.sign(basis.sum(axis=1, keepdims=True) + 1e-12)
    basis *= signs
    explained = singular**2 / np.sum(singular**2)
    return {"mean": mean, "basis": basis, "explained": explained}


def bspline_basis_torch(x: torch.Tensor, knots: torch.Tensor, degree: int) -> torch.Tensor:
    original = x.shape
    flat = x.reshape(-1)
    eps = 10 * torch.finfo(flat.dtype).eps
    xe = torch.clamp(flat, min=float(knots[0]), max=float(knots[-1]) - eps)
    basis = ((xe[:, None] >= knots[:-1]) & (xe[:, None] < knots[1:])).to(flat.dtype)
    for p in range(1, degree + 1):
        n_cols = len(knots) - p - 1
        left_den = knots[p : p + n_cols] - knots[:n_cols]
        right_den = knots[p + 1 : p + 1 + n_cols] - knots[1 : 1 + n_cols]
        left_num = xe[:, None] - knots[:n_cols]
        right_num = knots[p + 1 : p + 1 + n_cols] - xe[:, None]
        left = torch.where(
            left_den.abs() > 0,
            left_num / left_den.clamp_min(eps),
            torch.zeros_like(left_num),
        ) * basis[:, :n_cols]
        right = torch.where(
            right_den.abs() > 0,
            right_num / right_den.clamp_min(eps),
            torch.zeros_like(right_num),
        ) * basis[:, 1 : n_cols + 1]
        basis = left + right
    return basis.reshape(*original, basis.shape[-1])


class FunctionalEncoder(nn.Module):
    def __init__(self, n_basis: int, n_neurons: int, gram: np.ndarray):
        super().__init__()
        self.weight_coeffs = nn.Parameter(torch.empty(n_neurons, n_basis))
        nn.init.xavier_uniform_(self.weight_coeffs)
        self.register_buffer("gram", torch.tensor(gram, dtype=torch.float32))

    def forward(self, coeffs: torch.Tensor) -> torch.Tensor:
        return coeffs @ self.gram @ self.weight_coeffs.T

    def roughness(self) -> torch.Tensor:
        d2 = (
            self.weight_coeffs[:, 2:]
            - 2 * self.weight_coeffs[:, 1:-1]
            + self.weight_coeffs[:, :-2]
        )
        return d2.square().mean()


class KANBSplineLayer(nn.Module):
    def __init__(
        self,
        n_in: int,
        n_out: int,
        n_basis: int = KAN_EDGE_BASIS,
        degree: int = 3,
        grid_range: tuple[float, float] = (-2.5, 2.5),
    ):
        super().__init__()
        self.degree = degree
        self.norm = nn.LayerNorm(n_in)
        knots = open_uniform_knots(n_basis, degree, *grid_range)
        self.register_buffer("knots", torch.tensor(knots, dtype=torch.float32))
        self.base_weight = nn.Parameter(torch.empty(n_out, n_in))
        self.spline_scale = nn.Parameter(torch.ones(n_out, n_in))
        self.spline_coeff = nn.Parameter(torch.randn(n_out, n_in, n_basis) * 0.03)
        self.bias = nn.Parameter(torch.zeros(n_out))
        nn.init.kaiming_uniform_(self.base_weight, a=math.sqrt(5))

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = self.norm(x)
        base = F.silu(x) @ self.base_weight.T
        spline = torch.einsum(
            "big,oig,oi->bo",
            bspline_basis_torch(x, self.knots, self.degree),
            self.spline_coeff,
            self.spline_scale,
        )
        return base + spline + self.bias

    def spline_l1(self) -> torch.Tensor:
        return self.spline_coeff.abs().mean()


class LinearHead(nn.Module):
    def __init__(self, n_in: int, n_out: int):
        super().__init__()
        self.net = nn.Linear(n_in, n_out)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


class MLPHead(nn.Module):
    def __init__(self, n_in: int, n_out: int, hidden: int):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_in, hidden), nn.GELU(), nn.Linear(hidden, n_out)
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


class KANHead(nn.Module):
    def __init__(self, n_in: int, n_out: int, hidden: int):
        super().__init__()
        self.project = nn.Linear(n_in, hidden)
        self.kan = KANBSplineLayer(hidden, n_out)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.kan(self.project(x))


class FPCADecoder(nn.Module):
    def __init__(self, mean: np.ndarray, basis: np.ndarray):
        super().__init__()
        self.register_buffer("mean", torch.tensor(mean, dtype=torch.float32))
        self.register_buffer("basis", torch.tensor(basis, dtype=torch.float32))

    def forward(self, scores: torch.Tensor) -> torch.Tensor:
        return self.mean + scores @ self.basis


class FunctionalRainNetwork(nn.Module):
    def __init__(
        self,
        satellite_gram: np.ndarray,
        ground_gram: np.ndarray,
        fpca: dict[str, np.ndarray],
        head_kind: str,
        n_coords: int,
        kan_hidden: int,
        mlp_hidden: int,
    ):
        super().__init__()
        self.encoder = FunctionalEncoder(len(satellite_gram), N_FUNCTIONAL, satellite_gram)
        self.register_buffer("output_gram", torch.tensor(ground_gram, dtype=torch.float32))
        n_in = N_FUNCTIONAL + n_coords
        n_out = fpca["basis"].shape[0]
        if head_kind == "linear":
            self.head = LinearHead(n_in, n_out)
        elif head_kind == "mlp":
            self.head = MLPHead(n_in, n_out, mlp_hidden)
        elif head_kind == "kan":
            self.head = KANHead(n_in, n_out, kan_hidden)
        else:
            raise ValueError(head_kind)
        self.decoder = FPCADecoder(fpca["mean"], fpca["basis"])

    def forward(self, sat_coeffs: torch.Tensor, coords: torch.Tensor) -> torch.Tensor:
        z = self.encoder(sat_coeffs)
        features = torch.cat([z, coords], dim=1) if coords.shape[1] else z
        return self.decoder(self.head(features))


def count_params(module: nn.Module) -> int:
    return sum(p.numel() for p in module.parameters() if p.requires_grad)


def find_matched_widths(n_in: int, n_out: int) -> tuple[int, int, int, int, float]:
    best = None
    for hk in range(4, 33):
        pk = count_params(KANHead(n_in, n_out, hk))
        for hm in range(4, 257):
            pm = count_params(MLPHead(n_in, n_out, hm))
            diff = 100 * abs(pk - pm) / max(pk, pm)
            candidate = (diff, pk, hk, hm, pm)
            if best is None or candidate < best:
                best = candidate
    diff, pk, hk, hm, pm = best
    return hk, hm, pk, pm, diff


def functional_l2(pred: torch.Tensor, true: torch.Tensor, gram: torch.Tensor) -> torch.Tensor:
    d = pred - true
    return torch.einsum("bi,ij,bj->b", d, gram, d).mean()


def train_model(
    name: str,
    model: FunctionalRainNetwork,
    train_loader: DataLoader,
    val_tensors: tuple[torch.Tensor, torch.Tensor, torch.Tensor],
) -> tuple[FunctionalRainNetwork, pd.DataFrame]:
    optimizer = torch.optim.Adam(model.parameters(), lr=LEARNING_RATE, weight_decay=1e-5)
    sat_val, coord_val, y_val = val_tensors
    best_loss = float("inf")
    best_state = None
    wait = 0
    rows = []
    for epoch in range(1, MAX_EPOCHS + 1):
        model.train()
        train_sum = 0.0
        n_seen = 0
        for sat_b, coord_b, y_b in train_loader:
            optimizer.zero_grad()
            pred = model(sat_b, coord_b)
            loss = functional_l2(pred, y_b, model.output_gram)
            loss = loss + 1e-3 * model.encoder.roughness()
            if hasattr(model.head, "kan"):
                loss = loss + 1e-6 * model.head.kan.spline_l1()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
            train_sum += float(loss.detach()) * len(sat_b)
            n_seen += len(sat_b)
        model.eval()
        with torch.no_grad():
            val_loss = float(functional_l2(model(sat_val, coord_val), y_val, model.output_gram))
        rows.append(
            {"modelo": name, "epoca": epoch, "train_loss": train_sum / n_seen, "val_loss": val_loss}
        )
        if val_loss < best_loss - 1e-7:
            best_loss = val_loss
            best_state = copy.deepcopy(model.state_dict())
            wait = 0
        else:
            wait += 1
            if wait >= PATIENCE:
                break
    model.load_state_dict(best_state)
    return model, pd.DataFrame(rows)


def inverse_target(coeffs: np.ndarray, ground_design: np.ndarray) -> np.ndarray:
    sqrt_values = coeffs @ ground_design.T
    return np.clip(sqrt_values, 0, None) ** 2


def calculate_metrics(
    observed_mm: np.ndarray, target_curve_mm: np.ndarray, pred_mm: np.ndarray
) -> dict[str, float]:
    observed_error = pred_mm - observed_mm
    curve_error = pred_mm - target_curve_mm
    true_total = target_curve_mm.sum(axis=1)
    pred_total = pred_mm.sum(axis=1)
    return {
        "MAE_observado_mm": float(np.nanmean(np.abs(observed_error))),
        "RMSE_observado_mm": float(np.sqrt(np.nanmean(observed_error**2))),
        "RMSE_objetivo_funcional_mm": float(np.sqrt(np.mean(curve_error**2))),
        "sesgo_diario_observado_mm": float(np.nanmean(observed_error)),
        "sesgo_total_anual_mm": float(np.mean(pred_total - true_total)),
        "MAE_total_anual_mm": float(np.mean(np.abs(pred_total - true_total))),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    set_seed()
    torch.set_num_threads(1)

    sat = pd.read_csv(
        F2 / "entradas_satelitales_72_docente.csv.gz", dtype={"station_id": "string"}
    )
    ground = pd.read_csv(
        F2 / "salidas_terrestres_365_docente.csv.gz", dtype={"station_id": "string"}
    )
    key = ["station_id", "year"]
    if not sat[key].equals(ground[key]):
        raise RuntimeError("Las filas satelitales y terrestres no estan alineadas")

    sat_cols = [c for c in sat.columns if c.startswith("t")]
    ground_cols = [str(i) for i in range(1, 366)]
    partition = sat["particion"].to_numpy()
    train_mask = partition == "ENTRENAMIENTO"
    val_mask = partition == "VALIDACION"

    sat_raw = sat[sat_cols].to_numpy(float)
    ground_mm = ground[ground_cols].to_numpy(float)
    sat_mean = float(sat_raw[train_mask].mean())
    sat_std = float(sat_raw[train_mask].std())
    sat_scaled = (sat_raw - sat_mean) / sat_std
    ground_sqrt = np.sqrt(np.clip(ground_mm, 0, None))

    sat_knots, sat_gram = gram_matrix(N_SAT_BASIS, DEGREE)
    ground_knots, ground_gram = gram_matrix(N_GROUND_BASIS, DEGREE)
    sat_time = np.array([float(c[1:]) for c in sat_cols]) / 365.0
    ground_time = np.arange(1, 366) / 365.0
    sat_design = bspline_basis_numpy(sat_time, sat_knots, DEGREE)
    ground_design = bspline_basis_numpy(ground_time, ground_knots, DEGREE)
    sat_coeffs = curves_to_coefficients(sat_scaled, sat_design)
    ground_coeffs = curves_to_coefficients(ground_sqrt, ground_design)

    fpca_full = fit_weighted_fpca(ground_coeffs[train_mask], ground_gram)
    cumulative = np.cumsum(fpca_full["explained"])
    eligible_k = [k for k in range(4, 7) if cumulative[k - 1] >= 0.95]
    n_fpca = eligible_k[0] if eligible_k else 6
    fpca = {
        "mean": fpca_full["mean"],
        "basis": fpca_full["basis"][:n_fpca],
        "explained": fpca_full["explained"],
    }

    coords = sat[["latitude", "longitude", "elevation_m"]].to_numpy(float)
    coord_mean = coords[train_mask].mean(axis=0)
    coord_std = coords[train_mask].std(axis=0)
    coords_scaled = (coords - coord_mean) / coord_std

    x_sat = torch.tensor(sat_coeffs, dtype=torch.float32)
    x_geo = torch.tensor(coords_scaled, dtype=torch.float32)
    x_latlon = x_geo[:, :2]
    x_empty = torch.empty((len(sat), 0), dtype=torch.float32)
    y = torch.tensor(ground_coeffs, dtype=torch.float32)
    train_index = np.flatnonzero(train_mask)
    val_index = np.flatnonzero(val_mask)

    hk, hm, pk, pm, param_diff = find_matched_widths(
        N_FUNCTIONAL + 3, n_fpca
    )
    configs = [
        ("Lineal + lat/lon/elevación", "linear", "geo3"),
        ("MLP + lat/lon/elevación", "mlp", "geo3"),
        ("KAN + lat/lon/elevación", "kan", "geo3"),
        ("KAN + lat/lon", "kan", "latlon"),
        ("KAN sin ubicación", "kan", "none"),
    ]
    models = {}
    histories = []
    metrics = []
    predictions = {}

    for name, head, coordinate_kind in configs:
        # Misma inicializacion aleatoria y mismo orden de lotes para que la
        # comparacion entre cabezas no se confunda con una semilla diferente.
        set_seed(SEED)
        if coordinate_kind == "geo3":
            coord_tensor = x_geo
        elif coordinate_kind == "latlon":
            coord_tensor = x_latlon
        else:
            coord_tensor = x_empty
        model = FunctionalRainNetwork(
            sat_gram, ground_gram, fpca, head, coord_tensor.shape[1], hk, hm
        )
        dataset = TensorDataset(
            x_sat[train_index], coord_tensor[train_index], y[train_index]
        )
        generator = torch.Generator().manual_seed(SEED)
        loader = DataLoader(
            dataset, batch_size=BATCH_SIZE, shuffle=True, generator=generator
        )
        model, history = train_model(
            name,
            model,
            loader,
            (x_sat[val_index], coord_tensor[val_index], y[val_index]),
        )
        model.eval()
        with torch.no_grad():
            pred_coeff = model(x_sat[val_index], coord_tensor[val_index]).numpy()
        pred_mm = inverse_target(pred_coeff, ground_design)
        true_mm = ground_mm[val_index]
        target_curve_mm = inverse_target(ground_coeffs[val_index], ground_design)
        result = calculate_metrics(true_mm, target_curve_mm, pred_mm)
        result.update(
            {
                "modelo": name,
                "parametros_totales": count_params(model),
                "parametros_cabeza": count_params(model.head),
                "epocas": len(history),
                "mejor_val_L2": float(history["val_loss"].min()),
            }
        )
        metrics.append(result)
        histories.append(history)
        models[name] = model
        predictions[name] = pred_mm

    metrics_df = pd.DataFrame(metrics).sort_values("RMSE_observado_mm").reset_index(drop=True)
    history_df = pd.concat(histories, ignore_index=True)
    metrics_df.to_csv(OUT / "metricas_validacion.csv", index=False)
    history_df.to_csv(OUT / "historial_entrenamiento.csv", index=False)

    fpca_df = pd.DataFrame(
        {
            "componente": np.arange(1, len(fpca_full["explained"]) + 1),
            "varianza_explicada": fpca_full["explained"],
            "varianza_acumulada": cumulative,
        }
    )
    fpca_df.to_csv(OUT / "varianza_fpca_train.csv", index=False)

    fig, ax = plt.subplots(figsize=(10, 5))
    for name, group in history_df.groupby("modelo", sort=False):
        ax.plot(group["epoca"], group["val_loss"], label=name)
    ax.set(
        title="Pérdida funcional en validación",
        xlabel="Época",
        ylabel="L2 funcional en espacio raíz cuadrada",
    )
    ax.set_yscale("log")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT / "curvas_validacion.png", dpi=160)
    plt.close(fig)

    best_name = metrics_df.iloc[0]["modelo"]
    best_pred = predictions[best_name]
    observed_validation = ground_mm[val_index]
    functional_validation = inverse_target(ground_coeffs[val_index], ground_design)
    curve_rmse = np.sqrt(np.nanmean((best_pred - observed_validation) ** 2, axis=1))
    order = np.argsort(curve_rmse)
    chosen = [order[int(q * (len(order) - 1))] for q in (0.25, 0.50, 0.75)]
    days = np.arange(1, 366)
    fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
    val_meta = sat.iloc[val_index].reset_index(drop=True)
    for ax, idx in zip(axes, chosen):
        ax.scatter(
            days,
            observed_validation[idx],
            color="black",
            s=5,
            alpha=0.45,
            label="Terrestre observado",
        )
        ax.plot(
            days,
            functional_validation[idx],
            color="#64748b",
            lw=1.2,
            label="Objetivo funcional",
        )
        ax.plot(days, best_pred[idx], color="#7b2cbf", lw=1.5, label=best_name)
        row = val_meta.iloc[idx]
        ax.set_title(f"{row['station_name']} · {int(row['year'])} · RMSE={curve_rmse[idx]:.2f} mm")
        ax.set_ylabel("mm/día")
        ax.grid(alpha=0.2)
        ax.legend(fontsize=8)
    axes[-1].set_xlabel("Día del año")
    fig.tight_layout()
    fig.savefig(OUT / "ejemplos_reconstruccion_validacion.png", dpi=160)
    plt.close(fig)

    stations = sat.drop_duplicates("station_id")
    pixels = stations.drop_duplicates(["pixel_row", "pixel_col"])
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(
        pixels["pixel_center_lon"],
        pixels["pixel_center_lat"],
        s=18,
        color="#f59e0b",
        label="Centros de píxel",
        alpha=0.75,
    )
    ax.scatter(
        stations["longitude"],
        stations["latitude"],
        s=12,
        color="#2563eb",
        label="Estaciones",
        alpha=0.75,
    )
    ax.set(
        title="Estaciones emparejadas y píxeles satelitales",
        xlabel="Longitud",
        ylabel="Latitud",
    )
    ax.legend()
    ax.grid(alpha=0.2)
    ax.set_aspect("equal", adjustable="box")
    fig.tight_layout()
    fig.savefig(OUT / "mapa_emparejamiento.png", dpi=160)
    plt.close(fig)

    summary = {
        "estado": "piloto ejecutado; prueba no consultada",
        "pares_train": int(train_mask.sum()),
        "pares_validacion": int(val_mask.sum()),
        "pares_prueba_reservados": int((partition == "PRUEBA").sum()),
        "estaciones": int(sat["station_id"].nunique()),
        "pixeles": int(sat[["pixel_row", "pixel_col"]].drop_duplicates().shape[0]),
        "bases_satelite": N_SAT_BASIS,
        "bases_terrestre": N_GROUND_BASIS,
        "componentes_fpca": int(n_fpca),
        "varianza_fpca_acumulada": float(cumulative[n_fpca - 1]),
        "media_satelital_train_escala_original": sat_mean,
        "desviacion_satelital_train_escala_original": sat_std,
        "covariables_media_train_lat_lon_elevacion": coord_mean.tolist(),
        "covariables_desviacion_train_lat_lon_elevacion": coord_std.tolist(),
        "kan_hidden": int(hk),
        "mlp_hidden": int(hm),
        "parametros_cabeza_kan": int(pk),
        "parametros_cabeza_mlp": int(pm),
        "diferencia_parametros_pct": float(param_diff),
        "mejor_modelo_validacion": str(best_name),
        "limitaciones": [
            "la escala fisica y el posible techo 200 del satelite siguen pendientes",
            "las estaciones repetidas entre anos permiten evaluar transferencia temporal, no espacial",
            "una sola semilla no permite declarar superioridad entre arquitecturas",
        ],
    }
    (OUT / "resumen_piloto.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print("\nMétricas de validación:\n", metrics_df.to_string(index=False))


if __name__ == "__main__":
    main()

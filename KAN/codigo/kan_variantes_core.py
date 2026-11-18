"""Utilidades compartidas para experimentos KAN que solo usan entrenamiento y validación."""

from __future__ import annotations

import json
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

torch.set_default_dtype(torch.float64)
torch.set_num_threads(4)
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "salidas_R_funcional" / "para_python"
SEEDS = (2026, 2027, 2028)
N_ENCODER = 8
LAMBDA_W = 1e6
MAX_EPOCHS, PATIENCE, LR = 5000, 1000, 1e-2
DAYS = np.arange(1, 366, dtype=float)


def _read_problem():
    spec = json.loads((INPUT / "especificacion.json").read_text(encoding="utf-8"))
    pairs = pd.read_csv(INPUT / "pares.csv", dtype={"codigo": str})
    read_matrix = lambda name, **kwargs: pd.read_csv(INPUT / name, **kwargs).to_numpy(float)
    expected = {"entrenamiento": 2257, "validación": 672, "prueba": 370}
    counts = pairs["particion"].value_counts().to_dict()
    if len(pairs) != 3299 or counts != expected or not pairs["usar"].all():
        raise ValueError(f"Partición/insumos inesperados: {counts}")
    if not (pairs.loc[pairs.particion == "entrenamiento", "anio"].max() <= 1995
            and pairs.loc[pairs.particion == "validación", "anio"].between(1996, 2000).all()):
        raise ValueError("La partición temporal de entrenamiento/validación no coincide con el proyecto")
    train_rows = np.flatnonzero(pairs.particion.to_numpy() == "entrenamiento")
    valid_rows = np.flatnonzero(pairs.particion.to_numpy() == "validación")
    test_rows = np.flatnonzero(pairs.particion.to_numpy() == "prueba")
    eligible_rows = np.concatenate([train_rows, valid_rows])
    eligible_set = set(eligible_rows.tolist())
    n_train = len(train_rows)
    tau_e = np.asarray(spec["nodos_estaciones"], dtype=float)
    tau_s = np.asarray(spec["nodos_satelite"], dtype=float)
    c_sat = read_matrix("coef_satelite_raiz.csv")
    c_station = read_matrix(
        "coef_estacion_raiz.csv",
        skiprows=lambda row: row > 0 and row - 1 not in eligible_set,
    )
    c_station = c_station[np.searchsorted(np.sort(eligible_rows), eligible_rows)]
    w_e_r = read_matrix("W_e.csv", header=None)
    w_s_r = read_matrix("W_s.csv", header=None)
    r_s = read_matrix("R_s.csv", header=None)
    fpca = pd.read_csv(INPUT / "fpca_estaciones.csv")
    mu = fpca["mu"].to_numpy(float)
    basis_fpca = fpca.drop(columns="mu").to_numpy(float)
    if c_sat.shape[0] != len(pairs) or c_station.shape[0] != len(eligible_rows):
        raise ValueError("Las matrices de coeficientes no corresponden a las particiones")
    eligible_pairs = pairs.iloc[eligible_rows].reset_index(drop=True)
    test_pairs = pairs.iloc[test_rows].reset_index(drop=True)
    altitude = ((eligible_pairs["altitud"].to_numpy(float) - spec["altitud_media"]) / spec["altitud_de"])
    enso = eligible_pairs["ONI_anual"].to_numpy(float)
    test_altitude = ((test_pairs["altitud"].to_numpy(float) - spec["altitud_media"]) / spec["altitud_de"])
    test_enso = test_pairs["ONI_anual"].to_numpy(float)
    return {
        "spec": spec, "valid_pairs": pairs.iloc[valid_rows].reset_index(drop=True),
        "valid_global_rows": valid_rows, "test_pairs": test_pairs, "test_global_rows": test_rows,
        "test_c_sat": c_sat[test_rows],
        "test_extra": np.column_stack([test_altitude, test_enso]),
        "tau_e": tau_e, "tau_s": tau_s,
        "c_sat": c_sat[eligible_rows], "c_station": c_station,
        "mu": mu, "basis_fpca": basis_fpca,
        "w_e_r": w_e_r, "w_s_r": w_s_r, "w_s": w_s_r, "r_s": r_s,
        "altitude": altitude, "enso": enso,
        "train": np.arange(n_train),
        "valid": np.arange(n_train, len(eligible_rows)),
    }


DATA = _read_problem()
K_E = DATA["c_station"].shape[1]
K_S = DATA["c_sat"].shape[1]
M = DATA["basis_fpca"].shape[1]


def _basis_at_nodes(t, tau, order=4):
    t = np.asarray(t, dtype=float)
    basis = np.zeros((len(t), len(tau) - 1))
    for j in range(len(tau) - 1):
        basis[:, j] = (t >= tau[j]) & (t < tau[j + 1])
    basis[t == tau[-1], np.max(np.flatnonzero(tau[:-1] < tau[1:]))] = 1
    for degree in range(2, order + 1):
        next_basis = np.zeros((len(t), len(tau) - degree))
        for j in range(next_basis.shape[1]):
            d1 = tau[j + degree - 1] - tau[j]
            d2 = tau[j + degree] - tau[j + 1]
            if d1 > 0:
                next_basis[:, j] += (t - tau[j]) / d1 * basis[:, j]
            if d2 > 0:
                next_basis[:, j] += (tau[j + degree] - t) / d2 * basis[:, j + 1]
        basis = next_basis
    return basis


def _gram(tau, n_gl=4):
    unique = np.unique(tau)
    x, weights = np.polynomial.legendre.leggauss(n_gl)
    left, right = unique[:-1, None], unique[1:, None]
    tq = ((right - left) / 2 * x + (left + right) / 2).ravel()
    wq = ((right - left) / 2 * weights).ravel()
    basis = _basis_at_nodes(tq, tau)
    return basis.T @ (wq[:, None] * basis)


def _tensor(array):
    return torch.as_tensor(np.array(array, dtype=float, copy=True), device=DEVICE)


X_T = _tensor(DATA["c_sat"])
Y_T = _tensor(DATA["c_station"])
EXTRA_T = _tensor(np.column_stack([DATA["altitude"], DATA["enso"]]))
W_E_T = _tensor(_gram(DATA["tau_e"]))
W_S_T = _tensor(_gram(DATA["tau_s"]))
R_S_T = _tensor(DATA["r_s"])


class KANLayer(torch.nn.Module):
    def __init__(self, n_in: int, n_out: int, grid_intervals: int):
        super().__init__()
        if grid_intervals < 1:
            raise ValueError("grid_intervals debe ser positivo")
        self.grid_intervals = grid_intervals
        self.order = 3
        self.range = (-3.0, 3.0)
        step = (self.range[1] - self.range[0]) / grid_intervals
        self.register_buffer(
            "grid",
            torch.arange(-self.order, grid_intervals + self.order + 1, dtype=torch.float64) * step + self.range[0],
        )
        self.coef = torch.nn.Parameter(0.1 * torch.randn(n_in, n_out, grid_intervals + self.order))
        self.w_b = torch.nn.Parameter(torch.randn(n_in, n_out) / np.sqrt(n_in))
        self.w_s = torch.nn.Parameter(torch.ones(n_in, n_out))

    def bases(self, values):
        grid, values = self.grid, values.unsqueeze(-1)
        basis = ((values >= grid[:-1]) & (values < grid[1:])).to(values.dtype)
        for degree in range(1, self.order + 1):
            basis = (
                (values - grid[:-(degree + 1)]) / (grid[degree:-1] - grid[:-(degree + 1)]) * basis[..., :-1]
                + (grid[degree + 1:] - values) / (grid[degree + 1:] - grid[1:-degree]) * basis[..., 1:]
            )
        return basis

    def activations(self, values):
        spline = torch.einsum("nig,iog->nio", self.bases(values), self.coef)
        return self.w_b * torch.nn.functional.silu(values).unsqueeze(-1) + self.w_s * spline

    def forward(self, values):
        return self.activations(values).sum(1)


class KAN(torch.nn.Module):
    def __init__(self, widths: tuple[int, ...], grid_intervals: int):
        super().__init__()
        if len(widths) < 3 or widths[0] != N_ENCODER + 2 or widths[-1] != M:
            raise ValueError(f"Arquitectura KAN inválida: {widths}")
        self.layers = torch.nn.ModuleList(
            KANLayer(n_in, n_out, grid_intervals) for n_in, n_out in zip(widths[:-1], widths[1:])
        )

    def forward(self, values):
        for layer in self.layers:
            values = layer(values)
        return values

    def second_difference_penalty(self):
        penalties = []
        for layer in self.layers:
            effective_coef = layer.w_s.unsqueeze(-1) * layer.coef
            second_difference = torch.diff(effective_coef, n=2, dim=-1)
            penalties.append(second_difference.square().mean())
        return torch.stack(penalties).mean()

    def n_parameters(self):
        return int(sum(parameter.numel() for parameter in self.parameters()))


class FunctionalKAN(torch.nn.Module):
    def __init__(self, widths: tuple[int, ...], grid_intervals: int):
        super().__init__()
        self.encoder = torch.nn.Parameter(0.1 * torch.randn(N_ENCODER, K_S))
        self.head = KAN(widths, grid_intervals)
        self.register_buffer("mu", _tensor(DATA["mu"]))
        self.register_buffer("basis_fpca", _tensor(DATA["basis_fpca"]))

    def latent(self, c_sat):
        return c_sat @ W_S_T @ self.encoder.T

    def head_input(self, c_sat, extra):
        train_idx = torch.as_tensor(DATA["train"], device=DEVICE)
        train_latent = self.latent(X_T[train_idx])
        center, scale = train_latent.mean(0), train_latent.std(0)
        return torch.cat([(self.latent(c_sat) - center) / scale, extra], dim=1)

    def forward(self, c_sat, extra):
        scores = self.head(self.head_input(c_sat, extra))
        return self.mu + scores @ self.basis_fpca.T

    def encoder_roughness(self, r_s, w_s):
        return (
            torch.einsum("ki,ij,kj->k", self.encoder, r_s, self.encoder)
            / torch.einsum("ki,ij,kj->k", self.encoder, w_s, self.encoder)
        ).sum()


def _indices(name):
    return torch.as_tensor(DATA[name], dtype=torch.long, device=DEVICE)


def validation_loss(model, indices):
    c_sat = X_T[indices]
    extra = EXTRA_T[indices]
    target = Y_T[indices]
    difference = model(c_sat, extra) - target
    return torch.einsum("ni,ij,nj->n", difference, W_E_T, difference).mean()


def fit(
    widths: tuple[int, ...],
    grid_intervals: int,
    seed: int,
    lambda_second_difference: float = 0.0,
    *,
    lambda_w_encoder: float = LAMBDA_W,
    learning_rate: float = LR,
    max_epochs: int = MAX_EPOCHS,
    patience: int = PATIENCE,
):
    if lambda_second_difference < 0:
        raise ValueError("lambda_second_difference no puede ser negativo")
    if lambda_w_encoder < 0:
        raise ValueError("lambda_w_encoder no puede ser negativo")
    if learning_rate <= 0:
        raise ValueError("learning_rate debe ser positivo")
    if max_epochs < 1:
        raise ValueError("max_epochs debe ser al menos 1")
    if patience < 0:
        raise ValueError("patience no puede ser negativo")
    torch.manual_seed(seed)
    head = KAN(widths, grid_intervals)
    torch.manual_seed(seed)
    model = FunctionalKAN(widths, grid_intervals)
    model.head.load_state_dict(head.state_dict())
    model = model.to(DEVICE)
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    train_idx, valid_idx = _indices("train"), _indices("valid")
    x_train = X_T[train_idx]
    extra_all = EXTRA_T
    y_train = Y_T[train_idx]

    best_loss, best_epoch, best_state = np.inf, 0, None
    started = time.time()
    epoch = -1
    for epoch in range(max_epochs):
        model.train()
        optimizer.zero_grad()
        difference = model(x_train, extra_all[train_idx]) - y_train
        data_loss = torch.einsum("ni,ij,nj->n", difference, W_E_T, difference).mean()
        objective = data_loss + lambda_w_encoder * model.encoder_roughness(R_S_T, W_S_T)
        if lambda_second_difference:
            objective = objective + lambda_second_difference * model.head.second_difference_penalty()
        objective.backward()
        optimizer.step()
        model.eval()
        with torch.no_grad():
            current = float(validation_loss(model, valid_idx).item())
        if current < best_loss - 1e-10:
            best_loss, best_epoch = current, epoch
            best_state = {key: value.detach().clone() for key, value in model.state_dict().items()}
        elif epoch - best_epoch > patience:
            break
    if best_state is None or epoch < 0:
        raise RuntimeError("Entrenamiento sin un estado seleccionado por validación")
    model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        difference = model(x_train, extra_all[train_idx]) - y_train
        train_loss = torch.einsum("ni,ij,nj->n", difference, W_E_T, difference).mean().item()
    info = {
        "semilla": seed,
        "lambda_w_encoder": lambda_w_encoder,
        "lambda_segunda_diferencia": lambda_second_difference,
        "learning_rate": learning_rate,
        "max_epochs": max_epochs,
        "patience": patience,
        "epoca_optima": best_epoch + 1,
        "epocas_corridas": epoch + 1,
        "L2_entrenamiento": train_loss,
        "L2_validacion": best_loss,
        "segundos": round(time.time() - started, 1),
        "parametros_cabeza": model.head.n_parameters(),
        "parametros_totales": model.head.n_parameters() + N_ENCODER * K_S,
    }
    return model, info


def predict_validation(model):
    indices = _indices("valid")
    extra = EXTRA_T[indices]
    with torch.no_grad():
        return model(X_T[indices], extra).cpu().numpy()


def predict_test(model):
    test_x = _tensor(DATA["test_c_sat"])
    test_extra = _tensor(DATA["test_extra"])
    with torch.no_grad():
        return model(test_x, test_extra).cpu().numpy()


def save_test_run(output_dir: Path, rows: list[dict], predictions: list[dict], config: dict):
    output_dir.mkdir(parents=True, exist_ok=True)
    training = pd.DataFrame(rows)
    training.to_csv(output_dir / "entrenamientos_validacion.csv", index=False)
    frames = []
    pairs = DATA["test_pairs"]
    for model_name, seed, coefficients in predictions:
        frame = pd.DataFrame(coefficients, columns=[f"c{i + 1}" for i in range(K_E)])
        frame.insert(0, "semilla", seed)
        frame.insert(0, "modelo", model_name)
        frame.insert(0, "particion", "prueba")
        frame.insert(0, "anio", pairs["anio"].to_numpy())
        frame.insert(0, "codigo", pairs["codigo"].to_numpy())
        frame.insert(0, "fila", DATA["test_global_rows"] + 1)
        frames.append(frame)
    pred = pd.concat(frames, ignore_index=True)
    pred.to_csv(output_dir / "predicciones_prueba_raiz.csv", index=False, float_format="%.12g")
    checks = {
        "Todas las predicciones pertenecen a prueba": pred["particion"].eq("prueba").all(),
        "Cada corrida exporta 370 curvas": pred.groupby(["modelo", "semilla"]).size().eq(len(pairs)).all(),
        "Predicciones y pérdidas de validación son finitas":
            np.isfinite(pred.filter(regex=r"^c\d+$").to_numpy()).all()
            and np.isfinite(training[["L2_entrenamiento", "L2_validacion"]].to_numpy()).all(),
        "Tres semillas y épocas seleccionadas en validación":
            training["semilla"].nunique() == len(SEEDS) and training["epoca_optima"].gt(0).all(),
    }
    pd.DataFrame([{"prueba": name, "pasa": bool(ok)} for name, ok in checks.items()]).to_csv(
        output_dir / "P6_pruebas.csv", index=False
    )
    if not all(checks.values()):
        raise RuntimeError(f"Fallaron controles de salida de prueba: {checks}")
    config = {
        **config,
        "particion_predicha": "prueba (2001-2003)",
        "filas_prueba": int(len(pairs)),
        "semillas": SEEDS,
        "lambda_w_encoder_fijo": LAMBDA_W,
        "lambda_segunda_diferencia_fijo": 100.0,
        "seleccion_epocas": "solo L2 de validacion",
        "seleccion_hyperparametros_con_prueba": False,
        "evaluacion_de_prueba": "exploratoria: el conjunto ya se habia consultado para modelos base",
        "dispositivo": str(DEVICE),
    }
    (output_dir / "configuracion.json").write_text(
        json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    return pd.DataFrame(rows)


def save_run(output_dir: Path, config: dict, rows: list[dict], predictions: list[dict]):
    output_dir.mkdir(parents=True, exist_ok=True)
    training = pd.DataFrame(rows)
    training.to_csv(output_dir / "entrenamientos_validacion.csv", index=False)
    pred_frames = []
    pairs = DATA["valid_pairs"]
    for item in predictions:
        model_name, seed, coefficients = item
        frame = pd.DataFrame(coefficients, columns=[f"c{i + 1}" for i in range(K_E)])
        frame.insert(0, "semilla", seed)
        frame.insert(0, "modelo", model_name)
        frame.insert(0, "particion", "validación")
        frame.insert(0, "anio", pairs["anio"].to_numpy())
        frame.insert(0, "codigo", pairs["codigo"].to_numpy())
        frame.insert(0, "fila", DATA["valid_global_rows"] + 1)
        pred_frames.append(frame)
    pred = pd.concat(pred_frames, ignore_index=True)
    pred.to_csv(output_dir / "predicciones_validacion_raiz.csv", index=False, float_format="%.12g")
    checks = {
        "Todas las predicciones pertenecen a validación": pred["particion"].eq("validación").all(),
        "Cada corrida exporta 672 curvas": pred.groupby(["modelo", "semilla"]).size().eq(len(DATA["valid"])).all(),
        "Las predicciones y pérdidas son finitas": np.isfinite(pred.filter(regex=r"^c\d+$").to_numpy()).all()
        and np.isfinite(training[["L2_entrenamiento", "L2_validacion"]].to_numpy()).all(),
        "Tres semillas por variante": training.groupby("modelo")["semilla"].nunique().eq(len(SEEDS)).all(),
    }
    pd.DataFrame([{"prueba": name, "pasa": bool(ok)} for name, ok in checks.items()]).to_csv(
        output_dir / "P6_pruebas.csv", index=False
    )
    if not all(checks.values()):
        raise RuntimeError(f"Fallaron controles de salida: {checks}")
    config = {
        **config,
        "particion_evaluada": "validación (1996-2000)",
        "filas_validacion": int(len(DATA["valid"])),
        "semillas": SEEDS,
        "lambda_w_encoder_fijo": LAMBDA_W,
        "prueba_del_proyecto_evaluada": False,
        "dispositivo": str(DEVICE),
    }
    (output_dir / "configuracion.json").write_text(
        json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    return pd.DataFrame(rows)


def output_path(env_name: str, folder_name: str) -> Path:
    return Path(os.environ.get(env_name, str(ROOT / "KAN" / "salidas" / folder_name)))


def run_grid_sensitivity():
    widths = (N_ENCODER + 2, 4, M)
    rows, predictions = [], []
    for intervals in (3, 5, 7):
        for seed in SEEDS:
            model, info = fit(widths, intervals, seed)
            name = f"KAN G={intervals}"
            rows.append({"modelo": name, "G": intervals, "arquitectura": str(widths), **info})
            predictions.append((name, seed, predict_validation(model)))
            print(name, seed, f"L2_validacion={info['L2_validacion']:.6f}", flush=True)
            del model
    return save_run(
        output_path("LLUVIA_SALIDA_KAN_REJILLA", "sensibilidad_rejilla"),
        {
            "objetivo": "sensibilidad al numero de intervalos de la rejilla KAN",
            "arquitectura": widths,
            "valores_G": (3, 5, 7),
            "seleccion": "comparacion de L2 funcional medio de validacion entre tres semillas",
        },
        rows,
        predictions,
    )


def run_depth_comparison():
    architectures = {
        "KAN base (2 capas)": (N_ENCODER + 2, 4, M),
        "KAN profunda (3 capas)": (N_ENCODER + 2, 4, 3, M),
    }
    rows, predictions = [], []
    for name, widths in architectures.items():
        for seed in SEEDS:
            model, info = fit(widths, 5, seed)
            rows.append({"modelo": name, "G": 5, "arquitectura": str(widths), **info})
            predictions.append((name, seed, predict_validation(model)))
            print(name, seed, f"L2_validacion={info['L2_validacion']:.6f}", flush=True)
            del model
    return save_run(
        output_path("LLUVIA_SALIDA_KAN_PROFUNDIDAD", "comparacion_profundidad"),
        {
            "objetivo": "comparar profundidad con el mismo G y lambda_w",
            "arquitecturas": architectures,
            "capas_kan": {"KAN base (2 capas)": 2, "KAN profunda (3 capas)": 3},
            "seleccion": "comparacion descriptiva de validacion; no se usa el conjunto de prueba",
        },
        rows,
        predictions,
    )


def run_second_difference():
    widths = (N_ENCODER + 2, 4, M)
    candidates = (0.0, 0.1, 1.0, 10.0, 100.0)
    search_rows = []
    for penalty in candidates:
        model, info = fit(widths, 5, SEEDS[0], penalty)
        search_rows.append({"lambda_segunda_diferencia": penalty, **info})
        print(
            "D2",
            penalty,
            f"L2_validacion={info['L2_validacion']:.6f}",
            f"epoca={info['epoca_optima']}",
            flush=True,
        )
        del model
    selected = min(search_rows, key=lambda row: row["L2_validacion"])["lambda_segunda_diferencia"]
    rows, predictions = [], []
    configurations = [("KAN base (sin D2)", 0.0)]
    if selected != 0:
        configurations.append((f"KAN + D2 (lambda={selected:g})", selected))
    else:
        configurations = [("KAN base = seleccion D2", 0.0)]
    for name, penalty in configurations:
        for seed in SEEDS:
            model, info = fit(widths, 5, seed, penalty)
            rows.append({"modelo": name, "G": 5, "arquitectura": str(widths), **info})
            predictions.append((name, seed, predict_validation(model)))
            print(name, seed, f"L2_validacion={info['L2_validacion']:.6f}", flush=True)
            del model
    output_dir = output_path("LLUVIA_SALIDA_KAN_RUGOSIDAD", "regularizacion_segunda_diferencia")
    frame = save_run(
        output_dir,
        {
            "objetivo": "penalizar segundas diferencias de los coeficientes spline efectivos w_s * c",
            "arquitectura": widths,
            "G": 5,
            "candidatos_lambda": candidates,
            "lambda_seleccionado": selected,
            "seleccion": "menor L2 funcional de validacion con semilla 2026",
            "definicion_penalizacion": "media por capa de mean((diff(w_s * coef, n=2))**2)",
        },
        rows,
        predictions,
    )
    pd.DataFrame(search_rows).to_csv(output_dir / "busqueda_lambda_D2_validacion.csv", index=False)
    return frame


print(
    f"KAN variants core: {DEVICE}; train={len(DATA['train'])}, validation={len(DATA['valid'])}; "
    "test evaluation is an explicit post-fit step",
    flush=True,
)

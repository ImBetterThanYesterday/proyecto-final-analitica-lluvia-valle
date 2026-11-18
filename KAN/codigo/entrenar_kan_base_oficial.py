"""Entrena la KAN base definida por el enunciado del proyecto.

La cabeza contiene exclusivamente la parametrizacion solicitada:
    phi_ij(x) = w_b * SiLU(x) + w_s * sum_g c_g B_g(x)

No se aplica regularizacion adicional a los coeficientes de la cabeza KAN.
La unica penalizacion de rugosidad corresponde a las funciones w_k(t) del
codificador funcional, tal como exige el documento del proyecto.
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch


PROJECT_ROOT = Path(__file__).resolve().parents[2]
INPUT = PROJECT_ROOT / "salidas_R_funcional" / "para_python"
OUTPUT = PROJECT_ROOT / "KAN" / "salidas" / "kan_oficial"
ARCHITECTURE = OUTPUT / "arquitectura"
RESULTS = OUTPUT / "resultados"
FIGURES = OUTPUT / "graficos"

SEEDS = (2026, 2027, 2028)
N_FUNCTIONAL_NEURONS = 8
HIDDEN_WIDTH = 4
GRID_INTERVALS = 5
SPLINE_DEGREE = 3
GRID_RANGE = (-3.0, 3.0)
LEARNING_RATE = 1e-2
LAMBDA_ENCODER_ROUGHNESS = 1e6
MAX_EPOCHS = 5000
PATIENCE = 1000

torch.set_default_dtype(torch.float64)
torch.set_num_threads(4)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
if DEVICE.type != "cuda":
    raise RuntimeError("La ejecucion oficial requiere una GPU CUDA compatible")


def basis_at_nodes(t: np.ndarray, knots: np.ndarray, order: int = 4) -> np.ndarray:
    """Evalua la base B-spline definida por su vector de nodos extendido."""
    t = np.asarray(t, dtype=float)
    knots = np.asarray(knots, dtype=float)
    basis = np.zeros((len(t), len(knots) - 1))
    for j in range(len(knots) - 1):
        basis[:, j] = (t >= knots[j]) & (t < knots[j + 1])
    basis[t == knots[-1], np.max(np.flatnonzero(knots[:-1] < knots[1:]))] = 1
    for degree in range(2, order + 1):
        updated = np.zeros((len(t), len(knots) - degree))
        for j in range(updated.shape[1]):
            left_denominator = knots[j + degree - 1] - knots[j]
            right_denominator = knots[j + degree] - knots[j + 1]
            if left_denominator > 0:
                updated[:, j] += (t - knots[j]) / left_denominator * basis[:, j]
            if right_denominator > 0:
                updated[:, j] += (knots[j + degree] - t) / right_denominator * basis[:, j + 1]
        basis = updated
    return basis


def gram_matrix(knots: np.ndarray, n_gauss_legendre: int = 4) -> np.ndarray:
    """Calcula la matriz de Gram por cuadratura de Gauss-Legendre."""
    unique_knots = np.unique(knots)
    nodes, weights = np.polynomial.legendre.leggauss(n_gauss_legendre)
    left = unique_knots[:-1, None]
    right = unique_knots[1:, None]
    quadrature_nodes = (((right - left) / 2) * nodes + (left + right) / 2).ravel()
    quadrature_weights = (((right - left) / 2) * weights).ravel()
    basis = basis_at_nodes(quadrature_nodes, knots)
    return basis.T @ (quadrature_weights[:, None] * basis)


def load_problem() -> dict:
    specification = json.loads((INPUT / "especificacion.json").read_text(encoding="utf-8"))
    pairs = pd.read_csv(INPUT / "pares.csv", dtype={"codigo": str})
    expected = {"entrenamiento": 2257, "validación": 672, "prueba": 370}
    if pairs["particion"].value_counts().to_dict() != expected:
        raise ValueError("La particion temporal no coincide con la version consolidada")
    if not pairs["usar"].all():
        raise ValueError("Existen pares no elegibles en las matrices funcionales")

    train_global = np.flatnonzero(pairs["particion"].to_numpy() == "entrenamiento")
    validation_global = np.flatnonzero(pairs["particion"].to_numpy() == "validación")
    test_global = np.flatnonzero(pairs["particion"].to_numpy() == "prueba")
    eligible_global = np.concatenate([train_global, validation_global])
    eligible_set = set(eligible_global.tolist())

    read_matrix = lambda name, **kwargs: pd.read_csv(INPUT / name, **kwargs).to_numpy(float)
    satellite_coefficients_all = read_matrix("coef_satelite_raiz.csv")
    station_coefficients = read_matrix(
        "coef_estacion_raiz.csv",
        skiprows=lambda row: row > 0 and row - 1 not in eligible_set,
    )
    station_coefficients = station_coefficients[
        np.searchsorted(np.sort(eligible_global), eligible_global)
    ]
    eligible_pairs = pairs.iloc[eligible_global].reset_index(drop=True)
    altitude = (
        eligible_pairs["altitud"].to_numpy(float) - specification["altitud_media"]
    ) / specification["altitud_de"]
    oni = eligible_pairs["ONI_anual"].to_numpy(float)
    fpca = pd.read_csv(INPUT / "fpca_estaciones.csv")

    if len(test_global) != 370:
        raise ValueError("La particion de prueba fue modificada")
    if len(station_coefficients) != len(eligible_global):
        raise ValueError("Las curvas terrestres no coinciden con entrenamiento y validacion")

    return {
        "specification": specification,
        "pairs": pairs,
        "eligible_pairs": eligible_pairs,
        "train": np.arange(len(train_global)),
        "validation": np.arange(len(train_global), len(eligible_global)),
        "train_global": train_global,
        "validation_global": validation_global,
        "satellite": satellite_coefficients_all[eligible_global],
        "station": station_coefficients,
        "altitude": altitude,
        "oni": oni,
        "station_knots": np.asarray(specification["nodos_estaciones"], dtype=float),
        "satellite_knots": np.asarray(specification["nodos_satelite"], dtype=float),
        "station_roughness": read_matrix("R_s.csv", header=None),
        "fpca_mean": fpca["mu"].to_numpy(float),
        "fpca_basis": fpca.drop(columns="mu").to_numpy(float),
    }


DATA = load_problem()
W_STATION = gram_matrix(DATA["station_knots"])
W_SATELLITE = gram_matrix(DATA["satellite_knots"])
N_STATION_BASIS = DATA["station"].shape[1]
N_SATELLITE_BASIS = DATA["satellite"].shape[1]
N_FPCA = DATA["fpca_basis"].shape[1]
WIDTHS = (N_FUNCTIONAL_NEURONS + 2, HIDDEN_WIDTH, N_FPCA)


def tensor(values: np.ndarray) -> torch.Tensor:
    return torch.as_tensor(np.array(values, dtype=float, copy=True), device=DEVICE)


X = tensor(DATA["satellite"])
Y = tensor(DATA["station"])
EXTRA = tensor(np.column_stack([DATA["altitude"], DATA["oni"]]))
W_E = tensor(W_STATION)
W_S = tensor(W_SATELLITE)
R_S = tensor(DATA["station_roughness"])


class KANLayer(torch.nn.Module):
    """Capa KAN B-spline exactamente con la forma indicada en el enunciado."""

    def __init__(
        self,
        n_input: int,
        n_output: int,
        grid_intervals: int = GRID_INTERVALS,
        degree: int = SPLINE_DEGREE,
        grid_range: tuple[float, float] = GRID_RANGE,
    ):
        super().__init__()
        self.degree = degree
        self.grid_intervals = grid_intervals
        self.grid_range = grid_range
        step = (grid_range[1] - grid_range[0]) / grid_intervals
        grid = (
            torch.arange(
                -self.degree,
                grid_intervals + self.degree + 1,
                dtype=torch.float64,
            )
            * step
            + grid_range[0]
        )
        self.register_buffer("grid", grid)
        self.coefficients = torch.nn.Parameter(
            0.1 * torch.randn(n_input, n_output, grid_intervals + self.degree)
        )
        self.base_weight = torch.nn.Parameter(torch.randn(n_input, n_output) / math.sqrt(n_input))
        self.spline_weight = torch.nn.Parameter(torch.ones(n_input, n_output))

    def bases(self, values: torch.Tensor) -> torch.Tensor:
        values = values.unsqueeze(-1)
        basis = ((values >= self.grid[:-1]) & (values < self.grid[1:])).to(values.dtype)
        for degree in range(1, self.degree + 1):
            basis = (
                (values - self.grid[: -(degree + 1)])
                / (self.grid[degree:-1] - self.grid[: -(degree + 1)])
                * basis[..., :-1]
                + (self.grid[degree + 1 :] - values)
                / (self.grid[degree + 1 :] - self.grid[1:-degree])
                * basis[..., 1:]
            )
        return basis

    def edge_activations(self, values: torch.Tensor) -> torch.Tensor:
        spline = torch.einsum("nig,iog->nio", self.bases(values), self.coefficients)
        base = self.base_weight * torch.nn.functional.silu(values).unsqueeze(-1)
        return base + self.spline_weight * spline

    def forward(self, values: torch.Tensor) -> torch.Tensor:
        return self.edge_activations(values).sum(dim=1)


class KANHead(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.layers = torch.nn.ModuleList(
            [KANLayer(WIDTHS[0], WIDTHS[1]), KANLayer(WIDTHS[1], WIDTHS[2])]
        )

    def forward(self, values: torch.Tensor) -> torch.Tensor:
        for layer in self.layers:
            values = layer(values)
        return values

    def number_of_parameters(self) -> int:
        return int(sum(parameter.numel() for parameter in self.parameters()))


class FunctionalKAN(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.encoder = torch.nn.Parameter(
            0.1 * torch.randn(N_FUNCTIONAL_NEURONS, N_SATELLITE_BASIS)
        )
        self.head = KANHead()
        self.register_buffer("fpca_mean", tensor(DATA["fpca_mean"]))
        self.register_buffer("fpca_basis", tensor(DATA["fpca_basis"]))

    def latent(self, coefficients: torch.Tensor) -> torch.Tensor:
        return coefficients @ W_S @ self.encoder.T

    def head_input(self, coefficients: torch.Tensor, extra: torch.Tensor) -> torch.Tensor:
        train_index = torch.as_tensor(DATA["train"], device=DEVICE)
        train_latent = self.latent(X[train_index])
        center = train_latent.mean(dim=0)
        scale = train_latent.std(dim=0)
        standardized = (self.latent(coefficients) - center) / scale
        return torch.cat([standardized, extra], dim=1)

    def forward(self, coefficients: torch.Tensor, extra: torch.Tensor) -> torch.Tensor:
        scores = self.head(self.head_input(coefficients, extra))
        return self.fpca_mean + scores @ self.fpca_basis.T

    def encoder_roughness(self) -> torch.Tensor:
        numerator = torch.einsum("ki,ij,kj->k", self.encoder, R_S, self.encoder)
        denominator = torch.einsum("ki,ij,kj->k", self.encoder, W_S, self.encoder)
        return (numerator / denominator).sum()


def functional_loss(model: FunctionalKAN, indices: torch.Tensor) -> torch.Tensor:
    difference = model(X[indices], EXTRA[indices]) - Y[indices]
    return torch.einsum("ni,ij,nj->n", difference, W_E, difference).mean()


def fit(seed: int) -> tuple[FunctionalKAN, dict]:
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    # La cabeza se inicializa primero y se conserva. Luego se reinicia la misma
    # semilla para crear el codificador. Este orden replica exactamente el
    # control experimental usado al comparar la base con las extensiones: la
    # presencia o ausencia de una penalizacion no cambia los pesos iniciales.
    initial_head = KANHead().to(DEVICE)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    model = FunctionalKAN().to(DEVICE)
    model.head.load_state_dict(initial_head.state_dict())
    optimizer = torch.optim.Adam(model.parameters(), lr=LEARNING_RATE)
    train_index = torch.as_tensor(DATA["train"], dtype=torch.long, device=DEVICE)
    validation_index = torch.as_tensor(DATA["validation"], dtype=torch.long, device=DEVICE)
    best_loss = np.inf
    best_epoch = 0
    best_state = None
    started = time.time()

    for epoch in range(MAX_EPOCHS):
        model.train()
        optimizer.zero_grad()
        loss_data = functional_loss(model, train_index)
        objective = loss_data + LAMBDA_ENCODER_ROUGHNESS * model.encoder_roughness()
        objective.backward()
        optimizer.step()

        model.eval()
        with torch.no_grad():
            validation = float(functional_loss(model, validation_index).item())
        if validation < best_loss - 1e-10:
            best_loss = validation
            best_epoch = epoch
            best_state = {name: value.detach().clone() for name, value in model.state_dict().items()}
        elif epoch - best_epoch > PATIENCE:
            break

    if best_state is None:
        raise RuntimeError("No se obtuvo un estado valido durante el entrenamiento")
    model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        train_loss = float(functional_loss(model, train_index).item())
    information = {
        "semilla": seed,
        "epoca_optima": best_epoch + 1,
        "epocas_corridas": epoch + 1,
        "L2_funcional_entrenamiento": train_loss,
        "L2_funcional_validacion": best_loss,
        "segundos": round(time.time() - started, 2),
        "parametros_cabeza": model.head.number_of_parameters(),
        "parametros_totales": int(sum(parameter.numel() for parameter in model.parameters())),
    }
    return model, information


def validation_metrics(model: FunctionalKAN, seed: int) -> tuple[pd.DataFrame, float]:
    validation_index = torch.as_tensor(DATA["validation"], device=DEVICE)
    with torch.no_grad():
        predicted_coefficients = model(X[validation_index], EXTRA[validation_index]).cpu().numpy()

    days = np.arange(1, 366, dtype=float)
    basis = basis_at_nodes(days, DATA["station_knots"])
    observed_coefficients = pd.read_csv(INPUT / "coef_estacion_mm.csv").to_numpy(float)
    observed = np.maximum(observed_coefficients[DATA["validation_global"]] @ basis.T, 0)
    predicted = np.maximum(predicted_coefficients @ basis.T, 0) ** 2
    l2_relative = np.sqrt(
        np.sum((predicted - observed) ** 2, axis=1)
        / np.maximum(np.sum(observed**2, axis=1), 1e-12)
    )
    pairs = DATA["pairs"].iloc[DATA["validation_global"]].reset_index(drop=True)
    frame = pd.DataFrame(
        {
            "semilla": seed,
            "fila_global": DATA["validation_global"] + 1,
            "codigo": pairs["codigo"],
            "anio": pairs["anio"].astype(int),
            "L2_relativo": l2_relative,
        }
    )
    by_station = frame.groupby("codigo", as_index=False)["L2_relativo"].mean()
    return frame, float(by_station["L2_relativo"].mean())


def architecture_figure() -> None:
    fig, ax = plt.subplots(figsize=(13, 7))
    ax.axis("off")
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 11)
    input_y = np.linspace(1, 10, 10)
    hidden_y = np.linspace(2.2, 8.8, 4)
    output_y = np.linspace(1.5, 9.5, 6)
    for y1 in input_y:
        for y2 in hidden_y:
            ax.plot([2.2, 5.8], [y1, y2], color="#cbd5e1", lw=0.7)
    for y1 in hidden_y:
        for y2 in output_y:
            ax.plot([6.2, 9.8], [y1, y2], color="#cbd5e1", lw=0.7)
    names = [f"z{i}" for i in range(1, 9)] + ["altitud", "ONI"]
    for y, name in zip(input_y, names):
        ax.scatter(2, y, s=180, color="#2563eb", zorder=3)
        ax.text(1.65, y, name, ha="right", va="center", fontsize=9)
    for index, y in enumerate(hidden_y, 1):
        ax.scatter(6, y, s=240, color="#f59e0b", zorder=3)
        ax.text(6, y, f"h{index}", ha="center", va="center", fontsize=9, weight="bold")
    for index, y in enumerate(output_y, 1):
        ax.scatter(10, y, s=190, color="#16a34a", zorder=3)
        ax.text(10.35, y, f"score FPCA {index}", ha="left", va="center", fontsize=9)
    ax.text(2, 10.55, "10 entradas", ha="center", fontsize=13, weight="bold")
    ax.text(6, 9.45, "4 nodos ocultos", ha="center", fontsize=13, weight="bold")
    ax.text(10, 10.55, "6 salidas FPCA", ha="center", fontsize=13, weight="bold")
    ax.text(
        6,
        0.35,
        r"Aristas: $\phi_{ij}(x)=w_b\,SiLU(x)+w_s\sum_g c_gB_g(x)$",
        ha="center",
        fontsize=13,
        bbox={"boxstyle": "round,pad=.5", "fc": "#eff6ff", "ec": "#2563eb"},
    )
    ax.set_title("KAN base oficial [10, 4, 6] — sin regularización adicional en la cabeza", fontsize=16, weight="bold")
    fig.tight_layout()
    fig.savefig(FIGURES / "01_arquitectura_KAN_base_oficial.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def validation_figure(summary: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(summary["semilla"].astype(str), summary["L2_relativo_media"], color="#2563eb")
    ax.axhline(summary["L2_relativo_media"].mean(), color="#dc2626", ls="--", label="media")
    ax.set_xlabel("Semilla")
    ax.set_ylabel("L² relativo medio por estación")
    ax.set_title("KAN base oficial: estabilidad en validación 1996–2000")
    ax.grid(axis="y", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES / "02_resultado_tres_semillas.png", dpi=220)
    plt.close(fig)


def write_readme(summary: pd.DataFrame, training: pd.DataFrame) -> None:
    mean = summary["L2_relativo_media"].mean()
    deviation = summary["L2_relativo_media"].std(ddof=1)
    lines = [
        "# KAN base oficial",
        "",
        "Esta carpeta contiene la KAN base alineada con el enunciado del proyecto.",
        "La cabeza usa exclusivamente B-splines cúbicas y SiLU en sus aristas; no",
        "incorpora regularizadores adicionales sobre sus coeficientes.",
        "",
        "## Arquitectura",
        "",
        "```text",
        "Curva satelital -> codificador funcional (8 neuronas)",
        "-> [z1-z8, altitud, ONI] -> KAN [10,4,6] -> FPCA fijo -> curva terrestre",
        "```",
        "",
        "- Rejilla: G=5 intervalos sobre [-3,3].",
        "- B-splines cúbicas.",
        "- Learning rate: 0.01.",
        "- Rugosidad del codificador: lambda_w=1e6.",
        "- Semillas: 2026, 2027 y 2028.",
        "- Selección y parada temprana: validación 1996-2000.",
        "- Prueba 2001-2003: no utilizada.",
        "",
        "## Resultado",
        "",
        f"L² relativo medio por estación: **{mean:.6f} ± {deviation:.6f}**.",
        "",
        "| Semilla | L² relativo medio | Época óptima |",
        "|---:|---:|---:|",
    ]
    for row in summary.itertuples(index=False):
        epoch = int(training.loc[training["semilla"] == row.semilla, "epoca_optima"].iloc[0])
        lines.append(f"| {int(row.semilla)} | {row.L2_relativo_media:.6f} | {epoch} |")
    synthetic_summary_path = OUTPUT / "prueba_sintetica" / "resultados" / "resumen.json"
    if synthetic_summary_path.is_file():
        synthetic = json.loads(synthetic_summary_path.read_text(encoding="utf-8"))
        lines.extend(
            [
                "",
                "## Verificación sintética obligatoria",
                "",
                f"Función: `{synthetic['funcion']}` sobre entradas en [-2,2].",
                "",
                f"RMSE medio: **{synthetic['RMSE_medio']:.8f}**; RMSE máximo entre tres semillas: "
                f"**{synthetic['RMSE_maximo']:.8f}**.",
                "",
                f"Umbral del documento: RMSE <= {synthetic['umbral_RMSE']:.2f}. "
                f"Resultado: **{'CUMPLE' if synthetic['cumple_todas_las_semillas'] else 'NO CUMPLE'}**.",
            ]
        )
    lines.extend(
        [
            "",
            "## Alcance metodológico",
            "",
            "Los valores concretos de ancho, G, tasa de aprendizaje y lambda_w son",
            "hiperparámetros elegidos con validación. No son técnicas externas: son",
            "decisiones necesarias dentro de la arquitectura permitida por el documento.",
        ]
    )
    (OUTPUT / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    for folder in (ARCHITECTURE, RESULTS, FIGURES):
        folder.mkdir(parents=True, exist_ok=True)

    training_rows = []
    validation_rows = []
    summary_rows = []
    models = {}
    for seed in SEEDS:
        print(f"Entrenando KAN base oficial, semilla {seed} en {torch.cuda.get_device_name(0)}", flush=True)
        model, information = fit(seed)
        details, mean_l2 = validation_metrics(model, seed)
        information["L2_relativo_validacion"] = mean_l2
        training_rows.append(information)
        validation_rows.append(details)
        summary_rows.append({"semilla": seed, "L2_relativo_media": mean_l2})
        models[seed] = model
        torch.save(model.state_dict(), ARCHITECTURE / f"modelo_semilla_{seed}.pt")
        print(f"Semilla {seed}: L2 relativo={mean_l2:.6f}; epoca={information['epoca_optima']}", flush=True)

    training = pd.DataFrame(training_rows)
    validation = pd.concat(validation_rows, ignore_index=True)
    summary = pd.DataFrame(summary_rows)
    training.to_csv(RESULTS / "entrenamiento_tres_semillas.csv", index=False)
    validation.to_csv(RESULTS / "metricas_validacion_por_curva.csv", index=False)
    summary.to_csv(RESULTS / "L2_validacion_por_semilla.csv", index=False)

    head_parameters = models[SEEDS[0]].head.number_of_parameters()
    total_parameters = int(sum(parameter.numel() for parameter in models[SEEDS[0]].parameters()))
    configuration = {
        "nombre": "KAN base oficial",
        "arquitectura": list(WIDTHS),
        "neuronas_funcionales": N_FUNCTIONAL_NEURONS,
        "G_intervalos": GRID_INTERVALS,
        "grado_B_spline": SPLINE_DEGREE,
        "rango_rejilla": list(GRID_RANGE),
        "funcion_arista": "w_b*SiLU(x) + w_s*sum_g(c_g*B_g(x))",
        "regularizacion_cabeza": "ninguna",
        "lambda_w_codificador": LAMBDA_ENCODER_ROUGHNESS,
        "learning_rate": LEARNING_RATE,
        "max_epochs": MAX_EPOCHS,
        "patience": PATIENCE,
        "semillas": list(SEEDS),
        "parametros_cabeza": head_parameters,
        "parametros_totales": total_parameters,
        "particion_seleccion": "validacion 1996-2000",
        "prueba_2001_2003_usada": False,
        "dispositivo": str(DEVICE),
        "gpu": torch.cuda.get_device_name(0),
        "L2_relativo_media": float(summary["L2_relativo_media"].mean()),
        "L2_relativo_DE": float(summary["L2_relativo_media"].std(ddof=1)),
    }
    (ARCHITECTURE / "configuracion.json").write_text(
        json.dumps(configuration, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    pd.DataFrame(
        [
            {"capa": "entrada", "unidades": 10, "descripcion": "z1-z8, altitud, ONI"},
            {"capa": "KAN oculta", "unidades": 4, "descripcion": "B-splines cubicas + SiLU en aristas"},
            {"capa": "salida", "unidades": 6, "descripcion": "puntajes FPCA terrestres"},
        ]
    ).to_csv(ARCHITECTURE / "arquitectura.csv", index=False)

    architecture_figure()
    validation_figure(summary)
    write_readme(summary, training)
    controls = pd.DataFrame(
        [
            {"control": "GPU utilizada", "resultado": True, "detalle": torch.cuda.get_device_name(0)},
            {"control": "tres semillas", "resultado": len(summary) == 3, "detalle": str(SEEDS)},
            {"control": "solo validacion", "resultado": True, "detalle": "1996-2000"},
            {"control": "prueba cerrada", "resultado": True, "detalle": "2001-2003 no evaluada"},
            {"control": "regularizacion de cabeza ausente", "resultado": True, "detalle": "ninguna"},
            {"control": "parametros cabeza", "resultado": head_parameters == 640, "detalle": str(head_parameters)},
        ]
    )
    controls.to_csv(RESULTS / "controles.csv", index=False)
    print(json.dumps(configuration, indent=2, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()

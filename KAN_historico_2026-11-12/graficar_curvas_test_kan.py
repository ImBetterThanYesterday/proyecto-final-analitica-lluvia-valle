"""Superpone todas las curvas observadas y reconstruidas del conjunto de test."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from matplotlib.axes import Axes
from matplotlib.lines import Line2D

from kan_funcional_jacobo import (
    DEVICE,
    N_COMPONENTES_FPCA,
    RedFuncional,
    ajustar_fpca,
    cargar_datos,
    preparar_tensores,
    predecir_coef,
)


KAN_DIR = Path(__file__).resolve().parent
RESULTADOS = KAN_DIR / "resultados_reconstruccion"
GRAFICA = (
    RESULTADOS
    / "graficas_sustentacion"
    / "01_curvas_test_observadas_vs_kan.png"
)
GRAFICA_SIN_PROMEDIOS = (
    RESULTADOS
    / "graficas_sustentacion"
    / "04_curvas_test_sin_promedios.png"
)
GRAFICA_TEST_COMPLETO = (
    RESULTADOS
    / "graficas_sustentacion"
    / "05_curvas_test_todos_los_pares.png"
)
SEMILLA = 2026
PARES_POR_ANIO = 20


def _agregar_leyenda_curvas(ax: Axes) -> None:
    ax.legend(
        handles=[
            Line2D(
                [0],
                [0],
                color="#0057d9",
                linewidth=2.5,
                label="Estación observada",
            ),
            Line2D(
                [0],
                [0],
                color="#e31a1c",
                linewidth=2.5,
                label="Reconstrucción KAN",
            ),
        ]
    )


def generar_grafica_consolidada() -> Path:
    resumen_path = RESULTADOS / "resumen_reconstruccion.json"
    modelo_path = RESULTADOS / f"modelo_kan_seed{SEMILLA}.pt"
    metricas_path = RESULTADOS / "metricas_test_por_curva.csv"
    for path in (resumen_path, modelo_path, metricas_path):
        if not path.is_file():
            raise FileNotFoundError(f"No se encontró el insumo requerido: {path}")

    resumen = json.loads(resumen_path.read_text(encoding="utf-8"))
    if not resumen.get("prueba_utilizada", False):
        raise RuntimeError(
            "El modelo aún no se ha evaluado en test; primero debe ejecutarse "
            "la evaluación final con --evaluar-test."
        )

    datos = cargar_datos()
    tensores = preparar_tensores(datos)
    fpca = ajustar_fpca(
        tensores["y"][tensores["train"]],
        datos["W_est"],
        N_COMPONENTES_FPCA,
    )
    modelo = RedFuncional(
        "kan",
        datos["W_sat"],
        datos["phi_sat"],
        fpca["media"],
        fpca["direcciones"],
    ).to(DEVICE)
    estado = torch.load(modelo_path, map_location=DEVICE, weights_only=True)
    modelo.load_state_dict(estado)

    mascara_test = tensores["test"]
    coeficientes = predecir_coef(modelo, tensores, mascara_test)
    curvas_observadas = np.maximum(
        datos["y_mm_coef"][mascara_test] @ datos["phi_est"].T,
        0,
    )
    curvas_reconstruidas = np.maximum(
        coeficientes @ datos["phi_est"].T,
        0,
    ) ** 2
    ids_test = datos["ids"].loc[
        mascara_test,
        ["station_id", "station_name", "year"],
    ].reset_index(drop=True)
    metricas_test = pd.read_csv(metricas_path)
    if len(ids_test) != len(metricas_test):
        raise ValueError(
            "El número de filas de test no coincide con las métricas guardadas."
        )
    if not np.array_equal(
        ids_test["station_id"].astype(str).to_numpy(),
        metricas_test["station_id"].astype(str).to_numpy(),
    ) or not np.array_equal(ids_test["year"].to_numpy(), metricas_test["year"].to_numpy()):
        raise ValueError(
            "Los identificadores de las métricas no coinciden con el orden del test."
        )

    dias = np.arange(1, curvas_observadas.shape[1] + 1)
    n_curvas = len(curvas_observadas)
    observada_p10, observada_p90 = np.quantile(
        curvas_observadas, [0.1, 0.9], axis=0
    )
    predicha_p10, predicha_p90 = np.quantile(
        curvas_reconstruidas, [0.1, 0.9], axis=0
    )

    fig, ax = plt.subplots(figsize=(13, 7))
    ax.fill_between(
        dias,
        observada_p10,
        observada_p90,
        color="#555555",
        alpha=0.12,
        label="Observadas: percentiles 10–90",
    )
    ax.fill_between(
        dias,
        predicha_p10,
        predicha_p90,
        color="#2a9d8f",
        alpha=0.16,
        label="KAN: percentiles 10–90",
    )
    for indice in range(n_curvas):
        ax.plot(
            dias,
            curvas_observadas[indice],
            color="#444444",
            linewidth=0.55,
            alpha=0.055,
        )
        ax.plot(
            dias,
            curvas_reconstruidas[indice],
            color="#2a9d8f",
            linewidth=0.55,
            alpha=0.055,
        )
    ax.plot(
        dias,
        np.median(curvas_observadas, axis=0),
        color="black",
        linewidth=2.2,
        label="Mediana observada",
    )
    ax.plot(
        dias,
        np.median(curvas_reconstruidas, axis=0),
        color="#007f6e",
        linewidth=2.2,
        label="Mediana predicha por KAN",
    )
    ax.set(
        title=(
            f"Test: curvas observadas y reconstruidas "
            f"({n_curvas} pares, {ids_test['station_id'].nunique()} estaciones)"
        ),
        xlabel="Día del año",
        ylabel="Precipitación (mm/día)",
    )
    ax.grid(alpha=0.2)
    ax.legend(ncol=2)
    fig.tight_layout()
    GRAFICA.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(GRAFICA, dpi=180)
    plt.close(fig)

    rng = np.random.default_rng(SEMILLA)
    indices_muestra = []
    estaciones_seleccionadas = set()
    anios_test = ids_test["year"].to_numpy()
    estaciones_test = ids_test["station_id"].astype(str).to_numpy()
    for anio in np.sort(ids_test["year"].unique()):
        candidatos = np.flatnonzero(
            (anios_test == anio)
            & ~np.isin(estaciones_test, list(estaciones_seleccionadas))
        )
        if len(candidatos) < PARES_POR_ANIO:
            raise ValueError(
                f"No hay {PARES_POR_ANIO} estaciones distintas disponibles "
                f"para muestrear el año {anio}."
            )
        seleccion = rng.choice(candidatos, size=PARES_POR_ANIO, replace=False)
        indices_muestra.extend(seleccion.tolist())
        estaciones_seleccionadas.update(estaciones_test[seleccion])

    fig, ax = plt.subplots(figsize=(13, 7))
    for posicion, indice in enumerate(indices_muestra):
        ax.plot(
            dias,
            curvas_observadas[indice],
            color="#0057d9",
            linewidth=0.9,
            alpha=0.6,
        )
        ax.plot(
            dias,
            curvas_reconstruidas[indice],
            color="#e31a1c",
            linewidth=0.9,
            alpha=0.6,
        )
    ax.set(
        title=(
            f"Test: observadas y reconstrucciones KAN "
            f"({len(indices_muestra)} pares, {PARES_POR_ANIO} por año, "
            f"{len(estaciones_seleccionadas)} estaciones distintas)"
        ),
        xlabel="Día del año",
        ylabel="Precipitación (mm/día)",
    )
    ax.grid(alpha=0.2)
    _agregar_leyenda_curvas(ax)
    fig.tight_layout()
    GRAFICA_SIN_PROMEDIOS.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(GRAFICA_SIN_PROMEDIOS, dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(13, 7))
    for indice in range(n_curvas):
        ax.plot(
            dias,
            curvas_observadas[indice],
            color="#0057d9",
            linewidth=0.65,
            alpha=0.24,
        )
        ax.plot(
            dias,
            curvas_reconstruidas[indice],
            color="#e31a1c",
            linewidth=0.65,
            alpha=0.24,
        )
    ax.set(
        title=(
            f"Test: todas las curvas observadas y reconstruidas "
            f"({n_curvas} pares, {ids_test['station_id'].nunique()} estaciones)"
        ),
        xlabel="Día del año",
        ylabel="Precipitación (mm/día)",
    )
    ax.grid(alpha=0.2)
    _agregar_leyenda_curvas(ax)
    fig.tight_layout()
    fig.savefig(GRAFICA_TEST_COMPLETO, dpi=180)
    plt.close(fig)
    return GRAFICA


if __name__ == "__main__":
    print(f"Se generó una gráfica con todas las curvas de test en: {generar_grafica_consolidada()}")
    print(f"Se generó la gráfica sin promedios en: {GRAFICA_SIN_PROMEDIOS}")
    print(f"Se generó la gráfica con todos los pares de test en: {GRAFICA_TEST_COMPLETO}")

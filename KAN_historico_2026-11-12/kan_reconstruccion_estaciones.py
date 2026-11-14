"""Entrena una KAN para reconstruir curvas de estaciones desde datos satelitales."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch

from kan_funcional_jacobo import (
    KAN_BASES,
    N_COMPONENTES_FPCA,
    N_NEURONAS_FUNCIONALES,
    ajustar_fpca,
    cargar_datos,
    entrenar,
    metricas_curva,
    preparar_tensores,
    predecir_coef,
)


KAN_DIR = Path(__file__).resolve().parent
OUT = KAN_DIR / "resultados_reconstruccion"
SEED = 2026


def graficar_curva_aprendizaje(historia: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(historia["epoca"], historia["train_l2"], label="Entrenamiento")
    ax.plot(historia["epoca"], historia["val_l2"], label="Validación")
    ax.set(
        title="KAN: pérdida funcional por época",
        xlabel="Época",
        ylabel="Pérdida L² funcional (escala raíz)",
    )
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "curva_aprendizaje.png", dpi=160)
    plt.close(fig)


def graficar_curvas_ejemplo(
    ids_validacion: pd.DataFrame,
    metricas: pd.DataFrame,
    curvas_satelitales: np.ndarray,
    curvas_observadas: np.ndarray,
    curvas_reconstruidas: np.ndarray,
) -> pd.DataFrame:
    orden = np.argsort(metricas["l2_relativo"].to_numpy())
    indices = [orden[0], orden[len(orden) // 2], orden[-1]]
    etiquetas = ["Menor error", "Error mediano", "Mayor error"]

    fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    filas = []
    dias = np.arange(1, curvas_observadas.shape[1] + 1)
    for ax, indice, etiqueta in zip(axes, indices, etiquetas):
        identificador = ids_validacion.iloc[indice]
        error = float(metricas.iloc[indice]["l2_relativo"])
        ax.plot(dias, curvas_observadas[indice], label="Estación observada", color="black")
        ax.plot(
            dias,
            curvas_satelitales[indice],
            label="Satélite de entrada",
            color="#e76f51",
            alpha=0.8,
        )
        ax.plot(
            dias,
            curvas_reconstruidas[indice],
            label="Reconstrucción KAN",
            color="#2a9d8f",
        )
        ax.set_title(
            f"{etiqueta}: {identificador['station_id']} — {identificador['year']} "
            f"(L2 relativo={error:.3f})"
        )
        ax.set_ylabel("Precipitación (mm/día)")
        ax.grid(alpha=0.25)
        filas.append({
            "caso": etiqueta,
            "station_id": identificador["station_id"],
            "year": int(identificador["year"]),
            "l2_relativo": error,
        })

    axes[0].legend(ncol=3)
    axes[-1].set_xlabel("Día del año")
    fig.tight_layout()
    fig.savefig(OUT / "curvas_validacion_ejemplos.png", dpi=160)
    plt.close(fig)
    return pd.DataFrame(filas)


def graficar_error_por_estacion(metricas_por_curva: pd.DataFrame) -> pd.DataFrame:
    por_estacion = (
        metricas_por_curva.groupby("station_id", as_index=False)
        .agg(l2_relativo_medio=("l2_relativo", "mean"), pares_validacion=("year", "size"))
        .sort_values("l2_relativo_medio", ascending=False)
        .reset_index(drop=True)
    )
    por_estacion.to_csv(OUT / "l2_relativo_por_estacion.csv", index=False)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(por_estacion["l2_relativo_medio"], bins="auto", color="#457b9d", alpha=0.85)
    ax.axvline(
        por_estacion["l2_relativo_medio"].mean(),
        color="#e63946",
        linestyle="--",
        label=f"Media={por_estacion['l2_relativo_medio'].mean():.3f}",
    )
    ax.axvline(
        por_estacion["l2_relativo_medio"].median(),
        color="#2a9d8f",
        linestyle=":",
        label=f"Mediana={por_estacion['l2_relativo_medio'].median():.3f}",
    )
    ax.set(
        title="Distribución del error L² relativo medio por estación",
        xlabel="L² relativo medio por estación (menor es mejor)",
        ylabel="Número de estaciones",
    )
    ax.grid(axis="y", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "distribucion_l2_por_estacion.png", dpi=160)
    plt.close(fig)
    return por_estacion


def resumir_metricas(metricas: pd.DataFrame) -> dict:
    return {
        nombre: {
            "media": float(metricas[nombre].mean()),
            "desviacion_estandar": float(metricas[nombre].std()),
        }
        for nombre in metricas.columns
    }


def ejecutar_reconstruccion(evaluar_test: bool = False) -> dict:
    OUT.mkdir(exist_ok=True)
    datos = cargar_datos()
    tensores = preparar_tensores(datos)
    fpca = ajustar_fpca(
        tensores["y"][tensores["train"]],
        datos["W_est"],
        N_COMPONENTES_FPCA,
    )

    modelo, historia, mejor_l2_validacion, mejor_epoca = entrenar(
        "kan", SEED, datos, tensores, fpca
    )
    mascara_validacion = tensores["val"]
    coeficientes_predichos = predecir_coef(modelo, tensores, mascara_validacion)

    phi_est = datos["phi_est"]
    curvas_observadas = np.maximum(
        datos["y_mm_coef"][mascara_validacion] @ phi_est.T,
        0,
    )
    curvas_satelitales = np.maximum(
        datos["x_mm_coef"][mascara_validacion] @ datos["phi_sat"].T,
        0,
    )
    curvas_raiz_predichas = coeficientes_predichos @ phi_est.T
    curvas_reconstruidas = np.maximum(curvas_raiz_predichas, 0) ** 2
    metricas = metricas_curva(curvas_observadas, curvas_reconstruidas)
    ids_validacion = datos["ids"].loc[
        mascara_validacion,
        ["station_id", "station_name", "year", "fase_enso"],
    ].reset_index(drop=True)

    historia.to_csv(OUT / "historia_entrenamiento.csv", index=False)
    metricas_por_curva = pd.concat([ids_validacion, metricas], axis=1)
    metricas_por_curva.to_csv(OUT / "metricas_validacion_por_curva.csv", index=False)
    graficar_curva_aprendizaje(historia)
    casos_graficados = graficar_curvas_ejemplo(
        ids_validacion,
        metricas,
        curvas_satelitales,
        curvas_observadas,
        curvas_reconstruidas,
    )
    casos_graficados.to_csv(OUT / "casos_validacion_graficados.csv", index=False)
    metricas_estacion = graficar_error_por_estacion(metricas_por_curva)

    n_dias = curvas_reconstruidas.shape[1]
    columnas_dia = lambda prefijo: [
        f"{prefijo}_dia_{dia + 1}" for dia in range(n_dias)
    ]
    curvas_predicciones = pd.concat(
        [
            ids_validacion,
            pd.DataFrame(curvas_satelitales, columns=columnas_dia("satelite")),
            pd.DataFrame(curvas_observadas, columns=columnas_dia("observado")),
            pd.DataFrame(curvas_reconstruidas, columns=columnas_dia("reconstruido")),
        ],
        axis=1,
    )
    curvas_predicciones.to_csv(
        OUT / "curvas_validacion.csv.gz",
        index=False,
        compression="gzip",
    )
    torch.save(modelo.state_dict(), OUT / f"modelo_kan_seed{SEED}.pt")

    mascara_test = tensores["test"]
    ids_test = datos["ids"].loc[
        mascara_test,
        ["station_id", "station_name", "year", "fase_enso"],
    ].reset_index(drop=True)
    resumen_test = {
        "pares": int(mascara_test.sum()),
        "estaciones": int(ids_test["station_id"].nunique()),
        "evaluada": evaluar_test,
        "usada_para_seleccion_del_modelo": False,
    }
    if evaluar_test:
        coeficientes_test = predecir_coef(modelo, tensores, mascara_test)
        curvas_observadas_test = np.maximum(
            datos["y_mm_coef"][mascara_test] @ phi_est.T,
            0,
        )
        curvas_reconstruidas_test = np.maximum(
            coeficientes_test @ phi_est.T,
            0,
        ) ** 2
        metricas_test = metricas_curva(
            curvas_observadas_test,
            curvas_reconstruidas_test,
        )
        pd.concat([ids_test, metricas_test], axis=1).to_csv(
            OUT / "metricas_test_por_curva.csv",
            index=False,
        )
        estaciones_train = set(
            datos["ids"].loc[tensores["train"], "station_id"].astype(str)
        )
        estaciones_test = set(ids_test["station_id"].astype(str))
        resumen_test.update({
            "estaciones_compartidas_con_entrenamiento": len(
                estaciones_test & estaciones_train
            ),
            "estaciones_nuevas_vs_entrenamiento": len(
                estaciones_test - estaciones_train
            ),
            "metricas": resumir_metricas(metricas_test),
        })

    resumen = {
        "objetivo": "Reconstruir curvas de precipitación de estación con una KAN.",
        "entrada_modelo": [
            "curva satelital en escala raíz",
            "elevación",
            "ONI anual",
        ],
        "semilla": SEED,
        "filas": {
            "entrenamiento": int(tensores["train"].sum()),
            "validacion": int(tensores["val"].sum()),
            "prueba": int(tensores["test"].sum()),
        },
        "estaciones_validacion": int(ids_validacion["station_id"].nunique()),
        "prueba": resumen_test,
        "arquitectura": {
            "neuronas_funcionales": N_NEURONAS_FUNCIONALES,
            "componentes_fpca": N_COMPONENTES_FPCA,
            "varianza_fpca_acumulada": float(
                fpca["proporcion"][:N_COMPONENTES_FPCA].sum()
            ),
            "bases_spline_por_arista": KAN_BASES,
            "parametros_entrenables": int(
                sum(param.numel() for param in modelo.parameters() if param.requires_grad)
            ),
        },
        "entrenamiento": {
            "mejor_epoca": int(mejor_epoca),
            "perdida_funcional_validacion_raiz": float(mejor_l2_validacion),
        },
        "metricas_validacion": resumir_metricas(metricas),
        "metricas_por_estacion": {
            "n_estaciones": int(len(metricas_estacion)),
            "l2_relativo_medio": float(
                metricas_estacion["l2_relativo_medio"].mean()
            ),
            "l2_relativo_mediana": float(
                metricas_estacion["l2_relativo_medio"].median()
            ),
        },
        "graficas": [
            "curva_aprendizaje.png",
            "curvas_validacion_ejemplos.png",
            "distribucion_l2_por_estacion.png",
        ],
        "prueba_utilizada": evaluar_test,
        "emparejamiento_estacion_pixel_pendiente_confirmacion_docente": True,
    }
    (OUT / "resumen_reconstruccion.json").write_text(
        json.dumps(resumen, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    return resumen


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Entrena una KAN y evalúa la reconstrucción.")
    parser.add_argument(
        "--evaluar-test",
        action="store_true",
        help="Evalúa una sola vez en test, después de seleccionar el modelo con validación.",
    )
    args = parser.parse_args()
    print(json.dumps(
        ejecutar_reconstruccion(evaluar_test=args.evaluar_test),
        ensure_ascii=True,
        indent=2,
    ))

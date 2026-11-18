"""Compara hiperparametros KAN de uno en uno usando solo validacion."""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from kan_variantes_core import (
    DATA,
    DEVICE,
    K_E,
    LAMBDA_W,
    LR,
    M,
    MAX_EPOCHS,
    N_ENCODER,
    PATIENCE,
    ROOT,
    SEEDS,
    fit,
    predict_validation,
)


BASELINE = {
    "widths": (N_ENCODER + 2, 4, M),
    "grid_intervals": 5,
    "learning_rate": LR,
    "max_epochs": MAX_EPOCHS,
    "patience": PATIENCE,
    "lambda_w_encoder": LAMBDA_W,
    "lambda_second_difference": 0.0,
}


def build_cases():
    cases = []

    def add(parameter, value, **overrides):
        config = {**BASELINE, **overrides}
        cases.append({
            "parametro": parameter,
            "valor": str(value),
            **config,
        })

    for width in (3, 4, 5):
        add("ancho_oculto", width, widths=(N_ENCODER + 2, width, M))
    add("capas_kan", 2, widths=(N_ENCODER + 2, 4, M))
    add("capas_kan", 3, widths=(N_ENCODER + 2, 4, 3, M))
    for intervals in (3, 5, 7):
        add("intervalos_G", intervals, grid_intervals=intervals)
    for rate in (1e-3, 3e-3, 1e-2):
        add("tasa_aprendizaje", rate, learning_rate=rate)
    for wait in (200, 500, 1000):
        add("paciencia", wait, patience=wait)
    for cap in (250, 500, 1000, 5000):
        add("tope_epocas", cap, max_epochs=cap)
    for penalty in (1e6, 1e7, 1e8):
        add("lambda_w_encoder", penalty, lambda_w_encoder=penalty)
    for penalty in (0.0, 0.1, 1.0, 10.0, 100.0):
        add("lambda_D2", penalty, lambda_second_difference=penalty)
    return cases


def _cache_key(config, seed):
    return (
        tuple(config["widths"]),
        config["grid_intervals"],
        config["learning_rate"],
        config["max_epochs"],
        config["patience"],
        config["lambda_w_encoder"],
        config["lambda_second_difference"],
        seed,
    )


def main():
    cases = build_cases()
    cache = {}
    training_rows = []
    prediction_frames = []
    valid_pairs = DATA["valid_pairs"]
    coefficient_names = [f"c{i + 1}" for i in range(K_E)]

    for case_number, case in enumerate(cases, start=1):
        parameter, value = case["parametro"], case["valor"]
        scenario = f"{parameter}={value}"
        config = {key: case[key] for key in BASELINE}
        for seed in SEEDS:
            key = _cache_key(config, seed)
            if key not in cache:
                model, info = fit(
                    config["widths"],
                    config["grid_intervals"],
                    seed,
                    config["lambda_second_difference"],
                    lambda_w_encoder=config["lambda_w_encoder"],
                    learning_rate=config["learning_rate"],
                    max_epochs=config["max_epochs"],
                    patience=config["patience"],
                )
                prediction = predict_validation(model)
                cache[key] = (info, prediction)
                del model
                if DEVICE.type == "cuda":
                    torch.cuda.empty_cache()
            info, prediction = cache[key]
            training_rows.append({
                "parametro": parameter,
                "valor": value,
                "escenario": scenario,
                "arquitectura": str(config["widths"]),
                "G": config["grid_intervals"],
                **info,
            })
            frame = pd.DataFrame(prediction, columns=coefficient_names)
            frame.insert(0, "semilla", seed)
            frame.insert(0, "escenario", scenario)
            frame.insert(0, "valor", value)
            frame.insert(0, "parametro", parameter)
            frame.insert(0, "particion", "validación")
            frame.insert(0, "anio", valid_pairs["anio"].to_numpy())
            frame.insert(0, "codigo", valid_pairs["codigo"].to_numpy())
            frame.insert(0, "fila", DATA["valid_global_rows"] + 1)
            prediction_frames.append(frame)
        print(
            f"[{case_number}/{len(cases)}] {scenario}: "
            "tres semillas completas",
            flush=True,
        )

    output = Path(os.environ.get(
        "LLUVIA_SALIDA_KAN_PARAMETROS",
        str(ROOT / "KAN" / "salidas" / "analisis_parametros_kan"),
    ))
    output.mkdir(parents=True, exist_ok=True)
    training = pd.DataFrame(training_rows)
    predictions = pd.concat(prediction_frames, ignore_index=True)
    training.to_csv(output / "entrenamientos_validacion.csv", index=False)
    predictions.to_csv(
        output / "predicciones_validacion_raiz.csv",
        index=False,
        float_format="%.12g",
    )

    summary = (
        training.groupby(["parametro", "valor"], sort=False)
        .agg(
            L2_funcional_media=("L2_validacion", "mean"),
            L2_funcional_DE=("L2_validacion", "std"),
            epoca_optima_media=("epoca_optima", "mean"),
            epocas_corridas_media=("epocas_corridas", "mean"),
        )
        .reset_index()
    )
    summary.to_csv(output / "resumen_parametros_python.csv", index=False)

    grouped_predictions = predictions.groupby(["parametro", "valor", "semilla"])
    checks = {
        "Se incluyen solo predicciones de validación":
            predictions["particion"].eq("validación").all(),
        "Cada escenario y semilla contiene 672 curvas":
            grouped_predictions.size().eq(len(valid_pairs)).all(),
        "Se usaron tres semillas para cada escenario":
            training.groupby(["parametro", "valor"])["semilla"].nunique().eq(len(SEEDS)).all(),
        "Predicciones y pérdidas son finitas":
            np.isfinite(predictions[coefficient_names].to_numpy()).all()
            and np.isfinite(training[["L2_entrenamiento", "L2_validacion"]].to_numpy()).all(),
        "Las épocas óptimas fueron seleccionadas":
            training["epoca_optima"].gt(0).all(),
    }
    pd.DataFrame(
        [{"prueba": name, "pasa": bool(ok)} for name, ok in checks.items()]
    ).to_csv(output / "P6_pruebas.csv", index=False)
    if not all(checks.values()):
        raise RuntimeError(f"Fallaron controles del análisis de parámetros: {checks}")

    config_output = {
        "protocolo": "un parametro a la vez; tres semillas; seleccion y evaluacion en validacion",
        "semillas": SEEDS,
        "curvas_entrenamiento": len(DATA["train"]),
        "curvas_validacion": len(DATA["valid"]),
        "prueba_del_proyecto_usada": False,
        "cache_de_corridas_identicas": True,
        "escenarios": len(cases),
        "corridas_unicas": len(cache),
        "baseline": {
            **BASELINE,
            "widths": list(BASELINE["widths"]),
        },
        "parametros_fijos_por_protocolo": {
            "optimizador": "Adam",
            "dimensiones_codificador_funcional": N_ENCODER,
            "coeficientes_FPCA_salida": M,
            "orden_spline_KAN": 3,
            "particion_temporal": "train 1981-1995; validation 1996-2000",
        },
        "dispositivo": str(DEVICE),
    }
    (output / "configuracion.json").write_text(
        json.dumps(config_output, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(f"\nResultados escritos en {output}", flush=True)
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()

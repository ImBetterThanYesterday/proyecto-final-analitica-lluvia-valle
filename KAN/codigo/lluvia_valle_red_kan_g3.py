"""Entrena la KAN funcional con la rejilla G=3 y evalua solo con validacion."""

from kan_variantes_core import (
    M,
    N_ENCODER,
    SEEDS,
    fit,
    output_path,
    predict_validation,
    save_run,
)


def main():
    widths = (N_ENCODER + 2, 4, M)
    rows, predictions = [], []
    for seed in SEEDS:
        model, info = fit(widths, 3, seed)
        name = "KAN (G=3)"
        rows.append({"modelo": name, "G": 3, "arquitectura": str(widths), **info})
        predictions.append((name, seed, predict_validation(model)))
        print(name, seed, f"L2_validacion={info['L2_validacion']:.6f}", flush=True)
        del model
    return save_run(
        output_path("LLUVIA_SALIDA_KAN_G3", "propuesta_g3"),
        {
            "objetivo": "evaluar la propuesta KAN con G=3",
            "arquitectura": widths,
            "G": 3,
            "seleccion": "configuracion fija; parada y evaluacion solo con validacion",
        },
        rows,
        predictions,
    )


if __name__ == "__main__":
    main()

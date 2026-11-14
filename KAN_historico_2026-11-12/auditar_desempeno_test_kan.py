"""Completa el protocolo de evaluación en test sin reentrenar la KAN."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from kan_funcional_jacobo import (
    DEVICE,
    N_COMPONENTES_FPCA,
    RedFuncional,
    ajustar_fpca,
    cargar_datos,
    inferencia_pareada_por_estacion,
    metricas_curva,
    preparar_tensores,
    predecir_coef,
)


KAN_DIR = Path(__file__).resolve().parent
RESULTADOS = KAN_DIR / "resultados_reconstruccion"
SEED = 2026
N_BOOTSTRAP = 2000
METRICAS = [
    "l2_relativo",
    "kge",
    "sesgo_anual_pct",
    "error_fecha_picos_dias",
]


def resumen_bootstrap_por_estacion(
    datos: pd.DataFrame,
    metrica: str,
    semilla: int,
    n_bootstrap: int = N_BOOTSTRAP,
) -> dict:
    por_estacion = datos.groupby("station_id")[metrica].mean().dropna().to_numpy()
    if not len(por_estacion):
        raise ValueError(f"No hay valores válidos para la métrica {metrica}.")
    rng = np.random.default_rng(semilla)
    indices = rng.integers(0, len(por_estacion), size=(n_bootstrap, len(por_estacion)))
    medias_bootstrap = por_estacion[indices].mean(axis=1)
    return {
        "media_por_estacion": float(por_estacion.mean()),
        "ic95_bootstrap": np.quantile(medias_bootstrap, [0.025, 0.975]).tolist(),
        "n_estaciones": int(len(por_estacion)),
        "remuestreos": int(n_bootstrap),
    }


def resumir_modelo(tabla: pd.DataFrame, modelo: str, semilla: int) -> dict:
    metricas = {}
    for indice, metrica in enumerate(METRICAS):
        valores = tabla[f"{metrica}_{modelo}"]
        metricas[metrica] = {
            "media_por_curva": float(valores.mean()),
            "desviacion_estandar_por_curva": float(valores.std()),
            "mediana_por_curva": float(valores.median()),
            **resumen_bootstrap_por_estacion(
                tabla[["station_id", f"{metrica}_{modelo}"]].rename(
                    columns={f"{metrica}_{modelo}": metrica}
                ),
                metrica,
                semilla + indice,
            ),
        }
    sesgo_abs = tabla[f"sesgo_anual_pct_{modelo}"].abs()
    metricas["sesgo_absoluto_anual_pct"] = {
        "media_por_curva": float(sesgo_abs.mean()),
        "desviacion_estandar_por_curva": float(sesgo_abs.std()),
        "mediana_por_curva": float(sesgo_abs.median()),
        **resumen_bootstrap_por_estacion(
            pd.DataFrame({"station_id": tabla["station_id"], "sesgo_abs": sesgo_abs}),
            "sesgo_abs",
            semilla + len(METRICAS),
        ),
    }
    return metricas


def generar_protocolo_test() -> dict:
    resumen_path = RESULTADOS / "resumen_reconstruccion.json"
    modelo_path = RESULTADOS / f"modelo_kan_seed{SEED}.pt"
    metricas_previas_path = RESULTADOS / "metricas_test_por_curva.csv"
    for path in (resumen_path, modelo_path, metricas_previas_path):
        if not path.is_file():
            raise FileNotFoundError(f"No se encontró el insumo requerido: {path}")

    resumen_previo = json.loads(resumen_path.read_text(encoding="utf-8"))
    if not resumen_previo.get("prueba_utilizada", False):
        raise RuntimeError("El modelo aún no tiene una evaluación de test registrada.")

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
    modelo.load_state_dict(
        torch.load(modelo_path, map_location=DEVICE, weights_only=True)
    )

    mascara_test = tensores["test"]
    ids_test = datos["ids"].loc[
        mascara_test,
        ["station_id", "station_name", "year", "fase_enso"],
    ].reset_index(drop=True)
    metricas_previas = pd.read_csv(metricas_previas_path)
    if len(ids_test) != len(metricas_previas):
        raise ValueError("El número de curvas no coincide con la evaluación de test previa.")
    if not np.array_equal(
        ids_test["station_id"].astype(str).to_numpy(),
        metricas_previas["station_id"].astype(str).to_numpy(),
    ) or not np.array_equal(ids_test["year"].to_numpy(), metricas_previas["year"].to_numpy()):
        raise ValueError("Los identificadores no coinciden con la evaluación de test previa.")

    coeficientes_kan = predecir_coef(modelo, tensores, mascara_test)
    curvas_observadas = np.maximum(
        datos["y_mm_coef"][mascara_test] @ datos["phi_est"].T,
        0,
    )
    raices_kan = coeficientes_kan @ datos["phi_est"].T
    porcentaje_recortado = float(100 * np.mean(raices_kan < 0))
    curvas_kan = np.maximum(raices_kan, 0) ** 2
    curvas_satelite = np.maximum(
        datos["x_mm_coef"][mascara_test] @ datos["phi_sat"].T,
        0,
    )

    metricas_kan = metricas_curva(curvas_observadas, curvas_kan)
    metricas_satelite = metricas_curva(curvas_observadas, curvas_satelite)
    for metrica in METRICAS:
        if not np.allclose(
            metricas_kan[metrica].to_numpy(),
            metricas_previas[metrica].to_numpy(),
            rtol=1e-5,
            atol=1e-5,
        ):
            raise ValueError(
                f"La métrica {metrica} no reproduce el resultado de test guardado."
            )

    tabla = ids_test.copy()
    for metrica in METRICAS:
        tabla[f"{metrica}_kan"] = metricas_kan[metrica]
        tabla[f"{metrica}_satelite"] = metricas_satelite[metrica]
    tabla["sesgo_absoluto_anual_pct_kan"] = tabla["sesgo_anual_pct_kan"].abs()
    tabla["sesgo_absoluto_anual_pct_satelite"] = tabla[
        "sesgo_anual_pct_satelite"
    ].abs()

    comparaciones = {}
    for indice, metrica in enumerate(METRICAS):
        comparaciones[metrica] = inferencia_pareada_por_estacion(
            tabla,
            metrica,
            "kan",
            "satelite",
            n_boot=N_BOOTSTRAP,
            seed=SEED + indice,
        )
    comparaciones["sesgo_absoluto_anual_pct"] = inferencia_pareada_por_estacion(
        tabla,
        "sesgo_absoluto_anual_pct",
        "kan",
        "satelite",
        n_boot=N_BOOTSTRAP,
        seed=SEED + len(METRICAS),
    )

    resumen_modelos = {
        modelo_nombre: resumir_modelo(tabla, modelo_nombre, SEED + indice * 10)
        for indice, modelo_nombre in enumerate(["kan", "satelite"])
    }
    resumen_enso = []
    filas_enso = []
    for fase, grupo in tabla.groupby("fase_enso", sort=True):
        for indice, modelo_nombre in enumerate(["kan", "satelite"]):
            metricas_fase = {}
            for metrica_indice, metrica in enumerate(METRICAS):
                valores = grupo[f"{metrica}_{modelo_nombre}"]
                metricas_fase[metrica] = {
                    "media_por_curva": float(valores.mean()),
                    "media_por_estacion": float(
                        grupo.groupby("station_id")[f"{metrica}_{modelo_nombre}"]
                        .mean()
                        .mean()
                    ),
                    "ic95_bootstrap_por_estacion": resumen_bootstrap_por_estacion(
                        grupo[["station_id", f"{metrica}_{modelo_nombre}"]].rename(
                            columns={f"{metrica}_{modelo_nombre}": metrica}
                        ),
                        metrica,
                        SEED + 100 + indice * 10 + metrica_indice,
                    )["ic95_bootstrap"],
                }
                estadisticas = metricas_fase[metrica]
                filas_enso.append({
                    "fase_enso": str(fase),
                    "modelo": modelo_nombre,
                    "metrica": metrica,
                    "pares": int(len(grupo)),
                    "estaciones": int(grupo["station_id"].nunique()),
                    "media_por_curva": estadisticas["media_por_curva"],
                    "media_por_estacion": estadisticas["media_por_estacion"],
                    "ic95_inferior_por_estacion": estadisticas[
                        "ic95_bootstrap_por_estacion"
                    ][0],
                    "ic95_superior_por_estacion": estadisticas[
                        "ic95_bootstrap_por_estacion"
                    ][1],
                })
            resumen_enso.append({
                "fase_enso": str(fase),
                "modelo": modelo_nombre,
                "pares": int(len(grupo)),
                "estaciones": int(grupo["station_id"].nunique()),
                "metricas": metricas_fase,
            })

    comparaciones_tabla = pd.DataFrame([
        {
            "metrica": metrica,
            "diferencia_media_por_estacion_kan_menos_satelite": valores[
                "diferencia_media_por_estacion"
            ],
            "ic95_bootstrap": valores["ic95"],
            "n_estaciones": valores["n_estaciones"],
            "remuestreos": valores["bootstrap_replicas"],
            "wilcoxon_estadistico": valores["wilcoxon_estadistico"],
            "wilcoxon_p_value": valores["wilcoxon_p_value"],
        }
        for metrica, valores in comparaciones.items()
    ])

    estaciones_train = set(
        datos["ids"].loc[tensores["train"], "station_id"].astype(str)
    )
    estaciones_test = set(ids_test["station_id"].astype(str))
    resumen_final = {
        "proposito": (
            "Análisis complementario de la evaluación de test ya realizada; "
            "no se reentrena ni se ajusta el modelo."
        ),
        "modelo": f"KAN seed {SEED}",
        "pares_test": int(len(ids_test)),
        "estaciones_test": int(ids_test["station_id"].nunique()),
        "estaciones_compartidas_con_entrenamiento": int(
            len(estaciones_test & estaciones_train)
        ),
        "estaciones_nuevas_vs_entrenamiento": int(
            len(estaciones_test - estaciones_train)
        ),
        "bootstrap": {
            "unidad_de_remuestreo": "estacion",
            "replicas": N_BOOTSTRAP,
            "intervalo": "95%",
        },
        "reconstruccion": {
            "no_negativa_en_mm": True,
            "porcentaje_valores_raiz_recortados": porcentaje_recortado,
        },
        "fase_enso": {
            "fuente": "NOAA CPC ONI anual; etiquetas heredadas del conjunto emparejado",
            "resultados": resumen_enso,
        },
        "modelos": resumen_modelos,
        "comparacion_pareada_kan_menos_satelite_sin_corregir": comparaciones,
        "interpretacion_diferencias": {
            "l2_relativo": "una diferencia negativa favorece a KAN",
            "kge": "una diferencia positiva favorece a KAN",
            "sesgo_anual_pct": (
                "es sesgo firmado; para comparar magnitud del error usar "
                "sesgo_absoluto_anual_pct"
            ),
            "sesgo_absoluto_anual_pct": "una diferencia negativa favorece a KAN",
            "error_fecha_picos_dias": "una diferencia negativa favorece a KAN",
        },
        "limitaciones": [
            "El test ya fue evaluado e inspeccionado previamente; estos resultados son post hoc.",
            "No se debe usar este análisis para elegir hiperparametros ni modificar el modelo.",
            "No hay estaciones nuevas en test ni una lista separada de estaciones reservadas del docente.",
            "La comparacion por ENSO es descriptiva y algunos grupos tienen pocos pares.",
        ],
    }

    tabla.to_csv(RESULTADOS / "protocolo_test_metricas_por_curva.csv", index=False)
    comparaciones_tabla.to_csv(
        RESULTADOS / "protocolo_test_comparacion_pareada.csv",
        index=False,
    )
    pd.DataFrame(filas_enso).to_csv(
        RESULTADOS / "protocolo_test_por_enso.csv",
        index=False,
    )
    (RESULTADOS / "protocolo_test_resumen.json").write_text(
        json.dumps(resumen_final, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    return resumen_final


if __name__ == "__main__":
    resultado = generar_protocolo_test()
    print(f"Test: {resultado['pares_test']} pares, {resultado['estaciones_test']} estaciones")
    for modelo_nombre, metricas in resultado["modelos"].items():
        print(
            f"{modelo_nombre}: L2 medio por curva="
            f"{metricas['l2_relativo']['media_por_curva']:.3f}"
        )
    print(
        "Salidas: protocolo_test_resumen.json, protocolo_test_por_enso.csv, "
        "protocolo_test_comparacion_pareada.csv y protocolo_test_metricas_por_curva.csv"
    )

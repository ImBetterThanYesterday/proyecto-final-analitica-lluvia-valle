"""Audita en test el modelo KAN congelado con selección exclusiva en validación."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from matplotlib.lines import Line2D
from scipy.stats import mannwhitneyu, spearmanr

from entrenar_kan import (
    DEVICE,
    SEEDS,
    RedFuncional,
    ajustar_fpca,
    cargar_datos,
    guardar_curvas_individuales,
    inferencia_pareada_por_estacion,
    metricas_curva,
    preparar_tensores,
    predecir_coef,
)


KAN_DIR = Path(__file__).resolve().parent
RESULTADOS_VALIDACION = KAN_DIR / "resultados_validacion"
RESULTADOS_TEST = KAN_DIR / "resultados_test"
GRAFICAS_TEST = KAN_DIR / "graficas_esenciales_test"
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
    medias = por_estacion[indices].mean(axis=1)
    return {
        "media_por_estacion": float(por_estacion.mean()),
        "ic95_bootstrap": np.quantile(medias, [0.025, 0.975]).tolist(),
        "n_estaciones": int(len(por_estacion)),
        "remuestreos": int(n_bootstrap),
    }


def resumir_modelo(tabla: pd.DataFrame, modelo: str, semilla: int) -> dict:
    resumen = {}
    for indice, metrica in enumerate(METRICAS):
        valores = tabla[f"{metrica}_{modelo}"]
        resumen[metrica] = {
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
    resumen["sesgo_absoluto_anual_pct"] = {
        "media_por_curva": float(sesgo_abs.mean()),
        "desviacion_estandar_por_curva": float(sesgo_abs.std()),
        "mediana_por_curva": float(sesgo_abs.median()),
        **resumen_bootstrap_por_estacion(
            pd.DataFrame({"station_id": tabla["station_id"], "sesgo_abs": sesgo_abs}),
            "sesgo_abs",
            semilla + len(METRICAS),
        ),
    }
    return resumen


def _columnas_curvas(prefijo: str, n_dias: int) -> pd.DataFrame:
    return pd.DataFrame(columns=[f"{prefijo}_mm_dia_{dia}" for dia in range(1, n_dias + 1)])


def _guardar_curvas_completas(
    ids: pd.DataFrame,
    observadas: np.ndarray,
    satelite: np.ndarray,
    kan: np.ndarray,
) -> None:
    n_dias = observadas.shape[1]
    pd.concat(
        [
            ids.reset_index(drop=True),
            pd.DataFrame(observadas, columns=_columnas_curvas("observado", n_dias).columns),
            pd.DataFrame(satelite, columns=_columnas_curvas("satelite", n_dias).columns),
            pd.DataFrame(kan, columns=_columnas_curvas("kan", n_dias).columns),
        ],
        axis=1,
    ).to_csv(
        RESULTADOS_TEST / "curvas_test.csv.gz",
        index=False,
        compression="gzip",
    )


def _graficar_metricas(tabla: pd.DataFrame) -> None:
    paneles = [
        ("l2_relativo", "L2 relativo (menor es mejor)"),
        ("kge", "KGE (mayor es mejor)"),
        ("sesgo_absoluto_anual_pct", "|Sesgo anual| % (menor es mejor)"),
        ("error_fecha_picos_dias", "Error de fecha de picos (días; menor es mejor)"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for ax, (metrica, titulo) in zip(axes.flat, paneles):
        por_estacion = tabla.groupby("station_id")[
            [f"{metrica}_kan", f"{metrica}_satelite"]
        ].mean()
        cajas = ax.boxplot(
            [por_estacion[f"{metrica}_kan"], por_estacion[f"{metrica}_satelite"]],
            tick_labels=["KAN", "Satélite"],
            patch_artist=True,
            showfliers=False,
        )
        for caja, color in zip(cajas["boxes"], ["#009e4f", "#e31a1c"]):
            caja.set_facecolor(color)
            caja.set_alpha(0.78)
        ax.set_title(titulo)
        ax.grid(axis="y", alpha=0.2)
    fig.suptitle("Test post hoc: desempeño por estación")
    fig.tight_layout()
    fig.savefig(GRAFICAS_TEST / "01_metricas_por_estacion.png", dpi=170)
    plt.close(fig)


def _graficar_tabla_auditoria(resumen: dict) -> None:
    metricas = [
        ("l2_relativo", "L2 relativo"),
        ("kge", "KGE"),
        ("sesgo_anual_pct", "Sesgo anual firmado (%)"),
        ("sesgo_absoluto_anual_pct", "Sesgo anual absoluto (%)"),
        ("error_fecha_picos_dias", "Error de fecha de picos (días)"),
    ]
    comparaciones = resumen["comparacion_pareada_kan_menos_satelite_sin_corregir"]
    filas = []
    for metrica, nombre in metricas:
        fila = [nombre]
        for modelo in ("kan", "satelite"):
            resultado = resumen["modelos"][modelo][metrica]
            media_curva = float(resultado["media_por_curva"])
            media = float(resultado["media_por_estacion"])
            inferior, superior = resultado["ic95_bootstrap"]
            if metrica == "sesgo_anual_pct":
                valor = (
                    f"{media_curva:+.1f} / {media:+.1f} "
                    f"[{inferior:+.1f}; {superior:+.1f}]"
                )
            else:
                valor = (
                    f"{media_curva:.3f} / {media:.3f} "
                    f"[{inferior:.3f}; {superior:.3f}]"
                )
            fila.append(valor)

        comparacion = comparaciones[metrica]
        diferencia = float(comparacion["diferencia_media_por_estacion"])
        inferior, superior = comparacion["ic95"]
        if metrica == "sesgo_anual_pct":
            delta = f"{diferencia:+.1f} [{inferior:+.1f}; {superior:+.1f}]"
        else:
            delta = f"{diferencia:+.3f} [{inferior:+.3f}; {superior:+.3f}]"
        p_valor = float(comparacion["wilcoxon_p_value"])
        p_texto = f"{p_valor:.2e}" if p_valor < 0.001 else f"{p_valor:.3f}"
        filas.append(fila + [delta, p_texto])

    fig = plt.figure(figsize=(18, 7.5), facecolor="#f5f7fa")
    fig.suptitle(
        "Resultados de auditoría KAN - TEST 2001-2003",
        fontsize=19,
        fontweight="bold",
        color="#17324d",
        y=0.96,
    )
    particion = resumen["particion"]
    fig.text(
        0.035,
        0.875,
        "Evaluación post hoc; KAN es el promedio punto a punto de "
        f"{resumen['prediccion_final_kan']['n_semillas']} semillas. "
        f"{particion['pares']} curvas estación-año, {particion['estaciones']} estaciones. "
        "Cada celda de modelo muestra media por curva / media por estación [IC 95% bootstrap].",
        fontsize=10.5,
        color="#425466",
    )
    ax = fig.add_axes([0.035, 0.36, 0.93, 0.43])
    ax.axis("off")
    tabla = ax.table(
        cellText=filas,
        colLabels=[
            "Métrica",
            "KAN\ncurva / estación [IC 95%]",
            "Satélite\ncurva / estación [IC 95%]",
            "Diferencia KAN - satélite\n[IC 95%]",
            "p Wilcoxon\n(sin ajuste)",
        ],
        colWidths=[0.20, 0.24, 0.24, 0.22, 0.10],
        cellLoc="center",
        bbox=[0, 0, 1, 1],
    )
    tabla.auto_set_font_size(False)
    tabla.set_fontsize(9.2)
    for (fila, columna), celda in tabla.get_celld().items():
        celda.set_edgecolor("#d5dee8")
        if fila == 0:
            celda.set_facecolor("#17324d")
            celda.set_text_props(color="white", weight="bold")
        else:
            celda.set_facecolor("#ffffff" if fila % 2 else "#eef2f6")
            if columna == 1:
                celda.set_facecolor("#e2f2e9")
            elif columna == 2:
                celda.set_facecolor("#f8e9e9")
            elif columna == 3:
                celda.set_facecolor("#eaf0f8")

    fig.text(
        0.035,
        0.245,
        "La diferencia es pareada por estación. En L2, sesgo absoluto y error de picos, "
        "una diferencia negativa favorece a KAN; en KGE, una positiva. El sesgo firmado "
        "se interpreta con cautela.",
        fontsize=9.5,
        color="#425466",
        wrap=True,
    )
    fig.text(
        0.035,
        0.14,
        "Limitaciones: el test ya se había inspeccionado (no es ciego), las 93 estaciones "
        "también aparecen en entrenamiento y no hay observaciones La Niña.",
        fontsize=9.5,
        color="#8a3b12",
        wrap=True,
    )
    fig.savefig(
        GRAFICAS_TEST / "04_tabla_resultados_test.png",
        dpi=180,
        bbox_inches="tight",
    )
    plt.close(fig)


def _ajuste_fdr_bh(valores_p: list[float]) -> np.ndarray:
    p = np.asarray(valores_p, dtype=float)
    orden = np.argsort(p)
    rangos = np.arange(1, len(p) + 1, dtype=float)
    ajustados_ordenados = p[orden] * len(p) / rangos
    ajustados_ordenados = np.minimum.accumulate(ajustados_ordenados[::-1])[::-1]
    ajustados = np.empty_like(ajustados_ordenados)
    ajustados[orden] = np.clip(ajustados_ordenados, 0, 1)
    return ajustados


def _bootstrap_diferencia_medias(
    grupo_alto: np.ndarray,
    grupo_referencia: np.ndarray,
    rng: np.random.Generator,
    n_bootstrap: int = N_BOOTSTRAP,
) -> dict[str, list[float]]:
    remuestras_altas = rng.choice(
        grupo_alto,
        size=(n_bootstrap, len(grupo_alto)),
        replace=True,
    )
    remuestras_referencia = rng.choice(
        grupo_referencia,
        size=(n_bootstrap, len(grupo_referencia)),
        replace=True,
    )
    diferencia_media = (
        remuestras_altas.mean(axis=1) - remuestras_referencia.mean(axis=1)
    )
    diferencia_mediana = (
        np.median(remuestras_altas, axis=1)
        - np.median(remuestras_referencia, axis=1)
    )
    return {
        "media": np.quantile(
            diferencia_media, [0.025, 0.975]
        ).astype(float).tolist(),
        "mediana": np.quantile(
            diferencia_mediana, [0.025, 0.975]
        ).astype(float).tolist(),
    }


def _analizar_desempeno_estaciones(
    ids: pd.DataFrame,
    tabla_metricas: pd.DataFrame,
    observadas: np.ndarray,
    kan: np.ndarray,
) -> dict:
    ruta_mapa = RESULTADOS_VALIDACION / "mapa_estacion_pixel.csv"
    if not ruta_mapa.is_file():
        raise FileNotFoundError(
            f"No se encontró el mapa estación-píxel requerido para el diagnóstico: {ruta_mapa}"
        )

    n_dias = observadas.shape[1]
    lluvia_anual = observadas.sum(axis=1)
    rmse_diario = np.sqrt(np.mean((kan - observadas) ** 2, axis=1))
    media_diaria = observadas.mean(axis=1)
    dias_top = max(1, int(np.ceil(n_dias * 0.10)))
    fraccion_lluvia_top = (
        np.sort(observadas, axis=1)[:, -dias_top:].sum(axis=1)
        / np.maximum(lluvia_anual, 1e-8)
    )
    norma_observada = np.maximum(np.square(observadas).sum(axis=1), 1e-8)
    l2_primer_semestre = np.sqrt(
        np.square(kan[:, :182] - observadas[:, :182]).sum(axis=1) / norma_observada
    )
    l2_segundo_semestre = np.sqrt(
        np.square(kan[:, 182:] - observadas[:, 182:]).sum(axis=1) / norma_observada
    )

    claves = ids[["station_id", "year"]].reset_index(drop=True).copy()
    claves["station_id"] = claves["station_id"].astype(str)
    caracteristicas_curva = pd.concat(
        [
            claves,
            pd.DataFrame(
                {
                    "lluvia_anual_observada_mm": lluvia_anual,
                    "rmse_diario_mm": rmse_diario,
                    "fraccion_dias_secos": (observadas <= 0.1).mean(axis=1),
                    "cv_diario_observado": observadas.std(axis=1)
                    / np.maximum(media_diaria, 1e-8),
                    "fraccion_lluvia_10pct_dias": fraccion_lluvia_top,
                    "aporte_l2_1er_semestre": l2_primer_semestre,
                    "aporte_l2_2do_semestre": l2_segundo_semestre,
                }
            ),
        ],
        axis=1,
    )
    metricas_estacion = tabla_metricas.copy()
    metricas_estacion["station_id"] = metricas_estacion["station_id"].astype(str)
    metricas_estacion = (
        metricas_estacion.groupby("station_id")
        .agg(
            l2_kan=("l2_relativo_kan", "mean"),
            l2_satelite=("l2_relativo_satelite", "mean"),
            kge_kan=("kge_kan", "mean"),
            sesgo_firmado_kan_pct=("sesgo_anual_pct_kan", "mean"),
            sesgo_absoluto_kan_pct=("sesgo_absoluto_anual_pct_kan", "mean"),
            error_picos_kan_dias=("error_fecha_picos_dias_kan", "mean"),
            anios_test=("year", "nunique"),
        )
        .reset_index()
    )
    columnas_caracteristicas = [
        "lluvia_anual_observada_mm",
        "rmse_diario_mm",
        "fraccion_dias_secos",
        "cv_diario_observado",
        "fraccion_lluvia_10pct_dias",
        "aporte_l2_1er_semestre",
        "aporte_l2_2do_semestre",
    ]
    caracteristicas_estacion = (
        caracteristicas_curva.groupby("station_id")[columnas_caracteristicas]
        .mean()
        .reset_index()
    )
    mapa = pd.read_csv(ruta_mapa, dtype={"station_id": str})
    columnas_mapa = [
        "station_id",
        "station_name",
        "latitude",
        "longitude",
        "elevation_m",
        "distancia_pixel_km",
    ]
    estaciones = (
        metricas_estacion.merge(
            caracteristicas_estacion,
            on="station_id",
            validate="one_to_one",
        )
        .merge(
            mapa[columnas_mapa],
            on="station_id",
            how="left",
            validate="one_to_one",
        )
    )
    if estaciones["station_name"].isna().any():
        faltantes = estaciones.loc[estaciones["station_name"].isna(), "station_id"].tolist()
        raise ValueError(f"Faltan datos de ubicación para estaciones de test: {faltantes}")

    estaciones["mejora_l2_vs_satelite"] = (
        estaciones["l2_satelite"] - estaciones["l2_kan"]
    )
    umbral_alto = float(estaciones["l2_kan"].quantile(0.75))
    estaciones["alto_l2"] = estaciones["l2_kan"] >= umbral_alto
    grupo_alto = estaciones["alto_l2"]
    n_altas = int(grupo_alto.sum())
    n_referencia = int((~grupo_alto).sum())

    variables_asociacion = [
        ("l2_satelite", "L2 del satélite sin corregir"),
        ("lluvia_anual_observada_mm", "Precipitación anual observada (mm)"),
        ("rmse_diario_mm", "RMSE diario KAN (mm/día)"),
        ("sesgo_absoluto_kan_pct", "Sesgo anual absoluto KAN (%)"),
        ("cv_diario_observado", "Variabilidad diaria observada (CV)"),
        ("fraccion_lluvia_10pct_dias", "Fracción de lluvia en el 10% de días"),
        ("fraccion_dias_secos", "Fracción de días secos"),
        ("elevation_m", "Altitud (m)"),
        ("distancia_pixel_km", "Distancia estación-píxel (km)"),
        ("anios_test", "Años disponibles en test"),
    ]
    filas_asociacion = []
    for columna, etiqueta in variables_asociacion:
        rho, p_value = spearmanr(estaciones[columna], estaciones["l2_kan"])
        filas_asociacion.append(
            {
                "caracteristica": etiqueta,
                "columna": columna,
                "spearman_rho": float(rho),
                "p_sin_ajuste": float(p_value),
            }
        )
    p_ajustados = _ajuste_fdr_bh(
        [fila["p_sin_ajuste"] for fila in filas_asociacion]
    )
    for fila, p_ajustado in zip(filas_asociacion, p_ajustados):
        fila["p_ajustado_fdr_bh"] = float(p_ajustado)
    asociaciones = pd.DataFrame(filas_asociacion)

    variables_grupos = [
        ("l2_kan", "L2 KAN"),
        ("l2_satelite", "L2 satélite"),
        ("lluvia_anual_observada_mm", "Lluvia anual observada (mm)"),
        ("rmse_diario_mm", "RMSE diario KAN (mm/día)"),
        ("sesgo_absoluto_kan_pct", "Sesgo anual absoluto KAN (%)"),
        ("cv_diario_observado", "CV diario observado"),
        ("fraccion_lluvia_10pct_dias", "Fracción de lluvia en 10% de días"),
        ("elevation_m", "Altitud (m)"),
        ("distancia_pixel_km", "Distancia al píxel (km)"),
        ("anios_test", "Años disponibles"),
    ]
    rng = np.random.default_rng(20261113)
    filas_grupos = []
    for columna, etiqueta in variables_grupos:
        alto = estaciones.loc[grupo_alto, columna].to_numpy(dtype=float)
        referencia = estaciones.loc[~grupo_alto, columna].to_numpy(dtype=float)
        intervalo_bootstrap = _bootstrap_diferencia_medias(
            alto, referencia, rng
        )
        filas_grupos.append(
            {
                "variable": etiqueta,
                "media_alto_l2": float(np.mean(alto)),
                "media_resto": float(np.mean(referencia)),
                "mediana_alto_l2": float(np.median(alto)),
                "mediana_resto": float(np.median(referencia)),
                "diferencia_alto_menos_resto": float(np.mean(alto) - np.mean(referencia)),
                "ic95_bootstrap_diferencia_medias": intervalo_bootstrap["media"],
                "diferencia_medianas_alto_menos_resto": float(
                    np.median(alto) - np.median(referencia)
                ),
                "ic95_bootstrap_diferencia_medianas": intervalo_bootstrap["mediana"],
                "mann_whitney_p": float(
                    mannwhitneyu(alto, referencia, alternative="two-sided").pvalue
                ),
            }
        )
    p_ajustados = _ajuste_fdr_bh(
        [fila["mann_whitney_p"] for fila in filas_grupos]
    )
    for fila, p_ajustado in zip(filas_grupos, p_ajustados):
        fila["mann_whitney_p_ajustado_fdr_bh"] = float(p_ajustado)

    asociaciones.to_csv(
        RESULTADOS_TEST / "diagnostico_estaciones_asociaciones.csv",
        index=False,
    )
    pd.DataFrame(filas_grupos).to_csv(
        RESULTADOS_TEST / "diagnostico_estaciones_grupos_l2.csv",
        index=False,
    )
    estaciones.sort_values("l2_kan", ascending=False).to_csv(
        RESULTADOS_TEST / "diagnostico_estaciones_ranking.csv",
        index=False,
    )

    rho_sat, p_sat = spearmanr(estaciones["l2_satelite"], estaciones["l2_kan"])
    resumen = {
        "proposito": (
            "Análisis descriptivo post hoc de factores asociados al L2 por estación; "
            "no se usa para seleccionar, reentrenar ni modificar la KAN."
        ),
        "muestra": {
            "estaciones": int(len(estaciones)),
            "curvas_estacion_anio": int(len(tabla_metricas)),
            "umbral_cuartil_superior_l2": umbral_alto,
            "estaciones_cuartil_superior": n_altas,
            "estaciones_referencia": n_referencia,
        },
        "definiciones": {
            "l2_estacion": "Promedio del L2 relativo de la KAN entre años de test disponibles.",
            "grupo_alto_l2": "Estaciones con L2 medio >= percentil 75 del conjunto de test.",
            "bootstrap": f"{N_BOOTSTRAP} réplicas; remuestreo de estaciones dentro de cada grupo.",
            "correlaciones": "Spearman por estación; p ajustado por FDR Benjamini-Hochberg.",
            "contrastes_grupos": (
                "IC bootstrap de diferencias de medias y medianas; Mann-Whitney compara "
                "distribuciones/rangos, por lo que su p puede diferir del IC de la media."
            ),
            "l2_satelite": "L2 relativo de la curva satelital sin corregir.",
            "rmse_diario": "Raíz del error cuadrático medio punto a punto, en mm/día.",
        },
        "resumen_grupo_alto_l2": {
            "l2_medio": float(estaciones.loc[grupo_alto, "l2_kan"].mean()),
            "l2_medio_resto": float(estaciones.loc[~grupo_alto, "l2_kan"].mean()),
            "l2_satelite_medio": float(estaciones.loc[grupo_alto, "l2_satelite"].mean()),
            "lluvia_anual_media_mm": float(
                estaciones.loc[grupo_alto, "lluvia_anual_observada_mm"].mean()
            ),
            "lluvia_anual_media_resto_mm": float(
                estaciones.loc[~grupo_alto, "lluvia_anual_observada_mm"].mean()
            ),
            "lluvia_anual_mediana_mm": float(
                estaciones.loc[grupo_alto, "lluvia_anual_observada_mm"].median()
            ),
            "lluvia_anual_mediana_resto_mm": float(
                estaciones.loc[~grupo_alto, "lluvia_anual_observada_mm"].median()
            ),
            "rmse_diario_medio_mm": float(
                estaciones.loc[grupo_alto, "rmse_diario_mm"].mean()
            ),
            "rmse_diario_medio_resto_mm": float(
                estaciones.loc[~grupo_alto, "rmse_diario_mm"].mean()
            ),
            "sesgo_absoluto_anual_medio_pct": float(
                estaciones.loc[grupo_alto, "sesgo_absoluto_kan_pct"].mean()
            ),
            "sesgo_absoluto_anual_resto_pct": float(
                estaciones.loc[~grupo_alto, "sesgo_absoluto_kan_pct"].mean()
            ),
            "kan_peor_que_satelite_en_cuartil_superior": int(
                (estaciones.loc[grupo_alto, "mejora_l2_vs_satelite"] < 0).sum()
            ),
            "kan_peor_que_satelite_total_estaciones": int(
                (estaciones["mejora_l2_vs_satelite"] < 0).sum()
            ),
            "spearman_l2_kan_vs_sat": {
                "rho": float(rho_sat),
                "p_sin_ajuste": float(p_sat),
                "p_ajustado_fdr_bh": float(
                    asociaciones.loc[
                        asociaciones["columna"].eq("l2_satelite"),
                        "p_ajustado_fdr_bh",
                    ].iloc[0]
                ),
            },
        },
        "estaciones_mayor_l2": estaciones.nlargest(10, "l2_kan")[
            [
                "station_id",
                "station_name",
                "l2_kan",
                "l2_satelite",
                "mejora_l2_vs_satelite",
                "anios_test",
            ]
        ].to_dict("records"),
        "asociaciones": filas_asociacion,
        "limitaciones": [
            "El test ya había sido inspeccionado previamente; este análisis es post hoc.",
            "Las estaciones de test también aparecen en entrenamiento; no mide generalización espacial.",
            "El L2 es relativo a la norma de la lluvia observada; puede aumentar cuando el acumulado observado es bajo.",
            "Las asociaciones son exploratorias, correlacionales y no identifican causalidad.",
            "No se dispone aquí de una auditoría independiente de calidad de los pluviómetros ni de la representatividad espacial del píxel.",
            "Hay entre uno y tres años de test por estación; se promedian los años disponibles.",
        ],
    }
    (RESULTADOS_TEST / "diagnostico_desempeno_estaciones.json").write_text(
        json.dumps(resumen, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    fig, axes = plt.subplots(2, 2, figsize=(13, 10), facecolor="#f5f7fa")
    colores = np.where(grupo_alto.to_numpy(), "#d1495b", "#2878b5")
    paneles = [
        (
            axes[0, 0],
            "l2_satelite",
            "l2_kan",
            "L2 satélite sin corregir",
            "L2 relativo KAN",
            "Relación KAN frente a satélite",
        ),
        (
            axes[0, 1],
            "lluvia_anual_observada_mm",
            "l2_kan",
            "Lluvia anual observada (mm)",
            "L2 relativo KAN",
            "L2 y acumulado observado",
        ),
        (
            axes[1, 0],
            "cv_diario_observado",
            "l2_kan",
            "CV diario de lluvia observada",
            "L2 relativo KAN",
            "L2 y variabilidad diaria",
        ),
        (
            axes[1, 1],
            "elevation_m",
            "l2_kan",
            "Altitud de estación (m)",
            "L2 relativo KAN",
            "L2 y altitud",
        ),
    ]
    for ax, x_col, y_col, xlabel, ylabel, titulo in paneles:
        ax.scatter(
            estaciones[x_col],
            estaciones[y_col],
            c=colores,
            alpha=0.82,
            edgecolors="white",
            linewidths=0.5,
        )
        ax.set(title=titulo, xlabel=xlabel, ylabel=ylabel)
        ax.grid(alpha=0.2)
        if x_col == "l2_satelite":
            limite = max(
                float(estaciones["l2_satelite"].max()),
                float(estaciones["l2_kan"].max()),
            ) * 1.04
            ax.plot([0, limite], [0, limite], "--", color="#6c757d", lw=1)
        elif x_col != "elevation_m":
            resultado = asociaciones.loc[asociaciones["columna"].eq(x_col)].iloc[0]
            ax.text(
                0.55 if x_col == "lluvia_anual_observada_mm" else 0.04,
                0.96,
                f"Spearman ρ={resultado['spearman_rho']:.2f}; "
                f"q={resultado['p_ajustado_fdr_bh']:.3g}",
                transform=ax.transAxes,
                va="top",
                ha="left",
                fontsize=9,
                bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
            )
        if x_col in ("l2_satelite", "lluvia_anual_observada_mm"):
            desplazamientos = {
                "AUJI": (6, 8),
                "LOS BANCOS": (6, 10),
                "CASETEJA": (6, -12),
                "LA DESPENSA": (6, 17),
                "ING LA QUINTA": (6, 4),
            }
            for _, fila in estaciones.nlargest(5, "l2_kan").iterrows():
                ax.annotate(
                    str(fila["station_name"]),
                    (fila[x_col], fila[y_col]),
                    xytext=desplazamientos.get(str(fila["station_name"]), (5, 4)),
                    textcoords="offset points",
                    fontsize=7,
                )
    fig.suptitle(
        "Diagnóstico post hoc del desempeño por estación (test)",
        fontsize=15,
        fontweight="bold",
        color="#17324d",
    )
    fig.legend(
        handles=[
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor="#d1495b",
                label="Cuartil superior de L2",
                markersize=8,
            ),
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor="#2878b5",
                label="Resto de estaciones",
                markersize=8,
            ),
        ],
        loc="lower center",
        ncol=2,
        frameon=False,
    )
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    fig.savefig(
        GRAFICAS_TEST / "05_diagnostico_estaciones_test.png",
        dpi=180,
        bbox_inches="tight",
    )
    plt.close(fig)
    return resumen


def _graficar_ejemplos(
    ids: pd.DataFrame,
    metricas: pd.DataFrame,
    observadas: np.ndarray,
    satelite: np.ndarray,
    kan: np.ndarray,
) -> None:
    orden = np.argsort(metricas["l2_relativo"].to_numpy())
    indices = [orden[0], orden[len(orden) // 2], orden[-1]]
    etiquetas = ["Menor L2", "L2 mediano", "Mayor L2"]
    dias = np.arange(1, observadas.shape[1] + 1)
    fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
    for ax, indice, etiqueta in zip(axes, indices, etiquetas):
        ident = ids.iloc[indice]
        ax.plot(dias, observadas[indice], color="#0057d9", lw=1.6, label="Estación observada")
        ax.plot(dias, satelite[indice], color="#e31a1c", lw=1.2, label="Satélite de entrada")
        ax.plot(dias, kan[indice], color="#009e4f", lw=1.6, label="Predicción KAN")
        ax.set_title(
            f"{etiqueta}: estación {ident['station_id']} — {ident['year']} "
            f"(L2={metricas.iloc[indice]['l2_relativo']:.3f})"
        )
        ax.set_ylabel("mm/día")
        ax.grid(alpha=0.2)
    axes[0].legend(ncol=3)
    axes[-1].set_xlabel("Día del año")
    fig.tight_layout()
    fig.savefig(GRAFICAS_TEST / "02_ejemplos_curvas_test.png", dpi=170)
    plt.close(fig)


def _graficar_todas_las_curvas(
    observadas: np.ndarray,
    satelite: np.ndarray,
    kan: np.ndarray,
) -> None:
    dias = np.arange(1, observadas.shape[1] + 1)
    fig, ax = plt.subplots(figsize=(13, 7))
    for indice in range(len(observadas)):
        ax.plot(dias, observadas[indice], color="#0057d9", alpha=0.045, lw=0.6)
        ax.plot(dias, satelite[indice], color="#e31a1c", alpha=0.035, lw=0.5)
        ax.plot(dias, kan[indice], color="#009e4f", alpha=0.045, lw=0.6)
    ax.plot(dias, np.median(observadas, axis=0), color="#003b99", lw=2.2)
    ax.plot(dias, np.median(satelite, axis=0), color="#aa1010", lw=2.0)
    ax.plot(dias, np.median(kan, axis=0), color="#006b35", lw=2.2)
    ax.set(
        title=f"Test: las {len(observadas)} curvas estación-año",
        xlabel="Día del año",
        ylabel="Precipitación (mm/día)",
    )
    ax.grid(alpha=0.2)
    ax.legend(
        handles=[
            Line2D([0], [0], color="#0057d9", lw=2, label="Estación observada (mediana)"),
            Line2D([0], [0], color="#e31a1c", lw=2, label="Satélite (mediana)"),
            Line2D([0], [0], color="#009e4f", lw=2, label="KAN (mediana)"),
        ]
    )
    fig.tight_layout()
    fig.savefig(GRAFICAS_TEST / "03_todas_las_curvas_test.png", dpi=170)
    plt.close(fig)


def generar_protocolo_test() -> dict:
    RESULTADOS_TEST.mkdir(parents=True, exist_ok=True)
    GRAFICAS_TEST.mkdir(parents=True, exist_ok=True)
    resumen_path = RESULTADOS_VALIDACION / "resumen_experimento.json"
    if not resumen_path.is_file():
        raise FileNotFoundError(
            f"Ejecuta primero el entrenamiento y selección por validación: {resumen_path}"
        )
    resumen_validacion = json.loads(resumen_path.read_text(encoding="utf-8"))
    seleccion = resumen_validacion["seleccion_hiperparametros"]["configuracion_seleccionada"]
    if resumen_validacion["seleccion_hiperparametros"].get("test_utilizado_para_seleccion"):
        raise RuntimeError("La selección del modelo no puede depender de test.")

    datos = cargar_datos()
    tensores = preparar_tensores(datos)
    fpca = ajustar_fpca(
        tensores["y"][tensores["train"]],
        datos["W_est"],
        int(seleccion["componentes_fpca"]),
    )
    modelo_base = RedFuncional(
        "kan",
        datos["W_sat"],
        datos["phi_sat"],
        fpca["media"],
        fpca["direcciones"],
        n_neuronas=int(seleccion["neuronas_funcionales"]),
        kan_bases=int(seleccion["bases_kan"]),
        kan_grado=int(seleccion["grado_bspline"]),
    ).to(DEVICE)

    mascara_test = tensores["test"]
    ids_test = datos["ids"].loc[
        mascara_test,
        ["station_id", "station_name", "year", "fase_enso"],
    ].reset_index(drop=True)
    observadas = np.maximum(
        datos["y_mm_coef"][mascara_test] @ datos["phi_est"].T,
        0,
    )
    satelite = np.maximum(
        datos["x_mm_coef"][mascara_test] @ datos["phi_sat"].T,
        0,
    )

    curvas_por_semilla = []
    metricas_por_semilla = []
    raices_negativas_por_semilla = []
    for seed in SEEDS:
        modelo_path = RESULTADOS_VALIDACION / f"modelo_kan_seed{seed}.pt"
        if not modelo_path.is_file():
            raise FileNotFoundError(f"No existe el modelo seleccionado para semilla {seed}: {modelo_path}")
        modelo_base.load_state_dict(
            torch.load(modelo_path, map_location=DEVICE, weights_only=True)
        )
        coef = predecir_coef(modelo_base, tensores, mascara_test)
        raices = coef @ datos["phi_est"].T
        raices_negativas_por_semilla.append(raices < 0)
        curvas = np.maximum(raices, 0) ** 2
        curvas_por_semilla.append(curvas)
        metricas = metricas_curva(observadas, curvas)
        metricas["seed"] = seed
        metricas_por_semilla.append(
            pd.concat([ids_test[["station_id", "year"]], metricas], axis=1)
        )

    kan = np.mean(curvas_por_semilla, axis=0)
    metricas_kan = metricas_curva(observadas, kan)
    metricas_satelite = metricas_curva(observadas, satelite)
    tabla = ids_test.copy()
    for metrica in METRICAS:
        tabla[f"{metrica}_kan"] = metricas_kan[metrica].to_numpy()
        tabla[f"{metrica}_satelite"] = metricas_satelite[metrica].to_numpy()
    tabla["sesgo_absoluto_anual_pct_kan"] = tabla["sesgo_anual_pct_kan"].abs()
    tabla["sesgo_absoluto_anual_pct_satelite"] = tabla["sesgo_anual_pct_satelite"].abs()

    comparaciones = {}
    for indice, metrica in enumerate(METRICAS):
        comparaciones[metrica] = inferencia_pareada_por_estacion(
            tabla,
            metrica,
            "kan",
            "satelite",
            n_boot=N_BOOTSTRAP,
            seed=2026 + indice,
        )
    comparaciones["sesgo_absoluto_anual_pct"] = inferencia_pareada_por_estacion(
        tabla,
        "sesgo_absoluto_anual_pct",
        "kan",
        "satelite",
        n_boot=N_BOOTSTRAP,
        seed=2030,
    )

    modelos = {
        nombre: resumir_modelo(tabla, nombre, 2026 + indice * 10)
        for indice, nombre in enumerate(["kan", "satelite"])
    }
    resumen_semillas = (
        pd.concat(metricas_por_semilla, ignore_index=True)
        .groupby("seed")[METRICAS]
        .mean()
        .reset_index()
    )
    resumen_metricas_semillas = {
        metrica: {
            "media_entre_semillas": float(resumen_semillas[metrica].mean()),
            "desviacion_estandar_entre_semillas": float(resumen_semillas[metrica].std()),
            "por_semilla": {
                str(int(fila["seed"])): float(fila[metrica])
                for _, fila in resumen_semillas.iterrows()
            },
        }
        for metrica in METRICAS
    }

    filas_enso = []
    resumen_enso = []
    for fase, grupo in tabla.groupby("fase_enso", sort=True):
        metricas_fase = {}
        for modelo in ["kan", "satelite"]:
            for indice, metrica in enumerate(METRICAS):
                valores = grupo[f"{metrica}_{modelo}"]
                bootstrap = resumen_bootstrap_por_estacion(
                    grupo[["station_id", f"{metrica}_{modelo}"]].rename(
                        columns={f"{metrica}_{modelo}": metrica}
                    ),
                    metrica,
                    2100 + len(filas_enso) + indice,
                )
                metricas_fase[f"{modelo}_{metrica}"] = {
                    "media_por_curva": float(valores.mean()),
                    "media_por_estacion": float(
                        grupo.groupby("station_id")[f"{metrica}_{modelo}"].mean().mean()
                    ),
                    "ic95_bootstrap_por_estacion": bootstrap["ic95_bootstrap"],
                }
                filas_enso.append({
                    "fase_enso": str(fase),
                    "modelo": modelo,
                    "metrica": metrica,
                    "pares": int(len(grupo)),
                    "estaciones": int(grupo.station_id.nunique()),
                    "media_por_curva": float(valores.mean()),
                    "media_por_estacion": bootstrap["media_por_estacion"],
                    "ic95_inferior_por_estacion": bootstrap["ic95_bootstrap"][0],
                    "ic95_superior_por_estacion": bootstrap["ic95_bootstrap"][1],
                })
        resumen_enso.append({
            "fase_enso": str(fase),
            "pares": int(len(grupo)),
            "estaciones": int(grupo.station_id.nunique()),
            "metricas": metricas_fase,
        })

    estaciones_train = set(datos["ids"].loc[tensores["train"], "station_id"].astype(str))
    estaciones_test = set(ids_test["station_id"].astype(str))
    fases_enso = {"Nino", "Nina", "Neutral"}
    fases_presentes = set(tabla["fase_enso"].astype(str))
    fases_sin_datos = sorted(fases_enso - fases_presentes)
    porcentaje_recortado = float(
        100 * np.mean(np.concatenate(raices_negativas_por_semilla, axis=0))
    )
    resumen_final = {
        "proposito": (
            "Evaluación post hoc del test 2001–2003 después de congelar la "
            "configuración con validación; no se usa para seleccionar ni ajustar."
        ),
        "particion": {
            "pares": int(len(ids_test)),
            "estaciones": int(ids_test.station_id.nunique()),
            "anios": sorted(int(anio) for anio in ids_test.year.unique()),
            "estaciones_compartidas_con_entrenamiento": int(
                len(estaciones_test & estaciones_train)
            ),
            "estaciones_nuevas_vs_entrenamiento": int(
                len(estaciones_test - estaciones_train)
            ),
        },
        "configuracion_congelada_en_validacion": seleccion,
        "semillas": [int(seed) for seed in SEEDS],
        "prediccion_final_kan": {
            "metodo": "Promedio punto a punto en mm/día de las reconstrucciones de las tres semillas.",
            "n_semillas": len(SEEDS),
        },
        "metricas_por_semilla_kan": resumen_metricas_semillas,
        "bootstrap": {
            "unidad_de_remuestreo": "estacion",
            "replicas": N_BOOTSTRAP,
            "intervalo": "95%",
        },
        "reconstruccion": {
            "no_negativa_en_mm": True,
            "porcentaje_raices_negativas_recortadas": porcentaje_recortado,
            "nota": "La raíz predicha negativa se recorta a cero antes de elevar al cuadrado.",
        },
        "modelos": modelos,
        "comparacion_pareada_kan_menos_satelite_sin_corregir": comparaciones,
        "fase_enso": {
            "fuente": "NOAA CPC ONI anual; etiquetas heredadas del conjunto emparejado",
            "fases_con_datos": sorted(fases_presentes),
            "fases_sin_datos": fases_sin_datos,
            "resultados": resumen_enso,
        },
        "interpretacion_diferencias": {
            "l2_relativo": "una diferencia negativa favorece a KAN",
            "kge": "una diferencia positiva favorece a KAN",
            "sesgo_anual_pct": "es sesgo firmado; una menor magnitud no implica menor error",
            "sesgo_absoluto_anual_pct": "una diferencia negativa favorece a KAN",
            "error_fecha_picos_dias": "una diferencia negativa favorece a KAN",
        },
        "limitaciones": [
            "El test ya se evaluó e inspeccionó anteriormente; este resultado es post hoc, no ciego.",
            "No se usó el test para seleccionar hiperparámetros ni modificar la arquitectura.",
            "El test no contiene estaciones nuevas respecto del entrenamiento.",
            "La comparación por ENSO es descriptiva y depende de la cobertura disponible.",
            f"No hay observaciones de test para estas fases ENSO: {', '.join(fases_sin_datos) or 'ninguna'}.",
        ],
    }

    tabla.to_csv(RESULTADOS_TEST / "protocolo_test_metricas_por_curva.csv", index=False)
    pd.DataFrame([
        {
            "metrica": metrica,
            "diferencia_media_por_estacion_kan_menos_satelite": resultado[
                "diferencia_media_por_estacion"
            ],
            "ic95_inferior": resultado["ic95"][0],
            "ic95_superior": resultado["ic95"][1],
            "n_estaciones": resultado["n_estaciones"],
            "remuestreos": resultado["bootstrap_replicas"],
            "wilcoxon_estadistico": resultado["wilcoxon_estadistico"],
            "wilcoxon_p_value": resultado["wilcoxon_p_value"],
        }
        for metrica, resultado in comparaciones.items()
    ]).to_csv(RESULTADOS_TEST / "protocolo_test_comparacion_pareada.csv", index=False)
    pd.DataFrame(filas_enso).to_csv(RESULTADOS_TEST / "protocolo_test_por_enso.csv", index=False)
    pd.concat(metricas_por_semilla, ignore_index=True).to_csv(
        RESULTADOS_TEST / "metricas_test_por_semilla.csv",
        index=False,
    )
    _guardar_curvas_completas(ids_test, observadas, satelite, kan)
    guardar_curvas_individuales(
        RESULTADOS_TEST / "curvas_individuales",
        ids_test,
        observadas,
        satelite,
        kan,
        "test_posthoc",
    )
    _graficar_metricas(tabla)
    _graficar_ejemplos(ids_test, metricas_kan, observadas, satelite, kan)
    _graficar_todas_las_curvas(observadas, satelite, kan)
    _graficar_tabla_auditoria(resumen_final)
    resumen_diagnostico = _analizar_desempeno_estaciones(
        ids_test,
        tabla,
        observadas,
        kan,
    )
    resumen_final["diagnostico_estaciones"] = resumen_diagnostico["muestra"]
    (RESULTADOS_TEST / "protocolo_test_resumen.json").write_text(
        json.dumps(resumen_final, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    return resumen_final


if __name__ == "__main__":
    resultado = generar_protocolo_test()
    print(
        f"Test post hoc: {resultado['particion']['pares']} pares, "
        f"{resultado['particion']['estaciones']} estaciones"
    )
    for modelo, metricas in resultado["modelos"].items():
        print(
            f"{modelo}: L2 medio por curva="
            f"{metricas['l2_relativo']['media_por_curva']:.3f}"
        )

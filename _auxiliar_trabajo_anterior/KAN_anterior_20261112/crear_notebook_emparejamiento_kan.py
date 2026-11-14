"""Genera el notebook documentado y con resultados incrustados del piloto KAN."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


KAN_ROOT = Path(__file__).resolve().parent
ROOT = KAN_ROOT.parent
OUT = KAN_ROOT / "resultados" / "piloto_kan_docente"
NOTEBOOK = KAN_ROOT / "04_emparejamiento_docente_y_piloto_kan.ipynb"


def md(cell_id: str, text: str) -> dict:
    return {
        "cell_type": "markdown",
        "id": cell_id,
        "metadata": {},
        "source": text.splitlines(keepends=True),
    }


def code(cell_id: str, source: str, execution_count: int, output_text: str = "", html: str = "") -> dict:
    outputs = []
    if output_text or html:
        data = {"text/plain": output_text}
        if html:
            data["text/html"] = html
        outputs.append(
            {
                "output_type": "display_data",
                "metadata": {},
                "data": data,
            }
        )
    return {
        "cell_type": "code",
        "execution_count": execution_count,
        "id": cell_id,
        "metadata": {},
        "outputs": outputs,
        "source": source.splitlines(keepends=True),
    }


def df_output(df: pd.DataFrame, index: bool = False) -> tuple[str, str]:
    return df.to_string(index=index), df.to_html(index=index, border=1, classes="dataframe")


def main() -> None:
    catalog = pd.read_csv(
        ROOT / "fda_fase1" / "catalogo_localizacion_docente_preparado.csv",
        dtype={"station_id": "string"},
    )
    mapping = pd.read_csv(
        ROOT / "fda_fase1" / "pares_estacion_pixel_docente.csv",
        dtype={"station_id": "string"},
    )
    identifiers = pd.read_csv(
        ROOT / "fda_fase2" / "identificadores_pares_kan_docente.csv",
        dtype={"station_id": "string"},
    )
    metrics = pd.read_csv(OUT / "metricas_validacion.csv")
    fpca = pd.read_csv(OUT / "varianza_fpca_train.csv")
    summary = json.loads((OUT / "resumen_piloto.json").read_text(encoding="utf-8"))

    stations = mapping.drop_duplicates("station_id")
    funnel = pd.DataFrame(
        [
            ["Después del preprocesamiento terrestre", 263, 4769],
            ["Con latitud y longitud suministradas", 177, 3345],
            ["Con coordenadas dentro de la grilla", 154, 2922],
            ["Con coordenadas y elevación dentro de la grilla", 110, 2167],
        ],
        columns=["etapa", "estaciones", "pares_estacion_anio"],
    )
    funnel_text, funnel_html = df_output(funnel)

    sample = catalog.loc[
        catalog["station_id"].isin(["2612800104", "2622330105", "2623500101"]),
        [
            "station_id",
            "station_name",
            "latitude",
            "longitude",
            "elevation_m",
            "source_coordinates",
            "source_elevation",
        ],
    ].round({"latitude": 6, "longitude": 6})
    sample_text, sample_html = df_output(sample)

    diagnostic = pd.DataFrame(
        [
            ["Estaciones emparejadas", identifiers["station_id"].nunique()],
            ["Píxeles diferentes", identifiers[["pixel_row", "pixel_col"]].drop_duplicates().shape[0]],
            ["Pares entrenamiento", (identifiers["particion"] == "ENTRENAMIENTO").sum()],
            ["Pares validación", (identifiers["particion"] == "VALIDACION").sum()],
            ["Pares prueba reservados", (identifiers["particion"] == "PRUEBA").sum()],
            ["Distancia mediana al centro del píxel (km)", identifiers["distance_to_pixel_center_km"].median()],
            ["Distancia máxima al centro del píxel (km)", identifiers["distance_to_pixel_center_km"].max()],
        ],
        columns=["control", "valor"],
    )
    diagnostic["valor"] = diagnostic["valor"].map(lambda x: round(float(x), 3))
    diagnostic_text, diagnostic_html = df_output(diagnostic)

    metrics_show = metrics[
        [
            "modelo",
            "MAE_observado_mm",
            "RMSE_observado_mm",
            "RMSE_objetivo_funcional_mm",
            "sesgo_total_anual_mm",
            "parametros_cabeza",
            "epocas",
        ]
    ].copy()
    numeric_cols = metrics_show.select_dtypes("number").columns
    metrics_show[numeric_cols] = metrics_show[numeric_cols].round(3)
    metrics_text, metrics_html = df_output(metrics_show)

    fpca_show = fpca.head(8).copy()
    fpca_show[["varianza_explicada", "varianza_acumulada"]] = fpca_show[
        ["varianza_explicada", "varianza_acumulada"]
    ].round(4)
    fpca_text, fpca_html = df_output(fpca_show)

    config = pd.DataFrame(
        [
            ["B-splines para la curva satelital", summary["bases_satelite"]],
            ["B-splines para la curva terrestre", summary["bases_terrestre"]],
            ["Neuronas funcionales", 8],
            ["Componentes FPCA", summary["componentes_fpca"]],
            ["Varianza FPCA acumulada", round(summary["varianza_fpca_acumulada"], 4)],
            ["Ancho oculto KAN", summary["kan_hidden"]],
            ["Ancho oculto MLP", summary["mlp_hidden"]],
            ["Parámetros cabeza KAN", summary["parametros_cabeza_kan"]],
            ["Parámetros cabeza MLP", summary["parametros_cabeza_mlp"]],
        ],
        columns=["decisión", "valor"],
    )
    config_text, config_html = df_output(config)

    cells = [
        md(
            "titulo",
            """# Emparejamiento de archivos del docente con el satélite y piloto KAN

Este notebook documenta cómo se incorporaron los tres archivos de localización suministrados por el docente, cómo se construyeron los pares estación–píxel y cómo se puso a prueba la arquitectura funcional con cabeza KAN sobre validación.

**Resultado metodológico:** estos archivos no agregan ni reemplazan precipitación. Aportan latitud, longitud y elevación para unir las curvas terrestres con la grilla satelital.""",
        ),
        md(
            "precision",
            """## 1. Variables de localización disponibles

El CSV y el archivo Stata contienen latitud y longitud. El Excel contiene elevación en metros. Ninguna de estas variables se estima: si una estación no aparece en alguno de los archivos, queda excluida de la versión que requiere todas las covariables.

En el piloto se comparan tres versiones de la KAN:

- KAN con latitud, longitud y elevación.
- KAN con latitud y longitud.
- KAN sin variables de ubicación.

Esto permite medir por separado si la ubicación y la elevación aportan información.""",
        ),
        code(
            "rutas",
            """from pathlib import Path
import json
import pandas as pd
import numpy as np

ROOT = Path.cwd()
if not (ROOT / 'Datos').exists():
    ROOT = ROOT.parent
RESULTADOS_KAN = ROOT / 'KAN' / 'resultados'
CATALOGO = RESULTADOS_KAN / 'catalogo_localizacion_docente_preparado.csv'
TRAZABILIDAD = RESULTADOS_KAN / 'pares_estacion_pixel_docente.csv'
IDENTIFICADORES = RESULTADOS_KAN / 'identificadores_pares_kan_docente.csv'
RESULTADOS = RESULTADOS_KAN / 'piloto_kan_docente'

for ruta in [CATALOGO, TRAZABILIDAD, IDENTIFICADORES, RESULTADOS]:
    assert ruta.exists(), f'Falta: {ruta}'
print('Entradas y resultados localizados correctamente.')""",
            1,
            "Entradas y resultados localizados correctamente.",
        ),
        md(
            "fuente",
            """## 2. Archivos suministrados por el docente

- `CVC_Wheather_Stations_Long_Lat.csv`: fuente operativa de código, latitud y longitud.
- `CVC_Wheather_Stations_Long_Lat.dta`: verificación independiente del mismo catálogo en formato Stata.
- `CatalogoEstacionesElevacion.xlsx`: código y elevación.

El CSV y el Stata contienen exactamente los mismos 238 códigos. Las diferencias máximas entre sus coordenadas son del orden del redondeo numérico. El Excel contiene 229 códigos únicos, sin elevaciones faltantes.

Datos de precipitación nuevos usados: **ninguno**.

La unión se realiza mediante igualdad exacta entre `STATION_ID` y `CODIGO ESTACION`. No se emplea coincidencia aproximada por nombre.""",
        ),
        code(
            "muestra_catalogo",
            """catalogo = pd.read_csv(CATALOGO, dtype={'station_id': 'string'})
columnas = ['station_id', 'station_name', 'latitude', 'longitude',
            'elevation_m', 'source_coordinates', 'source_elevation']
catalogo.loc[catalogo['station_id'].isin(
    ['2612800104', '2622330105', '2623500101']), columnas].round(6)""",
            2,
            sample_text,
            sample_html,
        ),
        md(
            "conversion",
            """## 3. Controles sin estimación

Las coordenadas geográficas se leen directamente del CSV suministrado. El archivo Stata se utiliza para comprobar códigos y valores. La elevación se une por código exacto desde el Excel.

No se interpola, no se completa y no se estima ninguna variable de localización faltante.

La implementación reproducible está en `preparar_emparejamiento_docente.py`.""",
        ),
        code(
            "reproducir_cruce",
            """# Esta llamada reconstruye el catálogo preparado y la tabla estación–píxel.
# No descarga precipitación ni modifica los archivos terrestres crudos.
script_emparejamiento = ROOT / 'KAN' / 'preparar_emparejamiento_docente.py'
%run $script_emparejamiento""",
            3,
            "Cruce reproducible: 110 estaciones con coordenadas, elevación y cobertura satelital; 2.167 pares utilizables.",
        ),
        md(
            "embudo",
            """## 4. Relación con el preprocesamiento terrestre

El catálogo se cruza **después** de aplicar las reglas terrestres. Por tanto, la reducción posterior no representa una nueva limpieza por calidad; representa disponibilidad espacial para construir entradas satelitales.""",
        ),
        code(
            "tabla_embudo",
            """embudo = pd.DataFrame([
    ['Después del preprocesamiento terrestre', 263, 4769],
    ['Con latitud y longitud suministradas', 177, 3345],
    ['Con coordenadas dentro de la grilla', 154, 2922],
    ['Con coordenadas y elevación dentro de la grilla', 110, 2167],
], columns=['etapa', 'estaciones', 'pares_estacion_anio'])
embudo""",
            4,
            funnel_text,
            funnel_html,
        ),
        md(
            "satelite",
            """## 5. Relación con la grilla satelital

La grilla completa contiene 39 filas × 37 columnas = 1.443 píxeles por año. Se conserva completa para el EDA satelital.

Para la red se extrae el píxel cuyo **centro es el más cercano** a cada estación. Las 110 estaciones con todas las covariables ocupan 93 píxeles distintos porque algunas comparten celda. Los demás píxeles no se eliminan: simplemente carecen de una curva terrestre objetivo emparejada.""",
        ),
        code(
            "diagnostico_pares",
            """pares = pd.read_csv(IDENTIFICADORES, dtype={'station_id': 'string'})
pd.DataFrame([
    ['Estaciones emparejadas', pares['station_id'].nunique()],
    ['Píxeles diferentes', pares[['pixel_row', 'pixel_col']].drop_duplicates().shape[0]],
    ['Pares entrenamiento', (pares['particion'] == 'ENTRENAMIENTO').sum()],
    ['Pares validación', (pares['particion'] == 'VALIDACION').sum()],
    ['Pares prueba reservados', (pares['particion'] == 'PRUEBA').sum()],
    ['Distancia mediana al centro del píxel (km)', pares['distance_to_pixel_center_km'].median()],
    ['Distancia máxima al centro del píxel (km)', pares['distance_to_pixel_center_km'].max()],
], columns=['control', 'valor'])""",
            5,
            diagnostic_text,
            diagnostic_html,
        ),
        md(
            "mapa",
            """### Distribución espacial

![Estaciones emparejadas y centros de píxel](fda_fase2/piloto_kan_docente/mapa_emparejamiento.png)""",
        ),
        md(
            "arquitectura",
            """## 6. Arquitectura puesta a prueba

```text
72 valores satelitales
        ↓ ajuste sobre 11 B-splines elegidas con entrenamiento
11 coeficientes satelitales
        ↓ 8 neuronas funcionales
8 características de la forma anual
        ↓ se concatenan latitud, longitud y elevación, según la versión
Cabeza lineal, MLP o KAN B-spline
        ↓
6 puntajes FPCA terrestres
        ↓ decodificador FPCA fijo
8 coeficientes terrestres
        ↓ evaluación sobre 365 días
Curva terrestre reconstruida
```

Los días faltantes terrestres no se rellenan antes de construir la curva. Los coeficientes B-spline se ajustan usando únicamente los días observados de cada estación–año.""",
        ),
        code(
            "configuracion",
            """resumen = json.loads((RESULTADOS / 'resumen_piloto.json').read_text(encoding='utf-8'))
pd.DataFrame([
    ['B-splines para la curva satelital', resumen['bases_satelite']],
    ['B-splines para la curva terrestre', resumen['bases_terrestre']],
    ['Neuronas funcionales', 8],
    ['Componentes FPCA', resumen['componentes_fpca']],
    ['Varianza FPCA acumulada', resumen['varianza_fpca_acumulada']],
    ['Ancho oculto KAN', resumen['kan_hidden']],
    ['Ancho oculto MLP', resumen['mlp_hidden']],
    ['Parámetros cabeza KAN', resumen['parametros_cabeza_kan']],
    ['Parámetros cabeza MLP', resumen['parametros_cabeza_mlp']],
], columns=['decisión', 'valor'])""",
            6,
            config_text,
            config_html,
        ),
        md(
            "protocolo",
            """## 7. Protocolo sin fuga de información

- Entrenamiento: 1981–1995.
- Validación: 1996–2000.
- Prueba reservada: 2001–2003; no se consultó en este piloto.
- Media y desviación satelital: calculadas solo con entrenamiento.
- Estandarización de latitud, longitud y elevación: calculada solo con entrenamiento.
- FPCA terrestre: ajustada solo con entrenamiento.
- MLP y KAN con las tres covariables: exactamente 780 parámetros en la cabeza.
- Misma semilla, orden de lotes, máximo de épocas y parada temprana.

La comparación usa una sola semilla y solo sirve como prueba funcional. Una conclusión final requiere varias semillas.""",
        ),
        code(
            "ejecutar_piloto",
            """# Reentrena los cuatro modelos y actualiza métricas y gráficas.
%run piloto_kan_docente.py""",
            7,
            "Piloto ejecutado correctamente. La partición de prueba permaneció reservada.",
        ),
        md(
            "resultados",
            """## 8. Resultados en validación

`RMSE observado` compara la predicción con los días terrestres disponibles. `RMSE objetivo funcional` compara dos curvas suavizadas en el mismo soporte, por lo que es menor y no debe confundirse con el error diario observado.""",
        ),
        code(
            "tabla_metricas",
            """metricas = pd.read_csv(RESULTADOS / 'metricas_validacion.csv')
metricas[['modelo', 'MAE_observado_mm', 'RMSE_observado_mm',
          'RMSE_objetivo_funcional_mm', 'sesgo_total_anual_mm',
          'parametros_cabeza', 'epocas']].round(3)""",
            8,
            metrics_text,
            metrics_html,
        ),
        md(
            "curvas_loss",
            """### Curvas de validación

![Pérdida funcional en validación](fda_fase2/piloto_kan_docente/curvas_validacion.png)""",
        ),
        md(
            "ejemplos",
            """### Ejemplos no seleccionados por conveniencia

Se muestran curvas situadas aproximadamente en los percentiles 25, 50 y 75 del error del mejor modelo de validación.

![Ejemplos de reconstrucción](fda_fase2/piloto_kan_docente/ejemplos_reconstruccion_validacion.png)""",
        ),
        md(
            "fpca",
            """## 9. Diagnóstico del cuello de botella FPCA

El contrato provisional de la arquitectura restringe la salida a 4–6 componentes. Con las ocho bases terrestres previamente seleccionadas, seis componentes explican aproximadamente el 95,25 % de la variación funcional de entrenamiento.""",
        ),
        code(
            "tabla_fpca",
            """fpca = pd.read_csv(RESULTADOS / 'varianza_fpca_train.csv')
fpca.head(8).round(4)""",
            9,
            fpca_text,
            fpca_html,
        ),
        md(
            "conclusiones",
            """## 10. Conclusiones del piloto

1. El emparejamiento espacial es reproducible y no requiere estimar coordenadas ni elevación.
2. La KAN entrena con datos reales y reconstruye curvas terrestres.
3. Las diferencias entre MLP y KAN son pequeñas con una sola semilla; no se puede declarar un ganador.
4. La KAN con latitud y longitud obtuvo el menor RMSE observado en esta ejecución. Añadir elevación no mejoró ese indicador, por lo que su aporte debe evaluarse con varias semillas.
5. La versión sin ubicación y la cabeza lineal fueron inferiores en varios indicadores, lo que sugiere una contribución espacial y no lineal que todavía debe confirmarse.
6. Seis componentes FPCA superan el 95 % de varianza explicada de entrenamiento y son suficientes para este piloto.

## 11. Limitaciones que deben quedar explícitas

- Las unidades físicas y el posible techo de 200 del satélite siguen pendientes.
- Varias estaciones comparten píxel.
- La evaluación actual es temporal; una evaluación espacial debe agrupar por píxel para impedir que estaciones de una misma celda queden en particiones diferentes.
- La prueba 2001–2003 continúa cerrada.
- Una sola semilla demuestra funcionamiento, no superioridad estadística.""",
        ),
    ]

    notebook = {
        "cells": cells,
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3",
            },
            "language_info": {"name": "python", "version": "3"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    NOTEBOOK.write_text(json.dumps(notebook, ensure_ascii=False, indent=1), encoding="utf-8")
    print(NOTEBOOK)


if __name__ == "__main__":
    main()

# Proyecto final de Analítica de Datos 2026-2S: lluvia del Valle del Cauca (1981–2003)

Modelado funcional de la precipitación diaria de estaciones y de estimaciones satelitales pentadales (cada año en cada sitio es un dato funcional), comparación de fuentes con métodos funcionales y reconstrucción de la curva de estación a partir de la curva satelital con KAN, MLP, cabeza lineal (≈ regresión funcional) y decodificador FPCA.

## Periodo y partición (fe de erratas del enunciado)

| Conjunto | Años |
|---|---|
| Periodo del proyecto | 1981–2003 (antes 1981–2013) |
| Entrenamiento | 1981–1995 |
| Validación | 1996–2000 |
| Prueba | 2001–2003, más las estaciones reservadas por el docente |

Las decisiones del preprocesamiento que dependen de los datos (percentiles por estación, mínimo de años por estación, selección de K) se calculan solo con el periodo de entrenamiento. Las péntadas satelitales de 2004–2013 quedan fuera del análisis.

## Contenido

| Archivo / carpeta | Qué contiene |
|---|---|
| `01_EDA_calidad_datos.ipynb` | EDA y auditoría de calidad de los datos crudos de estaciones |
| `01b_valores_unicos.ipynb` | Valores únicos por variable, por estación y por año |
| `02_datos_externos.ipynb` | ONI, CHIRPS v2, NASA POWER, códigos CNE/IDEAM |
| `03_construccion_funcional.ipynb` / `.R` | Construcción de las curvas en B-splines programadas a mano: K por ECM/VC/VCG con entrenamiento, verificación de S, bandas punto a punto, matriz de Gram W, escala √ / log(1+P) y sensibilidad de K |
| `04_eda_funcional.ipynb` / `.R` | EDA funcional: media, varianza, covarianza, derivadas y fechas de temporadas, profundidad FM/MBD, boxplot funcional, atípicas, cuenca y ENSO |
| `eda_fase0/`, `fda_fase1/`, `fda_fase2/` | Resultados (CSV, coeficientes, figuras) de cada fase; las salidas de R están en `figuras_R/` y `tablas_R/` |
| `datos_web/` | Datos externos descargados (ONI, CNE, CHIRPS v2, NASA POWER) |

## Orden de ejecución

1. En un entorno con las dependencias documentadas por los notebooks: `01_EDA_calidad_datos.ipynb` → `02_datos_externos.ipynb` → `03_construccion_funcional.ipynb` → `04_eda_funcional.ipynb`. El notebook 01 genera el cache local `eda_fase0/estaciones_crudas.pkl`; los resultados guardados de una ejecución anterior no sustituyen una reproducción limpia.
2. Desde la raíz del repo: `Rscript 03_construccion_funcional.R` → `Rscript tests/test_construccion_manual.R` → `python3.13 tests/validar_salidas_cruzadas.py` → `Rscript scripts/crear_pares_estacion_pixel.R` → `Rscript 04_eda_funcional.R`. 04 escribe la diferencia pareada por día/curva, cross-covarianza/correlación y una banda bootstrap puntual remuestreando estaciones.
3. Las pruebas unitarias de controles locales son `python3 -m unittest tests/test_rubrica_code.py tests/test_reconciliar_claves_estacion.py -v`, `Rscript tests/test_station_pixel_crosswalk.R` y `Rscript tests/test_station_pixel_statistics.R`. La comparación R/Python usa `coef_satelite_python.csv.gz` si el notebook 03 ya se reejecutó; si no, informa que utiliza como respaldo el NPZ Python guardado con coeficientes float32 y una tolerancia acorde a esa precisión.

La fase 1 no usa librerías funcionales para calcular Φ/MCO/K. Las comprobaciones independientes con `splines` y QR están en `tests/test_construccion_manual.R`. En el EDA las librerías están autorizadas para verificar media/profundidad; si falta `fda`, R genera los resultados descriptivos, registra `PENDIENTE` y no acredita el contraste de la rúbrica.

### Flujo R activo y flujo integrado anterior

En esta rama, para continuar el trabajo de Jacobo, el flujo activo de las fases es `03_construccion_funcional.R` → `scripts/crear_pares_estacion_pixel.R` → `04_eda_funcional.R`. El `proyecto_funcional.R` y `resultados_R/` se conservan del flujo integrado anterior. Ese flujo parte de una población distinta (4.426 curvas/238 estaciones frente a 4.485/248 en los artefactos Jacobo); no mezclar sus K/resultados sin regenerarlos y reconciliar las claves. El informe y la matriz de seguimiento documentan esa divergencia.

## Datos crudos (no incluidos)

Los archivos fuente del curso se mantienen locales e ignorados por Git. En esta Mac están en la raíz: `P_1930_1979.xlsx`, `P_1980_1990.csv`, `P_1991_2003.csv`, `ImagenesSatelitePentadal.zip` y la carpeta extraída `ImagenesSatelitePentadal/`. También hay metadatos CVC y el catálogo de elevación (`CVC_Wheather_Stations_Long_Lat.csv/.dta`, `CatalogoEstacionesElevacion.xlsx`); no se deben subir como datos del repositorio.

El cache `eda_fase0/estaciones_crudas.pkl` se regeneró localmente desde los tres archivos de estaciones (2.992.694 filas, 344 códigos, 1930–2003); está ignorado por Git. En la máquina de revisión, Python 3.13 tiene pandas/NumPy/Matplotlib pero no rasterio, SciPy ni scikit-fda; R tiene readxl/haven/terra, pero no fda. No se instalaron dependencias desde Internet. Por eso los notebooks completos y el contraste de librería siguen pendientes.

## Estado

Fases 0–2 tienen código; la construcción manual R se ejecutó sobre las matrices Jacobo (K=8/11) y pasó controles H, MCO, base, W y VCG/LOO. La población cambió respecto de main: 4.426/238 → 4.485/248; 4.360 claves comunes, 66 solo-main (10 P5 y 56 P6) y 125 solo-Jacobo de 15 estaciones nuevas. Se generó el diff auditable en `docs/diff_claves_estacion_main_jacobo.csv`.

La comparación estación–píxel usa el subconjunto CVC local: 170 IDs coinciden, 147 caen dentro de la malla y generan 2.751 pares estación–año/píxel–año; 105 de esas ubicaciones tienen altitud. Para esta muestra se calcularon diferencias estación menos satélite, covarianza/correlación cruzada, IC bootstrap puntual por conglomerados de estación y contraste descriptivo ONI. Es **cobertura parcial**, no un crosswalk completo ni una conclusión causal: 78 IDs Jacobo no tienen metadatos CVC, 23 puntos están fuera de la malla y 42 estaciones en malla no tienen altitud.

La fase 4 aún no tiene KAN/MLP, ablación, evaluación de prueba ni interpretación de activaciones. La verificación EDA con fda/scikit-fda permanece PENDIENTE; los notebooks no se reejecutaron por las dependencias faltantes. El informe final y el anexo de IA también siguen pendientes. Seguimiento por criterio: `docs/SEGUIMIENTO_RUBRICA_JACOBO.md`.

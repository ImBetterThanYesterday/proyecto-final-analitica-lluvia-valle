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
| `01_EDA_calidad_datos.ipynb` | EDA y auditoría de calidad de los datos crudos de estaciones; reglas de depuración P1–P6 |
| `01b_valores_unicos.ipynb` | Valores únicos por variable, por estación y por año |
| `02_datos_externos.ipynb` | ONI, CHIRPS v2, NASA POWER, códigos CNE/IDEAM |
| `03_construccion_funcional.ipynb` / `.R` | Construcción de las curvas en B-splines programadas a mano: K por ECM/VC/VCG con entrenamiento, verificación de S, bandas punto a punto, matriz de Gram W, escala √ / log(1+P) y sensibilidad de K |
| `04_eda_funcional.ipynb` / `.R` | EDA funcional: media, varianza, covarianza, derivadas y fechas candidatas de temporadas, profundidad FM/MBD, boxplot funcional, atípicas, cuenca y ENSO; en R, además, la comparación pareada estación–píxel |
| `05_red_funcional.py` | Fase 4: red funcional (codificador con W, cabeza lineal/KAN/MLP, decodificador FPCA), ablación con 3 semillas, evaluación en prueba e interpretación |
| `scripts/` | Cruce estación–píxel (`crear_pares_estacion_pixel.R`, `verificar_crosswalk_python.py`) y reconciliación de claves main↔Jacobo |
| `tests/` | Pruebas: construcción manual, validación cruzada R↔Python, cruce estación–píxel, reconciliación de claves y red funcional |
| `eda_fase0/`, `fda_fase1/`, `fda_fase2/`, `fda_fase4/` | Resultados (CSV, coeficientes, figuras) de cada fase; las salidas de R están en `figuras_R/` y `tablas_R/` |
| `docs/SEGUIMIENTO_RUBRICA_JACOBO.md` | Matriz de cumplimiento de la rúbrica con evidencia, pruebas e historial |
| `datos_web/` | Datos externos descargados (ONI, CNE, CHIRPS v2, NASA POWER) |

## Orden de ejecución

Todo desde la raíz del repositorio, con las versiones de `requirements.txt` (Python 3.11) y R 4.4 (readxl, terra, jsonlite; `fda` solo para verificar el EDA).

1. Notebooks `01_EDA_calidad_datos` → `02_datos_externos` → `03_construccion_funcional` → `04_eda_funcional`. El 01 crea el cache local `eda_fase0/estaciones_crudas.pkl` desde los archivos fuente; el 03 escribe `fda_fase1/coef_satelite_python.csv.gz` para la comparación con R.
2. `Rscript 03_construccion_funcional.R` → `Rscript tests/test_construccion_manual.R` → `python tests/validar_salidas_cruzadas.py` y `Rscript tests/validar_salidas_cruzadas.R`.
3. `Rscript scripts/crear_pares_estacion_pixel.R` → `python scripts/verificar_crosswalk_python.py` → `Rscript 04_eda_funcional.R`.
4. `python 05_red_funcional.py` (≈5 min en CPU).
5. Pruebas: `python -m unittest tests/test_red_funcional.py tests/test_rubrica_code.py tests/test_reconciliar_claves_estacion.py`, `Rscript tests/test_station_pixel_crosswalk.R` y `Rscript tests/test_station_pixel_statistics.R`.

La fase 1 no usa librerías funcionales para calcular Φ, MCO ni K; las comprobaciones con `splines` y QR están en `tests/test_construccion_manual.R`. En el EDA las librerías (scikit-fda, fda) se usan para verificar media, profundidades y boxplot. En la fase 4 se usa PyTorch, autorizada por el enunciado; la KAN se programa con la forma φ(x) = w_b·SiLU(x) + w_s·Σ c_i B_i(x).

La comparación de coeficientes satelitales R↔Python se hace con la misma entrada (diferencia < 1e-8). El notebook ajusta con los valores float32 del TIFF y R con `curvas_satelite_mm.csv.gz` (5 decimales); esa diferencia se acota exactamente con la norma de y2cMap.

### Flujo R anterior

`proyecto_funcional.R` y `resultados_R/` se conservan del flujo integrado anterior, que parte de otra población (4.426 curvas/238 estaciones, K=6/11). El flujo activo es el de los archivos 03/04/05.

## Datos crudos (no incluidos)

Los archivos fuente del curso se mantienen locales e ignorados por Git y deben copiarse en la raíz: `P_1930_1979.xlsx`, `P_1980_1990.csv`, `P_1991_2003.csv`, `ImagenesSatelitePentadal.zip` (o la carpeta `ImagenesSatelitePentadal/`), y los metadatos `CVC_Wheather_Stations_Long_Lat.csv/.dta` y `CatalogoEstacionesElevacion.xlsx`.

## Estado

Seguimiento detallado por criterio en `docs/SEGUIMIENTO_RUBRICA_JACOBO.md`.

- **Fases 0–2.** Reproducidas desde los archivos fuente: el cache y las curvas regeneradas son idénticos a los versionados. K=8 (estaciones) y K=11 (satélite) elegidos solo con entrenamiento. R y Python coinciden en claves, K, coeficientes de ambas fuentes, W, profundidades, fechas y atípicas, y el EDA queda contrastado con scikit-fda y fda.
- **Población.** Respecto de main (4.426 curvas/238 estaciones) cambia a 4.485/248. Las 191 claves distintas se explican por regla en `docs/diff_claves_estacion_main_jacobo.csv`: el mínimo de años (P6) ahora se cuenta solo en entrenamiento (≥7 años en 1981–1995, antes ≥10 en 1981–2003), y 10 claves salen por cobertura (P5) con las reglas relativas a cada estación.
- **Fase 3.** 2.751 pares estación–año/píxel–año en 147 estaciones con coordenada CVC dentro de la malla (2.051 con altitud), verificados con la geotransformada del TIFF. Es cobertura parcial: 78 estaciones no tienen coordenada en el archivo CVC y 23 caen fuera de la malla.
- **Fase 4.** Red funcional entrenada con los pares que tienen altitud (1.349 de entrenamiento, 445 de validación y 257 de prueba). En prueba (2001–2003), las tres cabezas mejoran a la curva satelital sin corregir en error L², KGE, sesgo anual y fecha de los picos (IC bootstrap del 95 % por estación). La KAN no supera a la cabeza lineal de forma significativa. La equivalencia de la cabeza lineal con la regresión funcional se verifica a 4e-15.
- **Pendiente.** Estaciones reservadas (las debe entregar el docente), coordenadas faltantes, informe final y consolidación en un solo script R y uno Python.

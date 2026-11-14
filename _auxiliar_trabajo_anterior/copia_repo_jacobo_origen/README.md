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

1. Notebooks `01` → `02` → `03` → `04` (Python).
2. Scripts `03_construccion_funcional.R` → `04_eda_funcional.R`, desde la raíz del repositorio. Usan solo R base (`splines` y `fda` únicamente para verificar). El de la fase 1 lee las curvas preparadas por el notebook 03 y construye sus propios coeficientes; el de la fase 2 trabaja con esos coeficientes de R y al final comprueba que K, coeficientes, matriz W, profundidades, fechas y atípicas coinciden con los de Python.

## Datos crudos (no incluidos)

Los archivos `P_1930_1979.xlsx`, `P_1980_1990.csv`, `P_1991_2003.csv` y `ImagenesSatelitePentadal` los entrega el curso; hay que copiarlos en la raíz del repo para ejecutar los notebooks desde cero.

## Estado

Fases 0–2 listas, en Python y en R. El satélite se trabaja solo en los 1 353 píxeles de tierra (los 90 de océano vienen rellenados). Para la fase 4 ya quedan la matriz de Gram W (`fda_fase1/matrices_gram.npz`) y los coeficientes en escala √ y log(1+P) (`fda_fase1/coeficientes_transformados.npz`). La fase 3 (emparejamiento estación–píxel) depende de la tabla de coordenadas y altitud de las estaciones, pendiente de confirmar con el profesor. La fase 4 (modelos) está por hacer.

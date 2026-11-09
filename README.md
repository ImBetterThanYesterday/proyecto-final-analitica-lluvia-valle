# Proyecto final de Analítica de Datos 2026-2S: lluvia del Valle del Cauca (1981–2013)

Modelado funcional de la precipitación diaria de estaciones y de estimaciones satelitales pentadales (cada año en cada sitio es un dato funcional), comparación de fuentes con métodos funcionales y reconstrucción de la curva de estación a partir de la curva satelital con KAN, MLP, cabeza lineal (≈ regresión funcional) y decodificador FPCA.

## Contenido

| Archivo / carpeta | Qué contiene |
|---|---|
| `01_EDA_calidad_datos.ipynb` | EDA y auditoría de calidad de los datos crudos de estaciones |
| `01b_valores_unicos.ipynb` | Valores únicos por variable, por estación y por año |
| `02_datos_externos.ipynb` | ONI, CHIRPS v2, NASA POWER, códigos CNE/IDEAM |
| `03_construccion_funcional.ipynb` / `.R` | Representación en B-splines (K por validación cruzada), Python y R |
| `04_eda_funcional.ipynb` / `.R` | EDA funcional: media, covarianza, profundidad, boxplot funcional, atípicas |
| `eda_fase0/`, `fda_fase1/`, `fda_fase2/` | Resultados (CSV, coeficientes, figuras) de cada fase |
| `datos_web/` | Datos externos descargados (ONI, CNE, CHIRPS v2, NASA POWER) |

## Datos crudos (no incluidos)

Los archivos `P_1930_1979.xlsx`, `P_1980_1990.csv`, `P_1991_2003.csv` y `ImagenesSatelitePentadal` los entrega el curso; hay que copiarlos en la raíz del repo para ejecutar los notebooks desde cero.

## Estado

Fases 0–2 listas. La fase 3 (emparejamiento estación–píxel) depende de la tabla de coordenadas y altitud de las estaciones, pendiente de confirmar con el profesor. La fase 4 (modelos) está por hacer.

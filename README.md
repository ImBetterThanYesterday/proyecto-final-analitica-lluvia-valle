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
| `KAN/entrenar_kan.py` | Arquitectura funcional KAN, búsqueda escalonada de hiperparámetros solo con validación y reconstrucciones de validación |
| `KAN/auditar_kan.py` | Auditoría post hoc del test: métricas, bootstrap por estación, Wilcoxon pareado, ENSO y gráficas esenciales |
| `KAN/resultados_validacion/` | Configuración seleccionada, modelos de tres semillas, comparación con MLP/cabeza lineal, métricas, curvas estación-año y `tabla_resultados_validacion.png` |
| `KAN/resultados_test/` | Métricas, curvas completas, ranking por L2 y diagnóstico de asociaciones por estación |
| `KAN/graficas_esenciales_test/` | Comparación de métricas, ejemplos representativos, todas las curvas, `04_tabla_resultados_test.png` y `05_diagnostico_estaciones_test.png` |
| `KAN_historico_2026-11-12/` | Respaldo íntegro de los scripts y resultados anteriores a la reorganización |
| `datos_web/` | Datos externos descargados (ONI, CNE, CHIRPS v2, NASA POWER) |

## Orden de ejecución

1. Notebooks `01` → `02` → `03` → `04` (Python).
2. Scripts `03_construccion_funcional.R` → `04_eda_funcional.R`, desde la raíz del repositorio. Usan solo R base (`splines` y `fda` únicamente para verificar). El de la fase 1 lee las curvas preparadas por el notebook 03 y construye sus propios coeficientes; el de la fase 2 trabaja con esos coeficientes de R y al final comprueba que K, coeficientes, matriz W, profundidades, fechas y atípicas coinciden con los de Python.
3. Desde la raíz, `python KAN\entrenar_kan.py` compara en validación neuronas funcionales, componentes FPCA, bases y grado spline, rugosidad, tasa de aprendizaje, tamaño de lote y weight decay. La rugosidad es la segunda derivada media de las funciones de peso sobre el año normalizado a [0,1]. Usa tres semillas por configuración, selecciona una configuración parsimoniosa, guarda las curvas observada/satelital/KAN y crea una tabla gráfica de métricas y comparaciones inferenciales en `KAN\resultados_validacion\`.
4. Después de congelar la configuración, `python KAN\auditar_kan.py` evalúa una sola vez el test 2001–2003 y genera las métricas, bootstrap por estación (2.000 réplicas), Wilcoxon pareado, desglose ENSO, diagnóstico post hoc de factores asociados al L2 por estación y gráficos —incluidas las tablas resumen— en las carpetas de test. No reentrena ni ajusta con test.

## Datos crudos (no incluidos)

Los archivos `P_1930_1979.xlsx`, `P_1980_1990.csv`, `P_1991_2003.csv`, los catálogos de estaciones y `ImagenesSatelitePentadal` están centralizados en `Datos/`. Los notebooks leen los insumos desde esa carpeta.

## Estado

Fases 0–2 listas, en Python y en R. El satélite se trabaja solo en los 1 353 píxeles de tierra (los 90 de océano vienen rellenados). La KAN reconstruye la curva de estación con la curva satelital en escala raíz, elevación y ONI anual; FPCA y normalización se ajustan solo con entrenamiento. La selección por validación eligió 6 neuronas funcionales, M=4 (84,2 % de varianza FPCA), 8 bases cúbicas, λ=1e-4, tasa 0,002, lote 64 y weight decay 0; la KAN tiene 354 parámetros frente a 356 del MLP. La penalización de rugosidad se calcula como curvatura en el dominio anual normalizado [0,1]. En validación, el ensamble KAN obtuvo L2 relativo 0,394 y KGE 0,465. El ensamble post hoc de tres semillas obtuvo en test L2 0,435 por curva (media por estación 0,440; IC95 % 0,403–0,482), KGE 0,476, sesgo anual +14,7 % y error de picos 21,9 días; el satélite sin corregir obtuvo L2 0,538, KGE 0,360, sesgo +30,2 % y error de picos 26,0 días. El diagnóstico post hoc por estación encontró mayor L2 asociado con mayor L2 satelital (Spearman ρ=0,407; q=0,00013), menor acumulado anual observado (ρ=-0,429; q=0,000060), mayor variabilidad diaria (ρ=0,308; q=0,0044) y mayor concentración de lluvia en el 10 % de días (ρ=0,363; q=0,00070). El cuartil de L2 más alto tuvo RMSE diario medio 2,40 frente a 1,72 mm/día y sesgo anual absoluto 60 % frente a 18 %; altitud y distancia al píxel no mostraron asociación clara. Como el L2 es relativo y el sesgo comparte errores con esa métrica, estos son patrones correlacionales, no causas demostradas. El test cubre Neutral y El Niño, pero no La Niña; sus 93 estaciones ya están en entrenamiento, así que no es prueba ciega ni mide generalización espacial. La confirmación docente del emparejamiento estación–píxel y la lista de estaciones reservadas siguen pendientes.

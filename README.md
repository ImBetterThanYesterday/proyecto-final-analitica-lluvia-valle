# Lluvia del Valle del Cauca 1981–2003 · Proyecto final de Analítica de Datos (UAO, 2026-2S)

La lluvia del Valle del Cauca se modela como datos funcionales: cada año de cada sitio es una curva. Hay dos fuentes:
- estaciones con lluvia diaria;
- satélite con lluvia por péntadas.

El trabajo tiene tres partes: construir las curvas con B-splines, hacer el EDA funcional y reconstruir la curva de la estación a partir de la del satélite.

## Estructura

| Carpeta / archivo | Qué es |
|---|---|
| `lluvia_valle_datos.R` | **El código de datos, en un solo archivo.** Va de los crudos al dataset limpio: población, EDA, experimento de huecos, depuración y dataset final. |
| `salidas_R_datos/` | Lo que produce ese archivo. Tiene `lluvia_valle_limpia.rds` (el dataset final), además de `tablas/` y `graficos/`. |
| `DECISIONES.md` | Bitácora de decisiones (D-01 a D-36), con su justificación, y preguntas para el profesor. |
| `CONTEXTO_EDA.md` | Resumen para retomar el trabajo, también con Claude Code. |
| `datos_web/` | Fuentes externas: estaciones del Valle, límites del DANE, CNE, ONI, CHIRPS y NASA POWER. |
| `datos_curso/` | **Datos crudos del curso. No se suben a git:** hay que copiarlos aquí (ver abajo). |
| `enunciado/` | Enunciado y fe de erratas. |
| `material_clase/` | PDF y scripts de clase, con sus datos de ejemplo. Son la referencia de los métodos permitidos. |
| `archivo_flujo_anterior/` | Flujo anterior (notebooks, fases 0–4, `proyecto_funcional.R`, pruebas, informes). Usaba 248 estaciones sin filtrar por el Valle. Se conserva como referencia y para reutilizar código: ya no es el flujo activo. |

## Cómo se corre

1. Copiar en `datos_curso/`:
   - `P_1930_1979.xlsx`, `P_1980_1990.csv` y `P_1991_2003.csv`;
   - `ImagenesSatelitePentadal.zip` o la carpeta `ImagenesSatelitePentadal/`;
   - `CVC_Wheather_Stations_Long_Lat.csv`, `CVC_Wheather_Stations_Long_Lat.dta` y `CatalogoEstacionesElevacion.xlsx`.
2. Desde la raíz del repositorio:

   ```bash
   Rscript lluvia_valle_datos.R
   ```

   Usa R 4.5 con data.table, ggplot2, readxl y sf. Tarda unos minutos.

## Periodo y partición (fe de erratas)

| Conjunto | Años |
|---|---|
| Entrenamiento | 1981–1995 |
| Validación | 1996–2000 |
| Prueba | 2001–2003, más las estaciones que reserve el profesor |

Las reglas que dependen de los datos se calculan solo con entrenamiento.

## Por dónde vamos (12-nov-2026)

| Paso | Estado |
|---|---|
| 1. Ubicar los 344 códigos y quedarse con el Valle | ✅ 232 estaciones, todas dentro del polígono del DANE (D-31 a D-35) |
| 2. EDA de los datos crudos | ✅ secciones 1–6 de `lluvia_valle_datos.R` |
| 3. Reglas de depuración sacadas de los datos | ✅ experimento de huecos y reglas D-36 |
| 4. Dataset limpio de estaciones | ✅ 176 estaciones y 3 326 estación-años válidos, en `lluvia_valle_limpia.rds` |
| 5. Satélite: recortar al Valle, pasar a mm/día y revisar los años 1985, 1988 y 2000 | ⏳ siguiente |
| 6. Emparejar cada estación con su píxel | ⏳ |
| 7. Construcción funcional (B-splines propias; K por ECM/VC/VCG) con las 176 | ⏳ se reutiliza código del flujo anterior |
| 8. EDA funcional de cada fuente (media, varianza, profundidad, boxplot, ENSO) | ⏳ |
| 9. Regresión funcional y red (KAN/MLP) | ⏳ se reutiliza `archivo_flujo_anterior/05_red_funcional.py` |
| 10. Entrega: un archivo R, un archivo Python, el informe y el anexo de IA | ⏳ |

**Abierto con el profesor** (ver `DECISIONES.md`, §8):
- P-1: los 21 códigos sin ubicar;
- P-3: si se puede usar CHIRPS;
- P-4: qué hacer con el año 2000 del satélite, que cae en validación.

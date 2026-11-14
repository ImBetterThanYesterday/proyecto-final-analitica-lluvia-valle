# Lluvia del Valle del Cauca 1981–2003 · Proyecto final de Analítica de Datos (UAO, 2026-2S)

La lluvia del Valle del Cauca se modela como datos funcionales: cada año de cada sitio es una curva. Hay dos fuentes:
- estaciones con lluvia diaria;
- satélite con lluvia por péntadas.

El trabajo tiene tres partes: construir las curvas con B-splines, hacer el EDA funcional y reconstruir la curva de la estación a partir de la del satélite.

## Estructura

| Carpeta / archivo | Qué es |
|---|---|
| `lluvia_valle_datos.R` | **El código de datos, en un solo archivo.** Va de los crudos a los datasets limpios: población, EDA, experimento de huecos, depuración, y el satélite (fuente, recorte, EDA y comparación con las estaciones), los datos preparados y 27 pruebas automáticas. |
| `salidas_R_datos/` | Lo que produce ese archivo. Tiene `datos_preparados.rds` (las matrices listas para la construcción funcional), `lluvia_valle_limpia.rds` (estaciones) y `satelite_valle.rds` (satélite y pares estación–píxel), además de `tablas/` y `graficos/`. |
| `lluvia_valle_funcional.R` | **Toda la parte funcional, en un solo archivo:** emparejamiento, construcción de las curvas, EDA funcional, comparación de fuentes y regresión funcional. Solo lee `salidas_R_datos/datos_preparados.rds`. Los cambios de lo funcional se hacen aquí (D-48). |
| `salidas_R_funcional/` | Lo que produce ese archivo: `tablas/` y `graficos/` (prefijos `F1_` a `F6_`), `para_python/` (insumos de la red) y `resultados_funcional.rds`. |
| `lluvia_valle_red.py` | **La red neuronal funcional y la ablación de la cabeza (lineal, KAN y MLP), en Python.** Solo lee `salidas_R_funcional/para_python/`. Solo entrena y exporta las predicciones; R las evalúa en la sección 7 del funcional (D-54). Orden: `Rscript lluvia_valle_funcional.R` → `python lluvia_valle_red.py` (~15 min; `torch`, `numpy`, `pandas`) → `Rscript lluvia_valle_funcional.R`. |
| `salidas_python/` | Lo que produce: `predicciones_coef_raiz.csv` (lo lee R), `pesos_codificador.csv`, `tablas/` (prefijos `P1_` a `P6_`, construcción y entrenamiento) y `versiones.txt`. |
| `BITACORA_CAMBIOS.md` | **Qué cambió en cada sesión y por qué**, para ponerse al día. |
| `CONTROL_RUBRICA.md` | **Lista de lo que pide el enunciado y cómo vamos**, ítem por ítem. |
| `DECISIONES.md` | Bitácora de decisiones (D-01 a D-54), con su justificación, y preguntas para el profesor. |
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

   Usa R 4.5 con data.table, ggplot2, readxl, sf y terra. Tarda unos minutos.
3. La parte funcional (cerca de 1 minuto; usa además fda y splines, solo para verificar):

   ```bash
   Rscript lluvia_valle_funcional.R
   ```

## Periodo y partición (fe de erratas)

| Conjunto | Años |
|---|---|
| Entrenamiento | 1981–1995 |
| Validación | 1996–2000 |
| Prueba | 2001–2003 (el profesor ya no reserva estaciones) |

Las reglas que dependen de los datos se calculan solo con entrenamiento.

## Por dónde vamos (14-nov-2026)

El detalle criterio por criterio está en `CONTROL_RUBRICA.md`.

| Paso | Estado |
|---|---|
| 1. Ubicar los 344 códigos y quedarse con el Valle | ✅ 232 estaciones, todas dentro del polígono del DANE (D-31 a D-35) |
| 2. EDA de los datos crudos | ✅ secciones 1–6 de `lluvia_valle_datos.R` |
| 3. Reglas de depuración sacadas de los datos | ✅ R0 a R7: experimento de huecos, reglas D-36, R4 calculada (D-49) y R0 con el récord OMM (D-50) |
| 4. Dataset limpio de estaciones | ✅ 176 estaciones y 3 299 estación-años válidos, en `lluvia_valle_limpia.rds` |
| 5. Satélite: fuente, recorte al Valle, mm/día y EDA | ✅ es CHIRP solo-infrarrojo, confirmado (D-43); acierta el clima medio pero no los años Niño/Niña (D-44). Secciones 9–12 de `lluvia_valle_datos.R` |
| 6. Emparejar cada estación con su píxel | ✅ 176 estaciones en 148 píxeles (D-38, D-44) |
| 6b. Datos preparados para lo funcional, con pruebas | ✅ `datos_preparados.rds`: 3 299 curvas de estación, 15 663 de satélite y sus pares con el ONI (D-45); Univalle con 996 m (D-46) |
| 7. Construcción funcional (B-splines propias; K por ECM/VC/VCG) con las 176 | ✅ sección 2 de `lluvia_valle_funcional.R`: K = 8 / 11 (D-39, D-48); K = 8 confirmado frente a la fecha de los picos (D-51) |
| 8. EDA funcional de cada fuente (media, varianza, profundidad, boxplot, ENSO) | ✅ sección 3 |
| 8b. Comparación funcional estación vs. satélite (d(t), correlación cruzada, covarianzas, altitud, ENSO) | ✅ sección 4, con todos los años |
| 9. Regresión funcional función a función | ✅ sección 5: FPCA M = 6, λ por validación; en prueba baja el error L² relativo de 0,51 a 0,44 frente al satélite (D-48) |
| 9b. Red funcional (cabezas lineal, KAN y MLP) y ablación | ✅ `lluvia_valle_red.py` (D-52): la lineal reproduce a R; en prueba L² relativo 0,440 / 0,434 / 0,420 (lineal / KAN / MLP) frente a 0,513 del satélite |
| 9c. Interpretación (w_k, activaciones KAN, ajuste simbólico) y ENSO por cabeza | ⏳ |
| Años 1985, 1988 y 2000 del satélite | ✅ entran al modelo, por indicación del profesor y verificado con datos (D-47, D-48) |
| 10. Entrega: un archivo R, un archivo Python, el informe y el anexo de IA | ⏳ |

**Abierto con el profesor** (ver `DECISIONES.md`, §8):
- ✅ P-1: los 21 códigos sin ubicar se descartan;
- ✅ P-3: no se usa CHIRPS ni otra fuente; solo el satélite entregado;
- ✅ P-4: 1985, 1988 y 2000 del satélite entran al modelo; quitarlos solo se usa como sensibilidad (D-47, D-48);
- ✅ P-8: el satélite es CHIRP (confirmado).

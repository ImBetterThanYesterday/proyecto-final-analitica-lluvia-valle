# Seguimiento vivo de la rúbrica — integración Jacobo en rama María

**Norma aplicada:** resumen autocontenido de la Fe de Erratas nueva: 1981–2003; entrenamiento 1981–1995, validación 1996–2000, prueba 2001–2003 más estaciones reservadas; construcción funcional manual, sin librerías funcionales. La rúbrica completa conserva los pesos 1=10 %, 2=15 %, 3=15 %, 4.1–4.5=45 %, 5=10 %, 6=5 %.

**Rama local actual:** `maria-eda-funcional-construccion-funcional`. El commit merge `1e6f00d` integra la rama local `mejoras-rubrica-jacobo` (commit `5efceb0`); se usó la referencia Jacobo ya cacheada, sin `fetch`. El merge es local y todavía no está publicado; no se agregaron archivos crudos.

**Estados:** `PENDIENTE`, `EN PROGRESO`, `VERIFICADO`, `BLOQUEADO`. Un estado `VERIFICADO` exige una prueba/comando ejecutado y evidencia guardada; la existencia de código o salida vieja no basta.

## Matriz criterio por criterio

| Criterio / peso | Meta para aspirar a 5 | Hecho comprobado al iniciar | Falta / no hecho | Librerías y cumplimiento | Evidencia de cierre | Estado |
|---|---|---|---|---|---|---|
| 1 Construcción / 10 % | Φ/MCO manual, ECM/VC/VCG y K por fuente, S/H idempotente y traza K, bandas puntuales, dominio común. | R manual se ejecutó: K=8/11; 4.485/31.119 curvas; partición K solo entrenamiento; H, VC-LOO, W y base pasaron. Se quitó rm(list=ls()) y los verificadores splines/lsfit salieron del constructor. | El notebook 03 no se reejecutó; el CSV Python satelital nuevo no existe localmente. Cross-check R vs refs guardadas: estación 6,79e-13; satélite 5,25e-6 contra NPZ Python float32; W 8,35e-14/4,62e-14. Repetir tras reejecutar notebook 03. | Producción R no usa librerías funcionales; `splines`/QR se usan solo en `tests/test_construccion_manual.R`. | `Rscript tests/test_construccion_manual.R`; `python3.13 tests/validar_salidas_cruzadas.py`; reproducir notebook 03 con dependencias locales permitidas. | EN PROGRESO |
| 2 EDA / 15 % | Media verificada, varianza/DE/CV, covarianza/correlación, derivadas e interpretación de temporadas, FM/MBD, boxplot y claves, contexto. | R EDA ejecutado; MBD vs fuerza bruta error 1,7e-15. FM normalizado, FM_pdf separado, cuantiles con empates, fences inclusivos, fechas candidatas y ENSO descriptivo. | fda/scikit-fda no están instalados; Notebook Python no se pudo ejecutar y sus outputs guardados son antiguos. R registra media/boxplot como PENDIENTE. | EDA permite librerías; la verificación no se declara completa si no se ejecuta. | Reejecutar con fda/scikit-fda y verificar medias ambas fuentes, FM/MBD y fbplot con outputs nuevos. | EN PROGRESO |
| 3 Comparación / 15 % | Pares estación–píxel, regla/coordenadas, diferencias/bandas, covarianza/correlación, altitud y ENSO. | CVC CSV/DTA: 238 códigos y coordenadas concordantes; 170 IDs en Jacobo, 147 dentro de malla, 105 con altitud; 2.751 pares. Diferencia estación−satélite, crosscov/cor y bootstrap de estación 1000; 173/365 IC95 % puntuales excluyen cero. | 147/248 estaciones solamente; 78 sin CVC, 23 fuera de malla y 42 en malla sin altitud. No universalizar ni llamar sesgo causal. Falta sensibilidad/representatividad y revisión del criterio de nearest-pixel. | R `terra` valida el TIFF; fórmula explícita de centro/nearest, bootstrap por conglomerados. | Tablas de pares/cobertura/diferencias/IC/cov/cor y script reproducible; faltan cobertura integral y validación docente. | EN PROGRESO |
| 4.1 Depuración / 6 % | Reglas justificadas para imposibles, acumulados, ceros, huecos/cobertura y trazabilidad sin fuga. | Cache local regenerada: 2.992.694 filas, 344 códigos, 1930–2003. Diff exacto: 10 claves baseline-only por P5, 56 por P6; 125 claves Jacobo-only pertenecen a 15 IDs nuevos. | Notebook 01 no se pudo reejecutar completo por falta de rasterio; p99/flags P1–P6 requieren reproducción/sensibilidad. El pkl se generó localmente vía readxl+pandas, no se sube. | Sin nuevas descargas/instalaciones; fuente original intacta. | `docs/diff_claves_estacion_main_jacobo.csv`, cache local y reproducción completa de 01 cuando se autoricen/dispongan dependencias. | EN PROGRESO |
| 4.2 Red / 10 % | Gram, FPCA de entrenamiento, pérdida L²/rugosidad, equivalencia lineal ≤1e-4, no negatividad y recortes. | W y coeficientes transformados están preparados; no hay red completa. | Implementar encoder y cabezas lineal/KAN/MLP sin fuga. | PyTorch/pykan solo dentro del apartado autorizado de modelos. | Tests de Gram/FPCA/encoder/equivalencia y reporte de recortes. | PENDIENTE |
| 4.3 Ablación / 11 % | Modelos comparables, parámetros ±10 %, tres semillas, presupuesto equivalente, bootstrap e igualdad sintética KAN. | No hay ablación implementada. | Implementar protocolo y análisis de incertidumbre; no seleccionar con prueba. | Registrar librerías/versiones del apartado de modelos. | Matriz de parámetros, tres semillas, RMSE sintético y bootstrap reproducible. | PENDIENTE |
| 4.4 Evaluación / 8 % | Estaciones reservadas, baseline satelital, L²/KGE/sesgo/picos, bootstrap≥1000, Wilcoxon y CSV prueba. | La partición está escrita; K se elige con entrenamiento. | No hay red evaluada ni IDs de estaciones reservadas en el checkout. | Ninguna transformación/hyperparámetro aprende de validación/prueba. | CSV station_id/year/day/value_mm, métricas e intervalos reproducibles. | BLOQUEADO por estaciones reservadas/modelo |
| 4.5 Interpretación / 10 % | Pesos/ventanas KAN, activaciones, cabeza lineal, ajuste simbólico con R² y casos extremos. | Hay interpretación climatológica, no de un modelo. | Construir análisis de modelo luego de entrenar/evaluar. | Métodos del apartado de redes solamente. | Figuras/tablas y ajuste simbólico verificables. | PENDIENTE |
| 5 Informe / 10 % | PDF final ≤20 páginas, una columna, autores, fuentes, figuras interpretadas e informe IA/prompts. | Se generó un Word de auditoría/handoff para Claude; README y notebooks explican fases. Ninguno es el informe académico final. | Producir informe global y anexo IA con resultados validados. | Declarar herramienta y prompts reales; citar fuentes. | PDF/Word revisados por formato/páginas y anexos presentes. | PENDIENTE |
| 6 Scripts / 5 % | Un R y un Python finales, reproducibles, rutas/semillas/versiones, separados del informe. | Hay dos scripts R de fases y varios notebooks Python. | Consolidar los dos entregables y ejecutar desde ambiente limpio. | No instalar paquetes al vuelo; dependencias documentadas. | Registro de versiones/comandos y ejecución limpia con artefactos esperados. | PENDIENTE |

## Historial de cambios de estado

| Fecha | Criterio | Estado anterior → nuevo | Evidencia / comando | Nota |
|---|---|---|---|---|
| 2026-11-11 | Preparación de trabajo | — → EN PROGRESO | `git switch -c mejoras-rubrica-jacobo origin/actualizacion_jacobo`; `git status --short --branch` | Se usó la referencia local cacheada; no se consultó el servidor. |
| 2026-11-11 | Dependencias/datos | — → EN PROGRESO | Inventario local y regeneración del cache vía readxl + pandas 3.13. | Cache local: 2.992.694 filas, 344 IDs, 1930–2003; notebook completo sigue bloqueado por rasterio. |
| 2026-11-11 | Reconciliación de población | EN PROGRESO → EN PROGRESO | `python3 -m unittest tests/test_reconciliar_claves_estacion.py -v`; `python3 scripts/reconciliar_claves_estacion.py` | 4.360 claves comunes; 10 main-only P5 y 56 P6; 125 Jacobo-only corresponden a 15 IDs nuevos; CSV de diferencias generado. |
| 2026-11-11 | Crosswalk/comparación | EN PROGRESO → EN PROGRESO | `Rscript scripts/crear_pares_estacion_pixel.R`; `Rscript 04_eda_funcional.R` | 2.751 pares en 147 estaciones; 2.051 con altitud; cluster bootstrap puntual 1.000; 173/365 días IC95 % excluye cero. Parcial. |
| 2026-11-11 | Construcción manual | EN PROGRESO → EN PROGRESO | `Rscript tests/test_construccion_manual.R`; `python3.13 tests/validar_salidas_cruzadas.py` | K=8/11; base, QR/MCO, H/S, VCG/LOO y W OK; cross-check guardado: station 6,79e-13; sat 5,25e-6 (NPZ float32). Notebook 03 no rerun. |
| 2026-11-11 | Pruebas de fuente | — → VERIFICADO | 9 tests Python, crosswalk/bootstrap R y `parse()` de scripts | Pruebas locales pasan; no prueban redes ni holdout. |
| 2026-11-11 | EDA R/pares | EN PROGRESO → EN PROGRESO | `Rscript 04_eda_funcional.R` | EDA R y comparación parcial se ejecutan; fda y Python crosschecks quedan PENDIENTE explícito. |
| 2026-11-11 | Ambiente Python/R | — → BLOQUEADO para rerun notebooks | Chequeo de módulos con `python3.13`; `requireNamespace()` en R | Python 3.13 tiene pandas/numpy/matplotlib, pero falta rasterio/SciPy/scikit-fda/openpyxl; R tiene readxl/haven/terra, pero no fda. No se instalaron dependencias externas. |
| 2026-11-11 | Handoff de implementación | — → VERIFICADO | Pandoc + ZIP check de `INFORME_CAMBIOS_Y_ESTADO_RUBRICA_JACOBO.docx`; revisión de `PROMPT_IMPLEMENTAR_RUBRICA_JACOBO.md` | Word resume cambios/pendientes/pruebas para Claude; es contexto interno y no el informe final del criterio 5. |
| 2026-11-11 | Integración en rama María | — → VERIFICADO localmente | Merge `1e6f00d` (`mejoras-rubrica-jacobo` → `maria-eda-funcional-construccion-funcional`) | Se conservaron archivos previos de María y se eligieron los scripts completos 03/04 de Jacobo; `proyecto_funcional.R` queda como flujo integrado anterior. Merge aún no subido. |
| 2026-11-11 | Handoff Claude/Jacobo | — → VERIFICADO | `AUDITORIA_RUBRICA_RAMA_JACOBO.docx`; `INFORME_CAMBIOS_Y_ESTADO_RUBRICA_JACOBO.docx`; `PROMPT_IMPLEMENTAR_RUBRICA_JACOBO.md` | Informe final de implementación resume preexistente vs nuevo, pruebas, bloqueos y próximo paso; no reemplaza el informe final de curso. |
| 2026-11-11 | Informe posterior al merge | — → VERIFICADO localmente | `INFORME_FINAL_RAMA_MARIA_JACOBO.docx` validado con Pandoc y ZIP integrity | Resume merge 1e6f00d, estado de rúbrica y procedimiento de publicación; el remoto no se modificó. |

## Registro de ejecuciones

Añade cada prueba al cerrar un lote:

| Fecha | Lote | Comando | Resultado real | Artefacto/log |
|---|---|---|---|---|
| — | — | — | — | — |

## Bloqueos externos/locales

- Ampliar/justificar cobertura: 147/248 estaciones con coordenadas CVC dentro de la malla; 23 fuera, 78 sin CVC y 42 en malla sin altitud. No generalizar la comparación parcial.
- Ejecutar el notebook 04 actualizado para generar FM_pdf y salidas depth/boxplot nuevas; el Python local carece de rasterio/SciPy/scikit-fda.
- Reproducir los notebooks 01/03 desde crudos; el cache está regenerado, pero rasterio falta para el flujo completo del notebook.
- Conseguir del docente las coordenadas/altitudes e IDs de estaciones reservadas si el crosswalk local no cubre esos datos.
- Implementar 4.2–4.5 y no marcar criterios 3, 4.4 o 5 completos hasta cerrar la cobertura, holdout/modelos e informe.

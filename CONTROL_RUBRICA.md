# Control de la rúbrica

Lista de todo lo que pide el enunciado (`enunciado/Proyecto_Final_Analitica_FeErratas.docx`) y cómo vamos, criterio por criterio de la rúbrica. Se actualiza cada vez que se avanza.

**Estados:** ✅ hecho y verificado en el código · 🟡 hecho, falta redactarlo o un detalle · ⬜ pendiente · ♻️ existe en `archivo_flujo_anterior/05_red_funcional.py` con los insumos viejos (248 estaciones) y hay que reconectarlo a `salidas_R_funcional/para_python/`.

**Última actualización:** 14-nov-2026, commit `b1bd344` (R0, P-20 cerrada, 27 + 12 pruebas pasan en la Mac del grupo).

## Resumen

| Ítem | Peso | Estado | Qué falta |
|---|---|---|---|
| 1. Construcción de datos funcionales | 10 % | ✅ | Redactarlo en el informe |
| 2. EDA funcional | 15 % | ✅ | Redactarlo; ventanas de temporada en el Pacífico (ver abajo) |
| 3. Comparación funcional por fuente | 15 % | ✅ | Redactar la interpretación del sesgo |
| 4.1 Preprocesamiento y selección | 6 % | 🟡 | Decidir P-19 (registros 1930–1979) |
| 4.2 Red neuronal funcional | 10 % | ♻️ | La red en Python; la regresión de referencia en R ya cumple |
| 4.3 Ablación KAN, MLP y lineal | 11 % | ♻️ | Todo en Python |
| 4.4 Evaluación | 8 % | 🟡 | Hecha para la regresión; falta para las tres cabezas. Confirmar el CSV de estaciones reservadas |
| 4.5 Interpretación climatológica | 10 % | ♻️ | Todo en Python, más la lectura climatológica |
| 5. Informe PDF | 10 % | ⬜ | Todo |
| 6. Scripts R y Python | 5 % | 🟡 | Unir los dos R en uno; el Python |

Lo que pesa más y falta es lo de Python: 4.2 + 4.3 + 4.4 + 4.5 suman 39 % de la nota.

## 1. Construcción de datos funcionales (10 %) — sección 2 de `lluvia_valle_funcional.R`

| Lo que pide la rúbrica | Estado | Dónde se verifica |
|---|---|---|
| Φ y coeficientes por mínimos cuadrados, paso a paso, sin librerías | ✅ | Funciones `base_bspline` (Cox–de Boor), `ajuste_mco` y `construir`. Coincide con `splines::splineDesign` a 3e-16 (prueba 3) |
| K con ECM, VC y VCG: gráfica y valor elegido por fuente | ✅ | `F2_seleccion_K.png`, `F2_K_elegido.csv`: K = 8 (estaciones) y 11 (satélite), solo con entrenamiento (D-39) |
| K = 8 frente a la fecha de los picos | ✅ | `F2_K_picos_vs_referencia.csv`, `F2_K_que_gana_por_referencia.csv`: el K que mejora los picos depende de la referencia, así que se queda el de la VCG (D-51) |
| S idempotente y tr(S) = K | ✅ | `F2_controles_S.csv`, error < 1e-14 (prueba 4) |
| Covarianza residual y bandas punto a punto, interpretadas como punto a punto | ✅ | `F2_covarianza_residual.png`, `F2_bandas_punto_a_punto.png` |
| Ambas fuentes en el mismo dominio anual | ✅ | `F2_mismo_dominio.png`: día 1–365, péntadas en su día central |
| Tope: curvas con librerías → máximo 2 | ✅ no aplica | Las librerías solo verifican |

Para el informe: los días faltantes no se imputan; cada curva se ajusta con sus días observados y su propia S.

## 2. EDA funcional (15 %) — sección 3 de `lluvia_valle_funcional.R`

| Lo que pide la rúbrica | Estado | Dónde se verifica |
|---|---|---|
| Media paso a paso, verificada con la función de librería | ✅ | `F3_verificacion_media.csv`, contra `fda::mean.fd` (prueba 6) |
| Varianza y desviación estándar funcional | ✅ | `F3_media_y_dispersion.png`, `F3_dispersion_resumen.csv` |
| Covarianza y correlación de cada fuente | ✅ | `F3_covarianza_*.png`, `F3_correlacion_*.png`, también por vertiente y con anomalías |
| Derivadas e inicio y fin de las temporadas | 🟡 | `F3_derivada_media.png`, `F3_temporadas_por_curva.png`. Las ventanas suponen régimen bimodal y se aplican también al Pacífico; el fin de la 2.ª temporada no se puede fechar en el 62 % de las curvas. Hay que justificarlo o separar el Pacífico |
| Profundidad FM y MBD, boxplot funcional, atípicas por estación y año | ✅ | `F3_boxplot_funcional.png`, `F3_atipicas_*.csv`; MBD y atípicas coinciden con `fda::fbplot` (prueba 7). Opcional: verificar FM con `fda.usc::depth.FM` |
| Interpretación climatológica (bimodal, orografía, ENSO) | 🟡 | Escrita en el bloque "Qué dicen los datos"; falta pasarla al informe |
| Tope: subconjunto sin justificar → máximo 2 | ✅ no aplica | Se usan las 176 estaciones y los 23 años; representatividad en la sección 8 de `lluvia_valle_datos.R` |

Falta justificar por qué se modela la tasa en mm/día y no la lluvia acumulada que usa el ejemplo de clase `Script_AnalisisExploratorio.R`.

## 3. Comparación funcional por fuente (15 %) — secciones 1 y 4 de `lluvia_valle_funcional.R`

| Lo que pide la rúbrica | Estado | Dónde se verifica |
|---|---|---|
| Emparejamiento estación–píxel: coordenadas, regla y píxeles compartidos | ✅ | `F1_emparejamiento_estacion_pixel.csv`, `F1_pixeles_compartidos.csv`: 176 estaciones en 148 píxeles, píxel que contiene la estación (P-18) |
| d(t) con banda ± 2EE e intervalos donde no contiene el cero | ✅ | `F4_diferencia.png`, `F4_intervalos_diferencia.csv` |
| Correlación cruzada y su diagonal, interpretadas | ✅ | `F4_correlacion_cruzada.png`, `F4_diagonal_correlacion.png` |
| Comparación de las covarianzas | ✅ | `F4_covarianzas.png`, `F4_covarianzas_comparadas.csv` |
| Por altitud o subregión y por fase ENSO | ✅ | `F4_comparacion_por_grupo.csv`, `F4_diferencia_por_grupo.png`, `F4_sesgo_altitud.png` |
| Interpretación climatológica del sesgo | 🟡 | Base en D-44, sección 12 de datos y bloque de la sección 4; falta redactarla |
| Tope: sin emparejamiento documentado → máximo 1 | ✅ no aplica | |

El diagnóstico global no depende de usar K distinto por fuente; el veredicto día a día coincide en 281 de 365 días (`F4_sensibilidad_mismo_K.csv`). Opcional: ANOVA funcional o prueba t por permutación (scripts de clase) para ENSO y altitud.

## 4.1 Preprocesamiento y selección (6 %) — `lluvia_valle_datos.R`

| Lo que pide la rúbrica | Estado | Dónde se verifica |
|---|---|---|
| Regla explícita para valores imposibles y errores de digitación | ✅ | R0: negativos y > 1 825 mm (récord OMM), `2_R0_valores_imposibles.csv` (D-50). R2: extremos sin respaldo de las vecinas |
| Totales acumulados en un solo día y ceros que son faltantes | ✅ | R3 (totales mensuales) y R4 (rachas de ceros con la vecina y umbral de Bonferroni, D-49) |
| Cobertura mínima por estación y año, justificada | ✅ | R5 (≥ 90 % y hueco ≤ 15 días) y R6 (≥ 7 años de entrenamiento), con el experimento de huecos (sección 7) |
| Decisión justificada sobre agregar o eliminar datos | 🟡 | Sin imputación (D-13). **Falta P-19:** usar o no 1930–1979 para el FPCA. Recomendación: no usarlos (la fe de erratas fija 1981–2003 y el FPCA con 2 257 curvas ya explica 94 % con M = 6) |
| Trazabilidad: estaciones y pares antes y después de cada paso | ✅ | `8_trazabilidad.csv`: 323 → 220 → 176 estaciones, 3 299 pares |
| Nada usa datos de prueba | ✅ | Umbrales y experimentos solo con 1981–1995; prueba automática |

Pendiente menor: sensibilidad con y sin las 6 estaciones del cruce con el CNE y las 12 coordenadas a revisar.

## 4.2 Red neuronal funcional (10 %) — Python

| Lo que pide la rúbrica | Estado | Base |
|---|---|---|
| Codificador con w_k(t) en B-splines y z = A W c, Gram por cuadratura | ♻️ | `05_red_funcional.py` (8 neuronas). W ya viene por Gauss–Legendre en `para_python/W_s.csv` |
| Altitud y ONI agregados a z | ♻️ | En `para_python/pares.csv`, con la media y la DE de entrenamiento en `especificacion.json` |
| FPCA solo con entrenamiento; M justificado y varianza reportada | ✅ en R / ♻️ en Python | R: M = 6 (94,2 %), elegido en validación (`F5_fpca_varianza.csv`). Python debe usar el mismo M y la misma ψ (`fpca_estaciones.csv`) |
| Pérdida L² con W y penalización de rugosidad sobre w_k | ♻️ | `05_red_funcional.py`; R_s en `para_python/R_s.csv` |
| Demostración escrita y verificación ≤ 1e-4 de la equivalencia con la regresión funcional | ✅ en R / ♻️ en Python | R: 4e-15. Python: `verificar_equivalencia_lineal` |
| Reconstrucción no negativa con el % recortado | ✅ en R / ♻️ en Python | R reporta el % recortado (sección 5.5) |
| Tope: MLP sobre el vector de péntadas → máximo 1 | | No aplica si se usa el codificador funcional |

## 4.3 Ablación KAN, MLP y cabeza lineal (11 %) — Python

| Lo que pide la rúbrica | Estado | Base |
|---|---|---|
| Las tres cabezas con el mismo codificador y el mismo decodificador | ♻️ | Clase `RedFuncional` |
| Parámetros de KAN y MLP con `sum(p.numel())`, diferencia ≤ 10 % | ♻️ | `n_param` y `parametros_cabezas.csv`, con `assert` |
| Tres semillas, media y DE | ♻️ | `SEMILLAS = (2026, 2027, 2028)` |
| Mismo presupuesto: épocas, parada temprana con validación, pruebas de hiperparámetros | ♻️ | 6 000 épocas máx., paciencia 500; λ se elige con la cabeza lineal. Hay que tabular el presupuesto |
| KAN verificada con una función sintética (RMSE ≤ 1e-2) | ♻️ | `verificar_kan_sintetica` |
| Una cabeza supera a otra solo si el IC 95 % bootstrap de la diferencia excluye el 0 | ♻️ | Comparaciones pareadas por estación |
| Tope: hiperparámetros con prueba → máximo 1 | | Solo validación |

## 4.4 Evaluación (8 %)

| Lo que pide la rúbrica | Estado | Dónde |
|---|---|---|
| Partición 1981–1995 / 1996–2000 / 2001–2003 sin fuga | ✅ | Datos preparados; prueba de no fuga del funcional (prueba 9) |
| Estaciones reservadas | 🟡 | El profesor dijo que ya no van (P-10), pero la fe de erratas aún pide el CSV `station_id, year, day, value_mm`. **Confirmar con él** |
| Comparación con el satélite sin corregir | ✅ regresión / ♻️ red | `F5_metricas_prueba.csv` |
| Las cuatro métricas | ✅ regresión / ♻️ red | L² relativo 0,51 → 0,44; KGE 0,38 → 0,47; sesgo anual +26 % → +7 %; picos |
| IC bootstrap 95 % (≥ 1 000) y Wilcoxon pareado por estación | ✅ regresión / ♻️ red | 2 000 remuestreos, `F5_regresion_vs_satelite_pareado.csv` |
| Por fase ENSO con el ONI y su fuente | ✅ regresión / ♻️ red | `F5_metricas_por_fase_ENSO.csv`; ONI de NOAA CPC. P-16 (ONI anual = media de los 12 trimestres) sigue abierta |

## 4.5 Interpretación climatológica (10 %) — Python

| Lo que pide la rúbrica | Estado | Base |
|---|---|---|
| w_k(t) graficadas con la ventana del año de cada neurona | ♻️ | Figuras de `05_red_funcional.py` |
| Activaciones de la primera capa de la KAN en el rango de los datos | ♻️ | |
| Contraste con los coeficientes de la cabeza lineal | ♻️ | En R ya está β(s,t) (`F5_beta_superficie.png`) |
| Ajuste simbólico de al menos una activación con su R² | ♻️ | `ajuste_simbolico` |
| Lectura climatológica coherente con las gráficas | ⬜ | |
| Peores casos con sus curvas | ⬜ | |

## 5. Informe en PDF (10 %)

| Lo que pide la rúbrica | Estado |
|---|---|
| PDF a una columna, ≤ 20 páginas con referencias, todos los integrantes | ⬜ |
| Tablas y figuras numeradas, tituladas, citadas e interpretadas | ⬜ |
| Referencias de métodos y fuentes (estaciones, CHIRP, ONI, OMM) | 🟡 Fuentes listadas en `DECISIONES.md`; falta la cita formal de CHIRP |
| Anexo en Word con la herramienta de IA y los prompts | 🟡 **Hay que ir guardando los prompts** (se usó Claude) |
| Redacción propia | ⬜ |
| Topes: sin R Markdown, dos columnas ni solo Colab; IA no declarada → máximo 2 | |

## 6. Scripts en R y Python (5 %)

| Lo que pide la rúbrica | Estado |
|---|---|
| Un único R y un único Python | 🟡 Hoy hay dos R (datos y funcional); se pegan uno tras otro para la entrega (D-48). Falta el Python |
| Corren de principio a fin, rutas relativas, semillas fijas, versiones listadas | ✅ R: 27 pruebas en datos y 12 en funcional; `sessionInfo.txt`; fijan UTF-8 solos y leen el satélite por lotes (límite de archivos de macOS) |
| Solo procedimientos de clase; librerías extra solo donde se permite | 🟡 Construcción sin librerías; en Python solo PyTorch (pykan permitido) |
| Cada tabla y figura del informe ubicable en el código | 🟡 Prefijos `2_`…`15_` (datos) y `F1_`…`F6_` (funcional); falta en Python |

## Plan para cerrar (en orden)

1. **Red funcional en Python:** reconectar `05_red_funcional.py` a `para_python/`, con M = 6, ψ y escalas de R; un solo archivo `lluvia_valle_red.py`. Cubre 4.2, 4.3, 4.4 y 4.5.
2. **Decisiones abiertas:** P-19 (recomendado no usar 1930–1979), P-16 (ONI anual) y confirmar el CSV de estaciones reservadas.
3. **Texto respaldado:** tasa mm/día vs. acumulada; ventanas de temporada en el Pacífico; interpretación del sesgo.
4. **Unir los dos R** en un solo archivo de entrega.
5. **Informe PDF** y **anexo de IA**.

## Reglas generales que hay que cuidar

- Entrega: PDF llamado "Proyecto final AnalíticaDeDatos", más un único R y un único Python (y el repositorio en zip o rar).
- Sustentación de 40 minutos; el profesor puede preguntarle cualquier fase a cualquier integrante, y una respuesta mal dada resta 0,2.
- Solo scripts, procedimientos y metodologías vistos en clase (`material_clase/`).

## Pendientes con el profesor

Detalle en `DECISIONES.md`, §8.

- ✅ P-1, P-3, P-4, P-8, P-10, P-12 a P-15, P-17 y P-18 respondidas en la asesoría del 13-nov (D-47).
- ✅ P-20 (K de las estaciones) cerrada por el equipo con datos (D-51).
- ⬜ P-16 (definición del ONI anual) y P-19 (registros 1930–1979).
- ⬜ Confirmar que no hace falta el CSV de estaciones reservadas que todavía menciona la fe de erratas en 4.4.

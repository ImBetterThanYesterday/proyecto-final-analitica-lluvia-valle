# Parte funcional: dónde está cada ítem de la rúbrica

El orden lo fija el enunciado: construir las curvas, hacer el EDA funcional y comparar las fuentes. Después viene la red funcional. La regresión funcional función a función es la referencia lineal; la cabeza lineal de la red debe reproducirla.

Todo el código está en **`lluvia_valle_funcional.R`** (un solo archivo, D-48) y sus salidas en `salidas_R_funcional/`: figuras en `graficos/` y tablas en `tablas/`, con los prefijos `F1_` a `F6_` según la sección. La entrada es `salidas_R_datos/datos_preparados.rds`. Entran todos los años (D-47).

## 1. Construcción de datos funcionales (10 %): sección 2

| Ítem | Dónde | Resultado |
|---|---|---|
| Φ y coeficientes por mínimos cuadrados, sin librerías | Funciones `base_bspline`, `ajuste_mco` y `construir`, al inicio del archivo | Coincide con `splines::splineDesign` (3e-16) |
| K por ECM, VC y VCG para cada fuente | `F2_seleccion_K.png`, `F2_K_elegido.csv` | K = 8 en estaciones y K = 11 en satélite, solo con entrenamiento |
| Sensibilidad de K | `F2_sensibilidad_K_estaciones.csv` | Con K = 12, 16 o 24 los picos se mueven a lo sumo 9 días |
| S idempotente y tr(S) = K | `F2_controles_S.csv` | Error < 1e-14 en el peor ajuste |
| Covarianza residual y bandas punto a punto | `F2_covarianza_residual.png`, `F2_bandas_punto_a_punto.png` | Bandas con σ²I y σ²(t); son punto a punto |
| Ambas fuentes en el mismo dominio | `F2_mismo_dominio.png` | Día 1–365; las péntadas se ubican en su día central |

## 2. EDA funcional (15 %): sección 3

| Ítem | Dónde |
|---|---|
| Media paso a paso, verificada con `fda::mean.fd` | `F3_verificacion_media.csv`, `F3_media_y_dispersion.png` |
| Varianza y DE | `F3_media_y_dispersion.png`, `F3_dispersion_resumen.csv` |
| Covarianza y correlación de cada fuente | `F3_covarianza_*.png`, `F3_correlacion_*.png` |
| Covarianza por vertiente y con anomalías | `F3_covarianza_por_vertiente_y_anomalias.csv`, `F3_correlacion_anomalias.png` |
| Derivadas, inicio y fin de las temporadas | `F3_derivada_media.png`, `F3_temporadas_por_curva.png`, `F3_temporadas_fechas.csv` |
| FM, MBD y boxplot funcional; atípicas por estación y año | `F3_boxplot_funcional.png`, `F3_atipicas_*.csv`; las atípicas coinciden con `fda::fbplot` |
| Interpretación: régimen bimodal, orografía y ENSO | `F3_contexto.png`, `F3_contexto_vertiente_altitud_ENSO.csv` y el bloque "Qué dicen los datos" |
| Representatividad | Se usan las 176 estaciones y los 23 años, sin submuestrear |

## 3. Comparación funcional por fuente (15 %): secciones 1 y 4

| Ítem | Dónde |
|---|---|
| Emparejamiento documentado | Sección 1: `F1_emparejamiento_estacion_pixel.csv` y `F1_pixeles_compartidos.csv` |
| Vertiente de los píxeles, verificada | `F1_regla_vertiente_pixeles.csv` |
| d(t) ± 2EE e intervalos sin el cero | `F4_diferencia.png`, `F4_intervalos_diferencia.csv`; el EE se calcula por estación |
| Correlación cruzada y su diagonal | `F4_correlacion_cruzada.png`, `F4_diagonal_correlacion.png` |
| Covarianzas comparadas | `F4_covarianzas.png`, `F4_covarianzas_comparadas.csv` |
| Mismo K en las dos fuentes | `F4_sensibilidad_mismo_K.csv` |
| Por altitud, subregión y ENSO | `F4_diferencia_por_grupo.png`, `F4_comparacion_por_grupo.csv`, `F4_sesgo_altitud.png` |
| El error del satélite sigue al ENSO | `F1_cociente_satelite_estacion_vs_ONI.png`, `F1_satelite_vs_estaciones_por_anio.csv` |
| Interpretación del sesgo | Bloque "Qué dicen los datos" de la sección 4 |

## 4. Regresión funcional y paso a la red: sección 5

| Qué | Resultado |
|---|---|
| Orden de la raíz y el suavizado | `F5_orden_raiz_y_suavizado.csv`: suavizar en mm conserva el acumulado anual |
| FPCA de las estaciones, solo con entrenamiento | M = 6 (94,2 % de la varianza), elegido junto con λ = 3,2e8 en validación; `F5_fpca_varianza.csv` |
| β(s,t) = Σ β_m(s) ψ_m(t) con penalización ∫β'' | `F5_beta_superficie.png`, `F5_efectos.png` |
| Equivalencia entre la integral de β y los puntajes | Diferencia de 4e-15 (el criterio es ≤ 1e-4) |
| Las cuatro métricas con IC bootstrap | `F5_metricas_prueba.csv`: L² relativo 0,51 → 0,44; KGE 0,38 → 0,47; sesgo anual +26 % → +7 % |
| Wilcoxon pareado por estación | `F5_regresion_vs_satelite_pareado.csv` |
| Por fase ENSO | `F5_metricas_por_fase_ENSO.csv` (la prueba no tiene año Niña; la validación va como orientativa) |
| Sensibilidad sin 1985, 1988 y 2000 | `F5_sensibilidad_anios_marcados.csv`: quitarlos empeora la prueba |
| Insumos para Python | `para_python/`: pares, coeficientes en mm y en raíz, W, R, FPCA, coeficientes de la regresión y `especificacion.json` |

Las 12 pruebas automáticas están en la sección 6 (`F6_pruebas.csv`).

## 5. Lo que sigue

1. **Red funcional en Python** (`.py` único). Hay que adaptar `archivo_flujo_anterior/05_red_funcional.py` para que lea `salidas_R_funcional/para_python/`:
   - codificador z = A W c con 6 a 8 neuronas, más altitud y ONI;
   - decodificador FPCA congelado (μ y ψ de `fpca_estaciones.csv`);
   - cabezas lineal, KAN y MLP con parámetros igualados (±10 %) y 3 semillas;
   - evaluación con las mismas métricas de la sección 5.5.

   La cabeza lineal debe reproducir los coeficientes de `coeficientes_regresion.csv` o, al menos, la identidad de β(s,t).
2. **Entrega del R:** pegar `lluvia_valle_datos.R` y `lluvia_valle_funcional.R` en un solo archivo (D-48).
3. **Abierto con el profesor:** P-16 (definición del ONI anual) y P-19 (registros 1930–1979).
4. **Informe y anexo de IA.**

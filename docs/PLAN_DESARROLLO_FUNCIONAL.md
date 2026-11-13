# Parte funcional: dónde está cada ítem de la rúbrica

El orden lo fija el enunciado: construir las curvas, hacer el EDA funcional y comparar las fuentes. Después viene la red funcional. La regresión funcional función a función es la referencia lineal; la cabeza lineal de la red debe reproducirla.

Todo el código está en `funcional/` y sus salidas en `salidas_R_funcional/`. Las figuras están en `graficos/` y las tablas en `tablas/`.

## 1. Construcción de datos funcionales (10 %): sección 11

| Ítem | Dónde | Resultado |
|---|---|---|
| Φ y MCO paso a paso, sin librerías | `00_funciones.R` (`base_bspline`, `ajuste_mco`, `construir`) | La base coincide con `splines::splineDesign` (3e-16); se usa solo para verificar |
| K por ECM, VC y VCG para cada fuente | `11_seleccion_K.png`, `11_K_elegido.csv` | K = 8 en estaciones y K = 11 en satélite, solo con entrenamiento |
| S idempotente y tr(S) = K | `11_controles_S.csv` | Error < 1e-14 en el peor ajuste |
| Covarianza residual y bandas punto a punto | `11_covarianza_residual.png`, `11_bandas_punto_a_punto.png` | Bandas con σ²I y σ²(t); efecto de un hueco de 15 días |
| Ambas fuentes en el mismo dominio | `11_mismo_dominio.png` | Día 1–365; las péntadas se ubican en su día central |
| Matriz de Gram W | Salida de consola de la sección 11 | 1'W1 = 364; c'Wc coincide con el trapecio fino |

## 2. EDA funcional (15 %): sección 12

| Ítem | Dónde |
|---|---|
| Media paso a paso, verificada con `fda::mean.fd` | `12_verificacion_media.csv`, `12_media_y_dispersion.png` |
| Varianza y DE | `12_media_y_dispersion.png`, `12_dispersion_resumen.csv` |
| Covarianza y correlación de cada fuente | `12_covarianza_*.png`, `12_correlacion_*.png` |
| Derivadas, inicio y fin de las temporadas | `12_derivada_media.png`, `12_temporadas_por_curva.png`, `12_temporadas_fechas.csv` |
| FM, MBD y boxplot funcional; atípicas por estación y año | `12_boxplot_funcional.png`, `12_atipicas_*.csv`; las atípicas coinciden con `fda::fbplot` |
| Interpretación: régimen bimodal, orografía y ENSO | `12_contexto.png`, `12_contexto_vertiente_altitud_ENSO.csv` y el bloque "Qué dicen los datos" |

Se usan las 176 estaciones y todos los años de 1981–2003; la representatividad está en la sección 8 de `lluvia_valle_datos.R`.

## 3. Comparación funcional (15 %): secciones 10 y 13

| Ítem | Dónde |
|---|---|
| Emparejamiento documentado | Sección 10: `10_emparejamiento_estacion_pixel.csv` y `10_pixeles_compartidos.csv` (D-38) |
| d(t) ± 2EE e intervalos sin el cero | `13_diferencia.png`, `13_intervalos_diferencia.csv`; el EE se calcula por estación |
| Correlación cruzada y su diagonal | `13_correlacion_cruzada.png`, `13_diagonal_correlacion.png` |
| Covarianzas comparadas | `13_covarianzas.png`, `13_covarianzas_comparadas.csv` |
| Por altitud, subregión y ENSO | `13_diferencia_por_grupo.png`, `13_comparacion_por_grupo.csv`, `13_sesgo_altitud.png` |
| Interpretación del sesgo | Bloque "Qué dicen los datos" de la sección 13 |

## 4. Regresión funcional y paso a la red: sección 14

| Qué | Resultado |
|---|---|
| FPCA de las estaciones, solo con entrenamiento | M = 5 (91 % de la varianza), elegido junto con λ en validación |
| β(s,t) = Σ β_m(s) ψ_m(t) con penalización ∫β'' | `14_beta_superficie.png`, `14_efectos.png` |
| Equivalencia entre la integral de β y los puntajes | Diferencia de 5e-15 (el criterio es ≤ 1e-4) |
| Prueba 2001–2003 frente al satélite sin corregir | `14_metricas_prueba.csv`: RMSE 2,43 → 2,24; KGE 0,38 → 0,45 |
| Insumos para Python | `para_python/`: pares, coeficientes en mm y en raíz, W, R, FPCA, coeficientes de la regresión y `especificacion.json` |

## 5. Lo que sigue

1. **Red funcional en Python** (`.py` único). Hay que adaptar `archivo_flujo_anterior/05_red_funcional.py` para que lea `salidas_R_funcional/para_python/`:
   - codificador z = A W c con 6 a 8 neuronas, más altitud y ONI;
   - decodificador FPCA congelado (μ y ψ de `fpca_estaciones.csv`);
   - cabezas lineal, KAN y MLP con parámetros igualados (±10 %) y 3 semillas.

   La cabeza lineal debe reproducir los coeficientes de `coeficientes_regresion.csv` o, al menos, la identidad de β(s,t).
2. **Unir el `.R`** (D-41): `lluvia_valle_datos.R` + `funcional/00_funciones.R` + secciones 10 a 14, sin los `source(...)`.
3. **Pendiente en equipo:** P-4, si 1985, 1988 y 2000 del satélite entran al modelo (D-37/D-42 frente a D-44). Lo demás quedó resuelto por el profesor: no hay estaciones reservadas, los 21 códigos sin ubicar se descartan y no se usa CHIRPS.
4. **Informe y anexo de IA.**

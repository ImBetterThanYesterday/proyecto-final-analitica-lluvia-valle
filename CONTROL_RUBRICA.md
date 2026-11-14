# Control de la rúbrica

Lista de todo lo que pide el enunciado (`enunciado/Proyecto_Final_Analitica_FeErratas.docx`) y cómo vamos.
Se actualiza cada vez que se avanza.

**Estados:** ✅ hecho y verificado · 🟡 a medias · ⬜ pendiente · ♻️ existe en `archivo_flujo_anterior/` con la población vieja (248 estaciones) y hay que rehacerlo con las 176.

**Última actualización:** 13-nov-2026, después de cerrar la parte funcional en `lluvia_valle_funcional.R` (D-48).

## Resumen

| Ítem | Peso | Estado | Qué falta |
|---|---|---|---|
| 1. Construcción de datos funcionales | 10 % | ✅ | Sección 2 de `lluvia_valle_funcional.R`: K = 8 / 11, S verificada, bandas |
| 2. EDA funcional | 15 % | ✅ | Sección 3; falta pasar la interpretación al informe |
| 3. Comparación funcional por fuente | 15 % | ✅ | Secciones 1 (emparejamiento) y 4, con todos los años |
| 4.1 Preprocesamiento y selección | 6 % | ✅ | Solo decidir sobre los registros 1930–1979 (con el FPCA) |
| 4.2 Red neuronal funcional | 10 % | 🟡 | La regresión de referencia está (sección 5); falta la red en Python |
| 4.3 Ablación KAN, MLP y lineal | 11 % | ♻️ | Rehacer |
| 4.4 Evaluación | 8 % | 🟡 | Las cuatro métricas, IC bootstrap, Wilcoxon y fase ENSO hechos para la regresión (sección 5.5); falta para la red |
| 4.5 Interpretación climatológica | 10 % | ♻️ | Rehacer |
| 5. Informe PDF | 10 % | ⬜ | Todo |
| 6. Scripts R y Python | 5 % | 🟡 | El R son dos archivos (datos y funcional) que se pegan para la entrega (D-48); falta el Python |

## Antes de lo funcional: datos y EDA (lo que ya está)

- ✅ **Estaciones del Valle:** 232 con coordenada, todas dentro del polígono del DANE (D-31 a D-35).
- ✅ **EDA de los datos crudos de estaciones:** secciones 1–6 de `lluvia_valle_datos.R`.
- ✅ **Depuración con reglas sacadas de los datos (D-36):** 176 estaciones y 3 299 pares estación–año.
- ✅ **Fuente del satélite identificada:** CHIRP solo-infrarrojo (D-43), confirmado por el profesor. Se usa solo ese, sin rellenar con otra fuente.
- ✅ **EDA del satélite y comparación con las estaciones (D-44).**
- ✅ **Datos preparados para la construcción funcional:** `salidas_R_datos/datos_preparados.rds`.
- ✅ **Pruebas:** 27 comprobaciones automáticas dentro del script (sección 14), incluida la R0 (D-50).

## 1. Construcción de datos funcionales (10 %) — sección 2 de `lluvia_valle_funcional.R`

- ✅ Matriz de la base Φ (Cox–de Boor propio) y coeficientes por mínimos cuadrados, sin librerías (sección 2 de `lluvia_valle_funcional.R`). Coincide con `splines::splineDesign` (3e-16), que se usa solo para verificar.
- ✅ K por ECM, VC y VCG para cada fuente, solo con entrenamiento: K = 8 (estaciones) y K = 11 (satélite) (D-39).
- ✅ K = 8 se sostiene frente a la fecha de los picos: el K que la mejora depende de la referencia elegida (D-51, `F2_K_que_gana_por_referencia.csv`).
- ✅ S idempotente y tr(S) = K (error < 1e-14).
- ✅ Covarianza residual y bandas punto a punto.
- ✅ Ambas fuentes en el dominio 1–365 (péntadas en su día central).
- 🟡 Días faltantes: el ajuste usa solo los días observados (sin imputar). Hay que escribirlo en el informe.
- ⚠️ **Tope:** construir las curvas con librerías deja la nota en máximo 2.

## 2. Análisis exploratorio funcional (15 %) — sección 3 de `lluvia_valle_funcional.R`

- ✅ Media funcional paso a paso, verificada con `fda::mean.fd`.
- ✅ Varianza y desviación estándar funcional.
- ✅ Covarianza y correlación funcional de cada fuente.
- ✅ Derivadas, inicio y fin de las temporadas de lluvia.
- ✅ Profundidad (FM y MBD) y boxplot funcional, con las atípicas por estación y año (coinciden con `fda::fbplot`).
- 🟡 Interpretación climatológica (régimen bimodal, orografía, ENSO): está en los comentarios del código; falta redactarla en el informe.
- ✅ **Tope de representatividad:** las 176 estaciones representan a las 220 (sección 8 de `lluvia_valle_datos.R`).

## 3. Comparación funcional por fuente (15 %) — secciones 1 y 4 de `lluvia_valle_funcional.R`

- ✅ Emparejamiento estación–píxel documentado: 176 estaciones en 148 píxeles, 23 compartidos (D-38, D-44).
- ✅ Curva de diferencia d(t) con banda ± 2EE e intervalos donde no contiene el cero.
- ✅ Correlación cruzada funcional y su diagonal.
- ✅ Comparación de las covarianzas.
- ✅ Comparación por altitud, vertiente y fase ENSO.
- 🟡 Interpretación del sesgo satelital: base en D-44, en la sección 12 de `lluvia_valle_datos.R` y en la sección 4 de `lluvia_valle_funcional.R`; falta redactarla.
- ✅ Todos los años entran (D-47, D-48). 1985, 1988 y 2000 se muestran además aparte, como descripción.
- ✅ Verificado que el diagnóstico global (diferencia media, diagonal, cociente de varianzas) no sale de usar K distinto por fuente; el veredicto día a día coincide en 281 de 365 días, así que las fechas de días significativos llevan esa salvedad (`F4_sensibilidad_mismo_K.csv`).

## 4.1 Preprocesamiento y selección de datos (6 %)

- ✅ Regla explícita para valores imposibles y errores de digitación: R0 (negativos y valores sobre el récord OMM de 1 825 mm en 24 h, contados aunque den 0; D-50) y R2 (extremos sin respaldo de las vecinas).
- ✅ Totales acumulados en un solo día (R3) y ceros que son faltantes (R4, regla calculada con la vecina más cercana y umbral de Bonferroni: 30 rachas, D-49).
- ✅ Criterio de cobertura mínima por estación y año, justificado con el experimento de huecos (R5 y R6).
- 🟡 Decisión sobre agregar datos:
  - imputación: no se imputa, y está justificado (D-13);
  - registros fuera del periodo (1930–1979) para estimar el FPCA: **falta decidirlo y escribirlo**.
- ✅ Tabla de trazabilidad: estaciones y pares antes y después de cada paso (`tablas/8_trazabilidad.csv`).
- ✅ Ninguna decisión usa datos de prueba: los umbrales (p99), el experimento de huecos y el mínimo de años se calculan solo con 1981–1995. Hay una prueba automática.

## 4.2 Red neuronal funcional (10 %)

- ♻️ Codificador con w_k(t) en B-splines y z = A W c, con la matriz de Gram por cuadratura (6 a 8 neuronas).
- ✅ Datos de altitud de la estación (`$meta_estaciones`) y ONI del año (`$pares`) listos para z en `datos_preparados.rds`, con prueba automática. Univalle (`2622320201`) toma los 996 m de la estación del IDEAM en el mismo campus (D-46); falta agregarlos a z al construir la red.
- ♻️ Decodificador FPCA solo con entrenamiento; M (4 a 6) justificado y varianza explicada.
- ♻️ Pérdida L² funcional con W y penalización de rugosidad; se modela sobre √P o log(1+P).
- ♻️ Demostración y verificación numérica (≤ 1e-4) de que la cabeza lineal equivale a la regresión funcional.
- ♻️ Reconstrucción no negativa en mm, con el % de valores recortados.

## 4.3 Ablación KAN, MLP y cabeza lineal (11 %)

- ♻️ Las tres cabezas con el mismo codificador y el mismo decodificador.
- ♻️ Parámetros de KAN y MLP contados con `sum(p.numel())`, diferencia ≤ 10 %.
- ♻️ Tres semillas por configuración, con media y desviación estándar.
- ♻️ Mismo presupuesto de entrenamiento para todas las cabezas.
- ♻️ KAN verificada con una función sintética (RMSE ≤ 1e-2).
- ♻️ Conclusiones solo si el IC bootstrap del 95 % de la diferencia excluye el cero.
- ⚠️ **Tope:** hiperparámetros elegidos con el conjunto de prueba dejan la nota en máximo 1.

## 4.4 Evaluación (8 %)

- ✅ Partición 1981–1995 / 1996–2000 / 2001–2003 aplicada en los datos preparados.
- ✅ Estaciones reservadas por el docente: **ya no van** (el profesor cambió el plan). La prueba es 2001–2003.
- ♻️ Comparación contra la curva satelital sin corregir.
- 🟡 Las cuatro métricas (error L² relativo, KGE, sesgo del acumulado anual en % y error en la fecha de los dos picos): hechas para la regresión de referencia (`F5_metricas_prueba.csv`); falta para la red.
- 🟡 IC bootstrap del 95 % (2 000 remuestreos por estación) y Wilcoxon pareado por estación: hechos para la regresión (`F5_regresion_vs_satelite_pareado.csv`); falta para la red.
- 🟡 Resultados por fase ENSO con el ONI: hechos para la regresión (`F5_metricas_por_fase_ENSO.csv`). La prueba no tiene año Niña; se agrega la validación como orientativa. Falta para la red.

## 4.5 Interpretación climatológica (10 %)

- ♻️ w_k(t) graficadas, con la ventana del año que usa cada neurona.
- ♻️ Funciones de activación de la primera capa de la KAN.
- ♻️ Contraste con los coeficientes de la cabeza lineal.
- ♻️ Ajuste simbólico de al menos una función de activación, con su R².
- ⬜ Lectura climatológica coherente con las gráficas.
- ⬜ Análisis de los peores casos con sus curvas.

## 5. Informe en PDF (10 %)

- ⬜ PDF a una columna, máximo 20 páginas con referencias, con todos los integrantes.
- ⬜ Cada tabla y figura numerada, titulada, citada e interpretada.
- 🟡 Referencias para los métodos y las fuentes de datos. Las fuentes están en `DECISIONES.md`; el satélite necesita la cita de CHIRP (D-43).
- 🟡 Anexo en Word con la herramienta de IA y los prompts. **Hay que ir guardando los prompts.**
- ⬜ Redacción propia.
- ⚠️ No se acepta R Markdown, dos columnas ni solo Colab. IA no declarada: máximo 2.

## 6. Scripts en R y Python (5 %)

- 🟡 Un único archivo de R y un único archivo de Python. Hoy el R son dos: `lluvia_valle_datos.R` y `lluvia_valle_funcional.R`. El segundo solo lee `datos_preparados.rds`, así que para la entrega se pegan uno tras otro (D-48). Falta el Python.
- ✅ Los dos R corren de principio a fin, con rutas relativas, semilla fija y versiones listadas (`sessionInfo.txt` en cada carpeta de salidas). Pruebas automáticas: 26 en datos y 12 en funcional.
- 🟡 Solo procedimientos vistos en clase; librerías adicionales solo donde se permite. En el R de datos las librerías son de EDA; la base B-spline del experimento de huecos es propia.
- ⬜ Cada tabla y figura del informe ubicable en el código.

## Reglas generales que hay que cuidar

- Entrega: PDF llamado "Proyecto final AnalíticaDeDatos", más un único R y un único Python.
- Sustentación de 40 minutos; el profesor puede preguntarle cualquier fase a cualquier integrante.
- Solo scripts, procedimientos y metodologías vistos en clase (`material_clase/`).

## Pendientes con el profesor

Detalle en `DECISIONES.md`, §8.

- ✅ P-8 y P-3: el satélite es CHIRP y se usa solo el entregado.
- ✅ P-1: los 21 sin ubicar se descartan.
- ✅ P-13: solo el Valle del Cauca.
- ✅ P-4: no se eliminan años del satélite (D-47), ya aplicado en `lluvia_valle_funcional.R` (D-48).
- ✅ P-12, P-14, P-15, P-17 y P-18: fuentes del IDEAM, Univalle, tope de 200 mm (serie graficada), 29 de febrero y píxel de la estación, aceptados en la asesoría (D-47).
- P-16 (definición del ONI anual) y P-19 (registros 1930–1979): siguen abiertas.

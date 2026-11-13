# Control de la rúbrica

Lista de todo lo que pide el enunciado (`enunciado/Proyecto_Final_Analitica_FeErratas.docx`) y cómo vamos.
Se actualiza cada vez que se avanza.

**Estados:** ✅ hecho y verificado · 🟡 a medias · ⬜ pendiente · ♻️ existe en `archivo_flujo_anterior/` con la población vieja (248 estaciones) y hay que rehacerlo con las 176.

**Última actualización:** 13-nov-2026.

## Resumen

| Ítem | Peso | Estado | Qué falta |
|---|---|---|---|
| 1. Construcción de datos funcionales | 10 % | ♻️ | Rehacer con los datos preparados |
| 2. EDA funcional | 15 % | ♻️ | Rehacer con los datos preparados |
| 3. Comparación funcional por fuente | 15 % | 🟡 | El emparejamiento está; faltan las comparaciones funcionales |
| 4.1 Preprocesamiento y selección | 6 % | ✅ | Solo decidir sobre los registros fuera del periodo |
| 4.2 Red neuronal funcional | 10 % | ♻️ | Rehacer con los pares nuevos |
| 4.3 Ablación KAN, MLP y lineal | 11 % | ♻️ | Rehacer |
| 4.4 Evaluación | 8 % | ♻️ | Rehacer con la prueba 2001–2003 (ya no hay estaciones reservadas) |
| 4.5 Interpretación climatológica | 10 % | ♻️ | Rehacer |
| 5. Informe PDF | 10 % | ⬜ | Todo |
| 6. Scripts R y Python | 5 % | 🟡 | El R de datos está; falta lo funcional en R y el único archivo Python |

## Antes de lo funcional: datos y EDA (lo que ya está)

- ✅ **Estaciones del Valle:** 232 con coordenada, todas dentro del polígono del DANE (D-31 a D-35).
- ✅ **EDA de los datos crudos de estaciones:** secciones 1–6 de `lluvia_valle_datos.R`.
- ✅ **Depuración con reglas sacadas de los datos (D-36):** 176 estaciones y 3 326 pares estación–año.
- ✅ **Fuente del satélite identificada:** CHIRP solo-infrarrojo (D-37), confirmado por el profesor. Se usa solo ese, sin rellenar con otra fuente.
- ✅ **EDA del satélite y comparación con las estaciones (D-38).**
- ✅ **Datos preparados para la construcción funcional:** `salidas_R_datos/datos_preparados.rds`.
- ✅ **Pruebas:** 26 comprobaciones automáticas dentro del script (sección 14). Todas pasan.

## 1. Construcción de datos funcionales (10 %)

- ♻️ Matriz de la base Φ y coeficientes por mínimos cuadrados, paso a paso y sin librerías.
- ♻️ Selección de K (o de λ) con ECM, VC y VCG: gráfica y valor elegido para cada fuente.
- ♻️ Verificación numérica de la matriz de suavizado: S idempotente y tr(S) = K.
- ♻️ Covarianza residual y bandas de confianza punto a punto.
- 🟡 Curvas de ambas fuentes sobre el mismo dominio anual. Los datos ya están en el dominio 1–365 (estaciones por día, satélite por día central de la péntada); falta construir las curvas.
- ⚠️ **Tope:** construir las curvas con librerías deja la nota en máximo 2.

## 2. Análisis exploratorio funcional (15 %)

- ♻️ Media funcional paso a paso y verificada contra la función de librería.
- ♻️ Varianza y desviación estándar funcional.
- ♻️ Covarianza y correlación funcional de cada fuente.
- ♻️ Derivadas e interpretación del inicio y el fin de las temporadas de lluvia.
- ♻️ Profundidad (Fraiman–Muniz y MBD) y boxplot funcional, con las atípicas por estación y año.
- ⬜ Interpretación climatológica de cada resultado (régimen bimodal, orografía, ENSO).
- ✅ **Tope de representatividad:** las 176 estaciones representan a las 220 (misma altitud mediana, lluvia parecida; sección 8). Hay que citarlo en el informe.

## 3. Comparación funcional por fuente (15 %)

- ✅ Emparejamiento estación–píxel documentado: coordenadas, regla (el píxel que contiene la coordenada) y estaciones que comparten píxel (23 píxeles con más de una). Sección 12 y `tablas/12_pares_estacion_pixel.csv`.
- ⬜ Curva de diferencia d(t) = estación − satélite con banda d̄(t) ± 2EE(t) e intervalos donde la banda no contiene el cero.
- ⬜ Correlación cruzada funcional y su diagonal.
- ⬜ Comparación de la covarianza de cada fuente.
- 🟡 Comparación por altitud o subregión y por fase ENSO. El EDA ya muestra el efecto (D-38); falta hacerlo con las curvas.
- 🟡 Interpretación climatológica del sesgo satelital. La base está en D-38: el satélite aplana el espacio y no capta los años Niño y Niña.

## 4.1 Preprocesamiento y selección de datos (6 %)

- ✅ Regla explícita para valores imposibles y errores de digitación (R2: extremos sin respaldo de las vecinas).
- ✅ Totales acumulados en un solo día (R3) y ceros que son faltantes (R4, con revisión a mano).
- ✅ Criterio de cobertura mínima por estación y año, justificado con el experimento de huecos (R5 y R6).
- 🟡 Decisión sobre agregar datos:
  - imputación: no se imputa, y está justificado (D-13);
  - registros fuera del periodo (1930–1979) para estimar el FPCA: **falta decidirlo y escribirlo**.
- ✅ Tabla de trazabilidad: estaciones y pares antes y después de cada paso (`tablas/8_trazabilidad.csv`).
- ✅ Ninguna decisión usa datos de prueba: los umbrales (p99), el experimento de huecos y el mínimo de años se calculan solo con 1981–1995. Hay una prueba automática.

## 4.2 Red neuronal funcional (10 %)

- ♻️ Codificador con w_k(t) en B-splines y z = A W c, con la matriz de Gram por cuadratura (6 a 8 neuronas).
- ✅ Datos de altitud de la estación (`$meta_estaciones`) y ONI del año (`$pares`) listos para z en `datos_preparados.rds`, con prueba automática. Univalle (`2622320201`) toma los 996 m de la estación del IDEAM en el mismo campus (D-40); falta agregarlos a z al construir la red.
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
- ♻️ Las cuatro métricas: error L² relativo, KGE, sesgo del acumulado anual y error en la fecha de los dos picos.
- ♻️ IC bootstrap del 95 % (≥ 1000 remuestreos) y prueba de Wilcoxon pareada por estación.
- 🟡 Resultados por fase ENSO con el ONI y su fuente citada. El ONI está en `datos_web/`; falta aplicarlo a los resultados.

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
- 🟡 Referencias para los métodos y las fuentes de datos. Las fuentes están en `DECISIONES.md`; el satélite necesita la cita de CHIRP (D-37).
- 🟡 Anexo en Word con la herramienta de IA y los prompts. **Hay que ir guardando los prompts.**
- ⬜ Redacción propia.
- ⚠️ No se acepta R Markdown, dos columnas ni solo Colab. IA no declarada: máximo 2.

## 6. Scripts en R y Python (5 %)

- 🟡 Un único archivo de R y un único archivo de Python. Hoy: `lluvia_valle_datos.R` (datos). Falta sumarle lo funcional y armar el Python.
- ✅ El R corre de principio a fin, con rutas relativas, semilla fija y versiones listadas (`salidas_R_datos/sessionInfo.txt`).
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
- P-12: fuentes externas (ONI, CNE, DANE).
- ✅ P-13: solo el Valle del Cauca.
- P-4: el año 2000 del satélite (ya no es urgente, D-38).

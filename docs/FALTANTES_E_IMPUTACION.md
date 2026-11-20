# Datos faltantes: qué dice la literatura, cómo se rellenan y qué nos conviene

> **Actualización del 20-nov-2026 (D-61).** La recomendación de la sección 6 se revisó:
> - **R7 ahora rellena** los días faltantes de los años que pasan R5 con **kriging ordinario de las vecinas**, cada una escalada por la razón de medias mensuales. La razón: la OMM (2017, 4.6) dice que, para valores diarios, los días de antes y después «are rarely useful» y pone primero la interpolación espacial; y en nuestros datos rellenar con vecinas deja menos error en la curva que dejar que la B-spline cubra el hueco (0,30 frente a 0,55 con el máximo que acepta R5).
> - **R5 no cambia** (los límites de la OMM) y **R6 tiene dos conjuntos**. Desde D-63 el principal es el de **≥10 años** (la opción B de abajo, que el profesor prefirió), y el del 80 % de la OMM queda como sensibilidad.
> - Lo demás se mantiene: no se rellenan años sin datos ni huecos largos, y nunca con el satélite.
> - Detalle, fuentes y pruebas: `docs/FUNDAMENTOS_REGLAS.md`, secciones 2, 3 y 6.

Este documento responde tres preguntas de Gus (19-nov-2026):
1. ¿La regla de la OMM es la fuente oficial o hay más literatura?
2. ¿Cómo se rellenan (imputan) los datos faltantes, y cómo lo hace NOAA?
3. ¿Qué opciones tenemos, en qué se diferencian y qué nos tocaría hacer en cada una?

Las cifras son de nuestros datos: `salidas_R_datos/lluvia_valle_limpia.rds`, corrida del 19-nov-2026. Las reglas actuales están en `docs/FUNDAMENTOS_REGLAS.md`.

---

## 1. ¿Qué tan "oficial" es la OMM y qué más hay?

**La OMM es la fuente más alta en climatología.** Es la agencia de la ONU que fija las normas para los servicios meteorológicos nacionales, entre ellos el IDEAM. Su guía para datos faltantes es la *WMO Guidelines on the Calculation of Climate Normals* (WMO-No. 1203, 2017, sección 4.4). Pero:
- **es una guía para calcular normales**, es decir, promedios climáticos de 30 años, no para entrenar modelos;
- **cada país la adapta.** La misma guía dice que los países que usan reglas más estrictas pueden seguir usándolas.

Las otras fuentes encontradas muestran que **no hay un único umbral universal**: cada comunidad fija el suyo según para qué usa los datos.

| Fuente | Para qué | Regla de días dentro del año o el mes | Regla de años por estación |
|---|---|---|---|
| **OMM 2017** (WMO-No. 1203, 4.4) | Normales climatológicas | Mes válido si faltan ≤10 días y menos de 5 seguidos. En totales de lluvia, solo si esos días se estimaron | ≥80 % de los años del periodo, sin criterio de años seguidos |
| **OMM 1989** (WMO-TD 341) | Normales (versión anterior, más estricta) | ≤5 faltantes y ≤3 seguidos por mes | — |
| **NOAA, normales 1991–2020** | Normales de EE. UU. | Mes completo, aceptando acumulados de 2 y 3 días | "S" con ≥80 % de los años (24 de 30); "R" o "P" con ≥10 años. Según NOAA, la OMM pide 10 años como mínimo para cualquier normal |
| **ETCCDI** (climdex.pcic, RClimDex) | Índices de extremos climáticos | Mes válido con ≤3 días faltantes; año válido con ≤15 faltantes | — |
| **Kim y Ryu (2015)**, *Water Resources Management* | Rellenar series diarias | El relleno funciona hasta ~15 % de faltantes | — |
| Estudio de 550 pluviómetros en Malasia | Selección de estaciones | Se excluyen estaciones con >10 % de faltantes | — |
| Validación de CHIRPS en la Orinoquía (2025) | Validar satélite en Colombia | Se excluyen estaciones con >30 % de meses faltantes | — |

**Conclusión.** La OMM es la referencia de mayor jerarquía, y por eso la usamos. Pero su 80 % de años es para un propósito distinto al nuestro, y otras fuentes serias usan reglas más sueltas. Por eso conviene mostrar **cómo cambian los resultados con otra regla publicada**, en vez de tratar el 80 % como la única respuesta correcta.

---

## 2. Cómo rellena NOAA los huecos

Viene del documento de metodología de las normales 1991–2020 de NOAA y de su apéndice.

**1. NOAA no rellena días: rellena totales mensuales.** Para calcular un promedio mensual, primero arma el total de cada mes. Acepta acumulados de 2 y 3 días, y rellena solo los meses que faltan.

**2. El relleno es una regresión con estaciones vecinas:**
- se arma por separado para cada estación y cada mes del calendario (un modelo para los eneros, otro para los febreros…);
- las vecinas deben estar a menos de **500 km**, tener al menos **10 años** de ese mes en común con la estación y tener dato en todos los años que hay que rellenar;
- las vecinas se ordenan por su **índice de concordancia** con la estación, que mide a la vez correlación y amplitud;
- se empieza con la mejor vecina y se agrega una más mientras la concordancia mejore, hasta un máximo de **5**;
- es una regresión robusta, basada en la mediana de las desviaciones absolutas;
- si una estimación da negativa, se pone en 0.

**3. Cuándo no rellena:**
- si la estación tiene todos los totales de ese mes en 0;
- si no hay vecinas que cumplan;
- si la regresión no tiene solución única.

**4. Qué consigue.** Una estación con al menos 10 años, pero sin el 80 %, recibe la marca **"R" (representativa)** si se pudieron rellenar los años faltantes, o **"P" (provisional)** si no.

**¿Sirve para nosotros?** No directamente. NOAA rellena **totales mensuales** para sacar un promedio. Nosotros necesitamos la **forma diaria** de la curva de cada año, con sus dos picos, y un total mensual no la da. Lo que sí se puede copiar es la idea: un año que falta se puede reemplazar por una regresión con vecinas bien correlacionadas.

---

## 3. Métodos para rellenar lluvia diaria

**Antes que nada, un principio de la estadística** (Little y Rubin, *Statistical Analysis with Missing Data*): rellenar solo es válido si **la razón por la que falta el dato no depende del valor que falta**.
- **Una estación que no operaba** (todavía no existía o ya se había cerrado) no depende de cuánto llovió. Ahí no hay sesgo.
- **Un pluviómetro que nadie leyó** puede depender de la lluvia: por ejemplo, si no se lee en los aguaceros. Ahí rellenar puede meter sesgo.

| Familia | Cómo estima el día que falta | Ventaja | Problema | Fuentes |
|---|---|---|---|---|
| Razón normal / distancia inversa (IDW) | Promedio de las vecinas, ponderado por su lluvia media o por la distancia | Simple; muy usado en Colombia y Latinoamérica | En montaña, la vecina más cercana no es la más parecida | Estudios del Quindío, Chile, México y Argentina |
| Ponderar por correlación | Igual, pero el peso es la correlación con la vecina, no la distancia | Mejora a IDW en las comparaciones | Necesita años en común para calcular la correlación | Teegavarapu y Chandramouli (2005) |
| Regresión con vecinas | Regresión de la estación sobre sus mejores vecinas | Suele ganar en las comparaciones; es lo que usa NOAA | Aplana los extremos: la serie rellenada varía menos que la real | Wagner et al. (2012); NOAA |
| Kriging / cokriging | Interpolación espacial con el variograma | Da el error de cada estimación | La lluvia diaria tiene muchos ceros y no es normal | Revisiones comparativas |
| Imputación múltiple (MICE, Amelia) | Genera varias series plausibles y combina los resultados | Lleva la incertidumbre del relleno a los resultados | Más complejo; supone que la razón del faltante no depende de la lluvia | Estudio de Boyacá (Amelia), Carvalho et al. (2017) |
| Aprendizaje automático (redes, bosques) | Aprende la relación con las vecinas | A veces gana en precisión | Caja negra y riesgo de sobreajuste | Estudios de Irlanda y del Huila |
| **Completar curvas parciales** (datos funcionales) | Predice el tramo que falta de una curva a partir de su parte observada y de cómo varían todas las curvas (FPCA) | **Es del mismo marco que el proyecto** y da bandas para el tramo rellenado | Pide suficientes curvas completas para estimar la covarianza | **Kraus (2015)**; Yao, Müller y Wang (2005) |
| Con satélite o radar | El satélite estima el día que falta | La OMM lo menciona como opción | **En nuestro proyecto es circular:** el satélite es lo que el modelo usa para predecir. Rellenar la estación con satélite le enseñaría al modelo a copiarse a sí mismo | OMM 2017, 4.6 |

**Ninguna comparación encuentra un método que gane siempre.** Depende de la región, de la densidad de estaciones y de qué tanto se parecen las vecinas.

**Límites que pone la OMM al relleno** (2017, secciones 4.4.1 c y 4.6):
- solo para **huecos pequeños**, de hasta 10 días en un mes;
- un total mensual con días estimados vale solo si esos días no son 11 o más, ni 5 o más seguidos;
- **no se debe estimar lo que la estación nunca observó**;
- el relleno debe ser **sin sesgo**. Por ejemplo, no se puede usar el satélite solo para confirmar los días secos.

---

## 4. Qué pasa con nuestros datos

**¿Dónde se pierden los datos?** De 220 estaciones × 23 años = 5 060 años de estación posibles:

| Situación | Años de estación | ¿Se puede rellenar? |
|---|---|---|
| Válidos con la regla de la OMM | 3 279 | No hace falta |
| **La estación no operaba ese año** (sin ningún dato) | **1 117** | **No.** Sería reconstruir un año entero que la estación nunca midió; la OMM lo desaconseja |
| La estación operaba, pero el año no pasa R5 | 664 | En parte |
| — de ellos, con ≥90 % de los días | 332 | El hueco mayor tiene mediana de 14 días; solo 138 tienen huecos de 10 días o menos |
| — con 80 % a 90 % de los días | 106 | Huecos más largos |
| — con menos de 80 % | 226 | No razonable |

**Lo que muestra:**
1. **La mayor pérdida son años enteros en que la estación no existía.** Ningún método serio rellena eso desde la propia estación. El profesor dijo lo mismo en la revisión: reconstruir unos 1 000 días seguidos «eso no se lo reconstruye nadie».
2. **Los años que fallan R5 tienen huecos largos** (mediana de 14 días, y 31 días en el 10 % peor). La OMM permite estimar como máximo 4 días seguidos para que el mes cuente. Rellenarlos para que pasen la regla sería ir **más allá** de lo que la OMM acepta.
3. **Lo que la OMM sí acepta ya lo hacemos.** La curva B-spline estima los días faltantes de los meses que cumplen R5. Esa es la imputación permitida.

**¿Cuánto se ganaría rellenando?** Si se rellenaran los años con al menos 80 % de los días y se contaran como válidos, con la R6 del 80 %:
- quedarían **136 estaciones (12 pacíficas)** y 2 998 curvas, frente a 98 (7 pacíficas) y 2 066;
- pero serían huecos de 2 a 7 semanas rellenados con vecinas, fuera de la guía de la OMM.

---

## 5. Las opciones, lado a lado

| | A. OMM estricta (actual) | B. Mínimo de 10 años | C. Rellenar huecos largos con vecinas | D. Completar curvas parciales (FDA) |
|---|---|---|---|---|
| **Regla** | R5 de la OMM + ≥80 % de los años | R5 de la OMM + ≥10 años válidos | R5 relajada con relleno + ≥80 % de los años | R5 relajada con completación + ≥80 % de los años |
| **Estaciones (pacíficas)** | 98 (7) | 163 (22) | ~136 (~12) | ~136 (~12) |
| **Curvas** | 2 066 | 3 029 | ~2 998 | ~2 998 |
| **¿Rellena datos?** | No (solo lo que estima la base, dentro de R5) | No | Sí: días de huecos de 2 a 7 semanas | Sí: tramos de curva |
| **Respaldo** | OMM 2017, 4.4 | Práctica de NOAA para normales "P"; NOAA atribuye los 10 años mínimos a la OMM | NOAA (mensual), Teegavarapu y Chandramouli (2005); fuera del límite de la OMM | Kraus (2015); Yao, Müller y Wang (2005) |
| **A favor** | Norma de mayor jerarquía; responde al profesor | No inventa datos; recupera casi toda la zona pacífica | Recupera años casi completos | Del mismo marco del curso; da bandas de la parte rellenada |
| **En contra** | La zona pacífica queda con 7 estaciones | Hay estaciones con bloques de años vacíos, lo que el profesor cuestionó | Inventa datos y aplana los extremos; en montaña la vecina no siempre se parece | Más trabajo; también inventa datos |
| **Qué habría que hacer** | Nada: es lo que está | Cambiar el 19 por 10 en R6 y volver a correr (~25 min) | Programar el relleno, validarlo con días que sí existen, y volver a correr | Programar la completación (predicción de los puntajes FPCA), validarla y volver a correr |

---

## 6. Recomendación

1. **Dejar A como análisis principal.** Es la norma de mayor jerarquía, responde a lo que pidió el profesor y no inventa datos.
2. **Agregar B como sensibilidad.** Tiene respaldo publicado, no rellena nada y recupera la zona pacífica (22 estaciones). Cada curva es un año completo de una estación, así que un año que falta es solo una curva menos, no un hueco que haya que llenar.
   - Si el resultado es parecido con A y con B, la conclusión no depende de la regla.
   - Si cambia, sobre todo en la zona pacífica, eso es un hallazgo que se reporta.
3. **No rellenar huecos largos (C ni D) en el análisis principal.** La propia OMM lo desaconseja, la mayoría de lo que falta son años en que la estación no existía, y lo que la norma sí permite ya lo hace la curva.
   - D se puede mencionar como trabajo futuro: es la forma propia de datos funcionales de usar años incompletos.
4. **Nunca rellenar con el satélite:** sería circular.

---

## Referencias

- Carvalho et al. (2017) y Wagner et al. (2012): citados en las revisiones comparativas de métodos de relleno, por ejemplo *Cross Assessment of Twenty-One Different Methods for Missing Precipitation Data Estimation* (*Atmosphere*, 2020, 11(4), 389). No se revisó el texto original; hay que hacerlo antes de citarlos en el informe.
- Kim, J. y Ryu, J. H. (2015). Quantifying a threshold of missing values for gap filling processes in daily precipitation series. *Water Resources Management*, 29(11), 4173–4184.
- Kraus, D. (2015). Components and completion of partially observed functional data. *Journal of the Royal Statistical Society: Series B*, 77(4), 777–801.
- Little, R. J. A. y Rubin, D. B. (2019). *Statistical Analysis with Missing Data*, 3.ª ed. Wiley.
- NOAA NCEI (2021). *Normals Calculation Methodology 2020*, con el apéndice "Procedure for estimating monthly precipitation totals".
- Teegavarapu, R. S. V. y Chandramouli, V. (2005). Improved weighting methods, deterministic and stochastic data-driven models for estimation of missing precipitation records. *Journal of Hydrology*, 312, 191–206.
- WMO (2017). *WMO Guidelines on the Calculation of Climate Normals*. WMO-No. 1203 (secciones 4.4 y 4.6).
- Yao, F., Müller, H.-G. y Wang, J.-L. (2005). Functional data analysis for sparse longitudinal data. *Journal of the American Statistical Association*, 100(470), 577–590.
- Las demás fuentes (ETCCDI y climdex.pcic, el estudio de Malasia, la validación de la Orinoquía y los estudios de relleno de Colombia y Latinoamérica) están enlazadas en la conversación del 19-nov-2026 y en `docs/FUNDAMENTOS_REGLAS.md`.

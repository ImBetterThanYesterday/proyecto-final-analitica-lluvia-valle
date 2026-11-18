# Base de conocimiento: cómo trabaja y presenta el profesor

Este documento reúne lo que el profesor Johann Alexis Ospina-Galindez enseña, usa y espera ver en el análisis de datos funcionales, y cómo lo aplicamos al proyecto. Sirve para dos cosas:

- presentar el informe y la sustentación como él lo hace;
- defender cada decisión con el material del curso.

## 0. Fuentes revisadas

| Fuente | Qué es | Para qué sirve |
|---|---|---|
| Tesis de maestría (Univalle, 2016): *Modelo de regresión funcional de la precipitación pentadal, integrando datos de estaciones climáticas y de satélite* | El trabajo del profesor sobre **el mismo problema**: lluvia del Valle del Cauca, estaciones y satélite, 1993–2011 | Cómo plantea, grafica y concluye un estudio como el nuestro |
| `Introducción_ADF.pdf` y `ConstruccionDatosFuncionales.pdf` (semanas 10–11) | Teoría y construcción: mínimos cuadrados, K, bandas | Construcción (rúbrica 1) |
| `AnalisisExploratorioDatosFuncionales.pdf` y `Taller_Estudio_Independiente_ADF.pdf` (semanas 12–13) | EDA funcional con SO₂ y con `X_sat` y `Y_est` | EDA (rúbrica 2) |
| `InferenciaDatosFuncionales.pdf` (semana 14) | Prueba t funcional y ANOVA funcional | Comparación por fuente (rúbrica 3) |
| `RegresionFuncional.pdf` (semana 15) | Escalar-funcional, funcional-escalar, concurrente y funcional-funcional | Regresión (rúbrica 4) |
| Los cinco scripts de clase | El código de cada sesión | Qué funciones y figuras usa |
| Tres videos de clase (4 h 37 min) | Las sesiones grabadas | Lo que dice que espera; ver sección 7 |

---

## 1. La tesis del profesor, en resumen

**Problema.** Modelar la precipitación pentadal del Valle del Cauca entre 1993 y 2011 con 82 estaciones y el satélite. Los datos de satélite se los dieron FEWS NET y la Universidad de California en Santa Bárbara: es la familia CHIRP/CHIRPS, la misma que usamos nosotros.

**Pasos, en orden:**

1. **Datos faltantes:** los imputa con kriging ordinario año por año. Elige el semivariograma (exponencial, esférico, gaussiano o Matérn) por el menor ECM de validación cruzada.
2. **EDA de los datos crudos.** Compara las dos fuentes:
   - las estaciones son más variables;
   - el satélite es más homogéneo porque es "una estimación uniforme en toda la región".
   - Mide la correlación lineal estación–satélite por año (media de 0,33 a 0,67) y compara valle frente a montaña.
3. **Construcción.**
   - Trabaja la precipitación **acumulada** por pentadas (72 por año).
   - Usa B-splines **monótonas crecientes** (`smooth.monotone`), sin penalización.
   - Elige K por validación cruzada generalizada año por año: K = 10 en estaciones y K = 20 en satélite.
   - Muestra cada curva con su banda del 95 %.
4. **EDA funcional:** media y desviación estándar por año y fuente, por zona (valle o montaña).
5. **Regresión funcional concurrente:** precipitación de la estación(t) = β₀(t) + β₁(t)·satélite(t) + β₂(t)·altura + β₃(t)·altura² + ε(t).
   - Bandas de confianza con `y2cMap` y `SigmaE` (Ramsay y otros, 2011).
   - Bondad de ajuste con R²(t) y validación cruzada dejando una estación fuera.
6. **Correlación espacial de los residuos:**
   - traza-variograma de los residuos funcionales (Giraldo, 2009);
   - ajuste exponencial, esférico, gaussiano o Matérn por mínimos cuadrados;
   - el exponencial es el que más veces gana: 9 de los 19 años.
7. **Modelo final:** el concurrente con **mínimos cuadrados generalizados**, usando la matriz de covarianza espacial de los errores. En el apéndice G lo programa como un modelo funcional–funcional con β(s,t).
8. **Evaluación:**
   - R²(t) con y sin correlación espacial;
   - suma de cuadrados del error por validación cruzada, comparada con Mann–Whitney;
   - dispersión de las 82 curvas ajustadas frente a observadas;
   - 20 curvas observadas frente a estimadas;
   - mapas de la RECM por estación.

**Conclusiones de la tesis** (útiles para contrastar con las nuestras):

- El satélite **no capta la orografía**: la media de las estaciones es mayor en la montaña y la del satélite no cambia entre valle y montaña.
- La altura y la altura² **no aportan** en la regresión: sus bandas contienen el cero casi siempre.
- El satélite aporta más en los **años Niña** y en las temporadas húmedas (pentadas 18–30, marzo–mayo, y 54–66, septiembre–noviembre).
- Modelar la correlación espacial de los errores mejora la predicción en el 53–63 % de los años.
- Los residuos tienen tendencia al final del año, por la acumulación.

---

## 2. Cómo presenta: reglas que se repiten en todo el material

1. **Cada figura lleva una frase con números.** Nunca dice solo "se observa un patrón". Ejemplos de sus diapositivas:
   - "La media alcanza su máximo, cercano a 10,5, alrededor del registro 11,6."
   - "La correlación en la diagonal varía entre 0,10 y 0,62, con mediana de 0,54."
   - "Al final del año acumula 1 543,7 mm frente a 1 345,1 mm."
2. **Lo programado paso a paso se compara con la librería** y se dice la diferencia. Por ejemplo, la media calculada a mano frente a `mean.fd()` "coinciden: la diferencia máxima es del orden de 10⁻¹⁵".
3. **Las bandas son punto a punto, y lo advierte siempre:** "no garantizan la cobertura simultánea de toda la curva".
4. **El EDA responde seis preguntas** (diapositiva 4 de la semana 12):
   - ¿cuál es la curva típica?
   - ¿dónde se concentra la variabilidad?
   - ¿cómo se relacionan dos instantes?
   - ¿cómo se relacionan dos variables?
   - ¿a qué ritmo cambia?
   - ¿hay curvas atípicas?
5. **Los atípicos se confirman con el contexto:** "estos señalamientos deben confirmarse con el contexto de los datos".
6. **Primero mira dónde están las diferencias y luego decide globalmente.** En inferencia: curva t(t) o F(t), luego la prueba global y luego la región significativa con el criterio del máximo.
7. **Interpreta en el dominio:**
   - la derivada de la lluvia acumulada es "la intensidad de lluvia en mm por pentada";
   - la diagonal de la correlación cruzada es "la concordancia instante a instante entre ambas fuentes";
   - β(s,t) es "cuánto influye el valor de x en el instante s sobre la respuesta en el instante t".
8. **Separa por grupos con sentido físico:** valle y montaña, regiones climáticas (Canadá), fase ENSO (tesis).
9. **Los ejercicios de clase son el proyecto en pequeño.** Con `X_sat` y `Y_est` (10 sitios, 2011, acumuladas, K = 20) pide:
   - **EDA:** comparar media y mediana, la memoria temporal de cada fuente, `smooth.monotone` y la correlación cruzada de las intensidades.
   - **Inferencia (ejercicio 3):** "aplique la prueba F para muestras pareadas y compare con la banda d̄(t) ± 2 EE(t)".
   - **Regresión (ejercicio 3):** "ajuste el modelo concurrente Y_est(t) = β₀(t) + β₁(t) X_sat(t). ¿En qué parte del año el satélite explica mejor las mediciones en tierra?"

---

## 3. Catálogo de figuras que usa

| Etapa | Figura | Dónde la usa | En el proyecto |
|---|---|---|---|
| Construcción | Curvas con K = 4, 8, …, 24 | Semana 10–11 | `F2_seleccion_K` (criterios) |
| Construcción | ECM, VC y VCG frente a K | Semana 10–11, tesis 6.5 | `F2_seleccion_K` |
| Construcción | Cada curva con datos y banda del 95 % | Semana 10–11, tesis fig. 6.1 | `F2_bandas_punto_a_punto` |
| EDA | Curvas, media ± 2 DE y DE funcional | Semana 12, tesis figs. 6.2–6.3 | `F3_media_y_dispersion` |
| EDA | **CV funcional CV(t) = DE(t)/media(t)** | Semana 12 | `F3_media_y_dispersion`, 4.º panel (**nuevo**) |
| EDA | Covarianza y correlación (superficie o curvas de nivel) | Semana 12 | `F3_covarianza_*`, `F3_correlacion_*` |
| EDA | Correlación cruzada y su diagonal | Semana 12 | `F4_correlacion_cruzada`, `F4_diagonal_correlacion` |
| EDA | Derivada de la media, orden 2 frente a orden 4 | Semana 12 | `F3_derivada_media` |
| EDA | Diferencias por sitio y d̄(t) ± 2 EE | Semana 12 | `F4_diferencia` |
| EDA | Curvas ordenadas por profundidad; media, mediana y media recortada | Semana 12 | `F3_centralidad` |
| EDA | Boxplot funcional | Semana 12 | `F3_boxplot_funcional` |
| EDA | **Curvas acumuladas por fuente y d(t) acumulada** | Tesis, semana 12 | `F4_acumulada` (**nuevo**) |
| Inferencia | **t(t) con crítico puntual y crítico del máximo** | Semana 14 | `F4_prueba_t_pareada` (**nuevo**) |
| Inferencia | **F(t) del ANOVA funcional y tabla con gl × d̂** | Semana 14 | `F4_anova_funcional` (**nuevo**) |
| Regresión | OCV frente a log λ | Semana 15 | `F5_seleccion_M_lambda` |
| Regresión | β(t) con bandas ± 2 EE | Semana 15, tesis figs. 7.1, 9.1 y 9.2 | `F5_efectos` con bandas (**nuevo**), `C2_*`, `C6_*` |
| Regresión | **Panel del concurrente: x(t), y(t), R²(t), β₀, β₁ y correlación de residuos** | Semana 15 | `C6_concurrente_estilo_clase` (**nuevo**) |
| Regresión | β(s,t) en 3D y curvas de nivel con la diagonal | Semana 15 | `F5_beta_superficie` |
| Regresión | Ajuste en varias estaciones (datos, curva y modelo) | Semana 15, tesis fig. 10.3 | `F5_prueba`, **`F5_curvas_20_estaciones`** (**nuevo**) |
| Regresión | R²(t) | Semana 15, tesis fig. 10.1 | `F5_prueba`, `C3_R2_por_dia` |
| Regresión | **Prueba F por permutación (`Fperm.fd`)** | Semana 15 | `F5_permutacion_satelite` (**nuevo**) |
| Diagnóstico | **Residuos estudentizados y residuos frente a predichos** | Tesis figs. 9.3–9.4 | `F5_diagnostico` (**nuevo**) |
| Diagnóstico | **Traza-variograma de los residuos** | Tesis cap. 8 | `F5_diagnostico` (**nuevo**) |
| Evaluación | **Mapa del error por estación** | Tesis figs. 10.4–10.6 | `F5_mapa_error` (**nuevo**) |
| Datos, por año | **Precipitación acumulada por año, estación y satélite** | Tesis figs. 5.8 y 5.12 | `12_acumulada_por_anio` (**nuevo**, D-59) |
| Datos, por año | **Correlación estación–satélite de cada año** | Tesis cap. 5 | `12_correlacion_por_anio` (**nuevo**, D-59) |
| EDA, por año | **Curvas de cada año con promedio y desviación estándar** | Tesis figs. 6.2 y 6.3 | `F3_por_anio_estaciones`, `F3_por_anio_satelite` (**nuevo**, D-59) |
| Regresión, por año | **b₁(t) del concurrente ajustado en cada año, con banda** | Tesis fig. 9.2 | `C7_b1_por_anio` (**nuevo**, D-59) |

## 4. Funciones de R que enseña

- **Construcción:** `create.bspline.basis`, `smooth.basis`, `fdPar`, `eval.fd`, `smooth.monotone` y `optim.basis` (de fda.usc).
- **EDA:** `mean.fd`, `std.fd`, `var.fd`, `cor.fd`, `deriv.fd`, `fbplot`, `depth.FM` y `depth.mode` (de fda.usc).
- **Inferencia:**
  - `tperm.fd` para la t por permutación con criterio del máximo;
  - `Fperm.fd`;
  - el paquete `fdANOVA`, con las pruebas L², F, GPF, Fmax y la de Cuevas–Febrero–Fraiman.
- **Regresión:**
  - `fRegress`, con `y2cMap` y `SigmaE` para los errores estándar;
  - `fRegress.CV`;
  - `linmod` para el funcional–funcional;
  - `pca.fd`.

El proyecto programa estos métodos paso a paso y usa las funciones de librería para verificar, como pide el curso.

---

## 5. Qué aplicamos al proyecto y por qué

Cada cambio está en `lluvia_valle_funcional.R` (y en su copia `lluvia_valle.R`) o en `regresion_concurrente.R`, con su tabla, su gráfica y, cuando aplica, una prueba automática. Los resultados están en la sección 8 de este documento.

| Cambio | Fuente en el material | Por qué sirve |
|---|---|---|
| CV(t) en el EDA | Semana 12, diapositiva 12 | Compara la dispersión de dos fuentes con niveles distintos |
| Prueba t funcional pareada (L², F, permutación con criterio del máximo) | Semana 14, diapositivas 15–21 y ejercicio 3 | Convierte la banda d̄(t) ± 2 EE en una decisión global y controla el error en todo el año |
| ANOVA funcional del error del satélite por fase ENSO, vertiente y altitud | Semana 14, diapositivas 25–40 | Prueba con datos lo que antes solo se describía: si el sesgo del satélite depende del ENSO y del relieve |
| Vista acumulada | Tesis y semanas 12–13 | Muestra los resultados en la escala con que él trabaja, y justifica con datos por qué el modelo usa mm/día |
| Bandas de los efectos por bootstrap de estaciones | Semana 15 y tesis figs. 7.1 y 9.2 | Él siempre presenta β con banda; la nuestra respeta que los años de una estación no son independientes |
| Prueba F por permutación del aporte del satélite | Semana 15, diapositiva 15 (`Fperm.fd`) | Responde si el satélite aporta más allá de la altitud y el ONI |
| Residuos estudentizados, residuos frente a predichos y correlación de residuos | Tesis figs. 9.3 y 9.4; semana 15 | Diagnóstico del modelo como en su tesis |
| Traza-variograma de los residuos | Tesis cap. 8 | Mide si queda correlación espacial sin modelar, que fue su aporte principal |
| 20 curvas observadas frente a estimadas y mapa del error | Tesis figs. 10.3 a 10.6 | Así evaluó su modelo |
| Concurrente con bandas por bootstrap de estaciones y el panel de clase | Semana 15 y tesis cap. 7 | Es su modelo; las bandas por estación corrigen la independencia que asumen las sándwich |

## 6. Qué no aplicamos y por qué

| De su material | Por qué no | Qué hacemos en su lugar |
|---|---|---|
| Imputar faltantes con kriging | La clase y nuestras reglas no imputan: cada curva se ajusta con sus días observados, y solo se aceptan años con ≥ 90 % de días y ningún hueco > 15 días (R5) | El experimento de huecos (sección 7 del archivo de datos) muestra cuánto cambia la curva con huecos |
| Curvas acumuladas con `smooth.monotone` | El enunciado pide modelar sobre √P o log(1+P) y fechar los dos picos, que son rasgos de la intensidad. Además, la acumulada está dominada por el total que se va sumando (memoria larga) y oculta las fechas de las temporadas, que son una métrica de la rúbrica. La derivada de la acumulada es nuestra curva: no se pierde información | Se muestra la vista acumulada en el EDA (4.6) con la justificación en datos |
| Mínimos cuadrados generalizados con covarianza espacial | No se vio en clase y el enunciado pide la regresión como referencia de la red | Se mide la correlación espacial con el traza-variograma y se respeta con IC por bootstrap de estaciones |
| Altura² | En su tesis no aportó | La altitud entra lineal, como pide el enunciado |
| Base de Fourier para un año de lluvia | Él mismo lo advierte en el taller de repaso (ejercicio 3): "el valor final de un día no coincide necesariamente con el inicial del siguiente". En nuestros datos, enero y diciembre del mismo año solo se correlacionan 0,54 | B-splines, como en su tesis y en el enunciado |
| K por año | Nosotros elegimos K una vez con todo el entrenamiento, para que todas las curvas estén en la misma base (requisito de la comparación y de la red; Hall y Van Keilegom, citado en la semana 14) | — |

## 7. Lo que dice en los videos de clase

Las tres grabaciones se transcribieron automáticamente (66, 98 y 82 minutos). La calidad del audio es baja: las frases se leen a medias y aquí se resumen ideas, no citas textuales. Las clases son:
- construcción de datos funcionales y elección de K (semana 11);
- una práctica en R con SO₂ y PM2.5;
- el EDA funcional con SO₂ y con la lluvia `X_sat` / `Y_est` (semanas 12–13).

**Lo que dice del proyecto:**
- **El ejemplo de la lluvia es el proyecto en pequeño.** Presenta `X_sat` y `Y_est` como el ejemplo "asociado al proyecto final": precipitación acumulada por pentadas, el satélite como variable de entrada y la estación como salida.
- **La primera parte del proyecto** es "toda la construcción de los datos funcionales, toda la descripción funcional y toda la diferencia funcional". Con esos datos "se va a la arquitectura", es decir, a la red.
- **Los datos vienen con problemas:** valores negativos y faltantes. Hay que depurarlos antes de construir las curvas.
- **Advierte sobre la IA generativa:** pide "mucha atención" con su uso en el proyecto final; él mismo está haciendo pruebas con Claude. Esto refuerza dos cosas de la rúbrica: la redacción propia y el anexo con los prompts.
- **Lo que pesa es entender el método:** "si usted entiende el método, lo puede adaptar a cualquier aplicación". En la sustentación, cada integrante debe poder explicar cualquier fase.

**Lo que repite sobre el método:**
- **El supuesto de varianza constante del error es "muy fuerte" y en la realidad no se cumple:** hay que estimar la varianza. En el proyecto esto ya está hecho: la banda con σ²(t) y la covarianza residual de la sección 2.5.
- **Trabajar con curvas es más eficiente** que con escalares: por ejemplo, 686 curvas de píxeles en vez de 686 × 365 números, en memoria y en tiempo.
- **Lee cada figura en el dominio:**
  - en el SO₂, el pico de la mañana por los carros;
  - en la lluvia, la correlación satélite–estación baja en las pentadas secas de enero y febrero y sube a partir de la pentada 30.
- **De la correlación cruzada a la regresión.** Lo dice así: la correlación cruzada "me permite ver en qué instante hay menos correlación", y por eso luego "tenía que modelar eso con un modelo de regresión funcional: esa fue mi tesis". Es el mismo camino del proyecto: comparación (sección 4), luego regresión (sección 5).
- **La mediana funcional resiste las curvas atípicas y la media no.** Si una curva atípica entra en la media, "me daña la estimación". Por eso conviene mostrar la media, la mediana y la media recortada (`F3_centralidad`).

## 8. Resultados de lo aplicado

Corrida del 18-nov-2026: `lluvia_valle_datos.R` pasa 28 de 28 pruebas, `lluvia_valle_funcional.R` 16 de 16 y `regresion_concurrente.R` 10 de 10. El detalle está en D-57 y D-59 de `DECISIONES.md`.

| Pregunta | Prueba, como en clase | Resultado | Qué significa |
|---|---|---|---|
| ¿El satélite y la estación tienen la misma media? | t pareada: L², F y permutación con criterio del máximo; 176 estaciones | F: p = 0,003. Máximo: p = 0,0002. 129 días significativos | No. La diferencia se concentra en febrero–marzo (la estación llueve más) y en los picos de abril–junio y noviembre (el satélite exagera) |
| ¿El error del satélite depende del ENSO? | ANOVA funcional con los 23 años | F: p = 0,0003. GPF: p = 0,005. Fmax: p = 0,22 | Sí, con un efecto pequeño pero repartido en todo el año. Respalda dar el ONI al modelo |
| ¿Y de la vertiente y la altitud? | ANOVA funcional con las estaciones | p ≤ 0,005 | Sí: el satélite no ve la orografía, como concluye la tesis |
| ¿Qué muestra la vista acumulada? | d(t) acumulada ± 2 EE | Al 31-dic: −19 ± 71 mm | En el total anual los errores se compensan. Por eso el modelo trabaja en mm/día |
| ¿Aporta el satélite más allá de altitud y ONI? | F por permutación, como `Fperm.fd` | p = 0,003 en los 365 días | Sí |
| ¿Aporta la altitud? | Bandas por bootstrap de estaciones | La banda contiene el 0 todo el año | No, una vez se conoce el píxel. Igual que en la tesis |
| ¿Quedó correlación espacial? | Traza-variograma de los residuos | Rango práctico de unos 96 km, pepita del 17 % | Sí. Se respeta con IC por estación; mínimos cuadrados generalizados queda como trabajo futuro |
| ¿Eran correctas las bandas sándwich del concurrente? | Bootstrap de estaciones | Son 4 veces más angostas de lo que deberían | No. Se reportan las de bootstrap |

## 9. Su forma de graficar y qué metodologías suyas adoptamos (revisión del 18-nov-2026)

### 9.1 Cómo grafica (tesis, capítulos 5, 6 y 9)
- **Un panel por año.** Casi todo lo muestra así: los datos crudos (fig. 5.8), las curvas suavizadas con su promedio y desviación (figs. 6.2 y 6.3), las regiones (fig. 5.12) y hasta el coeficiente del satélite (fig. 9.2). Para él, el año es la unidad natural: la lluvia de un año depende del ENSO, y el panel por año lo deja ver sin una sola cifra.
- **Mismo eje en todos los paneles**, para comparar un año con otro de un vistazo.
- **Blanco y negro, con R base:**
  - curvas individuales en gris claro;
  - promedio en negro grueso;
  - desviación estándar en línea punteada;
  - bandas de confianza como sombra gris;
  - línea punteada horizontal en 0 en los coeficientes.
- **El año, en negrita, como título del panel**, y ejes en el lenguaje del problema: "Tiempo (día pentadal)" y "Precipitación acumulada (mm)".
- **Trabaja en lluvia acumulada.** Así las curvas siempre suben y se ven el total y las temporadas: donde la curva se aplana, es tiempo seco.

### 9.2 Lo que adoptamos de su estilo
- Las cuatro figuras por año del catálogo (sección 3), en su mismo estilo. A cada una le agregamos la fase ENSO en el título, porque es lo que explica las diferencias entre años.
- En la figura acumulada, el texto "satélite = x % de la estación" de cada año: es la "frase con números" de cada figura (regla 2.1).

### 9.3 ¿Usar sus metodologías para todas las decisiones?

Sí, donde el problema es el mismo. No, donde el enunciado del curso nos pide otra cosa: su tesis interpola en el espacio dentro de un mismo año, y nuestro proyecto debe **predecir años nuevos** (prueba 2001–2003) a partir de 1981–1995. Decisión por decisión:

| Su decisión en la tesis | La nuestra | Por qué |
|---|---|---|
| Graficar por año | **Adoptada** (cuatro figuras nuevas) | Es su forma de presentar y deja ver el ENSO |
| Elegir K por VCG | **Adoptada**, con una K para toda la muestra en lugar de una por año | Todas las curvas deben estar en la misma base para compararlas y para la red (sección 6) |
| Precipitación acumulada y `smooth.monotone` | **Como vista** (`F4_acumulada`, `12_acumulada_por_anio`); el modelo se queda en √P | El enunciado pide √P o log(1+P) y la fecha de los dos picos (sección 6) |
| Imputar faltantes con kriging | No | Las reglas R5 y R7 y el experimento de huecos: no hace falta imputar |
| Modelo concurrente | **Adoptado** (`regresion_concurrente.R`) | Es su modelo; queda como contraste del de función a función |
| Ajustar el modelo año por año | **Como descriptivo** (`C7_b1_por_anio`) | Un ajuste por año no puede predecir un año nuevo: necesita las estaciones de ese año. Sirve para ver cómo cambia el aporte del satélite con el ENSO, y confirma su conclusión (b₁ mayor en La Niña, r = −0,44) |
| Altura y altura² | Altura sí, altura² no | En su tesis ninguna de las dos aportó; en la nuestra la banda de la altitud contiene el 0 (D-57). La altitud se deja porque el enunciado la pide |
| Traza-variograma de los residuos | **Adoptado** | Mide la correlación espacial que queda |
| Mínimos cuadrados generalizados con covarianza espacial | No (trabajo futuro) | No se vio en clase; las bandas por bootstrap de estaciones respetan la dependencia |
| Validación cruzada dejando una estación fuera | No; se usa la partición por años | El enunciado fija entrenamiento, validación y prueba por años |
| Mann–Whitney para comparar errores | Wilcoxon **pareado** por estación | Es la misma familia de pruebas, en su versión pareada: los dos modelos se evalúan en las mismas estaciones |
| R²(t), 20 curvas y mapa del error | **Adoptados** (`F5_prueba`, `F5_curvas_20_estaciones`, `F5_mapa_error`) | Así evalúa él |
| Regiones valle y montaña | Zonas andina y pacífica, más franjas de altitud | Es el mismo criterio físico, aplicado a nuestro territorio |

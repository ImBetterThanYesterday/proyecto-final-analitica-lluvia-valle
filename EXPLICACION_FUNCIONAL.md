# Qué sale de la parte funcional y cómo se interpreta

Este documento sigue donde termina `EXPLICACION_DATOS.md`. Aquella explica el EDA de los datos diarios (`lluvia_valle_datos.R`). Esta recorre lo que viene después:

- la **reconstrucción**: pasar de 365 datos diarios a una curva;
- el **EDA funcional**;
- la **comparación** entre estaciones y satélite;
- la **regresión funcional**: la del equipo, de función a función, y la concurrente de María.

Para cada paso explica qué muestra cada tabla y gráfica, qué significa y qué decisión sale de ahí. Sirve para estudiar y para defender esta parte en la sustentación.

- **Dónde está el código:**
  - las secciones 1 a 6 de `lluvia_valle_funcional.R`;
  - el archivo unido `lluvia_valle.R` tiene las mismas secciones, idénticas salvo los encabezados de rúbrica, desde su línea 1535;
  - la regresión concurrente está en `regresion_concurrente.R`.
- **Las cifras** son de la corrida del 17-nov-2026:
  - `lluvia_valle_funcional.R` pasa sus 13 pruebas y `regresion_concurrente.R` las 8 suyas;
  - las tablas salen iguales a las que están en el repositorio.
- **Salidas:** gráficas en `salidas_R_funcional/graficos/` y tablas en `salidas_R_funcional/tablas/`, con prefijos F1 a F6. La concurrente deja las suyas en `salidas_R_concurrente/`, con prefijos C1 a C5.
- **Decisiones:** las D-xx están en `DECISIONES.md`.
- **Conceptos:** la guía de estudio (curva, base, matriz W, FPCA y demás) explica cada concepto con ejemplos.

> **Recordatorio.** "Cauca (interandina)" y "Pacífico" son las dos vertientes **dentro del Valle del Cauca**, según la cuenca de la CVC. No son departamentos.

## Resumen en una tabla

| Paso | Rúbrica | Qué se hace | Resultado clave | Sección del R |
|---|---|---|---|---|
| Emparejar | 3 | Cada estación toma el píxel del satélite que la contiene | 176 estaciones en 148 píxeles; 3 299 pares estación–año | 1 |
| Reconstruir | 1 (10 %) | B-splines cúbicas y mínimos cuadrados, sin librerías | K = 8 en estaciones y K = 11 en satélite, elegidos con VC y VCG | 2 |
| EDA funcional | 2 (15 %) | Media, varianza, covarianza, derivadas, profundidad y boxplot | Régimen bimodal; el 74 % de la varianza es diferencia entre sitios | 3 |
| Comparar fuentes | 3 (15 %) | d(t) con banda, correlación cruzada, por vertiente, altitud y ENSO | El satélite arranca tarde la primera temporada, exagera los picos y amortigua el ENSO | 4 |
| Regresión de función a función | 4.2 y 4.4 | β(s,t) con FPCA de entrenamiento y penalización | En prueba, el error L² relativo baja de 0,51 a 0,44 y el sesgo anual de +26 % a +7 % | 5 |
| Regresión concurrente | (contraste) | y(t) = b0(t) + b1(t)·x(t), como en el script de clase | El mismo error que la de función a función, pero peor en el sesgo y en el segundo pico | `regresion_concurrente.R` |

---

## 1. Datos y emparejamiento estación–píxel

**De dónde vienen los datos.** Todo sale de `salidas_R_datos/datos_preparados.rds`, que deja el archivo de datos. No se vuelve a leer el satélite.
- 3 299 curvas estación–año de 176 estaciones: 2 257 de entrenamiento (1981–1995), 672 de validación (1996–2000) y 370 de prueba (2001–2003). De ellas, 2 812 no tienen ningún día faltante.
- 15 663 curvas píxel–año: 681 píxeles por 23 años. Para el EDA se usan los 676 píxeles del Valle.

**Emparejamiento.** Cada estación toma la celda de 0,05° del satélite que contiene su coordenada, y una prueba comprueba que de verdad cae dentro.
- Las 176 estaciones caen en 148 píxeles.
- 23 píxeles tienen de 2 a 4 estaciones (51 estaciones en total), que comparten la misma curva satelital.
- El centro del píxel queda a 2,2 km de la estación en mediana, y nunca a más de 3,9 km.

**Vertiente de cada píxel.** Los píxeles no traen la cuenca de la CVC, así que hay que asignársela con alguna regla. Se compararon dos dejando una estación fuera cada vez:

| Regla | Acierta |
|---|---|
| Corte en la longitud −76,8° (la que se usaba) | 88,6 % |
| El mejor corte en longitud posible | 89,2 % |
| La vertiente de la estación más cercana | **97,7 %** |

Tabla: `F1_regla_vertiente_pixeles.csv`. El corte en longitud dejaba como "Cauca" a 20 de las 26 estaciones del Pacífico: las de las cuencas del Dagua, el Calima y el Garrapatas, sobre la cordillera Occidental. **Decisión:** se usa la estación más cercana. Así las dos fuentes tienen una sola definición de vertiente.

![El error del satélite sigue al ENSO](salidas_R_funcional/graficos/F1_cociente_satelite_estacion_vs_ONI.png)

**Qué muestra.** Cada punto es un año. El eje x es el ONI, y el eje y el cociente satélite/estación de la lluvia anual.

**Cifras.**
- La correlación es r = 0,68.
- El cociente medio es 0,99 en años Niña, 1,16 en neutros y 1,29 en años Niño.
- Los cinco cocientes más bajos son 1988, 1984, 2000, 1999 y 1996: cuatro de ellos son años Niña.

**Interpretación.** 1985, 1988 y 2000 no son archivos dañados del satélite. Son un error sistemático ligado al ENSO: en La Niña el satélite se queda corto y en El Niño se pasa. 1984 y 1999 se apartan igual y nunca se marcaron.

**Decisión.** No se excluye ningún año (D-47, lo pidió el profesor). El modelo recibe el ONI para que aprenda a corregir ese error.

---

## 2. Reconstrucción: de datos diarios a curvas (rúbrica 1)

### 2.1 La idea

![Datos discretos de una estación y su píxel](salidas_R_funcional/graficos/F2_datos_discretos.png)

**Qué muestra.**
- **A la izquierda**, la lluvia diaria de una estación en un año: muchos días en 0 y aguaceros sueltos.
- **A la derecha**, su píxel del satélite el mismo año: 72 péntadas (promedios de 5 días), cada una puesta en su día central.

Las dos fuentes tienen resoluciones distintas. La reconstrucción las pone a las dos como **curvas en el mismo dominio**, del día 1 al 365. Así se pueden sumar, restar, integrar y comparar.

**Cómo se hace.** La curva de cada año se escribe como una suma de K funciones base fijas, cada una multiplicada por un coeficiente:

x(t) = c₁ φ₁(t) + … + c_K φ_K(t)

- Las bases φ_k son **B-splines cúbicas**: polinomios de grado 3 por tramos, pegados sin saltos ni "codos".
- Toda la curva de un año queda guardada en K números, los coeficientes c.
- Las unidades son mm/día.

### 2.2 La base, programada sin librerías

![Las bases B-spline](salidas_R_funcional/graficos/F2_bases.png)

**Qué muestra.** Las 8 funciones base de las estaciones (izquierda) y las 11 del satélite (derecha). Cada una es una "colina" que sube desde 0, llega a un pico y vuelve a 0. Las líneas verticales son los nodos, donde se pegan los tramos. Con más bases las colinas son más estrechas y la curva puede seguir detalles más finos.

**Cómo se calcula.** Con la recursión de Cox–de Boor (funciones `nodos` y `base_en_nodos`): se empieza con funciones escalón y cada orden se arma promediando dos funciones del orden anterior.

**Verificaciones** (sección 2 y pruebas de la sección 6):
- La base propia coincide con `splines::splineDesign` con una diferencia de 3 × 10⁻¹⁶. La librería solo se usa para comprobar.
- En cada día las bases suman 1 y ninguna es negativa. Se comprueba con K = 6, 8, 11, 20 y 40.

**Con K = 8** (estaciones), los nodos interiores caen el 15 de marzo, el 27 de mayo, el 7 de agosto y el 19 de octubre. El año queda partido en 5 tramos.

### 2.3 Los coeficientes por mínimos cuadrados

Se arma la matriz Φ (365 × K), donde cada fila es un día y cada columna una base. Los coeficientes son los que dejan la curva lo más cerca posible de los datos:

ĉ = (ΦᵀΦ)⁻¹ Φᵀ y   y   ŷ = S y, con S = Φ (ΦᵀΦ)⁻¹ Φᵀ

**Verificación de S** (`F2_controles_S.csv`). La rúbrica pide comprobar tres propiedades:

| Fuente | K | Simetría, máx. ∣S − Sᵀ∣ | Idempotencia, máx. ∣SS − S∣ | ∣tr(S) − K∣ |
|---|---|---|---|---|
| Estaciones | 8 | 1 × 10⁻¹⁶ | 9 × 10⁻¹⁶ | 1 × 10⁻¹⁴ |
| Satélite | 11 | 8 × 10⁻¹⁷ | 2 × 10⁻¹⁶ | 2 × 10⁻¹⁵ |

Todas son cero hasta el error de redondeo del computador. También se cumplen en las curvas con días faltantes, que tienen cada una su propia S.

**Días faltantes.** No se imputan. Cada curva se ajusta solo con sus días observados: se quitan esas filas de Φ y de y. El archivo de datos ya garantiza que ningún año aceptado tiene un hueco de más de 15 días (regla R5).

### 2.4 Cuántas bases: elegir K

![Selección de K](salidas_R_funcional/graficos/F2_seleccion_K.png)

**Qué muestra.** Los tres criterios de clase en función de K, arriba para estaciones y abajo para satélite. Se calculan solo con las curvas completas de entrenamiento: 1 882 de estación y 10 215 de satélite.
- **ECM:** el error de ajuste. Siempre baja al subir K, así que no sirve para elegir.
- **VC:** el error al dejar un día fuera. Se calcula con el atajo de la diagonal de S, y se comprobó que da lo mismo que hacer los 365 ajustes.
- **VCG:** la validación cruzada generalizada.

**Cifras** (`F2_K_elegido.csv`).
- VC y VCG coinciden: **K = 8 en estaciones y K = 11 en satélite**.
- En el satélite el mínimo es nítido.
- En estaciones la VC es casi plana: cualquier K entre 5 y 17 queda a menos de 1 % del mínimo. La lluvia diaria es tan ruidosa que más bases solo persiguen aguaceros.

**Sensibilidad** (`F2_sensibilidad_K_estaciones.csv`). Como la VC es plana, se probó qué cambia con más bases:

| K | VC sobre el mínimo | Pico 1 (mediana) | Pico 2 (mediana) | Puntos negativos |
|---|---|---|---|---|
| **8** | 0 % | 18 abr | 25 oct | 2,0 % |
| 12 | 0,4 % | 27 abr | 28 oct | 2,1 % |
| 16 | 0,9 % | 23 abr | 29 oct | 3,3 % |
| 24 | 2,1 % | 18 abr | 24 oct | 5,4 % |

- La fecha de los picos se mueve a lo sumo 9 días.
- Eso es mucho menos que la dispersión entre curvas: el rango intercuartílico de la fecha del pico es de 29 a 45 días.
- Lo que sí sube con K son los valores negativos.

**P-20 (D-51).** Se probó también elegir K comparando la fecha del pico con una media móvil de los datos. El K que gana depende del ancho de la media móvil:
- con 15 o 31 días gana K = 20;
- con 61 días gana K = 10.

Las referencias difieren entre sí tanto como lo que se gana subiendo K, así que la fecha del pico no sirve para elegir K.

**Decisión.** K = 8 y K = 11, el mínimo de la VCG.

### 2.5 Residuos y bandas de confianza

![Bandas punto a punto](salidas_R_funcional/graficos/F2_bandas_punto_a_punto.png)

**Qué muestra.** La curva de una estación con su banda de ±2 errores estándar, su píxel con la suya, y una curva con 15 días faltantes (la marca roja).

**Cifras** (consola y `F2_varianza_error_estaciones.csv`).
- **El error no es igual todo el año.** Su varianza va de 26 (mm/día)² el 11 de agosto a 264 el 17 de octubre, la temporada húmeda. Por eso se dibujan dos bandas: con σ² constante y con σ²(t).
- **Los residuos de días seguidos casi no se correlacionan** (0,04 en general y 0,08 a un día; gráfica `F2_covarianza_residual.png`). Eso justifica la banda con σ²I. En el satélite la correlación es mayor (0,14), porque cada péntada ya promedia cinco días.
- Un hueco de 15 días ensancha la banda solo un 15 % en esos días.

**Interpretación.** Las bandas son **punto a punto**: cada día tiene 95 % de cobertura por separado, no la curva completa.

### 2.6 Valores negativos y efecto de borde

**Cifras** (`F2_valores_negativos.csv`).

| Fuente | Puntos negativos | Curvas con algún negativo | Mínimo |
|---|---|---|---|
| Estaciones | 1,9 % | 46,9 % | −12,0 mm/día |
| Satélite | 0 % | 1 % | −0,5 mm/día |

**Interpretación.** Los negativos aparecen en dos lugares:
- en la temporada seca de julio–agosto, donde la curva baja cerca de 0;
- en los bordes, el 1 de enero y el 31 de diciembre. Ahí la curva tiene datos de un solo lado y depende de un solo coeficiente: es el **efecto de borde** de una base que no es periódica.

**Decisión.** Para modelar, la curva se recorta en 0 (D-23). Se documenta como costo de la base.

### 2.7 ¿Por qué B-splines y no Fourier?

Fourier arma la curva con senos y cosenos, y obliga a que el 31 de diciembre empalme con el 1 de enero.

**Razones del código** (sección 2 y D-47):
- el enunciado exige B-splines para X(t) y para w_k(t) en el codificador de la red;
- el ejemplo de clase con esta misma lluvia usa B-splines cúbicas;
- una curva de un solo año no es periódica.

**Comprobación con los datos.** Se hizo el 16-nov-2026 con un script aparte y **todavía no está en el R**. Se ajustaron las 1 882 curvas completas de entrenamiento con las dos bases:

| | B-spline, K = 8 | Fourier, K = 7 |
|---|---|---|
| VC en todo el año | 105,2 | 105,0 |
| VC en los bordes, en el 25 % de curvas con enero y diciembre muy distintos | 161,2 | 165,7 |
| Curvas que bajan de 0 en los bordes | 34 % | 3 % |

- **Ajustan igual.** Fourier incluso tiene menos efecto de borde.
- **Pero un año no es periódico.** La lluvia del 1 al 15 de enero y la del 17 al 31 de diciembre del mismo año solo se correlacionan 0,54. En la climatología, el promedio de muchos años, la correlación sube a 0,97: ese es el caso en que Fourier es natural, como la temperatura de Canadá en clase.

**Conclusión.** La B-spline se queda porque el enunciado la exige, porque un año no es periódico y porque es local: un aguacero en abril casi no mueve la curva en octubre. El costo es el efecto de borde, que se maneja recortando en 0. Esta comprobación debería pasar al R para que quede respaldada en el código.

### 2.8 La matriz de Gram W

W = ∫ φ(t) φ(t)ᵀ dt permite calcular integrales con los coeficientes. Por ejemplo, ∫ x(t) y(t) dt = c_xᵀ W c_y. La red la necesita para z = A W c.
- Se calcula por cuadratura de Gauss–Legendre: 4 puntos en cada tramo entre nodos, lo que es exacto para estos polinomios.
- **Verificaciones:** 1ᵀW1 = 364, que es la longitud del dominio, y cᵀWc coincide con la regla del trapecio en una malla fina.

![Mismo dominio](salidas_R_funcional/graficos/F2_mismo_dominio.png)

**Qué muestra.** Una estación (K = 8) y su píxel (K = 11), ya como curvas, en el mismo dominio de 1 a 365. A partir de aquí las dos fuentes se pueden comparar día a día. Se ve también el efecto de borde: la curva de la estación baja a −2,6 mm/día el 31 de diciembre, y por eso se recorta en 0.

---

## 3. EDA funcional (rúbrica 2)

**Qué se usa.** Las 176 estaciones con sus 3 299 curvas de 1981–2003, sin submuestrear, y los 676 píxeles del Valle en los 23 años. La sección 8 del archivo de datos muestra que las estaciones incluidas representan a las 220 con datos. Esto evita el tope de la rúbrica por analizar solo un subconjunto.

### 3.1 Media, varianza y desviación estándar

![Media y dispersión](salidas_R_funcional/graficos/F3_media_y_dispersion.png)

**Qué muestra.** Las curvas en gris, la media en color y la banda de ±2 DE en rojo, para estaciones y para satélite. A la derecha, la desviación estándar funcional de las dos fuentes.

**Cómo se calcula.** La media es la base por el promedio de los coeficientes: x̄(t) = Φ(t) c̄. Se verifica contra `fda::mean.fd` con una diferencia de 3 × 10⁻¹⁵ (`F3_verificacion_media.csv`).

**Cifras** (`F3_extremos_de_la_media.csv`, `F3_dispersion_resumen.csv`).

| | Estaciones | Satélite |
|---|---|---|
| Primer máximo | 22 abr, 6,3 mm/día | 13 may, 10,7 mm/día |
| Segundo máximo | 28 oct, 7,2 mm/día | 5 nov, 12,3 mm/día |
| Mínimo de mitad de año | 1 ago, 2,2 mm/día | 20 jul, 4,6 mm/día |
| Coeficiente de variación medio | 0,96 | 0,61 |

**Interpretación.**
- **Las dos fuentes muestran el régimen bimodal:** dos temporadas de lluvia (abril–mayo y octubre–noviembre) y dos secas.
- **El satélite pone los picos de 1 a 3 semanas después y con casi el doble de magnitud.** Su media pesa por área, y el 43 % de los píxeles es de la vertiente Pacífica, donde llueven unos 4 000 mm/año. Por vertiente, el satélite supera a las estaciones en las dos: 1 735 frente a 1 551 mm/año en el Cauca, y 3 998 frente a 3 105 en el Pacífico.
- **La variación entre curvas es tan grande como la media** (CV de 0,96 en estaciones).
- **La DE de las estaciones** es alta en la segunda temporada y mínima el 27 de febrero. El máximo del 1 de enero (5,9) es efecto de borde, no un hecho del clima.

### 3.2 Covarianza y correlación

![Correlación de las anomalías](salidas_R_funcional/graficos/F3_correlacion_anomalias.png)

**Qué muestra.** La correlación entre la lluvia del día s y la del día t, después de restar a cada curva la media de su estación (anomalías). Las superficies sin restar están en `F3_covarianza_estaciones.png`, `F3_correlacion_estaciones.png` y sus equivalentes del satélite.

**Cifras** (`F3_covarianza_por_vertiente_y_anomalias.csv`).

| Grupo | 1.ª componente | Correlación abril–octubre | Correlación mínima |
|---|---|---|---|
| Estaciones, todas | 76 % | 0,74 | 0,39 |
| Estaciones del Cauca | 50 % | 0,54 | −0,14 |
| Estaciones, anomalías por estación | 29 % | 0,14 | −0,27 |
| Satélite, todos los píxeles | 88 % | 0,83 | 0,43 |

**Interpretación.**
- Con todas las estaciones juntas, todo se correlaciona positivamente y una sola componente explica el 76 %. Eso es sobre todo **el nivel de lluvia de cada sitio**: el 74 % de la varianza de las estaciones, y el 88 % de la del satélite, es diferencia *entre* sitios.
- Al quitar el nivel de cada estación queda la variación de un año a otro. Ahí la primera componente baja a 29 % y abril y octubre casi no se relacionan (0,14).
- **En palabras:** un sitio húmedo es húmedo todo el año, pero un año con una primera temporada lluviosa no anuncia una segunda temporada lluviosa.

### 3.3 Derivadas y temporadas de lluvia

![Derivada de la media](salidas_R_funcional/graficos/F3_derivada_media.png)

**Qué muestra.** La primera derivada de la media, es decir, la velocidad a la que cambia la lluvia.
- Donde es positiva, la lluvia sube; donde es negativa, baja.
- Donde cruza el 0 hay un máximo o un mínimo de la curva.

La derivada se calcula derivando la base: D x(t) = Σ c_k Dφ_k(t).

**Cifras** (`F3_temporadas_fechas.csv`, curvas medias).

| Momento | Estaciones del Cauca | Satélite en el Cauca |
|---|---|---|
| Inicio de la 1.ª temporada | 28 feb | 18 mar |
| Pico de la 1.ª temporada | 14 abr | 2 may |
| Fin de la 1.ª temporada | 14 jun | 11 jun |
| Inicio de la 2.ª temporada | 15 sep | 25 sep |
| Pico de la 2.ª temporada | 1 nov | 3 nov |

**Interpretación.**
- **El satélite fecha el inicio y el pico de la primera temporada 2 a 3 semanas tarde.**
- La gráfica `F3_temporadas_por_curva.png` hace lo mismo curva por curva. En el 62 % de las curvas de estación la lluvia sigue bajando hasta el borde del año, así que el fin de la segunda temporada no se puede fechar en ellas.

### 3.4 Profundidad, boxplot funcional y curvas atípicas

**Profundidad.** Mide qué tan "central" es una curva dentro del grupo. Se usan las dos de clase:
- **Fraiman–Muniz (FM)**, que promedia la profundidad de cada día;
- **MBD**, que cuenta cuántas veces la curva queda entre otras dos.

Las dos ordenan casi igual (Spearman 0,99 y 1,00). La curva más profunda, la mediana funcional, es **Betania 2002**. Gráfica: `F3_centralidad.png`.

![Boxplot funcional](salidas_R_funcional/graficos/F3_boxplot_funcional.png)

**Qué muestra.** El boxplot funcional:
- **región central**, en morado: el 50 % de las curvas más profundas;
- **bigotes**, en azul: 1,5 veces esa región;
- **curvas atípicas**, en rojo: las que se salen de los bigotes.

Arriba están las estaciones (todas, Cauca y Pacífico) y el satélite; abajo, las curvas medias de cada año. Las atípicas coinciden con las de `fda::fbplot`, y hay una prueba automática que lo verifica.

**Cifras** (`F3_atipicas_*.csv`, `F3_profundidad_anual.csv`).
- **Estaciones todas juntas:** 63 curvas atípicas. 62 son de cinco estaciones de la costa Pacífica (Bahía Málaga, Yurumanguí, San Francisco, Triana y El Tabano), con 7 000 a 12 600 mm/año.
- **Por vertiente:**
  - en el Cauca solo quedan 2 atípicas, Tenjo 1982 (4 482 mm/año, candidata a revisión) y Las Brisas 1981;
  - en el Pacífico quedan 44.
- **Satélite:** 330 atípicas. Todas son del Pacífico, y 229 tocan el tope de 200 mm por péntada.
- **Años:** ninguna curva media anual es atípica. Los años menos profundos de las estaciones son 1988, 1984, 1983 y 1999, tres de ellos Niña.

**Interpretación.** Las atípicas del Pacífico no son errores: son otro clima. Por eso el boxplot se separa por vertiente.

### 3.5 Contexto climático: vertiente, altitud y ENSO

![Contexto](salidas_R_funcional/graficos/F3_contexto.png)

**Qué muestra.** La media por vertiente, la de las estaciones del Cauca por franja de altitud y la de esas mismas estaciones por fase ENSO.

**Cifras** (`F3_contexto_vertiente_altitud_ENSO.csv`, `F3_picos_por_fase_ENSO.csv`).

| Grupo (estaciones del Cauca) | Lluvia anual |
|---|---|
| Menos de 1 200 m (fondo del valle) | 1 333 mm |
| 1 200 a 2 000 m | **1 715 mm** |
| Más de 2 000 m | 1 607 mm |
| Años neutros | 1 544 mm |
| Años Niña | 1 718 mm (+11 %) |
| Años Niño | 1 387 mm (−10 %) |

**Interpretación.**
- **Orografía:** las laderas medias son lo más lluvioso del Cauca, más que el fondo del valle y más que lo alto. Por eso la altitud entra al modelo.
- **La Niña se nota en la segunda temporada** (7,1 frente a 6,3 mm/día) y en una sequía de mitad de año más húmeda (2,4 frente a 1,4 mm/día).
- **El satélite no muestra la señal ENSO:** da menos lluvia en años Niña (1 664 mm) que en neutros (1 767). Si se quitan 1985, 1988 y 2000, el satélite pasa a dar 1 838 mm en La Niña. Es el mismo error ligado al ENSO de la sección 1, y la razón para darle el ONI al modelo.

### 3.6 Lista de la rúbrica 2

| La rúbrica pide | Dónde está |
|---|---|
| Media calculada paso a paso y verificada contra la librería | 3.1: `F3_verificacion_media.csv` (3 × 10⁻¹⁵ frente a `fda::mean.fd`) |
| Varianza y desviación estándar funcional | 3.1: `F3_media_y_dispersion.png`, `F3_dispersion_resumen.csv` |
| Covarianza y correlación de cada fuente | 3.2: superficies `F3_covarianza_*` y `F3_correlacion_*` |
| Derivadas e inicio y fin de las temporadas | 3.3: `F3_derivada_media.png`, `F3_temporadas_*.csv` |
| Profundidad FM y MBD, boxplot y atípicas por estación y año | 3.4: `F3_centralidad.png`, `F3_boxplot_funcional.png`, `F3_atipicas_*.csv` |
| Interpretación climatológica (bimodal, orografía, ENSO) | 3.1 a 3.5 |
| Sin subconjuntos sin justificar | Se usan todas las estaciones y todos los años |

---

## 4. Comparación funcional por fuente (rúbrica 3)

**Pares.** 3 299 pares estación–año y píxel–año: 176 estaciones en 148 píxeles, todos los años.

**El error estándar se calcula entre estaciones.** Primero se promedia cada estación y luego se compara entre estaciones. Así el EE sale 3,7 veces más grande que si cada curva fuera independiente, porque los años de una misma estación repiten su sesgo. Usar el EE ingenuo daría bandas falsamente estrechas.

### 4.1 La diferencia media d(t)

![Diferencia estación menos satélite](salidas_R_funcional/graficos/F4_diferencia.png)

**Qué muestra.**
- **A la izquierda**, la media de las estaciones y la de sus píxeles.
- **A la derecha**, la diferencia d(t) = estación − satélite, con su banda de ±2 EE. Las marcas rojas abajo son los días en que la banda excluye el 0. La línea roja discontinua es la diferencia solo en 1985, 1988 y 2000.

**Cifras** (`F4_intervalos_diferencia.csv`).

| Tramo | Días | Quién da más |
|---|---|---|
| 8 feb – 28 mar | 49 | La estación: el satélite arranca tarde la primera temporada |
| 15 abr – 9 jun | 56 | El satélite: exagera el primer pico |
| 28 oct – 18 nov | 22 | El satélite: exagera el segundo pico |
| 1–12 ene y 24–31 dic | 20 | Efecto de borde, no se interpreta |

**Interpretación.**
- **En el promedio anual casi no hay diferencia:** el satélite da 22 mm/año más. Pero el signo cambia con la temporada, y la banda excluye el 0 en 147 de los 365 días.
- Es el mismo retraso que mostraron las derivadas.
- En 1985, 1988 y 2000 el satélite pierde 337 mm/año frente a la estación; en el resto de años da 83 mm/año de más.

### 4.2 Correlación cruzada y covarianzas

![Correlación cruzada](salidas_R_funcional/graficos/F4_correlacion_cruzada.png)

**Qué muestra.** La correlación entre la estación en el día t y su píxel en el día s. La diagonal (s = t) dice qué tanto sigue el satélite a la estación el mismo día. La diagonal sola está en `F4_diagonal_correlacion.png`.

**Cifras.**
- **La diagonal** vale 0,62 en promedio, con un máximo de 0,78 a mediados de julio. La correlación máxima está en el mismo día (desfase mediano de −1 día).
- **Con anomalías**, restando la climatología de cada estación, la diagonal baja a 0,24: máximo 0,47 en febrero y apenas 0,07 en noviembre.
- **Covarianzas** (`F4_covarianzas.png`, `F4_covarianzas_comparadas.csv`):
  - la forma es casi igual: la correlación abril–octubre es 0,74 en estaciones y 0,72 en satélite, y la primera componente explica 76 % y 77 %;
  - la varianza del satélite es **una cuarta parte** de la de las estaciones (cociente mediano 0,24 por día).

**Interpretación.**
- Casi toda la correlación viene de que el satélite reconoce qué sitios son húmedos y cuáles secos. Sigue mal la variación de un año a otro, sobre todo en la segunda temporada.
- Un píxel de unos 5,5 km y una péntada de 5 días promedian los aguaceros locales. Por eso el satélite varía menos.

**¿Sale esto de usar K distinto?** No. Con el satélite también en K = 8 (`F4_sensibilidad_mismo_K.csv`), la diferencia media (−20 frente a −22 mm/año), la diagonal (0,64 frente a 0,62) y el cociente de varianzas (0,25 frente a 0,24) casi no cambian. El veredicto día a día sí es algo sensible: coincide en 281 de 365 días. Por eso las fechas exactas se reportan con esa salvedad.

### 4.3 Por vertiente, altitud y fase ENSO

![Sesgo según la altitud](salidas_R_funcional/graficos/F4_sesgo_altitud.png)

**Cifras** (`F4_comparacion_por_grupo.csv`, `F4_sesgo_por_estacion.csv`).

| Grupo | Estación (mm/año) | Satélite (mm/año) | Cociente satélite/estación | Días con la banda sin el 0 |
|---|---|---|---|---|
| Cauca | 1 551 | 1 669 | 1,08 | 263 |
| Pacífico | 3 105 | 2 827 | 0,91 | 46 |
| Años Niña | 1 926 | 1 731 | 0,90 | 273 |
| Años Niño | 1 570 | 1 833 | 1,17 | 269 |

**Interpretación.**
- **En el Cauca el satélite sobreestima un 8 %** (+109 mm/año).
- **En el Pacífico el sesgo depende del sitio:**
  - subestima a la mitad las estaciones más lluviosas (El Tabano, San Francisco y Bahía Málaga, cocientes de 0,43 a 0,50);
  - duplica la lluvia de los sitios secos a sotavento (Loboguerrero 2,6; La Despensa 2,3; Puerto Nuevo 2,3; Moravito 1,9).
  - Es decir, el satélite suaviza hacia la media regional: sube lo seco y baja lo muy húmedo.
- **Lluvia anual:** entre estaciones, la del píxel explica la de la estación con r = 0,79, pero dentro del Cauca solo con 0,42 (`F4_anual_estacion_vs_pixel.png`).
- **Altitud:** en el Cauca el cociente satélite/estación baja un 15 % por cada 1 000 m (p = 0,001). El satélite exagera en el fondo del valle y se queda corto en las laderas altas, donde la lluvia orográfica no se ve a 5 km.
- **ENSO:** el satélite subestima un 10 % en La Niña y sobreestima un 17 % en El Niño. Amortigua la señal.

**Decisión.** Esto justifica darle al modelo la **altitud** y el **ONI** junto con la curva del satélite, como pide el enunciado. La gráfica `F4_diferencia_por_grupo.png` muestra d(t) por grupo.

---

## 5. Regresión funcional de función a función (rúbrica 4.2 y 4.4)

### 5.1 Qué modelo es y por qué

Se predice **la curva de la estación completa** a partir de **la curva del satélite completa**, más la altitud y el ONI:

y_i(t) = α(t) + ∫ β(s,t) x_i(s) ds + γ(t)·altitud_i + δ(t)·ONI_i + error

- **β(s,t)** dice cuánto pesa la lluvia del satélite en el día s para predecir la estación en el día t. Es el modelo integral de clase.
- **Se resuelve en tres pasos:**
  1. Se hace un **FPCA** de las curvas de estación, solo con entrenamiento: y(t) ≈ μ(t) + Σ_m s_m ψ_m(t).
  2. Cada puntaje s_m se predice con la integral del satélite, que con coeficientes es z = C_x W, más la altitud y el ONI. Es una regresión de cresta en la que λ castiga la curvatura de β.
  3. Se arma β(s,t) = Σ_m β_m(s) ψ_m(t).
- **Es la referencia lineal de la red.** La cabeza lineal de la red tiene que reproducir este modelo (rúbrica 4.2).

**La escala raíz.** La lluvia es muy asimétrica, así que el modelo trabaja con la raíz de la curva en mm/día. Al final se eleva al cuadrado, lo que garantiza que la reconstrucción no es negativa.
- Recortar en 0 antes de sacar la raíz solo cambia el 2 % de los días de las estaciones.
- La raíz se reproyecta en la misma base con un error de 0,09 (estaciones) y 0,03 (satélite), en raíz de mm/día.

**Primero se suaviza y luego se saca la raíz** (`F5_orden_raiz_y_suavizado.csv`):

| Camino | Acumulado anual / dato diario (mediana) |
|---|---|
| A. Suavizar en mm, recortar en 0 y sacar la raíz (el que se usa) | **1,000** |
| B. Sacar la raíz al dato diario, suavizar y elevar al cuadrado | 0,329 |

El camino B reduce la lluvia anual a un tercio, porque el cuadrado de la media de las raíces no es la media. Como el sesgo del acumulado anual es una métrica de la rúbrica, se usa el camino A.

### 5.2 FPCA de las estaciones (solo entrenamiento)

![FPCA](salidas_R_funcional/graficos/F5_fpca.png)

**Qué muestra.** A la izquierda, la media de entrenamiento en escala raíz, con las dos temporadas. A la derecha, las 6 funciones propias ψ_m(t).

**Cifras** (`F5_fpca_varianza.csv`).

| Componente | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| Varianza explicada | 62 % | 13 % | 7 % | 5 % | 4 % | 3 % |
| Acumulada | 62 % | 75 % | 82 % | 87 % | 91 % | **94,2 %** |

- **La primera componente** es casi plana: es el nivel de lluvia de cada sitio. Las demás mueven las temporadas.
- **Se calcula sobre W¹ᐟ² Σ W¹ᐟ²**, para que la varianza se mida como integral y no como suma de coeficientes. Las funciones propias son ortonormales con W, y hay una prueba que lo verifica.

**P-19: ¿usar los registros de 1930–1979?** (D-56, `F5_P19_conteo.csv`).
- Habría 1 250 curvas válidas de 140 estaciones, 963 de ellas de los años setenta. Sumarían un 55 % a las de entrenamiento.
- Con ellas las cuatro primeras funciones propias no cambian (coseno 1,000) y la varianza explicada con M = 6 baja de 94,2 % a 93,2 %.
- **Decisión:** no se usan. Traen variabilidad de otra época, no estructura nueva, y la fe de erratas fija 1981–2003.

### 5.3 Elegir M y λ con validación

![Selección de M y lambda](salidas_R_funcional/graficos/F5_seleccion_M_lambda.png)

**Qué muestra.** El error en los años de validación (1996–2000) para cada número de componentes M y cada penalización λ.

**Regla.** Se toma el menor M cuyo mejor error de validación queda a menos de 0,5 % del mínimo.

**Resultado.** **M = 6** (94,2 % de la varianza; el enunciado pide M entre 4 y 6) y **λ = 3,2 × 10⁸**. Con ese λ el modelo usa unos 10 grados de libertad efectivos de 14 parámetros por componente. La prueba (2001–2003) no se usa para elegir nada.

### 5.4 La superficie β(s,t) y la equivalencia

![Superficie beta](salidas_R_funcional/graficos/F5_beta_superficie.png)

**Qué muestra.** β(s,t) en 3D y en curvas de nivel:
- en rojo, donde más lluvia del satélite en el día s sube la predicción del día t;
- en azul, donde la baja.

**Interpretación.**
- **Los valores altos siguen la diagonal:** la lluvia del satélite de una época predice sobre todo esa misma época de la estación.
- **Fuera de la diagonal hay contrastes entre temporadas.**

![Efectos](salidas_R_funcional/graficos/F5_efectos.png)

**Qué muestra.**
- **∫β(s,t) ds**, a la izquierda: cuánto pasa a la estación si el satélite sube igual todo el año. Lejos de los bordes va de 0,62 en diciembre a 1,03 el 1 de junio. El modelo traslada el nivel del satélite casi completo a mitad de año y lo amortigua al final.
- **El efecto del ONI δ(t)** resta lluvia todo el año, y sobre todo a mitad de año (−0,24 por unidad de ONI el 2 de julio). El Niño profundiza la sequía de junio–julio.
- **La altitud γ(t)** pesa poco una vez que se conoce el píxel.

**Equivalencia.** La predicción con la integral de β por cuadratura y la predicción por puntajes FPCA difieren en 4 × 10⁻¹⁵ (el criterio es 10⁻⁴). Es la identidad β(s,t) = Σ_m β_m(s) ψ_m(t) que la cabeza lineal de la red debe reproducir.

### 5.5 Evaluación en prueba (2001–2003)

![Prueba](salidas_R_funcional/graficos/F5_prueba.png)

**Qué muestra.**
- **Arriba**, tres curvas de prueba con la estación, el satélite y la regresión. La Despensa 2001 es un caso difícil: un sitio seco que el satélite ve húmedo.
- **Abajo**, la varianza explicada por día, el error por curva y el sesgo anual, comparados con el satélite sin corregir y con un modelo que solo usa altitud y ONI.

**Cifras** (`F5_metricas_prueba.csv`: 370 curvas de 137 estaciones; IC del 95 % por bootstrap de estaciones con 2 000 remuestreos).

| Modelo | L² relativo | KGE | Sesgo anual | Error pico 1 | Error pico 2 |
|---|---|---|---|---|---|
| Satélite sin corregir | 0,51 [0,47; 0,56] | 0,38 | +26,0 % | 23,1 días | 18,0 días |
| **Regresión funcional** | **0,44 [0,41; 0,46]** | **0,47** | **+7,2 %** | **18,1 días** | **17,4 días** |
| Solo altitud y ONI | 0,51 [0,47; 0,54] | 0,35 | +19,3 % | 21,7 días | 19,4 días |

**Pareado por estación** (`F5_regresion_vs_satelite_pareado.csv`):
- en cinco métricas el IC de la diferencia excluye el 0 y el Wilcoxon da p < 0,001;
- en el pico de la segunda temporada la ganancia es de medio día y **no es significativa** (IC de −1,2 a 0,08 días, p = 0,08);
- la regresión mejora el error en el 65 % de las curvas.

**Interpretación.**
- La regresión corrige sobre todo el **nivel** (sesgo de +26 % a +7 %) y la **fecha del primer pico** (5 días mejor).
- Sin el satélite, con solo altitud y ONI, el error queda igual que el del satélite crudo. **El satélite sí aporta.**

### 5.6 Por fase ENSO y sensibilidad

**Por fase** (`F5_metricas_por_fase_ENSO.csv`). La prueba solo tiene años neutros (2001 y 2003) y un Niño (2002), ningún año Niña. Por eso se muestra también la validación, de forma orientativa:

| Años | Satélite: L² / sesgo | Regresión: L² / sesgo |
|---|---|---|
| 1997, Niño (validación) | 0,61 / +29 % | **0,46 / −6 %** |
| 1999–2000, Niña (validación) | 0,41 / −0,5 % | 0,40 / **−11 %** |
| 2002, Niño (prueba) | 0,50 / +28 % | 0,39 / −14 % |

**Interpretación.**
- El modelo gana más en años Niño.
- En La Niña casi empata con el satélite y el sesgo empeora: baja el nivel del satélite como en un año normal, y el término del ONI no alcanza a devolver lo que el satélite ya perdió.
- Es el **punto débil de la referencia lineal**, y lo que una cabeza no lineal (KAN o MLP) podría mejorar.

**Sensibilidad a 1985, 1988 y 2000** (`F5_sensibilidad_anios_marcados.csv`). Sin esos años en el entrenamiento, la prueba empeora:
- L² relativo de 0,434 a 0,439 y sesgo de 7,4 % a 8,9 %;
- en 79 de 137 estaciones (Wilcoxon p = 0,003).

Confirma con datos la indicación del profesor de no quitarlos (D-47).

**No negatividad.** La reconstrucción en mm nunca es negativa. Solo el 0,02 % de los valores se recorta en 0.

### 5.7 Contraste con los métodos de clase (D-53)

Se compararon tres formas de hacer la regresión vista en clase (`F5_contraste_metodos_clase.csv`):

| Método | L² relativo | KGE | Sesgo anual |
|---|---|---|---|
| Nuestra regresión (penalización de rugosidad) | 0,44 | **0,47** | +7,2 % |
| Taller de clase: regresión sobre componentes principales (l = 4) | 0,44 | 0,44 | +6,9 % |
| `fda.usc::fregre.basis.fr` | 0,44 | 0,44 | +7,6 % |

**Interpretación.** Los tres llegan al mismo error. La nuestra sube el KGE unos 0,03 (IC que excluyen el 0), porque la penalización conserva mejor la forma y la amplitud que truncar en 4 componentes. **No es un método aparte:** es el modelo integral de clase. La altitud y el ONI no cambian la prueba, pero en validación bajan el RMSE de 2,42 a 2,29 mm/día, y por eso se dejan.

### 5.8 Lo que pasa a la red

La sección 5.7 del R exporta a `salidas_R_funcional/para_python/`:
- los coeficientes de estaciones y satélite, en mm y en raíz;
- las matrices W de las dos bases y la de rugosidad R;
- la media y las funciones propias del FPCA;
- los pares con su partición, altitud y ONI;
- los coeficientes de esta regresión, y un archivo con la especificación (K, M, λ).

La red (`lluvia_valle_red.py`) usa el mismo codificador (z = A W c) y el mismo decodificador FPCA. Su cabeza lineal reproduce esta regresión al decimal.

---

## 6. Regresión funcional concurrente (`regresion_concurrente.R`, de María)

### 6.1 Qué es "concurrente"

Es el ejercicio 3 del script de regresión funcional de clase. **Cada día t de la estación se explica solo con el mismo día t del satélite**:

- **Modelo A, el de clase:** y_i(t) = b0(t) + b1(t)·x_i(t) + error.
- **Modelo B:** agrega la altitud y el ONI, γ(t)·altitud_i + δ(t)·ONI_i.

**La diferencia con la sección 5:**
- En la de función a función, todo el año del satélite puede predecir cada día: β(s,t) es una superficie.
- En la concurrente solo cuenta la diagonal s = t: b1(t) es una curva.

Es un modelo más simple y más fácil de leer, porque dice en qué época del año el satélite sirve más.

**Cómo se ajusta.** b0 y b1 se escriben en una base B-spline de K_b funciones, con penalización de curvatura λ. El número de bases y λ se eligen con validación:
- modelo A: K_b = 12 y λ = 0;
- modelo B: K_b = 8 y λ = 10⁸.

Usa las mismas curvas en escala raíz, la misma partición y las mismas métricas que la sección 5. Pasa sus 8 pruebas, entre ellas que no hay fuga de datos y que las ecuaciones normales dan lo mismo que mínimos cuadrados directos. Además se contrastó con `fda::fRegress`: la diferencia del error de ajuste es del 0,4 %.

### 6.2 Qué dice b1(t)

![Coeficientes del modelo A](salidas_R_concurrente/figuras/C2_coeficientes_modeloA.png)

**Qué muestra.** El intercepto b0(t) y el efecto del satélite b1(t), con bandas de ±2 EE.
- Si b1 = 1, la estación sube lo mismo que el satélite.
- Si b1 < 1, sube menos.

**Cifras** (`C3_R2_y_b1_por_mes.csv`).
- b1 vale unos 0,68 en marzo y 0,73 en noviembre: en las temporadas de lluvia la estación sube menos que el satélite.
- b1 vale unos 1,14 en julio, en la temporada seca.

**Interpretación.** El satélite exagera la lluvia de las temporadas húmedas, y el modelo la amortigua. Es lo mismo que mostró la comparación (sección 4.1).

![R² por día](salidas_R_concurrente/figuras/C3_R2_por_dia.png)

**Qué muestra.** R²(t), qué parte de la variación de la estación explica el modelo cada día de los años de prueba.

**Cifras.**
- Es mayor a mitad de año: 0,59 en junio y 0,56 en julio.
- Es bajo en enero y marzo (0,17).
- En noviembre el satélite solo explica más que el modelo (0,42 frente a 0,35).

### 6.3 Resultados en prueba y comparación con la de función a función

| Modelo | L² relativo | KGE | Sesgo anual | Error pico 1 | Error pico 2 |
|---|---|---|---|---|---|
| Satélite sin corregir | 0,51 | 0,38 | +26,0 % | 23,1 días | 18,0 días |
| Concurrente A (clase) | 0,44 | 0,47 | +10,6 % | 18,9 días | **23,0 días** |
| Concurrente B (+ altitud y ONI) | 0,45 | 0,46 | +10,3 % | 20,8 días | 19,4 días |
| De función a función (sección 5) | 0,44 | 0,47 | **+7,2 %** | **18,1 días** | **17,4 días** |

Tablas: `C4_metricas_prueba.csv`, `C4_concurrente_vs_satelite_pareado.csv`, `F5_metricas_prueba.csv`.

**Interpretación.**
- **En error L² y KGE la concurrente empata con la de función a función.**
- **La de función a función gana en el sesgo anual** (+7 % frente a +11 %) y en los picos.
- **La concurrente A empeora el segundo pico** en 5 días frente al satélite sin corregir, y la diferencia es significativa (IC de 4,1 a 5,9 días). Al mirar solo el mismo día, no puede corregir un retraso de la temporada.
- **Conclusión:** la concurrente sirve como contraste de clase y para leer b1(t) y R²(t). La de función a función sigue siendo la referencia.

**Ojo con las bandas de b1(t).** Se calculan como si cada curva fuera independiente. Los años de una misma estación repiten su sesgo, y en la comparación (sección 4) eso multiplicó el EE por 3,7. Así que estas bandas probablemente son demasiado estrechas. Habría que calcularlas agrupando por estación, igual que el bootstrap de las métricas.

---

## 7. Pruebas automáticas de la parte funcional

Si una prueba falla, el script se detiene. Tabla: `F6_pruebas.csv`, con 13 de 13 pruebas que pasan.

| Prueba | Qué protege |
|---|---|
| El emparejamiento es el mismo del archivo de datos | Que las dos partes usen los mismos pares |
| Ningún año queda fuera del modelo | La indicación del profesor (D-47) |
| La base propia coincide con `splines::splineDesign` | La construcción sin librerías (rúbrica 1) |
| S es simétrica e idempotente, y tr(S) = K | Las propiedades que pide la rúbrica |
| W integra el dominio: 1ᵀW1 = 364 | Las integrales de la red |
| La media paso a paso coincide con `fda::mean.fd` | La verificación que pide la rúbrica 2 |
| Las atípicas coinciden con `fda::fbplot` | El boxplot programado |
| Las funciones propias son ortonormales con W | El FPCA |
| Sin fuga: alterar validación y prueba no cambia K, FPCA, escalas ni B | Que la prueba no se use para ajustar |
| Equivalencia de β(s,t) con los puntajes FPCA (≤ 10⁻⁴) | La rúbrica 4.2 |
| La reconstrucción en mm no es negativa | La rúbrica 4.2 |
| Bootstrap con al menos 1 000 remuestreos | Los IC de la rúbrica 4.4 |
| El contraste de clase está bien armado | Que la comparación con los métodos de clase sea justa |

---

## 8. Lo que falta o hay que cuidar

1. **Fourier frente a B-spline:** la comprobación de la sección 2.7 está en un script aparte. Hay que pasarla al R y a `DECISIONES.md` para que quede respaldada en el código.
2. **Bandas de b1(t) en la concurrente:** hay que calcularlas agrupando por estación (sección 6.3).
3. **El segundo pico:** ninguna regresión lo mejora de forma significativa, y la concurrente lo empeora. Hay que decirlo así en el informe.
4. **La Niña:** la prueba no tiene años Niña. Lo que se diga de La Niña sale de la validación y es orientativo.
5. **¿Qué R se edita?** Hay cuatro archivos: `lluvia_valle_datos.R`, `lluvia_valle_funcional.R`, el unido `lluvia_valle.R` y `regresion_concurrente.R`. Hay que decidir cuál es el definitivo para que las copias no se desincronicen.

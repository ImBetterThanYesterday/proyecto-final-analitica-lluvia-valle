# Fundamentos de las reglas y los métodos

Este documento responde la pregunta del profesor en la revisión por Teams: **«¿de dónde sacan esas reglas?»** Para cada regla de limpieza y cada método dice:
- **qué hace**;
- **de dónde sale**: la norma, el libro o el artículo que la respalda;
- **qué muestran nuestros datos**;
- **qué tomamos de lo que dijo el profesor y en qué fuimos críticos.**

Lo que dijo el profesor está en `docs/transcripcion_revision_profesor.txt`. Las cifras salen de `salidas_R_datos/tablas/` (corrida del 20-nov-2026, conjunto principal).

---

## 1. Resumen

| Regla o método | Qué hace | Fuente principal | Antes |
|---|---|---|---|
| R0 | Valor negativo o mayor que el récord mundial en 24 h → faltante | Récord de la OMM: 1 825 mm, Foc-Foc, 1966 | Igual |
| R1 | Se quita el 29 de febrero (dominio de 365 días) | Convención de curvas anuales (Ramsay y Silverman 2005) | Igual |
| R2 | Extremo sin respaldo de las vecinas → faltante | Control de consistencia espacial (Durre et al. 2010) | Igual |
| R3 | Total del mes anotado en un día → mes faltante | Acumulados sin marcar (Viney y Bates 2004; OMM 2017, 4.4.1 c) | Igual |
| **R4** | **Racha de ≥30 días en 0 que al menos 2 de sus 3 vecinas más cercanas contradicen → faltante** | **Regla de periodos secos de GSDR-QC (Lewis et al. 2021); q por temporada como la ventana de Durre et al. (2010)** | Decidía con una sola vecina y q de todo el año (D-61) |
| **R5** | **Año válido: en cada mes faltan ≤10 días y no hay 5 o más seguidos** | **OMM 2017 (WMO-No. 1203), 4.4.1** | ≥90 % de días y hueco ≤15 (propia) |
| **R6** | **Estación: ≥10 años válidos en 1981–2003. Sensibilidad: ≥80 % (19 de 23)** | **Mínimo de la OMM para cualquier promedio (WMO-No. 1203, sec. 3); NOAA, normales "representativas". Sensibilidad: OMM 2017, 4.4.2** | ≥7 años válidos en entrenamiento (propia) |
| **R7** | **Los faltantes de los años válidos se rellenan por kriging ordinario con las vecinas** | **OMM 2017, 4.6 (interpolación espacial); Cressie 1993; Terán-Chaves et al. 2021 (IDEAM)** | No se imputaba: la curva cubría los huecos (D-13) |
| Lluvia en mm/día, nunca acumulada | Se modela la intensidad diaria | Ramsay y Silverman 2005, cap. 6; el profesor | Había dos gráficas de lluvia acumulada |
| Comparar modelos por estación con IC bootstrap | El IC de la diferencia decide | Davison y Hinkley 1997, sec. 3.8 | Además había Wilcoxon |
| ANOVA funcional por grupos independientes | ¿El error del satélite cambia con ENSO, zona o altitud? | Cuevas et al. 2004; Zhang 2013; Górecki y Smaga 2015 | Igual |
| Sin prueba t ni F pareada | El sesgo del satélite lo responde la regresión | Ramsay y Silverman 2005, cap. 13–14 | Había t y F pareadas |

---

## 2. R5 y R6: cuántos datos faltantes se aceptan

### Lo que dijo el profesor
Las reglas del 90 % y de los 15 días eran subjetivas («si es subjetivo, no está bien, porque eso le quita el rigor científico»). Además, la serie de Bahía Málaga tenía bloques vacíos de años enteros, y pidió explicarlos.

### La norma: la guía de la OMM para datos faltantes
La Organización Meteorológica Mundial publica cómo tratar los datos faltantes al calcular valores climáticos (*WMO Guidelines on the Calculation of Climate Normals*, WMO-No. 1203, 2017, sección 4.4). Dice tres cosas que usamos.

**1. Un mes de un promedio** (4.4.1 a) no se calcula si se cumple cualquiera de estas dos condiciones:
- faltan 11 días o más;
- faltan 5 días o más seguidos.

La guía explica por qué cuentan los días seguidos: los valores diarios están autocorrelacionados, así que un hueco seguido pesa más que la misma cantidad de días sueltos.

**2. Un mes de un total**, como la lluvia (4.4.1 c), exige el mes completo. Hay una excepción: el mes se acepta si los días **estimados** cumplen el mismo límite, es decir, menos de 11 en total y menos de 5 seguidos.

**3. Para una estación** (4.4.2), se exigen datos válidos en al menos el 80 % de los años del periodo. La guía de 2017 **quitó** el criterio de años faltantes seguidos que traía la guía de 2011: un estudio de la OMM (2007) encontró que no reduce la incertidumbre.

### Por qué aplica a nuestras curvas
Los días faltantes de un año válido se **estiman** con las vecinas (R7, sección 6 de este documento). Por eso la regla que corresponde es la de totales con días estimados (punto 2):
- **R5:** un año vale si en **cada uno de sus 12 meses** faltan como máximo 10 días y no hay 5 o más seguidos.
- **R6 (D-63):** una estación entra si tiene **al menos 10 años válidos** en 1981–2003. Un año en que la estación no operaba cuenta como no válido. El 80 % de la OMM (19 de 23) queda como sensibilidad; por qué, más abajo.

Un valor que el control de calidad (R0 a R4) volvió faltante cuenta como faltante: lo dice la nota 2 de la sección 4.4.

### Lo que muestran nuestros datos
Con el experimento de huecos de la sección 7 se mide cuánto mueve la curva cada regla. Se borran días a años completos copiando los huecos reales de otros años. El error es cuánto se mueve la curva, dividido por lo que varía la estación de un año a otro (`7_reglas_comparadas.csv`):

| Regla | Años de entrenamiento que conserva | Error típico (mediana) | Error del peor 10 % | Error del peor 1 % |
|---|---|---|---|---|
| **OMM 2017 (R5)** | 2 227 | 0,03 | 0,12 | **0,34** |
| OMM 1989 (más estricta: 5 faltantes y 3 seguidos) | 2 191 | 0,03 | 0,11 | 0,30 |
| Regla propia anterior (≥90 %, hueco ≤15) | 2 361 | 0,04 | 0,21 | 0,64 |
| Sin regla (≥40 % de días) | 2 635 | 0,10 | 1,88 | 191,7 |

**Lectura:** la regla de la OMM conserva el 94 % de los años que conservaba la nuestra y deja **la mitad del error** en el peor 1 % (0,34 frente a 0,64). Funciona mejor porque limita los huecos **mes a mes**: la regla anterior aceptaba un hueco de 15 días seguidos, y la OMM no acepta más de 4. La regla de 1989 casi no cambia el resultado.

**R6: ¿80 % de los años o al menos 10?** Al principio adoptamos el 80 % de la OMM (D-60). Al comparar los dos criterios (D-61), el profesor dijo que le hace más sentido el de 10 años, bien fundamentado, porque no quita tantas estaciones del Pacífico. Es la regla principal desde D-63.

**El fundamento del mínimo de 10 años:**
1. **Es el mínimo de la OMM para cualquier promedio climatológico.** El Reglamento Técnico de la OMM, citado en la guía de 2017 (sección 3), define los *period averages* como promedios «computed for any period of at least ten years». El 80 % (4.4.2) es para las **normales estándar** de 30 años, que son otro producto: el promedio oficial de una estación.
2. **NOAA lo aplica igual, y rellena como nosotros.** En sus normales 1991–2020: «according to the WMO guidelines, a minimum of 10 years is required for calculating any normal». Una estación con ≥10 años que no llega al 80 % recibe la marca "R" (*Representative*) si sus huecos se rellenan con estaciones vecinas. Es nuestro caso: R7 rellena con kriging.
3. **En Colombia, para el mismo tipo de análisis.** Urrea, Ochoa y Mesa (2019), *Water Resources Research* 55:4149–4162, estudian el **ciclo anual diario** de la lluvia con estaciones del IDEAM y CHIRPS, en un dominio de 365 días como el nuestro. Verificado en el texto (sección 2): de 2 716 series diarias del IDEAM descartan 120 por baja calidad y **890 por tener menos de 20 años**; quedan 1 706, con registros entre 1970 y 2015 y «some missing data in many time series». No piden un porcentaje de años completos, sino un registro de 20 años en 46 (≈43 %); en nuestros 23 años esa proporción da unos 10. La equivalencia es una lectura nuestra.
4. **Nuestra unidad es la curva estación-año, no la normal de la estación.** El modelo aprende de curvas. Un año que falta es una curva menos, no un hueco que haya que rellenar. Thornton et al. (2021), con la estación-año como unidad, filtran años y **no** piden un mínimo de años por estación.

**Lo que muestran los datos** (`8_R6_umbrales.csv` y `8_representatividad.csv`):

| | ≥10 años (principal) | ≥15 años | ≥19 años (OMM, sensibilidad) | Las 220 del Valle |
|---|---|---|---|---|
| Estaciones | **163** | 139 | 100 | 220 |
| Zona pacífica | **22 (13 %)** | 15 | 7 (7 %) | 19 % |
| Curvas estación-año | **3 046** | 2 761 | 2 113 | — |
| Altitud mediana | 1 350 m | | 1 398 m | 1 330 m |
| Lluvia anual mediana | 1 430 mm | | 1 419 mm | 1 515 mm |

- Con ≥10 años el conjunto **se parece más a todo el Valle**, sobre todo en la zona pacífica.
- Las 163 estaciones tienen años de entrenamiento; 146 tienen validación, 136 prueba y 133 los tres.
- **La conclusión del modelo no depende de R6** (prueba 2001–2003, L² relativo): satélite 0,51 → regresión 0,44 con ≥10 años; 0,52 → 0,42 con la OMM. Con ≥10 años el satélite aporta más: sin él, el error sube a 0,50.
- Con la nueva R4 (D-61), el conjunto de la OMM tiene 100 estaciones y no 98: la R4 anterior les quitaba un año a Aují y Puerto Frazadas.

Se corre la sensibilidad con `CONJUNTO=r6_omm` (ver el README).

### Dónde fuimos críticos
- **Con los bloques de años vacíos.** El profesor pidió que no haya bloques de años vacíos. Con ≥10 años, algunas estaciones los tienen; Bahía Málaga, por ejemplo, entra con 15 años válidos.
  - El modelo nunca reconstruye un año que falta. Cada curva es un año de una estación, como en el ejemplo del clima de Canadá de Ramsay y Silverman. Un año ausente es simplemente una curva menos.
  - La OMM 2017 tampoco exige años seguidos, porque no reduce la incertidumbre (WMO 2007).
  - Donde sí importa es en los resultados **por estación**: una estación con pocos años tiene un promedio menos firme. Por eso esos resultados se leen con su número de años.
- **Con la jerarquía de las normas.** El 80 % es la regla de la OMM para normales estándar, y por eso la usamos primero. Pero la misma OMM fija 10 años como mínimo de cualquier promedio, y nuestro uso no es una normal. La sensibilidad con el 80 % muestra que la conclusión no cambia.
- **Con el tope de 200 mm del satélite.** Con ≥10 años entran las tres estaciones de la costa donde el satélite toca el tope (Bahía Málaga, Yurumanguí y San Francisco). Aun así, el satélite está en el tope en **33 de 218 806** péntadas estación-satélite (0,015 %). En las 950 péntadas en que la estación pasa de 200 mm, el satélite está en el tope solo en el 2 %: el tope no es la causa de que el satélite subestime los aguaceros (D-47). Con el conjunto de la OMM no hay ninguna.

---

## 3. R4: rachas de ceros y el ENSO

### Lo que dijo el profesor
En un año El Niño es normal que no llueva por semanas («estábamos en superniño... hace mucho tiempo no llovía»). Pidió revisar las rachas según la temporada (seca o lluviosa) y la fase ENSO.

### Qué busca la regla
**Ceros falsos:** meses en que el pluviómetro anotó 0 todos los días pero sí llovió (nadie lo leyó, estaba dañado o los días en blanco se pasaron como 0). Un 0 parece un dato válido, así que ninguna otra regla lo detecta.

### Las fuentes y cómo las usamos (D-61)

| Fuente | Qué hace con rachas largas sin lluvia | ¿Usa la fase ENSO? |
|---|---|---|
| **GSDR-QC** (Lewis et al. 2021; código `nclwater/intense-qc`, regla R13) | Ventanas de 15 días en 0: las marca si **todas** sus vecinas (hasta 10 a ≤50 km) llovieron ≥3 días y la racha supera la racha seca más larga de la zona | No |
| **GHCN-Daily** (Durre et al. 2010) | No revisa rachas de ceros («missing values and zeros are skipped»). Compara con vecinas solo valores altos, con una ventana de 29 días por época del año | No |
| **Terán-Chaves et al. (2021)**, datos diarios del IDEAM | Marca los meses con más días secos que el p99,5 de la estación, sin vecinas | No |
| **INQC** (paquete que usa el IDEAM, según el resumen del buscador de su nota técnica 2023-004) | `drywetlong`: rachas más largas de lo esperable en la propia estación, sin vecinas | No |

**Ninguna fuente trata distinto los años Niño:** una sequía real se protege porque las vecinas también estuvieron secas.

**La regla de GSDR, tal cual, no funciona en el Valle.** Solo usa vecinas que coinciden en más del 90 % de los días en si llovió o no; en el Valle la coincidencia mediana es 0,73 y solo 6 de 1 897 pares pasan. Sin vecinas, la regla borra toda racha más larga que la máxima de la zona, que son justo las sequías de El Niño. Por eso tomamos **su principio** (varias vecinas deben contradecir la racha) y la ventana estacional de Durre:
1. Para cada racha se toman las **3 vecinas más cercanas** a ≤25 km con datos en ≥80 % de la racha.
2. Para cada vecina, q = P(llueve en la estación | llueve en la vecina), con los días de entrenamiento del **mes central de la racha y sus dos vecinos** (q de la temporada).
3. Si la vecina llovió k días, p = (1 − q)^k. La vecina **contradice** la racha si p < 0,05 / 917 (Bonferroni).
4. La racha pasa a faltante si la contradicen **al menos 2 de las 3**.

El equipo de CHIRPS trata como sospechoso el 0 de una estación cuando **el satélite** estimó lluvia. Nosotros usamos las vecinas y no el satélite, porque el satélite es justo lo que queremos corregir.

### Lo que muestran nuestros datos
**Versiones de la regla** (`5b_rachas_ceros_versiones.csv`):

| Versión | Rachas a faltante |
|---|---|
| Vecina más cercana, q de todo el año (R4 anterior) | 30 |
| Vecina más cercana, q de la temporada | 18 |
| **Al menos 2 de las 3 vecinas, q de la temporada (R4 vigente)** | **10 (469 días)** |
| Las 3 vecinas (como exige GSDR) | 4 |

Que el número baje de 30 a 10 muestra que la regla anterior era frágil: dependía de una sola vecina. En las 10 rachas marcadas, al menos 2 vecinas a ≤20 km llovieron entre 8 y 44 días; por ejemplo, Loboguerrero marcó 0 durante 67 días mientras sus tres vecinas llovieron 39, 44 y 41 días (`5b_rachas_ceros_prueba.csv`). Como testigo independiente, el satélite CHIRP vio en esas rachas entre 0,54 y 1,45 veces su lluvia normal; es un testigo débil, porque casi no distingue sequías de 30 días en el Valle.

**Por temporada y fase ENSO.** La temporada se define con el ciclo de entrenamiento **de cada estación**: un mes es seco si su lluvia media está por debajo de la media mensual de esa estación (`5b_rachas_ceros_escenarios.csv`):

| Temporada | Fase | Rachas de ≥30 días en 0 | R4 las vuelve faltante | % | Días que llovió la vecina más cercana (mediana) |
|---|---|---|---|---|---|
| Seca | El Niño | 476 | 5 | 1,1 | 2 |
| Seca | Neutro | 343 | 3 | 0,9 | 3 |
| Seca | La Niña | 58 | 2 | 3,4 | 4 |
| Lluviosa | El Niño | 14 | 0 | 0 | 3,5 |
| Lluviosa | Neutro | 18 | 0 | 0 | 5,5 |
| Lluviosa | La Niña | 5 | 0 | 0 | 6 |

**Lectura:**
- **El profesor tiene razón en el escenario.** El 96 % de las rachas largas en 0 cae en temporada seca, y más de la mitad en temporada seca de un año El Niño. Son sequías reales, y R4 las deja: en las rachas de El Niño la vecina llovió 2 días de mediana.
- La regla no necesita conocer la fase: las vecinas la "traen" incluida.
- **Si se quitara R4**, entrarían 21 curvas más (2 134 en vez de 2 113) y una estación más (Juntas), pero con 10 meses secos que casi seguro no ocurrieron. El modelo aprendería que el satélite exagera en temporada seca por un error del pluviómetro.

---

## 4. R2 y R3: extremos y totales mensuales

- **R2.** Un valor mayor que 3 veces el p99 de la estación y de al menos 150 mm, sin respaldo de ninguna vecina, es casi seguro un error de digitación. Comparar un extremo con las vecinas es otro control de consistencia espacial (Durre et al. 2010).
- **R3.** En registros manuales, cuando nadie lee el pluviómetro varios días, el siguiente valor trae el **acumulado** de todos esos días sin marcarlo. Viney y Bates (2004) mostraron que esto ocurre en buena parte de las estaciones de la red australiana de alta calidad, y que distorsiona la frecuencia de días con lluvia y el largo de las rachas.
  - La OMM 2017 (4.4.1 c) dice que un acumulado solo sirve si se sabe qué días cubre.
  - Nuestra versión extrema de ese error es todo el mes anotado en un solo día; ahí el mes pasa a faltante.

---

## 5. Métodos

### 5.1 No se modela la lluvia acumulada
**Lo que dijo el profesor:** «no modelen lluvia acumulada»; lo había dicho en clase.

**La fuente:** una curva acumulada es **creciente**. Para que la estimación lo siga siendo, la base tiene que llevar esa restricción, con suavizado monótono (Ramsay y Silverman 2005, cap. 6). Además, cada error se arrastra a todos los días siguientes, que es lo que el profesor describió como «se arrastra en todo ese rollo».

**Lo que hicimos:**
- Siempre modelamos la intensidad en mm/día, y la regresión trabaja sobre su raíz.
- Quitamos las dos gráficas de lluvia acumulada: `12_acumulada_por_anio` y `F4_acumulada`, junto con la sección 4.6 funcional.
- El **sesgo del total anual** que pide la rúbrica es un número por curva, la suma del año (un "total" en el sentido de la OMM). No es una curva acumulada, así que se mantiene.

### 5.2 No se usa Wilcoxon
**Lo que dijo el profesor:** es estadística clásica y no tiene en cuenta la estructura funcional.

**Lo que hicimos:**
- Para comparar dos modelos se calcula la métrica de cada curva, se promedia por estación y se toma la diferencia entre modelos.
- El IC del 95 % de esa diferencia sale de un **bootstrap por conglomerados**: se remuestrea la estación entera, porque sus años no son independientes (Davison y Hinkley 1997, sec. 3.8).
- Un modelo gana solo si el IC excluye el 0. El Wilcoxon no agregaba nada al IC y se quitó.

**Dónde fuimos críticos:** la métrica (error L², KGE) ya resume la curva entera, así que comparar las métricas es legítimo. Lo que no aportaba era la prueba de rangos.

### 5.3 No se usan la prueba t ni la F pareadas
**Lo que dijo el profesor:** la t pareada es para el mismo individuo antes y después de un tratamiento, y aquí no aplica.

**Dónde fuimos críticos:**
- Una estación y su píxel miden el mismo sitio con dos instrumentos. En los libros de comparación de métodos de medición eso suele tratarse como datos pareados, así que no es obvio que la prueba estuviera mal planteada.
- Pero la pregunta que respondía («¿el satélite está sesgado?») **ya la responde la regresión**, y de forma más completa. Si el satélite midiera lo mismo que la estación, la regresión daría intercepto α(t) = 0 e integral de β igual a 1. Las bandas de esos efectos dicen en qué días no es así (Ramsay y Silverman 2005, cap. 13–14).
- Hay además una tensión que conviene aclarar con el profesor. El ejercicio 3 de su clase de inferencia funcional (semana 14) pide «aplique la prueba F para muestras pareadas» justamente a `X_sat` y `Y_est`, satélite contra estación (`docs/BASE_CONOCIMIENTO_PROFESOR.md`, sección 2).
- Por eso las quitamos: no por un error de aplicación, sino porque sobraban. Y queda como pregunta para la asesoría: ¿qué diferencia ve entre ese ejercicio y nuestro caso?
- La curva de la diferencia media estación − satélite con su banda (sección 4.1 funcional) se deja como **descripción** de la comparación por fuente que pide la rúbrica, no como prueba.

### 5.4 Sí se usa el ANOVA funcional
Responde si el error del satélite, d(t) = estación − satélite, cambia según la fase ENSO, la zona o la altitud. Son grupos **independientes**:
- para el ENSO, los años;
- para la zona y la altitud, las estaciones, con la curva media de cada una.

Es el caso para el que se hizo la prueba (Cuevas, Febrero y Fraiman 2004; Zhang 2013) y la que vimos en clase. Así se demuestra con datos el efecto del ENSO, como pidió el profesor («ustedes pueden hacer inferencia funcional, entonces eso lo tienen que demostrar»).

**Sobre "Welch":** el profesor sospechó de una prueba de Welch. No hay ANOVA de Welch. Lo que hay es la aproximación de dos momentos de Satterthwaite para la prueba tipo F (Zhang 2013). Los valores p que usamos para decidir son por **permutación** (GPF y F_max; Górecki y Smaga 2015), que no dependen de esa aproximación.

### 5.5 La tesis de 2016 no es guía
El profesor dijo que mostró la tesis y los dos papers para que **no** se hiciera lo mismo. Quitamos todo lo que se tomó de ella sin otra fuente. Lo que queda tiene respaldo propio:
- **Paneles por año:** son descriptivos y muestran la lluvia diaria en mm/día.
- **Traza-variograma de los residuos:** Giraldo, Delicado y Mateu (2011).
- **Ajuste del variograma por mínimos cuadrados ponderados:** Cressie (1985).

---

## 6. R7: relleno de los días faltantes (D-61)

### Qué se rellena y qué no
- **Sí:** los días faltantes de los años que pasaron R5 (≤10 por mes y menos de 5 seguidos), incluidos los que R0 o R2 volvieron faltantes.
- **No:** los años que no pasan R5, las rachas de R4 y los meses de R3 (sus años ya cayeron en R5), y los años en que la estación no operaba.
- **Nunca con el satélite:** es lo que el modelo usa para predecir; rellenar con él sería circular.
- Si un día no tiene ninguna vecina con dato, queda faltante y lo cubre la curva.

### Por qué con las vecinas y no con la curva
Antes (D-13) no se imputaba: la B-spline cubría el hueco con los días de antes y después. La OMM (2017, 4.6) dice que, para valores diarios, los valores de otros días «are rarely useful» para estimar un día faltante, y pone primero la **interpolación espacial**: datos de sitios vecinos. Pide además estimar solo huecos pequeños (≤10 días en un mes) y sin sesgo.

Nuestros datos lo confirman (`8_R7_por_que_vecinas.csv`, entrenamiento):

| Comparación | Correlación de la lluvia diaria (mediana) |
|---|---|
| Vecina a ≤5 km, mismo día | 0,54 |
| Vecina a 5–10 km, mismo día | 0,44 |
| Vecina a 10–15 km, mismo día | 0,37 |
| Vecina a 15–25 km, mismo día | 0,31 |
| Vecina a 25–50 km, mismo día | 0,26 |
| La misma estación con su día anterior | 0,12 |

Incluso una vecina a 25–50 km dice el doble sobre el día faltante que la propia estación el día antes. Y el parecido baja con la distancia, que es lo que aprovechan IDW y el kriging.

### Por qué kriging
Los dos métodos estiman el día como un promedio ponderado de las vecinas ese mismo día. **IDW** fija los pesos con una fórmula (1/distancia²) y un número de vecinas escogidos a mano. El **kriging ordinario** (Cressie 1993) **mide** con el semivariograma cuánto baja el parecido con la distancia, reparte el peso entre vecinas que dicen lo mismo y da el error de cada estimación. Tiene precedente aquí: el profesor rellena con kriging ordinario en su tesis, y Terán-Chaves et al. (2021) rellenan con kriging ordinario tridimensional los datos diarios del IDEAM.

**Dos ajustes, comprobados con los datos:**
1. **Cada vecina se escala por la razón de medias mensuales** (lluvia media de la estación en ese mes / la de la vecina, en entrenamiento): es la razón normal. El kriging ordinario supone media constante, y en el Valle una vecina a ≤25 km llueve en promedio 1,28 veces lo que la otra.
2. **El semivariograma sale de la lluvia diaria de 1981–1995**, no de la media anual como en la tesis. El modelo se elige por validación cruzada (menor ECM), como en la tesis (`8_R7_semivariograma_modelos.csv`): gana el **exponencial** (ECM 88,7, frente a 89,0 el esférico, 92,2 el Matérn y 92,8 el gaussiano). Su rango llega al tope de 300 km: hasta 100 km la semivarianza sube casi en línea recta (de 66 a 2,5 km a 86 a 37,5 km), con una pepita alta. Es decir, la lluvia de un día se parece poco aun entre vecinas cercanas, y el kriging les da pesos casi parejos.

**La prueba** (`8_R7_metodos_relleno.csv`): 1 932 años completos de entrenamiento a los que se les quitan días, con el patrón real de una curva que pasa R5 (2 días de mediana) o con el máximo que acepta R5 (10 por mes, 120 al año). Error de la curva frente a la del año completo, en fracción de la variación entre años de la estación:

| Método | Patrón real | Máximo de R5 | Casos en que deja menos error que IDW (máximo de R5) |
|---|---|---|---|
| Solo la curva (D-13) | 0,032 | 0,551 | 14 % |
| IDW con 3 vecinas | 0,018 | 0,321 | — |
| Kriging sin la razón | 0,020 | 0,344 | 42 % |
| **Kriging con la razón (R7)** | **0,018** | **0,300** | **59 %** |

**Lectura:** rellenar con vecinas es mucho mejor que dejar que la curva cubra el hueco. Entre los métodos con vecinas las diferencias son pequeñas; el kriging con la razón es el mejor y además sus pesos salen de los datos.

### Cuánto se rellena
- **Conjunto principal (≥10 años):** 719 días de 1 111 790 (0,065 %). Otros 36 días no tienen ninguna vecina con dato ese día: quedan faltantes y los cubre la curva.
- **Sensibilidad (OMM):** 366 días de 771 245 (0,047 %), sin faltantes que queden.
- En todo el Valle, 944 de los 980 faltantes de años válidos tienen vecina con dato (mediana: 10 vecinas).
- Cada día rellenado queda marcado (`8_R7_dias_rellenados.csv` y la matriz `rellenado` de los datos preparados), y R5 se mide siempre con los días **observados**.

---

## Referencias

- Cressie, N. (1985). Fitting variogram models by weighted least squares. *Mathematical Geology*, 17(5), 563–586.
- Cressie, N. (1993). *Statistics for Spatial Data*, ed. revisada. Wiley (kriging ordinario).
- Lewis, E. et al. (2021). Quality control of a global hourly rainfall dataset. *Environmental Modelling & Software*, 144, 105169. Código: github.com/nclwater/intense-qc.
- Lopes Martins, L. L. et al. (2022). Gap-filling of daily precipitation and streamflow time series: a method comparison at random and sequential gaps. *Hydrological Sciences Journal*, 68(1), 148–160.
- NOAA NCEI (2021). *Normals Calculation Methodology 2020*.
- Terán-Chaves, C. A., Duarte-Carvajalino, J. M. y Polo-Murcia, S. M. (2021). Quality control and filling of daily temperature and precipitation time series in Colombia. *Meteorologische Zeitschrift*, 30(6), 489–501.
- Thornton, P. E. et al. (2021). Gridded daily weather data for North America with comprehensive uncertainty quantification. *Scientific Data*, 8 (PMC8302764).
- Urrea, V., Ochoa, A. y Mesa, O. (2019). Seasonality of rainfall in Colombia. *Water Resources Research*, 55, 4149–4162.
- WMO (2017). *WMO Guidelines on the Calculation of Climate Normals*. WMO-No. 1203 (secciones 4.4 y 4.6).
- Cuevas, A., Febrero, M. y Fraiman, R. (2004). An anova test for functional data. *Computational Statistics & Data Analysis*, 47(1), 111–122.
- Davison, A. C. y Hinkley, D. V. (1997). *Bootstrap Methods and their Application*. Cambridge University Press (sec. 3.8, datos jerárquicos).
- Durre, I., Menne, M. J., Gleason, B. E., Houston, T. G. y Vose, R. S. (2010). Comprehensive automated quality assurance of daily surface observations. *Journal of Applied Meteorology and Climatology*, 49(8), 1615–1633.
- Funk, C. et al. (2015). The climate hazards infrared precipitation with stations — a new environmental record for monitoring extremes. *Scientific Data*, 2, 150066.
- Giraldo, R., Delicado, P. y Mateu, J. (2011). Ordinary kriging for function-valued spatial data. *Environmental and Ecological Statistics*, 18(3), 411–426.
- Górecki, T. y Smaga, Ł. (2015). A comparison of tests for the one-way ANOVA problem for functional data. *Computational Statistics*, 30(4), 987–1010.
- Ramsay, J. O. y Silverman, B. W. (2005). *Functional Data Analysis*, 2.ª ed. Springer (cap. 6, funciones con restricciones; cap. 13–14, modelos lineales funcionales).
- Viney, N. R. y Bates, B. C. (2004). It never rains on Sunday: the prevalence and implications of untagged multi-day rainfall accumulations in the Australian high quality data set. *International Journal of Climatology*, 24(9), 1171–1192.
- WMO (1989). *Calculation of Monthly and Annual 30-Year Standard Normals*. WCDMP-No. 10, WMO-TD No. 341.
- WMO (2007). *The Role of Climatological Normals in a Changing Climate*. WCDMP-No. 61, WMO-TD No. 1377.
- WMO (2017). *WMO Guidelines on the Calculation of Climate Normals*. WMO-No. 1203 (sección 4.4).
- Zhang, J.-T. (2013). *Analysis of Variance for Functional Data*. Chapman & Hall/CRC.

# Bitácora de decisiones — Proyecto final Analítica de Datos 2026-2S

Precipitación del Valle del Cauca 1981–2013: estaciones diarias vs. satélite pentadal, datos funcionales, regresión funcional y red KAN.

Este documento registra **cada decisión que tomamos, por qué y cómo la validamos**, para que cualquier integrante pueda defenderla en la sustentación. Cuando una decisión cambie, no se borra: se marca como reemplazada y se agrega la nueva.

**Estados:** ✅ Decidida · 🟡 Propuesta (falta confirmar con el profesor) · 🔁 Reemplazada

**Cómo agregar una decisión:** copiar la plantilla del final, darle el siguiente número (D-xx) y llenarla.

---

## Índice

0. [Fe de erratas y rama de trabajo](#0-fe-de-erratas-y-rama-de-trabajo-12-nov-2026)
1. [Alcance geográfico](#1-alcance-geográfico)
2. [Ubicación de las estaciones](#2-ubicación-de-las-estaciones)
3. [Calidad de los datos de estaciones](#3-calidad-de-los-datos-de-estaciones)
4. [Producto satelital](#4-producto-satelital)
5. [Fuentes externas](#5-fuentes-externas)
6. [Construcción funcional](#6-construcción-funcional)
7. [Regresión funcional y red](#7-regresión-funcional-y-red)
8. [Pendientes y preguntas al profesor](#8-pendientes-y-preguntas-al-profesor)

---

## 0. Fe de erratas y rama de trabajo (12-nov-2026)

### D-28 · Rige la fe de erratas del enunciado ✅
- **Decisión:** la norma vigente es `Proyecto_Final_Analitica_FeErratas.docx`. Cambia cuatro cosas: el periodo pasa a **1981–2003** (antes 1981–2013) en la introducción, el objetivo y la rúbrica del EDA, y la partición pasa a **entrenamiento 1981–1995 / validación 1996–2000 / prueba 2001–2003** más las estaciones reservadas por el docente.
- **Justificación:** según el correo del docente, la base de estaciones no se recuperó hasta 2013. Comparación párrafo a párrafo en `COMPARACION_WORD_FE_ERRATAS.md`.
- **Consecuencias:** reemplaza D-17 y D-24; las péntadas de 2004–2013 quedan fuera; ya no hay que preguntar por 2004–2013 (P-4).

### D-29 · Se trabaja sobre el flujo de la rama de María y Jacobo ✅
- **Decisión:** la base del trabajo es la rama `maria-eda-funcional-construccion-funcional` (fases 0–4 ya reproducidas con la fe de erratas). Nuestro aporte entra en la rama `gus-estaciones-valle`, creada a partir de ella.
- **Qué cambia frente a lo que habíamos hecho en `main`:**
  - **Depuración:** reglas P1–P6; el mínimo de años por estación (P6) se cuenta **solo en entrenamiento** (≥7 años válidos en 1981–1995, antes ≥10 en 1981–2003) y los umbrales relativos (percentiles) se calculan por estación con entrenamiento. Población: **248 estaciones / 4 485 pares** (antes 238 / 4 426). Reemplaza el umbral de D-11.
  - **K:** 8 para estaciones y 11 para satélite, elegidos solo con entrenamiento. Reemplaza D-21 (K = 6).
  - **Fase 3:** 2 751 pares estación–píxel en 147 estaciones; 78 estaciones sin coordenada en el archivo del curso y 23 fuera de la malla.
  - **Fase 4:** `05_red_funcional.py` (codificador con W, cabezas lineal/KAN/MLP, FPCA, ablación con 3 semillas). La equivalencia lineal ↔ regresión funcional se verifica a 4e‑15. Reemplaza nuestro `05_regresion_funcional.ipynb` (D-25 a D-27), que queda solo como referencia local.
- **Conflictos encontrados al integrar:** ninguno de archivos (nuestros archivos nuevos no existen en la rama; solo el `README.md` se tocó en ambos lados y se resolvió agregando la sección "Decisiones").

### D-30 · Qué aporta nuestro trabajo de estaciones a esa rama 🟡
- **Filtro del Valle (D-01):** la rama **no filtra por departamento**. De sus 248 estaciones, **167 son del Valle (3 122 pares)** y 81 son del Cauca o de departamento desconocido (1 363 pares). Falta aplicar el filtro en el flujo y rehacer las fases.
- **Coordenadas:** de las 78 estaciones sin coordenada en la rama, **27 se recuperan** con `datos_web/estaciones_valle_del_cauca.xlsx` (529 pares), más La Montañuela (`2624700111`) y Limones (`5421220101`), que entran por la nueva regla P6 y están en el catálogo CVC (18 pares; falta verificarlas con el criterio D-06).
- **Estaciones del Valle nuevas por la regla P6** (8): La Siberia, El Janeiro, El Jigual, La Leonera, Calabazas, Diguita (con coordenadas en el archivo del curso), La Montañuela y Limones. Tres del Valle salen por la nueva depuración: `2622100104` Peña Mona, `2624600101` Bolívar y `5420410112` La María.
- **Pendiente:** actualizar el Excel a la población de la rama y conectar sus coordenadas al cruce estación–píxel (`scripts/crear_pares_estacion_pixel.R`).

---

## 1. Alcance geográfico

### D-01 · Solo estaciones del departamento del Valle del Cauca ✅
- **Decisión:** el análisis usa únicamente estaciones ubicadas en el departamento del Valle del Cauca. Las de otros departamentos se excluyen y se reportan aparte.
- **Justificación:** el enunciado del proyecto delimita el estudio al Valle del Cauca. Los datos entregados incluyen estaciones del departamento del Cauca (Caloto, Corinto, Toribío, Puerto Tejada, Santander de Quilichao, Silvia, Suárez, Puracé, El Tambo, Popayán, etc.).
- **Cómo se aplicó:** se asignó el departamento por el municipio de cada estación (lista de los 42 municipios del Valle), usando los dos catálogos del curso y el catálogo oficial de la CVC. Los catálogos coinciden en el departamento de todas las estaciones que comparten.
- **Resultado (población de `main`; la de la rama actual está en D-30):** de 238 estaciones válidas → **162 del Valle** (3 092 pares estación-año), 37 del Cauca (712 pares) y 39 de departamento desconocido (622 pares).
- **Archivo:** `datos_web/estaciones_valle_del_cauca.xlsx` (hojas "Estaciones Valle" y "Excluidas").

### D-02 · Las estaciones de departamento desconocido se excluyen ✅
- **Decisión:** las 39 estaciones que no aparecen en ningún catálogo quedan fuera.
- **Justificación:** de ellas solo tenemos el código (los archivos de datos no traen nombre). Sin ubicación no se puede confirmar que sean del Valle ni emparejarlas con un píxel del satélite. Varios códigos (prefijos `2611`, `2612`) son de subcuencas que también cruzan el departamento del Cauca.
- **Cómo se validó:** se buscaron por código en el catálogo de la CVC (datos.gov.co) y en el CNE del IDEAM, y en la web. Ninguna apareció.
- **Reversible:** si el profesor o la CVC dan su ubicación y son del Valle, se incluyen.

### D-03 · Recortar el satélite al polígono del Valle 🟡
- **Decisión propuesta:** para el EDA del satélite y los promedios regionales, usar solo los píxeles dentro del límite del departamento.
- **Justificación:** el dominio entregado (lon −77.55 a −75.70, lat 3.05 a 5.00) incluye partes de Chocó, Cauca, Risaralda, Quindío y Tolima.
- **Falta:** descargar el polígono oficial del Valle (DANE/IGAC) y rehacer las fases afectadas.

---

## 2. Ubicación de las estaciones

### D-04 · Fuentes de coordenadas y orden de prioridad ✅
1. **Archivo del curso** `CVC_Wheather_Stations_Long_Lat.csv` (coordenadas en grados). Usado para 134 estaciones del Valle.
2. **Catálogo oficial de la CVC** en datos.gov.co (id `kfsz-ija2`): mismos códigos de 10 dígitos que nuestros datos, con coordenadas proyectadas.
3. **Catálogo Nacional de Estaciones del IDEAM (CNE)** en datos.gov.co (id `hp9r-jxuu`): fuente independiente para confirmar ubicación, altitud y fechas.
4. **Documentos oficiales** (POT, planes de manejo) cuando existe una medición GPS.
- **Justificación:** el archivo del curso solo tenía coordenadas de 166 de nuestras 238 estaciones. El catálogo de la CVC permite buscar por el mismo código; el CNE permite verificar con una segunda fuente.

### D-05 · Conversión de coordenadas del catálogo CVC ✅
- **Decisión:** las coordenadas planas del catálogo CVC se convierten a grados (WGS84) con `pyproj`: **EPSG:3115** (MAGNA-SIRGAS Colombia Oeste) para la mayoría y **EPSG:9377** (CTM12) para las que tienen Este > 3 000 000.
- **Cómo se validó:** en 191 estaciones presentes en el catálogo CVC y en el archivo del curso, la diferencia mediana entre ambas ubicaciones es **0.04 km**.

### D-06 · Criterio para considerar una ubicación "segura" ✅
Una estación que no está en el archivo del curso se acepta si pasa **al menos dos pruebas independientes**:
1. **Catálogo:** mismo código, nombre y municipio en el catálogo CVC, y la misma estación en el CNE con fecha de inicio anterior a nuestros datos y ubicación a ≤4 km.
2. **Altitud:** la elevación del terreno (SRTM 30 m, servicio OpenTopoData) en la coordenada coincide con la altitud registrada (Excel del curso o CNE).
3. **Vecinas:** la lluvia mensual de la estación se correlaciona más con las estaciones cercanas que con las lejanas.
- **Referencia para la prueba 3:** en todas las estaciones confirmadas, la correlación promedio de las anomalías mensuales (log(1+P)) baja con la distancia: 0.70 a ≤10 km, 0.61 a 10–20 km, 0.56 a 20–40 km, 0.51 a 40–80 km y 0.46 a >80 km. La correlación de Spearman entre r y distancia tiene mediana −0.43 (P10 −0.57, P90 −0.22). Una estación mal ubicada no seguiría este patrón.
- **Limitación:** la prueba 3 no distingue entre puntos a menos de ~10 km; para eso se usa la altitud.

### D-07 · Estaciones recuperadas con el catálogo CVC + CNE ✅
19 estaciones que estaban solo en el Excel de altitudes del curso: mismo código y nombre en el catálogo CVC, el CNE confirma ubicación (≤4 km) y altitud.

| Código | Estación | Municipio |
|---|---|---|
| 2610000202 | Zaragoza | Cartago |
| 2612710301 | La Diana | Florida |
| 2612720301 | Chambú | Palmira |
| 2612810201 | San Emigdio | Palmira |
| 2613600301 | La Magdalena | Buga |
| 2614100202 | Acueducto Tuluá | Tuluá |
| 2614100203 | Monteloro | Tuluá |
| 2614900201 | Miravalles | La Victoria |
| 2621500301 | Samarkanda | Jamundí |
| 2622100301 | San Antonio | Jamundí |
| 2622110201 | El Topacio | Cali |
| 2622330104 | El Descanso | Cali |
| 2622330107 | Edificio CVC | Cali |
| 2622400201 | La Teresita | Cali |
| 2622600201 | La Buitrera | Yumbo |
| 2622900301 | Vijes | Vijes |
| 2624700119 | La Calera | La Unión |
| 2624700301 | Sabanazo | La Unión |
| 5300000203 | Bahía Málaga | Buenaventura |

### D-08 · Estaciones verificadas una por una ✅
Estas no pasaban el criterio de D-07 a la primera (fechas que no cuadraban, solo coincidencia por código, o el CNE daba otra ubicación). Cada una se resolvió así:

| Código | Estación | Problema inicial | Cómo se resolvió | Ubicación usada |
|---|---|---|---|---|
| 2613600105 | Acueducto Buga | El archivo del curso dice que empezó en nov-1999, pero hay datos desde 1981 | El CNE registra ACUEDUCTO BUGA de la CVC (cód. 26090790) instalada en jul-1978 a ~2 km. Sus vecinas más parecidas son La Primavera y Mediacanoa (7–10 km, r≈0.71–0.73). La fecha de 1999 parece una reinstalación | Archivo del curso |
| 2614400201 | Barragán | El catálogo CVC dice que empezó en 2015, pero hay datos 1981–2000 | El CNE registra BARRAGAN (IDEAM, cód. 26105140) en Tuluá desde sep-1972, **en las mismas coordenadas**, a 2 902 m. La fecha de 2015 es la nueva estación automática en el mismo sitio. Vecinas más parecidas: Los Bancos, La Gitana, Alegrías (r≈0.70–0.72) | Catálogo CVC |
| 2612700201 | Ing. Manuelita | Solo coincidía por código | CNE: INGENIO MANUELITA (IDEAM) desde 1935, a 0.7 km. Terreno 1 034 m vs. 1 059 m registrados | Catálogo CVC |
| 2614100102 | El Placer | Solo coincidía por código | CNE: EL PLACER (IDEAM) desde 1971 en el mismo punto. Terreno 2 125 m vs. 2 132 m | Catálogo CVC |
| 2620000109 | Mediacanoa | Solo coincidía por código | CNE: MEDIACANOA-ALERTAS (IDEAM) desde 1979 a 0.3 km. Terreno 938 m. No está duplicada: la otra Mediacanoa (2610000118) no pasó el filtro de calidad | Catálogo CVC |
| 2620000202 | Garzonero | El CNE la ubicaba a 11 km | El CNE tiene **dos** Garzonero de la CVC; se había cruzado con la pluviométrica (terreno 1 654 m). La climatológica (cód. 26085140) está a ~2 km del punto CVC y tiene **exactamente** la altitud del Excel del curso (942 m). Terreno en el punto CVC: 953 m | Catálogo CVC |
| 2623300201 | Bosque Yotoco | El CNE la ubicaba a 6.5 km y con 140 m de diferencia | El punto CVC cae en la reserva Bosque de Yotoco con terreno a 1 601 m (Excel: 1 700 m). Las coordenadas del CNE están redondeadas al minuto y caen en el valle (964 m): el impreciso era el CNE | Catálogo CVC |
| 2612810203 | La Sirena | En el punto CVC el terreno está a 3 474 m, pero el Excel dice 2 605 m | El **plan de manejo del PNR del Nima (POT de Palmira)** da la ubicación medida con GPS Garmin (error ~1 m): Este 1108144.877, Norte 881317.939, a 2 675 m. Convertida: **3.52242, −76.104246**, terreno a **2 674 m**. El catálogo CVC la tenía corrida ~3.6 km y en **otro píxel del satélite**. Los datos empiezan en 1988, igual que la estación climatológica de la CVC | GPS del POT |
| 5421240201 | El Tábano | Registrada a 380 m, pero el terreno está a 874 m (punto CVC) y 559 m (punto CNE) | El CNE la registra 15/09/1981–15/02/1992; nuestros datos empiezan en sep-1981 y son continuos hasta ene-1992 (coinciden). Las coordenadas del CNE están redondeadas al minuto (±1 km); revisando el terreno cada 200 m dentro de ese cuadro, el mínimo es 390 m en el fondo del valle, **3.9458, −76.7063**. Lluvia de 7 400–9 500 mm/año, coherente con la vertiente Pacífica baja. Incertidumbre ~1 km, mismo píxel satelital que el punto CNE | Estimada con terreno |

- **Resultado:** las 162 estaciones del Valle tienen ubicación verificada (134 del archivo del curso + 28 verificadas).
- **Hallazgo para el informe:** el catálogo de la CVC tiene errores de ubicación de varios kilómetros en algunas estaciones (La Sirena, El Tábano). Por eso se cruzaron tres fuentes en vez de confiar en una sola.

### D-09 · Las 162 estaciones son el total correcto para nuestros datos ✅
- **Pregunta:** ¿cuántas estaciones de lluvia tiene el Valle? No hay un único total: los datos del curso son un subconjunto de la red de la CVC.
- **Cuadre con el catálogo oficial CVC:** 219 estaciones de lluvia del Valle (pluviográficas, pluviométricas y climatológicas) iniciadas antes de 2004:
  - 123 están en nuestras 162;
  - 42 vienen en los datos del curso pero no pasaron el filtro de calidad (mediana de 7.5 años con datos en 1981–2003);
  - 54 no vienen en los datos entregados por el curso.
- 31 de nuestras 162 no aparecen por código en el catálogo actual (vienen con coordenadas del archivo del curso; el catálogo las tiene con códigos nuevos o ya no las lista).
- El CNE cuenta 966 estaciones de lluvia en el Valle operando en algún año de 1981–2003 sumando todas las entidades (IDEAM, Cenicaña, cafeteros, particulares); de esas, 410 son de la CVC. No es comparable porque el curso solo entregó estaciones de la CVC.
- **Conclusión:** no se perdió ninguna estación por error nuestro.

### D-10 · Altitud de las estaciones 🟡
- **Decisión actual:** se usa la altitud del Excel del curso (`CatalogoEstacionesElevacion.xlsx`); para La Sirena y El Tábano, la verificada (2 675 m y 380 m).
- **Pendiente:** 40 estaciones del Valle no tienen altitud en el Excel (37 del archivo del curso y 3 verificadas: Ing. Manuelita, El Placer y Mediacanoa; para estas 3 el CNE da 1 059, 2 132 y 920–962 m). Propuesta: tomarla del CNE cuando exista o del modelo de elevación SRTM en su coordenada, después de validar SRTM con las estaciones que sí tienen altitud. Confirmar con el profesor.

---

## 3. Calidad de los datos de estaciones

### D-11 · Criterio de cobertura de un año-estación 🔁 (umbral de años reemplazado por D-29)
- **Decisión:** un año de una estación es válido si tiene ≥90 % de días con dato, ningún hueco mayor a 30 días, los 12 meses con ≥20 días, y la estación tiene ≥10 años válidos.
- **Justificación:** una curva anual con huecos largos no se puede suavizar sin inventar la forma de la temporada de lluvias; el umbral de 10 años asegura que cada estación aporte variabilidad interanual.
- **Resultado del embudo (todas las estaciones, antes del filtro del Valle):** 323 estaciones / 5 514 pares → 238 / 4 426. Detalle en `eda_fase0/embudo_depuracion_preliminar.csv`.

### D-12 · Días marcados como erróneos → faltantes ✅
- **Decisión:** se marcan como faltantes (no como cero): totales mensuales registrados en un solo día, rachas de ceros falsos (≥30 días en cero mientras la región llueve, contraste con p_lluvia ≥ 0.4) y años copiados entre estaciones.
- **Justificación:** son errores de digitación o de registro, no lluvia real; tratarlos como dato sesgaría las curvas.
- **Archivos:** `eda_fase0/dias_marcados.csv`, `totales_mensuales_un_dia.csv`, `rachas_ceros_30d_con_contraste_regional.csv`, `anios_copiados_entre_estaciones.csv`.

### D-13 · No se imputan datos de estaciones ✅
- **Decisión:** no se rellenan días faltantes con otras fuentes ni con modelos.
- **Justificación:** la base B‑spline se ajusta solo en los días observados, así que no hace falta imputar para construir la curva. Rellenar con satélite o CHIRPS contaminaría justamente lo que se quiere comparar (estación vs. satélite).

### D-14 · Dominio anual común ✅
- **Decisión:** día del año 1–365; el 29 de febrero se elimina en años bisiestos.
- **Justificación:** todas las curvas deben vivir en el mismo dominio para compararlas y promediarlas.

---

## 4. Producto satelital

### D-15 · Unidades del satélite: mm/día ✅
- **Decisión:** cada péntada se divide por su número real de días (3 a 6) y se ubica en su día central.
- **Justificación:** las péntadas no tienen la misma duración (la 6.ª de cada mes va del 26 al fin de mes). Sin dividir, la última péntada de los meses de 31 días aparece artificialmente más lluviosa.

### D-16 · Años defectuosos del satélite: 1985, 1988 y 2000 ✅
- **Decisión:** se excluyen del entrenamiento de los modelos y se reportan aparte.
- **Justificación:** en el EDA tienen déficit de lluvia anómalo frente a las estaciones y a CHIRPS v2. En CHIRPS v2 solo desaparece el déficit de 2000.

### D-17 · Satélite: 1981–2012 🔁 (reemplazada por D-28: 1981–2003)
- **Decisión:** 2013 se excluye porque solo tiene 24 péntadas.

---

## 5. Fuentes externas

### D-18 · Las fuentes externas solo sirven para comparar el satélite, no para reemplazar estaciones ✅
- **Decisión:** CHIRPS v2 y NASA POWER se usan como referencia para evaluar el producto satelital entregado. No se usan para rellenar ni reemplazar los datos de estaciones.
- **Justificación:**
  - **CHIRPS v2** incorpora pluviómetros (probablemente algunos de nuestras mismas estaciones), así que no es independiente para validar.
  - **NASA POWER** (MERRA-2) es independiente pero muy gruesa (0.5° × 0.625°).
- **Pendiente con el profesor:** si se permite reemplazar el satélite entregado por CHIRPS (ver §8).

### D-19 · ENSO con el índice ONI de NOAA ✅
- **Decisión:** las fases El Niño / La Niña / neutral se toman del ONI (NOAA CPC). Archivo `datos_web/ONI_NOAA_CPC_1981_2013.csv`.

---

## 6. Construcción funcional

### D-20 · Base B‑spline cúbica y escala mm/día ✅
- **Decisión:** B‑splines de orden 4 sobre el día del año, programadas paso a paso (Cox–de Boor) y verificadas contra `scipy`. Se ajustan en mm/día.
- **Justificación:** mm/día es la misma magnitud en ambas fuentes (tasa media). Si se transforma antes de suavizar (√ o log), las fuentes dejan de ser comparables porque la lluvia diaria es intermitente.

### D-21 · Número de bases K por VCG 🔁 (reemplazada por D-29: K = 8 / 11 solo con entrenamiento)
- **Decisión:** K = 6 para estaciones y K = 11 para satélite, elegidos por validación cruzada generalizada (VCG) en mm/día.
- **Justificación:** criterio usado en clase; no hace falta el mismo K en ambas fuentes porque se comparan evaluando en la misma rejilla de días. Detalle en `fda_fase1/resumen_seleccion_K.csv`.

---

## 7. Regresión funcional y red

### D-22 · Tipo de modelo: regresión función a función ✅
- **Decisión:** y(t) = α(t) + ∫β(s,t)x(s)ds, con x(s) la curva satelital y y(t) la de la estación.
- **Justificación:** el enunciado pide que la cabeza lineal de la red con decodificador FPCA sea equivalente a una regresión funcional función a función de rango bajo, y que se demuestre y verifique (diferencia ≤ 1e‑4).

### D-23 · Escala √P, aplicada después de suavizar ✅
- **Decisión:** se modela √P (lo pide el enunciado), pero la raíz se aplica a la curva ya suavizada en mm/día, no a la lluvia diaria.
- **Justificación y validación:** aplicar la raíz día a día subestima la tasa (con las estaciones: 1.8 mm/día frente a 5.4 mm/día reales). Volver a mm es elevar al cuadrado lo que se recorta en 0.

### D-24 · Partición entrenamiento / prueba 🔁 (reemplazada por D-28)
- **Propuesta:** prueba = años 1999, 2001, 2002, 2003 + 20 % de estaciones reservadas (semilla 42); 2000 fuera de todo; 1985 y 1988 fuera del entrenamiento.
- **Justificación:** permite evaluar años nuevos, estaciones nuevas y ambos. Archivo `fda_fase4_regresion/particion.csv`.

### D-25 · FPCA, Gram y penalización 🔁 (ver `05_red_funcional.py`, D-29)
- **FPCA** solo con entrenamiento; M = 4 componentes (≥90 % de varianza), verificado con `scikit-fda`.
- **Matriz de Gram** por Gauss–Legendre (exacta para B‑splines), verificada contra trapecio fino (diferencia 4e‑8).
- **Penalización:** rugosidad en s + un término pequeño de tamaño (P = R + δW), para que con λ grande el modelo pueda llegar a "no usar el satélite" (β = 0) y la validación cruzada pueda elegirlo.
- **λ** por validación cruzada dejando un año fuera (no una curva), porque las estaciones del mismo año comparten la entrada.

### D-26 · Verificación de la equivalencia red ↔ regresión 🔁 (ver D-29)
- Identidad red ↔ superficie β: diferencia máxima **1.5e‑14**.
- Red lineal entrenada por gradiente vs. solución cerrada: **1.3e‑10**.
- Ambas cumplen el criterio del enunciado (≤ 1e‑4). Archivo `fda_fase4_regresion/verificacion_equivalencia.csv`.

### D-27 · Entrada satelital provisional (regional) 🔁
- **Decisión anterior:** mientras no había coordenadas, la entrada x(s) fue la curva satelital promedio de la región por año.
- **Reemplazada por:** usar el píxel de cada estación con las coordenadas de §2. Falta correr de nuevo `05_regresion_funcional.ipynb` cambiando `curva_satelite_de()`.

---

### D-31 · Inventario completo de los 344 códigos y fuentes oficiales por código ✅
- **Qué se hizo:** se sacaron los códigos únicos (`STATION_ID`) de los tres archivos de lluvia: P_1930_1979 (264), P_1980_1990 (302) y P_1991_2003 (284) → **344 códigos únicos**. Se cruzaron con el CSV de lat/lon del curso (208 están) y el Excel de altitudes (186 están); 84 no están en ninguno.
- **Búsqueda por código en fuentes oficiales** (usan el mismo código de 10 dígitos):
  1. Catálogo de estaciones de la CVC en datos.gov.co (`kfsz-ija2`): 60 de los 136 sin coordenadas.
  2. Registros de estaciones automáticas de la CVC en datos.gov.co (`avya-p282`): `2612710202` = **CENICAÑA** (Florida).
  3. Geoportal GeoCVC, capa `Externos/Hidroclimatologia_Sec_Gestion_riesgo`: `2622320201` = **UNIVALLE** (Cali, 970 m, desde 1965).
  4. GeoCVC `Cambio_Climatico/Red_Hidroclimatologica` (250 estaciones): usa la codificación nueva, sin coincidencias.
  5. CNE del IDEAM: usa códigos propios; sirve para confirmar por nombre, no por código.
  6. Búsqueda web de los códigos: sin resultados útiles.
- **Resultado:** 272 códigos con ubicación oficial (226 Valle, 46 Cauca). Cenicaña y Univalle entran a "Estaciones Valle" como seguras → **164 estaciones del Valle con datos válidos 1981–2003 (3 133 pares)**.
- **Los 72 sin ninguna fuente** se clasificaron con dos pistas independientes:
  - **Cuenca por código:** los primeros 8 (o 7) dígitos identifican la cuenca; si todas las estaciones oficiales de esa cuenca están en municipios de un solo departamento, se asigna ese. No aplica a los prefijos del cauce principal del río Cauca (`2610000`, `2620000`), que cruzan varios departamentos. Ojo: `5300900` es Argelia **del Cauca** (río Hondo), no Argelia del Valle.
  - **Vecinas:** a qué estaciones ubicadas se parece su lluvia mensual (anomalías de log(1+P), 1930–2003). Validado con las 226 ubicadas que tienen datos suficientes: acierta el departamento en 214 de 226 (95 %).
  - **Resultado:** 41 Valle (30 con las dos pistas de acuerdo, 11 con una sola), 10 Cauca, 21 sin evidencia suficiente o con pistas contradictorias. Detalle en la hoja "Inventario 344 códigos" de `datos_web/estaciones_valle_del_cauca.xlsx`.
- **Uso propuesto:** las 22 del Valle inferidas que tienen años válidos pueden entrar al EDA y al FPCA (no necesitan ubicación), marcadas como "Valle por inferencia". No entran al emparejamiento con el satélite hasta tener coordenadas (pedirlas a la CVC o al profesor).

### D-32 · Estructura del código, recuperación por cruce con el CNE y uso de la pista de vecinas ✅
- **Orden de trabajo:** primero se completa nombre, latitud y longitud de los 344 códigos y se filtra por departamento (Valle). **No se filtra por años válidos en esta etapa**: ese filtro va en el EDA, con su justificación.
- **Qué significa el código (10 dígitos, codificación antigua de la CVC):** 1-2 zona hidrográfica (26 Cauca, 53 Pacífico, 54 San Juan); 3 margen del río Cauca (1 derecha/este, 2 izquierda/oeste); 4 tramo de sur a norte; 5-6 cuenca; 7-8 tipo de estación (01 pluviométrica/pluviográfica, 02 climatológica, 03 evaporimétrica, 04 limnigráfica); 9-10 consecutivo.
  - **Fuente:** el sistema de codificación de unidades hidrográficas que describe la CVC en el POMCA del río La Vieja (cap. 6), basado en la zonificación del IDEAM (área, zona, margen, cuencas numeradas de sur a norte).
  - **Comprobación con datos oficiales:** el dígito de margen coincide con el lado del río en 211 de 229 estaciones; las latitudes suben en orden con el dígito de tramo; los dígitos 7-8 coinciden con el tipo registrado en el catálogo CVC.
  - **Consecuencia:** el CSV de coordenadas del curso es el catálogo de la red pluviométrica (solo PM/PG = 01); por eso faltaban las climatológicas y evaporimétricas (02/03).
- **Recuperación por cruce con el CNE del IDEAM (1 228 estaciones de Valle y Cauca):** para cada código sin ubicación se buscaron estaciones del CNE que cumplan todo: miden lluvia; tipo compatible con los dígitos 7-8; operaban durante los años de datos del código; a ≤ 40 km de las estaciones oficiales de su cuenca (dígitos 1-6); no corresponden ya a otro código nuestro. Se identifica solo si además: el año de instalación coincide con el primer año de datos (±1), hay pocas candidatas compatibles, está a ≤ 15 km de su cuenca, ninguna otra estación nuestra apunta a la misma, y la lluvia es coherente con el sitio (ver abajo).
  - **Identificadas (7, todas del Valle):** 2612700202 = PALMIRA ICA (1930); 2612700203 = AEROPUERTO ALFONSO BONILLA ARAGÓN (1971); 2612710201 = CENTRAL CASTILLA, Pradera (1961); 2624000301 = MANUEL M. MALLARINO, Trujillo (1969); 2624700202 = CENT. ADMO. LA UNIÓN (1967); 5301400105 = ALTO ANCHICAYÁ (1967; verificar que no duplique a 5301400108); 5301420201 = EL QUEREMAL, Dagua (1966).
  - **Candidatas sin asignar (41):** hay una estación compatible pero no única o con evidencia débil; quedan anotadas en la columna "Candidata del CNE" para revisión manual, sin coordenadas asignadas.
- **Pista 2 (lluvia de vecinas): cómo se usa y cómo nos cubrimos.**
  - **Fundamento:** la lluvia está espacialmente correlacionada (las estaciones cercanas se parecen más); es el principio de las pruebas de homogeneidad y control de calidad con estaciones de referencia.
  - **Límites conocidos:** la orografía puede hacer que dos estaciones cercanas lluevan distinto (y lejanas parecido); no distingue bien a menos de ~10 km; falla cerca de los límites entre departamentos; en el Pacífico la red es dispersa y la correlación es baja.
  - **Reglas para cubrirnos:** (1) nunca se usa sola para dar coordenadas, solo como verificación o para el departamento; (2) se validó antes de usarla: acierta el departamento en 214 de 226 estaciones con ubicación oficial (95 %), y ese error se reporta; (3) se exige que coincida con otra pista independiente (código o CNE); (4) cada estación queda etiquetada con su nivel de evidencia; (5) en el EDA se corre un análisis de sensibilidad con y sin las estaciones de nivel 3.
- **Resultado (344 códigos):** Valle 267 (226 ubicación oficial + 7 identificadas por CNE + 34 solo departamento inferido); Cauca 56; desconocido 21. Hojas "Inventario 344 códigos" y "Valle (sin filtrar años)" de `datos_web/estaciones_valle_del_cauca.xlsx`.

## 8. Pendientes y preguntas al profesor

| # | Pregunta / pendiente | Estado |
|---|---|---|
| P-1 | Los 21 códigos sin ubicar y los 35 con departamento inferido se descartan y no se usan. Respuesta del profesor (13-nov-2026) | ✅ |
| P-2 | Resuelta (D-35): no se usa SRTM para la altitud | ✅ |
| P-3 | No: se usa solo el satélite entregado, sin rellenar ni reemplazar con otra fuente (CHIRPS u otra). Los problemas del satélite se tratan con el modelo. Respuesta del profesor (13-nov-2026) | ✅ |
| P-4 | ¿Cómo tratar el año 2000 del satélite, que cae en validación? Según D-38 no está dañado: es un año La Niña que el satélite subestima | 🟡 |
| P-9 | Resuelta por el enunciado: el ONI del año y la altitud se agregan al vector z de la red | ✅ |
| P-5 | Resuelta por la fe de erratas (D-28) | ✅ |
| P-6 | Aplicar el filtro del Valle (167 estaciones) y el recorte del satélite en el flujo de la rama y rehacer las fases | Por hacer |
| P-7 | Conectar las coordenadas verificadas al cruce estación–píxel de la rama | Por hacer |
| P-8 | Confirmado: el satélite es CHIRP (D-37). Respuesta del profesor (13-nov-2026) | ✅ |
| P-10 | Ya no hay estaciones reservadas: cambió el plan del curso. La prueba es solo 2001–2003 con nuestras estaciones. Respuesta del profesor (13-nov-2026) | ✅ |
| P-11 | Ya no aplica: el CSV de predicciones era para las estaciones reservadas | ✅ |
| P-12 | ¿Se aceptan las fuentes externas: ONI de la NOAA, altitudes y coordenadas del CNE del IDEAM, polígono del DANE? | ⬜ |
| P-13 | Confirmado: solo estaciones dentro del Valle del Cauca. Respuesta del profesor (13-nov-2026) | ✅ |

---

### D-33 · Alto Anchicayá, altitudes y coordenadas a revisar ✅
- **Decisión:** `5301400108` (ALTO - ANCHICAYA) y `5301400105` son estaciones distintas. El punto del CNE `0053100050` (3.5333, -76.8667, 630 m) se asocia al nombre y a la altitud (650 m) de `5301400108`, no a `5301400105`. Por eso `5301400105` pasa a nivel 3 (departamento inferido por cuenca 530140, sin coordenadas). La coordenada del CSV de `5301400108` queda marcada "a revisar" y no se corrige sin confirmación.
- **Justificación:** r mensual 0.85, menos de 1 % de meses idénticos y 24 % más lluvia en `5301400105`, así que no es un duplicado. El terreno SRTM es 403 m en el punto del CSV, 643 m en el del CNE, y la altitud registrada es 650 m.
- **Cómo se validó:** SRTM 30 m frente a la altitud oficial en 169 estaciones: diferencia mediana de 46 m, P75 133 m, P90 387 m. Se marcan 22 con diferencia mayor de 300 m (puede ser error de SRTM en terreno quebrado o coordenada corrida) en la columna "Coordenada a revisar".
- **Altitud faltante:** el Excel del curso trae altitud solo en 163 de las 232 (corregido en D-34; antes se contaban 169 por una etiqueta mal puesta en las 6 del cruce con el CNE) estaciones del Valle con coordenadas. 50 se completaron con el CNE (misma estación) y 13 con SRTM, marcadas como estimadas en "Fuente de la altitud".
- **Regla de la pista de vecinas (95 %, 214 de 226):** valida solo el método de asignar departamento. No confirma la ubicación de ninguna estación concreta. Esas estaciones quedan en nivel 3.
- **Resultado:** Valle sin filtrar años = 267 estaciones (232 con coordenadas, 35 solo con departamento). Archivo: `estaciones_valle_del_cauca.xlsx`.
- **Fecha y responsable:** 2026-11-12, Gus y Claude.

### D-34 · Conjunto final de estaciones del Valle y método de la altitud ✅
- **Decisión:** el análisis usa **232 estaciones**: las 226 con ubicación oficial y las 6 identificadas por cruce con el CNE. Se descartan las 35 del Valle con departamento inferido (no hay coordenada exacta que rescatar) y las 21 sin ubicar. Las 56 del Cauca siguen excluidas (D-31). Todavía no se filtra por años: eso se hace en el EDA.
- **Por qué se incluyen las 6 del cruce:** el año de inicio de los datos coincide con el año de instalación del CNE (diferencia de 0 o 1 años), el tipo coincide, están a 3–9 km de la cuenca que indica el código y son estaciones únicas y conocidas (Palmira ICA, Aeropuerto Alfonso Bonilla Aragón, Central Castilla, Manuel M. Mallarino, Cent. Admo. La Unión, El Queremal). El CNE redondea al minuto (±1 km), muy por debajo del píxel satelital (~5 km). Van marcadas como "Cruce con CNE" para el análisis de sensibilidad (con y sin ellas).
- **Altitud (orden de prioridad):**
  1. Excel del curso (`CatalogoEstacionesElevacion.xlsx`): 163 estaciones.
  2. CNE del IDEAM, cuando es la misma estación: 56 (50 oficiales + las 6 del cruce). Validación: en las 139 estaciones que están en ambas fuentes, la altitud es igual (diferencia mediana 0 m).
  3. Terreno SRTM 30 m en la coordenada (OpenTopoData), marcada como estimada: 13. Validación frente a la altitud oficial en 169 estaciones: diferencia mediana 46 m, P75 133 m, P90 387 m.
- **Control de coordenadas:** si la altitud difiere más de 300 m del terreno SRTM, la coordenada queda "a revisar" (22 estaciones). No se corrige sin una segunda fuente. Cotejo con el estudio pluviométrico de Palmira 2020: sus coordenadas CVC coinciden a 40 m con las del curso, así que vienen de la misma base y no confirman nada de forma independiente. Para La Sirena (`2612810203`), ese estudio da un punto a 3,6 km del nuestro, pero ahí el terreno está a 3.473 m frente a los 2.605 m registrados. Nuestro punto (GPS del POT del Nima, terreno a 2.674 m) es el coherente y se mantiene.
- **Resultado:** `estaciones_valle_finales.xlsx` / `.csv` (232 filas). El inventario completo sigue en `estaciones_valle_del_cauca.xlsx`.
- **Fecha y responsable:** 2026-11-12, Gus y Claude.

### D-35 · La altitud solo sale de fuentes oficiales; SRTM queda solo como control ✅
- **Decisión:** no se estima ninguna altitud con SRTM. El orden es: (1) Excel del curso, 163 estaciones; (2) CNE del IDEAM, 65 estaciones. Las 4 sin altitud oficial (Univalle, Los Llaneros, La Rosita, La Esperanza) quedan vacías. SRTM se usa solo para controlar coordenadas: si la altitud oficial difiere más de 300 m del terreno, la coordenada se marca "a revisar".
- **Cambio frente a D-34:** 9 de las 13 que se habían estimado con SRTM sí están en el CNE, con su propio código. Se identificaron por nombre, municipio igual o vecino, año de instalación compatible y altitud coherente con el terreno: Los Cuervos 2.620, Las Veraneras 1.440, Lloreda 930, La Cristalina 1.490, Venecia 1.532, Yurumanguí 21, La Siberia 1.017, San Fernando 990 y Cajones 2.472 m. La Esperanza no se asigna porque las del CNE están a 26–32 km.
- **Cajones (`2612200105`):** en la coordenada del CSV el terreno está a 3.970 m; en la del CNE (3.226417, -76.140806), a 2.388 m, coherente con su altitud de 2.472 m. La coordenada del CSV probablemente está mal; queda marcada con la del CNE como alternativa. Ahora son 23 coordenadas a revisar.
- **Qué es SRTM:** un modelo de elevación del terreno medido por radar desde el transbordador espacial de la NASA (año 2000), con celdas de 30 m. Se consulta en OpenTopoData con la latitud y la longitud.
- **Resultado:** `estaciones_valle_finales.xlsx` / `.csv` (232 estaciones). Resuelve P-2.
- **Fecha y responsable:** 2026-11-12, Gus y Claude.

### D-36 · Reglas de depuración sacadas del EDA de las 232 ✅
- **Población:** de las 232 (D-34), **220 tienen datos en 1981–2003**; las 12 restantes solo tienen datos fuera del periodo. Las 232 caen dentro del polígono DANE del Valle (`lluvia_valle_datos.R`, sección 1).
- **Qué dicen los datos** (`lluvia_valle_datos.R`, secciones 2–6):
  - Un faltante es una fila ausente. No hay `NA`, negativos, duplicados ni el valor imposible. Hay 5 144 valores con decimales: vienen de estaciones climatológicas (tipo 02), que leen a 0,1 mm. No son error.
  - La cobertura es bimodal: el 81 % de los estación-años tiene ≥99 % de días. Los huecos internos son el 2,7 % de la rejilla, con mediana de 3 días.
  - Hay 5 días >3·p99 y ≥150 mm, y las vecinas a ≤25 km registran como máximo 78 mm: Pardo 600, Heraclio Uribe 503, La Cascada 421, Costa Rica 283 y Venecia 283.
  - Hay 19 meses con un solo día de lluvia y ese valor supera el p99. Son posibles totales mensuales anotados en un día.
  - De 918 rachas de ceros de ≥30 días, solo 6 tienen vecinas lloviendo (≥40 %): 5 de Loboguerrero (cañón seco del Dagua) y 1 de La Despensa. El contraste por vertiente de D-12 sobremarcaba.
  - No hay años copiados entre estaciones.
- **Experimento de huecos** (`lluvia_valle_datos.R`, sección 7):
  - Se borran días de años completos y se mide cuánto cambia la curva B-spline, en unidades de la variabilidad interanual de la estación.
  - Un hueco continuo de 15 días mueve la curva 0,12; uno de 30 días, 0,23; uno de 60 días, 0,48 (K=8).
  - Con máscaras reales, el hueco más largo explica el error mejor que el % de cobertura: R² 0,58 frente a 0,46.
- **Reglas propuestas:**
  1. Un año es válido si tiene ≥90 % de días y ningún hueco mayor de 15 días. Con 30 días, el p99 del error pasa de 0,64 a 0,95 y solo se ganan 2 estaciones (D-11 pedía ≤30).
  2. Una estación entra si tiene ≥7 años válidos en entrenamiento: quedan **176 estaciones**, entre ellas las 6 del CNE y 12 con coordenada a revisar.
  3. Los 5 extremos sin respaldo de vecinas pasan a faltante.
  4. Los totales mensuales en un día pasan a faltante si el mes está casi vacío; si el mes está completo, deciden las vecinas.
  5. Las rachas de ceros no se borran por regla. Las 6 se revisaron a mano:
     - **Loboguerrero, oct–nov de 1995, pasa a faltante.** Es temporada de lluvias: la estación suele llover el 41 % de esos días, y Dagua, a 12 km en el mismo cañón, tuvo 12 días con lluvia (83 mm).
     - **Las otras 5 se mantienen.** La Despensa (ene–feb 1996) y Loboguerrero (1982, 1982–83, 1989 y 1991) caen en temporada seca o en El Niño, y la vecina del mismo cañón también estuvo casi seca.
  6. No se imputa (D-13).
- **Precisión de las reglas 3 y 4 (implementación):**
  - Un extremo pasa a faltante si ninguna vecina a ≤25 km supera P/3 ese día, el anterior o el siguiente.
  - Un mes sospechoso pasa entero a faltante si tiene menos de la mitad de sus días observados o si las vecinas llovieron en ≥40 % de sus días. Si el total se anotó en un día, los ceros del resto del mes también son falsos.
- **Confirmada por Gus el 2026-11-12. Resultado** (`lluvia_valle_datos.R`, sección 8):
  - **Días pasados a faltante:** 100. Son 5 extremos, 62 días de 13 meses (de los 19 sospechosos; en los otros 6 las vecinas también estaban secas) y 33 días de la racha de Loboguerrero. Casi todos caen en años que ya eran inválidos por sus huecos.
  - **Años válidos:** 3 482 de 3 943 (2 361 en entrenamiento, 722 en validación y 399 en prueba).
  - **Estaciones incluidas: 176.** Son 170 con ubicación oficial y 6 del CNE; 12 tienen la coordenada "a revisar". De ellas, 147 tienen algún año válido en validación, 137 en prueba y 136 en ambos.
  - **Conjunto final:** `salidas_R_datos/lluvia_valle_limpia.rds`, con 3 326 estación-años (2 278, 673 y 375 por partición) y 1 213 990 días, de los cuales 2 887 faltan dentro de años válidos.
  - **El `.rds` trae:**
    - `$diario`: solo incluidas y años válidos;
    - `$diario_todo`: las 220 estaciones;
    - `$anios`, `$estaciones`, `$marcas` y `$reglas`.
  - **Copias en CSV:** `salidas_R_datos/tablas/8_marcas_limpieza.csv`, `8_anios_validez.csv` y `8_estaciones_incluidas.csv`.
  - **Todo el proceso está en un solo archivo, `lluvia_valle_datos.R`:** población, EDA, experimento de huecos, depuración y dataset final. Se corre desde la raíz del repositorio.
- **Fecha y responsable:** 2026-11-12, Gus y Claude.

### D-37 · Fuente del satélite entregado: CHIRP (solo infrarrojo), versión anterior al 10-may-2013 ✅ (confirmado por el profesor)
- **Conclusión:** los 2 328 TIFF (`ppt198101` a `ppt201324`) son un recorte de **CHIRP**, el producto solo-satélite del Climate Hazards Group (UCSB/USGS), en una versión anterior a la que hoy publica el CHC. **No son CHIRPS** (no tienen la corrección con estaciones). Es una inferencia a partir de la evidencia; falta que el profesor la confirme (P-8).
- **Evidencia en los archivos:**
  - Malla de 0,05°, 72 péntadas por año desde 1981 y la misma escala de píxel que los archivos del CHC (`0.0500000007450581`).
  - 2 326 archivos creados el **6-may-2013** en unas 4 horas; la última péntada es la del 26–30 de abril de 2013. Los dos primeros (`ppt198101` y `ppt198102`) se volvieron a guardar el 13–14 de junio de 2013; sus valores no se distinguen del resto.
  - Escritos con convenciones de ArcGIS: `GCS_WGS_1984`, bloques de 128×128, enteros de 16 bits y sin dato = −9999.
  - Valores enteros de 0 a 200 mm por péntada, con **tope en 200** (8 196 valores, el 0,24 %). Casi nunca hay ceros (7 en 3,36 millones) y **hay datos sobre el mar**.
- **Comparación número por número** (ventana del Valle leída del servidor del CHC, 7 años; `datos_web/chirp_v2_valle_muestra_7anios.bin`):

  | Contra | r por péntada y píxel | r del mapa medio | Total anual frente al entregado | ¿Cubre el mar? | Tope |
  |---|---|---|---|---|---|
  | CHIRP v2.0 (solo satélite) | 0,84–0,92 | 0,996 | Igual (±2 %) en 1985, 1990, 2000 y 2012 | Sí | 250 |
  | CHIRPS v2.0 (con estaciones) | 0,77–0,87 | 0,965 | El entregado es 10–26 % menor | No | — |

  - Los archivos de CHIRP v2.0 de 1981–2012 que hoy están en el CHC se generaron el **10–11 de mayo de 2013**, cuatro días después de los nuestros. Por eso no coinciden valor a valor: los nuestros son de la corrida anterior.
- **¿Fueron modificados?** No hay señal de alteración manual: no hay capas repetidas ni constantes, todos los archivos tienen el mismo formato y los totales anuales siguen a CHIRP. Lo que sí hubo es procesamiento: recorte a la ventana, redondeo a enteros y tope de 200 mm.
- **Los años malos (D-16) vienen del producto, no de una modificación:** CHIRP v2.0 del CHC también da 1985 y 2000 secos (2 967 y 3 221 mm frente a 2 907 y 3 160 del entregado) y 1988 todavía más seco (2 327 mm). Las estaciones y CHIRPS dicen que esos años fueron normales o húmedos. Es un error del insumo infrarrojo que CHIRPS corrige con estaciones.
- **Consecuencias:**
  - El satélite es independiente de las estaciones, que es lo que necesita el proyecto (corregir el sesgo del satélite con las estaciones). Refuerza no reemplazarlo por CHIRPS (P-3, D-18).
  - Hay que pasar de mm por péntada a mm/día (D-15) y tener presente el tope de 200 al interpretar extremos.
  - Los 90 píxeles de mar se quitan con el recorte al polígono del Valle (D-03).
- **Fecha y responsable:** 2026-11-13, Gus y Claude.

### D-38 · EDA del satélite: qué acierta, qué no y los "años defectuosos" ✅ (revisa D-16)
- **Dónde está:** secciones 9–12 de `lluvia_valle_datos.R`. Salida: `salidas_R_datos/satelite_valle.rds`.
- **Preparación:**
  - Se usan los 676 píxeles con centro dentro del polígono del Valle (D-03) y las 1 656 péntadas de 1981–2003.
  - Cada péntada se divide por sus días reales para pasar a mm/día (D-15). Sin dividir, la 6.ª péntada de los meses de 31 días parece más lluviosa (49 mm frente a 37).
  - Cada estación se empareja con el píxel que contiene su coordenada (P-7): las 176 caen dentro de la ventana, en 148 píxeles distintos; 23 píxeles tienen más de una estación.
- **Qué dicen los datos:**
  - **Aplana el espacio.** El satélite va de ~1 100 a ~6 200 mm/año y las estaciones de ~800 a ~12 100. Sobreestima donde llueve poco y subestima donde llueve mucho (r = −0,68 entre el sesgo y la lluvia de la estación).
  - **Sigue poco a la estación péntada a péntada:** r ≈ 0,5 en conjunto y 0,42 de mediana por estación. Casi nunca marca cero; las estaciones, entre el 19 y el 26 % de las péntadas.
  - **Acierta el clima medio, no las anomalías:**

    | Componente | r satélite–estación |
    |---|---|
    | Dónde llueve más (media por estación) | 0,79 |
    | Cuándo llueve en el año (ciclo medio) | 0,78 |
    | Meses más o menos lluviosos de lo normal | 0,30 |
    | Años más o menos lluviosos de lo normal | 0,25 |

  - **La diferencia con las estaciones sigue al ENSO** (r = 0,74 con el ONI). El satélite se queda corto en La Niña y se pasa en El Niño. Su variación de un año a otro es la mitad (7 % frente a 14 %), y su relación con el ONI tiene el signo contrario al de las estaciones (+0,27 frente a −0,60).
  - **Años que más se apartan:** 1988, 1984, 2000 y 1999 (todos La Niña; el satélite queda 19–26 puntos por debajo) y 1992 (El Niño; 29 puntos por encima).
- **Cambio frente a D-16:** los años "defectuosos" no son archivos dañados. 1988 y 2000 son casos de este efecto sistemático, 1985 no se aparta especialmente de las estaciones (el satélite lo da seco y las estaciones también, algo menos) y faltaban 1984, 1992 y 1999. No hay razón para excluir años del entrenamiento: el modelo tiene que aprender justo esta corrección.
- **Consecuencias para el modelo:**
  - La curva del satélite aporta sobre todo el ciclo anual y la ubicación. Para corregir años Niño o Niña hace falta información que el satélite no trae: por eso el enunciado agrega el ONI del año al vector z de la red.
  - P-4 pierde peso: 2000 no es un año dañado sino un año La Niña, y en validación sirve para medir ese efecto.
  - Santa Teresa (`2612810102`) tiene r = 0 con su píxel y la coordenada "a revisar": entra al análisis de sensibilidad.
- **Fecha y responsable:** 2026-11-13, Gus y Claude.

### D-39 · Datos preparados para la construcción funcional, sin fuga y con pruebas ✅
- **Qué queda listo** (`lluvia_valle_datos.R`, secciones 13 y 14; salida `salidas_R_datos/datos_preparados.rds`):
  - `$estaciones$Y`: 3 326 curvas estación–año × 365 días, en mm/día, con `NA` en los faltantes (0,24 %).
  - `$satelite$X`: 15 663 curvas píxel–año × 72 péntadas, en mm/día, con `$satelite$t` = día central de cada péntada. Son 681 píxeles: los 676 del Valle y 5 que contienen una estación aunque su centro quede fuera del polígono.
  - `$pares`: une cada curva de estación con la curva del satélite de su píxel y su año (2 278 en entrenamiento, 673 en validación y 375 en prueba).
    Incluye el ONI anual y la fase ENSO de cada año (NOAA CPC, D-19), que van al vector z junto con la altitud de `$meta_estaciones`.
  - `$meta_estaciones`: coordenada, altitud, vertiente, lluvia anual y marcas para el análisis de sensibilidad (6 del cruce con el CNE, 12 con coordenada a revisar).
- **Sin fuga de información (rúbrica 4.1):** el p99 de cada estación, que es el umbral de las reglas R2 y R3, ahora se calcula solo con 1981–1995. Antes usaba todo el registro. El resultado no cambió: siguen siendo los mismos 100 días marcados.
- **Trazabilidad (rúbrica 4.1):** `tablas/8_trazabilidad.csv`.

  | Paso | Estaciones | Pares estación–año |
  |---|---|---|
  | Crudo, 1981–2003 | 323 | 5 514 |
  | Solo el Valle, con coordenada | 220 | 3 943 |
  | Reglas R1 a R4 (solo cambian días) | 220 | 3 943 |
  | R5: años válidos | 214 | 3 482 |
  | R6: ≥7 años válidos en entrenamiento | 176 | 3 326 |

- **Pruebas:** 26 comprobaciones automáticas dentro del script; si una falla, el script se detiene. Entre ellas:
  - cada valor final es igual al del archivo crudo, releído de nuevo (1 211 103 valores);
  - la matriz del satélite coincide con 12 TIFF releídos;
  - cada estación cae dentro de su píxel;
  - el p99 usa solo entrenamiento.
- **Reproducibilidad:** dos corridas completas producen archivos `.rds` idénticos byte a byte.
- **Pendientes que salen del enunciado:**
  - El ONI y la altitud entran al modelo por enunciado (se agregan al vector z). Univalle (`2622320201`) no tenía altitud oficial: resuelto en D-40.
  - Falta decidir si se usan registros fuera del periodo (1930–1979) para estimar el FPCA.
- **Control de la rúbrica:** `CONTROL_RUBRICA.md`.
- **Fecha y responsable:** 2026-11-13, Gus y Claude.

### D-40 · Altitud de Univalle: 996 m, de la estación del IDEAM en el mismo campus ✅ (excepción a D-35)
- **Problema:** la altitud entra al vector z de la red (enunciado), y Univalle (`2622320201`, 39 años de datos, en Cali) no tiene altitud ni en el Excel del curso ni en el CNE con su código.
- **Decisión:** usar 996 m, la altitud oficial de la estación del IDEAM "UNIVERSIDAD DEL VALLE – AUT" (CNE `0026055120`, 3.378, −76.534), que está en el mismo campus, a 1,6 km. Esa estación abrió en 2006 y no aporta lluvia; solo se toma su altitud.
- **Justificación:** sigue siendo una fuente oficial (D-35), no una estimación con SRTM. La alternativa, sacar a Univalle del modelo, perdería una estación con registro largo en Cali.
- **Cómo se validó:**
  - El SRTM en la coordenada de la CVC da 977 m: 19 m de diferencia. En las 28 estaciones a menos de 15 km, el SRTM difiere de la altitud oficial 15 m en mediana, así que la diferencia es normal.
  - Las vecinas son coherentes: Edificio CVC, a 4,8 km, está a 985 m.
  - Prueba automática nueva: todas las estaciones incluidas tienen altitud.
- **Dónde:** sección 1 de `lluvia_valle_datos.R` (columna `fuente_alt`). El CSV de `datos_web/` no se modifica.
- **Fecha y responsable:** 2026-11-13, Gus y Claude.

## Fuentes

- Enunciado del proyecto: `Proyecto_Final_Analitica.docx`.
- Archivos del curso: `CVC_Wheather_Stations_Long_Lat.csv`, `CatalogoEstacionesElevacion.xlsx`, `P_1930_1979.xlsx`, `P_1980_1990.csv`, `P_1991_2003.csv`, `ImagenesSatelitePentadal/`.
- Catálogo de estaciones hidroclimatológicas de la CVC (datos.gov.co, kfsz-ija2): https://www.datos.gov.co/Ambiente-y-Desarrollo-Sostenible/Cat-logo-de-estaciones-hidroclimatol-gicas-de-la-CVC/kfsz-ija2 — copia local: `datos_web/cvc_catalogo_datosgov_kfsz-ija2.csv`.
- Catálogo Nacional de Estaciones del IDEAM (datos.gov.co, hp9r-jxuu): https://www.datos.gov.co/Ambiente-y-Desarrollo-Sostenible/Cat-logo-Nacional-de-Estaciones-del-IDEAM/hp9r-jxuu — copia local: `datos_web/CNE_IDEAM_dominio_satelite.csv`.
- Plan de manejo PNR del Nima, POT de Palmira (GPS de La Sirena, págs. 44–45): https://www.palmira.gov.co/planeacion/wp-content/uploads/pot/potmoderno/PM%20PNR%20DEL%20%20NIMA.pdf
- Estudio pluviométrico de Palmira 2020: https://palmira.gov.co/wp-content/uploads/2021/11/Estudio_Pluviometrico_Palmira_V07.pdf
- Elevación SRTM 30 m vía OpenTopoData: https://api.opentopodata.org (ejemplo: `/v1/srtm30m?locations=3.52242,-76.104246`).
- Sistemas de coordenadas: https://epsg.io/3115 y https://epsg.io/9377
- ONI (NOAA CPC): https://origin.cpc.ncep.noaa.gov/products/analysis_monitoring/ensostuff/ONI_v5.php
- CHIRPS v2 (UCSB CHC) y NASA POWER (MERRA-2): ver `02_datos_externos.ipynb`.

---

## Plantilla para nuevas decisiones

```
### D-xx · Título corto [✅ | 🟡 | 🔁]
- **Decisión:**
- **Justificación:**
- **Cómo se validó:**
- **Resultado / archivo:**
- **Fecha y responsable:**
```

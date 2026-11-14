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
| P-4 | Resuelta: **no se eliminan años** del satélite (1985, 1988, 2000 incluidos); el modelo debe aprender a corregir también esos casos. Revisa D-37 y D-42 (ver D-47). Respuesta del profesor (13-nov-2026) | ✅ |
| P-5 | Resuelta por la fe de erratas (D-28) | ✅ |
| P-6 | Aplicar el filtro del Valle y el recorte del satélite y rehacer las fases | Resuelta (D-34, D-37 a D-40) |
| P-7 | Conectar las coordenadas verificadas al cruce estación–píxel | Resuelta (D-38; confirmado en D-44) |
| P-8 | Confirmado: el satélite es CHIRP (D-43). Respuesta del profesor (13-nov-2026) | ✅ |
| P-9 | Resuelta por el enunciado: el ONI del año y la altitud se agregan al vector z de la red | ✅ |
| P-10 | Ya no hay estaciones reservadas: cambió el plan del curso. La prueba es solo 2001–2003 con nuestras estaciones. Respuesta del profesor (13-nov-2026) | ✅ |
| P-11 | Ya no aplica: el CSV de predicciones era para las estaciones reservadas | ✅ |
| P-12 | Resuelta: se aceptan las fuentes externas oficiales, principalmente del IDEAM. Respuesta del profesor (13-nov-2026) | ✅ |
| P-13 | Confirmado: solo estaciones dentro del Valle del Cauca. Respuesta del profesor (13-nov-2026) | ✅ |
| P-14 | Resuelta: se acepta la altitud de Univalle (996 m) tomada del IDEAM, porque los datos del curso también vienen del IDEAM (D-46). Respuesta del profesor (13-nov-2026) | ✅ |
| P-15 | Respondida: hay que armar la serie y dibujarla para entender el tope. Hecho en la sección 12 (`12_serie_tope_200.png`): el tope afecta 33 de 238 224 péntadas comparables y solo 3 estaciones del Pacífico (D-47). Respuesta del profesor (13-nov-2026) | ✅ |
| P-16 | ¿Qué "ONI del año" espera: la media de los 12 trimestres del año (lo que se usa) o un trimestre concreto (p. ej. DEF)? | ⬜ |
| P-17 | Resuelta: está bien quitar el 29 de febrero. Respuesta del profesor (13-nov-2026) | ✅ |
| P-18 | Resuelta: el píxel que contiene la estación está bien; como el píxel es más grande que el punto, no hay que reescalar ni promediar. Respuesta del profesor (13-nov-2026) | ✅ |
| P-19 | ¿Se pueden usar los registros de 1930–1979 solo para estimar el FPCA, aunque la fe de erratas fije 1981–2003? | ⬜ |
| P-20 | Resuelta (D-51): K = 8 en las estaciones. La ventaja de K = 16 en la fecha de los picos dependía de la referencia elegida | ✅ |

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
  5. *(Reemplazado por D-49: R4 ahora es una regla calculada.)* Las rachas de ceros no se borran por regla. Las 6 se revisaron a mano:
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

### D-37 · Satélite recortado al Valle, en mm/día y con los años defectuosos marcados (cierra D-03 y aplica D-15/D-16)
- **Decisión:**
  - Se usan los píxeles con centro dentro del polígono DANE del Valle (676 de 1 443), más los 5 píxeles de estaciones cercanas al borde: 681 píxeles × 23 años = 15 663 curvas.
  - Cada péntada se divide por sus días reales (5, o de 3 a 6 en la última del mes) y se ubica en su día central del año de 365 días.
  - 1985, 1988 y 2000 se marcan como defectuosos. Se conservan en el EDA y se excluyen de la elección de K, de la comparación principal y del modelo.
- **Qué mostró el EDA** (`funcional/10_satelite_emparejamiento.R`):
  - Esos tres años son los más secos del satélite en el Valle: 2 219, 2 302 y 2 374 mm, frente a ~2 730 de mediana.
  - En los pares estación–píxel la diferencia se invierte en esos años: +336 mm/año a favor de la estación, frente a −83 en los demás.
  - **Los tres son años Niña.** Con ellos, el satélite no muestra la señal de La Niña; sin ellos, la recupera: Niña 1 949 > neutro 1 892 > Niño 1 875 mm (sección 12).
- **Tope de 200 mm/péntada:** 7,7 % de las curvas lo tocan alguna vez, todas en la vertiente Pacífica (al oeste de −76,8°).

### D-38 · Emparejamiento estación–píxel: la celda que contiene la estación (cierra P-7)
- **Regla:** cada estación se asigna a la celda del TIFF que contiene su coordenada de `estaciones_valle_finales.csv`.
- **Resultado:** las 176 estaciones caen en 148 píxeles. 23 píxeles tienen de 2 a 4 estaciones (51 estaciones en total, 16 de ellas en Cali), que comparten la misma curva satelital. El centro del píxel queda a 2,2 km en mediana (máximo 3,9 km).
- **Tablas:** `salidas_R_funcional/tablas/10_emparejamiento_estacion_pixel.csv` y `10_pixeles_compartidos.csv`.

### D-39 · K = 8 (estaciones) y K = 11 (satélite), con las 176 del Valle (confirma D-29)
- **Cómo se eligió:** ECM, VC y VCG con las curvas completas de entrenamiento: 1 899 de estación y 8 853 de satélite, sin 1985 ni 1988.
- **Resultado:**
  - VC y VCG coinciden en ambas fuentes.
  - En estaciones la VC es plana entre K = 5 y 17.
  - S es idempotente y tr(S) = K hasta 1e-14.
- **Valores negativos:** 47 % de las curvas de estación baja de 0 en algún día, casi siempre poco (mediana del mínimo −0,07 mm/día). Ocurre en julio–agosto y en los bordes de la base no periódica, y se resuelve recortando en 0 antes de la raíz (D-23).

### D-40 · Regresión funcional de referencia: FPCA con M = 5 y λ = 1e8, elegidos con validación
- **Modelo:** en la escala raíz,

  y(t) = μ(t) + α(t) + ∫β(s,t)x(s)ds + γ(t)·altitud + δ(t)·ONI,

  con β(s,t) = Σ_m β_m(s)ψ_m(t).
  - Las ψ_m y μ salen del FPCA de las estaciones, solo con entrenamiento.
  - Cada β_m(s) se representa en la base del satélite, con penalización ∫β''².
- **Elección de M y λ:** se eligen juntos con validación (1996–1999). El error deja de bajar desde M = 5, que explica el 91 % de la varianza; las componentes 6 a 8 tienen autovalores casi iguales, es decir, ruido.
- **Verificación:** la predicción por la integral de β(s,t) y la predicción por puntajes FPCA difieren en 5e-15. La cabeza lineal de la red debe reproducir esta identidad.
- **Prueba (2001–2003), frente al satélite sin corregir:**

  | Métrica | Satélite | Regresión |
  |---|---|---|
  | RMSE L2 (mm/día) | 2,43 | 2,24 |
  | KGE | 0,38 | 0,45 |
  | Error en el pico de la 1.ª temporada (días) | 23 | 18 |

  - La regresión mejora 60 % de las curvas, pero los intervalos se solapan.
  - Sin el satélite (solo altitud y ONI), el RMSE es 2,61.
  - En 2002 (Niño) el modelo corrige de más.
- **Univalle:** antes salía por no tener altitud oficial (17 pares); desde D-46 entra con 996 m. Cifras actualizadas al volver a correr (13-nov-2026): 375 curvas de prueba de 137 estaciones.

### D-41 · Organización del código funcional y paso a un solo archivo R (en curso)
- **Organización mientras se desarrolla:** los módulos están en `funcional/`:
  - `00_funciones.R` con las funciones comunes;
  - las secciones 10 a 14, cada una con su script.
  - Cada módulo lee el `.rds` del anterior desde `salidas_R_funcional/`, para no tener que correr todo cada vez.
- **Para la entrega:**
  - Se concatenan, después de `lluvia_valle_datos.R`, el contenido de `00_funciones.R` y luego las secciones 10 a 14, sin las líneas `source(...)`.
  - Los `readRDS` intermedios pueden quedarse, porque leen lo que el mismo archivo acaba de guardar.
  - Así queda un único `.R`, y el `.py` solo tiene la red (lee `salidas_R_funcional/para_python/`).
- **Reemplaza:** `01_datos_y_construccion_funcional.R` y `02_EDA_funcional.R`. Su código (Cox–de Boor, MCO, criterios de K, controles de S, momentos y profundidades) quedó en las secciones 11 y 12, con el satélite agregado.
- **Fecha y responsable:** 2026-11-13, Jacobo.

### D-42 · Sin estaciones reservadas, años defectuosos fuera y un `.py` autosuficiente
- **Estaciones reservadas:** no habrá. La prueba es 2001–2003 con las mismas estaciones; el informe no evalúa la generalización a estaciones nuevas.
- **Años defectuosos del satélite:** 1985, 1988 y 2000 siguen fuera de la elección de K, la comparación principal y el modelo, mientras se revisan (P-4).
- **Reparto R / Python:** el `.R` tiene la depuración, el EDA, la construcción funcional, el EDA funcional, la comparación y la regresión de referencia. El `.py` debe rehacer lo que necesita para la red: base B-spline, MCO, W, R, escala raíz y FPCA, a partir de las observaciones exportadas por R (lluvia diaria de la estación y péntadas del píxel de cada par). Los coeficientes calculados en R sirven solo para verificar que Python da lo mismo. Esto modifica D-41.
- **Fecha y responsable:** 2026-11-13, Jacobo.

> **Nota de la unión (13-nov-2026):** Jacobo (secciones funcionales) y Gus (datos) trabajaron en paralelo y ambos numeraron desde D-37. Se conservan D-37 a D-42 de Jacobo y las de Gus pasan a D-43 a D-46.

### D-43 · Fuente del satélite entregado: CHIRP (solo infrarrojo), versión anterior al 10-may-2013 ✅ (confirmado por el profesor)
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

### D-44 · EDA del satélite: qué acierta, qué no y los "años defectuosos" ✅ (revisa D-16; choca con D-37 y D-42, ver P-4)
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

### D-45 · Datos preparados para la construcción funcional, sin fuga y con pruebas ✅
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
  - El ONI y la altitud entran al modelo por enunciado (se agregan al vector z). Univalle (`2622320201`) no tenía altitud oficial: resuelto en D-46.
  - Falta decidir si se usan registros fuera del periodo (1930–1979) para estimar el FPCA.
- **Control de la rúbrica:** `CONTROL_RUBRICA.md`.
- **Fecha y responsable:** 2026-11-13, Gus y Claude.

### D-46 · Altitud de Univalle: 996 m, de la estación del IDEAM en el mismo campus ✅ (excepción a D-35; revisa D-40, donde Univalle quedaba fuera)
- **Problema:** la altitud entra al vector z de la red (enunciado), y Univalle (`2622320201`, 39 años de datos, en Cali) no tiene altitud ni en el Excel del curso ni en el CNE con su código.
- **Decisión:** usar 996 m, la altitud oficial de la estación del IDEAM "UNIVERSIDAD DEL VALLE – AUT" (CNE `0026055120`, 3.378, −76.534), que está en el mismo campus, a 1,6 km. Esa estación abrió en 2006 y no aporta lluvia; solo se toma su altitud.
- **Justificación:** sigue siendo una fuente oficial (D-35), no una estimación con SRTM. La alternativa, sacar a Univalle del modelo, perdería una estación con registro largo en Cali.
- **Cómo se validó:**
  - El SRTM en la coordenada de la CVC da 977 m: 19 m de diferencia. En las 28 estaciones a menos de 15 km, el SRTM difiere de la altitud oficial 15 m en mediana, así que la diferencia es normal.
  - Las vecinas son coherentes: Edificio CVC, a 4,8 km, está a 985 m.
  - Prueba automática nueva: todas las estaciones incluidas tienen altitud.
- **Dónde:** sección 1 de `lluvia_valle_datos.R` (columna `fuente_alt`). El CSV de `datos_web/` no se modifica.
- **Fecha y responsable:** 2026-11-13, Gus y Claude.

### D-47 · Respuestas de la asesoría del 13-nov-2026 ✅
- **No se eliminan años (P-4).** El profesor pidió no quitar 1985, 1988 ni 2000: el modelo tiene que aprender a corregir también los años Niño y Niña. Esto **revisa D-37 y D-42**. En la parte funcional, la constante `ANIOS_SAT_DEFECTUOSOS` debe dejar de excluir pares y quedar solo para un análisis de sensibilidad. La auditoría (`AUDITORIA_FUNCIONAL.md`, G1) ya había mostrado que incluirlos mejora la prueba.
- **29 de febrero (P-17), píxel de la estación (P-18), Univalle (P-14) y fuentes oficiales del IDEAM (P-12):** aceptados como están.
- **Tope de 200 mm (P-15).** Se armó la serie péntada a péntada (sección 12, `12_serie_tope_200.png` y `tablas/12_tope_200_por_estacion.csv`):
  - En los pares estación–píxel, el satélite toca el tope en solo 33 de 238 224 péntadas, en 3 estaciones del Pacífico (Bahía Málaga 17, Yurumanguí 9, San Francisco 7).
  - La estación pasa de 200 mm en 1 049 péntadas, y en solo el 2 % de ellas el satélite está en el tope.
  - **Conclusión:** el tope no es la causa principal de que el satélite subestime los aguaceros. El satélite casi nunca llega a esos valores por sí mismo. El tope se documenta y no se corrige.
- **Base de Fourier o B-splines.** Se revisó el material de clase:
  - El enunciado exige B-splines en el codificador (X(t) y w_k(t)).
  - El ejemplo de clase con esta misma lluvia (`Script_AnalisisExploratorio.R`, `X_sat` y `Y_est`) usa B-splines cúbicas.
  - Fourier se usa en clase (`Script_Regresion_Funcional.R`) para `CanadianWeather`, que son promedios climatológicos de muchos años y, por eso, de verdad periódicos.
  - Una curva de un solo año no es periódica: el 31 de diciembre de un año no tiene por qué empalmar con el 1 de enero del mismo año.
  - **Se mantienen las B-splines**, y el efecto de borde se documenta.
- **Principio para la sustentación:** toda regla debe estar respaldada por resultados en el código ("¿dónde lo verifican?"). Cada regla de D-36 tiene su tabla de evidencia y una prueba automática. La revisión manual de R4 (Loboguerrero 1995) está en la tabla `5b_revision_rachas.csv`.
- **Fecha y responsable:** 2026-11-13, Gus y Claude.

### D-48 · Cierre de la parte funcional: un solo archivo, todos los años y cada decisión verificada ✅ (revisa D-37, D-39, D-40, D-41 y D-42)
- **Un solo archivo funcional.** `lluvia_valle_funcional.R` reemplaza a `funcional/00` y `10` a `14` (quedan en el historial de git). Tiene cinco secciones: datos y emparejamiento, construcción, EDA funcional, comparación y regresión. Quien cambie algo de lo funcional lo cambia ahí.
  - **Entrada única:** `salidas_R_datos/datos_preparados.rds`. Ya no se vuelve a leer el satélite ni a emparejar: se quitó el procesamiento duplicado de la antigua sección 10, y las 26 pruebas del archivo de datos cubren todo lo que entra.
  - **Son dos archivos R mientras se trabaja** (`lluvia_valle_datos.R` y `lluvia_valle_funcional.R`), para no mezclar. La entrega pide un único R: se pegan uno tras otro al final, porque el segundo solo depende del `.rds` que deja el primero.
- **Todos los años entran (D-47).** Se quitó la exclusión de 1985, 1988 y 2000 en la elección de K, la comparación y la regresión. Verificación en el código:
  - el cociente satélite/estación de cada año sigue al ONI (r = 0,68); 1984 y 1999 se apartan tanto como los tres años marcados;
  - sin esos años la prueba empeora (L² relativo 0,439 frente a 0,434) en 79 de 137 estaciones (Wilcoxon p = 0,003).
- **Métricas de la rúbrica 4.4.** Error L² relativo, KGE, sesgo del acumulado anual en % y error en la fecha de los dos picos; IC del 95 % por bootstrap de estaciones (2 000 remuestreos); Wilcoxon pareado por estación; resultados por fase ENSO.
  - Prueba 2001–2003, regresión frente al satélite sin corregir: L² relativo 0,51 → 0,44; KGE 0,38 → 0,47; sesgo anual +26 % → +7 %; pico de la primera temporada 23 → 18 días. El IC de la diferencia excluye el cero en cinco métricas; en el pico de la segunda temporada la ganancia (medio día) no es significativa.
  - La prueba no tiene ningún año Niña. Por eso se reporta también la validación por fase, marcada como orientativa. En los años Niña 1999–2000 la regresión casi empata con el satélite y su sesgo empeora (−0,5 % → −11 %): es el punto débil de la referencia lineal.
- **Vertiente de los píxeles: la de la estación más cercana.** El corte en la longitud −76,8° solo acierta el 88,6 % al comprobarlo con las estaciones (ningún corte pasa de 89,2 %); la estación más cercana acierta el 97,7 % dejando una fuera. Queda una sola definición para las dos fuentes: la cuenca de la CVC.
- **K = 8 y K = 11 se mantienen (D-39),** ahora elegidos con todos los años de entrenamiento. Sensibilidad: con K = 12, 16 o 24 la fecha mediana de los picos se mueve a lo sumo 9 días y los valores negativos suben de 2 % a 5 %.
- **Comparación con el mismo K.** Con el satélite también en K = 8, la diferencia media, la diagonal de la correlación y el cociente de varianzas casi no cambian. Las diferencias entre fuentes no salen de usar K distinto.
- **Primero se suaviza en mm y luego se saca la raíz (D-23 se mantiene).** Ese orden conserva el acumulado anual (cociente 1,000); sacar la raíz al dato diario y luego suavizar lo reduce a un tercio (0,33).
- **EDA por vertiente y con anomalías.** El 74 % de la varianza de las estaciones es diferencia entre sitios. Con anomalías por estación, la primera componente baja de 76 % a 29 % y abril y octubre casi no se relacionan (0,14).
- **Regresión (revisa D-40):** M = 6 (94,2 % de la varianza, con los datos de D-49) y λ = 3,2e8, elegidos con validación con todos los años. La equivalencia β(s,t) frente a puntajes FPCA da 4e-15.
- **Pruebas:** 12 comprobaciones automáticas al final del archivo (`tablas/F6_pruebas.csv`); si una falla, el script se detiene.
- **Salidas:** `salidas_R_funcional/` con los prefijos `F1_` a `F6_`, `para_python/` y `resultados_funcional.rds`. `para_python/pares.csv` ya no lleva la columna `defectuoso`.
- **Fecha y responsable:** 2026-11-13, Gus y Claude.

### D-49 · R4 calculada con los datos: rachas de ceros sin respaldo de la vecina ✅ (revisa la R4 de D-36)
- **Problema:** la R4 anterior era una revisión a mano. Se revisaron las 6 rachas en las que las vecinas llovían más del 40 % de los días, y solo una (Loboguerrero, oct–nov de 1995) quedó escrita con sus fechas en el código. No había un criterio verificable ("¿dónde lo verifican?").
- **Regla nueva** (sección 5b de `lluvia_valle_datos.R`, tabla `5b_rachas_ceros_prueba.csv` y gráfica `5b_rachas_ceros_prueba.png`). Para cada racha de 30 días o más en 0:
  - se toma la vecina más cercana (a 25 km o menos) que tenga datos en al menos el 80 % de la racha, y se cuenta en cuántos días k llovió;
  - con los días de entrenamiento (1981–1995) en que ambas tienen dato, se calcula q = P(llueve en la estación | llueve en la vecina). No usa validación ni prueba;
  - p = (1 − q)^k es la probabilidad de que la estación marque 0 todos esos días por azar;
  - la racha pasa a faltante si p < 0,05 / 917 = 5,5e-5. Es la corrección de Bonferroni, porque se prueban 917 rachas a la vez.
- **Por qué este umbral:** p supone días independientes y la lluvia viene en rachas, así que p exagera lo raro. Un umbral exigente compensa eso. Con 0,001 se marcarían 72 rachas, con más riesgo de borrar sequías reales.
- **Resultado:** 30 de 917 rachas son falsas, en 27 estaciones (1 356 días). La mayoría de las rachas largas son sequías reales: la vecina también estuvo seca (El Niño de 1997, temporadas secas).
  - El caso más claro es Bella Vista 1982: 49 días en 0 mientras La Ceja (8,7 km) llovió 23 días (p = 7e-10).
  - **Loboguerrero oct–nov de 1995 vuelve a ser dato.** Dagua llovió 12 días, pero coincide con Loboguerrero solo la mitad de los días de lluvia (q = 0,50); p = 2e-4 no alcanza el umbral. Un enclave seco con 33 días sin lluvia en noviembre es posible.
  - **Loboguerrero dic-1982 a mar-1983 (67 días) sí pasa a faltante.** La revisión a mano la había dejado porque "Dagua también estuvo seca", pero Dagua solo tiene 4 días de datos en esas fechas. La regla usa la siguiente vecina con datos, Limones (16,7 km), que llovió 39 días (p = 5e-7).
- **Efecto en el conjunto:** cada racha deja un hueco de más de 15 días, así que su año deja de ser válido (R5). Quedan 176 estaciones y **3 299 pares estación–año** (2 257 de entrenamiento, 672 de validación y 370 de prueba), frente a 3 326. Las 26 pruebas de datos y las 12 funcionales pasan.
- **Efecto en lo funcional:** casi ninguno. K sigue en 8 y 11. M pasa de 5 a 6 (94,2 % de la varianza, dentro del rango 4–6 del enunciado). En prueba la regresión baja el error L² relativo de 0,51 a 0,44 frente al satélite.
- **Fecha y responsable:** 2026-11-14, Gus y Claude.

### D-50 · R0: valores físicamente imposibles a faltante, con umbral citado ✅ (completa D-36)
- **Problema:** la sección 2 contaba los negativos y el máximo, pero ninguna regla los convertía en faltante. La ausencia de valores imposibles se suponía, no se verificaba ("¿dónde lo verifican?").
- **Regla R0** (sección 2 de `lluvia_valle_datos.R`, tabla `2_R0_valores_imposibles.csv`): un valor diario P < 0 o P > 1 825 mm pasa a faltante, antes de R1. 1 825 mm es el récord mundial de lluvia en 24 h reconocido por la OMM (Foc-Foc, La Reunión, 7–8 de enero de 1966).
- **Resultado:** en las 220 estaciones del Valle R0 no quita nada (0 negativos; máximo 600 mm, Pardo). En el crudo completo detecta un solo valor, 859 375 mm, de una estación fuera del Valle. La trazabilidad lleva la fila "1b. R0" con 0 días que salen, y una prueba nueva de la sección 14 lo comprueba en la matriz final.
- **Efecto en el conjunto:** ninguno; 176 estaciones y 3 299 pares, como antes.
- **Fecha y responsable:** 2026-11-14, Gus y Claude.

### D-51 · K = 8 en las estaciones: la fecha de los picos no sirve para elegir K ✅ (cierra P-20; confirma D-39)
- **Problema:** la VC de las estaciones es casi plana entre K = 5 y 17. Una prueba del stash "fase A" mostraba que con K = 16 la fecha de los picos queda más cerca de una media móvil de 31 días de los datos crudos (6 y 5 días por curva frente a 12 y 10 con K = 8), con la VC solo 0,9 % peor.
- **Prueba** (sección 2 de `lluvia_valle_funcional.R`; tablas `F2_K_picos_vs_referencia.csv`, `F2_K_que_gana_por_referencia.csv` y `F2_diferencia_entre_referencias.csv`): 1 643 curvas completas interandinas de entrenamiento, con K = 8, 10, 12, 16 y 20, y tres referencias que no dependen de K (medias móviles de 15, 31 y 61 días).
- **Resultado:** el K que gana cambia con la referencia: K = 20 con 15 y 31 días, K = 10 con 61 días. Con 61 días, K = 16 queda peor que K = 8 (13 y 15 días frente a 10 y 9). Las referencias difieren entre sí 8 a 19 días en la fecha del pico, lo mismo que se gana subiendo K, y la fecha varía 32 a 36 días entre curvas (RIC). La prueba mide qué suavizado se parece a la referencia, no cuál ubica mejor el pico.
- **Decisión:** se mantiene K = 8, el mínimo de la VCG (criterio de clase). El script se detiene si algún día una misma K ganara con las tres referencias, porque entonces este argumento ya no valdría.
- **Consecuencias:** no cambia nada del análisis. El experimento de huecos de R5 sigue valiendo, porque se hizo con K = 8. K = 10 queda documentado como alternativa (VC 0,34 % peor, menos negativos), sin adoptarlo.
- **Fecha y responsable:** 2026-11-14, Gus y Claude.

### D-52 · Red neuronal funcional y ablación de la cabeza en Python, conectada a lo de R ✅ (rúbrica 4.2 y 4.3)
- **Por qué Python:** el enunciado permite PyTorch y pykan solo en este apartado, y son de Python. R ya tiene la versión lineal (regresión funcional, sección 5), así que Python lee lo que deja R en `salidas_R_funcional/para_python/` y no vuelve a construir nada.
- **Base:** se reutilizó `archivo_flujo_anterior/05_red_funcional.py` (KAN con la forma del enunciado, MLP, conteo de parámetros, semillas, prueba sintética), reconectado a los insumos actuales en un solo archivo: `lluvia_valle_red.py`. Salidas en `salidas_python/tablas` (prefijos P1 a P6).
- **Puente con R** (sección 1 y 3): W por cuadratura igual a la de R (8e-14); FPCA de entrenamiento igual (1e-14, salvo signo); las métricas de prueba de la regresión de R se reproducen; la red lineal cargada con los coeficientes de R da sus predicciones (3,6e-15). El decodificador usa el μ y las ψ de R (M = 6, 94,2 %), congelados.
- **Arquitectura:** z_k = a_kᵀ W_s c con 8 neuronas, más altitud estandarizada y ONI. z se estandariza con su media y DE en ENTRENAMIENTO (afín, no rompe la equivalencia lineal; deja las entradas de la KAN dentro de su rejilla). Cabezas: lineal (10→6), KAN [10,4,6] con G = 5 (640 parámetros) y MLP [10,37,6] (635; 0,8 % de diferencia).
- **Cambios frente al script viejo, con su razón:**
  - *Penalización de rugosidad relativa* Σ ∫w_k''² / ∫w_k²: con z estandarizado, la penalización absoluta a_kᵀRa_k se podía bajar encogiendo las w_k sin suavizarlas.
  - *Optimizador:* con lr = 3e-3 y paciencia 300 la cabeza lineal paraba en una pérdida de validación de 142, peor que mínimos cuadrados (137). Con lr = 1e-2, hasta 5 000 épocas y paciencia 1 000 llega al mínimo de entrenamiento de MCO, y con λ_w = 1e8 a 126,3 (la regresión de R: 126,9). Sin esto, la ventaja de las otras cabezas habría sido un artefacto.
  - *Rejilla de λ_w* {0, 1e6, 1e7, 1e8, 1e9}: la rugosidad relativa es del orden de 1e-6; la rejilla vieja (≤ 0,1) no penalizaba nada.
  - *M y ψ de R* en lugar de elegir M otra vez con "≥ 95 % de varianza": un solo criterio en todo el proyecto (validación, D-48).
- **Presupuesto igual para las tres cabezas:** Adam, lr 1e-2, lote completo, hasta 5 000 épocas, parada temprana con la pérdida L² de validación (paciencia 1 000), 5 pruebas de λ_w con la semilla 2026 y luego 3 semillas (2026-2028). `P5_presupuesto.csv`.
- **Verificaciones** (14 de 14 pasan, `P6_pruebas.csv`): equivalencia de la cabeza lineal con β(s,t) a 1,5e-14 en las 3 semillas (criterio 1e-4); KAN sintética con RMSE 2e-3 (criterio 1e-2); parámetros igualados; mismo codificador y decodificador; reconstrucción no negativa (0 % recortado).
- **Resultado** (prueba 2001-2003, 370 curvas, 137 estaciones): L² relativo satélite 0,513, lineal 0,440, KAN 0,434, MLP 0,420. El MLP supera a la lineal en L² y KGE por poco (IC 95 % pareado que apenas excluye el 0); la KAN no se distingue de la lineal en L² ni en KGE y reduce el |sesgo| 4 puntos; KAN y MLP no se distinguen salvo en el pico 2. En validación la ventaja no lineal es mayor (pérdida 100-108 frente a 126), pero las no lineales sobreajustan pronto (época ~100). Detalle en el bloque "Qué dicen los datos" de la sección 5.
- **Falta:** interpretación 4.5 (w_k(t), activaciones de la KAN, ajuste simbólico, peores casos) y las métricas por fase ENSO para las tres cabezas (4.4).
- **Fecha y responsable:** 2026-11-14, Gus y Claude.

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
- Récord mundial de lluvia en 24 h (R0, D-50): tabla de récords de la OMM, https://public.wmo.int/sites/default/files/2025-07/Table_Records_25Jul2025.pdf, y archivo de extremos https://wmo.asu.edu
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

# Antecedentes: cómo se ha modelado la lluvia del Valle del Cauca

En la revisión por Teams, el profesor preguntó: *"en un paper de modelación de lluvia diaria, ¿cuáles son los criterios de inclusión?"*. Este documento revisa los trabajos que usan datos parecidos a los nuestros (estaciones del IDEAM y la CVC, satélite de la familia CHIRP/CHIRPS) en el Valle del Cauca o junto a él. Para cada uno registra:

- qué datos usa;
- cómo trata los faltantes;
- qué modela;
- en qué se parece y en qué se diferencia de lo nuestro.

Búsqueda hecha el 19-nov-2026. La columna "Leído" dice si se leyó el texto completo o solo el resumen.

---

## 1. Resumen

| Trabajo | Zona | Datos | Criterio de inclusión y faltantes | Qué modela | Leído |
|---|---|---|---|---|---|
| **Ospina-Galindez, Giraldo y Andrade-Bejarano (2019)**, *J. Applied Statistics* 46(8). Base: tesis de Ospina (Univalle, 2016) | Valle del Cauca | 82 estaciones (CVC, 700–2 000 m) y CHIRP **pentadal**; 1993–2011 en la tesis | Tesis: los faltantes se **imputan** con kriging ordinario, año por año | Regresión **concurrente** estación ~ satélite, con errores correlacionados en el espacio (traza-variograma y mínimos cuadrados generalizados). En la tesis la lluvia es **acumulada** | Tesis completa; del artículo, solo el resumen |
| **Arango-Londoño y otros (2026)**, *Climate* 14(9), 188 | Valle del Cauca | 62 estaciones del IDEAM, CHIRPS **diario**, 2012–2020 | **No lo dice** | Modelo aditivo mixto funcional (FGAMM): curva **acumulada** anual de la estación ~ curva acumulada del satélite, con efectos de estación y año, y el SOI | Completo (versión preprint) |
| **Ocampo-Marulanda y otros (2021)**, *Data in Brief* 39 | Área metropolitana de Cali | 12 estaciones de la CVC y el IDEAM, **diarias**, 1969–2019 | **Hasta 6 % de faltantes** en el periodo; rellenan los índices de extremos con redes neuronales (NLPCA) | Índices de extremos (ETCCDI) | Completo |
| **CVC (2024)**, informe técnico del Grupo de Recursos Hídricos | Cuenca del río Cauca y Pacífico vallecaucano | De 486 series de la CVC quedan 136 estaciones **mensuales**, 1979–2023 | ≥20 años continuos con datos hasta 2008 o después; **máximo 40 % de faltantes** en el periodo común; atípicos contrastados con vecinas y con El Niño/La Niña; **imputación con redes neuronales** (NLPCA) | Tendencias (Mann-Kendall) | Métodos completos |
| **Cardona-Guerrero y otros (2014)**, *TecnoLógicas* 17(32) | Cuencas de los ríos Dagua y Cali | 18 estaciones (IDEAM, CVC, Cenicafé), **diarias**, de 21 a 56 años hasta 2010 | Inspección visual, control de calidad con RClimDex y **prueba de homogeneidad**: se descartan 2 de 20 estaciones por saltos bruscos; piden "un porcentaje bajo de datos faltantes", sin cifra | Tendencias de índices de extremos | Completo |
| **Puertas, Carvajal y Quintero (2011)**, *Dyna* 78(169) | Cuenca alta y media del Cauca | 108 estaciones **mensuales**, 1975–2006 | Según el buscador: **menos de 10 % de faltantes**, rellenados con el porcentaje mensual medio del total anual. **Falta verificarlo en el texto** | Tendencias (Mann-Kendall) | Solo el resumen |
| **Blanco y otros (2023)**, *Data in Brief* 50 | Departamento del Cauca | 12 estaciones del IDEAM, **diarias**, 2015–2021, y CHIRPS, ERA5-Land y MSWX | Se quitan errores de digitación y se verifica la continuidad; **no rellenan faltantes a propósito**, "para no introducir errores" en las métricas de validación | Reducción de escala con kriging; validación contra las estaciones | Completo |
| **Del Castillo-Gómez y otros (2023)**, *RBRH* 28 | Nariño | 46 estaciones del IDEAM, **mensuales**, 1983–2019 | Faltantes por estación de 0 % a 13 %; no fijan un umbral | Comparan cuatro métodos de relleno; gana la red neuronal (NLPCA) | Completo |
| **López-Bermeo y otros (2024)**, *J. South Am. Earth Sci.* 139 | Cuenca alta del Cauca | 30 pluviómetros, 1981–2018 | No se encontró | Corrección de sesgo de CHIRPS por cuantiles; empeora sobre 2 000 m | Solo el resumen |
| **Ocampo-Marulanda y otros (2022)**, *Ain Shams Eng. J.* 13 | Suroccidente de Colombia | 46 series **mensuales**, CHIRPS píxel contra punto | No se encontró | Validación de CHIRPS con componentes principales combinados; CHIRPS sobreestima en la costa Pacífica y subestima en las zonas andinas más lluviosas | Solo el resumen |

---

## 2. Los dos trabajos funcionales del Valle

### 2.1 Ospina-Galindez, Giraldo y Andrade-Bejarano (2019): el trabajo del profesor

Es el mismo problema que el nuestro: el Valle del Cauca, estaciones de la CVC, CHIRP pentadal y regresión funcional concurrente.

- **Datos:**
  - 82 estaciones entre 700 y 2 000 m, de 1993 a 2011 según la tesis.
  - Los faltantes se **imputan** con kriging ordinario año por año.
  - Excluye las estaciones con más del 20 % de faltantes: de 210 estaciones sacó 118 (tesis, cap. 5). El criterio es sobre todo el periodo, no mes a mes.
  - Rellena en diario y después suma a totales de 5 días (péntadas).
  - El texto da dos límites inferiores de altitud (700 y 900 m). Las estaciones van de 721 a 1 915 m.
- **Curvas:**
  - lluvia **acumulada** por pentadas, con B-splines monótonas;
  - K = 10 en estaciones y K = 20 en satélite.
- **Modelo:**
  - concurrente con la altitud;
  - los errores se correlacionan en el espacio (traza-variograma, Giraldo 2009);
  - se ajusta por mínimos cuadrados generalizados.
- **Hallazgos que podemos contrastar con los nuestros:**
  - el satélite no capta la orografía;
  - la altitud no aporta en la regresión;
  - el satélite aporta más en los años Niña;
  - los residuos tienen tendencia al final del año **por la acumulación**.

**Diferencia con lo nuestro:**

| | Ospina | Nosotros |
|---|---|---|
| Faltantes | Imputa con kriging ordinario; semivariograma de la media anual, por año | Imputa con kriging ordinario los huecos que acepta R5; semivariograma de la lluvia diaria de entrenamiento y vecinas escaladas por la razón de medias mensuales (D-61) |
| Curvas | Lluvia acumulada | mm/día, sin acumular |
| Partición | Dejando una estación fuera | Por años: entrenamiento, validación y prueba |
| Periodo | 1993–2011 | 1981–2003 |
| Estaciones | 82, de 700 a 2 000 m | 163 con ≥10 años válidos, de toda la altitud (100 con la regla del 80 % de la OMM, como sensibilidad) |

### 2.2 Arango-Londoño y otros (2026): FGAMM con CHIRPS diario

- **Datos:**
  - 62 estaciones del IDEAM, 2012–2020;
  - el satélite se extrae en el píxel de cada estación.
  - **El artículo no dice cómo se trataron los faltantes ni qué estaciones se descartaron.**
- **Curvas:**
  - lluvia **acumulada** del 1 de enero al 31 de diciembre;
  - se suaviza con P-splines (`pspline`).
  - Los autores dicen que acumulan porque la lluvia diaria es "demasiado variable" para modelarla.
- **Modelo:**
  - función sobre función, estación(t) = β₀(t) + ∫β₁(t,s)·satélite(s)ds + efectos de estación y de año;
  - se ajusta con `refund`.
  - La altitud, la latitud, la longitud y el SOI **no son significativos**. Según los autores, el efecto del año absorbe el ENSO.
- **Validación:**
  - 70 % de las estaciones para entrenar y 30 % para validar, repetido 50 veces;
  - error cuadrático medio de 0,68 mm/día.
  - La versión publicada advierte que la ventaja no se reproduce en otra base de datos.

**Lectura crítica:**

1. **El error de la curva acumulada no es el error diario.** Una curva acumulada es suave por construcción, así que su error en "mm/día" no se compara con el de un modelo de lluvia diaria. Es lo que advierte el profesor: la acumulación "arrastra" todo, y su propia tesis encontró residuos con tendencia al final del año.
2. **Sin un criterio de faltantes publicado**, no se puede reproducir qué entró al modelo. Eso es justo lo que el profesor nos pide justificar.
3. **Por qué su SOI no es significativo:** tienen un efecto por año, y cada año ya captura su fase ENSO. Nosotros no tenemos efecto de año y metemos el ONI directo, que es la alternativa que ellos mismos proponen como trabajo futuro.

**Probablemente estos dos son los "dos papers" que el profesor mostró para que no se siguieran:** los dos modelan la lluvia acumulada. Es una inferencia nuestra; hay que confirmarla con él.

---

## 3. Qué aprendemos para nuestras decisiones

### 3.1 Criterios de inclusión: lo que se usa en la región

| Fuente | Escala | Umbral de faltantes | ¿Rellena? |
|---|---|---|---|
| Ocampo-Marulanda y otros (2021), Cali | Diaria | ≤6 % en el periodo | Sí: índices, con NLPCA |
| Puertas y otros (2011), alto Cauca | Mensual | <10 % (por verificar) | Sí: proporción mensual |
| CVC (2024) | Mensual | ≤40 % en el periodo común, ≥20 años | Sí: NLPCA |
| Blanco y otros (2023), Cauca | Diaria | Continuidad, sin cifra | **No, a propósito** |
| Tesis de Ospina (2016), cap. 5 | Diaria para depurar y rellenar; pentadal para modelar | >20 % → fuera (quedaron 82 de 210 estaciones); además, solo de 700–900 a 2 000 m | Sí: kriging diario, y después suma a péntadas |
| IDEAM (2021), manual M-GDI-M-M005 | Diaria → mensual | Total mensual solo con ≥20 días con dato | No |
| **Nosotros (OMM 2017)** | **Diaria** | **Por mes: ≤10 faltantes y <5 seguidos; por estación: ≥10 años válidos (sensibilidad: ≥19 de 23, OMM)** | **Sí, solo esos huecos: kriging diario con las vecinas (D-61)** |

**Lectura:**

- **El único criterio mes a mes que encontramos es el del IDEAM** (≥20 días con dato para el total mensual), y es prácticamente el nuestro: ≤10 faltantes en un mes de 30 días. La OMM añade el límite de menos de 5 faltantes seguidos. Los estudios de la región usan un porcentaje sobre todo el periodo (6 %, 10 %, 20 %, 40 %), y eso permite que todos los faltantes se junten en un mismo mes.
- **Nuestra regla es más estricta y más explícita** porque viene de la norma de la OMM.
- En la práctica, nuestros datos quedan como los del estudio de Cali (≤6 %): en los años válidos falta el 0,05 % de los días.
- **Rellenar es lo común en la región**, sobre todo a escala **mensual**. Nosotros rellenamos solo los pocos días que acepta R5 (719 de 1 111 790 en el conjunto principal), con las vecinas y no con el satélite, porque para lluvia diaria los días de antes y después casi no informan (OMM 2017, 4.6).
- Blanco y otros (2023) no rellenan cuando las estaciones sirven para **validar el satélite**, que es nuestro caso, para no contaminar las métricas. Por eso no rellenamos con el satélite, y cada día rellenado queda marcado para poder calcular las métricas solo con días observados.

### 3.2 La regla de la vecina para los ceros (R4) tiene antecedente local

El informe de la CVC (2024) contrasta los atípicos "con registros de estaciones vecinas, considerando eventos históricos relevantes como El Niño y La Niña". Es la misma lógica de nuestra R4: la vecina decide, y el inventario por escenarios cruza la temporada con la fase ENSO. Cardona-Guerrero y otros (2014) agregan una prueba de homogeneidad, que nosotros no hacemos (ver la sección 4).

**¿Hay que tratar distinto los años Niño?** El Niño reduce la lluvia y el número de días con lluvia en el valle del río Cauca, sobre todo en las temporadas secas (Peña, Cortés y Montealegre, 2001; Cortés, 2010). La hoja metodológica del IDEAM (2025) pide tener en cuenta la fase ENSO al buscar atípicos. Lo probamos con los datos (prueba exploratoria, fuera del script):

| q que usa R4 = P(llueve en la estación \| llueve en la vecina) | Rachas falsas (de 909 con vecina) |
|---|---|
| Todos los días de entrenamiento (R4 anterior) | 30 |
| Solo años de la misma fase ENSO | 29 |
| Solo los mismos meses (mes central ± 1) | 18 |
| Misma fase y mismos meses | 17, pero 54 rachas quedan con menos de 30 días para estimar q |

- **La fase casi no cambia nada.** R4 no compara con lo que suele llover; compara con lo que llovió en la vecina esos mismos días. En una sequía Niño real la vecina también está seca: en las rachas de años Niño la vecina llovió una mediana de 2 días; en las 17 rachas falsas de años Niño, 11 días.
- **La temporada sí importa.** En los meses secos la lluvia es más local, así que el que llueva en la vecina predice peor que llueva en la estación (q baja de 0,62 a 0,58 en años Niño y de 0,57 a 0,52 en años neutros). Con la q de todo el año, R4 exagera lo raro de la racha y marca de más.
- **Aplicado (D-61):** q se estima con los mismos meses de la racha, no con todo el año ni por fase. Recoge lo que dice Cenicaña (el efecto del ENSO pasa por las temporadas secas) sin partir los datos en grupos muy pequeños.
- **Y con varias vecinas (D-61):** la racha pasa a faltante solo si la contradicen al menos 2 de sus 3 vecinas más cercanas, el principio de la regla de periodos secos de GSDR-QC (Lewis et al. 2021). Con eso R4 marca 10 rachas en vez de 30. Ninguna de las fuentes de control de calidad revisadas (GSDR-QC, GHCN-Daily, Terán-Chaves et al. 2021 con datos del IDEAM, INQC) usa la fase ENSO: todas protegen las sequías reales con las vecinas o con la propia climatología de la estación. Detalle en `docs/FUNDAMENTOS_REGLAS.md`, sección 3.

**El IDEAM y los acumulados.** El manual de control de calidad del IDEAM (2021, §7.2.1) dice: "Solo se aceptan hasta cuatro (4) días acumulados". Ordena también comparar con estaciones vecinas los valores de 135 mm o más, que es la capacidad del colector. Una racha de 30 ceros que termina en el total del mes no cumple esa regla; es lo que buscan nuestras R3 y R4.

### 3.3 Lluvia acumulada

- Los dos trabajos funcionales del Valle acumulan. El profesor nos pidió no hacerlo.
- Nuestra forma (mm/día, sin restricciones de forma) es la que él pidió.
- Además evita lo que encontró su propia tesis: residuos con tendencia al final del año.

### 3.4 Resultados que coinciden con los nuestros

| Hallazgo nuestro | Lo mismo en |
|---|---|
| El satélite sobreestima en la zona andina, +28,8 % anual en prueba | Arango-Londoño (2026): CHIRPS sobreestima en Cali, Dagua y El Dovio |
| La altitud aporta poco una vez está el satélite | Tesis de Ospina; Arango-Londoño (2026) |
| El error del satélite depende del ENSO | Tesis de Ospina: el satélite aporta más en años Niña |
| En el Pacífico el satélite no coincide con las pocas estaciones que quedan | Ocampo-Marulanda (2022): CHIRPS sobreestima en la costa Pacífica |

---

### 3.5 Escala: diaria o pentadal

- **La depuración siempre se hace en diario.** Así lo hacen el IDEAM (manual de 2021), la OMM y la tesis del profesor, que rellena en diario y suma a péntadas después. Las reglas R2 a R4 solo se pueden aplicar día a día.
- **CHIRP es pentadal de origen.** El CHIRPS diario se obtiene repartiendo la péntada entre los días con campos de un modelo (Funk y otros, 2015). En Colombia, CHIRPS reproduce bien la media y la estacionalidad en diario, mensual y anual, pero subestima la varianza diaria (Urrea y otros, 2016, 904 sitios; por verificar).
- **En nuestros datos** (prueba en `scratchpad`):
  - con K = 8, las curvas armadas desde datos diarios, pentadales o semanales difieren en 0,16 mm/día, frente a una dispersión de 2,85;
  - el 71 % de la varianza diaria está dentro de la péntada, donde el satélite no ve nada;
  - la correlación mediana estación–satélite es 0,22 en diario, 0,43 en péntadas y 0,67 por mes;
  - las 52 semanas del año cruzan el borde de una péntada.

## 4. Lo que estos trabajos hacen y nosotros no

1. **Prueba de homogeneidad** de las series (Cardona-Guerrero y otros, 2014; CVC, 2024, que cita a WMO-TD 1186). Detecta saltos por cambio de instrumento o de sitio. No la hacemos. Se puede agregar como verificación a las 98 estaciones.
2. **Correlación espacial de los errores** (Ospina y otros, 2019). Ya tenemos el traza-variograma de los residuos, pero no reajustamos con mínimos cuadrados generalizados.
3. **Validación dejando estaciones fuera** (tesis de Ospina; Arango-Londoño y otros, 2026). Nosotros validamos por años. Las dos responden preguntas distintas:
   - por años: ¿sirve el modelo en años que no vio?
   - por estaciones: ¿sirve en sitios que no vio?

---

## 5. Pendientes antes de citar en el informe

- **Urrea y otros (2016)** y **Peña y otros (2001):** confirmar autores y datos antes de citarlos.
- **Ospina y otros (2019):** conseguir el texto del artículo para confirmar el periodo y si las curvas son acumuladas. El acceso por la web está bloqueado (Taylor & Francis); se puede pedir al profesor.
- **Puertas y otros (2011):** verificar en el PDF el umbral del 10 % y el método de relleno.
- **López-Bermeo y otros (2024) y Ocampo-Marulanda y otros (2022):** solo se leyó el resumen; ScienceDirect bloquea la lectura automática.

---

## Referencias

**Normas y documentos técnicos de Colombia**

- IDEAM (2021). *Manual para el control de calidad de los datos meteorológicos en la preverificación, captura y verificación* (M-GDI-M-M005, versión 1, 25/03/2021). Subdirección de Meteorología. https://ideam.gov.co/sites/default/files/mapa-de-procesos/m-gdi-m-m005_manual_control_de_calidad_verificacion_v1.pdf
- IDEAM (2025). *Hoja metodológica del indicador Precipitación total* (GCI-OE-F002, versión 03, 10/07/2025). https://bart.ideam.gov.co/indiecosistemas/ind/clima/hm/HM_precipitacion.pdf
- Cortés Betancourt, E. (2010). Incidencia de "El Niño" y "La Niña" en el clima del valle del río Cauca. *Carta Trimestral Cenicaña*, 32(3–4), 6–7. https://www.cenicana.org/pdf_privado/carta_trimestral/ct2010/ct3y4_10/ct3y4_10_p6-7.pdf
- Peña, A. J., Cortés, E. y Montealegre, F. (2001). Incidencia de los fenómenos "El Niño" y "La Niña" sobre las condiciones climáticas en el valle del río Cauca, partes I y II. *Meteorología Colombiana*, 3, 103–108 y 119–128. (Citado por Cortés, 2010; no leído).
- Urrea, V., Ochoa, A. y Mesa, O. (2016). Validación de la base de datos de precipitación CHIRPS para Colombia a escala diaria, mensual y anual en el período 1981–2014. (Autores, año y evento por verificar; solo se leyó el resumen). https://www.researchgate.net/publication/310844678
- Funk, C. y otros (2015). The climate hazards infrared precipitation with stations: a new environmental record for monitoring extremes. *Scientific Data*, 2, 150066. https://doi.org/10.1038/sdata.2015.66
- Portillo y Medina (2012), precipitación diaria del Valle con un modelo lineal mixto generalizado, 1971–2009. Lo cita la tesis de Ospina (2016); no lo encontramos en la web.

**Estudios de la región**


- Arango-Londoño, D., Ortega-Lenis, D., Mazo-Lopera, M. A., Aparicio, J. S., Soto, D. y Moraga, P. (2026). Satellite-Based Daily Precipitation Bias Correction in a Tropical Mountainous Region Using Functional Generalized Additive Mixed Models: A Case Study in Valle del Cauca, Colombia. *Climate*, 14(9), 188. https://doi.org/10.3390/cli14090188 (preprint: https://www.preprints.org/manuscript/202605.1424)
- Blanco, K., Villamizar, S. R., Avila-Diaz, A., Marceló-Díaz, C., Santamaría, E. y Lesmes, M. C. (2023). Daily dataset of precipitation and temperature in the Department of Cauca, Colombia. *Data in Brief*, 50, 109542. https://doi.org/10.1016/j.dib.2023.109542
- Cardona-Guerrero, F., Ávila-Díaz, Á. J., Carvajal-Escobar, Y. y Jiménez-Escobar, H. (2014). Tendencias en las series de precipitación en dos cuencas torrenciales andinas del Valle del Cauca (Colombia). *TecnoLógicas*, 17(32), 85–95. http://www.scielo.org.co/scielo.php?script=sci_arttext&pid=S0123-77992014000100009
- CVC, Grupo de Recursos Hídricos (2024). *Análisis espaciotemporal de las tendencias de precipitación y temperatura en la cuenca del río Cauca y el Pacífico vallecaucano*. Mera Franco, J. D. y Soto Chávez, C. Y. https://ecopedia.cvc.gov.co/sites/default/files/archivosAdjuntos/Informe%20Tendencias%20Precipitacion%20Temperatura.pdf
- Del Castillo-Gómez, J. S., Canchala, T., Torres-López, W. A., Carvajal-Escobar, Y. y Ocampo-Marulanda, C. (2023). Estimation of monthly rainfall missing data in Southwestern Colombia: comparing different methods. *RBRH*, 28. https://doi.org/10.1590/2318-0331.282320230008
- López-Bermeo, C., Montoya, R. D., Caro-Lopera, F. J. y Dijón-Caro, L. F. (2024). Bias-corrected high-resolution precipitation datasets assessment over a tropical mountainous region in Colombia: a case of study in Upper Cauca River Basin. *Journal of South American Earth Sciences*, 139, 104870.
- Ocampo-Marulanda, C., Cerón, W. L., Avila-Diaz, A., Canchala, T., Alfonso-Morales, W., Kayano, M. T. y Torres, R. R. (2021). Missing data estimation in extreme rainfall indices for the Metropolitan area of Cali – Colombia: An approach based on artificial neural networks. *Data in Brief*, 39, 107592. https://doi.org/10.1016/j.dib.2021.107592
- Ocampo-Marulanda, C., Fernández-Álvarez, C., Cerón, W. L., Bedoya-Soto, J. M., Lomba-Romero, A., Cerda-Uribe, C. E. y Avila-Díaz, A. (2022). A spatiotemporal assessment of the high-resolution CHIRPS rainfall dataset in southwestern Colombia using combined principal component analysis. *Ain Shams Engineering Journal*, 13(5), 101739. https://doi.org/10.1016/j.asej.2022.101739
- Ospina-Galindez, J. A., Giraldo, R. y Andrade-Bejarano, M. (2019). Functional regression concurrent model with spatially correlated errors: application to rainfall ground validation. *Journal of Applied Statistics*, 46(8). https://doi.org/10.1080/02664763.2018.1544231
- Ospina Galíndez, J. A. (2016). *Modelo de regresión funcional de la precipitación pentadal, integrando datos de estaciones climáticas y de satélite*. Tesis de Maestría en Estadística, Universidad del Valle.
- Puertas Orozco, O. L., Carvajal Escobar, Y. y Quintero Angel, M. (2011). Estudio de tendencias de la precipitación mensual en la cuenca alta-media del río Cauca, Colombia. *Dyna*, 78(169), 112–120. https://revistas.unal.edu.co/index.php/dyna/article/view/21330

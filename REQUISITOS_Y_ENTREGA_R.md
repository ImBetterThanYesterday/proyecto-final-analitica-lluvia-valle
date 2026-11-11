# Requisitos textuales y entrega R — 9 de noviembre de 2026

> **Versión vigente, revisada nuevamente contra originales.** La auditoría
> independiente de esta revisión, sus diferencias y decisiones están en la
> [sección 8](#8-revisión-contra-originales-y-decisiones-corregidas).
> No interpretar los resultados de la primera intervención como prueba de
> equivalencia literal con todos los scripts docentes. La antigua guía Word es histórica.
> **Word vigente:** `INFORME_RECONSTRUCCION_EDA_ACTUAL.docx` (43 figuras, resultados R + diagramas).
> **Recorrido vigente:** BLOQUES 0–13 al inicio del único R, véase §10.
> **Actualización por fe de erratas:** norma vigente `Proyecto_Final_Analitica_FeErratas.docx`.
> Periodo **1981–2003**. Comparación íntegra: [COMPARACION_WORD_FE_ERRATAS.md](COMPARACION_WORD_FE_ERRATAS.md).
> Las cifras de §§3–6 son actuales; §8.4 conserva evidencia histórica explícita.

## 1. Fuente normativa: citas exactas antes de decidir el alcance

Se leyeron completos `Proyecto_Final_Analitica_FeErratas.docx` y `Proyecto_Final_Analitica.docx` extrayendo `word/document.xml`,
concatenando todos los runs de cada párrafo y recorriendo también las tablas.
Se citan **secciones y celdas de rúbrica**, no páginas de Word: la paginación
depende del renderizador y no se inventó un número de página.

| Sección del Word | Cita literal relevante | Consecuencia aplicada |
|---|---|---|
| Objetivos específicos, 1 | «Construir datos funcionales paso a paso (sin usar librerias).» | Φ, MCO y selección se programan en R estándar sin `fda`, `splines`, `smooth.basis` ni helpers de paquetes. |
| Objetivos específicos, 2 | «Describir el comportamiento de los datos funcionales (pueden usar librerias)» | La prohibición no se extiende al EDA. |
| Reglas generales | «Solo se pueden usar los scripts, procedimientos y metodologías vistos en clase. Las librerías adicionales se permiten en los apartados mencionados.» | No incorporar detectores de temporadas ni bandas ad hoc ajenos al material. |
| Rúbrica, criterio 1 | «Matriz de la base Φ y coeficientes por mínimos cuadrados programados paso a paso, sin librerías.» | Recursión propia y sistema normal resuelto explícitamente. |
| Rúbrica, criterio 1 | «Selección de K (o de λ) con ECM, VC y VCG: gráfica y valor elegido para cada fuente.» | Búsqueda completa repetida en R para ambas fuentes; no K fijado a mano. |
| Rúbrica, criterio 1 | «Verificación numérica de la matriz de suavizado: S idempotente y tr(S) = K, o ν = tr(S_λ) con penalización.» | Se eligió MCO; se verifica simetría, idempotencia y traza para todas las máscaras incompletas y la completa. |
| Rúbrica, criterio 1 | «Covarianza residual y bandas de confianza punto a punto, interpretadas como punto a punto.» | Segundo momento residual de clase y banda sigma²I de cada fuente, ejemplo identificado, sigma² de todas las curvas. |
| Rúbrica, criterio 1 | «Curvas de ambas fuentes (diaria y péntada) construidas sobre el mismo dominio anual.» | Dominio 1–365, observación diaria/pentadal diferenciadas. |
| Rúbrica, criterio 1, tope | «Construir las curvas con librerías: máximo 2.» | No se usa una librería de construcción ni siquiera para verificar dentro de esa fase. |
| Rúbrica, criterio 2 | «Media funcional calculada paso a paso y verificada contra la función de librería.» | Verificación real contra `fda::mean.fd` en ambas fuentes. |
| Rúbrica, criterio 2 | «Varianza y desviación estándar funcional.» / «Covarianza y correlación funcional de cada fuente.» | Tablas y figuras independientes por fuente. |
| Rúbrica, criterio 2 | «Derivadas e interpretación del inicio y el fin de las temporadas de lluvia.» | Derivada calculada y discutida; fechas meteorológicas no validadas, requisito parcialmente pendiente. |
| Rúbrica, criterio 2 | «Profundidad (Fraiman–Muniz y MBD) y boxplot funcional, con las curvas atípicas identificadas por estación y año.» | Ambas profundidades sobre todas las curvas del universo seleccionado, tablas de atípicas con claves. |
| Rúbrica, criterio 2 | «Interpretación climatológica de cada resultado (régimen bimodal, orografía, ENSO).» | Curvas de contexto por grupos descriptivos y ONI, límites de interpretación indicados abajo. |
| Rúbrica, criterio 2, tope | «Analizar solo un subconjunto de años o estaciones sin justificar que es representativo del periodo 1981–2003: máximo 2.» | Se incluyen todas las curvas preparadas de la ventana; todos los años presentes. Panel de estaciones incompleto y representatividad aún no demostrada. |
| Rúbrica, criterio 4.4 | «Partición 1981–1995 / 1996–2000 / 2001–2003 y estaciones reservadas por el docente, sin fuga de información.» | Actualiza la fase futura; el EDA global no constituye selección de modelo libre de fuga. |
| Reglas generales / criterio 6 | «Desarrollar el proyecto en R y Python.» / «Un único archivo de R y un único archivo de Python, separados del informe.» | Se entrega un principal R autónomo para estas fases según el alcance solicitado; no se declara completado el requisito global bilingüe. |
| Reglas generales / criterio 5 | «Si se usa alguna, deben entregarse en un documento de Word aparte los prompts utilizados y el nombre de la herramienta» | Esta intervención usó OpenCode, modelo `openai/gpt-6.1-sol`. Se creó un Word de estudio con resultados; el anexo separado de IA con prompts reales sigue pendiente. No se inventó un historial completo. |

No se desarrollaron las fases 3/4 (comparación pareada, regresión o redes).
La autorización PyTorch/pykan del apartado de red no autoriza su uso para
construcción. La restricción de no librerías **sí aplica a construcción**;
el permiso del EDA no la revoca.

## 2. Material de clase contrastado

Las páginas son físicas, incluyendo diapositivas repetidas:

- `../analitica/ConstruccionDatosFuncionales.pdf`: pp. 3–5 (MCO, Φ, S),
  20 (ECM de ajuste insuficiente), 25 y 30 (VC y atajo), 31 (VCG),
  33 (N^-1 Σ e_i e_i', no se confunde con covarianza centrada muestral),
  36 y 41 (varianza del ajuste y bandas puntuales), 46 (alternativa penalizada).
- `AnalisisExploratorioDatosFuncionales.pdf`: pp. 13 y 20 (coeficientes,
  media y verificación), 24 (varianza, DE directa y CV), 30–31 (covarianza y
  correlación), 43 y 48 (derivadas e interpretación), 53 (FM/MBD),
  57–60 (mediana, recorte, envolvente), 62 (atípicas), 17 y 37 (pares necesarios).
- `../analitica/Script_ConstrucciÃ³n_FDA.R`: 73–105, 215–305, 350–387.
  El script docente usa paquetes para bases y contraste; el Word del proyecto
  impone la programación manual de la construcción para esta entrega.
- `Script_AnalisisExploratorio.R`: 155–196, 253–289, 351–364, 445–518.
- `Script_t-Funcional.R`: 74–95 y 129–143, ajuste con observados, allí
  penalizado. Aquí se conserva MCO con K reducido y prueba de rango.

No se usaron métodos de inferencia ni regresión de los otros PDF/scripts:
pertenecen a fases fuera del alcance solicitado.

## 3. Decisiones implementadas y entradas heredadas

1. **Un solo script R secuencial** con funciones reutilizables. 03/04 son accesos
   compatibles. El EDA lee `resultados_R/construccion.rds`, nunca coeficientes
   históricos Python ni K constantes. `source()` ejecuta todo el análisis.
2. **Matrices observadas existentes como punto de entrada**, no se vuelve a
   limpiar desde archivos crudos ausentes. Se mantienen IDs y la máscara de
   observados en el RDS. El ajuste no modifica NA ni el dato original.
3. **B-splines cúbicos manuales**, nodos equiespaciados y dominio 1–365;
   selección descriptiva K por VCG entre 4:40 y 45:100 cada 5 en estaciones,
   4:40 en satélite. Estaciones: selección solo con 3.761 curvas completas,
   ajuste posterior de las 665 incompletas. Satélite: 33.189 completas (1981–2003).
   Se filtra por esquema explícito `STATION_ID,YEAR` / `anio,fila,col` antes de
   criterios y ajuste. Año inválido, claves duplicadas o ventana vacía se rechazan.
4. **ECM/VC/VCG exactos de clase**, calculados agregando cuadrados con Y'Y
   para evitar repetir productos de grandes matrices. Se contrasta contra
   residuos explícitos y VC dejando uno fuera en pruebas independientes.
5. **Bandas sigma²I** y segundo momento residual, por fuente. Se elimina del
   flujo vigente la banda estacional heredada sin receta de clase. Se exporta
   sigma² por curva, Φ observada, coeficientes y un ejemplo de banda por fuente.
   Las bandas no son simultáneas, ni intervalos de predicción, ni media ±2 DE.
6. **Máscara de tierra heredada**, transportada por claves del CSV histórico,
   sin usar sus coeficientes: 31.119 curvas en ventana (43.296 históricas). Su origen es la máscara de CHIRPS
   del notebook 02, no demuestra identidad del satélite entregado con CHIRPS.
7. **EDA completo del universo disponible**, con media, DE, CV, superficies,
   derivadas, FM y MBD inclusiva, medianas FM/MBD, recorte FM 25 %, envolvente
   50 % y factor 1,5. Tras esta revisión, `FM` usa la escala de `depth.FM()`
   del script docente; `FM_pdf` conserva la escala de la fórmula del PDF.
   Se promedian tiempos uniformes, como el código de `depth.FM`, no se reporta
   una integral sin normalizar. El recorte incluye `FM >= quantile(FM,.25)`;
   mediana y región central resuelven empates por orden original. Los límites
   de detección y los bigotes (envolvente no atípica) se exportan por separado.
8. **Retiro de fechas de temporadas heurísticas**: se conserva la derivada de
   clase y la interpretación de su signo. Máximos de pendiente no se llaman
   inicio meteorológico. No se incorporó otro detector sin respaldo docente.
9. **Grupos de contexto honestamente etiquetados**: `prefijo_26`/otro y
   corte de longitud −76,8 heredado, no cuencas confirmadas. Fases ENSO por
   promedio anual de ONI de `eda_fase0/oni_anual_fase.csv`, procedente del
   notebook 02 (NOAA CPC; copia local
   `datos_web/ONI_NOAA_CPC_1981_2013.csv` y `oni_noaa_cpc_completo.txt`).
   **Corrección tras leer el notebook 02:** su celda Markdown 0 sí registra
   `https://cpc.ncep.noaa.gov/data/indices/oni.ascii.txt` y noviembre de 2026.
   La celda de código lee una copia previa; falta evidencia de fecha exacta
   y versión de la descarga, no la URL documental.
   No equivalen automáticamente a episodios oficiales ni a pareamiento.
10. **fda instalado solo por necesidad de la rúbrica E01**, con sus 14
    dependencias faltantes, en `.Rlib/`; no se instala al ejecutar análisis.
    No se utiliza en la fase de construcción. No se instaló Python, scikit-fda
     ni herramientas nuevas de extracción de documentos.
11. **Ventana normativa única:** `ANIOS <- 1981:2003`. RDS con periodo declarado;
    EDA rechaza construcciones antiguas, valida años y exige ONI sin duplicados
    con los 23 años. Grupos/ONI parten de las curvas filtradas. Se exportan
    cobertura anual, combinaciones sitio–año ausentes y conteos del filtro.
12. **Presentación vigente (ampliada en §10):** preparación de funciones primero,
    `# CONSTRUCCIÓN FUNCIONAL ----` y `# EDA FUNCIONAL ----` después, con código
    secuencial directo. Source ejecuta; Cmd/Ctrl+Enter permite avanzar desde la
    primera línea. PDF/PNG más dispositivo interactivo, resumen en consola y
    `ver_tabla()` con View opcional. El dibujo se repite, el cálculo no.

## 4. Verificación efectuada y resultados

Comandos realmente ejecutados desde la raíz del proyecto:

```bash
Rscript test_fda_R.R
Rscript test_periodo_R.R
Rscript test_usabilidad_R.R
Rscript proyecto_funcional.R
Rscript verificar_flujo_R.R
```

La prueba añadida **falló antes de corregir** la aceptación de NA en
`profundidades()`: `rechaza(... matrix(c(NA,1,2,3),2)) is not TRUE`.
Ahora rechaza entradas no finitas/menos de dos curvas. También se observó
un fallo de lectura por doble cierre de conexión; se corrigió abriendo
explícitamente la conexión antes de `read.csv`. El verificador de grupos
falló por ausencia de la exportación antes de implementar ese componente.

Pruebas matemáticas que pasaron:

- MBD inclusiva contra definición exhaustiva por pares, empates y una malla
  de un solo punto; curvas idénticas dan MBD=1.
- Partición de unidad/no negatividad de la base; reproducción de un cúbico
  conocido; derivada analítica, incluidos extremos 1 y 365.
- MCO con huecos y rechazo de diseño no estimable, H y banda sigma²I.
- ECM/VC/VCG agregados contra residuos explícitos; atajo VC contra QR
  dejando un punto fuera cada vez.
- **12 controles integrales OK**, incluidos todos los coeficientes contra
  QR independiente en 4.426 + 33.189 curvas; IDs y máscaras observadas;
  momentos y señalamiento de atípicas exportados; conteos de grupos;
  medias verificadas con librería, periodo/cobertura y suma de curvas por ONI.
  Las cuatro entradas de datos y los tres Word comparados por MD5 permanecieron
  intactos (`integridad_entradas.csv`). Máscara filtrada contrastada por nombres,
  tipos, valores y orden de claves; se excluye solo el atributo row.names, pues
  numera filas distintas en dos CSV diferentes.
- Prueba temporal añadida **antes del código**: falló por función ausente; tras
  implementar pasan filtro de ambas fuentes, rechazo de esquemas/años inválidos,
  duplicados/ventana vacía, huecos interiores pese a min/max correctos, panel
  sitio–año, NA reales, ONI incompleto/duplicado y rechazo de RDS antiguo.
- Usabilidad: importación aislada del arnés, PDF y dispositivo visible PNG,
  construcción/EDA sintéticos reales, limpieza del dispositivo por error, tabla
  y ejecución secuencial/Source normal. No se afirma inspección visual de una GUI RStudio real.

| Resultado nuevo | Estaciones | Satélite en tierra |
|---|---:|---:|
| Curvas EDA | 4.426 | 31.119 |
| K por VC y VCG | 6 | 11 |
| Máximo error coeficientes contra QR (satélite completo) | 4,83e-13 | 8,53e-14 |
| Máximo error de los controles H | 1,07e-14 | 3,55e-15 |
| Media contra `fda::mean.fd` | 1,78e-14 | 9,41e-14 |
| Atípicas MBD globales | 170 | 0 |
| Puntos evaluados negativos | 24.361 | 408 |
| Mínimo ajustado, escala mm/día heredada | −11,6244 | −1,195546 |
| Correlación abril (105)–octubre (290) | 0,808861 | 0,824418 |

Las salidas nuevas son `resultados_R/*.csv`, ocho PDF de figuras (dos de
selección, dos de reconstrucción/bandas, dos EDA y dos de contexto), RDS y
`sessionInfo.txt`. El JSON `verificacion_fda_clase.json` conserva su fecha y
alcance histórico; **no se actualizó para hacer pasar evidencia Python por R**.

## 5. Lectura contextual de las salidas reales

- **Media y derivada**: en estaciones la derivada cambia de positiva a
  negativa entre días 114–115 y 316–317: máximos de la tasa anual alrededor
  de abril y noviembre. El mínimo intermedio está entre 220–221. En el
   satélite, máximos dominantes próximos a días 134–135 y 307–308, con
  oscilaciones adicionales pequeñas a principios/finales del año. Esto
  es compatible con dos épocas lluviosas, no prueba fechas exactas de inicio/fin.
- **Varianza y DE**: el máximo de DE terrestre está en el borde 365,
   mientras el satelital está en día 310. No se oculta el efecto de borde ni
  se atribuye automáticamente toda dispersión a clima: aquí se mezclan
  sitios, años y ajustes de curvas incompletas.
- **Covarianza/correlación**: la correlación positiva entre abril y octubre
  sugiere persistencia del orden relativo de curvas húmedas/seca entre
  épocas. No demuestra causalidad ni correlación cruzada entre fuentes.
- **Centralidad/atípicas**: las medianas y el recorte ofrecen resúmenes menos
  sensibles a las curvas más lluviosas. Las 170 atípicas son respecto de
  todas las estaciones seleccionadas; no se eliminó ninguna. Cero atípicas
  satelitales bajo esta regla no demuestra ausencia de problemas del producto.
- **Geografía**: media de tasa 11,279 en otros prefijos frente a 4,584 en
  prefijo 26; satélite 14,192 al oeste del corte −76,8 frente a 6,345 al
  este. Es compatible con un gradiente húmedo occidental, pero coordenadas
  y altitud de estaciones son necesarias para atribuirlo a orografía.
- **ONI**: las curvas por fase son descriptivas de la clasificación anual
  heredada. Diferencias entre las dos fuentes no se llaman sesgo ENSO:
  difieren sitios y ponderaciones, aunque ahora comparten periodo; no hay pares geográficos.
- **Residuos/bandas**: bajo errores iid la covarianza de los residuos es
   sigma²(I−H), no sigma²I. Por eso observar estructura residual no basta
   para refutar independencia: parte la induce el ajuste. Las bandas son
   condicionales a sigma²I y no corrigen sesgo. Las tasas negativas muestran
   el límite físico del MCO irrestricto; pasar controles no valida ese supuesto.

## 6. Límites reales y entrega pendiente

1. Crudos ausentes en esta raíz: no se reprodujo depuración desde cero.
   Se mantienen matrices preparadas y banderas históricas; no se han
   convertido umbrales exploratorios heredados en reglas autorizadas de clase.
2. Ventana vigente 1981–2003 cubierta en todos sus 23 años por ambas fuentes.
   No falta llegar a 2013: ese requisito fue corregido. Persisten 1.048
   combinaciones sin curva respecto de 238 estaciones × 23 años, y 2.871 NA
   diarios en 665 curvas incompletas. Se conservan, no se fabrican datos.
   El panel ausente es un diagnóstico sobre el catálogo observado, no prueba
   de que cada estación estuviera operativa toda la ventana ni sobre sitios
   completamente ausentes de los preparados. No acredita representatividad.
3. Unidad/producto y calendario exacto del satélite original siguen por
   confirmar. Los cálculos conservan la escala preparada sin mezclar CHIRPS
   externo ni rellenar años/meses con otro producto.
4. Inicio/fin meteorológicos y sensibilidad del ajuste siguen pendientes.
   Se reportan las curvas negativas, sin cambiar la solución MCO por recorte.
5. No hay comparación pareada ni altitud validada; no se entrenó red,
   no se implementó regresión adicional. La partición futura documentada se
   actualizó a 1981–1995 / 1996–2000 / 2001–2003; no se entrenaron modelos.
6. El PDF final, único Python global y anexo Word de uso de IA corresponden
   al cierre del proyecto completo. Esta entrega no afirma cumplirlos.
7. La guía Word y los árboles previos permanecen como material histórico;
   la versión vigente de instrucciones es README + este documento + aviso
   inicial de la guía Markdown. No se regeneraron Word/PDF documentales.

## 7. Continuidad para el coordinador

Ruta exclusiva trabajada:
`/Users/mariavalencia/Downloads/proyecto-final-analitica-lluvia-valle-main`.
No se tocó EduInteractiva ni se crearon commits. No se guardaron memorias
con el proyecto EduInteractiva: las herramientas de memoria disponibles
carecen de parámetro para seleccionar la ruta de otro proyecto. Este archivo
conserva requisitos, decisiones, correcciones y resultados para continuidad.

`skill_resolution=injected`: se aplicaron las reglas compactas recibidas.

## 8. Revisión contra originales y decisiones corregidas

### 8.1 Qué se leyó y qué significa «seguir la clase»

Se localizaron recursivamente los `.R` de esta carpeta y `../analitica/`;
se excluyeron únicamente los archivos internos de paquetes de `.Rlib/`.
Se leyeron íntegramente los nueve scripts docentes únicos de la tabla inferior,
además del principal, wrappers y verificadores del proyecto. Las cuatro copias
docentes de la raíz son idénticas byte a byte a las de `../analitica/` (SHA-256).
Los tres ZIP de `analitica` (`MapaCaliComuna`, `covariables`, `datos_casos`)
no contienen scripts `.R`. No se ejecutaron los originales que instalan paquetes,
borran el entorno o contienen rutas Windows.

| Script docente leído en `../analitica/` | Líneas | Finalidad efectiva | SHA-256 |
|---|---:|---|---|
| `Script_ConstrucciÃ³n_FDA.R` | 473 | MCO, H, selección y bandas; autoridad de construcción | `509e47bcb38a308d7f72fdd0b6b052e0ca698f6a8277460768e144ffed0770e5` |
| `Script_AnalisisExploratorio.R` | 530 | Momentos, derivadas, FM, recorte y boxplot | `c176836c021e24f1cc7b6eca1f6228fb5f67ba207a974a6ccdd98a2a8c993996` |
| `Script_t-Funcional.R` | 516 | Ajuste por días observados; inferencia no trasladada | `cd72dcc3fea9cbf17ce9a67db0cc95ea397151d147ded941b9142b619e334811` |
| `Script_anova_funcional.R` | 518 | Calidad de ICA, bordes y grupos; no convertir sus ceros/umbrales en reglas de lluvia | `315f1d539c9c591d938bebf5f127edaa6d0ccf670566668ad883412738c3b67c` |
| `Script_Regresion_Funcional.R` | 342 | Fourier, fRegress, fPCR, concurrente y linmod; fuera del alcance actual | `1f291aa75f4012673b618b00e1d7a5e891720f683e5d51ffe5bad116488a774a` |
| `Script_AFM.R` | 152 | Estandarización muestral, SVD, ponderación, RV; no se añadió AFM al FDA | `c55943f758a95f2d1dcbe9613e39c0687d82ccc50b0b85faa20524e0648bd79b` |
| `Script_AFM_Ejemplo2.R` | 268 | Estandarización poblacional, eigen y baricentro; no es intercambiable literalmente con el anterior | `fccfa52a24355eafed6f660b2d02295d61d026abb81df3d618f0b8b5c5884994` |
| `Practic_Geostatistics.R` | 118 | Meuse: exploración, logZn, variograma, kriging y CV | `b22d989a08e9c0e8f7e9f47f3a5dcf3b76a32c44f95dedf5ac2ae1dc9a8995aa` |
| `Practica_Lattice_Data_SEP291025-2.R` | 253 | NC: vecindades KNN/Rook/Queen, Moran/Geary y GLM; contexto espacial | `5fad7f7dd95081409572491080823c32d17d341e88b42fa576081e31113c8132` |

También se extrajeron por página física los PDF de construcción (49 páginas)
y EDA (66), y por párrafo concatenando runs el Word normativo y la guía Word.
Se leyeron los notebooks históricos **03 y 04 completos**, el **02 completo**
(calendario, máscara, ONI, catálogo y advertencias), el EDA original
`EDA_precipitacion_valle_cauca.ipynb` completo y se inspeccionaron las secciones
de preparación/calidad del `01_EDA_calidad_datos.ipynb`. Los otros cuatro
notebooks se localizaron y extrajeron, pero no se presentan como revisados
íntegramente. Extraer un archivo no es lo mismo que auditar todos sus métodos.

**Orden y notación.** La clase de construcción sigue Φ → coeficientes →
y2cMap → H → efecto de K → contraste con librería → VC/VCG → bandas →
penalización opcional (PDF p. 47). El principal agrupa funciones y ejecuta
selección de K → ajuste definitivo → H/residuos/banda → persistencia. No es
una copia del orden de exposición, pero los pasos matemáticos están mapeados.
En el proyecto `Y` es N×n y `C` es N×K; en clase se transponen: Y es n×N,
CoefEst es K×N y **C designa y2cMap**. Aquí `mapa` es el C de clase y
`Phi %*% mapa` es H/S. El EDA conserva media → dispersión → covarianza →
derivadas → profundidad; los cálculos pareados de clase se reservan hasta
tener pares. No se copiaron `rm(list=ls())`, instalaciones, `setwd()` ni GUI
Windows porque impedirían la ejecución reproducible en este equipo.

### 8.2 Discrepancias reales, correcciones y desviaciones explícitas

| Aspecto | Contraste directo | Decisión vigente |
|---|---|---|
| Librerías | Construcción docente carga fda/fda.usc y usa `create.bspline.basis`, `lsfit`, `smooth.basis` y `optim.basis`; el Word prohíbe construir con librerías | Φ/ajuste/criterios propios en R estándar. fda se consulta solo en EDA y pruebas de EDA. La restricción del proyecto prevalece. |
| Base | Construcción ilustra orden 2 y K=15; ejercicio línea 470 propone orden 4, y EDA 101–113 usa cúbicos para lluvia | Orden 4 respaldado. La recursión Cox–de Boor es implementación propia necesaria para hacer manual la base: no está desarrollada literalmente en los scripts/PDF revisados. No afirmar que el profesor entregó esa recursión; confirmar con él si exige un algoritmo manual concreto. |
| MCO/H | `solve(crossprod(Phi),t(Phi))`, SSE/(n−K), H simétrica/idempotente/traza K | Fórmulas coinciden; ajuste de huecos con filas observadas adapta t-funcional, cuyo ejemplo es penalizado. No se añade penalización. |
| Selección K | Clase usa 4:20 en SO2 y suma VCG por curva | Intervalos 4:40/45:100 y 4:40 se adaptan a n=365/72; ECM y VC promedios, VCG suma. Y'Y es reordenamiento algebraico, contrastado con residuos/LOO. La selección terrestre usa solo completas: posible sesgo de selección, no representatividad demostrada. |
| Residuos/bandas | PDF construcción pp. 33, 36, 41; script 350–387 | `crossprod(Res)/N_completas` es segundo momento residual, no covarianza centrada con N−1 ni Σ de errores sin más. EE usa sigma²I de cada curva y banda ±2EE; un ejemplo por fuente. No usar el segundo momento residual directamente como Σ en la banda: C(I−H)=0. |
| Media/varianza/covarianza | EDA 155–196 y 253–289 | Mismos coeficientes medios y N−1; DE directa no se sustituye por la aproximación suavizada `std.fd`. Media real verificada; pruebas adicionales contra `var.fd` y diagonal. |
| Derivadas | EDA 351–364, PDF p. 43; `deriv.fd()` reproyecta en la base original en fda 6.3.0 | Se conserva derivada analítica exacta. Se contrasta con `eval.fd(..., Lfdobj=1)`: `deriv.fd()` no da igualdad exacta en splines porque reproyecta. No confundir un fallo de ese contraste con un error de la fórmula propia. |
| Contexto de derivadas | Clase deriva lluvia acumulada en mm/péntada; aquí se hereda tasa mm/día | Nuestra derivada es cambio de tasa, mm/día². Positiva no significa automáticamente día de inicio; fechas no validadas. |
| **FM, corregida** | PDF p. 53 muestra d=1−abs(.5−F); `depth.FM` llama FM1 con scale=TRUE | `FM=2*promedio(d)−1` coincide en escala con el script; `FM_pdf=promedio(d)` deja visible la fórmula del PDF. F usa ≤, rango máximo en empates. No usar rangos medios para F. |
| **Recorte, corregido** | `depth.FM` selecciona profundidad ≥ cuantil .25; antes se tomaba ceil(.75N) | `centralidad_FM()` replica cuantil y empates; se exporta pertenencia. En estaciones entran 3.319 curvas, antes 3.320; satélite vigente 23.339 (histórico anterior: 32.472). |
| MBD con empates | Definición inclusiva del PDF vs fda 6.3.0, cuyo `fMBD` usa rangos medios | Se mantiene definición inclusiva exacta: 1−[L(L−1)+U(U−1)]/[N(N−1)]. Coincide con `fbplot` sin empates; **no** equivale al algoritmo de rangos medios con empates. Probada por enumeración de pares y con todas las curvas idénticas. |
| **Boxplot, corregido** | `fbplot` marca igualdad con los límites (≤/≥) y dibuja envolvente no atípica; antes solo < /> y límites | Se ajustan comparadores y se exportan/dibujan bigotes aparte de límites. Se contrasta `fbplot(depth=MBD_inclusiva)` sobre ambas fuentes completas. Si todas salen atípicas por región degenerada, bigotes NA: no inventar envolvente ni copiar el fallo de `apply` sobre cero curvas. |
| Datos preparados | Clase usa 10 sitios pareados X/Y de 72 péntadas; aquí miles de sitio-años no pareados, ahora en ventana común 1981–2003 | No reutilizar los datos de juguete ni igualar filas por posición. Sin coordenadas no hay correlación cruzada ni diferencias pareadas. |

El código de las librerías autorizadas se inspeccionó para resolver esas
diferencias, no para agregar métodos externos: `fda::fbplot`, `fda:::fMBD`
y `fda::deriv.fd` de la instalación local 6.3.0; `depth.FM` y `FM1` del código
de fda.usc (`R/depth.fdata.r` y `R/depth.mdata.R`, espejo CRAN en
`https://github.com/cran/fda.usc`, consultado 9/11/2026).
**fda.usc no está instalado:** se contrastó su implementación fuente, no se
afirma haber ejecutado su paquete. La prueba de FM usa ECDF independiente;
las de MBD y boxplot sí ejecutan fda local. El PDF y el script docente tienen
convenciones distintas que ahora están expresamente distinguidas.

### 8.3 Contexto de geoestadística y sesiones 3/4

Material localizado y efectivamente leído:

- `../analitica/Practic_Geostatistics.R`: EDA → transformación log10 →
  coordenadas → variograma empírico → modelo esférico/ajuste → grilla →
  kriging ordinario → varianza de predicción → LOO/RMSE. Usa **sp y gstat**.
  No se trasladan ni logZn ni rango 850 m ni kriging a las curvas de lluvia.
- `../analitica/Actividad_2.pdf`, **pp. 1–5**: CRS, distancia, estacionariedad,
  semivariograma de residuos y ajuste manual; **pp. 6–7**: ejemplo explícito
  de lluvia CHIRPS en Valle del Cauca, borde administrativo real, agregación
  **semanal**, distancias y kriging manual. Las excepciones geodata/chirps/terra
  pertenecen a ese taller, no autorizan construcción FDA con paquetes.
- `../analitica/Lecture_5_LatticeData (1).pdf`: portada dice **«Sesión 4»**;
  **p. 2** distingue geoestadística y datos de área, **pp. 8–20 y 25–28**
  vecindades/pesos, **pp. 30–38** Moran/Geary. El script Lattice aplica
  sf/spdep/tmap/RColorBrewer/dplyr a Carolina del Norte. Esa W espacial no
  es la matriz H de proyección temporal ni la Gram funcional del futuro modelo.
- `../analitica/Actividad_1.pdf`, pp. 1–3: dependencia espacial y GLS como
  contexto de los supuestos de varianza. No se añade GLS a esta entrega.
- `../analitica/IntroduciÃ³n al anÃ¡lisis de datos mapeados.pdf`, pp. 1–3:
  coordenadas, patrones de puntos, escala y efectos de borde; lectura contextual
  de esas páginas, no confundir patrones puntuales con kriging o FDA.

**Adaptación concreta al contexto de lluvia:** distinguir estación puntual de
píxel, no atribuir un corte por longitud a una cuenca oficial, no llamar
«Valle administrativo» a toda la caja de 39×37, y conservar como provisional
la máscara de tierra CHIRPS hasta confirmar georreferenciación/producto.
La covarianza funcional actual relaciona **días** dentro de una fuente; no
estima dependencia espacial entre estaciones. Miles de píxel-años no son miles
de réplicas necesariamente independientes. Los resúmenes son descriptivos,
no pruebas espaciales ni inferencia ENSO con independencia demostrada.

No hay archivos identificados inequívocamente como **«semana 3» de
geoestadística** en las dos carpetas buscadas. `S3_Taller1...` y `S4_Sesion1...`
son CNN, y `Taller 3.docx` trata Faster R-CNN: se descartaron como fuentes
espaciales. Sesión 4 sí está identificada; no se equipara automáticamente
«sesión» con «semana». Si se aludía a ejemplos adicionales de semanas 3/4,
falta el nombre o ruta precisa. No se rastrearon otras carpetas del Mac.

### 8.4 Evidencia histórica de la auditoría anterior (previa a fe de erratas)

**Los conteos y errores de esta subsección son históricos, no las salidas
vigentes.** Para resultados actuales usar §§3–6 y comparación de fe de erratas.

- Se añadió primero la prueba de la escala FM docente y **falló** con
  `max(abs(...$FM - (2 * fm_pdf - 1))) < 1e-12 is not TRUE`.
  Tras corregirla, pasan ECDF/empates, recorte por cuantil, límite exacto
  del boxplot, bigotes y las pruebas matemáticas anteriores.
- `Rscript test_fda_R.R`: **OK**, incluidos covarianza y derivada exacta
  contra fda, MBD sin empates y discrepancia conocida con empates.
  La reproyección de `deriv.fd` difiere 0,001450872 en el ejemplo aleatorio;
  no se relajó tolerancia para hacerla pasar: se cambió al contraste exacto.
- `Rscript verificar_flujo_R.R`: **OK**, rehace construcción+EDA y 10 controles,
  verifica las cuatro entradas intactas. Se regeneraron salidas R: K=6/11,
  170/0 atípicas, sin cambios en cobertura, coeficientes ni máscara.
  El EDA contrasta además los indicadores de atípicas contra `fbplot`
  usando la MBD inclusiva completa, no solo contra otra llamada propia.
- Contraste adicional sobre los objetos reales: error derivada analítica vs
  `eval.fd(Lfdobj=1)` **7,22e−16 / 7,99e−15**; vs `deriv.fd` reproyectada
  **0,02856774 / 0,06440525**. Covarianza vs `var.fd`: **4,27e−14** en ambas.
  Rangos FM en escala docente: **0,001062216–0,8916626** y
  **0,04352833–0,9554117**. Estas escalas sustituyen las anteriores del CSV.
- Actualizados: este documento de decisiones, trazabilidad, README y aviso
  de la guía Markdown. **`GUIA_RECONSTRUCCION_EDA.docx` fue inspeccionado y
  está obsoleto**: menciona secciones R antiguas, verificaciones Python y
  métodos históricos retirados. Se mantiene histórico; no se regeneró ni
  se presenta como Word actualizado. No es el Word normativo.
- Validación de trazabilidad: **31 filas**, ocho campos por fila, IDs únicos
  y nueve rutas de originales existentes. Última repetición de pruebas unitarias
  tras ajustar el caso de empate del recorte: **OK**. No se repitió el flujo
  pesado después porque los últimos cambios productivos fueron solo comentarios.
- SHA-256 Word normativo preservado: `5818aa2eb845c940589c872dad80cb2de1fe22c5e15c50b53517b578bd4e842c`.
  Guía Word histórica: `e143a545b85bc13fe2d68f5381cd91b702021d873e9d18ed847583bc73d8919a`.

**Conclusión verificable:** el núcleo MCO y los momentos siguen las fórmulas
docentes; esta revisión corrigió convenciones reales de profundidad/boxplot
y precisó una diferencia importante de derivación de librería. No se afirma
cumplimiento global ni copia literal de estructura, algoritmos o librerías.
En aquella auditoría seguían pendientes periodo completo/representatividad, fuentes/unidades,
coordenadas/altitud, fechas meteorológicas y confirmación docente de la
implementación manual de la base si exige una receta específica.

**Continuidad local:** código corregido y salidas regeneradas; próxima revisión
debe partir de esta sección y del CSV, no de las afirmaciones del Word histórico.
No hay nuevos paquetes instalados ni cambios en originales/notebooks/entradas.
El Git detectado está en la carpeta padre (la carpeta de este proyecto figura
sin seguimiento); no se creó commit. Consultas futuras de Git deben limitarse
con `-- .` para no mezclar el estado de otros trabajos de Descargas.

## 9. Decisión vigente por fe de erratas y continuidad

El correo docente, **según lo comunicado por la usuaria**, explica que no se
recuperó la base hasta 2013 y se trabajará con la serie histórica validada y
disponible hasta 2003. No se recibió aquí el archivo original del correo, su
fecha ni sus encabezados; no se inventa una cita textual. El inicio exacto
**1981** se obtuvo del Word nuevo, no se dedujo del correo.

Comparación completa de los 366 párrafos, incluidas tablas y rúbrica: cuatro
cambios textuales; véase `COMPARACION_WORD_FE_ERRATAS.md` para citas y revisión
de todos los criterios que permanecen. Las normas sobre construcción manual,
EDA con librería, métodos de clase, scripts únicos e informe siguen vigentes.

Se filtraron fuera 12.987 curvas satelitales de 2004–2012, sin alterar sus
archivos históricos. Estaciones sin exclusiones temporales; máscara en tierra
filtrada de 43.296 a 31.119. Se recalculó todo desde observaciones en R, no se
recortaron únicamente las tablas finales. K sigue 6/11 y atípicas 170/0, pero
medias, profundidades, momentos, ONI, grupos y negativos satelitales son nuevos.
Los 33.189 píxel–años son 1.443 × 23; los 31.119 de tierra son 1.353 × 23,
sin valores ni combinaciones ausentes en los preparados. Esto no certifica
integridad de los TIFF originales ni pertenencia a un borde administrativo.

Para el coordinador: continuar desde README, §§1–6 y 9 de este archivo y la
comparación normativa. No usar las cifras históricas de §8.4, guía Word,
notebooks o JSON Python como resultados actuales. `test_periodo_R.R` fue añadido;
`verificar_flujo_R.R` y `test_usabilidad_R.R` fueron actualizados sin relajar
tolerancias matemáticas. `proyecto_funcional.R` conserva las correcciones FM,
recorte con empates, MBD inclusiva y derivada analítica de la revisión anterior.
No hay paquetes adicionales, commits ni modificaciones a EduInteractiva.
La trazabilidad de esa revisión tenía **35 registros** (31 de métodos/auditoría anteriores
y N01–N03/U01 para norma, partición, comparación y ejecución visible), ocho
campos por registro e identificadores únicos.
`skill_resolution=injected`.

## 10. Recorrido pedagógico y Word vigente con resultados

**Petición aplicada:** ejecutar bloques de arriba abajo y ver su salida, no
limitarse a cargar definiciones y llamar al principal final. Se releen solo
los pasajes principales de construcción (73–105) y EDA (28–125): dimensiones,
observados, base, Phi, coeficientes, S y gráficos asociados. Se conserva ese
orden de aprendizaje, con construcción manual por la norma del proyecto.
No se atribuye la recursión propia ni cada línea al profesor.

### Diseño implementado

- `proyecto_funcional.R` es un **script secuencial normal**: preparación y
  definiciones completas primero; construcción y EDA después, con 13 pasos de
  cálculos directos. Ejecutar desde la primera línea con Cmd/Ctrl+Enter, o Source
  para todo. La ruta fija de esta entrega evita escribir salidas en otra carpeta.
- Filtro/observados → criterios → Phi/S → ajuste → máscara/RDS → momentos →
  dispersión → derivadas → profundidades → boxplot → contexto → diagnósticos.
  La selección de K calculada en el paso 3 se pasa al ajuste: no se repite.
- Environment conserva los objetos directos y los entornos `eda_est/eda_sat`
  con matrices/intermedios. `eda_momentos`, `eda_profundidad`, `eda_contexto`
  comparten fórmulas con `eda_fuente`; las entradas 03/04 ejecutan las secciones
  del mismo archivo mediante Source. El arnés separado `soporte_pruebas_R.R`
  importa solo la preparación para pruebas pequeñas. El usuario no lo necesita.
- Figuras individuales de observados, criterios, bases, bandas, residuos,
  media, varianza, DE, CV, covarianza, correlación, derivada, centralidad,
  boxplot, ONI y grupos: PNG180dpi y PDF más dispositivo interactivo. Mapas
  con barra de color numérica. Se conservan PDF compuestos compatibles;
  no son la única presentación. `ver_tabla()` abre como máximo la tabla pedida.
- No cambia ningún método estadístico, población, filtro, tolerancia ni escala.
  1981–2003 sigue antes de K/ajuste; construcción sigue sin paquetes.

### Documento para abrir y generación

**`INFORME_RECONSTRUCCION_EDA_ACTUAL.docx`**: 43 imágenes embebidas (36 figuras de resultados + 7 diagramas de decisión),
16 tablas (dos glosarios y 14 tablas de resultados). Incluye las figuras de
criterio por fuente (ECM, VC y VCG) y de dispersión (varianza, DE y CV), que
antes se generaban pero quedaban fuera del Word. Contexto/pregunta/objetivo,
fuentes/cobertura, vocabulario necesario, fórmulas/decisiones, valores extraídos
de CSV, figuras con pies e interpretaciones, conclusiones, limitaciones,
referencias y reproducción. Es guía de estas fases, no informe global completo.
La guía y generador antiguos se preservan como históricos.

`Rscript generar_diagramas.R` produce siete diagramas de decisión en R base
(flujo general, K, construcción, NA, centralidad, temporadas y comparación de
fuentes); cada uno muestra el camino elegido y los descartados con su razón.
`python3 generar_informe_actual.py` genera el Word y su Markdown desde los
CSV/PNG actuales, verifica ZIP/XML, imágenes, tablas y encabezados reales.
`resultados_R/manifiesto_informe_actual.json` registra SHA-256 de entradas del
informe y código. `resultados_R/manifiesto_figuras.csv` mapea cada figura del R
a su número en el Word y marca las compuestas duplicadas
(`seleccion_K_*`, `contexto_*`). No usa imágenes de la guía histórica ni inventa
prompts. El PDF requiere LibreOffice o Word (ausentes en este Mac), por lo que
no se incluye un PDF regenerado; si existe uno más reciente que el `.docx`,
PDFKit mediante `inspeccionar_informe.js` lo exporta a PNG y texto.
No se inspeccionó una sesión GUI de RStudio ni Word.

### Evidencia de esta intervención

1. `test_bloques_R.R` ejecuta selecciones de texto consecutivas del script real,
   desde la primera línea, sin precargar auxiliares, capturando la salida
   antes del siguiente paso. `verificacion_bloques.csv`: 14 registros con
   texto, duración y cantidad de objetos visibles; además exige selección
   previa al ajuste y media previa a profundidad, PNG y cifras esperadas.
2. `test_fda_R.R`, `test_periodo_R.R`, `test_usabilidad_R.R`: pasan. Se añadió
   cobertura PDF+PNG+pantalla y cierre del PNG ante error deliberado.
3. `verificar_flujo_R.R`: ejecuta Source normal y verifica objetos en la sesión;
   **12 controles OK**, K6/11, curvas4426/33189/31119, atípicas170/0 y entradas
   intactas. La modificación posterior de barras de color solo redibuja matrices
   exportadas; no precisa volver a seleccionar/ajustar todas las curvas.
4. Word: 43 imágenes/16 tablas/19 títulos H1; imágenes proporcionales y anchos
   de tablas/celdas explícitos. Metadatos satelitales, panel incompleto,
   negativos y temporadas pendientes declarados.
5. `python3 validar_informe_actual.py`: OK; cifras de cobertura y mínimos
   extraídas del Word contra CSV; imágenes embebidas byte a byte contra PNG;
   proporciones, hashes, entradas normativas/datos intactos y verificación de
   que toda figura exigida por la rúbrica esté incluida. Evidencia:
   `validacion_informe_actual.json`.

Trazabilidad ampliada: **38 registros** (se preservan los 35 y se agregan
U02–U04). Cambios exclusivamente en este proyecto externo; sin commits,
subdelegación ni cambios de código en EduInteractiva. `skill_resolution=injected`.

## 11. Corrección de usabilidad: script R normal

La versión vigente reemplaza el contenedor de bloques por código secuencial
directo. Primero están la ruta explícita de la entrega y las funciones
auxiliares completas; después las secciones de construcción y EDA. Se inicia
en la primera línea con Cmd/Ctrl+Enter, o Source ejecuta todo desde una sesión
nueva. Los objetos y salidas aparecen durante el recorrido. No se requiere
otra llamada de ejecución ni cargar previamente una sesión.

`soporte_pruebas_R.R` importa exclusivamente la preparación para las pruebas
matemáticas, temporales y sintéticas, sin recalcular los datos reales. Es una
herramienta del arnés, no una dependencia del archivo entregable al alumno.
Las entradas opcionales 03/04 ejecutan las secciones del mismo script con
Source; mantienen construcción y EDA por separado sin duplicar el pipeline.

Verificaciones reales de esta corrección:

- `test_fda_R.R`, `test_periodo_R.R`, `test_usabilidad_R.R`: OK.
- `test_bloques_R.R`: preparación y 13 selecciones normales consecutivas, OK;
  prueba selección anterior a ajuste y momentos anteriores a profundidad.
- `verificar_flujo_R.R`: Source normal, objetos visibles y 12 controles OK;
  K6/11, curvas4426/33189/31119, atípicas170/0; entradas intactas.
- `test_entradas_R.R`: cuatro procesos R --vanilla desde otra carpeta, OK:
  ejecución directa del script, Source normal con objetos en la sesión,
  entrada03 y entrada04. Todos los CSV deterministas son idénticos; ninguna
  salida se crea accidentalmente en la carpeta de partida. Los cuatro logs
  quedan en `resultados_R/prueba_*.log`.
- README, guía Markdown, trazabilidad U01/U02 y generador del Word actual se
  alinean con el uso secuencial. Los Word normativos y la guía Word histórica
  se conservan intactos. No se cambia ninguna fórmula ni resultado estadístico.

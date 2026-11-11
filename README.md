# Lluvia del Valle del Cauca: construcción y EDA funcional en R

**Entrada principal:** `proyecto_funcional.R`. Implementa únicamente construcción
y EDA, con los métodos de clase y los datos preparados existentes. La
construcción programa Φ y MCO sin librerías. El EDA permite `fda` y contrasta
la media manual contra `mean.fd`, como exige la rúbrica.

**Norma vigente:** `Proyecto_Final_Analitica_FeErratas.docx`, periodo exacto
**1981–2003 inclusive**. [Comparación íntegra y estado de requisitos](COMPARACION_WORD_FE_ERRATAS.md).
Se filtran ambas fuentes y la máscara antes de seleccionar K/ajustar; EDA,
ONI y grupos usan exclusivamente esa construcción. Los originales permanecen intactos.

## Ejecutar en RStudio (ver análisis, gráficos y tablas)

**Word vigente para abrir:** [INFORME_RECONSTRUCCION_EDA_ACTUAL.docx](INFORME_RECONSTRUCCION_EDA_ACTUAL.docx).
Incluye contexto, objetivos, glosario, métodos, resultados reales, **43 figuras** (36 de resultados + 7 diagramas de decisión),
16 tablas (incluidos dos glosarios), conclusiones y límites. Las figuras cubren
ECM, VC y VCG por fuente, varianza, DE y CV, residuos, covarianza/correlación,
derivadas, profundidades, boxplots y contexto ONI/grupos. La correspondencia
figura del R → número en el Word está en `resultados_R/manifiesto_figuras.csv`.
El PDF se obtiene exportando el `.docx` (requiere LibreOffice o Word; este Mac
no lo tiene instalado, por eso no se incluye un PDF regenerado).
La antigua `GUIA_RECONSTRUCCION_EDA.docx` se conserva como histórica.

1. Abra `proyecto_funcional.R` en el **editor** de RStudio.
2. Ejecute desde la **primera línea**, con **Cmd+Enter** (Mac) o **Ctrl+Enter**,
   avanzando de arriba abajo como en clase. La preparación fija la carpeta de
   esta entrega y define los auxiliares; RStudio ejecuta cada definición de función
   completa como una unidad. Después vienen **CONSTRUCCIÓN FUNCIONAL** y **EDA FUNCIONAL**.
3. También puede pulsar **Source** para ejecutar todo de una vez, desde una sesión
   nueva. No hace falta cargar nada antes ni escribir otra llamada en consola.
4. Los cálculos imprimen su tabla o resumen y dibujan sus figuras **antes** de
   continuar. Mire **Console** y **Plots**, con flechas de historial y **Zoom**.
   Los objetos quedan en **Environment**: `est`, `sat`, `seleccion_est`, `Phi_est`,
   `S_est`, `e`, `s`, `eda_est` y `eda_sat`. Los dos últimos son entornos que puede
   desplegar para inspeccionar `C`, `P`, `Y`, `mu`, `covar`, `pr` o `bx`.

| Sección | Cálculo y salida inmediata |
|---|---|
| Preparación | Ruta del proyecto y definiciones completas de auxiliares |
| 1 | Lectura, filtro 1981–2003, conteos y cobertura; aún no ajusta |
| 2 | Primeras observaciones y gráficos discretos de ambas fuentes |
| 3 | ECM/VC/VCG y mínimos; gráficos individuales por criterio/fuente |
| 4 | Phi, mapa, S y controles del ejemplo; bases por fuente |
| 5 | Ajustes con la selección ya calculada; pesos, controles, bandas y residuos |
| 6 | Máscara por claves, RDS y cobertura anual completa en CSV |
| 7 | Evaluación de curvas, momentos/media y contraste con `fda` |
| 8 | Varianza, DE, CV y mapas de covarianza/correlación con escala de color |
| 9 | Derivada de la media y tabla de cambio de tasa |
| 10 | Profundidades FM/MBD y centralidad robusta |
| 11 | Boxplots individuales y claves de atípicas |
| 12 | Medias y tamaños por ONI/grupos; contexto, no causalidad |
| 13 | Diagnósticos, sesión y ruta para regenerar el Word |

El archivo contiene código R secuencial: las asignaciones y cálculos están escritos
directamente después de los auxiliares. **Source ejecuta el análisis completo**;
las líneas o selecciones consecutivas permiten estudiarlo al paso. Los bucles y
las definiciones se ejecutan completos, como en cualquier script R de clase.

Cada figura se exporta en **PNG de 180 dpi y PDF**; el EDA ya no se presenta en
seis paneles minúsculos. Se conservan algunos PDF compuestos para compatibilidad.
No se abren decenas de ventanas View; una tabla completa a petición:

```r
ver_tabla("cobertura.csv")
ver_tabla("resumen_estaciones.csv")
ver_tabla("diagnosticos_EDA.csv")
# Solo consola, sin abrir el visor:
ver_tabla("seleccion_K_estaciones.csv", abrir = FALSE)
```

**Ejecutar todo desde consola** (equivale a Source):

```r
source("proyecto_funcional.R")
```

En sesiones interactivas las figuras se muestran por defecto. En batch se exportan
sin abrir dispositivos visibles. Ningún paquete se instala al analizar.

## Ejecutar con Rscript y verificar

Desde esta carpeta, sin cambiar rutas dentro de los scripts:

```bash
Rscript proyecto_funcional.R
Rscript test_usabilidad_R.R
Rscript test_bloques_R.R
Rscript test_fda_R.R
Rscript test_periodo_R.R
Rscript verificar_flujo_R.R
Rscript test_entradas_R.R
python3 generar_informe_actual.py
```

Con `Rscript` se ejecuta todo el script desde el principio,
sin dispositivo visible ni `View`: consulte consola, CSV y PDF en `resultados_R/`.
La prueba de usabilidad utiliza dispositivos PDF/PNG y datos sintéticos pequeños.
El verificador integral ejecuta Source normal y verifica los resultados reales;
el generador de Word solo consume los resultados ya calculados.
`test_entradas_R.R` comprueba ejecución directa, Source y entradas 03/04 en
procesos R limpios desde otra carpeta; exige CSV idénticos y objetos visibles
después de Source. Es una comprobación completa adicional, no un paso de uso.
Para repetir una fase: `Rscript 03_construccion_funcional.R` y después
`Rscript 04_eda_funcional.R` (reutiliza el RDS sin volver a ajustar).
Para la entrega de estas dos fases basta el **único** `proyecto_funcional.R`;
03 y 04 son entradas opcionales que ejecutan secciones del mismo archivo.
`soporte_pruebas_R.R` es exclusivo del arnés: importa definiciones para las pruebas
pequeñas sin ejecutar de nuevo el análisis; no se necesita para usar el script.

### Regenerar y comprobar el Word

Después de ejecutar el recorrido, las tablas y figuras son la fuente del Word;
el generador no vuelve a seleccionar K ni ajusta modelos. Desde esta carpeta:

```bash
Rscript generar_diagramas.R            # 7 diagramas de decisión (R base)
python3 generar_informe_actual.py
python3 validar_informe_actual.py
# PDF (opcional, requiere LibreOffice o Word instalado):
soffice --headless --convert-to pdf --outdir resultados_R INFORME_RECONSTRUCCION_EDA_ACTUAL.docx
# Inspección reproducible en macOS con PDFKit, sin instalar paquetes:
osascript -l JavaScript inspeccionar_informe.js "$PWD/resultados_R/INFORME_RECONSTRUCCION_EDA_ACTUAL.pdf" "$PWD/resultados_R"
```

`validar_informe_actual.py` contrasta cifras del Word con CSV, bytes de las **36**
imágenes con sus PNG, proporciones, manifiesto de hashes, siete entradas
preservadas y 38 registros de trazabilidad. Exige además que toda figura de la
rúbrica (ECM/VC/VCG, varianza/DE/CV, residuos por fuente) esté **incluida** en el
Word según `resultados_R/manifiesto_figuras.csv`. Si el PDF existe y es más
reciente que el `.docx`, valida también sus páginas; si no hay LibreOffice, lo
reporta como `no_regenerado_sin_LibreOffice` sin fallar.
El verificador integral vuelve a calcular el batch; no necesita repetirlo para
un cambio de redacción o de dibujo sin cambios estadísticos.

### Dependencia autorizada para EDA

Se instaló `fda` 6.3.0 y sus dependencias requeridas en `.Rlib/` (ignorada por
Git), porque la verificación contra la función de librería es un requisito
expreso. Ningún paquete se instala al ejecutar el análisis. En otro equipo:

```r
dir.create(".Rlib", showWarnings = FALSE)
install.packages("fda", lib = ".Rlib", repos = "https://cloud.r-project.org")
```

La construcción funciona sin esa dependencia. Si falta `fda`, el EDA guarda
explícitamente `PENDIENTE` en el contraste de librería; el verificador integral
lo rechaza. Versiones reales: `resultados_R/sessionInfo.txt`.

## Entradas preservadas y salidas

| Ruta | Función |
|---|---|
| `fda_fase1/curvas_estaciones_mm.csv.gz` | 4.426 curvas observadas preparadas, NA donde faltan datos; IDs y años. |
| `fda_fase1/curvas_satelite_mm.csv.gz` | 46.176 curvas observadas preparadas; tiempos pentadales en encabezados. |
| `fda_fase2/coef_satelite_tierra.csv.gz` | **Solo sus claves** definen la máscara de tierra heredada: 43.296 curvas. Sus coeficientes Python no se usan en el EDA vigente. |
| `eda_fase0/oni_anual_fase.csv` | Clasificación descriptiva por ONI anual heredado. |
| `resultados_R/` | Resultados nuevos y separados: selección K, Φ, coeficientes, máscaras observadas en RDS, sigma², bandas, momentos residuales, EDA, profundidades, atípicas, grupos/ONI, figuras y verificación. |
| `resultados_R/cobertura_anual_*.csv` | Los 23 años explícitos: curvas, sitios, ausencias, completas/incompletas y NA. |
| `resultados_R/sitio_anio_ausentes_*.csv` | Claves de combinaciones sin curva respecto del catálogo observado en la ventana. No implica que cada estación debiera existir todos los años. |
| `resultados_R/filtrado_periodo.csv`, `ONI_periodo.csv` | Conteos antes/después del filtro y clasificación ONI limitada a la ventana vigente. |
| `archivo_historico_notebooks/` | Ocho notebooks inspeccionados y preservados, incluida toda preparación y descarga previa. No son el flujo principal. |
| `REQUISITOS_Y_ENTREGA_R.md` | Citas exactas, fuentes, decisiones, resultados y pendientes de esta entrega. |
| `trazabilidad_metodos_clase.csv` | Métodos R01–R08/E01–E07 y funciones R vigentes. |
| `INFORME_RECONSTRUCCION_EDA_ACTUAL.docx` | **Word vigente** con resultados y figuras de R; abrir este. |
| `generar_informe_actual.py` | Generación reproducible desde CSV/PNG, validación ZIP/XML/tablas/imágenes y manifiesto SHA-256. Requiere pandoc. |
| `INFORME_RECONSTRUCCION_EDA_ACTUAL.md` | Fuente generada y legible del Word vigente; se regenera, editar el generador para cambios permanentes. |
| `GUIA_RECONSTRUCCION_EDA.md`, `.docx` | Guía y bitácora **históricas**, preservadas; no son el informe vigente. |
| `resultados_R/verificacion_bloques.csv` | Salida capturada de cada bloque antes de ejecutar el siguiente y objetos disponibles. |

`fda_fase1/`, `fda_fase2/`, `eda_fase0/` y `datos_web/` se conservan como
insumos/resultados históricos: el flujo no los sobrescribe. Las observaciones
y las evaluaciones ajustadas siguen separadas; NA no se cambia por cero.

## Estado verificado y límites

**Decisiones vigentes:** [REQUISITOS_Y_ENTREGA_R.md, fe de erratas y continuidad](REQUISITOS_Y_ENTREGA_R.md#9-decisión-vigente-por-fe-de-erratas-y-continuidad).
Incluye inventario de los nueve scripts docentes leídos con SHA-256, orden y
notación de clase, contexto de geoestadística/sesión 4 y desviaciones explícitas.
`GUIA_RECONSTRUCCION_EDA.docx` sigue siendo histórica. El Word vigente es
**`INFORME_RECONSTRUCCION_EDA_ACTUAL.docx`**, generado de la ejecución actual.

En la revisión se ajustaron la escala `FM` a `depth.FM` (se conserva `FM_pdf`),
el recorte por cuantil con sus empates y la igualdad en los límites del boxplot.
Ahora se distinguen límites de detección y bigotes no atípicos. La MBD sigue
la definición inclusiva, distinta de los rangos medios de fda cuando hay empates.
La derivada analítica se verifica con `eval.fd(Lfdobj=1)`; `deriv.fd` reproyecta
y no es idéntica. Los detalles y pruebas están en el documento de decisiones.

Se ejecutó el flujo R completo desde las matrices preparadas. K por VC y VCG:
**6 estaciones, 11 satélite**. Media contra `fda::mean.fd`: errores máximos
1,78e-14 y 9,41e-14. El verificador integral pasó **12 controles**, además de las
pruebas matemáticas/temporales, y comprobó integridad de cuatro entradas de datos
y los tres Word (ambas consignas y guía histórica). Atípicas MBD: **170/0**.
Se verificó el recorrido de 14 bloques con salidas intermedias y la visualización
PDF + PNG + dispositivo activo, incluyendo limpieza por error. El Word vigente
tiene 43 imágenes embebidas (36 figuras de resultados + 7 diagramas de decisión) y 16
tablas; la correspondencia figura a figura está en `manifiesto_figuras.csv`.
No se afirma inspección de la GUI real de RStudio.

**No equivale a reproducir la depuración desde crudos**: los Excel/CSV de
estaciones y TIFF originales no están en esta raíz. La selección preliminar,
la normalización temporal y las unidades satelitales se heredan y están
identificadas como tales. Construcción y EDA vigentes: **1981–2003** en ambas
fuentes. Se excluyen del cálculo 12.987 curvas satelitales de 2004–2012;
la máscara pasa de 43.296 a 31.119 claves. Los insumos históricos no se modifican.
Se analizan las 4.426 curvas de estaciones (3.761 completas, 665 incompletas) y
31.119 curvas satelitales terrestres (todas completas); la selección K satelital
y el ajuste incluyen los 33.189 píxel–años de la caja antes de aplicar la máscara.
Hay datos en cada uno de los 23 años, pero el panel terrestre no es balanceado:
238 estaciones observadas, 1.048 combinaciones estación–año ausentes y 2.871 NA
diarios en las curvas presentes. El satélite tiene 1.443 píxeles/año (1.353 en
tierra), sin combinaciones ni valores ausentes en las matrices preparadas.
Esta cobertura no demuestra representatividad espacial ni valida los crudos.
La necesidad de llegar a 2013 quedó sustituida por la fe de erratas; no es un
incumplimiento pendiente. El correo docente referido por la usuaria explica que
no se recuperó la base hasta 2013 y se usará la serie disponible validada hasta 2003.

La partición de la fase futura de modelado ahora es **1981–1995 / 1996–2000 /
2001–2003**, con estaciones reservadas por el docente. K actual es descriptivo
del EDA completo; no reutilizarlo como selección sin fuga para evaluar redes.
El EDA rechaza RDS históricos sin ventana o con periodo incorrecto y pide
reconstruir antes de continuar.

Se calculan derivadas y se interpretan aumentos/descensos; **inicio y fin
meteorológicos no están validados**. Se retiraron del flujo las fechas
heurísticas y la banda estacional sin respaldo localizado. Hay valores
negativos en el MCO sin restricciones: se reportan, no se recortan ocultamente.

La comparación pareada, regresión y redes quedan fuera de esta intervención.
El proyecto global exige también Python; este alcance solicitado solo cubre
las fases funcionales en R. El informe final y la declaración de uso de IA
exigidos por el Word todavía deben prepararse por el grupo.

Los verificadores `verificar_fda_clase.py` / `test_fda_clase.py` y su JSON
pertenecen al historial Python; su localizador admite notebooks archivados.
No sustituyen `verificar_flujo_R.R` ni acreditan la ejecución nueva.

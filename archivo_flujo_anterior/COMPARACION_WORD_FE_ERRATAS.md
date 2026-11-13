# Fe de erratas: comparación normativa completa y aplicación en R

## 1. Autoridad y método de comparación

Norma vigente: **`Proyecto_Final_Analitica_FeErratas.docx`**. Documento anterior:
`Proyecto_Final_Analitica.docx`. Se leyeron completos mediante Pandoc y extracción
de `word/document.xml`, uniendo los runs de cada párrafo, incluidas todas las
celdas de tablas. Ambos tienen **366 párrafos** (incluidos vacíos): **cuatro
párrafos con texto distinto y 362 con el mismo texto**. Se citan secciones y
celdas, sin inventar páginas de Word. Esta comparación es textual, no afirma
identidad de formato o metadatos.

El inventario ZIP no tiene partes añadidas/retiradas. Cambian `word/document.xml`,
`settings.xml`, `numbering.xml`, `styles.xml`, `fontTable.xml` y propiedades
`docProps/core.xml`, `app.xml`. El resto de partes es idéntico; no hay otros
cambios de texto normativo detectados en el cuerpo/tablas. Los cambios de
estilos/configuración/propiedades no se interpretan como nuevos requisitos.

SHA-256, originales preservados:

- Anterior: `5818aa2eb845c940589c872dad80cb2de1fe22c5e15c50b53517b578bd4e842c`.
- Fe de erratas: `0e91f06439d1f4d29daf2081a5d549891a620cdccaee46b04799da87bf247c50`.

**Correo docente referido por la usuaria:** no se recuperó la base hasta 2013;
se analizará la serie histórica validada/disponible hasta 2003. Es contexto
aportado por la usuaria, no una transcripción de un correo inspeccionado. No
se dispone aquí de fecha, remitente ni archivo del correo. El **inicio 1981**
proviene directamente de introducción/objetivo del Word nuevo.

## 2. Los cuatro cambios exactos

| Ubicación | Texto anterior | Texto vigente | Aplicación |
|---|---|---|---|
| Introducción, primera oración | «en el periodo 1981 a 2013» | «en el periodo 1981 a 2003» | Ambas fuentes: 1981–2003 inclusive, 23 años. |
| Objetivo general | «en el periodo 1981-2013» | «en el periodo 1981-2003» | Ventana del análisis y modelado. |
| Rúbrica 2, tope EDA | «Analizar solo un subconjunto de años o estaciones sin justificar que es representativo del periodo 1981–2013: máximo 2.» | «Analizar solo un subconjunto de años o estaciones sin justificar que es representativo del periodo 1981–2003: máximo 2.» | Se incluyen todas las curvas preparadas en ventana; se informa cobertura real y limitaciones de representatividad. |
| Rúbrica 4.4, primera exigencia | «Partición 1981–2005 / 2006–2009 / 2010–2013 y estaciones reservadas por el docente, sin fuga de información.» | «Partición 1981–1995 / 1996–2000 / 2001–2003 y estaciones reservadas por el docente, sin fuga de información.» | Nueva partición de entrenamiento/validación/prueba para la fase futura; no hay entrenamiento en esta entrega. |

No cambió solo el año final: **también cambió la partición de evaluación**.
El inicio no pasó a 1980, 1982 ni otro año. No corresponde seguir señalando
«falta llegar a 2013» como incumplimiento de la consigna vigente.

## 3. Requisitos revisados que permanecen iguales

### Reglas generales, objetivos y fuentes

- Reglas 1–5: condiciones del curso; PDF con nombre/asunto «Proyecto final
  AnalíticaDeDatos»; todos los integrantes; código único en R y Python,
  separado y reproducible; solo scripts/procedimientos/metodologías de clase,
  librerías adicionales exclusivamente donde se autoricen.
- Reglas 6–10: tablas/gráficos/indicadores e interpretación contextual; máximo
  20 páginas incluidas referencias, una columna; no manuscrito, solo Colab
  ni R Markdown; declaración de herramienta IA y prompts en Word aparte.
- Reglas 11–19: fecha/hora y plataforma institucional, no copia sin citar,
  asistencia, 40 minutos de presentación más preguntas, todos deben conocer
  todas las fases, consecuencias individuales de respuestas/asistencia,
  entrega de PDF y scripts antes de sustentación, sin versiones posteriores.
  La regla 16 sigue diciendo «reducción de -0.2» y la nota final de sustentación
  sigue diciendo que el punto falla si no se explica: esta diferencia interna
  ya existía y la fe de erratas no la resuelve.
- Secuencia: construcción → EDA → comparación funcional. Objetivos 1–7 sin
  cambios: construcción «paso a paso (sin usar librerias)»; EDA y comparación
  «pueden usar librerias»; KAN y referencias; interpretación de activaciones;
  red funcional/FPCA y ablación KAN/MLP/lineal.
- Fuentes sin cambios: estaciones diarias puntuales y satélite pentadal en
  malla regular; cada localización–año es una curva.

### Red y evaluación: especificaciones detalladas sin cambios adicionales

Codificador B-spline con **6–8 neuronas**, integral `z = A W c` y matriz de
Gram; altitud y ONI anual como covariables. Cabeza con **M=4–6** coordenadas;
FPCA de estaciones ajustado solo en entrenamiento y congelado, varianza
explicada reportada. Pérdida L² funcional con W y rugosidad de los pesos;
modelado en raíz de P o log(1+P), reconstrucción no negativa en mm. Cabezas
lineal, KAN B-spline/SiLU y MLP, con equivalencia función–función demostrada
y verificada. KAN/MLP con parámetros ±10 %, tres semillas, mismo codificador
y decodificador. Depuración, cobertura, justificación de datos agregados o
eliminados y trazabilidad de estaciones/pares por paso. Comparación con
satélite sin corregir; L², KGE, sesgo anual, fechas de dos picos, bootstrap y
ENSO; pesos/activaciones interpretables frente a cabeza lineal. **PyTorch y
pykan siguen autorizados en ese apartado**, no para construir las curvas R.

### Rúbrica completa: pesos, exigencias y topes

| Criterio | Exigencia que permanece (además de los cuatro cambios de §2) | Estado de esta entrega |
|---|---|---|
| 1, Construcción **10 %** | Φ/MCO manuales; ECM/VC/VCG con gráfico y elección; S idempotente y traza K o alternativa penalizada; covarianza residual/bandas puntuales; dominio común diario/pentadal. Construir con librería: máximo 2. | Implementado desde preparados con R estándar; no se reproduce limpieza cruda. Banda puntual condicional, un ejemplo identificado por fuente. |
| 2, EDA **15 %** | Media manual contrastada con librería; varianza/DE; covarianza/correlación; derivadas e inicio/fin de lluvias; FM/MBD, boxplot y atípicas con claves; climatología regional. | Momentos/derivadas/profundidades/tablas/figuras y contraste fda realizados. Inicio/fin meteorológicos y atribución orográfica no validados. Representatividad no demostrada por tener todos los años. |
| 3, Comparación **15 %** | Pares estación–píxel, coordenadas/regla/píxeles compartidos; diferencia y banda ±2EE y tramos sin cero; correlación cruzada/diagonal; covarianzas; altitud/subregión/ENSO; sesgo climatológico. Sin pares documentados: máximo 1. | Pendiente; grupos descriptivos no sustituyen pares. |
| 4.1, Depuración **6 %** | Imposibles/digitación, acumulados mensuales en un día, ceros falsos, faltantes/años incompletos; cobertura justificada; agregar/eliminar justificado; tabla de estaciones/pares por paso; sin información de prueba. | Filtro temporal y diagnóstico de preparados reproducibles; auditoría cruda/calidad y trazabilidad de limpieza completa pendientes. |
| 4.2, Red **10 %** | Gram por cuadratura, FPCA solo entrenamiento, pérdida/rugosidad; equivalencia lineal error ≤1e-4; no negatividad y % recorte. MLP de péntadas presentado como funcional: máximo 1. | Pendiente; no se presenta este EDA como red. |
| 4.3, Ablación **11 %** | KAN/MLP/lineal mismos bloques, parámetros con sum(p.numel()) ±10 %, tres semillas, presupuesto/validación iguales, KAN sintética RMSE ≤1e-2; superioridad solo si IC bootstrap 95 % de diferencia excluye cero. Hiperparámetros con prueba: máximo 1. | Pendiente. |
| 4.4, Evaluación **8 %** | Estaciones reservadas; baseline satelital; L² relativo/KGE/sesgo anual/error de dos picos; IC 95 % con ≥1000 bootstrap y Wilcoxon pareada por estación; ENSO con fuente ONI; CSV station_id,year,day,value_mm antes de sustentación. Fuga: 4.3/4.4 máximo 2. | Partición nueva documentada, sin evaluación predictiva. ONI anual actual es descriptivo, no validación de episodios oficiales. |
| 4.5, Interpretación **10 %** | Pesos w_k/ventanas y activaciones primera capa KAN, contraste lineal, ajuste simbólico de al menos una activación con R², climatología coherente y peores casos. Contradecir gráficas: máximo 2. | Pendiente de modelos. |
| 5, Informe **10 %** | PDF una columna/20 páginas con referencias e integrantes, tablas/figuras numeradas/tituladas/citadas/interpretadas, fuentes de métodos/datos/ONI, redacción propia, Word IA. Solo primeras 20 páginas; dos columnas/manuscrito/R Markdown/solo Colab no aceptados; IA no declarada máximo 2. | Informe final y anexo IA pendientes; estos Markdown son documentación técnica. |
| 6, Scripts **5 %** | Un archivo R y uno Python separados del informe, ejecución completa con rutas relativas/semillas/versiones, métodos de clase y localización de tablas/figuras. Repositorio y comprimido en plataforma. | Un R para construcción/EDA, con secciones separadas; Python/global/compresión final pendientes. Wrappers y tests son auxiliares. |

Se revisaron también la escala Excelente/Bueno/Aceptable/Insuficiente y el
desglose 4.1–4.5 de notas 0–5: sus textos no cambiaron. Ítem 4 sigue sumando
**45 %**, calculado `Σ(nota × peso)/45`, asignando nota inferior entre dos
descripciones. Los pesos globales siguen **10+15+15+45+10+5=100 %**.

## 4. Implementación y resultados vigentes

- Un mismo `proyecto_funcional.R`, con `# CONSTRUCCIÓN FUNCIONAL ----` y
  `# EDA FUNCIONAL ----`; ahora BLOQUES 0–13 de cálculo al principio, seleccionables
  con Cmd/Ctrl+Enter. El batch evalúa los mismos bloques; detalles en decisiones §10.
- Filtro por nombres de claves `STATION_ID,YEAR` / `anio,fila,col`, validado
  antes de K y ajuste, incluida máscara de tierra. No se filtran solo salidas.
- Satélite: **46.176 → 33.189** curvas; **12.987 excluidas** (2004–2012).
  Máscara: **43.296 → 31.119**; 12.177 claves fuera de ventana. Estaciones:
  **4.426 → 4.426**, 3.761 completas y 665 incompletas.
- Cobertura: 23 años con curvas en ambas fuentes; catálogo de **238 estaciones**,
  **1.048 combinaciones ausentes** respecto de 238×23 y **2.871 NA diarios** en
  curvas presentes. Estaciones anuales entre 155 y 218. Satélite **1.443 píxeles
  por año**, **1.353 en tierra**, sin combinaciones/valores ausentes en preparados.
  Estos paneles no prueban existencia operativa de cada estación todos los años
  ni detectan sitios totalmente excluidos por la preparación histórica.
- K por VC/VCG sigue **6 estaciones / 11 satélite**, nuevamente calculado.
  ECM mínimo en extremos de búsqueda (100/40); no se confunde con K elegido.
  Atípicas MBD **170 / 0** (estaciones/tierra). No se eliminan las atípicas.
- Media vs `fda::mean.fd`: **1,78e-14 / 9,41e-14**. Satélite cambia a mínimo
  **−1,195546**, **408** puntos negativos, máximo DE día **310** y correlación
  abril–octubre **0,824418**. Estaciones: mínimo −11,624358, 24.361 negativos,
  máximo DE día 365, correlación 0,808861. El MCO no se recorta ocultamente.
- Derivada media: máximos dominantes de estaciones entre días 114–115 y
  316–317; satélite entre 134–135 y 307–308. Compatible con régimen bimodal,
  no fechas meteorológicas certificadas de inicio/fin. Grupos longitudinales
  sugieren gradiente húmedo occidental; faltan coordenadas/altitud validadas
  para atribución orográfica. Las diferencias ONI no son sesgo pareado ENSO.
- RDS declara la ventana. EDA rechaza RDS antiguo/incompatible, valida años y
  ONI completo sin duplicados. K global actual describe el EDA; para modelos
  habrá que seleccionar transformaciones/FPCA/hiperparámetros solo con entrenamiento.

## 5. Verificación y continuidad

Ejecutados: `Rscript test_periodo_R.R` (falló antes de implementación y pasó
después), `Rscript test_fda_R.R`, `Rscript test_usabilidad_R.R`,
`Rscript proyecto_funcional.R` y `Rscript verificar_flujo_R.R`.
Último verificador: **12 controles OK**, coeficientes de todas las curvas
contra QR independiente, momentos, H, atípicas contra fbplot en EDA, claves,
filtro independiente, máscara, cobertura/ONI y cuatro datos+tres Word intactos
por MD5. Se mantuvieron tolerancias y fórmulas matemáticas previas.

Una comprobación inicial de máscara falló exclusivamente por `row.names`
distintos entre CSV; se diagnosticó con `all.equal`, y se contrastan ahora
nombres/tipos/valores/orden de cada columna, sin exigir índices de filas
históricas iguales. No se relajó igualdad de claves ni tolerancias numéricas.
La prueba UI verifica PDF y dispositivo PNG activo, no una sesión GUI RStudio.

Salidas actuales regeneradas en `resultados_R/`: PDF y PNG individuales, CSV de construcción,
EDA, grupos/ONI, filtro/cobertura/ausencias, RDS, sesión e integridad/verificación.
README tiene comandos RStudio para ver tablas/Plots y ejecutar fases separadas.
Decisiones metodológicas detalladas: `REQUISITOS_Y_ENTREGA_R.md` §§1–6 y 9.
§8.4, notebooks, JSON Python y la antigua guía Word son evidencia **histórica**.
**Word vigente: `INFORME_RECONSTRUCCION_EDA_ACTUAL.docx`**, 43 imágenes (36 figuras + 7 diagramas) y 16
tablas/glosarios, generado desde salidas R. Incluye las figuras de criterio
(ECM/VC/VCG) y dispersión (varianza/DE/CV) por fuente; el mapeo figura a figura
está en `resultados_R/manifiesto_figuras.csv`. El PDF se exporta desde el `.docx`
con LibreOffice/Word (no incluido: ausentes en este Mac).

Pendientes reales: datos crudos/calidad y representatividad, unidad/producto/
calendario satelital, coordenadas/altitud/pares, inicio/fin meteorológicos,
comparación y modelos/evaluación, informe PDF y anexo IA. No está pendiente
ampliar a 2013 por la consigna nueva. Métodos originales de clase conservados;
no se añadieron paquetes ni métodos externos. Se creó el Word vigente de estas
dos fases; informe global y anexo de prompts siguen pendientes. Sin commits.

Continuidad guardada localmente en este proyecto externo, sin guardar su
metodología bajo la memoria de EduInteractiva. `skill_resolution=injected`.

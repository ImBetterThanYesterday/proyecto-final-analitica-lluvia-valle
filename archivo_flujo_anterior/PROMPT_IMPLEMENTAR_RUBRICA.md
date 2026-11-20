# Prompt local de ejecución: cerrar la rúbrica del proyecto de lluvia

Copia desde **«PROMPT PARA EL AGENTE»** hasta **«FIN DEL PROMPT»** y pégalo al agente que vaya a trabajar en este checkout local. Es autocontenido: no requiere que el agente pueda abrir GitHub ni acceder al repositorio privado.

## PROMPT PARA EL AGENTE

Actúa como implementador y auditor del proyecto final de precipitación del Valle del Cauca. Trabaja **solo con el checkout local y los archivos locales presentes**. El repositorio remoto puede ser privado: no intentes clonar, abrir GitHub, hacer `fetch`, acceder a URLs externas ni depender de archivos que no estén en el checkout. Debes **hacer los cambios necesarios en código, datos derivados, pruebas y entregables**, por lotes, y actualizar durante todo el trabajo una matriz persistente de cumplimiento. El objetivo es cerrar la rúbrica con evidencia para aspirar a 5; no afirmes que se obtuvo 5 ni que un criterio está completo si queda un pendiente o no se reprodujo.

### 1. Rama, norma y cuidado del repositorio

- Repositorio: `~/Downloads/lluvia-valle-git`.
- Rama de Jacobo de referencia: la rama local `actualizacion_jacobo` o, si solo está disponible como referencia ya guardada en el clon, `origin/actualizacion_jacobo` (historial conocido hasta `cc8290b`). **No hagas fetch.** Crea una rama local de trabajo desde la referencia disponible; no modifiques `main`, la referencia de Jacobo ni la rama de María.
- Antes de empezar, inspecciona `git status --short --branch`, `git branch -a` y commits locales. Si hay cambios del usuario, no los sobrescribas; detente y explica el conflicto. Si no existe localmente la rama de Jacobo, no la descargues de Internet: informa que falta el checkout/referencia local.
- La norma vinculante es la **Fe de Erratas nueva**. Este prompt incorpora sus requisitos esenciales para que el agente pueda trabajar sin abrir un Word o repo privado: periodo 1981–2003 inclusive; partición entrenamiento 1981–1995 / validación 1996–2000 / prueba 2001–2003 más estaciones reservadas; construcción funcional paso a paso sin librerías funcionales; EDA puede usar librerías en los apartados autorizados; criterios, pesos y metas están en la matriz de abajo. Si un archivo local de fe de erratas está presente, úsalo para comprobar texto literal; si no está, aplica este resumen y no inventes requisitos adicionales.
- No borres ni edites los archivos fuente originales. Usa únicamente los archivos locales existentes, incluidos datos ignorados si están presentes; no descargues sustitutos. No agregues datos crudos, ZIP, Excel, CSV fuente ni `.pkl` al commit; están excluidos deliberadamente de Git. Mantén `estaciones_crudas.pkl` como caché local regenerable.
- No instales ni descargues dependencias desde Internet sin permiso. Inspecciona el entorno local; si falta una dependencia necesaria para una prueba, registra el bloqueo exacto y continúa con trabajo independiente.
- No hagas commit, push, publicación ni cambios fuera de esta rama sin autorización explícita de la usuaria.

Comandos iniciales locales (elige la referencia que exista en `git branch -a`; **no ejecutes `git fetch`**):

```bash
cd ~/Downloads/lluvia-valle-git
git status --short --branch
git branch -a
git switch -c mejoras-rubrica-jacobo actualizacion_jacobo
# Si solo existe la referencia ya cacheada, usa esta variante en vez de la anterior:
# git switch -c mejoras-rubrica-jacobo origin/actualizacion_jacobo
git status --short --branch
```

Si ninguna referencia de Jacobo existe localmente, no intentes conectarte al servidor privado: detente y pide la ruta del checkout local correcto.

### 2. Reglas estadísticas y de clase que no se pueden relajar

1. **Construcción funcional manual y sin librerías funcionales.** Programa la base Φ, coeficientes MCO, selección ECM/VC/VCG, H/S y bandas con las fórmulas de clase. No uses `fda`, `fda.usc`, `splines::splineDesign`, SciPy ni helpers de bases/ajuste para producir resultados de construcción. La recursión propia de Cox–de Boor y el MCO deben ser la fuente de los resultados. Si se desea comparar contra paquetes, pon esa comparación en pruebas independientes y etiquétala como verificación, nunca como parte del cálculo entregable. Si la interpretación literal de «sin librerías» afecta otra dependencia de Python/R, no la declares conforme: documenta el caso y pide aclaración antes de asumir.
2. **EDA** puede usar las librerías autorizadas, pero las cantidades exigidas deben calcularse paso a paso y contrastarse con la función de librería requerida. La comparación contra librería debe ser reproducible y fallar si el error excede tolerancia.
3. No uses validación, prueba ni estaciones reservadas para decidir limpieza, percentiles, escalas, K, FPCA, arquitectura o hiperparámetros. Toda transformación aprendida se ajusta en entrenamiento y se congela antes de evaluar.
4. No conviertas curvas no emparejadas, medias por región o diferencias por ONI en “sesgo satelital” ni causalidad. Mientras falten pares válidos, márcalo como descriptivo y deja criterio 3 bloqueado.
5. Conserva NA como faltantes, salvo regla aprobada y trazada. No imputes ceros ni elimines outliers sin motivo, claves y conteos antes/después.
6. Mantén tolerancias matemáticas estrictas. No cambies una prueba a advertencia ni aflojes tolerancias para hacerla pasar.

### 3. Actualiza primero esta matriz viva

Crea `docs/SEGUIMIENTO_RUBRICA.md` y conserva ahí la tabla siguiente. **Actualízala después de cada lote**, incluyendo evidencia/rutas y pruebas. Estados permitidos: `PENDIENTE`, `EN PROGRESO`, `VERIFICADO`, `BLOQUEADO`. Solo escribe `VERIFICADO` si se ejecutó el cálculo o test descrito y dejó evidencia reproducible. No uses casillas marcadas solo porque exista código.

| Criterio / peso | Meta concreta para aspirar a 5 | Hecho comprobado al iniciar | Falta / no hecho | Librerías: requisito y evidencia | Prueba/evidencia que cerrará el criterio | Estado inicial |
|---|---|---|---|---|---|---|
| 1 Construcción / 10 % | Φ y MCO paso a paso; ECM/VC/VCG con gráficos y K por fuente; S idempotente/traza K; covarianza residual y bandas puntuales; dominio común. | Cox–de Boor/MCO manuales; K reportado 8/11 en Jacobo; W y bandas preparadas. | Conciliar datos de origen; mover comparaciones de paquete fuera del constructor; probar igualdad R/Python de ambas fuentes y W con aserciones. | Sin paquetes funcionales en cómputo; fda/splines/SciPy solo en pruebas separadas, nunca como fuente del resultado. | Pruebas manuales contra fórmulas/QR independientes; K, coeficientes de ambas fuentes, H, W y conteos con `stopifnot`/tests. | EN PROGRESO |
| 2 EDA / 15 % | Media paso a paso + biblioteca; varianza/DE/CV; covarianza/correlación; derivada e inicio/fin interpretados; FM/MBD/boxplot/atípicas con claves; contexto climatológico. | Notebook contrasta media con scikit-fda; incluye momentos, depths, boxplot y ENSO. | Alinear escala FM, empates de centralidad/recorte y límites del boxplot; validar o rebajar fechas heurísticas; evitar inferir sesgo ENSO sin pares. | Bibliotecas permitidas en EDA; contraste obligatorio y reproducible para ambas fuentes. | Tests numéricos y de empates/bordes, claves de outliers, informe con interpretaciones limitadas por evidencia. | EN PROGRESO |
| 3 Comparación / 15 % | Pares estación–píxel, regla y coordenadas; diferencias/bandas, covarianza/correlación, geografía, ENSO y sesgo. | Jacobo documenta que faltan coordenadas oficiales. La tabla CVC local solo coincide con 170/248 IDs Jacobo. | Obtener crosswalk oficial, resolver IDs faltantes y validar píxel/CRS; ejecutar comparaciones pareadas. Sin pares: la rúbrica limita a 1/5. | No usar paquetes para inventar el emparejamiento; todo cruce documentado. | Tabla de pares con origen/ID/distancia/píxel, auditoría y métricas pareadas reproducibles. | BLOQUEADO |
| 4.1 Depuración / 6 % | Imposibles, digitación, totales mensuales, ceros, huecos/años incompletos; criterios y trazabilidad sin fuga. | Notebook define reglas P1–P6, cobertura, banderas y selección desde entrenamiento. | Regenerar desde crudos; revisar umbrales/sensibilidad; conciliar claves baseline/Jacobo y revisar cómo se marcan ceros sospechosos. | Las reglas deben estar justificadas por clase/fuente; no introducir detectores opacos. | Manifiesto de filas/estación-año por regla, hashes de entradas y revisión de ejemplos. | EN PROGRESO |
| 4.2 Red / 10 % | Gram, FPCA solo entrenamiento, pérdida L²/rugosidad, equivalencia lineal ≤1e-4, no negatividad y recortes. | W y coeficientes transformados son preparativos. | Implementar y probar encoder, FPCA y cabezas lineal/KAN/MLP; reportar % recorte. | PyTorch/pykan solo en el apartado autorizado de modelos. | Tests de integración, leakage, equivalencia lineal y reconstrucción. | PENDIENTE |
| 4.3 Ablación / 11 % | KAN/MLP/lineal comparables, parámetros ±10 %, tres semillas, presupuesto igual e IC bootstrap. | No hay ablación implementada. | Implementar experimento y prueba sintética KAN RMSE≤1e-2; no seleccionar con prueba. | Librerías solo las autorizadas en modelos; registrar versiones. | Tabla de parámetros, tres semillas, configuración y bootstrap reproducible. | PENDIENTE |
| 4.4 Evaluación / 8 % | Estaciones reservadas; baseline satelital; L²/KGE/sesgo/picos; bootstrap≥1000 y Wilcoxon pareada; CSV de prueba. | Partición escrita y K de Jacobo se limita a entrenamiento. | Conseguir estaciones reservadas; entrenar y evaluar modelo; producir métricas/intervalos/pares. | Herramientas permitidas en modelos, nunca para fuga de prueba. | CSV station_id/year/day/value_mm y reporte reproducible de holdout. | BLOQUEADO |
| 4.5 Interpretación / 10 % | Ventanas/pesos KAN, activaciones, cabeza lineal, una activación simbólica con R² y casos extremos. | Hay interpretación climática descriptiva, no del modelo. | Implementar interpretabilidad y resultados de modelos entrenados. | Ajustes de interpretación deben seguir la metodología permitida. | Figuras/tablas de activaciones y ajuste simbólico con R². | PENDIENTE |
| 5 Informe / 10 % | PDF final ≤20 páginas, una columna, autores, referencias, figuras/tablas interpretadas e informe separado de IA/prompts. | Hay README y notebooks; no equivalen a la entrega final. | Producir informe final y anexo IA usando resultados verificados. | Citar todas las fuentes y declarar herramienta/prompts reales. | PDF validado por número de páginas/formato + anexo Word. | PENDIENTE |
| 6 Scripts / 5 % | Un R y un Python finales, reproducibles, con rutas/semillas/versiones; separados del informe. | Hay dos scripts R de fases y notebooks Python. | Consolidar los entregables finales y un comando de ejecución limpio. | Sin instalación oculta al analizar; dependencias y versiones documentadas. | Ejecución completa desde ambiente limpio y lista de artefactos. | PENDIENTE |

Incluye una columna adicional o sección **“Historial de cambios de estado”** con fecha, criterio, estado anterior/nuevo, evidencia y comando ejecutado. Si un elemento está bloqueado por datos del docente, continúa las tareas independientes; no inventes ni marques completo lo bloqueado.

### 4. Ejecuta los lotes en este orden

**Lote A — Línea base y datos.** Revisa README, rubric snapshot de este prompt, notebooks, scripts y resultados dentro del checkout local. Regenera el cache ignorado `eda_fase0/estaciones_crudas.pkl` desde los archivos fuente locales y reproduce el notebook 01 si existen sus entradas y dependencias. Compara las claves estación–año de main (4.426/238) y Jacobo (4.485/248): documenta las 66 claves solo-main y 125 solo-Jacobo con regla P1–P6 por cada una. Revisa dependencias, unidades, TIFF, máscaras y banderas. No subas datos fuente/cache.

**Lote B — Construcción manual.** Corrige diferencias de entrada y lógica; garantiza que el cálculo productor no dependa de bibliotecas funcionales; separa verificadores. Reproduce Python y R en la misma población y bloquea con aserciones si no coinciden el periodo, K, filas, claves, coeficientes por ambas fuentes, matriz W y métricas matemáticas.

**Lote C — EDA.** Alinea FM/MBD, ties, recorte y boxplot con clase; crea casos de prueba sintéticos de empates/límites. Valida las fechas o decláralas heurísticas. Quita lenguaje de sesgo/causalidad que no esté respaldado por pares.

**Lote D — Geografía y comparación.** Busca el crosswalk solo en los archivos locales y revisa cobertura/matches de CVC/Catálogo disponibles. Si no hay station_id→lat/lon/altitud fiable en el checkout, deja criterio 3 BLOQUEADO, continúa otros criterios y anota el dato exacto que falta pedir al docente. No declares parejas basadas solo en prefijo ni coordenadas inventadas.

**Lote E — Modelos, ablation, holdout e interpretación.** Implementa 4.2–4.5 conforme a rúbrica y con pruebas contra fuga de información. No reutilices coeficientes, escalas, K, FPCA ni selección de hiperparámetros aprendidos con validación/prueba.

**Lote F — Informe y entregables.** Produce el informe global y anexo IA, consolida un R y un Python, repite la ejecución limpia y actualiza la matriz final. Deja cualquier criterio no probado como PENDIENTE/BLOQUEADO y explica la causa.

### 5. Contrato de salida por lote

Al cerrar cada lote:
- Actualiza `docs/SEGUIMIENTO_RUBRICA.md` con estado y evidencia.
- Entrega lista de archivos modificados, comandos ejecutados y sus resultados reales.
- Reporta criterios `VERIFICADO`, `PENDIENTE` y `BLOQUEADO` por separado.
- Explica discrepancias que cambien universo, K o conclusiones.
- No asignes/garantices nota 5; indica qué falta para poder aspirar a ella.
- Continúa automáticamente entre lotes; no esperes aprobación intermedia. Detente solo ante una decisión irreversible, riesgo de pérdida de datos o una ambigüedad normativa que no pueda resolverse con la fe de erratas. Al final, presenta la matriz y pregunta únicamente por los bloqueos que requieren información de la usuaria/docente.

**FIN DEL PROMPT**

# Contexto para continuar con el EDA (para Claude Code)

> Léelo completo antes de tocar código. Resume el trabajo de estaciones hecho el 12-nov-2026
> en la rama `gus-estaciones-valle`, que sale de `origin/maria-eda-funcional-construccion-funcional`.
> Las decisiones, con su justificación y su validación, están en `DECISIONES.md` (D-01 a D-35).
> Si algo de aquí contradice a `DECISIONES.md`, manda `DECISIONES.md`.

## 1. El proyecto en una línea

Proyecto final de Analítica de Datos (UAO, 2026-2S): modelar la lluvia del **Valle del Cauca** como
datos funcionales con dos fuentes, estaciones diarias (`P_1930_1979.xlsx`, `P_1980_1990.csv`,
`P_1991_2003.csv`) y satélite pentadal (`ImagenesSatelitePentadal/`). Se pasa por construcción
B-spline, EDA funcional y regresión función a función / red (KAN/MLP).

## 2. Normas que mandan (no negociables)

- **Fe de erratas** (`Proyecto_Final_Analitica_FeErratas.docx`): periodo **1981–2003**.
  Partición: entrenamiento 1981–1995, validación 1996–2000, prueba 2001–2003, más las estaciones
  que reserva el profesor. Sin fuga de información (D-28).
- **Solo métodos vistos en clase.** La construcción funcional (base Φ, MCO, selección de K con ECM/VC/VCG,
  verificación de S) se programa **sin librerías**. El EDA **sí** permite librerías.
- **Entrega:** un único archivo R y un único archivo Python, separados del informe. El flujo vigente de la rama
  es `proyecto_funcional.R` (ver `README.md` y `REQUISITOS_Y_ENTREGA_R.md`).
- **Anexo de IA:** hay que entregar en un Word aparte los prompts usados y el nombre de la herramienta.
  Guarda los prompts de esta etapa.
- Commits **sin** autoría de Claude/IA. Gus hace el push desde su terminal.

## 3. Estaciones: qué está listo

Archivo principal: **`datos_web/estaciones_valle_finales.csv`** (y `.xlsx`), con **232 estaciones** del Valle.

| Columna | Significado |
|---|---|
| `codigo` | Código CVC de 10 dígitos (es el `station_id` de los archivos de lluvia) |
| `nombre`, `municipio`, `lat`, `lon` | Ubicación, en WGS84 |
| `alt` | Altitud oficial en m; vacía en 4 estaciones (no se estima) |
| `nivel_c` | `Ubicación oficial` (226) o `Cruce con CNE` (6) |
| `fuente` | De dónde salió la coordenada |
| `metodo_alt`, `fuente_alt` | `Excel del curso` (163), `CNE IDEAM` (65) o `Sin altitud oficial` (4) |
| `srtm` | Terreno SRTM en la coordenada: **solo control de calidad**, no es altitud |
| `revisar` | Texto si la coordenada es dudosa (23 estaciones; altitud vs terreno > 300 m) |
| `años`, `n_años` | Años con datos, **sin filtrar** |

- Inventario completo de los 344 códigos (Valle, Cauca, sin ubicar, con niveles de evidencia):
  `datos_web/estaciones_valle_del_cauca.xlsx`.
- **Fuera del análisis:** 35 del Valle con departamento inferido y sin coordenadas (descartadas, D-34),
  56 del Cauca (fuera del área) y 21 sin ubicar (pendientes, P-1).
- Las 6 del cruce con el CNE (Palmira ICA, Aeropuerto Bonilla Aragón, Central Castilla, Manuel M. Mallarino,
  Cent. Admo. La Unión, El Queremal) se incluyen, pero el **análisis de sensibilidad** debe correrse con y sin ellas.
- Coordenadas "a revisar": no se corrigen sin una segunda fuente. Caso más claro: Cajones `2612200105`
  (el punto del CNE es coherente; el del CSV no).
- Estructura del código CVC (D-32): dígitos 1-2 zona (26 Cauca, 53 Pacífico, 54 San Juan), 3 margen del río Cauca
  (1 derecha/este, 2 izquierda/oeste), 4 tramo de sur a norte, 5-6 cuenca, 7-8 tipo (01 pluviométrica/pluviográfica,
  02 climatológica, 03 evaporimétrica), 9-10 consecutivo.

## 4. Qué sigue en el EDA (por hacer)

1. **Conectar las 232 estaciones al flujo** de la rama (`proyecto_funcional.R` / notebooks), en lugar de las listas anteriores (P-6, P-7).
2. **Filtro por años, justificado en el EDA:** cobertura por año-estación dentro de 1981–2003, criterio de año válido
   (D-11 y su reemplazo en D-29), cuántas estaciones sobreviven y por qué.
3. **Calidad de los datos:** faltantes por estación y año, valores imposibles o marcados como erróneos (D-12),
   duplicados. No se imputan estaciones sin justificarlo (D-13).
4. **Satélite:** recortar al polígono del Valle (`datos_web/limites/`, D-03), unidades mm/día (D-15),
   años defectuosos 1985, 1988 y 2000 (D-16; 2000 cae en validación: pregunta P-4).
5. **Emparejar cada estación con su píxel** usando `lat`/`lon` (P-7).
6. **EDA funcional por fuente:** media, varianza, covarianza y correlación funcionales; derivadas; profundidad FM y MBD;
   boxplot funcional con atípicas por estación y año; interpretación climatológica (régimen bimodal, orografía/altitud, ENSO con ONI).
7. **Representatividad:** si se usa un subconjunto de estaciones o años, demostrar que es representativo de 1981–2003.
8. **Sensibilidad:** con y sin las 6 del cruce con el CNE, y con y sin las 23 coordenadas a revisar.

## 5. Preguntas abiertas para el profesor (ver `DECISIONES.md`, sección 8)

- P-1: catálogo o coordenadas de los 21 códigos sin ubicar.
- P-3: ¿se puede reemplazar el satélite entregado por CHIRPS?
- P-4: ¿cómo tratar el año 2000 (satélite defectuoso) que cae en validación?

## 6. Código de datos (al 12-nov-2026)

Todo el entendimiento de los datos, la depuración y el dataset final están en **`lluvia_valle_datos.R`**,
en la raíz. Se corre de arriba abajo y deja `salidas_R_datos/lluvia_valle_limpia.rds`, más las tablas y los gráficos.
Las reglas están en D-36.

## 7. Archivos sin commit en la rama (histórico)

`DECISIONES.md`, `CONTEXTO_EDA.md`, `datos_web/estaciones_valle_del_cauca.xlsx`,
`datos_web/estaciones_valle_finales.{csv,xlsx}`, `datos_web/cne_valle_cauca_completo.csv`,
`datos_web/limites/`, `analisis_R/`. Borrar `datos_web/estaciones_valle_sin_coordenadas.xlsx` (reemplazado).

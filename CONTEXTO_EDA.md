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
  Partición: entrenamiento 1981–1995, validación 1996–2000, prueba 2001–2003
  (el profesor ya no reserva estaciones). Sin fuga de información (D-28).
- **Solo métodos vistos en clase.** La construcción funcional (base Φ, MCO, selección de K con ECM/VC/VCG,
  verificación de S) se programa **sin librerías**. El EDA **sí** permite librerías.
- **Entrega:** un único archivo R y un único archivo Python, separados del informe. El código de datos vigente
  es `lluvia_valle_datos.R`; lo funcional se le irá agregando. El flujo anterior (`proyecto_funcional.R`, notebooks,
  `REQUISITOS_Y_ENTREGA_R.md`) está en `archivo_flujo_anterior/`, solo como referencia.
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
  56 del Cauca (fuera del área) y 21 sin ubicar (descartados, P-1).
- Las 6 del cruce con el CNE (Palmira ICA, Aeropuerto Bonilla Aragón, Central Castilla, Manuel M. Mallarino,
  Cent. Admo. La Unión, El Queremal) se incluyen, pero el **análisis de sensibilidad** debe correrse con y sin ellas.
- Coordenadas "a revisar": no se corrigen sin una segunda fuente. Caso más claro: Cajones `2612200105`
  (el punto del CNE es coherente; el del CSV no).
- Estructura del código CVC (D-32): dígitos 1-2 zona (26 Cauca, 53 Pacífico, 54 San Juan), 3 margen del río Cauca
  (1 derecha/este, 2 izquierda/oeste), 4 tramo de sur a norte, 5-6 cuenca, 7-8 tipo (01 pluviométrica/pluviográfica,
  02 climatológica, 03 evaporimétrica), 9-10 consecutivo.

## 3b. Datos listos (13-nov-2026)

`lluvia_valle_datos.R` deja `salidas_R_datos/datos_preparados.rds`: matriz de estaciones (3 326 × 365), matriz del satélite
(15 663 × 72, mm/día) y los pares estación–píxel, con 26 pruebas automáticas (D-36, D-43 a D-45). El satélite es CHIRP solo-infrarrojo (D-43)
y no capta los años Niño/Niña (D-44). El control de lo que pide el enunciado está en `CONTROL_RUBRICA.md`.

## 4. Qué sigue

Ya están hechos: el filtro por años (D-36), la calidad de los datos (secciones 2–5 de `lluvia_valle_datos.R`) y la
representatividad de las 176 estaciones incluidas (sección 8). Sigue:

1. ~~Satélite: recorte, mm/día y EDA~~ — hecho (D-43, D-44).
2. ~~Emparejar cada estación con su píxel~~ — hecho (D-45).
3. ~~Construcción funcional, EDA funcional, comparación y regresión de referencia~~ — hechos por Jacobo y María en `funcional/10` a `14` (D-37 a D-42), ya con los datos nuevos.
4. **Acordar P-4:** si 1985, 1988 y 2000 del satélite entran al modelo (D-37/D-42 los sacan; D-44 dice que no están dañados).
5. **Sensibilidad:** con y sin las 6 del cruce con el CNE, y con y sin las 12 coordenadas a revisar que quedaron entre las 176.

## 5. Preguntas abiertas para el profesor (ver `DECISIONES.md`, sección 8)

- ✅ P-1: los 21 códigos sin ubicar se descartan.
- ✅ P-3: se usa solo el satélite entregado (CHIRP, confirmado), sin rellenar con otra fuente.
- P-4: ¿cómo tratar el año 2000 (satélite defectuoso) que cae en validación?

## 6. Organización del repositorio

Ver `README.md`:
- raíz: `lluvia_valle_datos.R`, `salidas_R_datos/`, `DECISIONES.md` y este archivo;
- `datos_curso/`: datos crudos, sin subir a git;
- `datos_web/`, `enunciado/`, `material_clase/`;
- `archivo_flujo_anterior/`: lo anterior.

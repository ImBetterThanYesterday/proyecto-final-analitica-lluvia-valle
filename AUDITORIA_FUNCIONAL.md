# Auditoría de la parte funcional (`funcional/10` a `14`)

Rama `gus-auditoria-funcional`, 13-nov-2026. Se revisó el código de Jacobo y María contra el enunciado (fe de erratas), la rúbrica y los datos cerrados (D-43 a D-46). Los puntos dudosos se probaron corriendo el código; los experimentos están al final.

## Estado después del cierre (13-nov-2026, D-48)

Todo lo funcional quedó en un solo archivo, `lluvia_valle_funcional.R`, que reemplaza a `funcional/10` a `14`. Qué pasó con cada hallazgo:

| Hallazgo | Estado | Dónde se ve |
|---|---|---|
| G1. Años 1985, 1988 y 2000 | ✅ Entran todos | Secciones 1 y 5.6; `F1_satelite_vs_estaciones_por_anio.csv`, `F5_sensibilidad_anios_marcados.csv` |
| G2. Métricas de la rúbrica | ✅ Agregadas | Sección 5.5; `F5_metricas_prueba.csv`, `F5_regresion_vs_satelite_pareado.csv`, `F5_metricas_por_fase_ENSO.csv` |
| G3. Código duplicado y seis archivos | ✅ Un archivo que lee `datos_preparados.rds` | Sección 1 |
| I1. EDA que mezcla dos climas | ✅ Por vertiente y con anomalías | Sección 3.2; `F3_covarianza_por_vertiente_y_anomalias.csv` |
| I2. K distinto por fuente | ✅ Verificado: no cambia las conclusiones | Sección 4.2; `F4_sensibilidad_mismo_K.csv` |
| I3. Orden de la raíz y el suavizado | ✅ Verificado: suavizar en mm conserva el acumulado anual | Sección 5.1; `F5_orden_raiz_y_suavizado.csv` |
| I4. K = 8 rígido | ✅ Sensibilidad hecha; se mantiene K = 8 | Sección 2; `F2_sensibilidad_K_estaciones.csv` |
| I5. Bordes | ✅ Documentado | Bloque "Qué dicen los datos" de la sección 2 |
| I6. Dos definiciones de vertiente | ✅ Una sola: estación más cercana | Sección 1; `F1_regla_vertiente_pixeles.csv` |
| Ventanas fijas de temporada | 🟡 Se mantienen; se reporta el % de picos en el borde de la ventana | `F2_sensibilidad_K_estaciones.csv` |
| Píxeles compartidos | 🟡 Documentado, sin corrección | `F1_pixeles_compartidos.csv` |
| Exportación a Python | ✅ `pares.csv` ya no lleva `defectuoso` | Sección 5.7 |

Lo que sigue describe la auditoría tal como se hizo, sobre el código anterior.

## Veredicto

**La base es sólida, pero hay que corregir varias cosas antes de pasar a la red.**

- La parte técnica de la construcción está bien hecha y verificada: base propia, mínimos cuadrados, S, W, FPCA y la equivalencia β(s,t).
- Los problemas están en las **decisiones de modelado y de evaluación**:
  - se quitan años sin una justificación sólida;
  - las métricas no son las que pide la rúbrica;
  - el EDA mezcla dos climas distintos;
  - el código del satélite está duplicado.

## Lo que está bien (se conserva)

| Qué | Evidencia |
|---|---|
| Base B-spline propia (Cox–de Boor) | Coincide con `splines::splineDesign` (3e-16) |
| Atajo de la VC con diag(S) | Da lo mismo que dejar cada día fuera (365 ajustes) |
| S simétrica, idempotente y tr(S) = K | Error < 1e-14, también con días faltantes |
| W por cuadratura de Gauss | 1'W1 = 364 exacto; c'Wc = trapecio fino |
| FPCA solo con entrenamiento | Coincide con `fda::pca.fd` |
| Equivalencia integral de β vs puntajes | 5e-15 (el criterio es ≤ 1e-4) |
| Sin fuga | K, μ, ψ, la escala de la altitud, M y λ salen de entrenamiento o validación; la prueba no se toca |
| EE por estación en d(t) y bootstrap por estación | Las curvas de una estación no se tratan como independientes |
| Modelo de referencia "solo altitud y ONI" | Muestra que el satélite sí aporta |
| Días faltantes | Cada curva se ajusta con sus días observados; no se imputa |

## Hallazgos, de más grave a menos grave

### Graves: afectan la nota o la validez

**G1. Quitar 1985, 1988 y 2000 no se sostiene.** El detalle está en la sección "Años del satélite".

**G2. Las métricas no son las de la rúbrica.** La rúbrica (4.4) pide:
- **error L² relativo** — el código reporta el RMSE absoluto;
- **sesgo del acumulado anual en %** — el código lo reporta en mm;
- **prueba de Wilcoxon pareada por estación** — no está.

Probado: con la regresión actual, el L² relativo baja de 0,52 a 0,45, el sesgo pasa de +26 % a +8 % y el Wilcoxon da p = 0,009. Hay que agregarlas, porque la red se va a evaluar con las mismas funciones.

**G3. Dos procesamientos del satélite y seis archivos R.**
- La sección 10 rehace lo que ya hacen las secciones 9 a 12 de `lluvia_valle_datos.R`: lectura, mm/día, recorte y emparejamiento.
- La sección 11 lee `lluvia_valle_limpia.rds` y no `datos_preparados.rds`, así que las 26 pruebas de los datos no cubren lo que entra a lo funcional.
- La entrega pide **un único archivo R**.
- Qué hacer: que lo funcional parta de `datos_preparados.rds`, borrar la sección 10 duplicada (dejar solo lo que no exista) y unir todo en un solo archivo.

### Importantes: metodológicas

**I1. El EDA funcional mezcla el Pacífico con el valle interandino.**
- La media, la varianza y la covarianza se calculan con todas las estaciones juntas. Unas tienen más de 7 000 mm/año y otras 1 000, y el Pacífico no tiene el régimen bimodal.
- El primer componente del FPCA explica el 61 %, y es "el nivel de lluvia de cada sitio": la covarianza la dominan las diferencias entre lugares, no la forma del año.
- Qué hacer: reportar el EDA **por vertiente** y además sobre **anomalías** (restando la curva media de cada estación). Las 6 fechas de temporada sí se calculan por vertiente, y eso está bien.

**I2. La comparación de fuentes usa K distintos: estación K = 8 y satélite K = 11.**
- La curva de la estación queda más suave que la del satélite.
- Así, la comparación de varianzas, covarianzas y correlación cruzada de la sección 13 en parte compara grados de suavizado, no fuentes.
- Qué hacer: en la sección 13, comparar con el mismo K para las dos fuentes, o al menos mostrar que las conclusiones no cambian.

**I3. Primero se suaviza en mm y después se saca la raíz.**
- El enunciado pide modelar sobre √P. El código ajusta la curva en mm, la recorta en 0, le saca raíz y la vuelve a proyectar (error de 0,09).
- La lluvia diaria tiene colas muy pesadas: el ajuste en mm y la elección de K quedan dominados por los aguaceros.
- La alternativa más limpia es sacar la raíz a los datos diarios y luego suavizar, eligiendo K en esa escala. Hay que probarla antes de decidir.

**I4. K = 8 para las estaciones es muy rígido para medir fechas de picos.**
- Con 8 bases hay un nodo cada ~73 días, y la VC es plana entre K = 5 y 17: el criterio no distingue.
- La métrica "error en la fecha de los dos picos" depende directamente de esa rigidez.
- Qué hacer: análisis de sensibilidad con K = 12 o 16 sobre las métricas de picos.

**I5. Bordes de la base no periódica.**
- El 47 % de las curvas de estación baja de 0 en algún día, sobre todo el 1 de enero, el 31 de diciembre y en julio–agosto.
- El satélite no tiene datos antes del día 3 ni después del 363, así que en los bordes la curva extrapola.
- Qué hacer: documentarlo. Una base periódica (Fourier) sería más natural para el ciclo anual, pero en clase se usan B-splines, así que basta con justificarlo y recortar en 0.

**I6. La vertiente se define distinto en cada fuente.** Para las estaciones sale del código CVC; para los píxeles, de un corte en la longitud −76,8°. Hay que unificar el criterio.

### Menores

- **Ventanas fijas de temporada** (días 60–200 y 230–365): asumen el régimen bimodal, que no se cumple en el Pacífico.
- **Píxeles compartidos:** 23 píxeles tienen 2 a 4 estaciones con la misma curva de satélite. El bootstrap por estación no lo corrige del todo; basta con mencionarlo.
- **Exportación a Python:** `para_python/pares.csv` lleva las columnas `defectuoso` y `usar`, así que la red heredaría la exclusión de años. Se corrige al resolver G1.

## Años del satélite (1985, 1988 y 2000): ¿se pueden quitar?

**Respuesta corta: no hay razón válida para quitarlos. Quitarlos empeora el modelo y debilita la evaluación.**

> **Actualización (asesoría del 13-nov-2026):** el profesor confirmó que **no se deben eliminar años**, porque el modelo tiene que aprender esos casos (D-47). Queda decidido.

### 1. El criterio para escogerlos no es coherente
Se eligieron por ser los tres años más secos del satélite. Pero las mismas notas de la sección 10 dicen que:
- 1985 no se aparta de las estaciones (cociente 1,03), porque ese año las estaciones también fueron secas;
- 1984 se aparta tanto como 1988 y no se marcó.

Así que la marca no corresponde a "el satélite falló".

### 2. No son archivos dañados: son años La Niña
- D-44 y las propias notas de Jacobo lo confirman.
- El satélite (CHIRP, solo infrarrojo) subestima la lluvia en La Niña. Ese es un error **sistemático**, y es justo lo que el modelo tiene que aprender a corregir; para eso está el ONI en z.
- Quitarlos le esconde al modelo los casos que más necesita ver: de los 4 años Niña de entrenamiento según el ONI anual (1984, 1985, 1988 y 1989) se pierden 2.

### 3. Quitar 2000 de validación arregla la evaluación a favor del modelo
- Si el modelo se elige sin el año en que el satélite falla, la validación se ve mejor de lo que es.
- Además, una regla así no se puede aplicar en prueba: no se puede decidir de antemano que un año "no cuenta".
- El profesor dijo que hay que trabajar con el satélite tal como viene y resolver los problemas que aparezcan.

### 4. Los datos lo confirman
Se volvió a correr la regresión de la sección 14 incluyendo esos años; todo lo demás quedó igual. Los dos modelos se evalúan en los mismos años de prueba (2001–2003).

| Métrica en prueba | Sin esos años (hoy) | Con esos años |
|---|---|---|
| RMSE L² (mm/día) | 2,24 | **2,21** |
| L² relativo | 0,45 | **0,44** |
| KGE | 0,45 | **0,46** |
| Sesgo anual | +8,5 % | **+7,5 %** |
| Curvas en que mejora al satélite | 60 % | **63 %** |
| Wilcoxon por estación vs satélite | p = 0,009 | **p = 0,0002** |
| Sesgo en 2003 (neutro) | +134 mm | **+45 mm** |

- Incluirlos mejora la prueba en los tres años (2001, 2002 y 2003) y en todas las métricas.
- La elección de M no cambia (M = 5) y λ sube un poco, de 1e8 a 3e8.

### Propuesta
1. Incluir los tres años en todo: elección de K, comparación, modelo y red.
2. Mantener la columna `defectuoso` solo para un **análisis de sensibilidad**: "sin esos años, el resultado es…".
3. Contarlo en el informe como hallazgo: el satélite falla en La Niña, y el modelo con el ONI lo corrige en parte.
4. Revisar D-37 y D-42 con Jacobo y cerrar P-4.

## Orden de trabajo sugerido

| # | Qué | Por qué primero |
|---|---|---|
| 1 | G1: incluir los años; sensibilidad sin ellos | Cambia todos los resultados de abajo |
| 2 | G2: L² relativo, sesgo en %, Wilcoxon | La red se evalúa igual; es rúbrica directa |
| 3 | G3: partir de `datos_preparados.rds` y unir en un solo R | Evita resultados que no cuadran entre los dos códigos |
| 4 | I3 e I4: probar la raíz antes de suavizar y la sensibilidad de K | Define la curva objetivo de la red |
| 5 | I1, I2 e I6: EDA por vertiente y anomalías, mismo K en la comparación, una sola definición de vertiente | Mejora la interpretación (15 % + 15 %) |
| 6 | I5 y menores: documentar | Informe |

## Cómo se hicieron las pruebas

- **G1 y G2:** se corrió una copia de `funcional/14_regresion_funcional.R` con una sola diferencia: incluir o no los pares de 1985, 1988 y 2000. Se agregaron el L² relativo, el sesgo en % y el Wilcoxon pareado (media por estación del L² relativo). Las salidas quedaron fuera del repositorio; `salidas_R_funcional/` no se tocó.
- **Lo demás** sale de leer el código y de sus propias salidas (`11_seleccion_K.csv`, `14_fpca_varianza.csv`, `12_temporadas_fechas.csv`).

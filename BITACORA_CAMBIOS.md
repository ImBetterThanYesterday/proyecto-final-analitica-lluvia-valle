# Bitácora de cambios

Qué cambió en cada sesión, en qué commit y por qué. El detalle de cada decisión está en `DECISIONES.md` (D-xx) y el estado frente a la rúbrica en `CONTROL_RUBRICA.md`.

## Sesión del 21-nov-2026 (rama `gus-base-conocimiento-profesor`)

| Archivo | Cambio | Por qué |
|---|---|---|
| `lluvia_valle_datos.R` | R7: Moran, revisión de pepita pura, validación de 10 alternativas (años impares para ajustar, pares para evaluar; regla del 1 % del proyecto 1), prueba por fase ENSO y prueba de un semivariograma por día | Validar cada decisión del kriging con el flujo de la clase (D-64) |
| `lluvia_valle_datos.R` | R7 no usa como vecinas las 7 estaciones sin ningún año válido | La Perdición (Calima) anota la lluvia de la semana en un solo día |
| Documentos | FUNDAMENTOS_REGLAS (sección 6: literatura y validación), DECISIONES (D-64), EXPLICACION_DATOS (8b) | Tesis del profesor, Russi y Larrahondo 2017, Lebel et al. 1987, Terán-Chaves 2021, Hofstra 2008, climatol |

Resultados en prueba iguales: satélite 0,51 → regresión 0,44 (principal).

## Sesión del 20-nov-2026 (rama `gus-base-conocimiento-profesor`)

Punto de partida: commit `f7ef057`. Detalle en D-61 y en `docs/FUNDAMENTOS_REGLAS.md` (secciones 2, 3 y 6).

| Archivo | Cambio | Por qué |
|---|---|---|
| `lluvia_valle_datos.R` | **R4:** decide con las 3 vecinas más cercanas (≥2 deben contradecir la racha) y q de la temporada; cita corregida a Lewis et al. 2021 | La versión con una sola vecina era frágil (30 rachas; ahora 10). Durre (2010) no revisa rachas de ceros |
| `lluvia_valle_datos.R` | **R7:** kriging ordinario de los faltantes de años válidos, con las vecinas escaladas por la razón de medias mensuales; semivariograma de entrenamiento elegido por validación cruzada (exponencial); prueba contra B-spline sola e IDW | OMM 2017, 4.6: para valores diarios sirve la interpolación espacial, no los días vecinos. En los datos, rellenar baja el error de la curva de 0,55 a 0,30 |
| Los cuatro scripts | Variable `CONJUNTO` (`principal` o `r6_omm`): R6 con ≥10 años (principal, D-63) o el 80 % de la OMM (sensibilidad), con salidas en carpetas `_r6_omm` | El profesor prefirió ≥10 años; la sensibilidad mide si la conclusión depende de R6 |
| `lluvia_valle_funcional.R` | "Curva completa" = 365 días observados; la banda del ejemplo de huecos se calcula sin los días rellenados | Con R7 ya no quedan faltantes en la matriz |
| `lluvia_valle.R` | Reconstruido con los cambios de los componentes | Sigue siendo la unión sin cambios de los dos `.R` |
| Documentos | FUNDAMENTOS_REGLAS, FALTANTES_E_IMPUTACION, ANTECEDENTES_VALLE, EXPLICACION_DATOS, EXPLICACION_FUNCIONAL (parcial), DECISIONES (D-13 → D-61, D-62, D-63), README, CONTROL_RUBRICA | Cifras nuevas: 163 estaciones y 3 046 curvas (principal); 100 y 2 113 (sensibilidad OMM) |
| `DECISIONES.md` | D-62: las estaciones se modelan en diario (el enunciado pide curvas diarias y de péntadas en el mismo dominio). D-63: R6 principal con ≥10 años | Revisión del profesor y fundamento en la OMM, NOAA y Urrea et al. 2019 |

Resultados en prueba (error L² relativo): principal (≥10 años), satélite 0,51 → regresión 0,44; sensibilidad (80 % OMM), 0,52 → 0,42. La conclusión no depende de R6.

## Sesión del 14-nov-2026 (rama `gus-cierre-funcional`)

Punto de partida: commit `fe34b7d` (parte funcional cerrada en un solo archivo, R4 calculada).

### 1. Correcciones rápidas — commit `46e3852`

| Archivo | Cambio | Por qué |
|---|---|---|
| `lluvia_valle_funcional.R` | "375 curvas" → **370**; título "1996-1999" → **1996-2000** | Eran errores de texto; los conteos salen de la tabla de pares |
| `lluvia_valle_funcional.R` | Línea del "mismo K": ahora dice que el diagnóstico global no cambia, pero el veredicto día a día coincide solo en **281 de 365 días** | Antes exageraba ("las diferencias no salen de usar K distinto") |
| `lluvia_valle_funcional.R` | **Prueba 9 reescrita:** prueba de no fuga real. Se reemplazan por ruido validación y prueba, se reajusta y la media, el FPCA, las escalas y B salen idénticos; con la prueba alterada se vuelven a elegir M y λ iguales | La prueba vieja solo revisaba las máscaras de partición (era tautológica) |
| Los dos `.R` | Fijan **UTF-8** solos al arrancar | En la Mac del grupo LANG está vacío y el script se detenía al comparar "validación" |
| `lluvia_valle_funcional.R` | Orden fijo de los grupos de altitud en las tablas F3 y F4 | El orden dependía del idioma de la máquina |
| `.md` | P-4 resuelta en README y CONTEXTO_EDA; CONTROL_RUBRICA y AUDITORIA sin referencias viejas; P-20 registrada | Estaban desactualizados |

### 2. Datos del curso localizados

`datos_curso/` nunca estuvo en git (está en `.gitignore`). Los archivos están en la rama de Federico, `origin/Arquitectura-KAN-Base`, carpeta `Datos/`. Se sacan sin cambiar de rama:

```bash
mkdir -p datos_curso
git archive origin/Arquitectura-KAN-Base Datos | tar -x --strip-components=1 -C datos_curso
find datos_curso -name '._*' -delete
```

⚠️ Esa rama sube los datos crudos al repositorio y quitó `datos_curso/` del `.gitignore`: **no fusionarla así**.
También se alineó el `main` local, que era un commit suelto: `git branch -f main origin/main`.

### 3. R0 y lectura del satélite — commits `23db96f` y `b1bd344`

| Cambio | Por qué | Resultado |
|---|---|---|
| **Regla R0** en `lluvia_valle_datos.R`: valor < 0 o > 1 825 mm → faltante (D-50) | La rúbrica 4.1 pide una regla explícita para valores imposibles; 1 825 mm es el récord mundial de lluvia en 24 h de la OMM | En el Valle no quita nada (máximo 600 mm); en el crudo completo detecta 859 375 mm de una estación de fuera. Fila "1b. R0" en la trazabilidad y prueba nueva: **27 pruebas** |
| El satélite se lee **por lotes de 200** archivos | Abrir los 2 328 GeoTIFF a la vez superaba el límite de archivos abiertos de macOS | El script corre completo en la Mac |
| **P-20 cerrada: K = 8** (D-51) | Con K = 16 la fecha del pico parecía mejorar, pero el K que gana cambia con la referencia (medias móviles de 15, 31 o 61 días) | Tablas `F2_K_picos_vs_referencia.csv` y `F2_K_que_gana_por_referencia.csv`; el script se detiene si algún día una sola K ganara con todas |

### 4. Red neuronal funcional y ablación — `lluvia_valle_red.py` (D-52)

Archivo nuevo en Python. Lee solo `salidas_R_funcional/para_python/`.

| Paso | Qué se hizo | Cómo se verificó |
|---|---|---|
| 1. Insumos | Lee pares, coeficientes en √P, W, R, FPCA y la especificación de R | W igual a la de R (8e-14); FPCA igual (1e-14); reproduce las métricas de prueba de la regresión de R |
| 2. Red | Codificador z = A W c (8 neuronas) + altitud + ONI; decodificador FPCA de R congelado (M = 6); cabezas lineal, KAN [10,4,6] y MLP [10,37,6] | — |
| 3. Lineal = regresión | Demostración escrita en el código | La red con los coeficientes de R da sus predicciones (3,6e-15); la identidad con β(s,t) a 1,5e-14 (criterio 1e-4) |
| 4. KAN y MLP | KAN probada con una función sintética; MLP con los mismos parámetros | RMSE 2e-3 (criterio 1e-2); 640 vs 635 parámetros (0,8 %) |
| 5. Ablación | Mismo presupuesto: 5 pruebas de λ_w por cabeza y 3 semillas | 14 de 14 pruebas pasan |

**Tres cambios frente al script viejo (`archivo_flujo_anterior/05_red_funcional.py`):**

1. **Rugosidad relativa** ∫w''²/∫w²: con z estandarizado, la penalización vieja se podía bajar encogiendo las w_k sin suavizarlas.
2. **Optimizador:** con lr = 3e-3 y paciencia 300, la cabeza lineal paraba antes de converger (pérdida de validación 142, peor que mínimos cuadrados). Ahora llega a 126, igual que la regresión de R. Sin eso, la ventaja de KAN y MLP habría sido falsa.
3. **Rejilla de λ_w** {0, 1e6, …, 1e9}: la vieja (≤ 0,1) no penalizaba nada.

**Resultado en prueba (2001–2003):** L² relativo satélite 0,513 → lineal 0,440, KAN 0,434, MLP 0,420. El MLP supera a la lineal por poco; la KAN no se distingue de la lineal en L² ni KGE, pero reduce el sesgo. La no linealidad se ve más en validación (que tiene el Niño de 1997 y las Niñas de 1999 y 2000) que en prueba.

### 5. Regresión contrastada con los métodos de clase — sección 5.5b del funcional (D-53)

El enunciado exige usar solo métodos vistos en clase. Nuestra regresión se comparó con las dos versiones de clase del modelo integral: la del taller de repaso (ejercicio 10, paso 3: puntajes de Y sobre puntajes de X) y `fda.usc::fregre.basis.fr()` (ejercicio 9). Con parámetros elegidos en validación, en prueba las tres dan el mismo L² relativo (0,43–0,44) y la nuestra tiene mejor KGE (+0,03, IC que excluye el 0). Se mantiene la nuestra. El script ahora usa `fda` y `fda.usc` y tiene **13 pruebas**. En la Mac hay que instalar `fda.usc` una vez: `Rscript -e 'install.packages("fda.usc", repos = "https://cloud.r-project.org")'`.

### 6. Python solo entrena; R evalúa — sección 7 del funcional (D-54)

El apartado de la red solo permite PyTorch y pykan. `lluvia_valle_red.py` ya no usa `scipy` ni `matplotlib`: entrena, verifica y exporta las predicciones (`salidas_python/predicciones_coef_raiz.csv`) y los pesos del codificador. La sección 7 nueva de `lluvia_valle_funcional.R` calcula las métricas, los IC, Wilcoxon y la tabla por fase ENSO con las mismas funciones de la regresión. Hallazgo: en años Niño (1997, 2002) KAN y MLP son mejores; en años Niña (1999-2000) la lineal es la mejor y las no lineales sobrestiman. Las tablas viejas de Python con métricas (`P1_regresion_R_reproducida`, `P5_ablacion_media_DE_semillas`, `P5_comparaciones_pareadas_prueba`, `P5_metricas_por_curva`) se reemplazan por las `F7_` de R.

### 7. Interpretación de la red (4.5) — sección 8 del funcional (D-55)

Python exporta además los pesos de la cabeza lineal, las activaciones de la primera capa de la KAN y los modelos (`salidas_python/modelos/`); los resultados no cambian (mismas predicciones, byte a byte). La sección 8 de R grafica las w_k, compara la β(s,t) de la red con la de la regresión, mide la importancia y la no linealidad de cada entrada, hace el ajuste simbólico y muestra los peores casos. Lo principal: KAN y MLP miran contrastes entre temporadas; la corrección satura (tanh) cuando el satélite marca mucha lluvia; la altitud es la entrada más no lineal; los peores casos son estaciones secas (fondo del valle del Cauca y cuenca del Dagua) que el píxel no representa.

### 8. Un solo archivo de R — `lluvia_valle.R`

Une `lluvia_valle_datos.R` (Parte A) y `lluvia_valle_funcional.R` (Parte B) sin cambiar el código: se verificó que, quitando lo agregado, quedan exactamente los dos archivos. Lo agregado: un encabezado con el mapa de la rúbrica, una marca "RÚBRICA x" antes de cada sección y, entre las partes, `rm(list = ls())` (vacía la memoria de R, no los archivos; la Parte B lee del disco lo que dejó la A). Probado aquí: Parte A hasta la sección 8 + Parte B completa en la misma sesión, 13 + 5 pruebas. Los dos archivos originales se mantienen hasta que el equipo decida cuál se edita (`RESUMEN_PARA_JACO.md`, decisión 4).

### 9. P-19 y P-16 cerradas, y corrida completa de `lluvia_valle.R` (D-56)

| Archivo | Cambio | Por qué |
|---|---|---|
| `lluvia_valle_funcional.R` y `lluvia_valle.R` | Nueva sección 5.2b: lee `P_1930_1979.xlsx`, aplica R5 a las 176 estaciones, ajusta las curvas y mide cuánto las explica el FPCA de entrenamiento y cuánto cambia si se agregan | P-19 necesitaba respaldo en el código. Resultado: 1 250 curvas válidas, explicadas al 89,8 %; agregarlas deja iguales las cuatro primeras funciones propias. No se usan |
| `DECISIONES.md` | D-56; P-16 y P-19 cerradas | ONI anual = media de los 12 trimestres; sin registros de 1930–1979 |
| `CONTROL_RUBRICA.md`, `RESUMEN_PARA_JACO.md` | 4.1 en ✅ | Ya no queda nada abierto en el preprocesamiento |

Corrida completa de `lluvia_valle.R` en Windows (R 4.4.2): pruebas 27 de 27, 13 de 13 y 5 de 5. Las salidas coinciden con las de la Mac salvo redondeo (diferencia máxima 1e-10), así que no se vuelven a subir; solo se agregan `F5_P19_conteo.csv` y `F5_P19_registros_1930_1979.csv`. Hace falta tener instalado `fda.usc`.

### Qué sigue

1. ~~4.5~~ — hecho (D-55).
2. ~~4.4 por fase ENSO~~ — hecho (D-54).
3. ~~P-16 y P-19~~: cerradas (D-56). Falta confirmar por escrito con el profesor que no hay CSV de estaciones reservadas.
4. Unir los dos R en uno para la entrega; informe PDF; anexo de IA con los prompts.

### Cómo correr todo

```bash
Rscript lluvia_valle.R              # todo R en uno (o los dos de abajo por separado)
Rscript lluvia_valle_datos.R        # ~8 min, 30 pruebas (necesita datos_curso/)
Rscript lluvia_valle_funcional.R    # ~1 min, 13 pruebas (necesita fda y fda.usc); deja para_python/
python lluvia_valle_red.py          # ~15 min, 13 pruebas (pip install torch numpy pandas)
Rscript lluvia_valle_funcional.R    # otra vez: secciones 7 (evaluación, 5 pruebas más) y 8 (interpretación)
```

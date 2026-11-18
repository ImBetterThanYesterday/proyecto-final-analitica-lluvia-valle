# Propuestas KAN de Federico (punto de partida)

Scripts seleccionados por sus resultados de validacion:

- `lluvia_valle_red_kan_rugosidad.py`: KAN + penalizacion D²; busca `lambda_D2` usando validacion y vuelve a entrenar la configuracion elegida con tres semillas.
- `lluvia_valle_red_kan_g3.py`: KAN `[10,4,6]` con rejilla fija `G=3`, la mejor rejilla en la comparacion exploratoria.
- `lluvia_valle_red_kan_l1.py` (penalización L1): se retiró de esta carpeta porque no se llevó a la validación cruzada por bloques ni quedó entre las variantes prometedoras.
- `kan_variantes_core.py`: componentes compartidos requeridos por los scripts D² y G=3. Lee los datos de `salidas_R_funcional/para_python/` y escribe en `KAN/salidas/`.
- `comparar_parametros_kan.py` y `evaluar_parametros_kan.R`: estudio controlado de parametros; entrena solo con entrenamiento, selecciona/parada con validacion y calcula en R el L² relativo normalizado.

## Comparación controlada de parámetros

Ejecuta una variable a la vez, manteniendo las demas en su valor de referencia, con las semillas 2026, 2027 y 2028:

| Parámetro | Valores comparados | Motivo |
|---|---|---|
| Ancho oculto | 3, 4, 5 | Comprobar capacidad de la capa oculta frente al ancho base 4. |
| Profundidad | 2 y 3 capas KAN | Aislar el efecto de añadir una capa; la profunda usa anchos `[10,4,3,6]`. |
| Rejilla `G` | 3, 5, 7 | Comparar resolución de las funciones spline. |
| Tasa de aprendizaje | 0,001; 0,003; 0,01 | Comprobar velocidad/estabilidad de Adam alrededor del valor base. |
| Paciencia | 200, 500, 1000 | Medir sensibilidad de la parada temprana. |
| Tope de épocas | 250, 500, 1000, 5000 | Es un máximo, no un número fijo de épocas; el mejor estado se elige con L² de validación. |
| `lambda_w` | 1e6, 1e7, 1e8 | Sensibilidad a la regularización de rugosidad del codificador. |
| `lambda_D2` | 0; 0,1; 1; 10; 100 | Sensibilidad a la penalización de segunda diferencia de las aristas. |

### Resultado de la corrida actual

Mejor valor observado dentro de cada eje, según el L² relativo normalizado por estación en validación (media ± DE de tres semillas):

| Eje | Mejor valor probado | L² relativo | Lectura |
|---|---:|---:|---|
| Ancho oculto | 5 | 0,406 ± 0,005 | Mejora frente al ancho 4 (0,423), con 888 parámetros totales frente a 728. |
| Capas KAN | 2 | 0,423 ± 0,005 | La alternativa de 3 capas dio 0,431 ± 0,029; no se justifica más profundidad. |
| Rejilla `G` | 3 | 0,415 ± 0,004 | Mejor que `G=5` (0,423) y `G=7` (0,430). |
| Tasa de aprendizaje | 0,003 | 0,416 ± 0,005 | `0,001` quedó prácticamente empatada (0,416 ± 0,006); `0,01` dio 0,423. |
| Paciencia | 200 | 0,423 ± 0,005 | Las tres paciencias dieron el mismo resultado; 200 redujo las épocas ejecutadas. |
| Tope de épocas | 250 | 0,423 ± 0,005 | Los mejores estados aparecieron en las épocas 107, 109 y 201; aumentar el tope no mejoró esta configuración. |
| `lambda_w` | 1e6 | 0,423 ± 0,005 | Los valores 1e7 y 1e8 empeoraron la validación. |
| `lambda_D2` | 100 | 0,403 ± 0,007 | Mejor valor probado; `lambda_D2=10` obtuvo 0,413 ± 0,003. |

Estos son mínimos por eje, **no una configuración conjunta ya probada**. Por ejemplo, el ancho 5 no se combinó con `G=3` y `lambda_D2=100`; esas interacciones requieren una corrida posterior en validación. Los resultados son para justificar y cribar decisiones, no una confirmación independiente.

Se mantienen fijos Adam, el codificador funcional de 8 dimensiones, los 6 coeficientes FPCA de salida, splines cúbicos y las particiones del proyecto. Son decisiones metodológicas/derivadas del flujo que se deben justificar, no cambiar inadvertidamente durante la comparación. Los escenarios idénticos se reutilizan para evitar entrenamientos duplicados.

Desde la raíz del proyecto:

```powershell
python -X utf8 .\KAN\codigo\comparar_parametros_kan.py
Rscript .\KAN\codigo\evaluar_parametros_kan.R
```

Los resultados quedan en `KAN/salidas/analisis_parametros_kan/`: entrenamientos, predicciones de validación, controles P6, resumen de L² funcional y métricas L² relativo en R. `mejor_valor_por_parametro_R.csv` muestra el mejor valor observado de cada eje. Esta es una comparación de cribado, no factorial: las interacciones entre parámetros no se estiman y el valor ganador de cada eje no debe combinarse automáticamente. El conjunto de prueba no se usa.

Desde la raiz del proyecto, ejecutar en PowerShell:

```powershell
$env:LLUVIA_SALIDA_KAN_RUGOSIDAD = "$PWD\KAN\salidas\propuesta_d2"
python -X utf8 .\KAN\codigo\lluvia_valle_red_kan_rugosidad.py

$env:LLUVIA_SALIDA_KAN_G3 = "$PWD\KAN\salidas\propuesta_g3"
python -X utf8 .\KAN\codigo\lluvia_valle_red_kan_g3.py

```

Los scripts seleccionan o paran usando entrenamiento/validacion y no evaluan la prueba. Se excluyo deliberadamente `lluvia_valle_red_kan_d2_prueba.py`: el conjunto de prueba ya fue consultado y esta reservado para reporte, no para seguir seleccionando propuestas.

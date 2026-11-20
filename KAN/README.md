# KAN: cabeza de la red funcional, MLP y auditoría de parámetros

Esta carpeta reúne todo el trabajo sobre la cabeza KAN de la red funcional: cómo se llegó a la **KAN oficial**, las variantes que se probaron, la comparación con el MLP (ablación), la interpretabilidad y la auditoría de cada parámetro.

**Regla común a todo:** solo se usan los años 1981–2000. La prueba 2001–2003 no se tocó para elegir nada. Los modelos se comparan con validación cruzada por 4 bloques de 5 años (1981–85, 86–90, 91–95, 96–00) y 3 semillas (2026, 2027, 2028). La métrica es el L² relativo de cada curva en mm/día, promediado por estación.

## Estructura

```
KAN/
├── README.md        este índice
├── codigo/          scripts de Python (y uno de R) que entrenan, analizan y grafican
├── informes/        un informe por etapa, en el orden en que se hicieron
└── salidas/         tablas resumen, figuras y modelos usados en las figuras
```

## El recorrido, en orden

| # | Etapa | Informe | Qué se concluyó |
|---|---|---|---|
| 0 | Propuestas de Federico | `informes/PROPUESTAS_FEDERICO_README.md` | Punto de partida: barrido de un parámetro a la vez con validación 1996–2000. Destacan G = 3 y la penalización D² |
| 1 | Experimentos complementarios (16-nov) | `informes/INFORME_EXPERIMENTOS_KAN.md` | Se combinan ancho, G y D², se prueba una KAN residual y se compara con MLP igualados. Con validación cruzada gana la KAN [10,4,6] G=5 D²=100 (0,429); la residual se descarta |
| 2 | Propuestas KAN nuevas (17-nov) | `informes/INFORME_PROPUESTAS_NUEVAS_KAN.md` | Estandarizar las neuronas ocultas corrige que la capa 2 caiga fuera de la rejilla. El promedio Federico + N6 mejora (0,415), pero es un ensamble |
| 3 | MLP y ablación oficial (17-nov) | `informes/INFORME_PROPUESTAS_MLP.md` | Se fija la **KAN oficial sin regularizar** (0,4397). 13 de 14 MLP del mismo tamaño la superan por poco: la interpretabilidad cuesta 0,007 (1,6 % del error) y la KAN conserva el 73 % de la ganancia no lineal sobre la cabeza lineal (0,4588) |
| 4 | Interpretabilidad de la KAN oficial (18-nov) | `informes/INFORME_INTERPRETABILIDAD_KAN_OFICIAL.md` | Gráficas G1–G8: las splines tienen forma (umbrales, campanas) y se identifica qué neurona lleva el nivel de lluvia, las temporadas y el ONI |
| 5 | Auditoría de parámetros (18-nov) | `informes/AUDITORIA_PARAMETROS_KAN.md` | 31 variantes de 13 parámetros, 372 ajustes. 10 de 13 parámetros quedan justificados. Candidatos a cambio: **ancho oculto 5** y **G = 3** |

## La KAN oficial

Script: `codigo/entrenar_kan_base_oficial.py` (de Federico; necesita GPU CUDA).

| Parámetro | Valor |
|---|---|
| Arquitectura de la cabeza | [10, 4, 6]: 8 neuronas funcionales + altitud + ONI → 4 ocultas → 6 puntajes FPCA |
| Rejilla de las splines | G = 5 intervalos, grado 3, rango [−3, 3], entradas estandarizadas |
| Regularización de la cabeza | ninguna (solo la rugosidad de las w_k(t), λ_w = 1e6, que pide el enunciado) |
| Entrenamiento | Adam, lr = 0,01, lote completo, hasta 5 000 épocas con paciencia 1 000 |
| Parámetros | 728 |
| L² relativo, validación cruzada | **0,4397** (DE entre semillas 0,0034) |

## Variantes más prometedoras

| Variante | L² CV | Por qué es prometedora | Por qué no se adopta todavía |
|---|---|---|---|
| Ancho oculto 5 | 0,4319 | Mejora en 3 de 4 bloques y es más estable | +22 % de parámetros; falta probarla combinada con G = 3 |
| G = 3 | 0,4342 | Mejora en los 4 bloques con menos parámetros (600) | Splines más rectas y 31 % de la capa 2 fuera de la rejilla |
| D² = 100 (KAN de Federico) | 0,4293 | La mejora más grande | Aplana las splines: la KAN pierde lo que la justifica |
| Federico + N6 (promedio de dos KAN) | 0,4152 | Gana en los 4 bloques, p < 0,001 | Es un ensamble de 1 520 parámetros, no una arquitectura |
| MLP M25 (3 × 14 + L2 = 200) | 0,4178 | El menor error de todos | Está regularizado: queda fuera de la ablación oficial |

El siguiente paso propuesto está en la sección 7 de la auditoría: una ronda de confirmación que combine ancho 5, G = 3 y el rango [−4, 4].

## Código

| Script | Etapa | Qué hace |
|---|---|---|
| `entrenar_kan_base_oficial.py` | oficial | Entrena la KAN oficial (escribe en `salidas/kan_oficial/`) |
| `kan_variantes_core.py` | todas | Datos, codificador, capa KAN, entrenamiento y métricas compartidos |
| `comparar_parametros_kan.py`, `evaluar_parametros_kan.R` | 0 | Barrido de Federico (un parámetro a la vez) y su métrica en R |
| `lluvia_valle_red_kan_rugosidad.py`, `lluvia_valle_red_kan_g3.py` | 0 | Las dos propuestas de Federico que siguieron adelante (D² y G = 3) |
| `experimentos_kan.py`, `analizar_experimentos_kan.py` | 1 | Combinaciones, KAN residual y referencias (validación 1996–2000) |
| `cv_bloques_kan.py`, `analizar_cv_bloques_kan.py` | 1 | Validación cruzada por bloques de los finalistas |
| `bsplines_kan.py` | 1–2 | Gráficas de las B-splines de cada capa |
| `cv_propuestas_nuevas_kan.py`, `analizar_propuestas_nuevas_kan.py`, `figuras_propuestas_nuevas_kan.py` | 2 | Propuestas N1–N7 y ensambles |
| `cv_mlp.py`, `analizar_mlp.py`, `figuras_mlp.py`, `analizar_ablacion_oficial.py` | 3 | Las 49 propuestas MLP (rondas r1–r5) y la ablación oficial |
| `interpretabilidad_kan.py`, `bspline_estilo_clase_kan.py` | 4 | Gráficas G1–G8 (`oficial` o `federico` como argumento) |
| `auditoria_parametros_kan.py`, `analizar_auditoria_kan.py` | 5 | Auditoría de parámetros y su verificación (paso 0) |
| `figura_unificada_kan.py` | 1–2 | Funciones de la capa 1: KAN base frente a D² = 100 (`bsplines/B2_...png`) |

Todos los scripts se corren desde `KAN/codigo/` con el entorno `.venv` del proyecto y leen los insumos de `salidas_R_funcional/para_python/`. Cada informe trae al final sus comandos.

## Salidas

- **En el repositorio:** tablas resumen (`*_resumen.csv`, `*_L2_por_bloque_y_semilla.csv`, contrastes), figuras (`fig1` a `fig9`, `auditoria_kan/fig_auditoria_parametros.png`, `interpretabilidad*/G*`, `bsplines/`) y los modelos que usan las figuras (`bsplines/*.pt`, `cv_nuevas/modelos/`).
- **Solo en local, fuera de Git** (ver `.gitignore`): las predicciones de cada entrenamiento (`corridas/`, `cv_bloques/`, `cv_mlp/`, `cv_nuevas/nv*.csv`, `auditoria_kan/corridas/`), los modelos intermedios y los logs, unos 190 MB. Los scripts de entrenamiento los regeneran; los de análisis los necesitan para recalcular las tablas.
- `interpretabilidad/` es la versión anterior, con la KAN D² = 100. Se conserva porque motivó las propuestas de la etapa 2.

## Nota sobre los insumos

Los entrenamientos se hicieron con `salidas_R_funcional/para_python/` antes de la decisión D-59 (corrección de la altitud de El Alcázar, de 17 a 977 m). Al recalcular las tablas con los insumos actuales, las diferencias son de 2 × 10⁻¹⁶, solo redondeo. Los modelos, sin embargo, vieron la altitud anterior en esa estación: un reentrenamiento final de la KAN oficial debe hacerse con los insumos actuales.

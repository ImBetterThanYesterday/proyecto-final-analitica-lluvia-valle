# Resumen del proyecto y decisiones pendientes (14-nov-2026)

Rama: `gus-cierre-funcional`. Detalle de cada cambio en `BITACORA_CAMBIOS.md`; de cada decisión en `DECISIONES.md` (D-01 a D-55); estado criterio por criterio en `CONTROL_RUBRICA.md`.

## 1. Qué está hecho

El análisis completo de la rúbrica (ítems 1 a 4.5) está programado, verificado y documentado. Falta redactar el informe.

| Ítem | Peso | Estado | Qué lo respalda |
|---|---|---|---|
| 1. Construcción de datos funcionales | 10 % | ✅ | B-splines paso a paso, K = 8 / 11 por ECM, VC y VCG (D-39, D-51) |
| 2. EDA funcional | 15 % | ✅ | Media verificada con `fda::mean.fd`, covarianzas, derivadas, profundidad, atípicas verificadas con `fbplot` |
| 3. Comparación por fuente | 15 % | ✅ | Emparejamiento documentado, d(t) con banda, correlación cruzada, por vertiente, altitud y ENSO |
| 4.1 Preprocesamiento | 6 % | 🟡 | Reglas R0 a R7, trazabilidad 323 → 220 → 176, sin fuga. **Falta justificar P-19** |
| 4.2 Red funcional | 10 % | ✅ | Gram por cuadratura, FPCA de entrenamiento (94,2 %), equivalencia demostrada y verificada (1,5e-14), no negativa |
| 4.3 Ablación | 11 % | ✅ | Mismo codificador y decodificador, KAN 640 vs MLP 635 parámetros, 3 semillas, mismo presupuesto, KAN verificada |
| 4.4 Evaluación | 8 % | ✅ | Cuatro métricas, IC bootstrap, Wilcoxon, por fase ENSO, contra el satélite sin corregir |
| 4.5 Interpretación | 10 % | ✅ | w_k(t), activaciones de la KAN, ajuste simbólico (tanh, R² 0,99), peores casos |
| 5. Informe | 10 % | ⬜ | Todo |
| 6. Scripts | 5 % | 🟡 | R unido en `lluvia_valle.R` (falta correrlo completo en la Mac); Python en un solo archivo |

Pruebas automáticas: 27 (datos) + 13 (funcional) + 5 (red en R) + 13 (Python), todas pasan.

## 2. Resultados clave (para el informe y la sustentación)

- **Prueba 2001-2003, L² relativo:** satélite sin corregir 0,513 → regresión 0,437 → red lineal 0,440, KAN 0,434, MLP 0,420. Sesgo anual del satélite +26 % → +7 % a +10 %.
- **La regresión es el método de clase:** frente a la versión del taller y a `fda.usc::fregre.basis.fr` da el mismo error y mejor KGE (D-53).
- **KAN y MLP ganan poco y depende del ENSO:** mejores en años Niño (1997, 2002), peores en años Niña (1999-2000). La prueba no tiene años Niña.
- **Qué aprendió la red:** KAN y MLP contrastan las dos temporadas de lluvia; la corrección satura cuando el satélite marca mucha lluvia; la altitud es la entrada más no lineal.
- **Peores casos:** estaciones secas del fondo del valle del Cauca y del Dagua; el píxel de 5 km no las representa.

## 3. Feedback

**Lo fuerte:**
- Cada decisión está respaldada por una tabla, una gráfica o una prueba en el código, como pidió el profesor.
- Python y R se verifican entre sí: Python reproduce la W, el FPCA y la regresión de R al decimal.
- Todo usa métodos de clase; la regresión se contrastó con las dos versiones vistas en clase.

**Lo que hay que cuidar:**
- Cada número del informe tiene que salir de una tabla del código (si no, el 4.2 baja a 1).
- La ventaja de KAN y MLP es pequeña: hay que presentarla como tal, no como un gran resultado.
- Cada integrante debe poder explicar cualquier fase (una respuesta mal dada resta 0,2). `CONTROL_RUBRICA.md` tiene la sección "Lo que hay que saber defender".

## 4. Decisiones para Jaco

| # | Decisión | Opciones | Recomendación |
|---|---|---|---|
| 1 | **P-19:** ¿usar los registros 1930-1979 para el FPCA? | (a) No. (b) Solo para el FPCA | **(a)**. La fe de erratas fija 1981-2003, el FPCA ya explica 94,2 %, y son datos de otra época sin depurar. Se respalda con un conteo en el código. **Es lo único que falta para el nivel 5 del 4.1** |
| 2 | **P-16:** ¿qué ONI anual? | (a) Media de los 12 trimestres (lo actual). (b) Un trimestre (DEF) | **(a)**, salvo que el profesor diga otra cosa |
| 3 | **Estaciones reservadas:** la fe de erratas aún pide un CSV, el profesor dijo que no van | Confirmar por escrito con el profesor | Preguntarle |
| 4 | **¿Qué R se entrega y cuál se edita?** Ahora hay tres: `lluvia_valle_datos.R`, `lluvia_valle_funcional.R` y el unido `lluvia_valle.R` (idéntico a los dos) | (a) Editar solo `lluvia_valle.R` y retirar los otros dos. (b) Seguir editando los dos y volver a unir al final | **(a)** después de correrlo completo en la Mac, para no tener copias que se desincronicen |
| 5 | **KAN propia o pykan** | (a) La nuestra (cumple la fórmula del enunciado y la prueba sintética). (b) Contrastar con pykan | **(a)**; (b) solo si sobra tiempo |
| 6 | **Detalles opcionales** de "como en clase" | Verificar FM con `fda.usc::depth.FM`; t funcional o ANOVA funcional para ENSO y altitud; justificar mm/día frente a acumulada; ventanas del Pacífico | Hacer al menos los dos primeros: son rápidos |
| 7 | **Pull request a `main`** | Cuándo fusionar `gus-cierre-funcional` | Cuando el equipo revise |

## 5. Lo que falta para entregar

1. Cerrar P-19 en el código (después de la decisión 1).
2. Correr `lluvia_valle.R` completo en la Mac (27 + 13 + 5 pruebas).
3. Informe PDF: una columna, ≤ 20 páginas, figuras y tablas numeradas, citadas e interpretadas.
4. Anexo en Word: herramienta de IA usada (Claude) y los prompts.
5. Guía de decisiones con diagramas para estudiar la sustentación.

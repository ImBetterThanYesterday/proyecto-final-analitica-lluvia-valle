"""Informe Word vigente desde CSV/PNG reales; no recalcula estadística.

Uso: python3 generar_informe_actual.py
Requiere pandoc. Para PDF/inspección: soffice --headless --convert-to pdf ...
Conserva intactos la guía histórica y ambos Word normativos.
"""
from pathlib import Path
import csv
import hashlib
import json
import subprocess
from datetime import datetime
from zipfile import ZipFile, ZIP_DEFLATED
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "resultados_R"
STEM = "INFORME_RECONSTRUCCION_EDA_ACTUAL"
W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
NS = {"w": W}
ET.register_namespace("w", W)
usados = set()
figuras = []
N_SECCIONES = 0


def leer(nombre):
    p = OUT / nombre
    usados.add(p)
    with p.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def num(x, dec=3):
    return f"{float(x):,.{dec}f}".replace(",", "@").replace(".", ",").replace("@", ".")


def entero(x):
    return num(x, 0)


def tabla(cabeceras, filas):
    def fila(xs):
        return "| " + " | ".join(str(x).replace("|", "/").replace("\n", " ") for x in xs) + " |"
    return "\n" + "\n".join([fila(cabeceras), fila(["---"] * len(cabeceras))] + [fila(f) for f in filas]) + "\n\n"


def figura(nombre, pie, ancho=6.1):
    p = OUT / (nombre + ".png")
    if not p.is_file():
        raise FileNotFoundError(f"Falta {p}; ejecute los bloques del R vigente.")
    usados.add(p)
    figuras.append(nombre)
    # Sin texto alternativo para que pandoc no duplique el pie explícito.
    return f"\n![](resultados_R/{nombre}.png){{width={ancho}in}}\n\n*Figura {len(figuras)}. {pie}*\n\n"


def pareja(nombres, pie):
    partes = []
    inicio = len(figuras) + 1
    for nombre in nombres:
        p = OUT / (nombre + ".png")
        if not p.is_file():
            raise FileNotFoundError(p)
        usados.add(p)
        figuras.append(nombre)
        partes.append(f"![](resultados_R/{nombre}.png){{width=3.05in}}")
    return "\n" + " ".join(partes) + f"\n\n*Figuras {inicio}–{len(figuras)}. {pie}*\n\n"


def pagina(titulo, texto):
    return f"# {titulo}\n\n{texto.strip()}\n\n"


SALTO = '\n```{=openxml}\n<w:p><w:r><w:br w:type="page"/></w:r></w:p>\n```\n\n'


def contenido():
    cobertura = leer("cobertura.csv")
    diag = {r["fuente"]: r for r in leer("diagnosticos_EDA.csv")}
    resumen = {n: leer(f"resumen_{n}.csv") for n in diag}
    controles = leer("verificacion_flujo.csv")
    bloques = leer("verificacion_bloques.csv")
    if len(controles) != 12 or any(r["estado"] != "OK" for r in controles):
        raise ValueError("Se requieren los 12 controles integrales en OK.")
    if len(bloques) != 14 or any(int(r["lineas_consola"]) == 0 for r in bloques):
        raise ValueError("Se requiere evidencia de salida de los 14 bloques.")
    if [int(r["curvas"]) for r in cobertura] != [4426, 33189, 31119]:
        raise ValueError("Cobertura inesperada: revisar antes de regenerar el informe.")
    if any(r["desde"] != "1981" or r["hasta"] != "2003" for r in cobertura):
        raise ValueError("Informe fuera del periodo normativo.")
    pags = []
    pags.append(pagina("Lluvia del Valle del Cauca", """
## Reconstrucción y análisis exploratorio funcional · 1981–2003

**Informe vigente y guía de estudio de las fases 1 y 2**

Este documento acompaña a `proyecto_funcional.R` y está escrito para que pueda
leerse **sin conocimientos previos** de estadística funcional. Explica, paso a
paso, qué problema se resolvió, con qué datos, qué se encontró, cómo se
construyeron las curvas y cómo se describieron. Sustituye como material vigente
de estas fases a la guía histórica; no sustituye la consigna del docente.

**Pregunta de trabajo.** ¿Cómo representar como curvas anuales las observaciones
diarias de estaciones y las observaciones pentadales satelitales, y qué revelan
sus formas, variabilidad y centralidad sobre el comportamiento de la lluvia
en la ventana común 1981–2003?

**Objetivo.** Reconstruir ambas fuentes con bases B-spline y mínimos cuadrados
programados en R estándar; después describir sus curvas mediante media,
dispersión, asociación entre tiempos, derivadas, profundidades y boxplot.

**Resultados principales.** K elegido por VCG: 6 en estaciones y 11 en satélite.
Se analizan 4.426 curvas estación–año y 31.119 curvas píxel–año en tierra;
el ajuste satelital incluye primero 33.189 curvas de la caja. El boxplot MBD
señala 170 y 0 curvas, respectivamente. La media se contrasta con `fda`.

**Advertencia de lectura.** Tener datos en los 23 años no equivale a tener
un panel completo ni una muestra representativa. Las unidades, identidad
y calendario del producto satelital siguen pendientes de confirmación.

**Cómo estudiar.** Cada sección tiene un título numerado; los diagramas de
decisión muestran el camino elegido y los caminos descartados con su razón.
Empiece por «El problema y las dos fases, desde cero», siga con «Datos» y los
glosarios, y después recorra la construcción y el EDA en orden. Las tablas
completas y figuras individuales están en `resultados_R/`.

**Alcance.** Guía de estas dos fases, no informe global terminado: no incluye
comparación pareada, regresión ni redes. La fe de erratas vigente fija 1981–2003;
no se mezcla esta ventana con los resultados históricos hasta 2012/2013.
"""))
    pags.append(pagina("0. El problema y las dos fases, desde cero", """
**Imagínese que no sabe nada del tema.** Tiene una pregunta muy concreta:
¿cómo se comporta la lluvia a lo largo de un año en el Valle del Cauca, y cómo
cambia ese comportamiento entre un año y otro?

**Qué se mide.**
- Una **estación** es un aparato en un punto del terreno que registra cuántos
  milímetros de lluvia caen cada día. De una estación y un año se obtiene una
  lista de hasta 365 valores diarios.
- Un **producto satelital** estima la lluvia sobre una cuadrícula (píxeles) a
  partir de sensores remotos. Aquí viene resumido en **pentadas** (bloques de
  unos cinco días), es decir, 72 valores por año y píxel, no 365.
- La **máscara de tierra** indica qué píxeles caen sobre continente; los demás
  son océano y no se usan para hablar de tierra.

**Qué se quiere lograr.** No basta con mirar listas de números. Queremos una
**curva anual** por cada estación–año y por cada píxel–año: una función suave
x(t) que diga cuánta lluvia cae, en promedio, en cada día del año. A eso se le
llama **dato funcional**: pasamos de observaciones dispersas a una función
estimada. Importante: la curva **no son datos nuevos medidos**; es una
representación construida con los datos que sí existen.

**Por qué hacerlo en dos fases.** Primero hay que **construir** las curvas
(fase de construcción funcional). Recién después tiene sentido **describirlas**
(fase de EDA, análisis exploratorio funcional): promedio, dispersión,
correlación entre días, pendientes, curvas atípicas.

**Qué SÍ se hace aquí.**
1. Filtrar ambas fuentes a la misma ventana: **1981–2003**.
2. Convertir cada serie anual en una curva suave, programando el método en R
   estándar (sin paquetes para construir).
3. Describir el conjunto de curvas con medidas funcionales y gráficos.
4. Contrastar la media calculada con la que da la librería `fda`, como control.

**Qué NO se hace (y por qué).** No se emparejan estaciones con píxeles (faltan
coordenadas y metadatos), no se estima un modelo causal, no se afirma que la
muestra sea representativa de todo el territorio y no se juzga la calidad del
satélite. Esas tareas requieren datos y fases posteriores. Confundirlas con lo
que sí se hizo sería el error más común al leer este informe.

**Cómo se decide.** Cada decisión metodológica se documenta en un **diagrama**:
se dibuja el camino elegido (verde) y, al lado, los caminos descartados (gris)
con la razón por la que no se usaron. Así se ve no solo *qué* se hizo, sino
*por qué* y *qué alternativas* existían.
""" + figura("diagrama_flujo_general", "Flujo completo de las dos fases. La ventana 1981–2003 se aplica antes de elegir K o ajustar, para que ambas fuentes hablen del mismo periodo.", 6.1)))
    filas = []
    for r in cobertura:
        nombre = r["fuente"]
        anual = leer(f"cobertura_anual_{nombre}.csv")
        ausentes = leer(f"sitio_anio_ausentes_{nombre}.csv")
        filas.append([nombre, entero(r["curvas"]), entero(sum(int(a["completas"]) for a in anual)),
                      entero(sum(int(a["incompletas"]) for a in anual)), entero(len(ausentes)),
                      entero(sum(int(a["valores_faltantes"]) for a in anual))])
    pags.append(pagina("1. Datos, procedencia y cobertura", """
**Estaciones.** `fda_fase1/curvas_estaciones_mm.csv.gz`: observaciones preparadas,
claves STATION_ID/YEAR y 365 columnas diarias. La unidad de trabajo es una
estación en un año, no una estación aislada. No se reconstruye aquí la limpieza
desde crudos ausentes de esta raíz.

**Satélite.** `fda_fase1/curvas_satelite_mm.csv.gz`: claves año/fila/columna y
72 tiempos heredados por año. Péntada es un intervalo de aproximadamente cinco
días; 72 columnas no prueban por sí solas el calendario ni su agregación exacta.

**Tierra.** De `fda_fase2/coef_satelite_tierra.csv.gz` se usan solo las claves de
la máscara heredada, NO sus coeficientes Python. La máscara fue construida con
referencia a CHIRPS, pero esto no identifica inequívocamente el producto recibido.

**ONI.** `eda_fase0/oni_anual_fase.csv` aporta la fase por ONI medio anual del
prototipo; el notebook 02 documenta NOAA CPC y una copia local. Falta confirmar
versión/fecha exacta de descarga. No es aquí una clasificación oficial de episodios.

**Tabla 1. Cobertura real después del filtro 1981–2003.** NA cuenta celdas;
ausencias cuenta combinaciones sitio–año sin curva, respecto del catálogo observado.
""" + tabla(["Fuente", "Curvas", "Completas", "Incompletas", "Ausencias", "NA"], filas) + """
Se observan 238 estaciones; 238 × 23 = 5.474 combinaciones posibles respecto
del catálogo, frente a 4.426 presentes. No se afirma que todas las estaciones
debieran existir todos los años. Las 665 curvas incompletas contienen 2.871 NA:
se ajustan solo donde hay observación, sin reemplazarlos por cero.

Satélite: 1.443 píxeles por año en la caja, 1.353 en tierra; la completitud de
estas matrices no verifica los datos crudos. Se excluyen 12.987 curvas satelitales
posteriores a 2003 antes de seleccionar K o ajustar.

**Qué significa «cobertura», en simple.** Para cada fuente contamos cuatro cosas:
cuántas curvas hay (filas), cuántas están completas, cuántas tienen algún hueco
y cuántas combinaciones sitio–año del catálogo no aparecen. Un valor **NA** es un
día que no se midió: no es cero de lluvia. En estaciones hay 238 sitios posibles
en 23 años; si todos existieran siempre habría 5.474 combinaciones, pero hay
4.426. Eso no es necesariamente un error: una estación puede no haber operado
en ciertos años. Lo dejamos declarado en vez de inventarlo.

**Qué se encontró, leyendo la tabla.** La fuente satelital tiene muchas más
curvas que la terrestre porque cada año aporta cientos de píxeles, no decenas de
estaciones. Las 665 curvas incompletas de estaciones concentran 2.871 NA; se
ajustan solo en los días observados (ver la decisión de NA más adelante). Antes
de seleccionar K se descartan 12.987 curvas satelitales que caen fuera de 2003,
para no mezclar la ventana normativa con datos posteriores.

**De dónde sale cada insumo.** `curvas_estaciones_mm.csv.gz` y
`curvas_satelite_mm.csv.gz` contienen las matrices de observaciones ya
preparadas; `coef_satelite_tierra.csv.gz` aporta **solo las claves** de la
máscara terrestre (no sus coeficientes de Python); `oni_anual_fase.csv` trae la
fase del índice Niño por año. No se reconstruye aquí la limpieza desde los datos
crudos, porque esos archivos no están en esta carpeta.
""" + figura("diagrama_decision_fuentes", "Decisión sobre comparar las dos fuentes. Sin pares estación–píxel verificados, la comparación se declara pendiente en lugar de estimar un sesgo que los datos no permiten.", 6.1)))
    glosario_a = [
        ("Dato discreto", "Una observación en un tiempo concreto; las celdas de las matrices de entrada."),
        ("NA", "Marca de un dato faltante. No significa cero de lluvia y no se convierte en cero para ajustar."),
        ("Dato funcional / curva", "Una función x(t) que permite evaluar el fenómeno a lo largo del año. Es una representación estimada, no nuevos datos medidos."),
        ("Estación / píxel", "Punto de medición terrestre / celda espacial del producto satelital. No se suponen pares entre fuentes."),
        ("Péntada", "Intervalo cercano a cinco días. Duraciones, fechas y agregación del producto deben verificarse."),
        ("Base B-spline", "Colección de funciones locales que se combinan para representar una curva flexible."),
        ("Orden / grado", "Orden 4 significa polinomios de grado 3 por tramos; no significa K=4."),
        ("K", "Número de funciones base y de coeficientes por curva. Mayor K permite más detalle y puede sobreajustar."),
        ("Phi (Φ)", "Matriz n × K: en cada fila se evalúan las K bases en un tiempo observado."),
        ("Coeficiente", "Peso de una función base en la curva; los K pesos resumen una curva, pero no son observaciones de lluvia."),
        ("MCO", "Mínimos cuadrados ordinarios: elige pesos que minimizan la suma de residuos al cuadrado en los puntos observados."),
        ("Residuo", "Observación menos valor ajustado en ese mismo tiempo. No se calcula residuo donde faltaba el dato."),
        ("S / H", "Matriz suavizadora o hat: transforma observaciones en ajustes. En este MCO es una proyección, S²=S y tr(S)=K."),
        ("ECM", "Error cuadrático medio de ajuste en los datos usados. Su descenso al aumentar K no prueba mejor generalización."),
        ("VC", "Validación cruzada dejando un tiempo fuera. Aquí se calcula con el atajo por la diagonal de S."),
        ("VCG", "Validación cruzada generalizada: aproxima VC usando la traza de S. Su mínimo decide K en este trabajo."),
        ("Ventana común", "Periodo compartido por ambas fuentes: 1981–2003. Todo cálculo se restringe a él."),
        ("Nodo (interior)", "Punto donde se unen los tramos de una B-spline; controla dónde puede cambiar la forma de la curva."),
        ("Rango completo", "La matriz Phi tiene columnas linealmente independientes; sin esto el sistema de MCO no tiene solución única."),
        ("Sobreajuste", "Cuando un modelo es tan flexible que copia el ruido de los datos usados y generaliza peor."),
    ]
    pags.append(pagina("2. Glosario · de observaciones a curvas", tabla(["Término", "Significado para explicar el proyecto"], glosario_a) + """
**Ejemplo de traducción.** «K=6» quiere decir que cada curva de estación se
representa mediante seis funciones base y seis pesos. No quiere decir seis
estaciones, seis días ni seis años.
"""))
    glosario_b = [
        ("Banda puntual / EE / IC", "Curva ajustada ±2 errores estándar (EE). IC significa intervalo de confianza. Bajo sigma²I (varianza residual constante e independencia), la banda es puntual, no simultánea."),
        ("Media funcional", "Promedio de curvas evaluadas en un mismo día; se obtiene también promediando coeficientes."),
        ("Varianza / DE", "Dispersión entre curvas en cada día. DE es la raíz de la varianza y conserva la unidad de la curva."),
        ("CV", "Coeficiente de variación: DE/media. No tiene unidad; se vuelve inestable cerca de media cero."),
        ("Covarianza", "Cómo cambian conjuntamente los valores de dos días entre curvas. Depende de la escala y tiene unidad al cuadrado."),
        ("Correlación", "Covarianza normalizada entre −1 y 1, sin unidad; aquí relaciona días dentro de una fuente, no estaciones con píxeles."),
        ("Derivada", "Pendiente instantánea. Al derivar una tasa en mm/día respecto de días se obtiene cambio de tasa en mm/día²."),
        ("Profundidad", "Medida de centralidad de una curva respecto del conjunto; valor alto significa más central bajo esa definición."),
        ("FM", "Fraiman–Muniz: combina centralidad punto a punto. Se conserva escala del script depth.FM y también la escala FM_pdf."),
        ("MBD", "Profundidad de banda modificada: fracción de tiempos en que la curva cae entre pares de curvas; aquí incluye empates en el borde."),
        ("Mediana funcional", "Curva observada en el conjunto ajustado con profundidad máxima. Mediana FM y MBD pueden ser distintas."),
        ("Media recortada", "Media de curvas con FM al menos igual al cuantil 25 %. Incluye empates, por lo que no siempre retiene exactamente 75 %."),
        ("Boxplot funcional", "Envolvente del 50 % más profundo por MBD; límites expandidos 1,5 veces su amplitud y bigotes de las no atípicas."),
        ("Atípica", "Curva que toca o rebasa algún límite de detección según la regla implementada de clase. No equivale automáticamente a error."),
        ("ONI / ENSO", "ONI: Índice Niño Oceánico, usado como contexto de ENSO, El Niño–Oscilación del Sur. Una fase anual heredada no prueba causa local."),
        ("Banda puntual vs simultánea", "La puntual vale para un día concreto; la simultánea cubriría todos los días a la vez. Aquí solo se dibuja la puntual."),
        ("Empate", "Dos curvas con la misma profundidad. Afecta cómo se elige la mediana y cómo se aplica el recorte."),
    ]
    pags.append(pagina("3. Glosario · cómo interpretar el EDA", tabla(["Término", "Significado y cautela"], glosario_b)))
    pags.append(pagina("4. De los puntos a una representación anual", """
**La idea central, sin fórmulas.** Cada fila de la matriz es una serie en el
tiempo. Si la dibujamos tal cual, aparece una línea quebrada con huecos y ruido.
Queremos una curva suave que resuma su forma. Ese paso —de puntos a función— es
el corazón de la construcción funcional.

**Bloques 1–2.** Primero se filtran ambas fuentes a 1981–2003. Las filas son
curvas y las columnas tiempos, al revés de la notación matricial de varios
ejemplos docentes. Los identificadores (estación y año, o año/fila/columna)
siempre viajan con sus filas, para no perder a qué unidad pertenece cada curva.

**La malla de evaluación.** Evaluamos todas las curvas en los días 1 a 365, para
que se puedan comparar entre sí aunque unas tengan 365 días observados (estación)
y otras 72 pentadas (satélite). Esa malla no inventa datos: solo es el conjunto
de días donde se evalúa la curva estimada. Tampoco convierte 72 pentadas en 365
observaciones reales ni resuelve por sí misma los años bisiestos o los metadatos
satelitales pendientes; son límites que se declaran.
""" + figura("observados_estaciones", "Cuatro curvas estación–año, identificadas en la leyenda. Las irregularidades son datos discretos antes del ajuste.")
    + figura("observados_satelite", "Cuatro curvas píxel–año de la caja. Son ejemplos gráficos, no el subconjunto usado para estimar: el cálculo usa todas las filas elegibles.")
    + "**Lectura prudente.** Las dos figuras no representan los mismos lugares. No permiten inferir sesgo satelital ni equivalencia de unidades.\n"))
    seleccion = {n: leer(f"seleccion_K_{n}.csv") for n in ("estaciones", "satelite")}
    filas = []
    for n, datos in seleccion.items():
        for criterio in ("ECM", "VC", "VCG"):
            mejor = min(datos, key=lambda r: float(r[criterio]))
            filas.append([n, criterio, mejor["K"], num(mejor[criterio])])
    pags.append(pagina("5. Elegir K sin confundir ajuste con validación", """
**Qué es K, en palabras simples.** K es cuántas piezas (funciones base) usamos
para armar cada curva. Con K pequeño la curva es rígida; con K grande es tan
flexible que puede «copiar» el ruido. Elegir K es decidir cuánta flexibilidad
permitimos. Ni poco ni demasiado: hay que validarlo.

**La trampa del ajuste.** El ECM mide qué tan cerca pasa la curva de los puntos
usados. Si uno elige K bajando el ECM, siempre gana el K más grande y se
sobreajusta. Por eso se usan criterios que **castigan la flexibilidad** (VC y
VCG) y no solo el ajuste.

**Bloque 3.** Se examina K=4…40 y 45…100 de cinco en cinco en estaciones;
K=4…40 en satélite. Se selecciona con las 3.761 curvas completas de estaciones
y las 33.189 curvas satelitales de la ventana. Las incompletas se ajustan después.

**Qué es «validar dejando un tiempo fuera».** VC quita un día, ajusta con los
demás y predice el día quitado; VCG aproxima esa idea sin reajustar una vez por
cada día, usando la traza de la matriz suavizadora. El K que mejor predice (no
el que mejor ajusta) es el que se adopta.

**Tabla 2. Mínimos reales.** ECM y VC se promedian por curva y tiempo;
VCG conserva la suma sobre curvas del script docente. No compare sus magnitudes
entre fuentes como si fueran el mismo universo o una medida de sesgo.
""" + tabla(["Fuente", "Criterio", "K mínimo", "Valor"], filas)
    + pareja(["criterio_ECM_estaciones", "criterio_ECM_satelite"], "ECM por fuente. La línea vertical marca su mínimo; por sí solo no decide, porque decrece al aumentar K.")
    + pareja(["criterio_VC_estaciones", "criterio_VC_satelite"], "VC rápidas por fuente. Penaliza la flexibilidad y su mínimo suele quedar cerca del K adoptado.")
    + pareja(["criterio_VCG_estaciones", "criterio_VCG_satelite"], "VCG por fuente. Su mínimo es el K adoptado: 6 en estaciones y 11 en satélite. Las tres gráficas cubren la exigencia ECM, VC y VCG.")
    + figura("diagrama_decision_K", "Decisión de K: se descarta minimizar el ECM de ajuste (siempre baja) y fijar K a ojo; se usa VCG por su atajo de la traza de S. Resultado: 6 estaciones, 11 satélite.", 6.1)))
    phi = leer("Phi_estaciones.csv")
    filas = [[r["tiempo"]] + [num(r[f"X{i}"] if f"X{i}" in r else r[list(r)[i]], 4) for i in range(1, 7)] for r in phi[:4]]
    pags.append(pagina("6. Base, matriz Phi y coeficientes", """
**Qué es una B-spline, con una analogía.** Imagine una regla flexible partida en
tramos que se unen suavemente en puntos llamados **nodos**. Cada tramo es un
polinomio cúbico (grado 3). Las **funciones base** son piezas locales que, al
sumarlas con ciertos pesos (los **coeficientes**), producen la curva. Es como
armar una figura con K piezas móviles: los coeficientes dicen cuánto pesa cada
pieza, y K dice cuántas piezas hay (la flexibilidad).

**Por qué cúbica y con nodos equiespaciados.** El orden 4 (grado 3) es un
equilibrio clásico entre suavidad y flexibilidad; nodos igualmente espaciados
evitan introducir supuestos sobre dónde «debería» cambiar la lluvia. No se
optimiza la posición de los nodos: eso añadiría decisiones que la consigna no
pide y que complicarían la reproducibilidad.

**Bloques 4–5.** Se adoptan B-splines cúbicos (orden 4), nodos interiores
equiespaciados y dominio 1–365. La recursión de Cox-de Boor se programa en R
estándar. No se emplea `fda`, `splines` ni `smooth.basis` en construcción.

**De la matriz a los coeficientes.** Para una curva: `c = solve(t(Phi) %*% Phi, t(Phi)) %*% y`.
Con huecos se restringen Phi e y a los tiempos observados; se exige rango
completo y más observaciones que coeficientes. En el script C es N × K. Si la
matriz Phi no tuviera rango completo, el sistema no tendría solución única; por
eso se verifica antes de resolver.

**Qué se decidió y qué se evitó.** Programar la base y el MCO es requisito de la
consigna; usar un paquete habría ocultado el método. Los diagramas resumen esa
decisión y, aparte, cómo se tratan los datos faltantes.
""" + figura("bases_estaciones", "Seis funciones locales de la base terrestre. La curva ajustada es su combinación ponderada.", 5.8)
    + figura("bases_satelite", "Once funciones satelitales; más K no implica que el satélite sea más preciso.", 5.8)
    + "**Tabla 3. Primeras cuatro filas de Phi de estaciones (sin unidad).**\n"
    + tabla(["Día"] + [f"B{i}" for i in range(1, 7)], filas)
    + figura("diagrama_decision_construccion", "Decisión de construcción: se descarta usar paquetes o cambiar el grado; se programan la B-spline cúbica y el MCO por ecuaciones normales.", 6.1)
    + figura("diagrama_decision_NA", "Decisión sobre faltantes: no se rellena con cero ni se borra la curva; se ajusta solo en los días observados, con su propia máscara.", 6.1)))
    filas = []
    for n in ("estaciones", "satelite"):
        for r in leer(f"controles_H_{n}.csv"):
            filas.append([n, r["control"], f'{float(r["error_max"]):.2e}'])
    pags.append(pagina("7. Ajuste MCO y banda puntual", """
**Qué hace el MCO.** Busca los pesos que hacen que la curva pase lo más cerca
posible de los puntos observados, midiendo «cerca» con la suma de residuos al
cuadrado. Es una proyección: convierte las observaciones en valores ajustados.

**Qué es la matriz S.** `S` (o «hat») aplica esa proyección: dado un vector de
observaciones, devuelve el vector ajustado. Dos propiedades útiles la delatan:
es simétrica e idempotente (aplicarla dos veces da lo mismo que una), y su traza
es igual a K. Comprobarlas es un control de que la implementación es correcta,
no una propiedad del clima.

`S = Phi %*% solve(t(Phi) %*% Phi, t(Phi))`. Se comprueban simetría,
idempotencia y traza=K en la malla completa y en cada máscara incompleta.
Los errores numéricos son próximos a la precisión de cálculo, no errores climáticos.

**Qué es la banda.** La banda gris alrededor de la curva es la curva ajustada
más/menos dos errores estándar. Es una banda **puntual** (vale día por día) y
supone varianza residual constante e independencia. Dibujarla no demuestra que
esos supuestos se cumplan; solo visualiza la incertidumbre bajo ese modelo.

**Tabla 4. Máximo error de cada propiedad entre máscaras.**
""" + tabla(["Fuente", "Control", "Error máximo"], filas)
    + figura("curva_banda_residuos_estaciones", "Primera curva completa terrestre; identidad en banda_ejemplo_estaciones.csv. Puntos observados, línea azul ajustada y banda ±2 EE.", 5.8)
    + figura("curva_banda_residuos_satelite", "Análogo para el primer píxel–año completo. Banda bajo sigma²I: puntual, no simultánea ni intervalo de predicción.", 5.8)
    + "No se recortan valores negativos para embellecer el ajuste: el MCO no impone no negatividad.\n"))
    filas = []
    for n in ("estaciones", "satelite"):
        coef = leer(f"coef_{n}.csv")
        r = coef[0]
        clave = " / ".join(r[k] for k in (("STATION_ID", "YEAR") if n == "estaciones" else ("anio", "fila", "col")))
        filas.append([n, clave, num(r["c1"]), num(r["c2"]), num(r["sigma2"]), r["n_observado"]])
    pags.append(pagina("8. Residuos: qué incertidumbre describen", """
**Qué es un residuo, en simple.** Es la diferencia, en un día concreto, entre lo
que se midió y lo que la curva estima. Residuos pequeños indican que la curva
pasa cerca de los datos; no indican que la curva sea «verdadera». Donde no hubo
dato, no se calcula residuo.

Un residuo es el dato observado menos el valor de su propia curva ajustada.
La varianza residual se estima como `sigma² = suma(residuo²)/(n_observado−K)`.
Para el ejemplo, `EE(t)² = sigma² × suma((Phi(t) %*% mapa)²)`.
La banda usa independencia y varianza residual constante; no se afirma que esas
suposiciones hayan quedado verificadas por dibujarla.

El segundo momento residual de clase es `Sigma_res = t(Res) %*% Res / N`,
solo con curvas completas. No se presenta como covarianza muestral centrada
con divisor N−1. Estructura fuera de la diagonal advierte dependencia residual.
""" + pareja(["residual_estaciones", "residual_satelite"], "Segundo momento residual por fuente, en escala heredada al cuadrado. Los colores de cada panel tienen su propia escala.")
    + "**Tabla 5. Primera curva de cada fuente: extracto de coeficientes (escala heredada), sigma² y número de observados.**\n"
    + tabla(["Fuente", "Clave", "c1", "c2", "sigma²", "n obs."], filas)
    + """
Los archivos `coef_*.csv` guardan todos los pesos y claves; `banda_ejemplo_*.csv`
guarda la identidad exacta del ejemplo, curva, EE y límites en 365 días.
El RDS conserva además la máscara observada. Las evaluaciones en días sin dato
son reconstrucciones, no imputaciones que se escriban sobre la matriz original.
"""))
    filas = []
    for n, datos in resumen.items():
        for d in (15, 105, 200, 290, 365):
            r = datos[d-1]
            filas.append([n.replace("satelite_tierra", "satélite tierra"), d, num(r["media"]), num(r["de"]), num(r["cv"])])
    pags.append(pagina("9. Media funcional y contraste de librería", """
**Qué es la media funcional, en simple.** Es la curva promedio: para cada día se
promedian los valores de todas las curvas. Da la forma «típica» del año. Ojo:
cada estación–año pesa igual, sin importar cuántos días tenga observados; por
eso una curva con muchos huecos no domina el promedio.

**Cómo se calcula y por qué coincide.** Promediar los coeficientes y luego
evaluar la base es equivalente a evaluar cada curva y promediar en cada día,
porque la evaluación es lineal. El script hace ambas cosas y comprueba que
coinciden; esa es una verificación de implementación.

**Bloque 7.** La media se calcula como `mu = P %*% colMeans(C)` y se comprueba
contra el promedio de las curvas evaluadas. Solo en EDA se crea un objeto `fda`
con los coeficientes ya calculados y se contrasta con `fda::mean.fd`.
No se vuelven a construir las curvas mediante un paquete.
""" + figura("media_estaciones", "Media de las 4.426 curvas terrestres; cada estación–año recibe el mismo peso.", 5.8)
    + figura("media_satelite_tierra", "Media de las 31.119 curvas en tierra. Ambas formas permiten discutir dos épocas de tasas altas, no demostrar concordancia espacial.", 5.8)
    + "**Tabla 6. Días ilustrativos, no promedios mensuales.** Media y DE en mm/día heredados; CV sin unidad.\n"
    + tabla(["Fuente", "Día", "Media", "DE", "CV"], filas)))
    filas = []
    for n, datos in resumen.items():
        for variable in ("varianza", "de", "cv"):
            valido = [r for r in datos if r[variable]]
            r = max(valido, key=lambda r: float(r[variable]))
            filas.append([n, variable, r["dia"], num(r[variable])])
    pags.append(pagina("10. Varianza, DE y variación relativa", """
**Qué miden, en simple.** En cada día, la **varianza** dice cuánto se separan las
curvas entre sí (unidades al cuadrado). La **DE** es su raíz y vuelve a las
unidades de la lluvia, así que es más fácil de leer. El **CV** divide la DE por
la media y quita las unidades, pero se dispara cuando la media es casi cero.

**La confusión que hay que evitar.** Esto no es la incertidumbre de la media ni
el error de una medición: es **cuánto varían las curvas entre sí** ese día. Una
DE grande puede deberse a que hay sitios muy lluviosos junto a otros secos, no a
que el instrumento sea impreciso.

**Bloque 8.** Se usa varianza muestral de las curvas reconstruidas (divisor N−1),
DE como raíz directa y CV=DE/media si la media no es casi cero.
Esta es variabilidad entre curvas, no EE de una curva ni incertidumbre de la media.
`std.fd` suavizada es una alternativa de clase, no se declara idéntica a esta DE.
""" + pareja(["varianza_estaciones", "varianza_satelite_tierra"], "Varianza funcional por fuente, en (mm/día)² heredados. Es variabilidad entre curvas, no EE de una curva.")
    + pareja(["de_estaciones", "de_satelite_tierra"], "DE por fuente: raíz directa de la varianza. En estaciones el máximo en el día 365 exige cautela por el borde y la mezcla de sitios, años y máscaras de observación.")
    + pareja(["cv_estaciones", "cv_satelite_tierra"], "Coeficiente de variación CV=DE/media por fuente (adimensional). Una media pequeña puede producir CV alto sin una gran DE absoluta.")
    + "**Tabla 7. Máximos diarios reales.** Varianza en (mm/día)² heredados, DE en mm/día heredados, CV adimensional.\n"
    + tabla(["Fuente", "Variable", "Día máximo", "Máximo"], filas)
    + "Las tres medidas describen dispersión entre curvas; no son incertidumbre de la media ni error satelital.\n"))
    filas = [[n, num(r["cor_abril_octubre"], 6)] for n, r in diag.items()]
    pags.append(pagina("11. Covarianza y correlación entre días", """
**Qué significan, en simple.** La covarianza pregunta: cuando un día es lluvioso
en una curva, ¿tiende a serlo también otro día en la misma curva? La correlación
es la misma idea pero sin unidades, entre −1 y 1. Aquí se comparan **días dentro
de una misma fuente**; nunca estaciones contra píxeles.

**Cómo leer el mapa de calor.** Una franja roja indica que esos dos días suelen
moverse juntos entre curvas; una azul, que van en sentidos opuestos. La diagonal
siempre es 1 (cada día consigo mismo). Asociación no es causalidad.

Se centran los coeficientes: `Cc = C − media_coef`.
`Cov = P %*% (t(Cc) %*% Cc/(N−1)) %*% t(P)`.
Su diagonal debe coincidir con la varianza diaria. La correlación divide por
las DE de los dos días. Cada celda compara días s y t dentro de una fuente.
""" + pareja(["covarianza_estaciones", "covarianza_satelite_tierra"], "Covarianza terrestre y satelital. Unidad heredada al cuadrado; los rangos de color son propios de cada panel.")
    + pareja(["correlacion_estaciones", "correlacion_satelite_tierra"], "Correlación entre días; escala fijada de −1 a 1. Asociación temporal entre curvas, no correlación entre las dos fuentes.")
    + tabla(["Fuente", "Correlación día 105–290"], filas)
    + "La asociación positiva sugiere que curvas con tasas altas en un periodo también tienden a tenerlas en otro. No demuestra causalidad, orografía ni sesgo de medición.\n"))
    filas = []
    for n, datos in resumen.items():
        cambios = [f'{i}–{i+1}' for i in range(1, 365)
                   if float(datos[i-1]["derivada_media"]) > 0 >= float(datos[i]["derivada_media"])]
        filas.append([n, ", ".join(cambios)])
    pags.append(pagina("12. Derivadas y temporadas: qué sí se puede decir", """
**Qué es una derivada, en simple.** Es la pendiente de la curva: cuánto sube o
baja la tasa de lluvia de un día al siguiente. Si la pendiente pasa de positiva
a negativa, la curva alcanza un máximo local: un candidato a «época lluviosa».
Ojo: se deriva una **tasa** (mm/día) medida en días, así que las unidades de la
derivada son mm/día². Es una tasa de cambio, no una cantidad de lluvia.

**Por qué se usa la derivada analítica.** Al derivar la función base (un
polinomio) se obtiene la derivada exacta de la curva, sin aproximar con restas.
Las diferencias finitas son una aproximación que depende del paso elegido, y
`deriv.fd` vuelve a proyectar la curva, por lo que no da exactamente lo mismo.
La implementación se contrasta con `fda::eval.fd(..., Lfdobj=1)`.

**Bloque 9.** La derivada analítica se obtiene con la derivada de la base y los
coeficientes medios, sin diferencias finitas. Se verifica en pruebas contra
`fda::eval.fd(..., Lfdobj=1)`. `deriv.fd` reproyecta y no es idéntica.

Pendiente positiva: aumenta la tasa media; negativa: disminuye. Cruce de positivo
a negativo: candidato a máximo local. Aquí se deriva una tasa, no lluvia acumulada.
""" + figura("derivada_estaciones", "Cambios de tasa en estaciones. Los cruces no son, por sí solos, fechas meteorológicas de inicio o fin.", 5.8)
    + figura("derivada_satelite_tierra", "Satélite: además de máximos dominantes hay oscilaciones menores. La flexibilidad y los bordes afectan la lectura.", 5.8)
    + "**Tabla 8. Intervalos de días con cambio positivo → no positivo de la derivada de la media.**\n"
    + tabla(["Fuente", "Intervalos de máximos locales candidatos"], filas)
    + "El patrón dominante es compatible con dos épocas lluviosas. Sigue pendiente justificar y validar una regla docente de inicio/fin; no se inventa un detector.\n"
    + figura("diagrama_decision_temporadas", "Decisión sobre temporadas: se descartan diferencias finitas y deriv.fd como definición; se usa la derivada analítica. La regla de inicio/fin queda pendiente, no inventada.", 6.1)))
    filas = []
    for n in diag:
        pr = leer(f"profundidades_{n}.csv")
        claves = ("STATION_ID", "YEAR") if n == "estaciones" else ("anio", "fila", "col")
        for criterio in ("FM", "MBD"):
            r = max(pr, key=lambda r: float(r[criterio]))
            filas.append([n, criterio, " / ".join(r[k] for k in claves), num(r[criterio], 6)])
    pags.append(pagina("13. Profundidad, mediana y media recortada", """
**Qué es «profundidad», con una imagen.** Piense en muchas curvas dibujadas una
sobre otra. La profundidad mide qué tan «central» es una curva: una curva muy
profunda está rodeada por las demás; una poco profunda vive en el borde del
conjunto. La curva **mediana** es la más profunda. Es una curva real del
conjunto, no un promedio.

**Por qué hay dos profundidades.** FM (Fraiman–Muniz) combina la centralidad día
a día. MBD (banda modificada) cuenta la fracción de días en que la curva cae
entre pares de otras curvas. Pueden elegir medianas distintas, sobre todo con
empates; por eso se reportan ambas en lugar de fingir que hay una sola verdad.

**La media recortada.** Es una media que deja fuera las curvas menos profundas
(el 25 % menos central). Sirve de puente entre la media (sensible a extremos) y
la mediana (robusta). Aquí se aplica con empates, así que no siempre retiene
exactamente el 75 %.

**Bloque 10.** FM agrega centralidad univariada; MBD compara cada curva con
bandas entre pares. Se calculan con todas las curvas de cada universo EDA en
la malla diaria. No se usa una muestra de 200 curvas para estimar profundidad.

La mediana es la primera curva de profundidad máxima; el recorte FM incluye
los empates del cuantil 25 %. MBD usa bandas inclusivas: los empates no se
tratan con rangos medios. Se conserva la escala docente FM y la columna FM_pdf.
""" + figura("EDA_estaciones", "Media, mediana FM y media recortada terrestre. Su separación describe asimetría o influencia de curvas grandes; no las elimina.", 5.8)
    + figura("EDA_satelite_tierra", "Centralidad satelital. Una mediana es una curva concreta del conjunto, no necesariamente el promedio puntual.", 5.8)
    + "**Tabla 9. Curvas medianas identificadas por claves y profundidad.**\n"
    + tabla(["Fuente", "Criterio", "Clave mediana", "Profundidad"], filas)
    + figura("diagrama_decision_centralidad", "Decisión sobre cómo resumir la curva típica: se reportan media, mediana FM, mediana MBD y media recortada, más boxplot; no se usa solo la media.", 6.1)))
    atip = leer("atipicas_estaciones.csv")
    filas = [[r["STATION_ID"], r["YEAR"], num(r["FM"], 5), num(r["MBD"], 5)] for r in atip[:4]]
    pags.append(pagina("14. Boxplot funcional y curvas atípicas", """
**Qué es un boxplot funcional, en simple.** Es la versión «por curvas» del
diagrama de caja clásico. Se toma el 50 % más profundo del conjunto y se dibuja
su envolvente como una banda central. Esa banda muestra dónde vive la mayoría
«típica». Alrededor, unos límites ampliados marcan cuándo una curva se aleja
demasiado; las curvas que los tocan o rebasan se señalan.

**Qué significa «atípica» aquí.** Solo que una curva cumple la regla de la clase
para ser señalada. No es un diagnóstico de error ni una lista de datos malos:
puede ser un extremo real, un ajuste difícil o algo a revisar. Por eso se
conservan e identifican, no se eliminan.

**Bloque 11.** Se ordena por MBD, se toma el 50 % más profundo y se construye
su envolvente. Los límites se expanden 1,5 veces la amplitud central. Una curva
se señala si toca o rebasa un límite en algún día, conforme al código de clase.
Los bigotes son la envolvente de las curvas no señaladas, no esos límites.
""" + figura("boxplot_estaciones", f"Estaciones: {len(atip)} curvas señaladas sobre 4.426. Se conservan y se identifican; pueden reflejar extremos, ajustes o problemas a revisar.", 5.8)
    + figura("boxplot_satelite_tierra", "Satélite terrestre: 0 señaladas bajo esta regla global. No significa ausencia de extremos locales ni calidad validada.", 5.8)
    + "**Tabla 10. Primeras cuatro atípicas, no ranking de severidad.** Lista completa de 170 en `atipicas_estaciones.csv`.\n"
    + tabla(["Estación", "Año", "FM", "MBD"], filas)))
    filas = []
    for n in diag:
        ens = leer(f"media_ONI_anual_{n}.csv")
        for r in ens:
            if r["dia"] == "1":
                filas.append([n, r["fase"], entero(r["n_curvas"])])
    pags.append(pagina("15. Contexto por ONI anual", """
**Qué es el ONI, en simple.** El Índice Niño Oceánico resume si el Pacífico está
en una fase cálida (El Niño), fría (La Niña) o neutra. Aquí se usa una fase
**anual heredada** como contexto: agrupamos los años por fase y comparamos las
curvas de cada grupo.

**Qué se puede y qué no se puede concluir.** Se puede describir si, en estos
datos, los años de una fase tienen medias distintas. No se puede afirmar que el
ENSO cause esas diferencias: los grupos mezclan sitios y años, no se ajustan
otras variables y los años no son independientes. Una asociación descriptiva no
es un efecto causal.

**Bloque 12.** Se enlaza la fase heredada por año a cada curva reconstruida.
Los tamaños se cuentan por curvas, no por episodios ENSO independientes.
El promedio mezcla sitios y años y no ajusta covariables ni dependencia espacial.
""" + figura("ONI_estaciones", "Medias terrestres por fase ONI anual. Las diferencias son descriptivas y pueden estar influidas por composición del panel.", 5.8)
    + figura("ONI_satelite_tierra", "Medias satelitales en tierra por la misma regla anual. No es un efecto causal estimado de ENSO.", 5.8)
    + "**Tabla 11. Tamaños de los grupos usados en cada media.**\n"
    + tabla(["Fuente", "Fase anual heredada", "Curvas"], filas)
    + "Para estudiar episodios oficiales se necesitaría verificar la serie y la definición temporal pertinente; aquí no se añade esa inferencia.\n"))
    filas = []
    for n in diag:
        for r in leer(f"media_grupos_descriptivos_{n}.csv"):
            if r["dia"] == "1":
                filas.append([n, r["grupo"], entero(r["n_curvas"])])
    pags.append(pagina("16. Grupos descriptivos y lectura territorial", """
**De dónde salen los grupos, en simple.** No hay una división geográfica oficial
en estos archivos, así que se usan dos reglas heredadas del prototipo: para
estaciones, el prefijo del identificador (26 u otro); para satélite, un corte de
longitud aproximado. Son formas de **agrupar para mirar**, no cuencas ni zonas
orográficas verificadas.

**Por qué se rotulan con cautela.** Comparemos curvas de dos «grupos» no prueba
que la diferencia se deba al relieve: los grupos mezclan altitudes, suelos y
años distintos. La orografía es una hipótesis razonable para explorar, no una
explicación demostrada por estas gráficas.

Estaciones: prefijo del identificador 26 frente a otro prefijo. Satélite: corte
de longitud heredado −76,8. Estas reglas se conservan del prototipo y se
rotulan por lo que son: no delimitan cuencas ni zonas orográficas verificadas.
""" + figura("grupos_estaciones", "Comparación por prefijo de estación; antes de interpretar por relieve se necesitan localización y metadatos.", 5.8)
    + figura("grupos_satelite_tierra", "Comparación por corte longitudinal aproximado en la máscara terrestre. No se han emparejado estaciones con píxeles.", 5.8)
    + "**Tabla 12. Curvas por grupo.**\n"
    + tabla(["Fuente", "Regla de grupo", "Curvas"], filas)
    + "La orografía es una hipótesis contextual razonable para explorar, no una explicación demostrada por estos grupos. No se estima sesgo entre fuentes sin pares válidos.\n"))
    filas = []
    for n in diag:
        lib = leer(f"verificacion_media_libreria_{n}.csv")[0]
        r = diag[n]
        filas.append([n, f'{float(lib["error_max"]):.2e}', entero(r["puntos_negativos"]), num(r["minimo"], 4)])
    pags.append(pagina("17. Conclusiones, límites y reproducción", """
**El recorrido en una frase.** Convertimos listas de lluvia en curvas anuales
suaves, elegimos su flexibilidad con validación, y describimos esas curvas con
promedios, dispersión, asociación, pendientes y profundidades; en el camino
dejamos explícito qué alternativas se descartaron y por qué.

**Conclusión.** La construcción manual entrega curvas y coeficientes trazables
en 1981–2003. Las medias son compatibles con dos épocas anuales de tasas altas;
la dispersión, las medianas y el boxplot complementan esa síntesis. Ni estos
resúmenes ni ONI establecen causalidad, representatividad o concordancia espacial.

**Tabla 13. Control de media y límites físicos del ajuste.** Mínimos y puntos
negativos se refieren a la malla evaluada, no a datos de lluvia negativos observados.
""" + tabla(["Fuente", "Error vs fda", "Puntos <0", "Mínimo ajustado"], filas) + """
**Pendientes sustantivos.** Confirmar producto, unidad y calendario satelital;
recuperar crudos y metadatos; justificar representatividad y sesgo de cobertura;
validar temporadas de lluvia y su sensibilidad a K/bordes. Los NA no se borran
por ajustar curvas. MCO sin restricción permite negativos. Media ±2 DE no es IC.

**Reproducir y estudiar.** En RStudio, abrir `proyecto_funcional.R` y ejecutar
desde la primera línea con Cmd/Ctrl+Enter, avanzando de arriba abajo. Primero
se definen los auxiliares completos; luego construcción y EDA calculan y muestran.
Source ejecuta todo desde una sesión nueva. Inspeccionar `Phi_est`, `S_est`, `e`, `s`, `eda_est`
y `eda_sat` en Environment. Ver una tabla completa: `ver_tabla("cobertura.csv")`.
Batch: `Rscript proyecto_funcional.R`. Regenerar Word: `python3 generar_informe_actual.py`.

**Evidencia.** `test_bloques_R.R`: preparación y 13 pasos secuenciales, con salida al paso.
`test_fda_R.R`, `test_periodo_R.R`, `test_usabilidad_R.R` y `verificar_flujo_R.R`:
matemática, periodo, dispositivos y 12 controles integrales; entradas preservadas.
El manifiesto del informe guarda hashes de sus CSV/PNG y del código. Esto no
reemplaza una inspección de la GUI real de RStudio.

**Fuentes.** Fe de erratas: objetivos y rúbrica 1/2/4.4; scripts docentes
`../analitica/Script_ConstrucciÃ³n_FDA.R` (73–105, 215–305, 350–387) y
`Script_AnalisisExploratorio.R` (155–196, 253–289, 351–364, 445–518);
PDF Construcción pp. 3–5, 25–36 y PDF EDA pp. 20–62 (páginas físicas).
Para el uso de observados: `Script_t-Funcional.R` 129–143, allí penalizado.
Procedencia ONI documentada: https://cpc.ncep.noaa.gov/data/indices/oni.ascii.txt.
Trazabilidad completa en `REQUISITOS_Y_ENTREGA_R.md` y CSV de métodos.

**Autoría y alcance.** Adaptación propia apoyada por OpenCode, modelo
`openai/gpt-6.1-sol`; no se atribuye cada línea al docente. Debe completarse el
anexo Word de IA con los prompts efectivamente conservados; no se inventa un
historial ni se declara terminado ese requisito global. Esta guía de estudio
no es el informe final completo ni implementa las fases posteriores.
"""))
    global N_SECCIONES
    N_SECCIONES = len(pags)
    assert N_SECCIONES >= 18
    return SALTO.join(pags)


def formato(path):
    with ZipFile(path) as z:
        partes = {n: z.read(n) for n in z.namelist()}
    estilos = ET.fromstring(partes["word/styles.xml"])
    for estilo in estilos.findall("w:style", NS):
        ident = estilo.get(f"{{{W}}}styleId")
        if ident in ("Normal", "BodyText", "FirstParagraph", "Compact", "Table"):
            rp = estilo.find("w:rPr", NS)
            if rp is None:
                rp = ET.SubElement(estilo, f"{{{W}}}rPr")
            fuentes = ET.SubElement(rp, f"{{{W}}}rFonts")
            for k in ("ascii", "hAnsi", "cs"):
                fuentes.set(f"{{{W}}}{k}", "Calibri")
            sz = rp.find("w:sz", NS)
            if sz is None:
                sz = ET.SubElement(rp, f"{{{W}}}sz")
            sz.set(f"{{{W}}}val", "20")
    partes["word/styles.xml"] = ET.tostring(estilos, encoding="utf-8", xml_declaration=True)
    doc = ET.fromstring(partes["word/document.xml"])
    for sec in doc.findall(".//w:sectPr", NS):
        referencia = ET.Element(f"{{{W}}}footerReference")
        referencia.set(f"{{{W}}}type", "default")
        referencia.set("{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id", "rIdInformeFooter")
        sec.insert(0, referencia)
        tam = sec.find("w:pgSz", NS)
        if tam is None:
            tam = ET.SubElement(sec, f"{{{W}}}pgSz")
        tam.set(f"{{{W}}}w", "11906"); tam.set(f"{{{W}}}h", "16838")
        marg = sec.find("w:pgMar", NS)
        if marg is None:
            marg = ET.SubElement(sec, f"{{{W}}}pgMar")
        for k in ("top", "bottom", "left", "right"):
            marg.set(f"{{{W}}}{k}", "1000")
        marg.set(f"{{{W}}}footer", "400")
    pie = ET.fromstring(f'''<w:ftr xmlns:w="{W}"><w:p><w:pPr><w:jc w:val="center"/></w:pPr>
        <w:r><w:rPr><w:sz w:val="16"/></w:rPr><w:t xml:space="preserve">Guía de estudio · Construcción y EDA · </w:t></w:r>
        <w:fldSimple w:instr="PAGE"><w:r><w:t>1</w:t></w:r></w:fldSimple></w:p></w:ftr>''')
    partes["word/footer-informe.xml"] = ET.tostring(pie, encoding="utf-8", xml_declaration=True)
    rels = ET.fromstring(partes["word/_rels/document.xml.rels"])
    ET.SubElement(rels, "{http://schemas.openxmlformats.org/package/2006/relationships}Relationship",
        Id="rIdInformeFooter", Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/footer",
        Target="footer-informe.xml")
    ET.register_namespace("", "http://schemas.openxmlformats.org/package/2006/relationships")
    partes["word/_rels/document.xml.rels"] = ET.tostring(rels, encoding="utf-8", xml_declaration=True)
    tipos = ET.fromstring(partes["[Content_Types].xml"])
    ET.SubElement(tipos, "{http://schemas.openxmlformats.org/package/2006/content-types}Override",
        PartName="/word/footer-informe.xml", ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.footer+xml")
    ET.register_namespace("", "http://schemas.openxmlformats.org/package/2006/content-types")
    partes["[Content_Types].xml"] = ET.tostring(tipos, encoding="utf-8", xml_declaration=True)
    for t in doc.findall(".//w:tbl", NS):
        grid = t.find("w:tblGrid", NS)
        cols = list(grid)
        # Todas las tablas tienen ancho explícito y cada celda coincide con grid.
        anchos = [int(c.get(f"{{{W}}}w")) for c in cols]
        total = sum(anchos)
        anchos = [round(a / total * 9800) for a in anchos]
        if "Término" in "".join(t.find("w:tr", NS).itertext()):
            anchos = [2300, 7500]
        anchos[-1] += 9800 - sum(anchos)
        for c, ancho in zip(cols, anchos):
            c.set(f"{{{W}}}w", str(ancho))
        pr = t.find("w:tblPr", NS)
        ancho_tabla = pr.find("w:tblW", NS)
        ancho_tabla.set(f"{{{W}}}type", "dxa"); ancho_tabla.set(f"{{{W}}}w", "9800")
        for fila in t.findall("w:tr", NS):
            fp = fila.find("w:trPr", NS)
            if fp is None:
                fp = ET.SubElement(fila, f"{{{W}}}trPr")
            ET.SubElement(fp, f"{{{W}}}cantSplit")
            for celda, ancho in zip(fila.findall("w:tc", NS), anchos):
                cp = celda.find("w:tcPr", NS)
                cw = cp.find("w:tcW", NS)
                if cw is None:
                    cw = ET.SubElement(cp, f"{{{W}}}tcW")
                cw.set(f"{{{W}}}type", "dxa"); cw.set(f"{{{W}}}w", str(ancho))
    partes["word/document.xml"] = ET.tostring(doc, encoding="utf-8", xml_declaration=True)
    with ZipFile(path, "w", ZIP_DEFLATED) as z:
        for nombre, datos in partes.items():
            z.writestr(nombre, datos)


def main():
    texto = contenido()
    md = ROOT / (STEM + ".md")
    md.write_text(texto, encoding="utf-8")
    destino = ROOT / (STEM + ".docx")
    subprocess.run(["pandoc", str(md), "--from=markdown", "--resource-path=" + str(ROOT),
                    "-M", "lang=es-CO", "-o", str(destino)], check=True)
    formato(destino)
    with ZipFile(destino) as z:
        assert z.testzip() is None
        for n in z.namelist():
            if n.endswith(".xml"):
                ET.fromstring(z.read(n))
        doc = ET.fromstring(z.read("word/document.xml"))
        imagenes = [n for n in z.namelist() if n.startswith("word/media/")]
        assert len(imagenes) == len(figuras)
        assert len(doc.findall(".//w:tbl", NS)) == 16  # 13 numeradas + correlación + 2 glosarios
        headings = [p for p in doc.findall(".//w:p", NS)
                    if p.find('w:pPr/w:pStyle[@w:val="Heading1"]', NS) is not None]
        assert len(headings) == N_SECCIONES
    usados.update([ROOT / "proyecto_funcional.R", Path(__file__),
                   ROOT / "Proyecto_Final_Analitica_FeErratas.docx"])
    manifiesto = {"generado": datetime.now().astimezone().isoformat(), "documento": destino.name,
                  "figuras": len(figuras), "tablas": 16, "secciones": N_SECCIONES,
                  "archivos": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in sorted(usados)}}
    (OUT / "manifiesto_informe_actual.json").write_text(json.dumps(manifiesto, indent=2, ensure_ascii=False) + "\n")
    # Correspondencia explícita figura R -> archivo -> número de figura en el Word.
    import csv as _csv
    incluidas = {n: i for i, n in enumerate(figuras, 1)}
    pngs = sorted(p.stem for p in OUT.glob("*.png")
                  if not p.name.startswith("pagina-") and p.stem != STEM)
    filas_fig = []
    for nombre in pngs:
        num = incluidas.get(nombre)
        if num:
            estado = "incluida"
        elif nombre.startswith("seleccion_K_") or nombre.startswith("contexto_"):
            estado = "compuesta_duplicada"
        else:
            estado = "solo_exportada"
        filas_fig.append([nombre, f"resultados_R/{nombre}.png", num or "", estado])
    with (OUT / "manifiesto_figuras.csv").open("w", encoding="utf-8", newline="") as f:
        w = _csv.writer(f)
        w.writerow(["figura_R", "archivo", "figura_word", "estado"])
        w.writerows(filas_fig)
    print(f"Generado {destino}: {len(figuras)} imágenes, 16 tablas, {N_SECCIONES} secciones; ZIP/XML verificados.")


if __name__ == "__main__":
    main()

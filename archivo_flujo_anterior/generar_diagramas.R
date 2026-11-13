# Diagramas de decisión del proyecto de lluvia (R base, sin paquetes).
# Uso: Rscript generar_diagramas.R [carpeta_salida]
# Cada diagrama explica el camino elegido y los caminos descartados con su razón.
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else file.path(getwd(), "resultados_R")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

VERDE <- "#d8f0d8"; GRIS <- "#ececec"; AMAR <- "#fdf3c8"; AZUL <- "#dbe8fb"
abrir <- function(nombre, alto = 5.2) {
  png(file.path(outdir, paste0(nombre, ".png")),
      width = 1342, height = round(220 * alto), res = 220, pointsize = 11)
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new(); plot.window(xlim = c(0, 100), ylim = c(0, 100))
}
cerrar <- function() dev.off()
caja <- function(x, y, w, h, txt, fill = "#ffffff", border = "black", lty = 1,
                 cex = 0.72, font = 1, col = "black") {
  rect(x - w/2, y - h/2, x + w/2, y + h/2, col = fill, border = border,
       lty = lty, lwd = 1.5)
  text(x, y, txt, cex = cex, font = font, col = col)
}
flecha <- function(x0, y0, x1, y1, col = "black", lty = 1, lwd = 1.5) {
  arrows(x0, y0, x1, y1, length = 0.07, col = col, lty = lty, lwd = lwd)
}
titulo <- function(t, sub = NULL) {
  text(50, 96, t, cex = 0.98, font = 2)
  if (!is.null(sub)) text(50, 92, sub, cex = 0.68, col = "#444444")
}
leyenda <- function() {
  text(3, 2.5, "verde = elegido   ·   gris = descartado   ·   amarillo = pendiente",
       adj = 0, cex = 0.62, col = "#333333")
}

# 1. Flujo general de las dos fases -----------------------------------------
abrir("diagrama_flujo_general", 6.0)
titulo("Flujo completo: de los datos a la interpretación",
       "Dos fases: construcción funcional y análisis exploratorio funcional")
caja(50, 84, 74, 8, "Datos preparados: matrices de observaciones diarias/pentadales\n(estaciones, satélite y máscara de tierra)", AZUL, cex = 0.7)
flecha(50, 80, 50, 74)
caja(50, 68, 74, 9, "FASE 1 · Construcción funcional: B-spline cúbica + MCO propio\n→ coeficientes por curva y curva evaluada en 365 días", VERDE, cex = 0.7)
flecha(50, 63.5, 50, 57.5)
caja(50, 51, 74, 9, "FASE 2 · EDA funcional: media, dispersión, asociación,\nderivadas, profundidades, boxplot y contexto", VERDE, cex = 0.7)
flecha(50, 46.5, 50, 40.5)
caja(50, 34, 74, 9, "Interpretación prudente: formas, épocas lluviosas y límites\n(no causalidad, no sesgo satelital, no muestra representativa)", AMAR, cex = 0.7)
text(50, 20, "Decisión transversal: ventana común 1981–2003 en ambas fuentes; nunca se mezcla con resultados históricos hasta 2012/2013.",
     cex = 0.64, col = "#333333")
leyenda(); cerrar()

# 2. Cómo elegir K -----------------------------------------------------------
abrir("diagrama_decision_K", 5.8)
titulo("Decisión 1 · ¿Cómo elegir el número de bases K?",
       "K = flexibilidad de la curva; debe decidirse por validación, no por ajuste")
caja(50, 86, 40, 7, "Se prueban varios K y hay que escoger uno", AZUL, cex = 0.72)
# rama izquierda: ECM
flecha(40, 82.5, 22, 77); caja(19, 70, 30, 13,
  "Minimizar ECM de ajuste\nDESCARTADO: siempre baja al\nsubir K ⇒ sobreajuste", GRIS, lty = 2, cex = 0.64)
# rama central: VC
flecha(50, 82.5, 50, 77); caja(50, 70, 30, 13,
  "VC dejando un tiempo fuera\nválida pero cara de calcular;\nse usa su atajo por diagonal de S", VERDE, cex = 0.64)
# rama derecha: fijar a ojo
flecha(60, 82.5, 78, 77); caja(81, 70, 30, 13,
  "Fijar K a ojo/experto\nDESCARTADO: no reproducible\nni justificable en la rúbrica", GRIS, lty = 2, cex = 0.64)
flecha(50, 63.5, 50, 58)
caja(50, 51, 52, 10, "VCG (validación cruzada generalizada):\naproxima la VC con la traza de S y penaliza la flexibilidad", VERDE, cex = 0.68)
flecha(50, 46, 50, 40.5)
caja(50, 33, 52, 10, "ELEGIDO: K por mínimo de VCG\n→ 6 en estaciones y 11 en satélite", VERDE, font = 2, cex = 0.72)
text(50, 20, "Por qué VCG y no VC exacta: el atajo por la diagonal de S evita reajustar una vez por cada tiempo,\nsin cambiar el criterio de escoger el K que menos sobreajusta.",
     cex = 0.62, col = "#333333")
leyenda(); cerrar()

# 3. Cómo construir las curvas ----------------------------------------------
abrir("diagrama_decision_construccion", 5.8)
titulo("Decisión 2 · ¿Cómo construir cada curva anual?",
       "La consigna exige programar la construcción, no llamar un paquete")
caja(50, 86, 52, 7, "Representar una curva suave a partir de datos discretos", AZUL, cex = 0.7)
flecha(40, 82.5, 22, 77); caja(19, 70, 30, 13,
  "fda::smooth.basis,\nsplines::bs o lm()\nDESCARTADO en construcción:\ncaja negra del paquete", GRIS, lty = 2, cex = 0.62)
flecha(50, 82.5, 50, 77); caja(50, 70, 32, 13,
  "B-spline cúbica (orden 4),\nnodos equiespaciados y\nrecursión de Cox–de Boor propia", VERDE, cex = 0.62)
flecha(60, 82.5, 78, 77); caja(81, 70, 30, 13,
  "Grado/orden distinto\nDESCARTADO: más flexibilidad\nlocal y bordes inestables", GRIS, lty = 2, cex = 0.62)
flecha(50, 63.5, 50, 58)
caja(50, 51, 56, 10, "MCO por ecuaciones normales:\nc = solve(t(Phi)Phi, t(Phi)y)", VERDE, cex = 0.68)
flecha(50, 46, 50, 40.5)
caja(50, 33, 56, 10, "ELEGIDO: base + MCO propios, con\nsimetría, idempotencia y traza=K verificadas", VERDE, font = 2, cex = 0.66)
text(50, 19, "Se usa fda SOLO para contrastar la media ya calculada (no para construir).\nSi la matriz Phi tiene rango completo, no hace falta QR/SVD; se comprueba numéricamente.",
     cex = 0.62, col = "#333333")
leyenda(); cerrar()

# 4. Qué hacer con los NA ----------------------------------------------------
abrir("diagrama_decision_NA", 5.6)
titulo("Decisión 3 · ¿Qué hacer con los datos faltantes (NA)?",
       "Un NA no es «no llovió»: es un dato que no se midió")
caja(50, 86, 46, 7, "Algunas curvas tienen días sin observación", AZUL, cex = 0.7)
flecha(38, 82.5, 20, 77); caja(18, 70, 30, 13,
  "Rellenar con 0\nDESCARTADO: inventa\n«no llovió» y sesga a la baja", GRIS, lty = 2, cex = 0.62)
flecha(50, 82.5, 50, 77); caja(50, 70, 30, 13,
  "Ajustar solo en los\ntiempos observados (máscara),\nexigiendo n>K y rango completo", VERDE, cex = 0.62)
flecha(62, 82.5, 80, 77); caja(82, 70, 30, 13,
  "Borrar toda la curva\nDESCARTADO: se pierden\n665 curvas con parte observada", GRIS, lty = 2, cex = 0.62)
flecha(50, 63.5, 50, 58)
caja(50, 51, 64, 10, "Cada curva conserva su propia máscara de días observados;\nla evaluación en días sin dato es reconstrucción, no imputación", VERDE, cex = 0.66)
flecha(50, 46, 50, 40.5)
caja(50, 33, 64, 10, "ELEGIDO: ajustar con máscara\n→ 4.426 curvas estación–año y 31.119 píxel–año en tierra", VERDE, font = 2, cex = 0.66)
text(50, 19, "Los 2.871 NA de las 665 curvas incompletas no se borran ni se reemplazan por cero:\nse ajusta donde hay dato y se declara que fuera de eso es una curva reconstruida.",
     cex = 0.62, col = "#333333")
leyenda(); cerrar()

# 5. Resumir la curva típica -------------------------------------------------
abrir("diagrama_decision_centralidad", 6.0)
titulo("Decisión 4 · ¿Cómo resumir la «curva típica»?",
       "Cada resumen responde a una pregunta distinta; se reportan varios")
caja(50, 88, 48, 6, "Resumir muchas curvas en pocas curvas guía", AZUL, cex = 0.7)
flecha(30, 85, 18, 80); caja(18, 74, 28, 10,
  "Media funcional\nELEGIDA: referencia, pero\nsensible a extremos", VERDE, cex = 0.62)
flecha(45, 85, 42, 80); caja(42, 74, 28, 10,
  "Mediana FM\nELEGIDA: curva observada\ncentral por profundidad", VERDE, cex = 0.62)
flecha(58, 85, 62, 80); caja(62, 74, 28, 10,
  "Mediana MBD\nELEGIDA: banda entre pares;\npuede diferir de FM", VERDE, cex = 0.62)
flecha(72, 85, 82, 80); caja(82, 74, 28, 10,
  "Media recortada FM 25 %\nELEGIDA: robusta\ny conserva empates", VERDE, cex = 0.62)
flecha(50, 88, 30, 76)
text(30, 63, "DESCARTADO: usar solo la media\noculta la influencia de curvas extremas", cex = 0.6, col = "#555555")
flecha(50, 60, 50, 55)
caja(50, 47, 68, 10, "Boxplot funcional MBD: envolvente del 50 % más profundo y\nlímites expandidos 1,5× la amplitud central para señalar atípicas", VERDE, cex = 0.64)
flecha(50, 42, 50, 37)
caja(50, 30, 68, 9, "ELEGIDO: reportar media, medianas, recortada y boxplot\n→ 170 atípicas en estaciones y 0 en satélite (regla global)", VERDE, font = 2, cex = 0.66)
text(50, 18, "«Atípica» = toca o rebasa un límite de la regla de clase; NO equivale a error ni a dato malo.\nLas atípicas se conservan e identifican, no se eliminan.",
     cex = 0.62, col = "#333333")
leyenda(); cerrar()

# 6. Temporadas y máximos ----------------------------------------------------
abrir("diagrama_decision_temporadas", 5.7)
titulo("Decisión 5 · ¿Cómo detectar épocas y máximos?",
       "Se deriva una tasa (mm/día), no una lluvia acumulada")
caja(50, 87, 46, 7, "¿Cuándo sube o baja la tasa media anual?", AZUL, cex = 0.7)
flecha(38, 83.5, 20, 78); caja(18, 71, 30, 13,
  "Diferencias finitas\nDESCARTADO: aproxima y\ndepende del paso elegido", GRIS, lty = 2, cex = 0.62)
flecha(50, 83.5, 50, 78); caja(50, 71, 30, 13,
  "Derivada analítica de la base\n+ coeficientes medios\n(se contrasta con fda, Lfdobj=1)", VERDE, cex = 0.62)
flecha(62, 83.5, 80, 78); caja(82, 71, 30, 13,
  "fda::deriv.fd como definición\nDESCARTADO: reproyecta\nla curva y no es idéntica", GRIS, lty = 2, cex = 0.62)
flecha(50, 64.5, 50, 59)
caja(50, 52, 60, 10, "Cruce de derivada positiva→no positiva:\ncandidato a máximo local (no fecha meteorológica)", VERDE, cex = 0.66)
flecha(50, 47, 50, 41.5)
caja(50, 34, 60, 10, "ELEGIDO: describir el patrón dominante\n(compatible con dos épocas lluviosas)", VERDE, font = 2, cex = 0.66)
text(50, 20, "PENDIENTE: una regla docente de inicio/fin de temporada. No se inventa un detector\nsin justificarlo; los cruces de la derivada no bastan como fechas oficiales.",
     cex = 0.62, col = "#333333")
leyenda(); cerrar()

# 7. Comparar las dos fuentes ------------------------------------------------
abrir("diagrama_decision_fuentes", 5.5)
titulo("Decisión 6 · ¿Cómo comparar estaciones y satélite?",
       "Sin pares válidos no se estima sesgo entre fuentes")
caja(50, 86, 50, 7, "Se tienen dos fuentes con coberturas distintas", AZUL, cex = 0.7)
flecha(38, 82.5, 20, 77); caja(19, 70, 32, 13,
  "Emparejar estación–píxel por cercanía\nDESCARTADO: faltan coordenadas,\naltitud y pares válidos", GRIS, lty = 2, cex = 0.62)
flecha(50, 82.5, 50, 77); caja(50, 70, 32, 13,
  "Comparar solo resúmenes agregados\nELEGIDO como descripción:\nformas y dispersión por separado", VERDE, cex = 0.62)
flecha(62, 82.5, 80, 77); caja(81, 70, 32, 13,
  "Modelo conjunto/regresión\nFUERA DE ALCANCE: es una\nfase posterior del proyecto", AMAR, lty = 3, cex = 0.62)
flecha(50, 63.5, 50, 58)
caja(50, 51, 64, 10, "Se evita afirmar sesgo, concordancia espacial o equivalencia de unidades\nporque no hay pares ni metadatos verificados", VERDE, cex = 0.66)
flecha(50, 46, 50, 40.5)
caja(50, 33, 64, 10, "ELEGIDO: interpretar cada fuente por separado\ny declarar la comparación como pendiente", VERDE, font = 2, cex = 0.66)
text(50, 19, "Péndulo opuesto: concluir que el satélite «mide mal» comparando promedios de lugares distintos.\nEso confunde composición del panel con sesgo de medición.",
     cex = 0.62, col = "#333333")
leyenda(); cerrar()

cat("Diagramas generados en:", outdir, "\n")

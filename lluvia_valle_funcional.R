# =============================================================================
# LLUVIA EN EL VALLE DEL CAUCA 1981-2003 — PARTE FUNCIONAL
# Proyecto final de Analítica de Datos (UAO, 2026-2S)
# =============================================================================
# Este archivo contiene TODA la parte funcional en R, en orden:
#   1. Datos preparados y emparejamiento estación-píxel
#   2. Construcción de datos funcionales            (rúbrica 1)
#   3. Análisis exploratorio de datos funcionales   (rúbrica 2)
#   4. Comparación funcional por fuente             (rúbrica 3)
#   5. Regresión funcional función a función        (rúbrica 4.2 y 4.4, referencia de la red)
#
# Entrada única: salidas_R_datos/datos_preparados.rds, que produce
#   lluvia_valle_datos.R (estaciones depuradas, satélite en mm/día, pares,
#   altitud y ONI). Aquí no se vuelve a leer ni a limpiar ningún dato crudo:
#   así las 26 pruebas del archivo de datos cubren todo lo que entra.
# Salidas: salidas_R_funcional/tablas, /graficos y /para_python (insumos de la red).
#
# Cómo se corre (desde la raíz del repositorio):
#   Rscript lluvia_valle_datos.R        # una vez, deja datos_preparados.rds
#   Rscript lluvia_valle_funcional.R    # este archivo, ~1 minuto
#
# Reglas del proyecto que se cumplen aquí:
#   - La base B-spline (Cox-de Boor), los mínimos cuadrados, S, W, la media, la
#     covarianza, las derivadas, las profundidades y el FPCA están programados
#     paso a paso. Las librerías fda y splines solo se usan para VERIFICAR.
#   - Toda elección (K, M, lambda, escalas) sale de entrenamiento (1981-1995) o
#     de validación (1996-2000). La prueba (2001-2003) solo se usa para evaluar.
#   - No se elimina ningún año (D-47, asesoría del 13-nov-2026).
#   - Cada decisión tiene su verificación numérica en el código: una tabla, una
#     gráfica o un stopifnot(). El bloque "QUÉ DICEN LOS DATOS" de cada sección
#     resume lo que salió.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
theme_set(theme_minimal(base_size = 11))

ANIOS     <- 1981:2003
PARTICION <- function(anio) fcase(anio <= 1995, "entrenamiento", anio <= 2000, "validación", default = "prueba")
RANGO     <- c(1, 365)
ORDEN     <- 4L
T_DIA     <- 1:365
# Años en que el satélite se queda más corto (La Niña). NO se excluyen (D-47): solo se usan
# para el análisis de sensibilidad "con y sin esos años".
ANIOS_MARCADOS <- c(1985L, 1988L, 2000L)

SALIDA_F <- "salidas_R_funcional"
dir.create(file.path(SALIDA_F, "tablas"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(SALIDA_F, "graficos"), showWarnings = FALSE, recursive = TRUE)
tabla_f   <- function(x, nombre) { fwrite(x, file.path(SALIDA_F, "tablas", paste0(nombre, ".csv"))); invisible(x) }
grafico_f <- function(p, nombre, ancho = 9, alto = 5) {
  ggsave(file.path(SALIDA_F, "graficos", paste0(nombre, ".png")), p, width = ancho, height = alto, dpi = 110, bg = "white")
  if (interactive()) print(p)
  invisible(p)
}
grafico_base <- function(nombre, dibujar, ancho = 11, alto = 5) {
  png(file.path(SALIDA_F, "graficos", paste0(nombre, ".png")), width = ancho, height = alto, units = "in", res = 110)
  dibujar(); dev.off()
  if (interactive()) dibujar()
  invisible(NULL)
}
titulo <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
fecha_dia <- function(d) format(as.Date(round(d) - 1, origin = "2001-01-01"), "%d-%b")

nodos <- function(K, orden = ORDEN) {
  interiores <- seq(RANGO[1], RANGO[2], length.out = K - orden + 2L)
  c(rep(RANGO[1], orden), interiores[-c(1L, length(interiores))], rep(RANGO[2], orden))
}
base_en_nodos <- function(t, tau, orden) {
  B <- outer(t, seq_len(length(tau) - 1L), function(x, j) (x >= tau[j]) & (x < tau[j + 1L])) * 1
  B[t == max(tau), max(which(diff(tau) > 0))] <- 1
  if (orden > 1L) for (m in 2:orden) {
    S <- matrix(0, length(t), length(tau) - m)
    for (j in seq_len(ncol(S))) {
      d1 <- tau[j + m - 1L] - tau[j]; d2 <- tau[j + m] - tau[j + 1L]
      if (d1 > 0) S[, j] <- S[, j] + (t - tau[j]) / d1 * B[, j]
      if (d2 > 0) S[, j] <- S[, j] + (tau[j + m] - t) / d2 * B[, j + 1L]
    }
    B <- S
  }
  B
}
base_bspline <- function(t, K, orden = ORDEN) base_en_nodos(t, nodos(K, orden), orden)
deriv_base <- function(t, K) {
  tau <- nodos(K); B <- base_en_nodos(t, tau, ORDEN - 1L); D <- matrix(0, length(t), K)
  for (j in seq_len(K)) {
    a <- tau[j + ORDEN - 1L] - tau[j]; b <- tau[j + ORDEN] - tau[j + 1L]
    if (a > 0) D[, j] <- D[, j] + (ORDEN - 1L) * B[, j] / a
    if (b > 0) D[, j] <- D[, j] - (ORDEN - 1L) * B[, j + 1L] / b
  }
  D
}

derivada_base <- function(t, K, n = 1L, tau = nodos(K), m = ORDEN) {
  if (n == 0L) return(base_en_nodos(t, tau, m))
  B <- derivada_base(t, K, n - 1L, tau, m - 1L); D <- matrix(0, length(t), length(tau) - m)
  for (j in seq_len(ncol(D))) {
    a <- tau[j + m - 1L] - tau[j]; b <- tau[j + m] - tau[j + 1L]
    if (a > 0) D[, j] <- D[, j] + (m - 1L) * B[, j] / a
    if (b > 0) D[, j] <- D[, j] - (m - 1L) * B[, j + 1L] / b
  }
  D
}
gram_derivada <- function(K, n = 2L) {
  tau <- unique(nodos(K)); a <- head(tau, -1); b <- tail(tau, -1)
  tq <- as.vector(outer(GL_x, seq_along(a), function(x, k) (b[k] - a[k]) / 2 * x + (a[k] + b[k]) / 2))
  wq <- as.vector(outer(GL_w, seq_along(a), function(w, k) (b[k] - a[k]) / 2 * w))
  P <- derivada_base(tq, K, n); crossprod(P, wq * P)
}

ajuste_mco <- function(y, Phi) {
  obs <- !is.na(y); P <- Phi[obs, , drop = FALSE]
  if (nrow(P) <= ncol(P) || qr(P)$rank != ncol(P)) stop("Curva no estimable con este K")
  mapa <- solve(crossprod(P), t(P))
  coef <- drop(mapa %*% y[obs]); H <- P %*% mapa
  res <- y[obs] - drop(P %*% coef)
  list(coef = coef, H = H, mapa = mapa, obs = obs, sigma2 = sum(res^2) / (nrow(P) - ncol(P)))
}
control_S <- function(H, K) c(simetria = max(abs(H - t(H))), idempotencia = max(abs(H %*% H - H)), traza_menos_K = abs(sum(diag(H)) - K))

criterios_K <- function(Y, t, Ks) {
  G <- crossprod(Y); n <- ncol(Y); N <- nrow(Y)
  rbindlist(lapply(Ks, function(K) {
    P <- base_bspline(t, K); mapa <- solve(crossprod(P), t(P)); h <- rowSums(P * t(mapa))
    T0 <- G %*% t(mapa)
    sse <- pmax(0, diag(G) - 2 * rowSums(P * T0) + rowSums((P %*% (mapa %*% T0)) * P))
    data.table(K = K, ECM = sum(sse) / (N * n), VC = sum(sse / (1 - h)^2) / (N * n), VCG = n * sum(sse) / (n - K)^2 / N)
  }))
}

construir <- function(Y, t, K) {
  Phi <- base_bspline(t, K); completas <- rowSums(is.na(Y)) == 0L
  C <- matrix(NA_real_, nrow(Y), K); sigma2 <- rep(NA_real_, nrow(Y))
  mapa <- solve(crossprod(Phi), t(Phi)); H <- Phi %*% mapa
  C[completas, ] <- Y[completas, , drop = FALSE] %*% t(mapa)
  Res <- Y[completas, , drop = FALSE] - C[completas, , drop = FALSE] %*% t(Phi)
  sigma2[completas] <- rowSums(Res^2) / (length(t) - K)
  controles <- control_S(H, K)
  for (i in which(!completas)) {
    a <- ajuste_mco(Y[i, ], Phi); C[i, ] <- a$coef; sigma2[i] <- a$sigma2
    controles <- pmax(controles, control_S(a$H, K))
  }
  colnames(C) <- paste0("c", seq_len(K))
  list(coef = C, K = K, Phi = Phi, mapa = mapa, H = H, sigma2 = sigma2, Res = Res, completas = completas,
       controles = controles, n_obs = rowSums(!is.na(Y)))
}

GL_x <- c(-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053)
GL_w <- c(0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454)
gram <- function(K) {
  tau <- unique(nodos(K)); a <- head(tau, -1); b <- tail(tau, -1)
  tq <- as.vector(outer(GL_x, seq_along(a), function(x, k) (b[k] - a[k]) / 2 * x + (a[k] + b[k]) / 2))
  wq <- as.vector(outer(GL_w, seq_along(a), function(w, k) (b[k] - a[k]) / 2 * w))
  P <- base_bspline(tq, K); crossprod(P, wq * P)
}

profundidades <- function(Y) {
  N <- nrow(Y); fm <- mbd <- numeric(N)
  for (j in seq_len(ncol(Y))) {
    r <- rank(Y[, j], ties.method = "max"); L <- rank(Y[, j], ties.method = "min") - 1; U <- N - r
    fm  <- fm + 1 - abs(.5 - r / N)
    mbd <- mbd + 1 - (L * (L - 1) + U * (U - 1)) / (N * (N - 1))
  }
  list(FM = fm / ncol(Y), MBD = mbd / ncol(Y))
}
boxplot_funcional <- function(Y, prof, factor = 1.5) {
  centro <- order(-prof, seq_along(prof))[seq_len(ceiling(nrow(Y) / 2))]
  lo <- apply(Y[centro, , drop = FALSE], 2, min); hi <- apply(Y[centro, , drop = FALSE], 2, max)
  inf <- lo - factor * (hi - lo); sup <- hi + factor * (hi - lo)
  atip <- apply(Y, 1, function(y) any(y < inf | y > sup))
  list(lo = lo, hi = hi, inf = inf, sup = sup, atip = atip, mediana = Y[which.max(prof), ],
       bigote_inf = apply(Y[!atip, , drop = FALSE], 2, min), bigote_sup = apply(Y[!atip, , drop = FALSE], 2, max))
}
dibujar_boxplot <- function(Y, bx, titulo, t = T_DIA) {
  plot(t, bx$mediana, type = "n", ylim = range(c(0, bx$bigote_inf, bx$bigote_sup, Y[bx$atip, ])), xlab = "Día del año", ylab = "mm/día", main = titulo)
  if (any(bx$atip)) matlines(t, t(Y[bx$atip, , drop = FALSE]), col = "#cc333330", lty = 1)
  polygon(c(t, rev(t)), c(bx$lo, rev(bx$hi)), col = "#b39ddbcc", border = NA)
  lines(t, bx$bigote_inf, col = 4); lines(t, bx$bigote_sup, col = 4)
  lines(t, bx$mediana, lwd = 2)
  legend("topleft", c("Región central 50 %", "Mediana (MBD)", "Bigotes", "Atípicas"),
         col = c("#b39ddb", 1, 4, "#cc3333"), lty = c(NA, 1, 1, 1), pch = c(15, NA, NA, NA), cex = .75, bty = "n")
}
mapa_calor <- function(t, z, titulo, limites = range(z, na.rm = TRUE), paleta = "Blue-Red 3", xlab = "Día s", ylab = "Día t") {
  colores <- hcl.colors(64, paleta)
  layout(matrix(1:2, 1), widths = c(5, 1)); par(mar = c(4, 4, 3, 1))
  image(t, t, z, zlim = limites, col = colores, main = titulo, xlab = xlab, ylab = ylab)
  par(mar = c(4, 0, 3, 4)); escala <- seq(limites[1], limites[2], length.out = 65)
  image(c(0, 1), escala, matrix(rep(escala, each = 2), nrow = 2), col = colores, zlim = limites, axes = FALSE, xlab = "", ylab = "")
  axis(4, las = 1, cex.axis = .8); box(); layout(1)
}

# =============================================================================
# 1. DATOS PREPARADOS Y EMPAREJAMIENTO ESTACIÓN-PÍXEL  — rúbrica 3 (emparejamiento)
# =============================================================================
# Lee datos_preparados.rds y deja las matrices con las que trabaja todo el
# archivo: Y (estación-año x 365 días, mm/día, NA = día faltante) y X (píxel-año
# x 72 péntadas, mm/día, cada péntada en su día central). Documenta el
# emparejamiento (regla, distancias y píxeles compartidos) y verifica con datos
# dos decisiones: la vertiente de los píxeles y el trato de 1985, 1988 y 2000.
# =============================================================================
titulo("1. Datos preparados y emparejamiento estación-píxel")

dp  <- readRDS("salidas_R_datos/datos_preparados.rds")
oni <- fread("datos_web/oni_anual_fase.csv")
Y <- dp$estaciones$Y; X <- dp$satelite$X; T_SAT <- dp$satelite$t
storage.mode(Y) <- "double"
Y_dia <- Y; X_pen <- X      # copias de los datos discretos (las secciones 4 y 5 reutilizan los nombres Y y X para curvas)
stopifnot(ncol(Y) == 365, ncol(X) == 72, !anyNA(X), all(Y >= 0, na.rm = TRUE), all(X >= 0),
          identical(dp$estaciones$id$fila, seq_len(nrow(Y))), identical(dp$satelite$id$fila, seq_len(nrow(X))),
          all(PARTICION(dp$estaciones$id$anio) == dp$estaciones$id$particion))

est <- copy(dp$meta_estaciones)
id_est <- dp$estaciones$id[, .(codigo, anio, particion, n_obs, hueco_max)]
id_est[, curva_id := paste(codigo, anio, sep = "_")]

# --- Emparejamiento: cada estación toma la celda de la malla que la contiene ---
px <- dp$pixeles[, .(celda = pixel, lon, lat, en_valle, con_estacion)]
emparejamiento <- px[, .(celda, lon_px = lon, lat_px = lat, px_en_valle = en_valle)][
  est[, .(codigo, nombre, vertiente, altitud, lat, lon, celda = pixel, dist_km = km_al_centro)], on = "celda"]
emparejamiento[, estaciones_en_el_pixel := .N, by = celda]
stopifnot(!anyNA(emparejamiento$lon_px), abs(emparejamiento$lon - emparejamiento$lon_px) <= 0.025 + 1e-6,
          abs(emparejamiento$lat - emparejamiento$lat_px) <= 0.025 + 1e-6)   # la estación cae DENTRO de su celda de 0,05°
tabla_f(emparejamiento[order(celda, codigo)], "F1_emparejamiento_estacion_pixel")
compartidos <- emparejamiento[estaciones_en_el_pixel > 1, .(estaciones = .N, codigos = paste(codigo, collapse = " / "),
                                                            nombres = paste(nombre, collapse = " / ")), by = celda]
tabla_f(compartidos, "F1_pixeles_compartidos")
cat("Estaciones:", nrow(emparejamiento), "| píxeles distintos:", uniqueN(emparejamiento$celda),
    "| píxeles con más de una estación:", nrow(compartidos), "(", sum(compartidos$estaciones), "estaciones ) |",
    "distancia al centro del píxel: mediana", round(median(emparejamiento$dist_km), 2), "km, máx", round(max(emparejamiento$dist_km), 2), "km |",
    "píxel con centro fuera del polígono del Valle:", emparejamiento[px_en_valle == FALSE, .N], "\n")
print(compartidos[order(-estaciones)][1:5])

# --- Vertiente de cada píxel: la de la estación más cercana (regla verificada) ---
# La vertiente de las estaciones viene del código CVC (cuenca hidrográfica). Los
# píxeles no traen ese dato. Se comparan dos reglas dejando una estación fuera
# cada vez y prediciendo su vertiente: (a) un corte en la longitud; (b) la
# vertiente de la estación más cercana.
haversine <- function(la1, lo1, la2, lo2) 12742 * asin(sqrt(sin((la2 - la1) * pi / 360)^2 + cos(la1 * pi / 180) * cos(la2 * pi / 180) * sin((lo2 - lo1) * pi / 360)^2))
D_ee <- outer(seq_len(nrow(est)), seq_len(nrow(est)), function(i, j) haversine(est$lat[i], est$lon[i], est$lat[j], est$lon[j])); diag(D_ee) <- Inf
cortes <- seq(-77.1, -76.4, by = .05)
acierto_corte <- sapply(cortes, function(k) mean((est$lon < k) == (est$vertiente == "Pacífico")))
acierto_vecina <- mean(est$vertiente[apply(D_ee, 1, which.min)] == est$vertiente)
regla_vertiente <- data.table(regla = c(paste0("corte en longitud ", cortes[which.max(acierto_corte)], "° (el mejor corte)"), "corte en longitud -76.8° (el que se usaba)", "estación más cercana (dejando una fuera)"),
                              pct_acierto = round(100 * c(max(acierto_corte), acierto_corte[abs(cortes + 76.8) < 1e-9], acierto_vecina), 1))
print(tabla_f(regla_vertiente, "F1_regla_vertiente_pixeles"))
D_pe <- outer(seq_len(nrow(px)), seq_len(nrow(est)), function(i, j) haversine(px$lat[i], px$lon[i], est$lat[j], est$lon[j]))
px[, `:=`(vertiente = est$vertiente[apply(D_pe, 1, which.min)], km_estacion_cercana = round(apply(D_pe, 1, min), 1))]
stopifnot(acierto_vecina > max(acierto_corte))
cat("Píxeles por vertiente (estación más cercana):\n"); print(px[, .(pixeles = .N, km_mediana_a_la_estacion = median(km_estacion_cercana)), by = vertiente])

id_sat <- px[dp$satelite$id[, .(celda = pixel, anio, particion)], on = "celda"]
stopifnot(identical(id_sat$celda, dp$satelite$id$pixel), dp$satelite$dias[12] == 3)
dias_reales <- matrix(dp$satelite$dias, nrow(X), 72, byrow = TRUE); dias_reales[id_sat$anio %% 4 == 0, 12] <- 4
id_sat[, `:=`(n_tope = rowSums(X * dias_reales >= 199.99), marcado = anio %in% ANIOS_MARCADOS, curva_id = paste(anio, celda, sep = "_"))]
cat("Curvas de estación:", nrow(Y), "(", uniqueN(id_est$codigo), "estaciones ) | curvas de satélite:", nrow(X), "=", uniqueN(id_sat$anio), "años x", uniqueN(id_sat$celda), "píxeles |",
    "con alguna péntada en el tope de 200 mm:", sum(id_sat$n_tope > 0), sprintf("(%.1f %%)", 100 * mean(id_sat$n_tope > 0)), "\n")
print(id_est[, .(curvas = .N, completas = sum(n_obs == 365), estaciones = uniqueN(codigo)), by = particion])

# --- 1985, 1988 y 2000: ¿archivos dañados o error que sigue al ENSO? ---
# Se compara, año por año, la lluvia anual de cada píxel con la de su estación
# (datos discretos, sin suavizar) y se cruza con el ONI.
anual <- dp$pares[, .(codigo, anio, particion, ONI, fase = fase_enso, P_est = rowMeans(Y[fila_est, ], na.rm = TRUE) * 365, P_sat = rowMeans(X[fila_sat, ]) * 365)]
por_anio <- anual[, .(pares = .N, P_est_mediana = round(median(P_est)), P_sat_mediana = round(median(P_sat)), cociente_mediano = round(median(P_sat / P_est), 3),
                      ONI = ONI[1], fase = fase[1], particion = particion[1]), by = anio][order(anio)]
por_anio[, marcado := anio %in% ANIOS_MARCADOS]
print(tabla_f(por_anio, "F1_satelite_vs_estaciones_por_anio"))
r_oni <- por_anio[, cor(ONI, cociente_mediano)]; r_oni_ent <- por_anio[particion == "entrenamiento", cor(ONI, cociente_mediano)]
cat(sprintf("Correlación entre el ONI del año y el cociente satélite/estación: r = %.2f con los 23 años (%.2f solo con entrenamiento)\n", r_oni, r_oni_ent))
cat("Cociente mediano por fase ENSO:\n"); print(por_anio[, .(anios = .N, cociente = round(mean(cociente_mediano), 3), anios_lista = paste(anio, collapse = " ")), by = fase])
cat("Años con el cociente más bajo:", paste(por_anio[order(cociente_mediano)][1:5, sprintf("%d (%.2f, %s)", anio, cociente_mediano, fase)], collapse = ", "), "\n")
grafico_f(ggplot(por_anio, aes(ONI, cociente_mediano)) + geom_hline(yintercept = 1, linetype = 3) +
            geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey50", linewidth = .5) +
            geom_text(aes(colour = fase, label = anio), size = 3.2, fontface = ifelse(por_anio$marcado, "bold", "plain")) +
            scale_colour_manual(values = c(Neutro = "grey30", `Niña` = "steelblue", `Niño` = "firebrick")) +
            labs(x = "ONI anual", y = "Mediana de satélite / estación (lluvia anual)", colour = NULL,
                 title = "El error del satélite sigue al ENSO", subtitle = sprintf("r = %.2f; en negrilla 1985, 1988 y 2000", r_oni)),
          "F1_cociente_satelite_estacion_vs_ONI", ancho = 8, alto = 5)

# QUÉ DICEN LOS DATOS
#   Datos: 3 299 curvas estación-año de 176 estaciones (2 257 de entrenamiento,
#   672 de validación y 370 de prueba; 2 812 sin ningún día faltante) y 15 663
#   curvas píxel-año (681 píxeles x 23 años). Todo viene de datos_preparados.rds.
#   Emparejamiento: cada estación toma la celda de 0,05° que la contiene (se
#   comprueba que la coordenada cae dentro). Las 176 estaciones caen en 148
#   píxeles; 23 píxeles tienen 2 a 4 estaciones (51 estaciones), que comparten la
#   misma curva satelital. El centro del píxel queda a 2,2 km en mediana y nunca
#   a más de 3,9 km. Cinco píxeles tienen el centro fuera del polígono del Valle,
#   pero su estación está dentro.
#   Vertiente de los píxeles: antes se usaba un corte en la longitud -76,8°. Al
#   comprobarlo con las estaciones, ese corte solo acierta el 88,6 % y deja como
#   "Cauca" a 20 de las 26 estaciones que la CVC pone en la vertiente Pacífica
#   (cuencas del Dagua, Calima y Garrapatas, sobre la cordillera Occidental).
#   Ningún corte en longitud pasa de 89,2 %. La regla "vertiente de la estación
#   más cercana" acierta el 97,7 % dejando una estación fuera cada vez, así que
#   se usa esa: una sola definición (la cuenca de la CVC) para las dos fuentes.
#   1985, 1988 y 2000: no son archivos dañados. El cociente satélite/estación de
#   cada año sigue al ONI (r = 0,68): en años Niña el satélite se queda corto
#   (cociente medio 0,99) y en años Niño se pasa (1,29); en neutros da 1,16. Los
#   cinco cocientes más bajos son de 1988, 1984, 2000, 1999 y 1996: cuatro Niña.
#   1984 y 1999 se apartan tanto como los tres años marcados y nunca se marcaron.
#   Es un error sistemático del satélite ligado al ENSO, justo lo que el modelo
#   debe aprender con el ONI. Por eso no se excluye ningún año (D-47).


# =============================================================================
# 2. CONSTRUCCIÓN DE DATOS FUNCIONALES (ESTACIONES Y SATÉLITE)  — rúbrica 1
# =============================================================================
# Base B-spline cúbica (Cox-de Boor) y mínimos cuadrados programados sin
# librerías; K por ECM, VC y VCG solo con entrenamiento (1981-1995, todos los
# años); verificación de S; covarianza residual; bandas punto a punto; matriz
# de Gram W. Ambas fuentes quedan en el mismo dominio 1-365.
# Días faltantes: cada curva se ajusta solo con sus días observados (no se imputa).
titulo("2. Construcción de datos funcionales")

completas_est <- rowSums(is.na(Y)) == 0
stopifnot(all(completas_est == (id_est$n_obs == 365L)))

grafico_base("F2_datos_discretos", function() {
  par(mfrow = c(1, 2))
  i <- which(completas_est & id_est$codigo == "2622330101")[1]
  plot(T_DIA, Y[i, ], type = "h", col = "grey40", xlab = "Día del año", ylab = "mm/día",
       main = paste("Estación", id_est$codigo[i], id_est$anio[i], "(diario)"))
  j <- which(id_sat$celda == emparejamiento[codigo == "2622330101", celda] & id_sat$anio == id_est$anio[i])
  plot(T_SAT, X[j, ], type = "b", pch = 16, cex = .6, col = "darkorange", xlab = "Día del año", ylab = "mm/día",
       main = paste("Su píxel, mismo año (péntada en su día central)"))
})

for (K in c(6, 8, 11, 20, 40)) {
  P <- base_bspline(T_DIA, K)
  stopifnot(all(is.finite(P)), min(P) >= -1e-12, max(abs(rowSums(P) - 1)) < 1e-10)
}
ref_splines <- splines::splineDesign(nodos(11), T_DIA, ord = ORDEN, outer.ok = TRUE)
cat("Base propia vs splines::splineDesign (K = 11), diferencia máxima:", format(max(abs(base_bspline(T_DIA, 11) - ref_splines)), digits = 2), "\n")

i_ej <- which(completas_est)[1]; y1 <- Y[i_ej, ]
P8 <- base_bspline(T_DIA, 8); a8 <- ajuste_mco(y1, P8)
vc_atajo <- mean(((y1 - a8$H %*% y1) / (1 - diag(a8$H)))^2)
vc_bruta <- mean(sapply(1:365, function(j) (y1[j] - P8[j, ] %*% solve(crossprod(P8[-j, ]), crossprod(P8[-j, ], y1[-j])))^2))
cat(sprintf("VC por el atajo de diag(S) = %.8f | VC dejando un día fuera (365 ajustes) = %.8f\n", vc_atajo, vc_bruta))

ent_est <- completas_est & id_est$particion == "entrenamiento"
ent_sat <- id_sat$particion == "entrenamiento"
sel_est <- criterios_K(Y[ent_est, ], T_DIA, c(4:40, seq(45, 100, 5)))
sel_sat <- criterios_K(X[ent_sat, ], T_SAT, 4:40)
sel <- rbind(sel_est[, fuente := "estaciones"], sel_sat[, fuente := "satélite"])
tabla_f(sel, "F2_seleccion_K")
resumen_K <- sel[, .(K_VC = K[which.min(VC)], K_VCG = K[which.min(VCG)],
                     K_VC_1pct = paste(range(K[VC <= 1.01 * min(VC)]), collapse = "-")), by = fuente]
resumen_K[, curvas_entrenamiento := c(sum(ent_est), sum(ent_sat))]
print(tabla_f(resumen_K, "F2_K_elegido"))
K_EST <- resumen_K[fuente == "estaciones", K_VCG]; K_SAT <- resumen_K[fuente == "satélite", K_VCG]

# Sensibilidad de K en las estaciones. La VC es casi plana alrededor del mínimo, así que
# se comprueba qué cambia con una base más flexible: las fechas de los picos (métrica de
# la rúbrica 4.4) y los valores negativos de la curva. Solo con curvas de entrenamiento.
sens_K <- rbindlist(lapply(sort(unique(c(K_EST, 12L, 16L, 24L))), function(K) {
  P <- base_bspline(T_DIA, K); Fk <- Y[ent_est, ] %*% t(solve(crossprod(P), t(P))) %*% t(P)
  p1 <- 59 + apply(Fk[, 60:200], 1, which.max); p2 <- 229 + apply(Fk[, 230:350], 1, which.max)
  data.table(K = K, VC_sobre_el_minimo_pct = round(100 * (sel_est$VC[sel_est$K == K] / min(sel_est$VC) - 1), 2),
             pico_1_mediana = fecha_dia(median(p1)), pico_1_RIC_dias = IQR(p1), pico_2_mediana = fecha_dia(median(p2)), pico_2_RIC_dias = IQR(p2),
             pct_picos_en_el_borde_de_la_ventana = round(100 * mean(p1 %in% c(60, 200) | p2 %in% c(230, 350)), 1),
             pct_puntos_negativos = round(100 * mean(Fk < 0), 2))
}))
print(tabla_f(sens_K, "F2_sensibilidad_K_estaciones"))

sel_largo <- melt(sel, id.vars = c("fuente", "K"), variable.name = "criterio")
marcas_K <- sel_largo[criterio != "ECM", .SD[which.min(value)], by = .(fuente, criterio)]
grafico_f(ggplot(sel_largo[K <= 60], aes(K, value)) + geom_line(colour = "steelblue") + geom_point(size = .8, colour = "steelblue") +
            geom_vline(data = marcas_K, aes(xintercept = K), colour = "red", linetype = 2) +
            geom_text(data = marcas_K, aes(K, Inf, label = paste("K =", K)), vjust = 1.5, hjust = -.1, colour = "red", size = 3) +
            facet_wrap(fuente ~ criterio, scales = "free_y", ncol = 3) + scale_y_log10() +
            labs(x = "Número de bases K", y = "Criterio (log)", title = "Selección de K con curvas de entrenamiento"),
          "F2_seleccion_K", ancho = 11, alto = 6)

e <- construir(Y, T_DIA, K_EST)
s <- construir(X, T_SAT, K_SAT)
controles <- rbind(data.table(fuente = "estaciones", K = K_EST, t(e$controles)), data.table(fuente = "satélite", K = K_SAT, t(s$controles)))
print(tabla_f(controles, "F2_controles_S"))
stopifnot(all(controles[, .(simetria, idempotencia, traza_menos_K)] < 1e-8))

Phi_e <- base_bspline(T_DIA, K_EST); Phi_s365 <- base_bspline(T_DIA, K_SAT)
F_est <- e$coef %*% t(Phi_e); F_sat <- s$coef %*% t(Phi_s365)
negativos <- data.table(fuente = c("estaciones", "satélite"),
                        pct_puntos_negativos = round(100 * c(mean(F_est < 0), mean(F_sat < 0)), 2),
                        pct_curvas_con_negativos = round(100 * c(mean(apply(F_est < 0, 1, any)), mean(apply(F_sat < 0, 1, any))), 1),
                        minimo_mm_dia = round(c(min(F_est), min(F_sat)), 2))
print(tabla_f(negativos, "F2_valores_negativos"))

grafico_base("F2_bases", function() {
  par(mfrow = c(1, 2))
  for (K in c(K_EST, K_SAT)) {
    matplot(T_DIA, base_bspline(T_DIA, K), type = "l", lty = 1, xlab = "Día del año", ylab = expression(phi[k](t)),
            main = paste("Base B-spline cúbica, K =", K, if (K == K_EST) "(estaciones)" else "(satélite)"))
    abline(v = unique(nodos(K)), lty = 3, col = "grey60")
  }
})

Sigma_res_est <- crossprod(e$Res) / nrow(e$Res); Sigma_res_sat <- crossprod(s$Res) / nrow(s$Res)
Cor_res_est <- cov2cor(Sigma_res_est); Cor_res_sat <- cov2cor(Sigma_res_sat)
v_t <- diag(Sigma_res_est) / (1 - diag(e$H)); forma <- v_t / mean(v_t)
cat(sprintf("|Correlación residual| fuera de la diagonal: estaciones %.3f (rezago de 1 día %.3f) | satélite %.3f (rezago de 1 péntada %.3f)\n",
            mean(abs(Cor_res_est[upper.tri(Cor_res_est)])), mean(diag(Cor_res_est[-1, -365])),
            mean(abs(Cor_res_sat[upper.tri(Cor_res_sat)])), mean(diag(Cor_res_sat[-1, -72]))))
cat(sprintf("Varianza del error por día (estaciones): mínima %.1f el %s, máxima %.1f el %s\n",
            min(v_t), fecha_dia(which.min(v_t)), max(v_t), fecha_dia(which.max(v_t))))
tabla_f(data.table(dia = T_DIA, varianza_error = v_t), "F2_varianza_error_estaciones")

grafico_base("F2_covarianza_residual", function() {
  par(mfrow = c(1, 3))
  image(T_DIA, T_DIA, Cor_res_est, zlim = c(-.3, .3), col = hcl.colors(30, "RdBu", rev = TRUE),
        main = "Correlación residual: estaciones", xlab = "Día s", ylab = "Día t")
  image(T_SAT, T_SAT, Cor_res_sat, zlim = c(-1, 1), col = hcl.colors(30, "RdBu", rev = TRUE),
        main = "Correlación residual: satélite", xlab = "Día s", ylab = "Día t")
  plot(T_DIA, v_t, type = "l", xlab = "Día del año", ylab = "(mm/día)^2", main = "Varianza del error por día (estaciones)")
}, ancho = 14, alto = 4.8)

ee_banda <- function(Phi_eval, mapa, Sigma) { M <- Phi_eval %*% mapa; sqrt(rowSums((M %*% Sigma) * M)) }
ee_iid <- ee_banda(Phi_e, e$mapa, e$sigma2[i_ej] * diag(365))
ee_het <- ee_banda(Phi_e, e$mapa, diag(e$sigma2[i_ej] * forma))
racha <- function(m) { r <- rle(m); max(c(0, r$lengths[r$values])) }
inc <- which(!e$completas); ig <- inc[which.max(sapply(inc, function(k) racha(is.na(Y[k, ]))))]
ag <- ajuste_mco(Y[ig, ], Phi_e)
ee_g <- ee_banda(Phi_e, ag$mapa, ag$sigma2 * diag(sum(ag$obs)))
ee_g_completa <- ee_banda(Phi_e, e$mapa, ag$sigma2 * diag(365))
js <- which(id_sat$celda == emparejamiento[codigo == id_est$codigo[i_ej], celda] & id_sat$anio == id_est$anio[i_ej])
ee_s <- ee_banda(Phi_s365, s$mapa, s$sigma2[js] * diag(72))
cat(sprintf("Curva %s con hueco de %d días: el EE en los días faltantes es %.2f veces el de la misma curva sin hueco\n",
            id_est$curva_id[ig], sum(!ag$obs), mean((ee_g / ee_g_completa)[!ag$obs])))
cat(sprintf("EE medio de la banda: estación %.2f mm/día (sigma %.1f) | su píxel %.2f mm/día (sigma %.1f)\n",
            mean(ee_iid), sqrt(e$sigma2[i_ej]), mean(ee_s), sqrt(s$sigma2[js])))
tabla_f(data.table(curva_id = id_est$curva_id[i_ej], dia = T_DIA, curva = F_est[i_ej, ], ee_iid, ee_het,
                   inferior = F_est[i_ej, ] - 2 * ee_iid, superior = F_est[i_ej, ] + 2 * ee_iid), "F2_banda_ejemplo_estacion")

grafico_base("F2_bandas_punto_a_punto", function() {
  par(mfrow = c(1, 3))
  yl <- range(c(F_est[i_ej, ] - 2 * ee_het, Y[i_ej, ]), na.rm = TRUE)
  plot(T_DIA, Y[i_ej, ], pch = 16, cex = .4, col = "grey50", ylim = yl, xlab = "Día del año", ylab = "mm/día",
       main = paste("Estación", id_est$curva_id[i_ej]))
  lines(T_DIA, F_est[i_ej, ], lwd = 2, col = 4)
  lines(T_DIA, F_est[i_ej, ] + 2 * ee_iid, lty = 2, col = 4); lines(T_DIA, F_est[i_ej, ] - 2 * ee_iid, lty = 2, col = 4)
  lines(T_DIA, F_est[i_ej, ] + 2 * ee_het, lty = 3); lines(T_DIA, F_est[i_ej, ] - 2 * ee_het, lty = 3)
  legend("topright", c("Curva", "±2 EE, sigma^2 I", "±2 EE, sigma^2(t)"), lty = 1:3, col = c(4, 4, 1), bty = "n", cex = .8)
  cs <- F_sat[js, ]
  plot(T_SAT, X[js, ], pch = 16, cex = .7, col = "grey50", ylim = range(c(X[js, ], cs - 2 * ee_s, cs + 2 * ee_s)),
       xlab = "Día del año", ylab = "mm/día", main = paste("Su píxel, año", id_sat$anio[js]))
  lines(T_DIA, cs, lwd = 2, col = "darkorange"); lines(T_DIA, cs + 2 * ee_s, lty = 2); lines(T_DIA, cs - 2 * ee_s, lty = 2)
  plot(T_DIA, Y[ig, ], pch = 16, cex = .4, col = "grey50", ylim = c(-5, max(Y[ig, ], na.rm = TRUE)), xlab = "Día del año",
       ylab = "mm/día", main = paste("Estación", id_est$curva_id[ig], "-", sum(!ag$obs), "días faltantes"))
  lines(T_DIA, F_est[ig, ], lwd = 2, col = 4)
  lines(T_DIA, F_est[ig, ] + 2 * ee_g, lty = 2, col = 4); lines(T_DIA, F_est[ig, ] - 2 * ee_g, lty = 2, col = 4)
  rug(T_DIA[!ag$obs], col = 2, lwd = 2)
}, ancho = 15, alto = 4.8)

W_est <- gram(K_EST); W_sat <- gram(K_SAT)
for (w in list(list("estaciones", W_est), list("satélite", W_sat))) {
  W <- w[[2]]
  stopifnot(max(abs(W - t(W))) < 1e-10, min(eigen(W, symmetric = TRUE)$values) > 0, abs(sum(W) - 364) < 1e-8)
  cat(sprintf("W %s (%d x %d): simétrica, definida positiva (autovalor mínimo %.3f), 1'W1 = %.8f = longitud del dominio 364\n",
              w[[1]], nrow(W), ncol(W), min(eigen(W, symmetric = TRUE)$values), sum(W)))
}
t_fino <- seq(1, 365, by = .01); x_fino <- drop(base_bspline(t_fino, K_EST) %*% e$coef[i_ej, ])
cat(sprintf("Integral de x(t)^2: c'Wc = %.4f | trapecio con paso 0,01 = %.4f\n",
            drop(e$coef[i_ej, ] %*% W_est %*% e$coef[i_ej, ]), sum(diff(t_fino) * (head(x_fino, -1)^2 + tail(x_fino, -1)^2) / 2)))

grafico_base("F2_mismo_dominio", function() {
  par(mfrow = c(1, 1))
  plot(T_DIA, F_est[i_ej, ], type = "l", lwd = 2, col = 4, ylim = range(c(F_est[i_ej, ], F_sat[js, ], 0)),
       xlab = "Día del año (1-365)", ylab = "mm/día", main = paste("Estación", id_est$curva_id[i_ej], "y su píxel en el mismo dominio"))
  lines(T_DIA, F_sat[js, ], lwd = 2, col = "darkorange")
  points(T_SAT, X[js, ], pch = 16, cex = .6, col = "darkorange")
  legend("topright", c(paste("Estación, K =", K_EST), paste("Satélite, K =", K_SAT)), col = c(4, "darkorange"), lwd = 2, bty = "n")
}, ancho = 9, alto = 4.5)

# QUÉ DICEN LOS DATOS
#   La base propia coincide con splines::splineDesign (3e-16) y el atajo de la VC
#   con diag(S) da lo mismo que dejar cada día fuera (365 ajustes).
#   K se elige solo con entrenamiento (1 882 curvas completas de estación y
#   10 215 de satélite, con todos los años): VC y VCG coinciden, K = 8 en
#   estaciones y K = 11 en satélite. El ECM de ajuste baja siempre con K, por eso
#   no sirve para elegir. En el satélite el mínimo es nítido. En estaciones la VC
#   es casi plana entre K = 5 y 17 (a menos de 1 % del mínimo): la lluvia diaria
#   es tan ruidosa que más bases solo persiguen aguaceros.
#   Sensibilidad de K (estaciones): con K = 12, 16 o 24 la fecha mediana de los
#   picos se mueve a lo sumo 9 días (18 a 27 de abril; 24 a 29 de octubre), mucho
#   menos que la dispersión entre curvas (29 a 45 días), y los valores negativos
#   suben de 2,0 % a 5,4 %. Una base más flexible no cambia las conclusiones y
#   empeora la VC, así que se mantiene K = 8.
#   S es simétrica, idempotente y tr(S) = K hasta 1e-14 en el peor de los ajustes
#   (también en las curvas con días faltantes, que usan su propia S).
#   El error de las estaciones no es homogéneo: su varianza va de 26 (mm/día)^2 en
#   agosto a 264 en octubre, la temporada húmeda. La correlación residual entre
#   días es casi nula (0,04; 0,08 a un día), así que la banda con sigma^2 I es
#   razonable. En el satélite es mayor (0,14): la péntada ya promedia cinco días.
#   Las bandas son punto a punto: cada día tiene 95 % de cobertura por separado,
#   no la curva completa. Un hueco de 15 días ensancha la banda solo 15 % en esos
#   días.
#   Valores negativos: 47 % de las curvas de estación baja de 0 en algún día (2 %
#   de los puntos), en la temporada seca de julio-agosto y en los bordes 1-ene y
#   31-dic (efecto de borde de la base no periódica). En el satélite casi no
#   ocurre (1 % de las curvas). Para modelar, la curva se recorta en 0 (D-23).
#   Por qué B-splines y no Fourier: una curva de un solo año no es periódica (el
#   31-dic no tiene por qué empalmar con el 1-ene del mismo año), el ejemplo de
#   clase con esta lluvia usa B-splines cúbicas y el enunciado las exige en la
#   red. El costo es el efecto de borde, que queda documentado.
#   W integra exactamente: 1'W1 = 364 y c'Wc coincide con el trapecio fino.

o <- list(id_est = id_est, id_sat = id_sat, coef_est = e$coef, coef_sat = s$coef, K_est = K_EST, K_sat = K_SAT,
             sigma2_est = e$sigma2, sigma2_sat = s$sigma2, n_obs_est = e$n_obs, W_est = W_est, W_sat = W_sat,
             t_sat = T_SAT, seleccion_K = sel, controles_S = controles, Sigma_res_est = Sigma_res_est, Sigma_res_sat = Sigma_res_sat,
             emparejamiento = emparejamiento, estaciones = est)

# =============================================================================
# 3. ANÁLISIS EXPLORATORIO DE DATOS FUNCIONALES (CADA FUENTE)  — rúbrica 2
# =============================================================================
# Media paso a paso (verificada con fda), varianza, desviación estándar,
# covarianza y correlación, derivadas y temporadas, profundidad FM y MBD,
# boxplot funcional con atípicas por estación y año, y contexto climatológico
# (vertiente, altitud y ENSO). Se usan todas las estaciones incluidas y todos
# los años 1981-2003; el satélite se analiza en los píxeles del Valle.
titulo("3. EDA funcional")

Phi_e <- base_bspline(T_DIA, o$K_est); Phi_s <- base_bspline(T_DIA, o$K_sat)
valle_s <- o$id_sat$en_valle
C_e <- o$coef_est; C_s <- o$coef_sat[valle_s, ]
id_e <- est[, .(codigo, nombre, vertiente, altitud)][o$id_est, on = "codigo"]
id_s <- o$id_sat[valle_s]
Ye <- C_e %*% t(Phi_e); Xs <- C_s %*% t(Phi_s)
cat("Estaciones:", nrow(Ye), "curvas,", uniqueN(id_e$codigo), "estaciones, años", paste(range(id_e$anio), collapse = "-"),
    "| satélite (Valle):", nrow(Xs), "curvas,", uniqueN(id_s$celda), "píxeles\n")
print(id_e[, .N, by = vertiente])

titulo("3.1 Media funcional")
media_e <- drop(Phi_e %*% colMeans(C_e)); media_s <- drop(Phi_s %*% colMeans(C_s))
fd_e <- fda::fd(t(C_e), fda::create.bspline.basis(RANGO, nbasis = o$K_est, norder = ORDEN))
fd_s <- fda::fd(t(C_s), fda::create.bspline.basis(RANGO, nbasis = o$K_sat, norder = ORDEN))
verif_media <- data.table(fuente = c("estaciones", "satélite"),
                          vs_curvas_evaluadas = c(max(abs(media_e - colMeans(Ye))), max(abs(media_s - colMeans(Xs)))),
                          vs_fda_mean.fd = c(max(abs(media_e - fda::eval.fd(T_DIA, fda::mean.fd(fd_e)))),
                                             max(abs(media_s - fda::eval.fd(T_DIA, fda::mean.fd(fd_s))))))
print(tabla_f(verif_media, "F3_verificacion_media"))
stopifnot(verif_media$vs_fda_mean.fd < 1e-10)
cat(sprintf("Lluvia media: estaciones %.0f mm/año | satélite en el Valle %.0f mm/año\n", mean(media_e) * 365, mean(media_s) * 365))

titulo("3.2 Varianza, desviación estándar, covarianza y correlación")
momentos <- function(C, Phi) {
  Cc <- sweep(C, 2, colMeans(C)); Sc <- crossprod(Cc) / (nrow(C) - 1)
  covar <- Phi %*% Sc %*% t(Phi)
  list(var = diag(covar), de = sqrt(diag(covar)), cov = covar, cor = cov2cor(covar))
}
me <- momentos(C_e, Phi_e); ms <- momentos(C_s, Phi_s)
cat("Covarianza Phi Sigma_c Phi' vs cov() de las curvas evaluadas (estaciones):", format(max(abs(me$cov - cov(Ye))), digits = 2), "\n")
resumen_disp <- data.table(fuente = c("estaciones", "satélite"),
                           dia_DE_max = c(which.max(me$de), which.max(ms$de)), DE_max = c(max(me$de), max(ms$de)),
                           dia_DE_min = c(which.min(me$de), which.min(ms$de)), DE_min = c(min(me$de), min(ms$de)),
                           CV_medio = c(mean(me$de / media_e), mean(ms$de / media_s)),
                           cor_abr_oct = c(me$cor[105, 290], ms$cor[105, 290]), cor_abr_jul = c(me$cor[105, 200], ms$cor[105, 200]),
                           cor_ene_dic = c(me$cor[15, 350], ms$cor[15, 350]))
resumen_disp[, `:=`(fecha_DE_max = fecha_dia(dia_DE_max), fecha_DE_min = fecha_dia(dia_DE_min))]
print(tabla_f(resumen_disp, "F3_dispersion_resumen"))

# ¿Qué mide la covarianza con todas las estaciones juntas? Se separa por vertiente y se
# repite con ANOMALÍAS (cada curva menos la curva media de su estación o de su píxel),
# que quitan las diferencias entre sitios y dejan la variación de un año a otro.
anom_coef <- function(C, g) { g <- as.character(g); m <- rowsum(C, g) / as.vector(table(g)[rownames(rowsum(C, g))]); C - m[g, ] }
espectro <- function(C, Phi, nombre) {
  W <- gram(ncol(C)); ev <- eigen(W, symmetric = TRUE); W12 <- ev$vectors %*% diag(sqrt(ev$values)) %*% t(ev$vectors)
  l <- eigen(W12 %*% cov(C) %*% W12, symmetric = TRUE, only.values = TRUE)$values; m <- momentos(C, Phi)
  data.table(grupo = nombre, curvas = nrow(C), DE_media_mm_dia = round(mean(m$de), 2), pct_varianza_1a_componente = round(100 * l[1] / sum(l), 1),
             pct_varianza_2a_componente = round(100 * l[2] / sum(l), 1), cor_abr_oct = round(m$cor[105, 290], 2), cor_abr_jul = round(m$cor[105, 200], 2), cor_minima = round(min(m$cor), 2))
}
Ca_e <- anom_coef(C_e, id_e$codigo); Ca_s <- anom_coef(C_s, id_s$celda)
cov_grupos <- rbind(espectro(C_e, Phi_e, "Estaciones: todas"), espectro(C_e[id_e$vertiente != "Pacífico", ], Phi_e, "Estaciones: Cauca"),
                    espectro(C_e[id_e$vertiente == "Pacífico", ], Phi_e, "Estaciones: Pacífico"), espectro(Ca_e, Phi_e, "Estaciones: anomalías por estación"),
                    espectro(C_s, Phi_s, "Satélite: todos los píxeles"), espectro(C_s[id_s$vertiente != "Pacífico", ], Phi_s, "Satélite: Cauca"),
                    espectro(C_s[id_s$vertiente == "Pacífico", ], Phi_s, "Satélite: Pacífico"), espectro(Ca_s, Phi_s, "Satélite: anomalías por píxel"))
print(tabla_f(cov_grupos, "F3_covarianza_por_vertiente_y_anomalias"))
ma_e <- momentos(Ca_e, Phi_e); ma_s <- momentos(Ca_s, Phi_s)
cat(sprintf("Parte de la varianza que es diferencia ENTRE sitios: estaciones %.0f %% | satélite %.0f %%\n",
            100 * (1 - sum(ma_e$var) / sum(me$var)), 100 * (1 - sum(ma_s$var) / sum(ms$var))))
grafico_base("F3_correlacion_anomalias", function() {
  layout(matrix(1:3, 1), widths = c(5, 5, 1)); par(mar = c(4, 4, 3, 1)); colores <- hcl.colors(64, "Blue-Red 3")
  image(T_DIA, T_DIA, ma_e$cor, zlim = c(-1, 1), col = colores, main = "Correlación de las anomalías: estaciones", xlab = "Día s", ylab = "Día t")
  image(T_DIA, T_DIA, ma_s$cor, zlim = c(-1, 1), col = colores, main = "Correlación de las anomalías: satélite", xlab = "Día s", ylab = "Día t")
  par(mar = c(4, 0, 3, 4)); escala <- seq(-1, 1, length.out = 65)
  image(c(0, 1), escala, matrix(rep(escala, each = 2), nrow = 2), col = colores, zlim = c(-1, 1), axes = FALSE, xlab = "", ylab = ""); axis(4, las = 1, cex.axis = .8); box(); layout(1)
}, ancho = 13, alto = 5.5)

grafico_base("F3_media_y_dispersion", function() {
  par(mfrow = c(1, 3)); set.seed(1)
  for (f in list(list(Ye, media_e, me, "Estaciones", 4), list(Xs, media_s, ms, "Satélite (Valle)", "darkorange"))) {
    matplot(T_DIA, t(f[[1]][sample(nrow(f[[1]]), 250), ]), type = "l", lty = 1, col = "grey85", ylim = c(-2, 22),
            xlab = "Día del año", ylab = "mm/día", main = paste(f[[4]], ": media ± 2 DE"))
    lines(T_DIA, f[[2]], lwd = 3, col = f[[5]])
    lines(T_DIA, f[[2]] + 2 * f[[3]]$de, lty = 2, lwd = 2, col = 2); lines(T_DIA, pmax(f[[2]] - 2 * f[[3]]$de, -2), lty = 2, lwd = 2, col = 2)
  }
  plot(T_DIA, me$de, type = "l", lwd = 2, col = 4, ylim = range(c(0, me$de, ms$de)), xlab = "Día del año", ylab = "mm/día",
       main = "Desviación estándar funcional")
  lines(T_DIA, ms$de, lwd = 2, col = "darkorange")
  legend("topleft", c("Estaciones", "Satélite"), col = c(4, "darkorange"), lwd = 2, bty = "n")
}, ancho = 15, alto = 4.8)
for (f in list(list(me, "estaciones"), list(ms, "satelite"))) {
  grafico_base(paste0("F3_covarianza_", f[[2]]), function() mapa_calor(T_DIA, f[[1]]$cov, paste("Covarianza funcional:", f[[2]]), paleta = "YlOrRd"), ancho = 7, alto = 5.5)
  grafico_base(paste0("F3_correlacion_", f[[2]]), function() mapa_calor(T_DIA, f[[1]]$cor, paste("Correlación funcional:", f[[2]]), c(-1, 1)), ancho = 7, alto = 5.5)
}

titulo("3.3 Derivadas y temporadas de lluvia")
DPhi_e <- deriv_base(T_DIA, o$K_est); DPhi_s <- deriv_base(T_DIA, o$K_sat)
ref_d <- splines::splineDesign(nodos(o$K_est), T_DIA, ord = ORDEN, derivs = rep(1, 365), outer.ok = TRUE)
cat("Derivada de la base propia vs splines::splineDesign:", format(max(abs(DPhi_e - ref_d)), digits = 2), "\n")
dmu_e <- drop(DPhi_e %*% colMeans(C_e)); dmu_s <- drop(DPhi_s %*% colMeans(C_s))
cruces <- function(media, dmu, fuente) {
  k <- which(diff(sign(dmu)) != 0)
  data.table(fuente, dia = k + 1, fecha = fecha_dia(k + 1), tipo = fifelse(dmu[k] > 0, "máximo", "mínimo"), mm_dia = round(media[k + 1], 2))
}
extremos_media <- rbind(cruces(media_e, dmu_e, "estaciones"), cruces(media_s, dmu_s, "satélite"))
print(tabla_f(extremos_media, "F3_extremos_de_la_media"))

temporadas <- function(Y, D) {
  t(sapply(seq_len(nrow(Y)), function(i) {
    y <- Y[i, ]; d <- D[i, ]
    p1 <- 59 + which.max(y[60:200]); p2 <- 229 + which.max(y[230:365])
    r <- c(inicio_1 = if (p1 > 21) 19 + which.max(d[20:(p1 - 1)]) else NA_integer_, pico_1 = p1, fin_1 = p1 + which.min(d[(p1 + 1):229]),
           inicio_2 = 200 + which.max(d[201:(p2 - 1)]), pico_2 = p2, fin_2 = if (p2 < 349) p2 + which.min(d[(p2 + 1):350]) else NA_integer_)
    r[c("inicio_1", "fin_2")][r[c("inicio_1", "fin_2")] %in% c(20, 350)] <- NA
    r
  }))
}
tem_e <- temporadas(Ye, C_e %*% t(DPhi_e)); tem_s <- temporadas(Xs, C_s %*% t(DPhi_s))
fechas_tem <- rbind(data.table(fuente = "estaciones", vertiente = id_e$vertiente, tem_e), data.table(fuente = "satélite", vertiente = id_s$vertiente, tem_s))
fechas_tem <- melt(fechas_tem, id.vars = c("fuente", "vertiente"), variable.name = "momento", value.name = "dia")
resumen_tem <- fechas_tem[, .(pct_no_identificable = round(100 * mean(is.na(dia)), 1), dia_mediano = median(dia, na.rm = TRUE), p25 = quantile(dia, .25, na.rm = TRUE), p75 = quantile(dia, .75, na.rm = TRUE)),
                          by = .(fuente, vertiente, momento)]
resumen_tem[, `:=`(fecha_mediana = fecha_dia(dia_mediano), rango_intercuartil = paste(fecha_dia(p25), "a", fecha_dia(p75)))]
print(tabla_f(dcast(resumen_tem, momento ~ fuente + vertiente, value.var = "fecha_mediana"), "F3_temporadas_fechas"))
tabla_f(resumen_tem, "F3_temporadas_detalle")

grafico_base("F3_derivada_media", function() {
  par(mfrow = c(1, 2))
  for (f in list(list(media_e, dmu_e, "Estaciones", 4), list(media_s, dmu_s, "Satélite (Valle)", "darkorange"))) {
    plot(T_DIA, f[[2]], type = "l", lwd = 2, col = f[[4]], xlab = "Día del año", ylab = "mm/día por día", main = paste(f[[3]], ": derivada de la media"))
    abline(h = 0, lty = 2, col = "grey40")
    k <- which(diff(sign(f[[2]])) != 0) + 1; abline(v = k, lty = 3, col = "grey60")
    axis(3, at = k, labels = fecha_dia(k), cex.axis = .7, las = 2, tick = FALSE, line = -.8)
  }
}, ancho = 12, alto = 4.5)
grafico_f(ggplot(fechas_tem[!is.na(dia)], aes(dia, momento, fill = fuente)) + geom_boxplot(outlier.size = .2, alpha = .7) +
            facet_wrap(~vertiente) + scale_x_continuous(breaks = c(1, 60, 121, 182, 244, 305, 365), labels = fecha_dia(c(1, 60, 121, 182, 244, 305, 365))) +
            labs(x = NULL, y = NULL, fill = NULL, title = "Inicio, pico y fin de cada temporada en cada curva (derivada)"),
          "F3_temporadas_por_curva", ancho = 11, alto = 5)

titulo("3.4 Profundidad, mediana y media recortada")
Yp <- Ye[1:40, ]; pares <- combn(40, 2)
mbd_bruta <- sapply(1:40, function(i) mean(apply(pares, 2, function(p) mean(pmin(Yp[p[1], ], Yp[p[2], ]) <= Yp[i, ] & Yp[i, ] <= pmax(Yp[p[1], ], Yp[p[2], ])))))
cat("MBD por rangos vs definición con todas las bandas de dos curvas:", format(max(abs(profundidades(Yp)$MBD - mbd_bruta)), digits = 2), "\n")
pr_e <- profundidades(Ye); pr_s <- profundidades(Xs)
centralidad <- function(Y, pr) {
  dentro <- pr$FM >= quantile(pr$FM, .25)
  list(mediana_FM = Y[which.max(pr$FM), ], mediana_MBD = Y[which.max(pr$MBD), ], recortada_25 = colMeans(Y[dentro, , drop = FALSE]))
}
ce_e <- centralidad(Ye, pr_e); ce_s <- centralidad(Xs, pr_s)
cat("Curva más profunda (MBD), estaciones:", id_e$curva_id[which.max(pr_e$MBD)], id_e$nombre[which.max(pr_e$MBD)],
    "| (FM):", id_e$curva_id[which.max(pr_e$FM)], "| satélite (MBD): año", id_s$anio[which.max(pr_s$MBD)], "píxel", id_s$celda[which.max(pr_s$MBD)], "\n")
cat(sprintf("Correlación de rangos FM vs MBD: estaciones %.2f | satélite %.2f\n",
            cor(pr_e$FM, pr_e$MBD, method = "spearman"), cor(pr_s$FM, pr_s$MBD, method = "spearman")))
grafico_base("F3_centralidad", function() {
  par(mfrow = c(1, 2))
  for (f in list(list(media_e, ce_e, "Estaciones"), list(media_s, ce_s, "Satélite (Valle)"))) {
    matplot(T_DIA, cbind(f[[1]], f[[2]]$mediana_FM, f[[2]]$mediana_MBD, f[[2]]$recortada_25), type = "l", lty = c(1, 2, 3, 1),
            col = c(4, 2, 1, 3), lwd = 2, xlab = "Día del año", ylab = "mm/día", main = paste("Centralidad:", f[[3]]))
    legend("topright", c("Media", "Mediana FM", "Mediana MBD", "Media recortada 25 % (FM)"), col = c(4, 2, 1, 3), lty = c(1, 2, 3, 1), lwd = 2, bty = "n", cex = .75)
  }
}, ancho = 12, alto = 4.5)

titulo("3.5 Boxplot funcional y curvas atípicas")
bx_e <- boxplot_funcional(Ye, pr_e$MBD); bx_s <- boxplot_funcional(Xs, pr_s$MBD)
grupos_v <- split(seq_len(nrow(Ye)), id_e$vertiente)
bx_v <- lapply(grupos_v, function(i) boxplot_funcional(Ye[i, ], profundidades(Ye[i, ])$MBD))
fb_e <- fda::fbplot(t(Ye), T_DIA, depth = pr_e$MBD, plot = FALSE)
fb_s <- fda::fbplot(t(Xs), T_DIA, depth = pr_s$MBD, plot = FALSE)
verif_fb <- c(estaciones = setequal(as.integer(fb_e$outpoint), which(bx_e$atip)), satelite = setequal(as.integer(fb_s$outpoint), which(bx_s$atip)))
cat("Atípicas iguales a fda::fbplot:", paste(names(verif_fb), verif_fb, collapse = " | "), "\n")
stopifnot(verif_fb)

atip_e <- id_e[bx_e$atip, .(codigo, nombre, anio, vertiente, altitud, particion)]
atip_e[, `:=`(MBD = round(pr_e$MBD[bx_e$atip], 3), P_anual_mm = round(rowMeans(Ye[bx_e$atip, , drop = FALSE]) * 365), dia_max = apply(Ye[bx_e$atip, , drop = FALSE], 1, which.max))]
atip_v <- rbindlist(lapply(names(grupos_v), function(v) id_e[grupos_v[[v]][bx_v[[v]]$atip], .(codigo, nombre, anio, vertiente)]))
tabla_f(atip_e[order(-P_anual_mm)], "F3_atipicas_estaciones_todas")
tabla_f(atip_v, "F3_atipicas_estaciones_por_vertiente")
cat(sprintf("Atípicas, todas las estaciones juntas: %d de %d curvas (%d estaciones; %d del Pacífico)\n",
            nrow(atip_e), nrow(Ye), uniqueN(atip_e$codigo), atip_e[vertiente == "Pacífico", .N]))
print(head(atip_e[, .(curvas = .N, anios = paste(sort(anio), collapse = " "), P_anual_mediana = median(P_anual_mm)), by = .(codigo, nombre, vertiente)][order(-curvas)], 12))
cat("Atípicas dentro de cada vertiente:\n"); print(atip_v[, .(curvas = .N, estaciones = uniqueN(codigo)), by = vertiente])
print(atip_v[vertiente != "Pacífico", .N, by = .(codigo, nombre, anio)])
atip_s <- id_s[bx_s$atip]
cat(sprintf("Atípicas del satélite: %d de %d curvas píxel-año (%.1f %%)\n", nrow(atip_s), nrow(Xs), 100 * mean(bx_s$atip)))
print(atip_s[, .N, by = .(vertiente)]); print(head(atip_s[, .N, by = anio][order(-N)], 6))
tabla_f(atip_s[, .(anio, celda, lon, lat, vertiente, n_tope)], "F3_atipicas_satelite")

anios_e <- sort(unique(id_e$anio))
media_anual <- function(Y, g, a) t(sapply(a, function(x) colMeans(Y[g == x, , drop = FALSE])))
Ya_e <- media_anual(Ye[id_e$vertiente != "Pacífico", ], id_e$anio[id_e$vertiente != "Pacífico"], anios_e)
Ya_s <- media_anual(Xs, id_s$anio, ANIOS)
pa_e <- profundidades(Ya_e)$MBD; pa_s <- profundidades(Ya_s)$MBD
bxa_e <- boxplot_funcional(Ya_e, pa_e); bxa_s <- boxplot_funcional(Ya_s, pa_s)
prof_anual <- oni[data.table(anio = anios_e, MBD_estaciones_Cauca = round(pa_e, 3), atipico_estaciones = bxa_e$atip,
                             MBD_satelite = round(pa_s, 3), atipico_satelite = bxa_s$atip,
                             P_estaciones_Cauca = round(rowMeans(Ya_e) * 365), P_satelite = round(rowMeans(Ya_s) * 365)), on = "anio"]
print(tabla_f(prof_anual[order(MBD_estaciones_Cauca)], "F3_profundidad_anual"))
cat(sprintf("Correlación (Spearman) entre la profundidad anual de estaciones y de satélite: %.2f\n", cor(pa_e, pa_s, method = "spearman")))

grafico_base("F3_boxplot_funcional", function() {
  par(mfrow = c(2, 3))
  dibujar_boxplot(Ye, bx_e, paste("Estaciones, todas:", sum(bx_e$atip), "atípicas"))
  for (v in names(bx_v)) dibujar_boxplot(Ye[grupos_v[[v]], ], bx_v[[v]], paste("Estaciones,", v, ":", sum(bx_v[[v]]$atip), "atípicas"))
  dibujar_boxplot(Xs, bx_s, paste("Satélite (Valle):", sum(bx_s$atip), "atípicas"))
  dibujar_boxplot(Ya_e, bxa_e, "Curvas medias anuales: estaciones del Cauca")
  dibujar_boxplot(Ya_s, bxa_s, "Curvas medias anuales: satélite")
}, ancho = 15, alto = 9)

titulo("3.6 Contexto: vertiente, altitud y ENSO")
id_e[, banda_altitud := cut(altitud, c(-Inf, 1200, 2000, Inf), labels = c("< 1200 m", "1200-2000 m", "> 2000 m"))]
id_e <- oni[, .(anio, fase)][id_e, on = "anio"]; id_s <- oni[, .(anio, fase)][id_s, on = "anio"]
grupo_media <- function(Y, g) { g <- as.character(g); k <- !is.na(g); sapply(split(seq_len(nrow(Y))[k], g[k]), function(i) colMeans(Y[i, , drop = FALSE])) }
m_vert_e <- grupo_media(Ye, id_e$vertiente); m_vert_s <- grupo_media(Xs, id_s$vertiente)
m_alt <- grupo_media(Ye[id_e$vertiente != "Pacífico", ], id_e$banda_altitud[id_e$vertiente != "Pacífico"])
m_enso_e <- grupo_media(Ye[id_e$vertiente != "Pacífico", ], id_e$fase[id_e$vertiente != "Pacífico"])
m_enso_s <- grupo_media(Xs[id_s$vertiente != "Pacífico", ], id_s$fase[id_s$vertiente != "Pacífico"])
k_sd <- id_s$vertiente != "Pacífico" & !id_s$marcado
m_enso_s_sd <- grupo_media(Xs[k_sd, ], id_s$fase[k_sd])
contexto <- rbind(
  data.table(fuente = "estaciones", agrupacion = "vertiente", grupo = colnames(m_vert_e), P_anual = colMeans(m_vert_e) * 365,
             curvas = as.vector(table(id_e$vertiente)[colnames(m_vert_e)])),
  data.table(fuente = "satélite", agrupacion = "vertiente", grupo = colnames(m_vert_s), P_anual = colMeans(m_vert_s) * 365,
             curvas = as.vector(table(id_s$vertiente)[colnames(m_vert_s)])),
  data.table(fuente = "estaciones (Cauca)", agrupacion = "altitud", grupo = colnames(m_alt), P_anual = colMeans(m_alt) * 365,
             curvas = as.vector(table(id_e[vertiente != "Pacífico", banda_altitud])[colnames(m_alt)])),
  data.table(fuente = "estaciones (Cauca)", agrupacion = "ENSO", grupo = colnames(m_enso_e), P_anual = colMeans(m_enso_e) * 365,
             curvas = as.vector(table(id_e[vertiente != "Pacífico", fase])[colnames(m_enso_e)])),
  data.table(fuente = "satélite (Cauca)", agrupacion = "ENSO", grupo = colnames(m_enso_s), P_anual = colMeans(m_enso_s) * 365,
             curvas = as.vector(table(id_s[vertiente != "Pacífico", fase])[colnames(m_enso_s)])),
  data.table(fuente = "satélite (Cauca) sin 1985, 1988 y 2000", agrupacion = "ENSO", grupo = colnames(m_enso_s_sd), P_anual = colMeans(m_enso_s_sd) * 365,
             curvas = as.vector(table(id_s[k_sd, fase])[colnames(m_enso_s_sd)])))
cat("Años marcados (1985, 1988 y 2000) y su fase ENSO:\n"); print(oni[anio %in% ANIOS_MARCADOS])
contexto[, P_anual := round(P_anual)]
print(tabla_f(contexto, "F3_contexto_vertiente_altitud_ENSO"))
picos_enso <- rbindlist(lapply(colnames(m_enso_e), function(f) data.table(fase = f, pico_1 = fecha_dia(59 + which.max(m_enso_e[60:200, f])),
  mm_pico_1 = round(max(m_enso_e[60:200, f]), 2), pico_2 = fecha_dia(229 + which.max(m_enso_e[230:365, f])), mm_pico_2 = round(max(m_enso_e[230:365, f]), 2),
  minimo_jul_ago = round(min(m_enso_e[160:260, f]), 2))))
print(tabla_f(picos_enso, "F3_picos_por_fase_ENSO"))

grafico_base("F3_contexto", function() {
  par(mfrow = c(1, 3))
  matplot(T_DIA, cbind(m_vert_e, m_vert_s), type = "l", lty = c(1, 1, 2, 2), lwd = 2, col = c(4, "darkgreen", 4, "darkgreen"),
          xlab = "Día del año", ylab = "mm/día", main = "Media por vertiente")
  legend("topleft", c(paste("Estaciones,", colnames(m_vert_e)), paste("Satélite,", colnames(m_vert_s))), col = c(4, "darkgreen", 4, "darkgreen"),
         lty = c(1, 1, 2, 2), lwd = 2, bty = "n", cex = .75)
  matplot(T_DIA, m_alt, type = "l", lty = 1, lwd = 2, col = c("goldenrod", "darkorange3", "brown"), xlab = "Día del año", ylab = "mm/día",
          main = "Estaciones del Cauca por altitud")
  legend("topleft", colnames(m_alt), col = c("goldenrod", "darkorange3", "brown"), lwd = 2, bty = "n")
  matplot(T_DIA, m_enso_e, type = "l", lty = 1, lwd = 2, col = c(1, 4, 2)[match(colnames(m_enso_e), c("Neutro", "Niña", "Niño"))],
          xlab = "Día del año", ylab = "mm/día", main = "Estaciones del Cauca por fase ENSO")
  legend("topleft", colnames(m_enso_e), col = c(1, 4, 2)[match(colnames(m_enso_e), c("Neutro", "Niña", "Niño"))], lwd = 2, bty = "n")
}, ancho = 15, alto = 4.8)

# QUÉ DICEN LOS DATOS
#   Representatividad: se usan las 176 estaciones incluidas y sus 3 299 curvas de
#   1981-2003, sin submuestrear (la sección 8 del archivo de datos muestra que
#   las incluidas representan a las 220 con datos), y los 676 píxeles del Valle
#   en los 23 años.
#   Media: Phi c_barra coincide con fda::mean.fd (5e-15). Ambas fuentes muestran
#   el régimen bimodal: estaciones con máximos el 22-abr (6,3 mm/día) y el 28-oct
#   (7,2), mínimos el 25-ene y el 1-ago (2,2). El satélite pone los máximos 1 a 3
#   semanas después (13-may y 5-nov) y con casi el doble de magnitud, porque su
#   media pesa por área y la vertiente Pacífica (43 % de los píxeles) llueve
#   4 000 mm/año. Por vertiente, el satélite supera a las estaciones en ambas:
#   1 735 frente a 1 551 mm en el Cauca y 3 998 frente a 3 105 en el Pacífico.
#   Dispersión: la DE de las estaciones es alta en la segunda temporada y mínima
#   el 27-feb; el máximo del 1-ene es efecto de borde. El CV medio es 0,96: la
#   variación entre curvas es tan grande como la media.
#   Qué mide la covarianza: con todas las estaciones juntas la correlación es
#   positiva todo el año (mínimo 0,39) y la primera componente explica 76 %. Eso
#   es sobre todo "el nivel de lluvia de cada sitio": el 74 % de la varianza de
#   las estaciones (88 % en el satélite) es diferencia ENTRE sitios. Dentro del
#   Cauca la primera componente baja a 50 % y la correlación abril-octubre a
#   0,54. Con anomalías (cada curva menos la media de su estación) queda la
#   variación de un año a otro: la primera componente explica solo 29 %, abril y
#   octubre casi no se relacionan (0,14) y aparecen correlaciones negativas
#   (-0,27). Es decir: un sitio húmedo lo es todo el año, pero un año con una
#   primera temporada lluviosa no anuncia una segunda temporada lluviosa.
#   Derivadas y temporadas: en el Cauca la primera temporada empieza a finales de
#   febrero, tiene su pico a mediados de abril y termina a mediados de junio; la
#   segunda empieza a mediados de septiembre y tiene su pico el 1-nov. En 62 % de
#   las curvas de estación la lluvia sigue bajando hasta el borde del año y el
#   fin de la segunda temporada no se puede fechar. El satélite fecha el inicio y
#   el pico de la primera temporada 2 a 3 semanas tarde.
#   Profundidad: FM y MBD ordenan casi igual (Spearman 0,99 y 1,00). La curva más
#   profunda es Betania 2002. Boxplot funcional (mismas atípicas que fda::fbplot):
#   con todas las estaciones juntas hay 63 atípicas, 62 de cinco estaciones de la
#   costa Pacífica (Bahía Málaga, Yurumanguí, San Francisco, Triana y El Tabano,
#   7 000 a 12 600 mm/año): no son errores sino otro clima. Por eso se separa por
#   vertiente: en el Cauca solo quedan Tenjo 1982 (4 482 mm, candidata a revisión)
#   y Las Brisas 1981. Las 330 atípicas del satélite son todas del Pacífico y 229
#   tocan el tope de 200 mm. Ningún año es atípico en las curvas medias anuales;
#   los menos profundos de las estaciones son 1988, 1984, 1983 y 1999 (tres Niña).
#   Orografía: en el Cauca la franja de 1 200-2 000 m es la más lluviosa
#   (1 715 mm), por encima de > 2 000 m (1 607) y del fondo del valle (1 333).
#   ENSO: en las estaciones del Cauca La Niña aporta +11 % (1 718 mm) y El Niño
#   -10 % (1 387) frente a neutro (1 544); la Niña se nota en la segunda temporada
#   (7,1 frente a 6,3 mm/día) y en una sequía de mitad de año más húmeda (2,4
#   frente a 1,4). El satélite NO muestra esa señal: da menos lluvia en años Niña
#   (1 664) que en neutros (1 767). Es el mismo error ligado al ENSO de la
#   sección 1, y es la razón para darle el ONI al modelo.

eda <- list(media_e = media_e, media_s = media_s, me = me, ms = ms, pr_e = pr_e, pr_s = pr_s, atip_e = atip_e, atip_s = atip_s,
             temporadas = resumen_tem, prof_anual = prof_anual, contexto = contexto)

# =============================================================================
# 4. COMPARACIÓN FUNCIONAL POR FUENTE (ESTACIÓN VS. SU PÍXEL)  — rúbrica 3
# =============================================================================
# Pares estación-año / píxel-año con el emparejamiento de la sección 1.
# Curva de diferencia d(t) = estación - satélite con banda d_barra(t) ± 2EE(t);
# el EE se calcula entre estaciones (las curvas de una misma estación no son
# independientes). Correlación cruzada funcional y su diagonal, covarianzas de
# cada fuente, comparación por vertiente, altitud y fase ENSO. Entran todos los
# años; 1985, 1988 y 2000 se muestran además aparte, como descripción.
titulo("4. Comparación funcional por fuente")

Phi_e <- base_bspline(T_DIA, o$K_est); Phi_s <- base_bspline(T_DIA, o$K_sat)

pares <- o$emparejamiento[, .(codigo, nombre, vertiente, altitud, celda, dist_km, estaciones_en_el_pixel)][o$id_est[, .(fila_est = .I, codigo, anio, particion, curva_id)], on = "codigo"]
pares <- o$id_sat[, .(fila_sat = .I, anio, celda, marcado)][pares, on = .(anio, celda)]
pares <- oni[, .(anio, ONI_anual, fase)][pares, on = "anio"]
pares[, banda_altitud := cut(altitud, c(-Inf, 1200, 2000, Inf), labels = c("< 1200 m", "1200-2000 m", "> 2000 m"))]
stopifnot(!anyNA(pares$fila_sat), nrow(pares) == nrow(o$id_est))
cat("Pares estación-año / píxel-año:", nrow(pares), "en", uniqueN(pares$codigo), "estaciones y", uniqueN(pares$celda), "píxeles |",
    "en 1985, 1988 y 2000:", pares[marcado == TRUE, .N], "\n")
print(pares[, .(pares = .N, estaciones = uniqueN(codigo)), by = .(particion, marcado)][order(particion, marcado)])
tabla_f(pares[, .(codigo, nombre, anio, particion, celda, dist_km, estaciones_en_el_pixel, vertiente, altitud, ONI_anual, fase, marcado)], "F4_pares_estacion_pixel")

Y <- o$coef_est[pares$fila_est, ] %*% t(Phi_e)
X <- o$coef_sat[pares$fila_sat, ] %*% t(Phi_s)
D <- Y - X
ok <- rep(TRUE, nrow(pares))   # todos los años entran (D-47)
marc <- pares$marcado

banda_diferencia <- function(D, grupo) {
  medias <- rowsum(D, grupo) / as.vector(table(grupo)[rownames(rowsum(D, grupo))])
  n <- nrow(medias)
  list(d = colMeans(medias), ee = apply(medias, 2, sd) / sqrt(n), n = n, ee_ingenuo = apply(D, 2, sd) / sqrt(nrow(D)))
}
intervalos <- function(inf, sup) {
  signo <- ifelse(inf > 0, "estación > satélite", ifelse(sup < 0, "satélite > estación", NA))
  r <- rle(signo); fin <- cumsum(r$lengths); ini <- fin - r$lengths + 1
  data.table(desde = fecha_dia(ini), hasta = fecha_dia(fin), dias = r$lengths, signo = r$values)[!is.na(signo)]
}

b <- banda_diferencia(D[ok, ], pares$codigo[ok])
cat(sprintf("EE por estación / EE tratando cada curva como independiente: %.2f en promedio\n", mean(b$ee / b$ee_ingenuo)))
int_d <- intervalos(b$d - 2 * b$ee, b$d + 2 * b$ee)
print(tabla_f(int_d, "F4_intervalos_diferencia"))
cat(sprintf("d(t) medio: %.2f mm/día (%.0f mm/año); días en que la banda excluye el 0: %d de 365\n",
            mean(b$d), mean(b$d) * 365, sum(int_d$dias)))
b_def <- banda_diferencia(D[marc, ], pares$codigo[marc]); b_resto <- banda_diferencia(D[!marc, ], pares$codigo[!marc])
cat(sprintf("En 1985, 1988 y 2000: d(t) medio %.2f mm/día (%.0f mm/año) | en el resto: %.2f mm/día (%.0f mm/año)\n", mean(b_def$d), mean(b_def$d) * 365, mean(b_resto$d), mean(b_resto$d) * 365))
tabla_f(data.table(dia = T_DIA, d = b$d, ee = b$ee, inferior = b$d - 2 * b$ee, superior = b$d + 2 * b$ee, d_1985_1988_2000 = b_def$d, d_resto = b_resto$d), "F4_diferencia_por_dia")

grafico_base("F4_diferencia", function() {
  par(mfrow = c(1, 2))
  plot(T_DIA, colMeans(Y[ok, ]), type = "l", lwd = 2, col = 4, ylim = range(c(0, colMeans(Y[ok, ]), colMeans(X[ok, ]))),
       xlab = "Día del año", ylab = "mm/día", main = "Media de los pares: estación y su píxel")
  lines(T_DIA, colMeans(X[ok, ]), lwd = 2, col = "darkorange")
  legend("bottomleft", c("Estaciones", "Satélite (su píxel)"), col = c(4, "darkorange"), lwd = 2, bty = "n")
  plot(T_DIA, b$d, type = "n", ylim = range(c(b$d - 2 * b$ee, b$d + 2 * b$ee, b_def$d, 0)), xlab = "Día del año", ylab = "mm/día",
       main = expression(paste(bar(d), "(t) = estación - satélite  ±2 EE(t)")))
  polygon(c(T_DIA, rev(T_DIA)), c(b$d - 2 * b$ee, rev(b$d + 2 * b$ee)), col = adjustcolor("steelblue", .3), border = NA)
  lines(T_DIA, b$d, lwd = 2); lines(T_DIA, b_def$d, lty = 2, col = "firebrick"); abline(h = 0, lty = 3)
  excl <- (b$d - 2 * b$ee > 0) | (b$d + 2 * b$ee < 0); rug(T_DIA[excl], col = "firebrick")
  legend("bottomleft", c("Todos los años (banda por estación)", "Solo 1985, 1988 y 2000", "Días en que la banda excluye el 0"),
         col = c(1, "firebrick", "firebrick"), lty = c(1, 2, 1), lwd = c(2, 1, 1), bty = "n", cex = .8)
}, ancho = 13, alto = 4.8)

titulo("4.2 Correlación cruzada y covarianzas")
Ce <- o$coef_est[pares$fila_est[ok], ]; Cs <- o$coef_sat[pares$fila_sat[ok], ]
Cec <- sweep(Ce, 2, colMeans(Ce)); Csc <- sweep(Cs, 2, colMeans(Cs)); n <- nrow(Ce)
Cov_YY <- Phi_e %*% (crossprod(Cec) / (n - 1)) %*% t(Phi_e)
Cov_XX <- Phi_s %*% (crossprod(Csc) / (n - 1)) %*% t(Phi_s)
Cov_YX <- Phi_e %*% (crossprod(Cec, Csc) / (n - 1)) %*% t(Phi_s)
R_YX <- Cov_YX / sqrt(outer(diag(Cov_YY), diag(Cov_XX)))
cat("Covarianza cruzada por coeficientes vs cov() de las curvas evaluadas:", format(max(abs(Cov_YX - cov(Y[ok, ], X[ok, ]))), digits = 2), "\n")
diag_r <- diag(R_YX)
cat(sprintf("Diagonal r_XY(t,t): media %.2f, máxima %.2f el %s, mínima %.2f el %s\n",
            mean(diag_r), max(diag_r), fecha_dia(which.max(diag_r)), min(diag_r), fecha_dia(which.min(diag_r))))
mejor_s <- apply(R_YX, 1, which.max)
cat(sprintf("Para cada día t de la estación, el día s del satélite más correlacionado está a %+.1f días en mediana\n", median(mejor_s - T_DIA)))

anom <- function(M, g) M - (rowsum(M, g) / as.vector(table(g)[rownames(rowsum(M, g))]))[as.character(g), ]
Ya <- anom(Y[ok, ], pares$codigo[ok]); Xa <- anom(X[ok, ], pares$codigo[ok])
R_anom <- diag(cov(Ya, Xa)) / sqrt(apply(Ya, 2, var) * apply(Xa, 2, var))
cat(sprintf("Diagonal con anomalías (restando la climatología de cada estación): media %.2f, máxima %.2f el %s, mínima %.2f el %s\n",
            mean(R_anom), max(R_anom), fecha_dia(which.max(R_anom)), min(R_anom), fecha_dia(which.min(R_anom))))
tabla_f(data.table(dia = T_DIA, r_diagonal = diag_r, r_diagonal_anomalias = R_anom, var_estacion = diag(Cov_YY), var_satelite = diag(Cov_XX)), "F4_diagonal_y_varianzas")

cov_comp <- data.table(fuente = c("estación", "satélite"), varianza_total = c(sum(diag(Cov_YY)), sum(diag(Cov_XX))),
                       DE_media = c(mean(sqrt(diag(Cov_YY))), mean(sqrt(diag(Cov_XX)))),
                       cor_abr_oct = c(cov2cor(Cov_YY)[105, 290], cov2cor(Cov_XX)[105, 290]),
                       cor_abr_jul = c(cov2cor(Cov_YY)[105, 200], cov2cor(Cov_XX)[105, 200]),
                       var_explicada_1a_componente = c(eigen(Cov_YY, symmetric = TRUE, only.values = TRUE)$values[1] / sum(diag(Cov_YY)),
                                                       eigen(Cov_XX, symmetric = TRUE, only.values = TRUE)$values[1] / sum(diag(Cov_XX))))
print(tabla_f(cov_comp, "F4_covarianzas_comparadas"))
cat(sprintf("Cociente de varianzas satélite / estación por día: mediana %.2f (rango %.2f a %.2f)\n",
            median(diag(Cov_XX) / diag(Cov_YY)), min(diag(Cov_XX) / diag(Cov_YY)), max(diag(Cov_XX) / diag(Cov_YY))))

# ¿Las diferencias entre fuentes salen de usar K distinto (8 y 11)? Se repite la
# comparación con el satélite construido con el mismo K de las estaciones.
X_k <- construir(X_pen[pares$fila_sat, ], T_SAT, o$K_est)$coef %*% t(Phi_e); b_k <- banda_diferencia(Y - X_k, pares$codigo)
excluye <- function(bb) sum((bb$d - 2 * bb$ee > 0) | (bb$d + 2 * bb$ee < 0))
mismo_K <- data.table(K_estacion = o$K_est, K_satelite = c(o$K_sat, o$K_est),
                      d_medio_mm_anio = round(c(mean(b$d), mean(b_k$d)) * 365), dias_banda_sin_cero = c(excluye(b), excluye(b_k)),
                      dias_con_el_mismo_veredicto = c(365L, sum(((b$d - 2 * b$ee > 0) - (b$d + 2 * b$ee < 0)) == ((b_k$d - 2 * b_k$ee > 0) - (b_k$d + 2 * b_k$ee < 0)))),
                      r_diagonal_media = round(c(mean(diag_r), mean(diag(cor(Y, X_k)))), 3),
                      cociente_varianzas_mediana = round(c(median(diag(Cov_XX) / diag(Cov_YY)), median(apply(X_k, 2, var) / apply(Y, 2, var))), 3))
print(tabla_f(mismo_K, "F4_sensibilidad_mismo_K"))

grafico_base("F4_correlacion_cruzada", function() {
  par(mfrow = c(1, 1)); mapa_calor(T_DIA, R_YX, "Correlación cruzada r(s,t): satélite en s, estación en t", c(-1, 1), xlab = "Día t (estación)", ylab = "Día s (satélite)")
}, ancho = 7, alto = 5.5)
grafico_base("F4_diagonal_correlacion", function() {
  plot(T_DIA, diag_r, type = "l", lwd = 2, ylim = c(0, 1), xlab = "Día del año", ylab = "r(t, t)", main = "Diagonal de la correlación cruzada")
  lines(T_DIA, R_anom, lwd = 2, col = "firebrick")
  legend("bottomleft", c("Curvas", "Anomalías (sin la climatología de cada estación)"), col = c(1, "firebrick"), lwd = 2, bty = "n")
}, ancho = 9, alto = 4.5)
lim_cov <- range(c(Cov_YY, Cov_XX))
grafico_base("F4_covarianzas", function() {
  layout(matrix(1:3, 1)); par(mar = c(4, 4, 3, 1))
  image(T_DIA, T_DIA, Cov_YY, zlim = lim_cov, col = hcl.colors(64, "YlOrRd", rev = TRUE), main = "Covarianza: estaciones", xlab = "Día s", ylab = "Día t")
  image(T_DIA, T_DIA, Cov_XX, zlim = lim_cov, col = hcl.colors(64, "YlOrRd", rev = TRUE), main = "Covarianza: su píxel", xlab = "Día s", ylab = "Día t")
  lim_d <- max(abs(Cov_XX - Cov_YY)) * c(-1, 1)
  image(T_DIA, T_DIA, Cov_XX - Cov_YY, zlim = lim_d, col = hcl.colors(64, "Blue-Red 3"), main = "Satélite - estación", xlab = "Día s", ylab = "Día t")
  layout(1)
}, ancho = 15, alto = 4.8)

titulo("4.3 Por vertiente, altitud y fase ENSO")
por_grupo <- function(var) {
  g <- as.character(pares[[var]]); k <- ok & !is.na(g)
  rbindlist(lapply(sort(unique(g[k])), function(v) {
    i <- k & g == v; bb <- banda_diferencia(D[i, , drop = FALSE], pares$codigo[i])
    data.table(agrupacion = var, grupo = v, pares = sum(i), estaciones = bb$n,
               P_estacion = mean(Y[i, ]) * 365, P_satelite = mean(X[i, ]) * 365, sesgo_mm_anio = mean(bb$d) * 365,
               cociente = mean(X[i, ]) / mean(Y[i, ]), dias_banda_sin_cero = sum((bb$d - 2 * bb$ee > 0) | (bb$d + 2 * bb$ee < 0)),
               r_diagonal_media = mean(diag(cor(Y[i, ], X[i, ]))))
  }))
}
grupos <- rbind(por_grupo("vertiente"), por_grupo("banda_altitud"), por_grupo("fase"))
grupos[, c("P_estacion", "P_satelite", "sesgo_mm_anio") := lapply(.SD, round), .SDcols = c("P_estacion", "P_satelite", "sesgo_mm_anio")]
grupos[, c("cociente", "r_diagonal_media") := lapply(.SD, round, 2), .SDcols = c("cociente", "r_diagonal_media")]
print(tabla_f(grupos, "F4_comparacion_por_grupo"))

pares[, `:=`(P_est = rowMeans(Y) * 365, P_sat = rowMeans(X) * 365)]
por_est <- pares[ok, .(P_est = mean(P_est), P_sat = mean(P_sat), altitud = altitud[1], vertiente = vertiente[1], nombre = nombre[1]), by = codigo]
por_est[, cociente := P_sat / P_est]
ajuste_alt <- lm(log(cociente) ~ altitud, data = por_est[vertiente != "Pacífico"])
cat(sprintf("Estaciones del Cauca: log(satélite/estación) cambia %.3f por cada 1 000 m de altitud (p = %.3g); r = %.2f\n",
            1000 * coef(ajuste_alt)[2], summary(ajuste_alt)$coefficients[2, 4], cor(por_est[vertiente != "Pacífico" & !is.na(altitud), .(altitud, log(cociente))])[1, 2]))
cat(sprintf("Lluvia anual estación vs su píxel (promedio por estación): r = %.2f; en el Cauca r = %.2f\n",
            cor(por_est$P_est, por_est$P_sat), por_est[vertiente != "Pacífico", cor(P_est, P_sat)]))
tabla_f(por_est[order(cociente)], "F4_sesgo_por_estacion")
print(rbind(head(por_est[order(cociente)], 4), tail(por_est[order(cociente)], 4)))

grafico_f(ggplot(por_est, aes(P_est, P_sat, colour = vertiente)) + geom_abline(linetype = 2) + geom_point() +
            scale_x_log10() + scale_y_log10() + labs(x = "Estación (mm/año, log)", y = "Su píxel (mm/año, log)", colour = NULL,
                                                   title = "Lluvia anual media: estación frente a su píxel (1981-2003)"),
          "F4_anual_estacion_vs_pixel", ancho = 8, alto = 6)
grafico_f(ggplot(por_est[!is.na(altitud)], aes(altitud, cociente, colour = vertiente)) + geom_hline(yintercept = 1, linetype = 2) +
            geom_point() + geom_smooth(method = "lm", formula = y ~ x, se = FALSE) + scale_y_log10() +
            labs(x = "Altitud de la estación (m)", y = "Satélite / estación (log)", colour = NULL, title = "Sesgo relativo del satélite según la altitud"),
          "F4_sesgo_altitud", ancho = 8, alto = 5)

curvas_d <- rbindlist(lapply(c("vertiente", "banda_altitud", "fase"), function(var) {
  g <- as.character(pares[[var]]); k <- ok & !is.na(g)
  rbindlist(lapply(sort(unique(g[k])), function(v) { i <- k & g == v; bb <- banda_diferencia(D[i, , drop = FALSE], pares$codigo[i])
    data.table(agrupacion = var, grupo = v, dia = T_DIA, d = bb$d, inf = bb$d - 2 * bb$ee, sup = bb$d + 2 * bb$ee) }))
}))
grafico_f(ggplot(curvas_d, aes(dia, d, colour = grupo, fill = grupo)) + geom_hline(yintercept = 0, linetype = 3) +
            geom_ribbon(aes(ymin = inf, ymax = sup), alpha = .15, colour = NA) + geom_line(linewidth = .8) +
            facet_wrap(~agrupacion, ncol = 3, labeller = as_labeller(c(vertiente = "Vertiente", banda_altitud = "Altitud", fase = "Fase ENSO"))) + labs(x = "Día del año", y = "d(t), mm/día", colour = NULL, fill = NULL,
                                                      title = "Diferencia estación - satélite por grupo (±2 EE por estación)"),
          "F4_diferencia_por_grupo", ancho = 14, alto = 4.8)

# QUÉ DICEN LOS DATOS
#   Pares: 3 299 estación-año / píxel-año, 176 estaciones en 148 píxeles, todos
#   los años. El EE se calcula entre estaciones (promedio de cada estación
#   primero): es 3,7 veces el que saldría tratando cada curva como independiente,
#   porque los años de una misma estación repiten su sesgo.
#   d(t): en el promedio anual casi no hay diferencia (el satélite da 22 mm/año
#   más), pero el signo cambia con la temporada. La estación supera al satélite
#   del 8-feb al 28-mar (inicio de la primera temporada: el satélite la arranca
#   tarde); el satélite supera a la estación del 15-abr al 9-jun y del 28-oct al
#   18-nov, es decir, exagera los picos de las dos temporadas. Es el mismo
#   retraso que mostraron las derivadas. La banda excluye el 0 en 147 de 365
#   días. Los tramos del 1 al 12 de enero y del 24 al 31 de diciembre son efecto
#   de borde de las bases. En 1985, 1988 y 2000 el satélite pierde 337 mm/año
#   frente a la estación; en el resto de años da 83 mm/año de más.
#   Correlación cruzada: la diagonal r(t,t) vale 0,62 en promedio (0,78 a
#   mediados de julio) y la correlación máxima está en el mismo día (desfase
#   mediano de -1 día). Casi toda viene de que el satélite reconoce qué sitios
#   son húmedos y cuáles secos: con anomalías (restando la climatología de cada
#   estación) baja a 0,24, con un máximo de 0,47 en febrero y 0,07 en noviembre.
#   El satélite sigue mal la variación de un año a otro, sobre todo en la segunda
#   temporada.
#   Covarianzas: la forma es casi igual (correlación abril-octubre 0,74 y 0,72;
#   la primera componente explica 75 % y 77 %), pero la varianza del satélite es
#   una cuarta parte de la de las estaciones (mediana 0,24 por día): el píxel de
#   ~5,5 km y la péntada promedian los aguaceros locales.
#   Mismo K: al construir el satélite también con K = 8, la diferencia media
#   (-20 frente a -22 mm/año), la diagonal (0,64 frente a 0,62) y el cociente de
#   varianzas (0,25 frente a 0,24) casi no cambian, y el veredicto de la banda coincide en 281
#   de 365 días. Las diferencias entre fuentes no salen de usar K distinto.
#   Subregión: en el Cauca el satélite sobreestima 8 % (+109 mm/año) y la banda
#   excluye el 0 en 263 días. En la vertiente Pacífica el sesgo depende del
#   sitio. El satélite subestima a la mitad las estaciones más lluviosas (El
#   Tabano, San Francisco y Bahía Málaga, cociente 0,43 a 0,50) y duplica la
#   lluvia de los sitios secos a sotavento (Loboguerrero 2,6; La Despensa 2,3;
#   Puerto Nuevo 2,3; Moravito 1,9). El satélite suaviza hacia la media regional:
#   sube lo seco y baja lo muy húmedo. Entre estaciones, la lluvia anual del píxel
#   explica la de la estación con r = 0,79, pero dentro del Cauca solo con 0,42.
#   Altitud: en el Cauca el cociente satélite/estación baja 15 % por cada 1 000 m
#   (p = 0,001): el satélite exagera en el fondo del valle y se queda corto en
#   las laderas altas, donde la lluvia orográfica no se ve a 5 km de resolución.
#   ENSO: en años Niña las estaciones reciben 1 926 mm y el satélite 1 731
#   (subestima 10 %); en años Niño 1 570 frente a 1 833 (sobreestima 17 %). El
#   satélite amortigua la señal ENSO. Esto justifica dar al modelo la altitud y
#   el ONI junto con la curva satelital, como pide el enunciado.

comparacion <- list(pares = pares, d = b, d_marcados = b_def, intervalos = int_d, R_YX = R_YX, Cov_YY = Cov_YY, Cov_XX = Cov_XX,
             grupos = grupos, sesgo_estacion = por_est)

# =============================================================================
# 5. REGRESIÓN FUNCIONAL FUNCIÓN A FUNCIÓN (SATÉLITE -> ESTACIÓN)  — rúbrica 4.2 y 4.4
# =============================================================================
# Modelo: y_i(t) = mu(t) + alfa(t) + integral beta(s,t) x_i(s) ds
#                  + gamma(t) altitud_i + delta(t) ONI_i + e_i(t)
# sobre la raíz de la lluvia (curvas en mm recortadas en 0 y luego con raíz, D-23).
# beta(s,t) = sum_m beta_m(s) psi_m(t): psi_m son las funciones propias del
# FPCA de las estaciones (solo entrenamiento) y beta_m(s) vive en la base
# B-spline del satélite con penalización de rugosidad lambda * integral beta_m''^2.
# Es el modelo de rango bajo al que equivale la cabeza lineal de la red
# funcional. M y lambda se eligen con validación (1996-2000) y se evalúa en
# prueba (2001-2003) frente al satélite sin corregir. Entran todos los años.
titulo("5. Regresión funcional")
set.seed(20261113)

K_e <- o$K_est; K_s <- o$K_sat
Phi_e <- base_bspline(T_DIA, K_e); Phi_s <- base_bspline(T_DIA, K_s)
W_e <- o$W_est; W_s <- o$W_sat
R_s <- gram_derivada(K_s, 2L)
ref_d2 <- splines::splineDesign(nodos(K_s), T_DIA, ord = ORDEN, derivs = rep(2, 365), outer.ok = TRUE)
cat("Segunda derivada de la base propia vs splines::splineDesign:", format(max(abs(derivada_base(T_DIA, K_s, 2L) - ref_d2)), digits = 2), "\n")

titulo("5.1 Escala raíz")
a_raiz <- function(C, Phi) {
  F <- C %*% t(Phi); G <- sqrt(pmax(F, 0)); Cr <- G %*% t(solve(crossprod(Phi), t(Phi)))
  list(coef = Cr, pct_recortado = 100 * mean(F < 0), rmse_reproyeccion = sqrt(mean((Cr %*% t(Phi) - G)^2)))
}
re <- a_raiz(o$coef_est[pares$fila_est, ], Phi_e); rs <- a_raiz(o$coef_sat[pares$fila_sat, ], Phi_s)
print(data.table(fuente = c("estaciones", "satélite"), pct_dias_recortados_en_0 = round(c(re$pct_recortado, rs$pct_recortado), 2),
                 rmse_reproyeccion_raiz = signif(c(re$rmse_reproyeccion, rs$rmse_reproyeccion), 3)))

# ¿Suavizar en mm y luego sacar la raíz, o sacar la raíz al dato diario y luego suavizar?
# Se comprueba cuál conserva el acumulado anual (métrica de la rúbrica 4.4), con las
# curvas completas de entrenamiento y el mismo K.
i_cmp <- which(id_est$n_obs == 365L & id_est$particion == "entrenamiento")
mapa_e <- solve(crossprod(Phi_e), t(Phi_e)); obs_mm <- rowMeans(Y_dia[i_cmp, ])
camino_A <- pmax(Y_dia[i_cmp, ] %*% t(mapa_e) %*% t(Phi_e), 0)
camino_B <- (sqrt(Y_dia[i_cmp, ]) %*% t(mapa_e) %*% t(Phi_e))^2
orden_raiz <- data.table(camino = c("A. Suavizar en mm, recortar en 0 y luego raíz (el que se usa)", "B. Raíz al dato diario, suavizar y elevar al cuadrado"),
                         curvas = length(i_cmp), cociente_acumulado_anual_mediano = round(c(median(rowMeans(camino_A) / obs_mm), median(rowMeans(camino_B) / obs_mm)), 3),
                         cociente_p05 = round(c(quantile(rowMeans(camino_A) / obs_mm, .05), quantile(rowMeans(camino_B) / obs_mm, .05)), 3),
                         cociente_p95 = round(c(quantile(rowMeans(camino_A) / obs_mm, .95), quantile(rowMeans(camino_B) / obs_mm, .95)), 3))
print(tabla_f(orden_raiz, "F5_orden_raiz_y_suavizado"))
stopifnot(abs(orden_raiz$cociente_acumulado_anual_mediano[1] - 1) < .02)
Y_mm <- pmax(o$coef_est[pares$fila_est, ] %*% t(Phi_e), 0)
X_mm <- pmax(o$coef_sat[pares$fila_sat, ] %*% t(Phi_s), 0)

pares[, usar := !is.na(altitud)]   # todos los años entran (D-47)
cat("Pares fuera del modelo por no tener altitud:", pares[is.na(altitud), .N], "\n")
print(pares[usar == TRUE, .(pares = .N, estaciones = uniqueN(codigo), anios = paste(range(anio), collapse = "-")), by = particion])
ent <- pares$usar & pares$particion == "entrenamiento"
val <- pares$usar & pares$particion == "validación"
pru <- pares$usar & pares$particion == "prueba"

titulo("5.2 FPCA de las estaciones (solo entrenamiento)")
ev <- eigen(W_e, symmetric = TRUE); W12 <- ev$vectors %*% diag(sqrt(ev$values)) %*% t(ev$vectors); Wm12 <- ev$vectors %*% diag(1 / sqrt(ev$values)) %*% t(ev$vectors)
mu_c <- colMeans(re$coef[ent, ])
Cc <- sweep(re$coef[ent, ], 2, mu_c)
eg <- eigen(W12 %*% (crossprod(Cc) / (sum(ent) - 1)) %*% W12, symmetric = TRUE)
psi_c <- Wm12 %*% eg$vectors
var_fpca <- data.table(m = seq_len(K_e), autovalor = eg$values, pct = 100 * eg$values / sum(eg$values))[, acumulado := cumsum(pct)]
print(tabla_f(var_fpca, "F5_fpca_varianza"))
stopifnot(max(abs(t(psi_c) %*% W_e %*% psi_c - diag(K_e))) < 1e-10)
fd_r <- fda::fd(t(re$coef[ent, ]), fda::create.bspline.basis(RANGO, nbasis = K_e, norder = ORDEN))
pca_fda <- fda::pca.fd(fd_r, nharm = 4, centerfns = TRUE)
psi_fda <- fda::eval.fd(T_DIA, pca_fda$harmonics)
cat(sprintf("Contra fda::pca.fd (4 primeras): autovalores x N/(N-1) a %.1e relativo, funciones propias a %.1e (salvo signo); fda integra W numéricamente\n",
            max(abs(pca_fda$values[1:4] * sum(ent) / (sum(ent) - 1) - eg$values[1:4]) / eg$values[1:4]),
            max(abs(abs(psi_fda) - abs(Phi_e %*% psi_c[, 1:4])))))
S_todos <- sweep(re$coef, 2, mu_c) %*% W_e %*% psi_c

titulo("5.3 Ajuste y elección de M y lambda con validación")
alt_m <- mean(pares$altitud[ent]); alt_s <- sd(pares$altitud[ent])
Zx <- rs$coef %*% W_s
diseno <- function(i, con_sat = TRUE) cbind(1, if (con_sat) Zx[i, , drop = FALSE], (pares$altitud[i] - alt_m) / alt_s, pares$ONI_anual[i])
penal <- function(con_sat = TRUE) { P <- matrix(0, 3 + con_sat * K_s, 3 + con_sat * K_s); if (con_sat) P[1 + 1:K_s, 1 + 1:K_s] <- R_s; P }
ajustar_B <- function(lambda, M, con_sat = TRUE) { Z <- diseno(which(ent), con_sat); solve(crossprod(Z) + lambda * penal(con_sat), crossprod(Z, S_todos[ent, 1:M, drop = FALSE])) }
predecir_raiz <- function(B, i, con_sat = TRUE) { M <- ncol(B); sweep((diseno(i, con_sat) %*% B) %*% t(Phi_e %*% psi_c[, 1:M, drop = FALSE]), 2, drop(Phi_e %*% mu_c), "+") }
a_mm <- function(R) pmax(R, 0)^2
rmse_l2 <- function(P, i) sqrt(rowMeans((P - Y_mm[i, ])^2))

lambdas <- 10^seq(-2, 10, by = .5)
sel_lambda <- rbindlist(lapply(seq_len(K_e), function(M) rbindlist(lapply(lambdas, function(l) {
  B <- ajustar_B(l, M)
  data.table(M = M, lambda = l, rmse_entrenamiento = mean(rmse_l2(a_mm(predecir_raiz(B, which(ent))), which(ent))),
             rmse_validacion = mean(rmse_l2(a_mm(predecir_raiz(B, which(val))), which(val))))
}))))
tabla_f(sel_lambda, "F5_seleccion_M_lambda")
mejor_por_M <- sel_lambda[, .SD[which.min(rmse_validacion)], by = M]
print(mejor_por_M)
minimo <- min(mejor_por_M$rmse_validacion)
M_sel <- mejor_por_M[rmse_validacion <= minimo * 1.005, min(M)]
LAMBDA <- mejor_por_M[M == M_sel, lambda]; M <- M_sel
psi_c <- psi_c[, 1:M, drop = FALSE]; psi_t <- Phi_e %*% psi_c; S <- S_todos[, 1:M, drop = FALSE]
cat(sprintf("Elegidos: M = %d (%.1f %% de la varianza; el menor M a menos de 0,5 %% del mejor RMSE de validación) y lambda = %s\n",
            M, var_fpca$acumulado[M], format(LAMBDA, scientific = TRUE)))
grafico_f(ggplot(sel_lambda, aes(lambda, rmse_validacion, colour = factor(M))) + geom_line() +
            geom_point(data = mejor_por_M, size = 2) + geom_vline(xintercept = LAMBDA, linetype = 2, colour = "red") + scale_x_log10() +
            labs(x = expression(lambda), y = "RMSE L2 medio en validación (mm/día)", colour = "M",
                 title = "Elección de M y lambda con los años de validación (1996-1999)"),
          "F5_seleccion_M_lambda", ancho = 9, alto = 5)

B <- ajustar_B(LAMBDA, M); B_sin_sat <- ajustar_B(0, M, con_sat = FALSE)
gl_df <- function(Z, P, l) sum(diag(solve(crossprod(Z) + l * P, crossprod(Z))))
cat(sprintf("Grados de libertad efectivos por componente: %.1f de %d parámetros\n", gl_df(diseno(which(ent)), penal(), LAMBDA), ncol(diseno(1))))

titulo("5.4 Superficie beta(s,t) y verificación de la equivalencia")
beta_m <- B[1 + 1:K_s, , drop = FALSE]
BETA <- Phi_s %*% beta_m %*% t(psi_t)
alfa_t <- drop(psi_t %*% B[1, ]); gamma_t <- drop(psi_t %*% B[K_s + 2, ]); delta_t <- drop(psi_t %*% B[K_s + 3, ])
tau_s <- unique(nodos(K_s)); a <- head(tau_s, -1); b <- tail(tau_s, -1)
sq <- as.vector(outer(GL_x, seq_along(a), function(x, k) (b[k] - a[k]) / 2 * x + (a[k] + b[k]) / 2))
wq <- as.vector(outer(GL_w, seq_along(a), function(w, k) (b[k] - a[k]) / 2 * w))
BETA_q <- base_bspline(sq, K_s) %*% beta_m %*% t(psi_t)
i_pru <- which(pru)
x_q <- rs$coef[i_pru, ] %*% t(base_bspline(sq, K_s))
pred_superficie <- sweep(x_q %*% (wq * BETA_q), 2, drop(Phi_e %*% mu_c) + alfa_t, "+") +
  outer((pares$altitud[i_pru] - alt_m) / alt_s, gamma_t) + outer(pares$ONI_anual[i_pru], delta_t)
pred_puntajes <- predecir_raiz(B, i_pru)
dif_equiv <- max(abs(pred_superficie - pred_puntajes))
cat(sprintf("Predicción con la integral de beta(s,t) (cuadratura) vs predicción por puntajes FPCA: diferencia máxima %.2e (criterio <= 1e-4)\n", dif_equiv))
stopifnot(dif_equiv <= 1e-4)

grafico_base("F5_fpca", function() {
  par(mfrow = c(1, 2))
  plot(T_DIA, drop(Phi_e %*% mu_c), type = "l", lwd = 2, xlab = "Día del año", ylab = "raíz(mm/día)", main = "Media de entrenamiento (escala raíz)")
  matplot(T_DIA, psi_t, type = "l", lty = 1, lwd = 2, col = 1:M, xlab = "Día del año", ylab = expression(psi[m](t)),
          main = paste("Funciones propias:", paste(sprintf("%.0f %%", var_fpca$pct[1:M]), collapse = ", ")))
  abline(h = 0, lty = 3)
}, ancho = 12, alto = 4.5)
grafico_base("F5_beta_superficie", function() {
  par(mfrow = c(1, 2))
  sem <- seq(1, 365, by = 7)
  persp(sem, sem, BETA[sem, sem], theta = 35, phi = 30, col = "lightblue", border = "grey50", xlab = "s (satélite)", ylab = "t (estación)",
        zlab = "beta(s,t)", main = expression(hat(beta)(s, t)), ticktype = "simple")
  lim <- max(abs(BETA)) * c(-1, 1)
  image(T_DIA, T_DIA, BETA, zlim = lim, col = hcl.colors(64, "Blue-Red 3"), xlab = "s: día del satélite", ylab = "t: día de la estación", main = "Curvas de nivel")
  contour(T_DIA, T_DIA, BETA, add = TRUE, nlevels = 10, drawlabels = FALSE, col = "grey30"); abline(0, 1, lty = 2)
}, ancho = 13, alto = 5.5)
grafico_base("F5_efectos", function() {
  par(mfrow = c(1, 3))
  plot(T_DIA, colSums(BETA), type = "l", lwd = 2, xlab = "Día t de la estación", ylab = "", main = expression(integral(beta(s, t) * ds)))
  abline(h = 0, lty = 3)
  plot(T_DIA, gamma_t, type = "l", lwd = 2, col = "brown", xlab = "Día del año", ylab = "raíz(mm/día) por 1 DE", main = "Efecto de la altitud gamma(t)")
  abline(h = 0, lty = 3)
  plot(T_DIA, delta_t, type = "l", lwd = 2, col = 4, xlab = "Día del año", ylab = "raíz(mm/día) por unidad de ONI", main = "Efecto del ONI delta(t)")
  abline(h = 0, lty = 3)
}, ancho = 14, alto = 4.5)
cat(sprintf("Integral de beta sobre s: máxima el %s, mínima el %s | efecto ONI más negativo el %s (%.3f) | efecto altitud medio %.3f\n",
            fecha_dia(which.max(colSums(BETA))), fecha_dia(which.min(colSums(BETA))), fecha_dia(which.min(delta_t)), min(delta_t), mean(gamma_t)))

titulo("5.5 Evaluación en prueba (2001-2003) con las métricas de la rúbrica")
# Las cuatro métricas de la rúbrica 4.4, por curva:
#   - error L² relativo  = ||ŷ - y|| / ||y||   (y además el error integrado L² en mm/día)
#   - KGE
#   - sesgo del acumulado anual en %
#   - error en la fecha de los dos picos (días)
# IC del 95 % por bootstrap de ESTACIONES (2 000 remuestreos): las curvas de una
# misma estación no son independientes. Prueba de Wilcoxon pareada por estación.
kge <- function(P, Yo) {
  r <- sapply(seq_len(nrow(P)), function(k) cor(P[k, ], Yo[k, ]))
  1 - sqrt((r - 1)^2 + (apply(P, 1, sd) / apply(Yo, 1, sd) - 1)^2 + (rowMeans(P) / rowMeans(Yo) - 1)^2)
}
pico <- function(M, ventana) apply(M[, ventana, drop = FALSE], 1, which.max) + min(ventana) - 1
metricas_de <- function(P, i, modelo) {
  Yo <- Y_mm[i, , drop = FALSE]
  data.table(modelo = modelo, codigo = pares$codigo[i], anio = pares$anio[i], fase = pares$fase[i], particion = pares$particion[i],
             l2_rel = sqrt(rowSums((P - Yo)^2) / rowSums(Yo^2)), rmse = sqrt(rowMeans((P - Yo)^2)), kge = kge(P, Yo),
             sesgo_pct = 100 * (rowMeans(P) - rowMeans(Yo)) / rowMeans(Yo), sesgo_anual = (rowMeans(P) - rowMeans(Yo)) * 365,
             err_pico_1 = abs(pico(P, 60:200) - pico(Yo, 60:200)), err_pico_2 = abs(pico(P, 230:350) - pico(Yo, 230:350)))
}
N_BOOT <- 2000L
boot_ic <- function(v, g) { m <- tapply(v, g, mean); c(mean(m), quantile(replicate(N_BOOT, mean(sample(m, replace = TRUE))), c(.025, .975))) }
boot_media <- function(v, g, dec = 2) { b <- boot_ic(v, g); sprintf(paste0("%.", dec, "f [%.", dec, "f, %.", dec, "f]"), b[1], b[2], b[3]) }
resumir <- function(mc) mc[, .(curvas = .N, estaciones = uniqueN(codigo), L2_relativo = boot_media(l2_rel, codigo), RMSE_L2_mm_dia = boot_media(rmse, codigo),
                               KGE = boot_media(kge, codigo), sesgo_anual_pct = boot_media(sesgo_pct, codigo, 1),
                               error_pico_1_dias = boot_media(err_pico_1, codigo, 1), error_pico_2_dias = boot_media(err_pico_2, codigo, 1))]

modelos <- list(`Satélite sin corregir` = X_mm[i_pru, ], `Regresión funcional` = a_mm(pred_puntajes),
                `Solo altitud y ONI` = a_mm(predecir_raiz(B_sin_sat, i_pru, FALSE)))
recorte_pred <- 100 * mean(pred_puntajes < 0)
metricas_curva <- rbindlist(lapply(names(modelos), function(nm) metricas_de(modelos[[nm]], i_pru, nm)))
metricas <- metricas_curva[, resumir(.SD), by = modelo]
print(tabla_f(metricas, "F5_metricas_prueba"))
tabla_f(metricas_curva, "F5_metricas_prueba_por_curva")
cat(sprintf("Reconstrucción en mm no negativa: %.2f %% de los valores predichos en la escala raíz se recortan en 0\n", recorte_pred))

# ¿La regresión supera al satélite sin corregir? Diferencia pareada por estación.
ancho <- dcast(metricas_curva[modelo != "Solo altitud y ONI"], codigo + anio ~ modelo, value.var = c("l2_rel", "rmse", "kge"))
por_estacion <- metricas_curva[, lapply(.SD, mean), by = .(modelo, codigo), .SDcols = c("l2_rel", "rmse", "kge", "sesgo_pct", "err_pico_1", "err_pico_2")]
pe <- dcast(por_estacion, codigo ~ modelo, value.var = c("l2_rel", "rmse", "kge", "sesgo_pct", "err_pico_1", "err_pico_2"))
comparar <- function(v, mejor_si_baja = TRUE, absoluto = FALSE) {
  a <- pe[[paste0(v, "_Regresión funcional")]]; s <- pe[[paste0(v, "_Satélite sin corregir")]]
  if (absoluto) { a <- abs(a); s <- abs(s) }
  d <- a - s; q <- quantile(replicate(N_BOOT, mean(sample(d, replace = TRUE))), c(.025, .975))
  data.table(metrica = v, satelite = mean(s), regresion = mean(a), diferencia = mean(d), ic_inf = q[1], ic_sup = q[2],
             p_wilcoxon = wilcox.test(a, s, paired = TRUE)$p.value, ic_excluye_0 = q[1] > 0 | q[2] < 0,
             mejora = if (mejor_si_baja) mean(d) < 0 else mean(d) > 0)
}
pareada <- rbind(comparar("l2_rel"), comparar("rmse"), comparar("kge", FALSE), comparar("sesgo_pct", absoluto = TRUE), comparar("err_pico_1"), comparar("err_pico_2"))
pareada[metrica == "sesgo_pct", metrica := "|sesgo_pct|"]
pareada[, c("satelite", "regresion", "diferencia", "ic_inf", "ic_sup") := lapply(.SD, round, 3), .SDcols = c("satelite", "regresion", "diferencia", "ic_inf", "ic_sup")]
pareada[, p_wilcoxon := signif(p_wilcoxon, 2)]
cat("Regresión frente al satélite sin corregir, pareado por estación (", nrow(pe), "estaciones; IC 95 % bootstrap de la diferencia y Wilcoxon):\n")
print(tabla_f(pareada, "F5_regresion_vs_satelite_pareado"))
cat(sprintf("La regresión baja el error L² relativo frente al satélite sin corregir en %.0f %% de las curvas de prueba\n",
            100 * mean(ancho$`l2_rel_Regresión funcional` < ancho$`l2_rel_Satélite sin corregir`)))

# Por fase ENSO. La prueba (2001-2003) no tiene ningún año Niña; por eso se agrega
# la validación (1996-2000), que sí tiene Niña y Niño. Esa parte es ORIENTATIVA:
# la validación se usó para elegir M y lambda (dos números), no para ajustar B.
i_val <- which(val)
mc_val <- rbind(metricas_de(X_mm[i_val, ], i_val, "Satélite sin corregir"), metricas_de(a_mm(predecir_raiz(B, i_val)), i_val, "Regresión funcional"))
por_fase <- rbind(metricas_curva[modelo != "Solo altitud y ONI"], mc_val)[, resumir(.SD), by = .(particion, fase, modelo)]
por_fase[, anios := sapply(seq_len(.N), function(k) paste(sort(unique(pares[particion == por_fase$particion[k] & fase == por_fase$fase[k] & usar, anio])), collapse = " "))]
setcolorder(por_fase, c("particion", "fase", "anios")); setorder(por_fase, -particion, fase, modelo)
print(tabla_f(por_fase, "F5_metricas_por_fase_ENSO"))
por_anio_pru <- metricas_curva[, .(curvas = .N, L2_relativo = round(mean(l2_rel), 3), KGE = round(mean(kge), 2), sesgo_anual_pct = round(mean(sesgo_pct), 1)), by = .(anio, fase, modelo)]
print(tabla_f(por_anio_pru[order(anio, modelo)], "F5_metricas_prueba_por_anio"))

titulo("5.6 Sensibilidad: ¿qué pasa si se quitan 1985, 1988 y 2000?")
# No se quitan (D-47). Esta sección solo mide cuánto cambia la prueba si el modelo
# NO viera esos años: mismo FPCA y mismo M; lambda se vuelve a elegir en validación sin 2000.
ent_todos <- ent; val_todos <- val
ent <- ent_todos & !pares$marcado; val <- val_todos & !pares$marcado
rmse_val_s <- sapply(lambdas, function(l) mean(rmse_l2(a_mm(predecir_raiz(ajustar_B(l, M), which(val))), which(val))))
LAMBDA_S <- lambdas[which.min(rmse_val_s)]; B_sens <- ajustar_B(LAMBDA_S, M)
ent <- ent_todos; val <- val_todos
mc_sens <- metricas_de(a_mm(predecir_raiz(B_sens, i_pru)), i_pru, "Regresión sin 1985, 1988 y 2000")
sens <- rbind(metricas_curva[modelo == "Regresión funcional"][, modelo := "Regresión con todos los años (la que se usa)"], mc_sens)
sensibilidad <- sens[, .(pares_entrenamiento = 0L, lambda = 0, L2_relativo = mean(l2_rel), RMSE_L2 = mean(rmse), KGE = mean(kge), sesgo_anual_pct = mean(sesgo_pct),
                         error_pico_1 = mean(err_pico_1), error_pico_2 = mean(err_pico_2)), by = modelo]
sensibilidad[, `:=`(pares_entrenamiento = c(sum(ent_todos), sum(ent_todos & !pares$marcado)), lambda = c(LAMBDA, LAMBDA_S))]
sens_anio <- dcast(sens[, .(L2 = round(mean(l2_rel), 3), sesgo = round(mean(sesgo_pct), 1)), by = .(modelo, anio)], anio ~ modelo, value.var = c("L2", "sesgo"))
ps <- merge(sens[modelo %like% "todos", .(a = mean(l2_rel)), by = codigo], mc_sens[, .(s = mean(l2_rel)), by = codigo], by = "codigo")
sensibilidad[, c("L2_relativo", "RMSE_L2", "KGE", "sesgo_anual_pct", "error_pico_1", "error_pico_2") := lapply(.SD, round, 3), .SDcols = c("L2_relativo", "RMSE_L2", "KGE", "sesgo_anual_pct", "error_pico_1", "error_pico_2")]
print(tabla_f(sensibilidad, "F5_sensibilidad_anios_marcados"))
print(sens_anio)
cat(sprintf("Con todos los años el error L² relativo en prueba es menor en %d de %d estaciones (Wilcoxon pareado p = %.3g)\n",
            sum(ps$a < ps$s), nrow(ps), wilcox.test(ps$a, ps$s, paired = TRUE)$p.value))

R2_t <- function(P, Yo) 1 - colSums((Yo - P)^2) / colSums(sweep(Yo, 2, colMeans(Yo))^2)
grafico_base("F5_prueba", function() {
  par(mfrow = c(2, 3))
  ej <- i_pru[order(metricas_curva[modelo == "Regresión funcional", rmse])[c(round(length(i_pru) * c(.1, .5, .9)))]]
  for (k in ej) {
    j <- match(k, i_pru)
    matplot(T_DIA, cbind(Y_mm[k, ], modelos[[1]][j, ], modelos[[2]][j, ]), type = "l", lty = c(1, 2, 1), lwd = 2, col = c(1, "darkorange", 4),
            xlab = "Día del año", ylab = "mm/día", main = paste(pares$nombre[k], pares$anio[k]))
    legend("topleft", c("Estación", "Satélite", "Regresión"), col = c(1, "darkorange", 4), lty = c(1, 2, 1), lwd = 2, bty = "n", cex = .8)
  }
  plot(T_DIA, R2_t(modelos[[2]], Y_mm[i_pru, ]), type = "l", lwd = 2, col = 4, ylim = c(-1, 1), xlab = "Día del año", ylab = expression(R^2*(t)),
       main = "Varianza explicada por día (prueba)")
  lines(T_DIA, R2_t(modelos[[1]], Y_mm[i_pru, ]), lwd = 2, col = "darkorange"); abline(h = 0, lty = 3)
  legend("bottomleft", c("Regresión", "Satélite sin corregir"), col = c(4, "darkorange"), lwd = 2, bty = "n")
  boxplot(l2_rel ~ modelo, metricas_curva, ylab = "Error L² relativo", xlab = "", main = "Error por curva (prueba)", cex.axis = .8)
  boxplot(sesgo_pct ~ modelo, metricas_curva, ylab = "% del acumulado anual", xlab = "", main = "Sesgo anual (prueba)", cex.axis = .8); abline(h = 0, lty = 3)
}, ancho = 15, alto = 8)

titulo("5.7 Exportación para la red funcional (Python)")
dir_py <- file.path(SALIDA_F, "para_python"); dir.create(dir_py, showWarnings = FALSE)
fwrite(pares[, .(codigo, anio, particion, usar, celda, altitud, ONI_anual, fase, vertiente)], file.path(dir_py, "pares.csv"))
fwrite(as.data.table(re$coef), file.path(dir_py, "coef_estacion_raiz.csv")); fwrite(as.data.table(rs$coef), file.path(dir_py, "coef_satelite_raiz.csv"))
fwrite(as.data.table(o$coef_est[pares$fila_est, ]), file.path(dir_py, "coef_estacion_mm.csv")); fwrite(as.data.table(o$coef_sat[pares$fila_sat, ]), file.path(dir_py, "coef_satelite_mm.csv"))
for (nm in c("W_e", "W_s", "R_s")) fwrite(as.data.table(get(nm)), file.path(dir_py, paste0(nm, ".csv")), col.names = FALSE)
fwrite(data.table(mu = mu_c, psi_c), file.path(dir_py, "fpca_estaciones.csv"))
fwrite(as.data.table(B), file.path(dir_py, "coeficientes_regresion.csv"), col.names = FALSE)
writeLines(c(sprintf('{"rango": [1, 365], "orden": %d, "K_estaciones": %d, "K_satelite": %d, "M": %d, "lambda": %g,', ORDEN, K_e, K_s, M, LAMBDA),
             sprintf('"nodos_estaciones": [%s], "nodos_satelite": [%s],', paste(nodos(K_e), collapse = ", "), paste(nodos(K_s), collapse = ", ")),
             sprintf('"altitud_media": %.6f, "altitud_de": %.6f, "escala": "raiz de la curva en mm/dia recortada en 0"}', alt_m, alt_s)),
           file.path(dir_py, "especificacion.json"))
cat("Archivos para Python en", dir_py, ":", paste(list.files(dir_py), collapse = ", "), "\n")

# QUÉ DICEN LOS DATOS
#   Escala: al recortar en 0 solo cambia el 2 % de los días de las estaciones y
#   casi ningún día del satélite; la raíz se reproyecta en la misma base con un
#   error de 0,09 (estaciones) y 0,03 (satélite) en raíz de mm/día.
#   Orden de la raíz y el suavizado: suavizar en mm y luego sacar la raíz conserva
#   el acumulado anual (cociente mediano 1,000 frente al dato diario). Sacar la
#   raíz al dato diario y luego suavizar lo reduce a un tercio (0,33), porque la
#   lluvia diaria es intermitente y el cuadrado de la media de raíces no es la
#   media. Como el sesgo del acumulado anual es una métrica de la rúbrica, se
#   usa el primer camino.
#   Datos del modelo: 3 299 pares con todos los años (2 257 de entrenamiento, 672
#   de validación, 370 de prueba). Univalle entra con 996 m (D-46).
#   FPCA (solo entrenamiento): la primera componente explica 62 % (el nivel de
#   lluvia de cada sitio), la segunda 13 % y la tercera 7 %. M y lambda se eligen
#   juntos con validación: el error deja de bajar desde M = 6 (94,2 % de la
#   varianza; el enunciado pide M entre 4 y 6), con lambda = 3,2e8 (10 grados de
#   libertad efectivos de 14).
#   Equivalencia: la predicción con la integral de beta(s,t) por cuadratura y la
#   predicción por puntajes FPCA difieren en 4e-15 (criterio 1e-4). Es la
#   identidad beta(s,t) = sum_m beta_m(s) psi_m(t) que debe reproducir la cabeza
#   lineal de la red.
#   Interpretación: lejos de los bordes, la integral de beta sobre s va de 0,62
#   (diciembre) a 1,03 (1-jun): el modelo traslada el nivel del satélite casi
#   completo a mitad de año y lo amortigua al final. El ONI resta lluvia todo el
#   año y sobre todo a mitad de año (-0,24 raíz(mm/día) por unidad el 2-jul): El Niño profundiza la sequía de
#   junio-julio. La altitud pesa poco una vez que se conoce el píxel.
#   Prueba (2001-2003, 375 curvas de 137 estaciones, IC 95 % por bootstrap de
#   estaciones con 2 000 remuestreos). Frente al satélite sin corregir, la
#   regresión baja el error L² relativo de 0,51 a 0,44, sube el KGE de 0,38 a
#   0,47, baja el sesgo del acumulado anual de +26 % a +7 % y acerca el pico de
#   la primera temporada de 23 a 18 días. Pareado por estación, el IC de la
#   diferencia excluye el cero y el Wilcoxon da p < 0,001 en cinco métricas; en
#   el pico de la segunda temporada la ganancia es de medio día y no es
#   significativa (IC con el cero, p = 0,08). Mejora en 65 % de las curvas. Sin
#   el satélite (solo altitud y ONI) el error queda como el del satélite crudo
#   (0,51) y lejos del 0,44 de la regresión: el satélite sí aporta.
#   Por fase ENSO: la prueba solo tiene años neutros (2001, 2003) y un Niño
#   (2002); no tiene Niña. En 2002 el modelo corrige de más (sesgo -14 % frente a
#   +28 % del satélite). En validación (orientativa) el Niño de 1997 es donde más
#   gana (L² 0,61 a 0,46). En los años Niña 1999-2000 casi empata con el
#   satélite (0,41 a 0,40) y el sesgo empeora (de -0,5 % a -11 %): el modelo baja
#   el nivel del satélite como en un año normal, y el término del ONI no alcanza
#   a devolver lo que el satélite ya perdió en La Niña. Es el punto débil de la
#   referencia lineal y lo que una cabeza no lineal podría mejorar.
#   Sensibilidad: si el modelo no viera 1985, 1988 y 2000, la prueba empeora (L²
#   0,439 frente a 0,434; sesgo 8,9 % frente a 7,4 %) en 79 de 137 estaciones
#   (Wilcoxon p = 0,003). Quitar esos años no ayuda: confirma con
#   datos la indicación del profesor (D-47).
#   La reconstrucción en mm es no negativa y solo 0,02 % de los valores se
#   recorta en 0.
#   Esta regresión es la referencia lineal. La red funcional (codificador con W,
#   decodificador FPCA con estas psi_m y cabezas lineal, KAN y MLP) sigue en el
#   archivo de Python con los insumos de para_python/.

regresion <- list(M = M, lambda = LAMBDA, B = B, mu_c = mu_c, psi_c = psi_c, var_fpca = var_fpca, BETA = BETA, metricas = metricas,
             metricas_curva = metricas_curva, dif_equivalencia = dif_equiv, alt_m = alt_m, alt_s = alt_s)

# =============================================================================
# 6. PRUEBAS Y CIERRE
# =============================================================================
titulo("6. Pruebas de la parte funcional")
pruebas <- list(
  "El emparejamiento es el mismo del archivo de datos" = identical(pares[order(fila_est), fila_sat], dp$pares[order(fila_est), fila_sat]) && identical(pares[order(fila_est), codigo], dp$pares[order(fila_est), codigo]),
  "Ningún año queda fuera del modelo" = all(pares$usar) && setequal(unique(pares$anio[ent]), 1981:1995) && setequal(unique(pares$anio[val]), 1996:2000) && setequal(unique(pares$anio[pru]), 2001:2003),
  "La base propia coincide con splines::splineDesign" = max(abs(base_bspline(T_DIA, 11) - ref_splines)) < 1e-12,
  "S simétrica, idempotente y tr(S) = K en las dos fuentes" = all(controles[, .(simetria, idempotencia, traza_menos_K)] < 1e-8),
  "W integra el dominio: 1'W1 = 364" = abs(sum(W_e) - 364) < 1e-8 && abs(sum(W_s) - 364) < 1e-8,
  "La media paso a paso coincide con fda::mean.fd" = all(verif_media$vs_fda_mean.fd < 1e-10),
  "Las atípicas coinciden con fda::fbplot" = all(verif_fb),
  "Funciones propias ortonormales con W" = max(abs(t(psi_c) %*% W_e %*% psi_c - diag(M))) < 1e-10,
  "K, FPCA y escalas solo con entrenamiento; M y lambda con validación" = !any(ent & pares$anio > 1995) && !any(val & (pares$anio < 1996 | pares$anio > 2000)),
  "Equivalencia beta(s,t) vs puntajes FPCA (<= 1e-4)" = dif_equiv <= 1e-4,
  "La reconstrucción en mm no es negativa" = min(modelos[["Regresión funcional"]]) >= 0,
  "Bootstrap con al menos 1000 remuestreos" = N_BOOT >= 1000
)
res_pruebas <- data.table(prueba = names(pruebas), pasa = unlist(pruebas))
print(tabla_f(res_pruebas, "F6_pruebas"))
cat(sprintf("\nPruebas: %d de %d pasan\n", sum(res_pruebas$pasa), nrow(res_pruebas)))
stopifnot(res_pruebas$pasa)
saveRDS(list(construccion = o, eda = eda, comparacion = comparacion, regresion = regresion, pruebas = res_pruebas),
        file.path(SALIDA_F, "resultados_funcional.rds"), compress = "xz")
writeLines(capture.output(sessionInfo()), file.path(SALIDA_F, "sessionInfo.txt"))

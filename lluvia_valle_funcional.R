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
#   así las 30 pruebas del archivo de datos cubren todo lo que entra.
# Salidas: salidas_R_funcional/tablas, /graficos y /para_python (insumos de la red).
#
# Cómo se corre (desde la raíz del repositorio):
#   Rscript lluvia_valle_datos.R        # una vez, deja datos_preparados.rds
#   Rscript lluvia_valle_funcional.R    # este archivo, ~1 minuto (fija UTF-8 solo)
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

# Locale: si R arranca sin UTF-8 (en la Mac del grupo LANG está vacío), las tildes del
# código ("validación") no coinciden con las de los datos y el script se detiene en la
# primera comparación de texto. Se fija UTF-8 aquí, antes de cualquier cadena con tilde.
if (!isTRUE(l10n_info()$`UTF-8`)) {
  for (loc in c("en_US.UTF-8", "es_CO.UTF-8", "es_ES.UTF-8", "C.UTF-8"))
    if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break
  if (!isTRUE(l10n_info()$`UTF-8`)) stop("No hay un locale UTF-8 disponible; corra con LC_ALL=en_US.UTF-8 Rscript ...")
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(fda)        # herramientas de clase: verificación (media, fbplot) y bases del contraste 5.5b
  library(fda.usc)    # herramientas de clase: fdata y fregre.basis.fr (contraste 5.5b)
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
COB_MIN_R5 <- 0.90
HUECO_MAX_R5 <- 15L

# Conjunto de datos (D-61, D-63): principal (R6 con ≥10 años válidos) o r6_omm (sensibilidad,
# R6 de la OMM, 19 de 23 años). Se elige con la variable de entorno CONJUNTO, igual que en la parte de datos.
CONJUNTO_F <- Sys.getenv("CONJUNTO", "principal")
if (!CONJUNTO_F %in% c("principal", "r6_omm")) stop("CONJUNTO debe ser 'principal' o 'r6_omm'")
SUFIJO_F <- if (CONJUNTO_F == "principal") "" else paste0("_", CONJUNTO_F)
SALIDA_F <- paste0("salidas_R_funcional", SUFIJO_F)
SALIDA_PY <- paste0("salidas_python", SUFIJO_F)
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

dp  <- readRDS(file.path(paste0("salidas_R_datos", SUFIJO_F), "datos_preparados.rds"))
stopifnot(identical(dp$conjunto, CONJUNTO_F))
cat("Conjunto:", CONJUNTO_F, "| curvas de estación:", nrow(dp$estaciones$Y), "| días rellenados por kriging (R7):", sum(dp$estaciones$rellenado), "\n")
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
# La zona es geografía, no calidad de datos: se usan TODAS las estaciones del Valle con
# coordenada oficial (232, no solo las 98 que entran al modelo). Si se usaran solo las
# incluidas, con 7 estaciones pacíficas muchos píxeles de la costa quedarían como "Cauca".
haversine <- function(la1, lo1, la2, lo2) 12742 * asin(sqrt(sin((la2 - la1) * pi / 360)^2 + cos(la1 * pi / 180) * cos(la2 * pi / 180) * sin((lo2 - lo1) * pi / 360)^2))
est_zona <- fread("datos_web/estaciones_valle_finales.csv", colClasses = list(character = "codigo"))[!is.na(lat) & !is.na(lon)]
est_zona[, vertiente := fifelse(substr(codigo, 1, 2) == "26", "Cauca (interandina)", "Pacífico")]
stopifnot(all(est$codigo %in% est_zona$codigo), all(est_zona[match(est$codigo, codigo), vertiente] == est$vertiente))
D_ee <- outer(seq_len(nrow(est_zona)), seq_len(nrow(est_zona)), function(i, j) haversine(est_zona$lat[i], est_zona$lon[i], est_zona$lat[j], est_zona$lon[j])); diag(D_ee) <- Inf
cortes <- seq(-77.1, -76.4, by = .05)
acierto_corte <- sapply(cortes, function(k) mean((est_zona$lon < k) == (est_zona$vertiente == "Pacífico")))
acierto_vecina <- mean(est_zona$vertiente[apply(D_ee, 1, which.min)] == est_zona$vertiente)
regla_vertiente <- data.table(regla = c(paste0("corte en longitud ", cortes[which.max(acierto_corte)], "° (el mejor corte)"), "corte en longitud -76.8° (el que se usaba)", "estación más cercana (dejando una fuera)"),
                              pct_acierto = round(100 * c(max(acierto_corte), acierto_corte[abs(cortes + 76.8) < 1e-9], acierto_vecina), 1))
print(tabla_f(regla_vertiente, "F1_regla_vertiente_pixeles"))
D_pe <- outer(seq_len(nrow(px)), seq_len(nrow(est_zona)), function(i, j) haversine(px$lat[i], px$lon[i], est_zona$lat[j], est_zona$lon[j]))
px[, `:=`(vertiente = est_zona$vertiente[apply(D_pe, 1, which.min)], km_estacion_cercana = round(apply(D_pe, 1, min), 1))]
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
#   Las cifras y su interpretación están en EXPLICACION_FUNCIONAL.md, sección 1 (corrida del
#   19-nov-2026, con las 98 estaciones de la regla de la OMM); las tablas de esta
#   sección las imprimen las líneas de arriba.


# =============================================================================
# 2. CONSTRUCCIÓN DE DATOS FUNCIONALES (ESTACIONES Y SATÉLITE)  — rúbrica 1
# =============================================================================
# Base B-spline cúbica (Cox-de Boor) y mínimos cuadrados programados sin
# librerías; K por ECM, VC y VCG solo con entrenamiento (1981-1995, todos los
# años); verificación de S; covarianza residual; bandas punto a punto; matriz
# de Gram W. Ambas fuentes quedan en el mismo dominio 1-365.
# Días faltantes: cada curva se ajusta solo con sus días observados (no se imputa).
titulo("2. Construcción de datos funcionales")

# "Completa" = los 365 días OBSERVADOS: los días que rellenó el kriging (R7, D-61) no cuentan
completas_est <- rowSums(is.na(Y)) == 0 & rowSums(dp$estaciones$rellenado) == 0
stopifnot(all(completas_est == (id_est$n_obs == 365L)))

grafico_base("F2_datos_discretos", function() {
  par(mfrow = c(1, 2))
  ej_cod <- names(sort(-table(id_est$codigo[completas_est])))[1]   # la estación con más años completos
  i <- which(completas_est & id_est$codigo == ej_cod)[1]
  plot(T_DIA, Y[i, ], type = "h", col = "grey40", xlab = "Día del año", ylab = "mm/día",
       main = paste("Estación", id_est$codigo[i], id_est$anio[i], "(diario)"))
  j <- which(id_sat$celda == emparejamiento[codigo == ej_cod, celda] & id_sat$anio == id_est$anio[i])
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

# P-20 (D-51): con K mayor la fecha del pico se acerca a una media móvil de 31 días de
# los datos crudos (prueba del stash "fase A"). ¿Es una mejora real o depende de la
# referencia? Se mide contra tres medias móviles (15, 31 y 61 días; ninguna depende de K)
# en las curvas completas de la vertiente interandina (régimen bimodal), solo de
# entrenamiento. Si el K que gana cambia con la referencia, la fecha del pico no sirve
# para elegir K y se queda el de la VCG.
vert_curva <- est$vertiente[match(id_est$codigo, est$codigo)]
cur_p20 <- ent_est & vert_curva != "Pacífico"
media_movil <- function(M, w) t(apply(M, 1, function(y) as.numeric(stats::filter(y, rep(1 / w, w), sides = 2))))
pico_en <- function(M, v) apply(M[, v, drop = FALSE], 1, which.max) + min(v) - 1
V1 <- 60:200; V2 <- 230:350; ANCHOS_REF <- c(15L, 31L, 61L)
picos_ref <- lapply(ANCHOS_REF, function(w) { R <- media_movil(Y[cur_p20, ], w); cbind(pico_en(R, V1), pico_en(R, V2)) })
K_picos <- rbindlist(lapply(c(K_EST, 10L, 12L, 16L, 20L), function(k) {
  P <- base_bspline(T_DIA, k); Fk <- Y[cur_p20, ] %*% t(P %*% solve(crossprod(P), t(P)))
  pk <- cbind(pico_en(Fk, V1), pico_en(Fk, V2))
  rbindlist(lapply(seq_along(ANCHOS_REF), function(r) data.table(K = k, referencia_dias = ANCHOS_REF[r],
    err_pico_1 = median(abs(pk[, 1] - picos_ref[[r]][, 1])), err_pico_2 = median(abs(pk[, 2] - picos_ref[[r]][, 2])))))
}))
K_picos[, err_total := err_pico_1 + err_pico_2]
gana_por_ref <- K_picos[, .(K_que_gana = K[which.min(err_total)]), by = referencia_dias]
dif_refs <- rbindlist(lapply(list(c(1, 2), c(1, 3), c(2, 3)), function(ab) data.table(
  referencias = paste(ANCHOS_REF[ab], collapse = " vs "),
  pico_1 = median(abs(picos_ref[[ab[1]]][, 1] - picos_ref[[ab[2]]][, 1])), pico_2 = median(abs(picos_ref[[ab[1]]][, 2] - picos_ref[[ab[2]]][, 2])))))
print(tabla_f(dcast(K_picos, K ~ referencia_dias, value.var = c("err_pico_1", "err_pico_2")), "F2_K_picos_vs_referencia"))
print(tabla_f(gana_por_ref, "F2_K_que_gana_por_referencia"))
print(tabla_f(dif_refs, "F2_diferencia_entre_referencias"))
cat(sprintf("Curvas: %d. Dispersión de la fecha del pico entre curvas (RIC, ref. 31 días): %d y %d días\n",
            sum(cur_p20), as.integer(IQR(picos_ref[[2]][, 1])), as.integer(IQR(picos_ref[[2]][, 2]))))
# La decisión depende de este resultado: si una misma K ganara con las tres referencias,
# el argumento de D-51 no valdría y el script se detiene para revisarlo.
stopifnot(uniqueN(gana_por_ref$K_que_gana) > 1)

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
# Curva con el hueco más largo que rellenó el kriging (R7): se ajusta SIN esos días para ver
# cuánto se abriría la banda si no se hubieran rellenado
Rell <- dp$estaciones$rellenado
inc <- which(rowSums(Rell) > 0); ig <- inc[which.max(sapply(inc, function(k) racha(Rell[k, ])))]
y_ig <- replace(Y[ig, ], Rell[ig, ], NA)
ag <- ajuste_mco(y_ig, Phi_e)
ee_g <- ee_banda(Phi_e, ag$mapa, ag$sigma2 * diag(sum(ag$obs)))
ee_g_completa <- ee_banda(Phi_e, e$mapa, ag$sigma2 * diag(365))
js <- which(id_sat$celda == emparejamiento[codigo == id_est$codigo[i_ej], celda] & id_sat$anio == id_est$anio[i_ej])
ee_s <- ee_banda(Phi_s365, s$mapa, s$sigma2[js] * diag(72))
cat(sprintf("Curva %s con %d días rellenados por kriging: sin rellenarlos, el EE en esos días sería %.2f veces el de la curva completa\n",
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
  plot(T_DIA, y_ig, pch = 16, cex = .4, col = "grey50", ylim = c(-5, max(Y[ig, ], na.rm = TRUE)), xlab = "Día del año",
       ylab = "mm/día", main = paste("Estación", id_est$curva_id[ig], "-", sum(!ag$obs), "días rellenados (banda sin ellos)"))
  points(T_DIA[!ag$obs], Y[ig, !ag$obs], pch = 4, col = 2)
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
#   Las cifras y su interpretación están en EXPLICACION_FUNCIONAL.md, sección 2 (corrida del
#   19-nov-2026, con las 98 estaciones de la regla de la OMM); las tablas de esta
#   sección las imprimen las líneas de arriba.

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
  par(mfrow = c(1, 4)); set.seed(1)
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
  # Coeficiente de variación funcional CV(t) = DE(t) / media(t) (clase, semanas 12-13): compara la
  # dispersión de fuentes con escalas distintas.
  plot(T_DIA, me$de / media_e, type = "l", lwd = 2, col = 4, ylim = range(c(0, me$de / media_e, ms$de / media_s)), xlab = "Día del año", ylab = "CV(t)",
       main = "Coeficiente de variación funcional")
  lines(T_DIA, ms$de / media_s, lwd = 2, col = "darkorange")
}, ancho = 19, alto = 4.8)
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

# Curvas por año (datos funcionales, en mm/día): un panel por año con las
# curvas en gris, el promedio en negro y la desviación estándar punteada. Todos los paneles
# tienen el mismo eje, para que los años se comparen a simple vista; el título lleva la fase ENSO.
por_anio <- function(M, anios, fases, cada, titulo_fig) {
  lim <- c(0, quantile(M, .995))
  par(mfrow = c(6, 4), mar = c(3, 3.2, 2, .6), mgp = c(1.9, .5, 0), oma = c(0, 0, 3, 0))
  for (a in ANIOS) {
    k <- anios == a
    matplot(T_DIA, t(M[k, , drop = FALSE]), type = "l", lty = 1, col = "grey80", ylim = lim, xaxt = "n", xlab = "Día del año",
            ylab = "Precipitación (mm/día)", main = sprintf("%d · %s", a, fases[which(k)[1]]), cex.axis = .75, cex.lab = .8)
    axis(1, at = c(1, 91, 182, 274, 365), labels = c("ene", "abr", "jul", "oct", "dic"), cex.axis = .75)
    lines(T_DIA, colMeans(M[k, , drop = FALSE]), lwd = 2); lines(T_DIA, apply(M[k, , drop = FALSE], 2, sd), lwd = 1.5, lty = 2)
  }
  plot.new(); legend("center", c(cada, "Promedio", "Desviación estándar"), col = c("grey70", "black", "black"),
                     lty = c(1, 1, 2), lwd = c(1, 2, 1.5), bty = "n", cex = 1)
  mtext(titulo_fig, outer = TRUE, cex = 1.1, font = 2)
}
grafico_base("F3_por_anio_estaciones", function() por_anio(Ye, id_e$anio, id_e$fase, "Cada estación",
             "Datos funcionales de las estaciones: promedio y desviación estándar por año, 1981–2003"), ancho = 11, alto = 15)
grafico_base("F3_por_anio_satelite", function() por_anio(Xs, id_s$anio, id_s$fase, "Cada píxel del Valle",
             "Datos funcionales del satélite: promedio y desviación estándar por año, 1981–2003"), ancho = 11, alto = 15)
resumen_anio <- oni[data.table(anio = ANIOS, estaciones = sapply(ANIOS, function(a) sum(id_e$anio == a)),
                               P_estaciones = sapply(ANIOS, function(a) round(mean(Ye[id_e$anio == a, ]) * 365)),
                               P_satelite = sapply(ANIOS, function(a) round(mean(Xs[id_s$anio == a, ]) * 365)),
                               DE_estaciones = sapply(ANIOS, function(a) round(mean(apply(Ye[id_e$anio == a, , drop = FALSE], 2, sd)), 2)),
                               DE_satelite = sapply(ANIOS, function(a) round(mean(apply(Xs[id_s$anio == a, , drop = FALSE], 2, sd)), 2))), on = "anio"]
print(tabla_f(resumen_anio, "F3_por_anio_resumen"))
cat(sprintf("Lluvia anual media vs ONI (23 años): estaciones r = %.2f | satélite r = %.2f\n",
            cor(resumen_anio$P_estaciones, resumen_anio$ONI_anual), cor(resumen_anio$P_satelite, resumen_anio$ONI_anual)))
# Orden de los grupos fijo (niveles del factor u orden "radix", que no depende del locale):
# así las tablas salen con las mismas filas en cualquier máquina.
niveles_de <- function(g) if (is.factor(g)) levels(droplevels(g)) else sort(unique(g[!is.na(g)]), method = "radix")
grupo_media <- function(Y, g) {
  k <- !is.na(g); g <- factor(as.character(g), levels = niveles_de(g))
  sapply(split(seq_len(nrow(Y))[k], g[k], drop = TRUE), function(i) colMeans(Y[i, , drop = FALSE]))
}
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
#   Las cifras y su interpretación están en EXPLICACION_FUNCIONAL.md, sección 3.6 (corrida del
#   19-nov-2026, con las 98 estaciones de la regla de la OMM); las tablas de esta
#   sección las imprimen las líneas de arriba.

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
  rbindlist(lapply(intersect(niveles_de(pares[[var]]), unique(g[k])), function(v) {
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
  rbindlist(lapply(intersect(niveles_de(pares[[var]]), unique(g[k])), function(v) { i <- k & g == v; bb <- banda_diferencia(D[i, , drop = FALSE], pares$codigo[i])
    data.table(agrupacion = var, grupo = v, dia = T_DIA, d = bb$d, inf = bb$d - 2 * bb$ee, sup = bb$d + 2 * bb$ee) }))
}))
grafico_f(ggplot(curvas_d, aes(dia, d, colour = grupo, fill = grupo)) + geom_hline(yintercept = 0, linetype = 3) +
            geom_ribbon(aes(ymin = inf, ymax = sup), alpha = .15, colour = NA) + geom_line(linewidth = .8) +
            facet_wrap(~agrupacion, ncol = 3, labeller = as_labeller(c(vertiente = "Vertiente", banda_altitud = "Altitud", fase = "Fase ENSO"))) + labs(x = "Día del año", y = "d(t), mm/día", colour = NULL, fill = NULL,
                                                      title = "Diferencia estación - satélite por grupo (±2 EE por estación)"),
          "F4_diferencia_por_grupo", ancho = 14, alto = 4.8)

titulo("4.4 ANOVA funcional del error del satélite (clase, semana 14)")
# Ya no se hace la prueba t (ni la F) funcional PAREADA estación vs satélite. Esa prueba es
# para el mismo individuo medido antes y después de un tratamiento; aquí la pregunta "¿el
# satélite está sesgado?" la responde la regresión, con las bandas de sus efectos (5.4):
# si el satélite midiera lo mismo, alfa(t) = 0 y la integral de beta sobre s = 1
# (regresión de comparación de métodos; Ramsay y Silverman 2005, cap. 13-14).
# Aquí queda la pregunta de grupos INDEPENDIENTES: ¿cambia el error del satélite,
# d(t) = estación - satélite, según la fase ENSO, la vertiente o la altitud?
# Pruebas de la clase: F puntual F_n(t), tipo F global (Zhang 2013, con la aproximación de
# dos momentos de Satterthwaite), GPF = ∫ F_n(t) dt y F_max = sup F_n(t), calibradas por
# permutación de las etiquetas de grupo (Cuevas, Febrero y Fraiman 2004; Górecki y Smaga
# 2015). Unidades independientes: para el ENSO, los AÑOS (curva media de d de cada año; la
# fase es del año, no de la estación); para vertiente y altitud, las ESTACIONES (curva media
# de d de cada estación, porque sus años no son independientes).
B_PERM <- 5000L
w_trap <- c(.5, rep(1, 363), .5)                         # regla del trapecio en la malla diaria
trazas_zhang <- function(G, N, w = w_trap) {             # G: covarianza evaluada en la malla; N: sus grados de libertad
  DG <- G * w; t1 <- sum(diag(DG)); t2 <- sum(DG * t(DG))  # tr(D G) y tr((D G)^2)
  t2_ins <- N^2 / ((N - 1) * (N + 2)) * (t2 - t1^2 / N)
  t1c_ins <- N * (N + 1) / ((N - 1) * (N + 2)) * (t1^2 - 2 * t2 / (N + 1))
  list(tr = t1, beta = t2_ins / t1, d = t1c_ins / t2_ins)
}
D_est <- rowsum(D, pares$codigo) / as.vector(table(pares$codigo)[rownames(rowsum(D, pares$codigo))])
vert_est <- pares[, .(vertiente = vertiente[1], banda_altitud = banda_altitud[1]), by = codigo][match(rownames(D_est), codigo)]
anova_funcional <- function(Xg, g, B_perm = B_PERM) {
  g <- factor(g); k <- nlevels(g); n <- nrow(Xg)
  estad <- function(gg) {
    nk <- tabulate(gg, k); med <- rowsum(Xg, gg, reorder = TRUE) / nk; mg <- colMeans(Xg)
    SSR <- colSums(nk * sweep(med, 2, mg)^2); SSE <- colSums((Xg - med[gg, , drop = FALSE])^2)
    list(Ft = (SSR / (k - 1)) / (SSE / (n - k)), SSRn = sum(w_trap * SSR), SSEn = sum(w_trap * SSE), med = med)
  }
  o_ <- estad(as.integer(g)); E <- Xg - o_$med[as.integer(g), , drop = FALSE]; z <- trazas_zhang(crossprod(E) / (n - k), n - k)
  Fn <- (o_$SSRn / (k - 1)) / (o_$SSEn / (n - k))
  perm <- replicate(B_perm, { e <- estad(sample(as.integer(g))); c(GPF = sum(w_trap * e$Ft), Fmax = max(e$Ft), Fn = (e$SSRn / (k - 1)) / (e$SSEn / (n - k))) })
  Ft_perm <- replicate(min(B_perm, 2000L), estad(sample(as.integer(g)))$Ft)
  list(k = k, n = n, niveles = levels(g), Ft = o_$Ft, med = o_$med, Fn = Fn, d = z$d, p_F = pf(Fn, (k - 1) * z$d, (n - k) * z$d, lower.tail = FALSE),
       GPF = sum(w_trap * o_$Ft), Fmax = max(o_$Ft), p_GPF = (1 + sum(perm["GPF", ] >= sum(w_trap * o_$Ft))) / (B_perm + 1),
       p_Fmax = (1 + sum(perm["Fmax", ] >= max(o_$Ft))) / (B_perm + 1), p_F_perm = (1 + sum(perm["Fn", ] >= Fn)) / (B_perm + 1),
       c_alfa = apply(Ft_perm, 1, quantile, .95), c_max = unname(quantile(perm["Fmax", ], .95)), SSEn = o_$SSEn, tr = z$tr)
}
set.seed(20261118)
D_anio <- rowsum(D, pares$anio) / as.vector(table(pares$anio)[rownames(rowsum(D, pares$anio))])
fase_anio <- oni[match(as.integer(rownames(D_anio)), anio), fase]
cauca_est <- vert_est$vertiente != "Pacífico"
av <- setNames(list(anova_funcional(D_anio, fase_anio), anova_funcional(D_est, vert_est$vertiente),
                    anova_funcional(D_est[cauca_est, ], vert_est$banda_altitud[cauca_est])),
               c(sprintf("Fase ENSO (%d años)", nrow(D_anio)), sprintf("Vertiente (%d estaciones)", nrow(D_est)),
                 sprintf("Altitud en el Cauca (%d estaciones)", sum(cauca_est))))
tabla_anova <- rbindlist(lapply(names(av), function(nm) with(av[[nm]], data.table(
  factor = nm, grupos = paste(niveles, collapse = " / "), k = k, n = n, F_n = round(Fn, 2), d = round(d, 2), p_F = signif(p_F, 2),
  p_F_permutacion = signif(p_F_perm, 2), GPF = round(GPF), p_GPF = signif(p_GPF, 2), F_max = round(Fmax, 2), p_F_max = signif(p_Fmax, 2),
  dias_F_puntual = sum(Ft > c_alfa), dias_criterio_maximo = sum(Ft > c_max)))))
print(tabla_f(tabla_anova, "F4_anova_funcional"))
# Comprobación de la clase (ejercicio 1): SSE_n = (n - k) tr(γ̂); y con k = 2 la F puntual es t² de dos muestras.
chequeo_anova <- c(SSE = max(sapply(av, function(a) abs(a$SSEn - (a$n - a$k) * a$tr) / a$SSEn)),
                   t2 = { a <- av[[2]]; gv <- factor(vert_est$vertiente); x1 <- D_est[gv == levels(gv)[1], 100]; x2 <- D_est[gv == levels(gv)[2], 100]
                          abs(a$Ft[100] - t.test(x1, x2, var.equal = TRUE)$statistic^2) / a$Ft[100] })
cat(sprintf("Comprobaciones: SSE_n = (n-k) tr(γ̂) a %.1e relativo; con k = 2, F(t) = t² a %.1e relativo\n", chequeo_anova[1], chequeo_anova[2]))
grafico_base("F4_anova_funcional", function() {
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
  cols <- c("grey30", "steelblue", "firebrick", "darkgreen")
  for (nm in names(av)) with(av[[nm]], {
    matplot(T_DIA, t(med), type = "l", lty = 1, lwd = 2, col = cols[seq_len(k)], xlab = "Día del año", ylab = "d(t), mm/día", main = nm)
    abline(h = 0, lty = 3); legend("bottomleft", niveles, col = cols[seq_len(k)], lwd = 2, bty = "n", cex = .8)
  })
  for (nm in names(av)) with(av[[nm]], {
    plot(T_DIA, Ft, type = "l", lwd = 2, ylim = c(0, max(Ft, c_max) * 1.05), xlab = "Día del año", ylab = "F(t)",
         main = sprintf("F puntual: GPF p = %.3f, Fmax p = %.3f", p_GPF, p_Fmax))
    lines(T_DIA, c_alfa, lty = 2, col = "steelblue", lwd = 2); abline(h = c_max, lty = 2, col = "firebrick", lwd = 2)
  })
}, ancho = 15, alto = 8)

# QUÉ DICEN LOS DATOS
#   Las cifras y su interpretación están en EXPLICACION_FUNCIONAL.md, sección 4.4 (corrida del
#   19-nov-2026, con las 98 estaciones de la regla de la OMM); las tablas de esta
#   sección las imprimen las líneas de arriba.

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
# El FPCA es una función de los coeficientes y de las filas que lo ajustan: así la prueba
# de no fuga (sección 6) lo puede volver a correr con los datos de fuera alterados.
fpca_ent <- function(coef, filas) {
  mu <- colMeans(coef[filas, ]); Cc <- sweep(coef[filas, ], 2, mu)
  eg <- eigen(W12 %*% (crossprod(Cc) / (sum(filas) - 1)) %*% W12, symmetric = TRUE)
  list(mu = mu, values = eg$values, psi = Wm12 %*% eg$vectors)
}
eg <- fpca_ent(re$coef, ent); mu_c <- eg$mu; psi_c <- eg$psi; psi_todas <- psi_c
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

titulo("5.2b P-19: registros de 1930-1979 y estabilidad del FPCA")
h30 <- rbindlist(lapply(readxl::excel_sheets("datos_curso/P_1930_1979.xlsx"), function(s) as.data.table(readxl::read_excel("datos_curso/P_1930_1979.xlsx", sheet = s))))
h30 <- h30[as.character(STATION_ID) %in% unique(id_est$codigo) & !(MONTH == 2 & DAY == 29) & !is.na(VALUE) & VALUE >= 0 & VALUE <= 1825,
           .(codigo = as.character(STATION_ID), anio = as.integer(YEAR), dia = yday(as.IDate(sprintf("%d-%02d-%02d", YEAR, MONTH, DAY))) - (YEAR %% 4 == 0 & MONTH > 2), P = VALUE)]
h30 <- unique(h30, by = c("codigo", "anio", "dia"))
anios30 <- h30[, .(n_obs = .N, hueco_max = { r <- rle(!(1:365 %in% dia)); if (any(r$values)) max(r$lengths[r$values]) else 0L }), by = .(codigo, anio)]
anios30[, valido := n_obs >= COB_MIN_R5 * 365 & hueco_max <= HUECO_MAX_R5]
h30 <- h30[anios30[valido == TRUE, .(codigo, anio)], on = .(codigo, anio)]
h30[, fila := .GRP, by = .(codigo, anio)]
Y30 <- matrix(NA_real_, max(h30$fila), 365); Y30[cbind(h30$fila, h30$dia)] <- h30$P
C30 <- a_raiz(construir(Y30, T_DIA, K_e)$coef, Phi_e)$coef
eg_30 <- fpca_ent(rbind(re$coef[ent, ], C30), rep(TRUE, sum(ent) + nrow(C30)))
explicada <- function(C, mu, psi, M = 6L) {
  D <- sweep(C, 2, mu); R <- D - (D %*% W_e %*% psi[, 1:M]) %*% t(psi[, 1:M])
  100 * (1 - sum(rowSums((R %*% W_e) * R)) / sum(rowSums((D %*% W_e) * D)))
}
M_ref <- 6L
p19 <- data.table(
  curvas = c("Entrenamiento 1981-1995", "Validación 1996-2000", "Prueba 2001-2003", "Registros 1930-1979 válidos"),
  n = c(sum(ent), sum(val), sum(pru), nrow(C30)),
  pct_explicado_FPCA_entrenamiento = round(c(explicada(re$coef[ent, ], mu_c, psi_c, M_ref), explicada(re$coef[val, ], mu_c, psi_c, M_ref),
                                             explicada(re$coef[pru, ], mu_c, psi_c, M_ref), explicada(C30, mu_c, psi_c, M_ref)), 1))
print(tabla_f(p19, "F5_P19_registros_1930_1979"))
print(tabla_f(anios30[, .(estaciones = uniqueN(codigo), estacion_anios = .N, validos = sum(valido), estaciones_con_validos = uniqueN(codigo[valido]),
                          validos_1970_1979 = sum(valido & anio >= 1970)), by = NULL], "F5_P19_conteo"))
cos_p19 <- svd(t(psi_c[, 1:M_ref]) %*% W_e %*% eg_30$psi[, 1:M_ref])$d
d_mu <- eg_30$mu - mu_c
cat(sprintf("FPCA con 1930-1979 agregados: %.1f %% con M = %d (sin ellos %.1f %%) | cosenos principales entre subespacios: %s | cambio de la media: %.1f %% de su norma\n",
            100 * sum(eg_30$values[1:M_ref]) / sum(eg_30$values), M_ref, 100 * sum(eg$values[1:M_ref]) / sum(eg$values),
            paste(round(cos_p19, 3), collapse = ", "), 100 * sqrt(drop(t(d_mu) %*% W_e %*% d_mu) / drop(t(mu_c) %*% W_e %*% mu_c))))
rm(h30, Y30)

# QUÉ DICEN LOS DATOS (P-19)
#   Las cifras y su interpretación están en EXPLICACION_FUNCIONAL.md, sección 5.2b (corrida del
#   19-nov-2026, con las 98 estaciones de la regla de la OMM); las tablas de esta
#   sección las imprimen las líneas de arriba.

titulo("5.3 Ajuste y elección de M y lambda con validación")
alt_m <- mean(pares$altitud[ent]); alt_s <- sd(pares$altitud[ent])
Zx <- rs$coef %*% W_s
# Los datos entran como argumentos con valor por defecto (los reales) para que la prueba
# de no fuga pueda reajustar con datos alterados sin tocar los objetos del análisis.
diseno <- function(i, con_sat = TRUE, zx = Zx, alt = pares$altitud, oni = pares$ONI_anual, am = alt_m, ds = alt_s)
  cbind(1, if (con_sat) zx[i, , drop = FALSE], (alt[i] - am) / ds, oni[i])
penal <- function(con_sat = TRUE) { P <- matrix(0, 3 + con_sat * K_s, 3 + con_sat * K_s); if (con_sat) P[1 + 1:K_s, 1 + 1:K_s] <- R_s; P }
ajustar_B <- function(lambda, M, con_sat = TRUE, S = S_todos, ...) {
  Z <- diseno(which(ent), con_sat, ...); solve(crossprod(Z) + lambda * penal(con_sat), crossprod(Z, S[ent, 1:M, drop = FALSE]))
}
predecir_raiz <- function(B, i, con_sat = TRUE, psi = psi_c, ...) {
  M <- ncol(B); sweep((diseno(i, con_sat, ...) %*% B) %*% t(Phi_e %*% psi[, 1:M, drop = FALSE]), 2, drop(Phi_e %*% mu_c), "+")
}
a_mm <- function(R) pmax(R, 0)^2
rmse_l2 <- function(P, i, Y = Y_mm) sqrt(rowMeans((P - Y[i, ])^2))

lambdas <- 10^seq(-2, 10, by = .5)
tabla_seleccion <- function(S = S_todos, Y = Y_mm, ...) rbindlist(lapply(seq_len(K_e), function(M) rbindlist(lapply(lambdas, function(l) {
  B <- ajustar_B(l, M, S = S, ...)
  data.table(M = M, lambda = l, rmse_entrenamiento = mean(rmse_l2(a_mm(predecir_raiz(B, which(ent), psi = psi_todas, ...)), which(ent), Y)),
             rmse_validacion = mean(rmse_l2(a_mm(predecir_raiz(B, which(val), psi = psi_todas, ...)), which(val), Y)))
}))))
# Regla: el menor M cuyo mejor RMSE de validación queda a menos de 0,5 % del mínimo.
elegir_M_lambda <- function(sel) {
  mejor <- sel[, .SD[which.min(rmse_validacion)], by = M]
  M_sel <- mejor[rmse_validacion <= min(rmse_validacion) * 1.005, min(M)]
  list(mejor = mejor, M = M_sel, lambda = mejor[M == M_sel, lambda])
}
sel_lambda <- tabla_seleccion()
tabla_f(sel_lambda, "F5_seleccion_M_lambda")
eleccion <- elegir_M_lambda(sel_lambda); mejor_por_M <- eleccion$mejor
print(mejor_por_M)
LAMBDA <- eleccion$lambda; M <- eleccion$M
psi_c <- psi_c[, 1:M, drop = FALSE]; psi_t <- Phi_e %*% psi_c; S <- S_todos[, 1:M, drop = FALSE]
cat(sprintf("Elegidos: M = %d (%.1f %% de la varianza; el menor M a menos de 0,5 %% del mejor RMSE de validación) y lambda = %s\n",
            M, var_fpca$acumulado[M], format(LAMBDA, scientific = TRUE)))
grafico_f(ggplot(sel_lambda, aes(lambda, rmse_validacion, colour = factor(M))) + geom_line() +
            geom_point(data = mejor_por_M, size = 2) + geom_vline(xintercept = LAMBDA, linetype = 2, colour = "red") + scale_x_log10() +
            labs(x = expression(lambda), y = "RMSE L2 medio en validación (mm/día)", colour = "M",
                 title = "Elección de M y lambda con los años de validación (1996-2000)"),
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
# Bandas del 95 % de los efectos funcionales. La clase (semana 15) presenta cada coeficiente
# con su banda y lee dónde excluye el 0 (Ramsay y Silverman 2005, sec. 13.5 y 14.6).
# Aquí la banda sale de un bootstrap de ESTACIONES: se remuestrean estaciones de entrenamiento
# (con todos sus años) y se reajusta B con el mismo FPCA, M y lambda. Así la banda respeta que
# los años de una estación no son independientes (la banda sándwich de clase los trata como tales).
B_BOOT_EF <- 500L
filas_est <- split(which(ent), pares$codigo[ent])
efectos_de <- function(Bm) cbind(int_beta = colSums(Phi_s %*% Bm[1 + 1:K_s, , drop = FALSE] %*% t(psi_t)),
                                 gamma = drop(psi_t %*% Bm[K_s + 2, ]), delta = drop(psi_t %*% Bm[K_s + 3, ]))
set.seed(20261119)
boot_ef <- replicate(B_BOOT_EF, {
  f <- unlist(filas_est[sample(length(filas_est), replace = TRUE)], use.names = FALSE)
  Z <- diseno(f); efectos_de(solve(crossprod(Z) + LAMBDA * penal(), crossprod(Z, S[f, , drop = FALSE])))
})
ef_obs <- efectos_de(B); ef_lo <- apply(boot_ef, c(1, 2), quantile, .025); ef_hi <- apply(boot_ef, c(1, 2), quantile, .975)
lejos_borde <- 16:350
bandas_ef <- data.table(efecto = c("Integral de beta sobre s", "Altitud gamma(t)", "ONI delta(t)"),
                        dias_banda_excluye_0 = colSums(ef_lo > 0 | ef_hi < 0),
                        dias_banda_excluye_1 = c(sum(ef_lo[, 1] > 1 | ef_hi[, 1] < 1), NA, NA),
                        dias_banda_excluye_0_lejos_del_borde = colSums((ef_lo > 0 | ef_hi < 0)[lejos_borde, ]),
                        minimo = round(apply(ef_obs, 2, min), 3), maximo = round(apply(ef_obs, 2, max), 3),
                        ancho_medio_banda = round(colMeans(ef_hi - ef_lo), 3))
print(tabla_f(bandas_ef, "F5_efectos_bandas"))
tabla_f(data.table(dia = T_DIA, ef_obs, inf = ef_lo, sup = ef_hi), "F5_efectos_bandas_por_dia")
grafico_base("F5_efectos", function() {
  par(mfrow = c(1, 3))
  titulos <- list(expression(integral(beta(s, t) * ds)), "Efecto de la altitud gamma(t)", "Efecto del ONI delta(t)")
  ejes <- c("", "raíz(mm/día) por 1 DE", "raíz(mm/día) por unidad de ONI"); cols <- c("black", "brown", 4)
  for (j in 1:3) {
    plot(T_DIA, ef_obs[, j], type = "n", ylim = range(c(ef_lo[, j], ef_hi[, j], 0)), xlab = "Día del año", ylab = ejes[j], main = titulos[[j]])
    polygon(c(T_DIA, rev(T_DIA)), c(ef_lo[, j], rev(ef_hi[, j])), col = adjustcolor(cols[j], .2), border = NA)
    lines(T_DIA, ef_obs[, j], lwd = 2, col = cols[j]); abline(h = c(0, if (j == 1) 1), lty = 3)
  }
  mtext(sprintf("Bandas del 95 %% por bootstrap de estaciones (%d remuestreos)", B_BOOT_EF), side = 1, line = -1.2, outer = TRUE, cex = .75)
}, ancho = 14, alto = 4.8)
cat(sprintf("Integral de beta sobre s: máxima el %s, mínima el %s | efecto ONI más negativo el %s (%.3f) | efecto altitud medio %.3f\n",
            fecha_dia(which.max(colSums(BETA))), fecha_dia(which.min(colSums(BETA))), fecha_dia(which.min(delta_t)), min(delta_t), mean(gamma_t)))

titulo("5.5 Evaluación en prueba (2001-2003) con las métricas de la rúbrica")
# Las cuatro métricas de la rúbrica 4.4, por curva:
#   - error L² relativo  = ||ŷ - y|| / ||y||   (y además el error integrado L² en mm/día)
#   - KGE
#   - sesgo del acumulado anual en %
#   - error en la fecha de los dos picos (días)
# IC del 95 % por bootstrap de ESTACIONES (2 000 remuestreos): las curvas de una
# misma estación no son independientes, así que se remuestrea la estación entera
# (bootstrap por conglomerados; Davison y Hinkley 1997, sec. 3.8). Para comparar dos
# modelos se calcula la diferencia de la métrica en cada estación y su IC con el mismo
# bootstrap: un modelo gana solo si el IC excluye el 0. No se usa Wilcoxon: es una
# prueba clásica de rangos que no aporta sobre el IC y no mira la estructura funcional.
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

# ¿La regresión supera al satélite sin corregir? Diferencia de la métrica en cada estación.
ancho <- dcast(metricas_curva[modelo != "Solo altitud y ONI"], codigo + anio ~ modelo, value.var = c("l2_rel", "rmse", "kge"))
por_estacion <- metricas_curva[, lapply(.SD, mean), by = .(modelo, codigo), .SDcols = c("l2_rel", "rmse", "kge", "sesgo_pct", "err_pico_1", "err_pico_2")]
pe <- dcast(por_estacion, codigo ~ modelo, value.var = c("l2_rel", "rmse", "kge", "sesgo_pct", "err_pico_1", "err_pico_2"))
comparar <- function(v, mejor_si_baja = TRUE, absoluto = FALSE) {
  a <- pe[[paste0(v, "_Regresión funcional")]]; s <- pe[[paste0(v, "_Satélite sin corregir")]]
  if (absoluto) { a <- abs(a); s <- abs(s) }
  d <- a - s; q <- quantile(replicate(N_BOOT, mean(sample(d, replace = TRUE))), c(.025, .975))
  data.table(metrica = v, satelite = mean(s), regresion = mean(a), diferencia = mean(d), ic_inf = q[1], ic_sup = q[2],
             estaciones_que_mejoran = sum(if (mejor_si_baja) d < 0 else d > 0), ic_excluye_0 = q[1] > 0 | q[2] < 0,
             mejora = if (mejor_si_baja) mean(d) < 0 else mean(d) > 0)
}
pareada <- rbind(comparar("l2_rel"), comparar("rmse"), comparar("kge", FALSE), comparar("sesgo_pct", absoluto = TRUE), comparar("err_pico_1"), comparar("err_pico_2"))
pareada[metrica == "sesgo_pct", metrica := "|sesgo_pct|"]
pareada[, c("satelite", "regresion", "diferencia", "ic_inf", "ic_sup") := lapply(.SD, round, 3), .SDcols = c("satelite", "regresion", "diferencia", "ic_inf", "ic_sup")]
cat("Regresión frente al satélite sin corregir, por estación (", nrow(pe), "estaciones; IC 95 % bootstrap de la diferencia):\n")
print(tabla_f(pareada, "F5_regresion_vs_satelite_por_estacion"))
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

titulo("5.5b Contraste con los métodos de regresión vistos en clase")
# La regla del curso: solo métodos vistos en clase. Nuestra regresión (puntajes FPCA de la
# estación sobre ∫φ_j(s)X(s)ds con penalización de rugosidad) es el modelo integral
# funcional-funcional, pero en clase se ajustó de otras dos formas:
#   (a) Taller de repaso, ejercicio 10, paso 3 (Wang et al., 2016): FPCA también de X y
#       regresión de los puntajes de Y sobre los l primeros puntajes de X, B = (ΞᵀΞ)⁻¹ΞᵀH.
#   (b) Taller, ejercicio 9: fda.usc::fregre.basis.fr(), superficie β(s,t) en base tensorial
#       con penalizaciones lambda.s y lambda.t.
# Se ajustan las tres con las mismas curvas √P, el mismo FPCA de la estación (M) y la misma
# regla: parámetros (l, λ_s, λ_t) elegidos con el RMSE de VALIDACIÓN. fregre.basis.fr no
# admite covariables escalares, así que la comparación justa es sin altitud ni ONI; la del
# taller se hace con y sin ellas. Se comparan en prueba con las mismas métricas y bootstrap.
iE_c <- which(ent); iV_c <- which(val)
H_c <- S_todos[, 1:M, drop = FALSE]
curva_de_puntajes <- function(Sc) sweep(Sc %*% t(Phi_e %*% psi_c[, 1:M, drop = FALSE]), 2, drop(Phi_e %*% mu_c), "+")
rmse_val_c <- function(R) mean(rmse_l2(a_mm(R[iV_c, , drop = FALSE]), iV_c))
extras_c <- cbind((pares$altitud - alt_m) / alt_s, pares$ONI_anual)

# (a) FPCA del satélite (solo entrenamiento) y regresión puntajes sobre puntajes
ev_s <- eigen(W_s, symmetric = TRUE)
W12_s <- ev_s$vectors %*% diag(sqrt(ev_s$values)) %*% t(ev_s$vectors); Wm12_s <- ev_s$vectors %*% diag(1 / sqrt(ev_s$values)) %*% t(ev_s$vectors)
mu_x <- colMeans(rs$coef[ent, ])
eg_x <- eigen(W12_s %*% (crossprod(sweep(rs$coef[ent, ], 2, mu_x)) / (sum(ent) - 1)) %*% W12_s, symmetric = TRUE)
psi_x <- Wm12_s %*% eg_x$vectors
Xi <- sweep(rs$coef, 2, mu_x) %*% W_s %*% psi_x                      # puntajes de X en todas las curvas
ajuste_taller <- function(l, con_extras) {
  Z <- cbind(1, Xi[, seq_len(l), drop = FALSE], if (con_extras) extras_c)
  B <- solve(crossprod(Z[iE_c, ]), crossprod(Z[iE_c, ], H_c[iE_c, ]))
  list(B = B, Z = Z, R = curva_de_puntajes(Z %*% B))
}
sel_taller <- rbindlist(lapply(c(TRUE, FALSE), function(ce) rbindlist(lapply(seq_len(K_s), function(l)
  data.table(metodo = "taller PC", con_altitud_ONI = ce, l = l, rmse_validacion = rmse_val_c(ajuste_taller(l, ce)$R))))))
l_con <- sel_taller[con_altitud_ONI == TRUE][which.min(rmse_validacion), l]
l_sin <- sel_taller[con_altitud_ONI == FALSE][which.min(rmse_validacion), l]
# El taller dice "lm() solo para verificar": la B de mínimos cuadrados debe ser la de lm()
aj_t <- ajuste_taller(l_con, TRUE)
dif_lm_taller <- max(abs(aj_t$B - coef(lm(H_c[iE_c, ] ~ aj_t$Z[iE_c, -1]))))

# Nuestra regresión sin altitud ni ONI (λ elegido igual, en validación)
ajuste_nuestra_sin <- function(l) {
  Z <- cbind(1, Zx); P <- matrix(0, ncol(Z), ncol(Z)); P[-1, -1] <- R_s
  curva_de_puntajes(Z %*% solve(crossprod(Z[iE_c, ]) + l * P, crossprod(Z[iE_c, ], H_c[iE_c, ])))
}
sel_nuestra <- data.table(metodo = "nuestra", con_altitud_ONI = FALSE, lambda = lambdas, rmse_validacion = sapply(lambdas, function(l) rmse_val_c(ajuste_nuestra_sin(l))))
LAMBDA_SIN <- sel_nuestra[which.min(rmse_validacion), lambda]

# (b) fda.usc::fregre.basis.fr con las mismas bases B-spline (nodos equiespaciados, como las nuestras)
X_fd <- fdata(rs$coef %*% t(Phi_s), argvals = T_DIA); Y_fd <- fdata(re$coef %*% t(Phi_e), argvals = T_DIA)
base_s_fda <- create.bspline.basis(RANGO, K_s); base_t_fda <- create.bspline.basis(RANGO, K_e)
ajuste_fr <- function(ls, lt) fregre.basis.fr(X_fd[iE_c], Y_fd[iE_c], basis.s = base_s_fda, basis.t = base_t_fda, lambda.s = ls, lambda.t = lt)
rej_fr <- CJ(lambda_s = 10^seq(-2, 14, 2), lambda_t = 10^c(2, 6, 10, 14))
sel_fr <- rej_fr[, .(metodo = "fregre.basis.fr", con_altitud_ONI = FALSE, rmse_validacion = rmse_val_c(predict(ajuste_fr(lambda_s, lambda_t), X_fd)$data)), by = .(lambda_s, lambda_t)]
mejor_fr <- sel_fr[which.min(rmse_validacion)]
fit_fr <- ajuste_fr(mejor_fr$lambda_s, mejor_fr$lambda_t)
stopifnot(mejor_fr$lambda_s < max(rej_fr$lambda_s))   # el óptimo no queda en el borde de la rejilla
tabla_f(rbindlist(list(sel_taller, sel_nuestra, sel_fr), fill = TRUE), "F5_contraste_seleccion_validacion")

modelos_clase <- list(
  `Satélite sin corregir` = X_mm[i_pru, ],
  `Nuestra regresión (con altitud y ONI)` = a_mm(pred_puntajes),
  `Taller PC (con altitud y ONI)` = a_mm(aj_t$R[i_pru, ]),
  `Nuestra regresión (sin altitud ni ONI)` = a_mm(ajuste_nuestra_sin(LAMBDA_SIN)[i_pru, ]),
  `Taller PC (sin altitud ni ONI)` = a_mm(ajuste_taller(l_sin, FALSE)$R[i_pru, ]),
  `fda.usc fregre.basis.fr (sin altitud ni ONI)` = a_mm(predict(fit_fr, X_fd[i_pru])$data))
mc_clase <- rbindlist(lapply(names(modelos_clase), function(nm) metricas_de(modelos_clase[[nm]], i_pru, nm)))
contraste <- mc_clase[, resumir(.SD), by = modelo]
contraste[, parametros := c("", sprintf("M = %d, λ = %.1e", M, LAMBDA), sprintf("M = %d, l = %d", M, l_con),
                            sprintf("M = %d, λ = %.1e", M, LAMBDA_SIN), sprintf("M = %d, l = %d", M, l_sin),
                            sprintf("λs = %.0e, λt = %.0e", mejor_fr$lambda_s, mejor_fr$lambda_t))]
print(tabla_f(contraste, "F5_contraste_metodos_clase"))

# Diferencia de la métrica en cada estación con IC 95 % bootstrap, mismas reglas que 5.5
pe_c <- mc_clase[, .(l2_rel = mean(l2_rel), kge = mean(kge)), by = .(modelo, codigo)]
dif_par <- function(a, b, v) {
  x <- pe_c[modelo == a][order(codigo)][[v]]; y <- pe_c[modelo == b][order(codigo)][[v]]; d <- x - y
  q <- quantile(replicate(N_BOOT, mean(sample(d, replace = TRUE))), c(.025, .975))
  data.table(comparacion = paste(a, "−", b), metrica = v, diferencia = round(mean(d), 4), ic_inf = round(q[1], 4), ic_sup = round(q[2], 4),
             ic_excluye_0 = q[1] > 0 | q[2] < 0)
}
pares_clase <- list(c("Nuestra regresión (con altitud y ONI)", "Taller PC (con altitud y ONI)"),
                    c("Nuestra regresión (sin altitud ni ONI)", "Taller PC (sin altitud ni ONI)"),
                    c("Nuestra regresión (sin altitud ni ONI)", "fda.usc fregre.basis.fr (sin altitud ni ONI)"),
                    c("Taller PC (sin altitud ni ONI)", "fda.usc fregre.basis.fr (sin altitud ni ONI)"),
                    c("Nuestra regresión (con altitud y ONI)", "Nuestra regresión (sin altitud ni ONI)"))
contraste_par <- rbindlist(lapply(pares_clase, function(p) rbind(dif_par(p[1], p[2], "l2_rel"), dif_par(p[1], p[2], "kge"))))
print(tabla_f(contraste_par, "F5_contraste_por_estacion"))

# QUÉ DICEN LOS DATOS (5.5b)
#   Las cifras y su interpretación están en EXPLICACION_FUNCIONAL.md, sección 5.5b (corrida del
#   19-nov-2026, con las 98 estaciones de la regla de la OMM); las tablas de esta
#   sección las imprimen las líneas de arriba.

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
q_sens <- quantile(replicate(N_BOOT, mean(sample(ps$a - ps$s, replace = TRUE))), c(.025, .975))
cat(sprintf("Con todos los años el error L² relativo en prueba es menor en %d de %d estaciones (diferencia media %.4f, IC 95 %% [%.4f, %.4f])\n",
            sum(ps$a < ps$s), nrow(ps), mean(ps$a - ps$s), q_sens[1], q_sens[2]))

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

titulo("5.6b Inferencia y diagnóstico de residuos")
# (a) Prueba F por permutación del aporte del satélite, a la manera de Fperm.fd (clase, semana 15):
#     F(t) = Σ_i (ŷ_i(t) - ȳ̂(t))² / Σ_i (y_i(t) - ŷ_i(t))², en la escala raíz y en entrenamiento.
#     H0: la curva del satélite no aporta nada más allá de la altitud y el ONI. Se permutan las
#     curvas del satélite ENTRE estaciones del mismo año (se conservan el año y su fase ENSO) y
#     se reajusta B con el mismo FPCA, M y lambda. Se reporta el criterio del máximo y el integrado.
Y_raiz <- re$coef %*% t(Phi_e); fil_ent <- which(ent)
F_t_de <- function(P, filas) colSums(sweep(P, 2, colMeans(P))^2) / colSums((Y_raiz[filas, ] - P)^2)
F_obs <- F_t_de(predecir_raiz(B, fil_ent), fil_ent)
anio_ent <- pares$anio[fil_ent]
perm_dentro_anio <- function() { nuevo <- fil_ent; for (a in unique(anio_ent)) { k <- which(anio_ent == a); nuevo[k] <- fil_ent[k][sample.int(length(k))] }; nuevo }
B_FPERM <- 300L
set.seed(20261120)
F_perm <- replicate(B_FPERM, { zxp <- Zx; zxp[fil_ent, ] <- Zx[perm_dentro_anio(), ]
  F_t_de(predecir_raiz(ajustar_B(LAMBDA, M, zx = zxp), fil_ent, zx = zxp), fil_ent) })
c_alfa_F <- apply(F_perm, 1, quantile, .95); c_max_F <- quantile(apply(F_perm, 2, max), .95)
fperm <- data.table(estadistico = c("F máximo", "F integrado (media sobre el año)"), observado = round(c(max(F_obs), mean(F_obs)), 3),
                    critico_95 = round(c(c_max_F, quantile(colMeans(F_perm), .95)), 3),
                    p_permutacion = signif(c((1 + sum(apply(F_perm, 2, max) >= max(F_obs))) / (B_FPERM + 1), (1 + sum(colMeans(F_perm) >= mean(F_obs))) / (B_FPERM + 1)), 2),
                    dias_sobre_el_critico_puntual = c(sum(F_obs > c_alfa_F), NA), permutaciones = B_FPERM)
print(tabla_f(fperm, "F5_permutacion_satelite"))

# (b) Diagnóstico de residuos (clase, correlación de los residuos; Ramsay y Silverman 2005, cap. 16):
#     residuos estudentizados por día con la DE de entrenamiento, residuos frente a predichos y
#     correlación de los residuos entre días.
res_ent <- Y_raiz[fil_ent, ] - predecir_raiz(B, fil_ent); sig_t <- apply(res_ent, 2, sd)
res_pru <- Y_raiz[i_pru, ] - pred_puntajes; stud_pru <- sweep(res_pru, 2, sig_t, "/")
cor_res <- cor(res_ent)
diag_res <- data.table(medida = c("Residuos estudentizados de prueba fuera de ±2 (%)", "Residuos estudentizados de prueba fuera de ±3 (%)",
                                  "Residuo medio de prueba (raíz de mm/día)", "Correlación de los residuos a 1 día", "Correlación de los residuos a 30 días",
                                  "Correlación de los residuos a 90 días", "Correlación residuo-predicho en prueba"),
                       valor = round(c(100 * mean(abs(stud_pru) > 2), 100 * mean(abs(stud_pru) > 3), mean(res_pru),
                                       mean(cor_res[cbind(1:364, 2:365)]), mean(cor_res[cbind(1:335, 31:365)]), mean(cor_res[cbind(1:275, 91:365)]),
                                       cor(as.vector(res_pru), as.vector(pred_puntajes))), 3))
print(tabla_f(diag_res, "F5_diagnostico_residuos"))

# (c) Traza-variograma de los residuos (Giraldo, Delicado y Mateu 2011): γ(h) = ½ E||e_i - e_j||²
#     para pares de estaciones del mismo año a distancia h. Si crece con h, los residuos de
#     estaciones cercanas se parecen: hay correlación espacial que el modelo no recoge. Se ajusta
#     el modelo exponencial, γ(h) = c0 + σ²(1 - exp(-h/φ)), por mínimos cuadrados ponderados (Cressie 1985).
coord_par <- est[match(pares$codigo[fil_ent], codigo), .(lat, lon)]
tv <- rbindlist(lapply(sort(unique(anio_ent)), function(a) {
  k <- which(anio_ent == a); E <- res_ent[k, , drop = FALSE]
  G <- E %*% (t(E) * w_trap); nn <- diag(G); D2 <- outer(nn, nn, "+") - 2 * G
  H <- outer(seq_along(k), seq_along(k), function(i, j) haversine(coord_par$lat[k[i]], coord_par$lon[k[i]], coord_par$lat[k[j]], coord_par$lon[k[j]]))
  u <- which(upper.tri(D2), arr.ind = TRUE); data.table(anio = a, h = H[u], g = D2[u] / 2)
}))
tv[, bin := cut(h, c(seq(0, 100, by = 10), 150, 250), include.lowest = TRUE)]
tv_bin <- tv[, .(h = mean(h), gamma = mean(g), pares = .N), by = bin][order(h)]
obj_tv <- function(p) with(tv_bin, sum(pares * (gamma - (p[1] + p[2] * (1 - exp(-h / p[3]))))^2))
fit_tv <- optim(c(min(tv_bin$gamma), diff(range(tv_bin$gamma)), 30), obj_tv, method = "L-BFGS-B", lower = c(0, 0, 1), upper = c(Inf, Inf, 500))$par
tabla_tv <- tv_bin[, .(bin = as.character(bin), h_km = round(h, 1), traza_variograma = round(gamma, 2), pares)]
print(tabla_f(tabla_tv, "F5_traza_variograma_residuos"))
cat(sprintf("Traza-variograma exponencial de los residuos: pepita c0 = %.1f, meseta parcial σ² = %.1f, φ = %.1f km (rango práctico %.0f km); la pepita es %.0f %% de la meseta\n",
            fit_tv[1], fit_tv[2], fit_tv[3], 3 * fit_tv[3], 100 * fit_tv[1] / (fit_tv[1] + fit_tv[2])))
tabla_f(data.table(parametro = c("pepita_c0", "meseta_parcial", "phi_km", "rango_practico_km", "pct_pepita"),
                   valor = round(c(fit_tv, 3 * fit_tv[3], 100 * fit_tv[1] / sum(fit_tv[1:2])), 2)), "F5_traza_variograma_ajuste")

grafico_base("F5_diagnostico", function() {
  par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))
  smoothScatter(as.vector(pred_puntajes), as.vector(res_pru), nrpoints = 0, xlab = "Predicho (raíz de mm/día)", ylab = "Residuo",
                main = "Residuos frente a predichos (prueba)"); abline(h = 0, lty = 2, col = "red")
  q <- apply(stud_pru, 2, quantile, c(.05, .25, .5, .75, .95))
  matplot(T_DIA, t(q), type = "l", lty = c(3, 2, 1, 2, 3), lwd = c(1, 1, 2, 1, 1), col = 1, ylim = c(-4, 4), xlab = "Día del año",
          ylab = "Residuo estudentizado", main = "Residuos estudentizados por día (prueba): cuantiles 5-95 %"); abline(h = c(-2, 2), lty = 3, col = "red")
  contour(T_DIA, T_DIA, cor_res, nlevels = 10, xlab = "Día s", ylab = "Día t", main = "Correlación de los residuos (entrenamiento)")
  with(tv_bin, plot(h, gamma, pch = 19, cex = sqrt(pares) / max(sqrt(pares)) * 2 + .5, ylim = c(0, max(gamma) * 1.1),
                    xlab = "Distancia entre estaciones (km)", ylab = "Traza-variograma", main = "Traza-variograma de los residuos"))
  hh <- seq(0, 250, length.out = 200); lines(hh, fit_tv[1] + fit_tv[2] * (1 - exp(-hh / fit_tv[3])), col = "red", lwd = 2)
}, ancho = 12, alto = 9)

grafico_base("F5_permutacion_satelite", function() {
  matplot(T_DIA, F_perm, type = "l", lty = 1, col = "grey85", ylim = c(0, max(F_obs, F_perm) * 1.05), xlab = "Día del año", ylab = "F(t)",
          main = sprintf("Aporte del satélite: F(t) observado frente a %d permutaciones dentro de cada año", B_FPERM))
  lines(T_DIA, F_obs, lwd = 2); lines(T_DIA, c_alfa_F, lty = 2, col = "steelblue", lwd = 2); abline(h = c_max_F, lty = 2, col = "firebrick", lwd = 2)
  legend("topright", c("Observado", "Permutaciones", "Crítico puntual 95 %", "Crítico del máximo 95 %"), col = c(1, "grey70", "steelblue", "firebrick"),
         lty = c(1, 1, 2, 2), lwd = 2, bty = "n", cex = .85)
}, ancho = 10, alto = 5)

# (d) Curvas observadas y ajustadas de 20 estaciones de un año de prueba.
#     Se toma 2002 (año Niño) y 20 estaciones repartidas por cuantiles del error.
k02 <- i_pru[pares$anio[i_pru] == 2002]
err02 <- metricas_curva[modelo == "Regresión funcional" & anio == 2002, l2_rel]
sel20 <- k02[order(err02)][round(seq(1, length(k02), length.out = 20))]
grafico_base("F5_curvas_20_estaciones", function() {
  par(mfrow = c(4, 5), mar = c(2.5, 2.5, 2.5, .5), mgp = c(1.5, .5, 0))
  for (k in sel20) {
    j <- match(k, i_pru)
    matplot(T_DIA, cbind(Y_mm[k, ], modelos[[1]][j, ], modelos[[2]][j, ]), type = "l", lty = c(1, 2, 1), lwd = c(2, 1, 2),
            col = c(1, "darkorange", 4), xlab = "", ylab = "mm/día", main = sprintf("%s (L² %.2f)", substr(pares$nombre[k], 1, 18), err02[match(k, k02)]), cex.main = .85)
  }
}, ancho = 15, alto = 11)

# (e) Mapa del error por estación en prueba.
mapa_err <- merge(metricas_curva[modelo != "Solo altitud y ONI", .(l2 = mean(l2_rel)), by = .(codigo, modelo)][, dcast(.SD, codigo ~ modelo, value.var = "l2")],
                  est[, .(codigo, nombre, lat, lon, vertiente)], by = "codigo")
setnames(mapa_err, c("Regresión funcional", "Satélite sin corregir"), c("l2_reg", "l2_sat"))
mapa_err[, mejora := l2_reg - l2_sat]
tabla_f(mapa_err[order(mejora)], "F5_error_por_estacion_prueba")
cat(sprintf("Prueba, por estación: la regresión baja el error L² relativo en %d de %d estaciones; error medio %.2f en el Cauca y %.2f en el Pacífico\n",
            sum(mapa_err$mejora < 0), nrow(mapa_err), mapa_err[vertiente != "Pacífico", mean(l2_reg)], mapa_err[vertiente == "Pacífico", mean(l2_reg)]))
valle <- sf::st_read("datos_web/limites/deptos_valle_y_vecinos_DANE_2020.gpkg", quiet = TRUE)
grafico_f(ggplot() + geom_sf(data = valle, fill = "grey97", colour = "grey60") +
            geom_point(data = mapa_err, aes(lon, lat, colour = mejora, size = l2_reg)) +
            scale_colour_gradient2(low = "steelblue", mid = "grey80", high = "firebrick", midpoint = 0) + scale_size(range = c(1, 4.5)) +
            coord_sf(xlim = range(mapa_err$lon) + c(-.15, .15), ylim = range(mapa_err$lat) + c(-.15, .15)) +
            labs(x = NULL, y = NULL, colour = "Regresión - satélite\n(error L², < 0 = mejora)", size = "Error L² relativo\nde la regresión",
                 title = "Error por estación en prueba (2001-2003)"),
          "F5_mapa_error", ancho = 8, alto = 7.5)

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
#   Las cifras y su interpretación están en EXPLICACION_FUNCIONAL.md, sección 5.7 (corrida del
#   19-nov-2026, con las 98 estaciones de la regla de la OMM); las tablas de esta
#   sección las imprimen las líneas de arriba.

regresion <- list(M = M, lambda = LAMBDA, B = B, mu_c = mu_c, psi_c = psi_c, var_fpca = var_fpca, BETA = BETA, metricas = metricas,
             metricas_curva = metricas_curva, dif_equivalencia = dif_equiv, alt_m = alt_m, alt_s = alt_s)

# =============================================================================
# 6. PRUEBAS Y CIERRE
# =============================================================================
titulo("6. Pruebas de la parte funcional")
# Prueba de no fuga. (a) Se reemplaza por ruido TODO lo que no es entrenamiento (curvas
# de estación, píxel, altitud y ONI de validación y prueba) y se reajusta con las mismas
# funciones: media, autovalores, escalas de altitud y B tienen que salir idénticos.
# (b) Se altera solo la prueba y se vuelve a elegir M y lambda: tienen que salir iguales.
# (c) K se eligió con ent_est / ent_sat: se revisa que esas filas sean solo de 1981-1995.
set.seed(7)
ruido_filas <- function(A, filas) { A[filas, ] <- matrix(rnorm(sum(filas) * ncol(A), sd = 3 * sd(A)), sum(filas)); A }
alterar <- function(filas) {
  coef <- ruido_filas(re$coef, filas)
  list(coef = coef, Y = abs(ruido_filas(Y_mm, filas)), zx = ruido_filas(Zx, filas),
       alt = replace(pares$altitud, filas, runif(sum(filas), 0, 4500)), oni = replace(pares$ONI_anual, filas, rnorm(sum(filas), sd = 2)))
}
d <- alterar(!ent)
fp_r <- fpca_ent(d$coef, ent); am_r <- mean(d$alt[ent]); ds_r <- sd(d$alt[ent])
S_r <- sweep(d$coef, 2, fp_r$mu) %*% W_e %*% fp_r$psi
B_r <- ajustar_B(LAMBDA, M, S = S_r, zx = d$zx, alt = d$alt, oni = d$oni, am = am_r, ds = ds_r)
no_fuga_ent <- c(media = max(abs(fp_r$mu - mu_c)), autovalores = max(abs(fp_r$values - var_fpca$autovalor)),
                 escalas = max(abs(c(am_r, ds_r) - c(alt_m, alt_s))), B = max(abs(B_r - B)))
d <- alterar(pru)
S_p <- sweep(d$coef, 2, mu_c) %*% W_e %*% psi_todas
el_p <- elegir_M_lambda(tabla_seleccion(S = S_p, Y = d$Y, zx = d$zx, alt = d$alt, oni = d$oni))
no_fuga_sel <- el_p$M == M && el_p$lambda == LAMBDA
no_fuga_K <- all(id_est$anio[ent_est] <= 1995) && all(id_sat$anio[ent_sat] <= 1995)
cat(sprintf("No fuga: con validación y prueba reemplazadas por ruido, diferencia máxima en media %.1e, autovalores %.1e, escalas %.1e y B %.1e
",
            no_fuga_ent[1], no_fuga_ent[2], no_fuga_ent[3], no_fuga_ent[4]))
cat(sprintf("No fuga: con la prueba alterada se vuelve a elegir M = %d y lambda = %s (original M = %d, lambda = %s)
",
            el_p$M, format(el_p$lambda, scientific = TRUE), M, format(LAMBDA, scientific = TRUE)))
rm(d)
pruebas <- list(
  "El emparejamiento es el mismo del archivo de datos" = identical(pares[order(fila_est), fila_sat], dp$pares[order(fila_est), fila_sat]) && identical(pares[order(fila_est), codigo], dp$pares[order(fila_est), codigo]),
  "Ningún año queda fuera del modelo" = all(pares$usar) && setequal(unique(pares$anio[ent]), 1981:1995) && setequal(unique(pares$anio[val]), 1996:2000) && setequal(unique(pares$anio[pru]), 2001:2003),
  "La base propia coincide con splines::splineDesign" = max(abs(base_bspline(T_DIA, 11) - ref_splines)) < 1e-12,
  "S simétrica, idempotente y tr(S) = K en las dos fuentes" = all(controles[, .(simetria, idempotencia, traza_menos_K)] < 1e-8),
  "W integra el dominio: 1'W1 = 364" = abs(sum(W_e) - 364) < 1e-8 && abs(sum(W_s) - 364) < 1e-8,
  "La media paso a paso coincide con fda::mean.fd" = all(verif_media$vs_fda_mean.fd < 1e-10),
  "Las atípicas coinciden con fda::fbplot" = all(verif_fb),
  "Funciones propias ortonormales con W" = max(abs(t(psi_c) %*% W_e %*% psi_c - diag(M))) < 1e-10,
  "Sin fuga: K, FPCA, escalas y B no cambian si se alteran validación y prueba; M y lambda no cambian si se altera la prueba" =
    no_fuga_K && all(no_fuga_ent < 1e-10) && no_fuga_sel,
  "Equivalencia beta(s,t) vs puntajes FPCA (<= 1e-4)" = dif_equiv <= 1e-4,
  "La reconstrucción en mm no es negativa" = min(modelos[["Regresión funcional"]]) >= 0,
  "Bootstrap con al menos 1000 remuestreos" = N_BOOT >= 1000,
  "ANOVA funcional: SSE_n = (n-k) tr(γ̂) y, con k = 2, F(t) = t² (clase, ejercicio 1)" = all(chequeo_anova < 1e-8),
  "Bandas por bootstrap de estaciones y permutación del satélite con réplicas suficientes" = B_BOOT_EF >= 500 && B_FPERM >= 200 && B_PERM >= 1000,
  "Contraste con los métodos de clase: la B del taller es la de lm() y λ de fregre.basis.fr no queda en el borde" =
    dif_lm_taller < 1e-8 && mejor_fr$lambda_s < max(rej_fr$lambda_s)
)
res_pruebas <- data.table(prueba = names(pruebas), pasa = unlist(pruebas))
print(tabla_f(res_pruebas, "F6_pruebas"))
cat(sprintf("\nPruebas: %d de %d pasan\n", sum(res_pruebas$pasa), nrow(res_pruebas)))
stopifnot(res_pruebas$pasa)
saveRDS(list(construccion = o, eda = eda, comparacion = comparacion, regresion = regresion, pruebas = res_pruebas),
        file.path(SALIDA_F, "resultados_funcional.rds"), compress = "xz")
writeLines(capture.output(sessionInfo()), file.path(SALIDA_F, "sessionInfo.txt"))

titulo("7. Evaluación de la red funcional (predicciones de lluvia_valle_red.py)")
# Python entrena la red (PyTorch, lo único permitido en ese apartado) y exporta los coeficientes
# √P predichos para todas las curvas por cada cabeza y semilla. Aquí se evalúan con las MISMAS
# funciones de la sección 5.5 (métricas y bootstrap por estación), así la regresión y
# las tres cabezas se miden exactamente igual y todo el análisis usa herramientas de clase.
# Orden de ejecución: este archivo → python lluvia_valle_red.py → este archivo otra vez.
ARCH_PRED <- file.path(SALIDA_PY, "predicciones_coef_raiz.csv")
if (!file.exists(ARCH_PRED)) {
  cat("No está", ARCH_PRED, "→ corra python lluvia_valle_red.py y vuelva a correr este archivo para evaluar la red.\n")
} else {
  set.seed(20261114)
  pr <- fread(ARCH_PRED, colClasses = list(character = "codigo"))
  cols_c <- paste0("c", seq_len(K_e))
  # Control 1: las filas son los mismos pares, en el mismo orden
  ok_filas <- pr[, all(codigo == pares$codigo[fila] & anio == pares$anio[fila] & particion == pares$particion[fila])]
  # Control 2: la regresión de R reconstruida en Python da las mismas curvas que aquí
  rp <- pr[modelo == "regresión R (en Python)"][order(fila)]
  dif_py_R <- max(abs(as.matrix(rp[, ..cols_c]) %*% t(Phi_e) - predecir_raiz(B, rp$fila)))
  cat(sprintf("Regresión de R hecha en Python vs en R: diferencia máxima %.1e (√(mm/día))\n", dif_py_R))

  cabezas <- pr[modelo %in% c("lineal", "kan", "mlp")]
  i_vp <- which(val | pru)
  mc_red <- rbindlist(c(
    list(metricas_de(X_mm[i_vp, ], i_vp, "satélite sin corregir")[, semilla := 0L],
         metricas_de(a_mm(predecir_raiz(B, i_vp)), i_vp, "regresión (R)")[, semilla := 0L]),
    lapply(split(cabezas[fila %in% i_vp], by = c("modelo", "semilla")), function(d) {
      d <- d[order(fila)]
      metricas_de(a_mm(as.matrix(d[, ..cols_c]) %*% t(Phi_e)), d$fila, d$modelo[1])[, semilla := d$semilla[1]]
    })))
  tabla_f(mc_red, "F7_metricas_red_por_curva")
  vars_m <- c("l2_rel", "rmse", "kge", "sesgo_pct", "err_pico_1", "err_pico_2")

  # (a) Media ± DE entre semillas: en cada semilla, media de las medias por estación (como en 5.5)
  por_semilla <- mc_red[, lapply(.SD, mean), by = .(particion, modelo, semilla, codigo), .SDcols = vars_m][
    , lapply(.SD, mean), by = .(particion, modelo, semilla), .SDcols = vars_m]
  ablacion <- por_semilla[, c(.(semillas = .N), unlist(lapply(.SD, function(v) list(media = round(mean(v), 3), DE = round(sd(v), 3))), recursive = FALSE)),
                          by = .(particion, modelo), .SDcols = vars_m]
  print(tabla_f(ablacion, "F7_ablacion_media_DE_semillas"))

  # (b) IC 95 % bootstrap por estación, con las métricas de cada curva promediadas en las 3 semillas
  mc_prom <- mc_red[, lapply(.SD, mean), by = .(particion, modelo, codigo, anio, fase), .SDcols = vars_m]
  ic_red <- mc_prom[, resumir(.SD), by = .(particion, modelo)]
  print(tabla_f(ic_red, "F7_ablacion_IC_bootstrap"))

  # (c) Comparaciones por estación en PRUEBA: una cabeza supera a otra solo si el IC 95 %
  #     bootstrap de la diferencia por estación excluye el 0 (rúbrica 4.3).
  pe_red <- mc_prom[particion == "prueba", lapply(.SD, mean), by = .(modelo, codigo), .SDcols = vars_m]
  comparar_red <- function(a, b, v) {
    x <- pe_red[modelo == a][order(codigo)][[v]]; y <- pe_red[modelo == b][order(codigo)][[v]]
    if (v == "sesgo_pct") { x <- abs(x); y <- abs(y) }
    d <- x - y; q <- quantile(replicate(N_BOOT, mean(sample(d, replace = TRUE))), c(.025, .975))
    excl <- q[1] > 0 | q[2] < 0; mejor_a <- if (v == "kge") mean(d) > 0 else mean(d) < 0
    data.table(comparacion = paste(a, "−", b), metrica = if (v == "sesgo_pct") "|sesgo_pct|" else v,
               media_a = round(mean(x), 3), media_b = round(mean(y), 3), diferencia = round(mean(d), 4), ic_inf = round(q[1], 4), ic_sup = round(q[2], 4),
               estaciones_a_mejor = sum(if (v == "kge") d > 0 else d < 0), ic_excluye_0 = excl,
               gana = if (!excl) "ninguna (IC con 0)" else if (mejor_a) a else b)
  }
  pares_red <- list(c("kan", "lineal"), c("mlp", "lineal"), c("kan", "mlp"), c("lineal", "regresión (R)"),
                    c("lineal", "satélite sin corregir"), c("kan", "satélite sin corregir"), c("mlp", "satélite sin corregir"))
  comp_red <- rbindlist(lapply(pares_red, function(p) rbindlist(lapply(vars_m, function(v) comparar_red(p[1], p[2], v)))))
  print(tabla_f(comp_red, "F7_comparaciones_por_estacion_prueba"))

  # (d) Por fase ENSO (rúbrica 4.4). La prueba no tiene años Niña: la validación se agrega como
  #     ORIENTATIVA (con ella se eligieron λ_w y la parada temprana, no los pesos).
  fase_red <- mc_prom[, resumir(.SD), by = .(particion, fase, modelo)]
  fase_red[, anios := sapply(seq_len(.N), function(k) paste(sort(unique(pares[particion == fase_red$particion[k] & fase == fase_red$fase[k], anio])), collapse = " "))]
  setcolorder(fase_red, c("particion", "fase", "anios")); setorder(fase_red, -particion, fase, modelo)
  print(tabla_f(fase_red, "F7_metricas_por_fase_ENSO"))

  pruebas_red <- list(
    "Las predicciones de Python corresponden a los mismos pares y en el mismo orden" = ok_filas,
    "La regresión de R hecha en Python coincide con la de R (< 1e-8)" = dif_py_R < 1e-8,
    "Tres semillas por cabeza y todas las curvas en cada una" = cabezas[, .N, by = .(modelo, semilla)][, uniqueN(semilla) == 3 && all(N == nrow(pares)), by = modelo][, all(V1)],
    "Métricas finitas en todas las curvas de validación y prueba" = mc_red[, all(is.finite(l2_rel) & is.finite(kge))],
    "Bootstrap con al menos 1000 remuestreos" = N_BOOT >= 1000)
  res_red <- data.table(prueba = names(pruebas_red), pasa = unlist(pruebas_red))
  print(tabla_f(res_red, "F7_pruebas_red"))
  cat(sprintf("\nPruebas de la evaluación de la red: %d de %d pasan\n", sum(res_red$pasa), nrow(res_red)))
  stopifnot(res_red$pasa)
}

# QUÉ DICEN LOS DATOS (sección 7)
#   Las cifras y su interpretación están en EXPLICACION_FUNCIONAL.md, sección 7 (corrida del
#   19-nov-2026, con las 98 estaciones de la regla de la OMM); las tablas de esta
#   sección las imprimen las líneas de arriba.

titulo("8. Interpretación de la red (rúbrica 4.5)")
# Se grafica y se lee lo que aprendió la red. Python solo exporta números (pesos del codificador,
# de la cabeza lineal y activaciones de la primera capa de la KAN en el rango de los datos);
# todo el análisis se hace aquí. Se usa la semilla 2026 para las figuras y las tres semillas
# para ver qué tan estable es cada lectura.
ARCH_INT <- file.path(SALIDA_PY, c("pesos_codificador.csv", "cabeza_lineal_pesos.csv", "kan_activaciones_capa1.csv", "kan_importancia_aristas.csv"))
if (!all(file.exists(ARCH_INT)) || !exists("mc_prom")) {
  cat("Faltan las salidas de lluvia_valle_red.py → la sección 8 se omite.\n")
} else {
  SEM <- 2026L
  meses <- data.table(dia = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335), mes = c("E", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"))
  cod <- fread(ARCH_INT[1]); lin <- fread(ARCH_INT[2]); act <- fread(ARCH_INT[3]); imp <- fread(ARCH_INT[4])

  # ---- 8.1 Las w_k(s): qué ventana del año mira cada neurona ----------------------------------
  # w_k(s) = Σ_j a_kj φ_j(s) con la base del satélite. La escala de w_k no importa (z se
  # estandariza), así que se normaliza a max|w_k| = 1. La "ventana" son los días con |w_k| ≥ 0,5.
  segmentos <- function(dias) {
    if (!length(dias)) return("")
    corte <- c(0, which(diff(dias) > 1), length(dias))
    paste(sapply(seq_len(length(corte) - 1), function(k) paste(fecha_dia(dias[corte[k] + 1]), "a", fecha_dia(dias[corte[k + 1]]))), collapse = "; ")
  }
  W_k <- cod[, {
    w <- drop(Phi_s %*% unlist(.SD)); w_n <- w / max(abs(w))
    .(dia = T_DIA, w = w, w_norm = w_n)
  }, by = .(cabeza, semilla, neurona), .SDcols = paste0("a", seq_len(K_s))]
  ventanas <- W_k[, .(dia_max_abs = which.max(abs(w_norm)), signo_en_max = sign(w_norm[which.max(abs(w_norm))]),
                     ventana = segmentos(dia[abs(w_norm) >= .5]), dias_en_ventana = sum(abs(w_norm) >= .5)), by = .(cabeza, semilla, neurona)]
  ventanas[, fecha_max := fecha_dia(dia_max_abs)]
  tabla_f(ventanas, "F8_ventanas_w_k")
  grafico_f(ggplot(W_k[semilla == SEM], aes(dia, w_norm, colour = cabeza)) + geom_hline(yintercept = 0, colour = "grey60") + geom_hline(yintercept = c(-.5, .5), linetype = 3, colour = "grey60") +
              geom_line(linewidth = .7) + facet_wrap(~ paste("neurona", neurona), ncol = 4) +
              scale_x_continuous(breaks = meses$dia, labels = meses$mes) +
              labs(x = "día del año del satélite (s)", y = "w_k(s) normalizada", colour = "cabeza",
                   title = "Pesos del codificador w_k(s): qué parte del año del satélite usa cada neurona (semilla 2026)",
                   subtitle = "Ventana = días con |w_k| ≥ 0,5 (líneas punteadas)"), "F8_pesos_w_k", ancho = 12, alto = 6)

  # ---- 8.2 Cabeza lineal: su β(s,t) frente a la regresión de R y el peso de cada entrada -----
  # La red lineal es la regresión funcional (sección 3 de Python): β(s,t) = Σ_k w_k(s) Σ_m (B_mk/σ_k) ψ_m(t).
  Bl <- dcast(lin[semilla == SEM & entrada != "sesgo"], componente ~ entrada, value.var = "peso")
  ent_nom <- c(paste0("z", 1:8), "altitud", "ONI")
  Bm <- as.matrix(Bl[order(componente), ..ent_nom])                     # M × 10
  cl <- cod[cabeza == "lineal" & semilla == SEM][order(neurona)]
  A_lin <- as.matrix(cl[, paste0("a", seq_len(K_s)), with = FALSE])
  BETA_red <- Phi_s %*% t(A_lin) %*% diag(1 / cl$de_z) %*% t(Bm[, 1:8]) %*% t(psi_t)
  cor_beta <- cor(as.vector(BETA_red), as.vector(BETA))
  # Estabilidad: la β(s,t) de cada semilla frente a la de R y entre semillas. Las neuronas sueltas no
  # son identificables (con cabeza lineal solo importa el espacio que generan las w_k: se pueden
  # rotar sin cambiar la predicción), pero β(s,t), que es lo que predice, sí debe ser estable.
  beta_de <- function(s) {
    B_s <- as.matrix(dcast(lin[semilla == s & entrada != "sesgo"], componente ~ entrada, value.var = "peso")[order(componente), ..ent_nom])
    c_s <- cod[cabeza == "lineal" & semilla == s][order(neurona)]
    Phi_s %*% t(as.matrix(c_s[, paste0("a", seq_len(K_s)), with = FALSE])) %*% diag(1 / c_s$de_z) %*% t(B_s[, 1:8]) %*% t(psi_t)
  }
  betas <- lapply(sort(unique(lin$semilla)), beta_de)
  estab_beta <- data.table(semilla = sort(unique(lin$semilla)), cor_con_regresion_R = sapply(betas, function(b) cor(as.vector(b), as.vector(BETA))),
                           cor_con_semilla_2026 = sapply(betas, function(b) cor(as.vector(b), as.vector(betas[[1]]))))
  print(tabla_f(estab_beta, "F8_estabilidad_beta_lineal"))
  cat(sprintf("β(s,t) de la red lineal (semilla %d) frente a la de la regresión de R: correlación %.3f\n", SEM, cor_beta))
  grafico_base("F8_beta_red_vs_regresion", function() {
    par(mfrow = c(1, 2)); lim <- max(abs(c(BETA, BETA_red))) * c(-1, 1)
    image(T_DIA, T_DIA, BETA, zlim = lim, col = hcl.colors(64, "Blue-Red 3"), xlab = "s: día del satélite", ylab = "t: día de la estación", main = "Regresión funcional (R)")
    image(T_DIA, T_DIA, BETA_red, zlim = lim, col = hcl.colors(64, "Blue-Red 3"), xlab = "s: día del satélite", ylab = "t: día de la estación",
          main = sprintf("Red con cabeza lineal (correlación %.2f)", cor_beta))
  }, ancho = 12, alto = 5)
  # Importancia de cada entrada en la cabeza lineal: norma L² del cambio en la curva √P cuando la
  # entrada sube 1 DE, ||Σ_m B_mi ψ_m(t)||; en la KAN: Σ_j DE de φ_ij(x) sobre entrenamiento.
  imp_lin <- rbindlist(lapply(unique(lin$semilla), function(s) {
    Bs <- as.matrix(dcast(lin[semilla == s & entrada != "sesgo"], componente ~ entrada, value.var = "peso")[order(componente), ..ent_nom])
    data.table(semilla = s, entrada = ent_nom, importancia = sqrt(colSums((psi_t %*% Bs)^2) / 365))
  }))
  imp_kan <- imp[, .(importancia = sum(de_phi_entrenamiento)), by = .(semilla, entrada)]
  importancia <- merge(imp_lin[, .(imp_lineal_media = mean(importancia), imp_lineal_DE = sd(importancia)), by = entrada],
                       imp_kan[, .(imp_kan_media = mean(importancia), imp_kan_DE = sd(importancia)), by = entrada], by = "entrada")
  importancia[, `:=`(rango_lineal = frank(-imp_lineal_media), rango_kan = frank(-imp_kan_media))]
  importancia <- merge(importancia, ventanas[semilla == SEM & cabeza == "lineal", .(entrada = paste0("z", neurona), ventana_lineal = ventana)], by = "entrada", all.x = TRUE)
  importancia <- merge(importancia, ventanas[semilla == SEM & cabeza == "kan", .(entrada = paste0("z", neurona), ventana_kan = ventana)], by = "entrada", all.x = TRUE)
  setorder(importancia, rango_kan)
  print(tabla_f(importancia, "F8_importancia_entradas"))

  # ---- 8.3 Activaciones de la primera capa de la KAN en el rango de los datos ------------------
  # No linealidad de cada arista = 1 − R² de una recta ajustada a φ_ij(x) en el rango p1–p99.
  aristas <- act[, .(no_linealidad = 1 - summary(lm(phi ~ x))$r.squared, rango_phi = diff(range(phi))), by = .(semilla, entrada, oculta)]
  aristas <- merge(aristas, imp, by = c("semilla", "entrada", "oculta"))
  tabla_f(aristas, "F8_kan_aristas")
  a_sem <- aristas[semilla == SEM]
  grafico_f(ggplot(merge(act[semilla == SEM], a_sem[, .(entrada, oculta, de_phi_entrenamiento)], by = c("entrada", "oculta")),
                   aes(x, phi, colour = de_phi_entrenamiento)) + geom_line(linewidth = .8) +
              facet_grid(factor(paste("oculta", oculta)) ~ factor(entrada, levels = ent_nom), scales = "free") +
              scale_colour_viridis_c(name = "importancia\n(DE de φ)") +
              labs(x = "entrada estandarizada (rango p1–p99 de entrenamiento)", y = "φ(x)",
                   title = "Activaciones de la primera capa de la KAN (semilla 2026): cada panel es una arista entrada → neurona oculta"),
            "F8_kan_activaciones", ancho = 15, alto = 7)

  # ---- 8.4 Ajuste simbólico de las activaciones más importantes --------------------------------
  # Candidatas de menor a mayor complejidad; se elige la más simple con R² ≥ 0,99.
  candidatas <- list(
    lineal = list(f = phi ~ x, ini = NULL),
    cuadratica = list(f = phi ~ x + I(x^2), ini = NULL),
    exponencial = list(f = phi ~ a + b * exp(c * x), ini = list(c(a = 0, b = 1, c = .5), c(a = 0, b = -1, c = -.5), c(a = 0, b = 1, c = -.5))),
    tanh = list(f = phi ~ a + b * tanh(c * x + d), ini = list(c(a = 0, b = 1, c = 1, d = 0), c(a = 0, b = -1, c = 1, d = 0), c(a = 0, b = 1, c = 3, d = 1))),
    cubica = list(f = phi ~ x + I(x^2) + I(x^3), ini = NULL),
    seno = list(f = phi ~ a + b * sin(c * x + d), ini = list(c(a = 0, b = 1, c = 1, d = 0), c(a = 0, b = 1, c = 2, d = 1))))
  ajusta <- function(d, nom) {
    cand <- candidatas[[nom]]
    m <- if (is.null(cand$ini)) lm(cand$f, d) else {
      fits <- lapply(cand$ini, function(s0) tryCatch(suppressWarnings(nls(cand$f, d, start = s0, control = nls.control(maxiter = 500, warnOnly = TRUE))), error = function(e) NULL))
      fits <- Filter(Negate(is.null), fits); if (!length(fits)) return(NULL); fits[[which.min(sapply(fits, deviance))]] }
    r2 <- 1 - sum(residuals(m)^2) / sum((d$phi - mean(d$phi))^2)
    data.table(forma = nom, R2 = r2, coeficientes = paste(sprintf("%s = %.3f", names(coef(m)), coef(m)), collapse = ", "))
  }
  top <- a_sem[order(-de_phi_entrenamiento)][1:4]
  simbolico <- rbindlist(lapply(seq_len(nrow(top)), function(k) {
    d <- act[semilla == SEM & entrada == top$entrada[k] & oculta == top$oculta[k]]
    r <- rbindlist(lapply(names(candidatas), function(nm) ajusta(d, nm)))
    r[, `:=`(entrada = top$entrada[k], oculta = top$oculta[k], importancia = top$de_phi_entrenamiento[k], elegida = FALSE)]
    r[which(R2 >= .99)[1], elegida := TRUE]
    if (!any(r$elegida)) r[which.max(R2), elegida := TRUE]
    r
  }))
  setcolorder(simbolico, c("entrada", "oculta", "importancia", "forma", "R2", "elegida"))
  print(tabla_f(simbolico[, R2 := round(R2, 4)], "F8_ajuste_simbolico"))

  # ---- 8.5 Peores casos en prueba ----------------------------------------------------------------
  # Curvas de prueba con mayor error L² relativo medio entre las tres cabezas (promedio de semillas).
  peor <- mc_prom[particion == "prueba" & modelo %in% c("lineal", "kan", "mlp"), .(l2_cabezas = mean(l2_rel)), by = .(codigo, anio, fase)][order(-l2_cabezas)][1:6]
  peor <- merge(peor, est[, .(codigo, nombre, vertiente, altitud)], by = "codigo")[order(-l2_cabezas)]
  peor[, fila := sapply(seq_len(.N), function(k) which(pares$codigo == codigo[k] & pares$anio == anio[k]))]
  peor[, `:=`(P_anual_obs = rowSums(Y_mm[fila, , drop = FALSE]), P_anual_sat = rowSums(X_mm[fila, , drop = FALSE]))]
  curva_cabeza <- function(cab, f) colMeans(do.call(rbind, lapply(split(pr[modelo == cab & fila == f], by = "semilla"), function(d) a_mm(as.matrix(d[, ..cols_c]) %*% t(Phi_e)))))
  peor[, P_anual_kan := sapply(fila, function(f) sum(curva_cabeza("kan", f)))]
  print(tabla_f(peor[, .(codigo, nombre, anio, fase, vertiente, altitud, l2_cabezas = round(l2_cabezas, 2),
                         P_anual_obs = round(P_anual_obs), P_anual_sat = round(P_anual_sat), P_anual_kan = round(P_anual_kan))], "F8_peores_casos"))
  dib <- rbindlist(lapply(seq_len(nrow(peor)), function(k) {
    f <- peor$fila[k]; et <- sprintf("%s %d (%s, %s, %.0f m)", peor$nombre[k], peor$anio[k], peor$fase[k], peor$vertiente[k], peor$altitud[k])
    rbind(data.table(caso = et, dia = T_DIA, serie = "estación", mm = Y_mm[f, ]), data.table(caso = et, dia = T_DIA, serie = "satélite", mm = X_mm[f, ]),
          rbindlist(lapply(c("lineal", "kan", "mlp"), function(cb) data.table(caso = et, dia = T_DIA, serie = cb, mm = curva_cabeza(cb, f)))))
  }))
  grafico_f(ggplot(dib, aes(dia, mm, colour = serie, linewidth = serie == "estación")) + geom_line() + facet_wrap(~ caso, ncol = 2, scales = "free_y") +
              scale_linewidth_manual(values = c(.5, 1.1), guide = "none") + scale_x_continuous(breaks = meses$dia, labels = meses$mes) +
              labs(x = NULL, y = "mm/día", colour = NULL, title = "Los seis peores casos de prueba (mayor error L² relativo de las tres cabezas)"),
            "F8_peores_casos", ancho = 12, alto = 9)
}

# QUÉ DICEN LOS DATOS (sección 8, rúbrica 4.5) — semilla 2026, con la estabilidad entre semillas
#   Las cifras y su interpretación están en EXPLICACION_FUNCIONAL.md, sección 8 (corrida del
#   19-nov-2026, con las 98 estaciones de la regla de la OMM); las tablas de esta
#   sección las imprimen las líneas de arriba.

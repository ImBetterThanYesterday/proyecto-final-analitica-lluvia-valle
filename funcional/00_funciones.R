# =============================================================================
# FUNCIONES COMUNES DE LA PARTE FUNCIONAL (secciones 10 a 14)
# =============================================================================
# Base B-spline (Cox-de Boor), derivada, mínimos cuadrados, criterios de K,
# matriz de Gram, profundidades y boxplot funcional, programados sin librerías.
# Las librerías (fda) solo se usan para verificar.
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
ANIOS_SAT_DEFECTUOSOS <- c(1985L, 1988L, 2000L)

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

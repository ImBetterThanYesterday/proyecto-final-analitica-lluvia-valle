#####################################################################################
# Proyecto final Analítica de Datos 2026-2S
# Fase 2: Análisis exploratorio de datos funcionales (versión en R)
#
# Reproduce en R lo central del notebook 04_eda_funcional.ipynb:
#   - media funcional paso a paso (c_barra' phi) y verificación con fda::mean.fd (si fda está instalado)
#   - varianza, desviación estándar, covarianza Phi Sigma_c Phi' y correlación
#   - derivadas de la base programadas a mano y fechas de las temporadas
#   - profundidad Fraiman-Muniz y MBD (fórmula de rangos) y boxplot funcional (Sun y Genton)
#   - comparación numérica con los resultados de Python
#
# Entradas: fda_fase1/base_especificacion.json, fda_fase1/coef_estaciones.csv,
#           fda_fase2/coef_satelite_tierra.csv.gz, eda_fase0/oni_anual_fase.csv
#####################################################################################

rm(list = ls())
# setwd("ruta/a/proyecto3")
dir_fig <- "fda_fase2/figuras_R"; dir.create(dir_fig, showWarnings = FALSE, recursive = TRUE)

RANGO <- c(1, 365); ORDEN <- 4; K_EST <- 6; K_SAT <- 11     # ver fda_fase1/base_especificacion.json
t <- 1:365

#------------------------------------------------------------------------------------
# Base B-spline y su derivada (Cox-de Boor programado a mano)
#------------------------------------------------------------------------------------
nodos <- function(K, orden = ORDEN) {
  int <- seq(RANGO[1], RANGO[2], length.out = K - orden + 2); int <- int[-c(1, length(int))]
  c(rep(RANGO[1], orden), int, rep(RANGO[2], orden))
}
base_en_nodos <- function(t, tau, m) {
  B <- matrix(0, length(t), length(tau) - 1)
  for (j in 1:(length(tau) - 1)) B[, j] <- as.numeric(t >= tau[j] & t < tau[j + 1])
  B[t == tau[length(tau)], max(which(tau[-length(tau)] < tau[-1]))] <- 1
  if (m >= 2) for (k in 2:m) {
    Bn <- matrix(0, length(t), length(tau) - k)
    for (j in 1:(length(tau) - k)) {
      d1 <- tau[j + k - 1] - tau[j]; d2 <- tau[j + k] - tau[j + 1]
      Bn[, j] <- (if (d1 > 0) (t - tau[j]) / d1 * B[, j] else 0) + (if (d2 > 0) (tau[j + k] - t) / d2 * B[, j + 1] else 0)
    }
    B <- Bn
  }
  B
}
base_bspline <- function(t, K) base_en_nodos(t, nodos(K), ORDEN)
deriv_base <- function(t, K) {
  tau <- nodos(K); B <- base_en_nodos(t, tau, ORDEN - 1); D <- matrix(0, length(t), K)
  for (j in 1:K) {
    d1 <- tau[j + ORDEN - 1] - tau[j]; d2 <- tau[j + ORDEN] - tau[j + 1]
    D[, j] <- (ORDEN - 1) * ((if (d1 > 0) B[, j] / d1 else 0) - (if (d2 > 0) B[, j + 1] / d2 else 0))
  }
  D
}

#------------------------------------------------------------------------------------
# Datos funcionales (coeficientes de la fase 1)
#------------------------------------------------------------------------------------
ce <- read.csv("fda_fase1/coef_estaciones.csv"); C_est <- as.matrix(ce[, -(1:2)])
cs <- read.csv(gzfile("fda_fase2/coef_satelite_tierra.csv.gz")); C_sat <- as.matrix(cs[, -(1:3)])
cuenca <- ifelse(substr(ce$STATION_ID, 1, 2) == "26", "Cauca", "Pacifico")
Phi_e <- base_bspline(t, K_EST); Phi_s <- base_bspline(t, K_SAT)
Ye <- C_est %*% t(Phi_e); Xs <- C_sat %*% t(Phi_s)
cat("Estaciones:", dim(Ye), "| Satélite:", dim(Xs), "\n")

#------------------------------------------------------------------------------------
# 1. Media funcional paso a paso y verificación
#------------------------------------------------------------------------------------
media_e <- as.vector(Phi_e %*% colMeans(C_est)); media_s <- as.vector(Phi_s %*% colMeans(C_sat))
cat("Media paso a paso vs media de curvas evaluadas:", max(abs(media_e - colMeans(Ye))), "\n")
if (requireNamespace("fda", quietly = TRUE)) {
  base_fda <- fda::create.bspline.basis(RANGO, nbasis = K_EST, norder = ORDEN)
  fd_e <- fda::fd(t(C_est), base_fda)
  cat("Media paso a paso vs fda::mean.fd:", max(abs(media_e - fda::eval.fd(t, fda::mean.fd(fd_e)))), "\n")
} else cat("(fda no instalado: la verificación contra mean.fd se hace en el notebook de Python con scikit-fda)\n")

#------------------------------------------------------------------------------------
# 2. Varianza, DE y CV; 3. Covarianza y correlación
#------------------------------------------------------------------------------------
var_e <- colSums(sweep(Ye, 2, media_e)^2) / (nrow(Ye) - 1); de_e <- sqrt(var_e)
var_s <- colSums(sweep(Xs, 2, media_s)^2) / (nrow(Xs) - 1); de_s <- sqrt(var_s)
cov_paso <- function(C, Phi) { Cc <- sweep(C, 2, colMeans(C)); Phi %*% (crossprod(Cc) / (nrow(C) - 1)) %*% t(Phi) }
cov_e <- cov_paso(C_est, Phi_e); cov_s <- cov_paso(C_sat, Phi_s)
cor_e <- cov2cor(cov_e); cor_s <- cov2cor(cov_s)
cat("diag(Cov) - Var:", max(abs(diag(cov_e) - var_e)), "| Cov vs cov() de R:", max(abs(cov_e - cov(Ye))), "\n")
cat(sprintf("Correlación abril-octubre: estaciones %.2f, satélite %.2f\n", cor_e[105, 290], cor_s[105, 290]))

pdf(file.path(dir_fig, "media_de_covarianza.pdf"), width = 12, height = 7)
par(mfrow = c(2, 3))
matplot(t, t(Ye[sample(nrow(Ye), 200), ]), type = "l", lty = 1, col = "grey80", xlab = "Día", ylab = "mm/día", main = "Estaciones: media ± 2 DE")
lines(t, media_e, lwd = 3, col = 4); lines(t, media_e + 2 * de_e, lty = 2, col = 2, lwd = 2); lines(t, media_e - 2 * de_e, lty = 2, col = 2, lwd = 2)
image(t, t, cov_e, col = hcl.colors(30), main = "Covarianza estaciones", xlab = "s", ylab = "t")
image(t, t, cor_e, zlim = c(-1, 1), col = hcl.colors(30, "Blue-Red 3"), main = "Correlación estaciones", xlab = "s", ylab = "t")
matplot(t, t(Xs[sample(nrow(Xs), 200), ]), type = "l", lty = 1, col = "grey80", xlab = "Día", ylab = "mm/día", main = "Satélite: media ± 2 DE")
lines(t, media_s, lwd = 3, col = "darkorange"); lines(t, media_s + 2 * de_s, lty = 2, col = 2, lwd = 2); lines(t, media_s - 2 * de_s, lty = 2, col = 2, lwd = 2)
image(t, t, cov_s, col = hcl.colors(30), main = "Covarianza satélite", xlab = "s", ylab = "t")
image(t, t, cor_s, zlim = c(-1, 1), col = hcl.colors(30, "Blue-Red 3"), main = "Correlación satélite", xlab = "s", ylab = "t")
dev.off()

#------------------------------------------------------------------------------------
# 4. Derivadas y fechas de las temporadas
#------------------------------------------------------------------------------------
DPhi_e <- deriv_base(t, K_EST); DPhi_s <- deriv_base(t, K_SAT)
if (requireNamespace("splines", quietly = TRUE)) {
  ref <- splines::splineDesign(nodos(K_EST), t, ord = ORDEN, derivs = rep(1, length(t)), outer.ok = TRUE)
  cat("Derivada de la base propia vs splines::splineDesign:", max(abs(DPhi_e - ref)), "\n")
}
DYe <- C_est %*% t(DPhi_e); DXs <- C_sat %*% t(DPhi_s)
temporadas <- function(Y, D) {
  t(sapply(seq_len(nrow(Y)), function(i) {
    y <- Y[i, ]; d <- D[i, ]
    p1 <- which(t >= 60 & t <= 200)[which.max(y[t >= 60 & t <= 200])]
    p2 <- which(t >= 230)[which.max(y[t >= 230])]
    c(inicio_1 = which.max(d[1:(p1 - 1)]), pico_1 = p1,
      fin_1 = p1 + which.min(d[(p1 + 1):230]), inicio_2 = 183 + which.max(d[184:(p2 - 1)]),
      pico_2 = p2, fin_2 = if (p2 < 365) p2 + which.min(d[(p2 + 1):365]) else NA)   # pico en el borde: fin no definido
  }))
}
tem_e <- temporadas(Ye, DYe)
print(apply(tem_e, 2, median, na.rm = TRUE))
py_tem <- read.csv("fda_fase2/profundidad_y_fechas_estaciones.csv")[, c("inicio_1", "pico_1", "fin_1", "inicio_2", "pico_2", "fin_2")]
cat("Fechas de temporadas: coincidencia con Python:", mean(tem_e == as.matrix(py_tem), na.rm = TRUE), "| faltantes iguales:", identical(is.na(tem_e), is.na(as.matrix(py_tem))), "\n")

#------------------------------------------------------------------------------------
# 5. Profundidad FM y MBD (fórmula de rangos) y boxplot funcional
#------------------------------------------------------------------------------------
profundidades <- function(Y) {
  N <- nrow(Y); r <- apply(Y, 2, function(v) rank(v, ties.method = "max"))
  list(FM = rowMeans(1 - abs(0.5 - r / N)), MBD = rowMeans(((r - 1) * (N - r) + (N - 1)) / choose(N, 2)))
}
# Verificación de la fórmula de rangos contra la definición por pares (40 curvas)
Yp <- Ye[1:40, ]; pares <- combn(40, 2)
bruta <- sapply(1:40, function(i) mean(apply(pares, 2, function(p)
  mean(pmin(Yp[p[1], ], Yp[p[2], ]) <= Yp[i, ] & Yp[i, ] <= pmax(Yp[p[1], ], Yp[p[2], ])))))
cat("MBD rangos vs fuerza bruta:", max(abs(profundidades(Yp)$MBD - bruta)), "\n")

pr_e <- profundidades(Ye)
py <- read.csv("fda_fase2/profundidad_y_fechas_estaciones.csv")
cat("FM R vs Python:", max(abs(pr_e$FM - py$FM)), "| MBD R vs Python:", max(abs(pr_e$MBD - py$MBD)), "\n")

boxplot_funcional <- function(Y, factor = 1.5) {
  prof <- profundidades(Y)$MBD; centro <- order(-prof)[1:ceiling(nrow(Y) / 2)]
  lo <- apply(Y[centro, , drop = FALSE], 2, min); hi <- apply(Y[centro, , drop = FALSE], 2, max)
  inf <- lo - factor * (hi - lo); sup <- hi + factor * (hi - lo)
  atip <- apply(Y, 1, function(y) any(y < inf | y > sup))
  list(prof = prof, lo = lo, hi = hi, inf = inf, sup = sup, atip = atip, mediana = Y[which.max(prof), ])
}
bx_ec <- boxplot_funcional(Ye[cuenca == "Cauca", ])
at <- ce[cuenca == "Cauca", 1:2][bx_ec$atip, ]
cat("Curvas atípicas (estaciones del Cauca):", nrow(at), "\n"); print(at)
py_at <- read.csv("fda_fase2/atipicas_estaciones_cauca.csv")
cat("¿Mismas atípicas que en Python?", setequal(paste(at$STATION_ID, at$YEAR), paste(py_at$STATION_ID, py_at$YEAR)), "\n")

# Años atípicos con las curvas medias anuales
oni <- read.csv("eda_fase0/oni_anual_fase.csv", fileEncoding = "UTF-8")
anios <- sort(unique(ce$YEAR[cuenca == "Cauca"]))
Ye_anual <- t(sapply(anios, function(a) colMeans(Ye[cuenca == "Cauca" & ce$YEAR == a, , drop = FALSE])))
bx_ya <- boxplot_funcional(Ye_anual)
print(head(merge(data.frame(anio = anios, MBD = round(bx_ya$prof, 3), atipico = bx_ya$atip), oni)[order(bx_ya$prof), ], 6))

pdf(file.path(dir_fig, "boxplot_funcional.pdf"), width = 12, height = 4.5)
par(mfrow = c(1, 2))
dib <- function(Y, bx, titulo) {
  plot(t, bx$mediana, type = "n", ylim = range(c(bx$inf, bx$sup)), xlab = "Día del año", ylab = "mm/día", main = titulo)
  polygon(c(t, rev(t)), c(bx$lo, rev(bx$hi)), col = "#b39ddb", border = NA)
  lines(t, bx$sup, col = 4, lwd = 2); lines(t, bx$inf, col = 4, lwd = 2); lines(t, bx$mediana, lwd = 2)
  if (any(bx$atip)) matlines(t, t(Y[bx$atip, , drop = FALSE]), col = 2, lty = 1)
}
dib(Ye[cuenca == "Cauca", ], bx_ec, paste("Estaciones del Cauca:", sum(bx_ec$atip), "atípicas"))
dib(Ye_anual, bx_ya, "Curvas medias anuales (estaciones del Cauca)")
dev.off()

#------------------------------------------------------------------------------------
# 6. ENSO: lluvia media anual por fase (estaciones del Cauca)
#------------------------------------------------------------------------------------
fase <- oni$fase[match(ce$YEAR, oni$anio)]
print(round(tapply(rowMeans(Ye[cuenca == "Cauca", ]) * 365, fase[cuenca == "Cauca"], mean)))

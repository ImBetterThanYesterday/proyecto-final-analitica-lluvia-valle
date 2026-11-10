#####################################################################################
# Proyecto final Analítica de Datos 2026-2S
# Fase 1: Construcción de datos funcionales paso a paso (versión en R)
#
# Reproduce en R base (sin fda) lo del notebook 03_construccion_funcional.ipynb:
#   - base B-spline por la recursión de Cox-de Boor programada a mano
#   - coeficientes por mínimos cuadrados: c = (Phi'Phi)^-1 Phi' y
#   - verificación de H: simétrica, idempotente y traza = K
#   - selección de K con ECM, VC (atajo y dejar-uno-fuera) y VCG por fuente
#   - sigma^2, covarianza residual y bandas de confianza punto a punto
#   - comparación numérica con los resultados de Python
#
# Entradas (generadas por el notebook 03, ver fda_fase1/):
#   curvas_estaciones_mm.csv.gz : STATION_ID, YEAR, días 1..365 (mm/día, NA = faltante)
#   curvas_satelite_mm.csv.gz   : anio, fila, col, 72 péntadas (mm/día)
#   coef_estaciones.csv, resumen_seleccion_K.csv : resultados de Python para verificar
#####################################################################################

rm(list = ls())
library(splines)          # paquete base de R; solo se usa para VERIFICAR la base propia

# setwd("ruta/a/proyecto3")   # ajustar si es necesario
dir_in  <- "fda_fase1"
dir_fig <- file.path(dir_in, "figuras_R"); dir.create(dir_fig, showWarnings = FALSE)

RANGO <- c(1, 365)
ORDEN <- 4

#------------------------------------------------------------------------------------
# 1. Datos
#------------------------------------------------------------------------------------
est <- read.csv(gzfile(file.path(dir_in, "curvas_estaciones_mm.csv.gz")), check.names = FALSE)
sat <- read.csv(gzfile(file.path(dir_in, "curvas_satelite_mm.csv.gz")), check.names = FALSE)

Y     <- as.matrix(est[, -(1:2)])            # curvas x 365
X     <- as.matrix(sat[, -(1:3)])            # curvas x 72
t_est <- 1:365
t_sat <- as.numeric(sub("t", "", colnames(sat)[-(1:3)]))
completas <- rowSums(is.na(Y)) == 0
cat("Estaciones:", nrow(Y), "curvas (completas:", sum(completas), ") | Satélite:", nrow(X), "curvas\n")

#------------------------------------------------------------------------------------
# 2. Base B-spline paso a paso (Cox-de Boor)
#------------------------------------------------------------------------------------
nodos <- function(rango, K, orden = ORDEN) {
  interiores <- seq(rango[1], rango[2], length.out = K - orden + 2)
  interiores <- interiores[-c(1, length(interiores))]
  c(rep(rango[1], orden), interiores, rep(rango[2], orden))
}

base_bspline <- function(t, rango, K, orden = ORDEN) {
  tau <- nodos(rango, K, orden)
  B <- matrix(0, length(t), length(tau) - 1)
  for (j in 1:(length(tau) - 1)) B[, j] <- as.numeric(t >= tau[j] & t < tau[j + 1])
  ultimo <- max(which(tau[-length(tau)] < tau[-1]))
  B[t == rango[2], ultimo] <- 1                       # incluye el extremo derecho
  for (m in 2:orden) {
    Bn <- matrix(0, length(t), length(tau) - m)
    for (j in 1:(length(tau) - m)) {
      d1 <- tau[j + m - 1] - tau[j]; d2 <- tau[j + m] - tau[j + 1]
      izq <- if (d1 > 0) (t - tau[j]) / d1 * B[, j] else 0
      der <- if (d2 > 0) (tau[j + m] - t) / d2 * B[, j + 1] else 0
      Bn[, j] <- izq + der
    }
    B <- Bn
  }
  B
}

# Verificación contra splines::splineDesign (solo comprobación)
for (K in c(6, 11, 17, 40)) {
  ref <- splineDesign(nodos(RANGO, K), t_est, ord = ORDEN, outer.ok = TRUE)
  cat(sprintf("K=%2d: max|Phi_propia - Phi_splines| = %.1e\n", K, max(abs(base_bspline(t_est, RANGO, K) - ref))))
}

#------------------------------------------------------------------------------------
# 3. Mínimos cuadrados paso a paso y propiedades de H
#------------------------------------------------------------------------------------
verificar_H <- function(H, K, nombre) {
  cat(sprintf("%-22s simétrica %.1e | idempotente %.1e | traza %.10f (K = %d)\n",
              nombre, max(abs(H - t(H))), max(abs(H %*% H - H)), sum(diag(H)), K))
}
Phi  <- base_bspline(t_est, RANGO, 11)
y2c  <- solve(crossprod(Phi), t(Phi))         # C = (Phi'Phi)^-1 Phi'  (y2cMap)
H    <- Phi %*% y2c
verificar_H(H, 11, "Estaciones, K = 11")
i_ej <- which(completas)[1]
cat("Coeficientes vs lsfit:", max(abs(y2c %*% Y[i_ej, ] - lsfit(Phi, Y[i_ej, ], intercept = FALSE)$coef)), "\n")

#------------------------------------------------------------------------------------
# 4. Selección de K: ECM, VC y VCG
#------------------------------------------------------------------------------------
# Atajo de VC vs dejar-uno-fuera (una curva, K = 6)
y1 <- Y[i_ej, ]; P <- base_bspline(t_est, RANGO, 6); H6 <- P %*% solve(crossprod(P), t(P))
vc_atajo <- mean(((y1 - H6 %*% y1) / (1 - diag(H6)))^2)
pred <- sapply(1:365, function(j) P[j, ] %*% solve(crossprod(P[-j, ]), t(P[-j, ]) %*% y1[-j]))
cat(sprintf("VC atajo = %.8f | VC dejar-uno-fuera = %.8f\n", vc_atajo, mean((y1 - pred)^2)))

criterios <- function(Ymat, t, Ks) {
  n <- length(t)
  do.call(rbind, lapply(Ks, function(K) {
    P <- base_bspline(t, RANGO, K); H <- P %*% solve(crossprod(P), t(P)); h <- diag(H)
    E <- Ymat - Ymat %*% t(H)
    data.frame(K = K, ECM = mean(E^2), VC = mean(sweep(E, 2, 1 - h, "/")^2),
               VCG = sum(n * rowSums(E^2) / (n - K)^2))
  }))
}
# Selección de K solo con curvas de entrenamiento (1981-1995): sin fuga de información
ent_est <- completas & est$YEAR <= 1995
ent_sat <- sat$anio <= 1995
sel_est <- criterios(Y[ent_est, ], t_est, c(4:40, seq(45, 100, 5)))
sel_sat <- criterios(X[ent_sat, ], t_sat, 4:40)
K_EST <- sel_est$K[which.min(sel_est$VCG)]; K_SAT <- sel_sat$K[which.min(sel_sat$VCG)]
cat("K por VC  — estaciones:", sel_est$K[which.min(sel_est$VC)], "| satélite:", sel_sat$K[which.min(sel_sat$VC)], "\n")
cat("K por VCG — estaciones:", K_EST, "| satélite:", K_SAT, "\n")

pdf(file.path(dir_fig, "seleccion_K.pdf"), width = 11, height = 6)
par(mfrow = c(2, 3))
for (d in list(list(sel_est, "Estaciones (diario)"), list(sel_sat, "Satélite (péntada)"))) {
  for (crit in c("ECM", "VC", "VCG")) {
    plot(d[[1]]$K, d[[1]][[crit]], type = "b", pch = 19, cex = .6, xlab = "Número de bases (K)", ylab = crit,
         main = paste(d[[2]], "-", crit))
    if (crit != "ECM") abline(v = d[[1]]$K[which.min(d[[1]][[crit]])], lty = 3, col = "red")
  }
}
dev.off()

#------------------------------------------------------------------------------------
# 5. Construcción de todas las curvas con el K elegido
#------------------------------------------------------------------------------------
Phi_est <- base_bspline(t_est, RANGO, K_EST)
coef_est <- t(apply(Y, 1, function(y) {
  m <- !is.na(y); P <- Phi_est[m, ]; solve(crossprod(P), t(P) %*% y[m])
}))
H_est <- Phi_est %*% solve(crossprod(Phi_est), t(Phi_est)); verificar_H(H_est, K_EST, "Estaciones (K elegido)")
# Curva con faltantes: H propia sobre los días observados
i_gap <- which(!completas)[1]; m <- !is.na(Y[i_gap, ]); Pm <- Phi_est[m, ]
verificar_H(Pm %*% solve(crossprod(Pm), t(Pm)), K_EST, "Estación con faltantes")

Phi_sat <- base_bspline(t_sat, RANGO, K_SAT)
C_sat <- solve(crossprod(Phi_sat), t(Phi_sat)); H_sat <- Phi_sat %*% C_sat
verificar_H(H_sat, K_SAT, "Satélite")
coef_sat <- X %*% t(C_sat)

#------------------------------------------------------------------------------------
# 6. Covarianza residual y bandas punto a punto
#------------------------------------------------------------------------------------
Yc <- Y[completas, ]
Res <- Yc - Yc %*% t(H_est)
Cor_res <- cor(Res)
cat(sprintf("Correlación residual media fuera de la diagonal: %.3f | entre días consecutivos: %.3f\n",
            mean(abs(Cor_res[upper.tri(Cor_res)])), mean(diag(Cor_res[-1, -365]))))
v_t <- colMeans(Res^2) / (1 - diag(H_est)); forma <- v_t / mean(v_t)

y <- Y[i_ej, ]; C_est <- solve(crossprod(Phi_est), t(Phi_est))
yhat <- Phi_est %*% C_est %*% y
sigma2 <- sum((y - yhat)^2) / (365 - K_EST)
M <- Phi_est %*% C_est
ee_iid <- sqrt(sigma2 * rowSums(M^2))                           # sigma^2 Phi C C' Phi'
ee_het <- sqrt(rowSums((M %*% diag(sigma2 * forma)) * M))       # Phi C diag(sigma^2(t)) C' Phi'

pdf(file.path(dir_fig, "bandas_y_covarianza.pdf"), width = 11, height = 4)
par(mfrow = c(1, 3))
plot(t_est, y, pch = 19, cex = .3, col = "grey50", xlab = "Día del año", ylab = "mm/día",
     main = paste("Estación", est$STATION_ID[i_ej], est$YEAR[i_ej]), ylim = c(-5, 40))
lines(t_est, yhat, lwd = 2, col = 4)
lines(t_est, yhat + 2 * ee_iid, lty = 2, col = 4); lines(t_est, yhat - 2 * ee_iid, lty = 2, col = 4)
lines(t_est, yhat + 2 * ee_het, lty = 3); lines(t_est, yhat - 2 * ee_het, lty = 3)
legend("topright", c("Curva", "±2 EE (sigma^2 I)", "±2 EE (sigma^2(t))"), lty = c(1, 2, 3), col = c(4, 4, 1), bty = "n")
image(t_est, t_est, Cor_res, zlim = c(-.3, .3), col = hcl.colors(30, "RdBu", rev = TRUE),
      main = "Correlación residual (estaciones)", xlab = "Día", ylab = "Día")
plot(t_est, v_t, type = "l", xlab = "Día del año", ylab = "Varianza", main = "Varianza del error por día")
dev.off()

#------------------------------------------------------------------------------------
# 7. Verificación cruzada con Python
#------------------------------------------------------------------------------------
py_coef <- as.matrix(read.csv(file.path(dir_in, "coef_estaciones.csv"))[, -(1:2)])
py_K    <- read.csv(file.path(dir_in, "resumen_seleccion_K.csv"))
cat("max|coef_R - coef_Python| (estaciones):", format(max(abs(coef_est - py_coef)), digits = 3), "\n")
cat("K elegido por VCG — R:", K_EST, K_SAT, "| Python:", py_K$K.por.VCG[1:2], "\n")

write.csv(data.frame(STATION_ID = est$STATION_ID, YEAR = est$YEAR, coef_est),
          file.path(dir_in, "coef_estaciones_R.csv"), row.names = FALSE)

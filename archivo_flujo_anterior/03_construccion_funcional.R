#####################################################################################
# Proyecto final Analítica de Datos 2026-2S — lluvia del Valle del Cauca (1981–2003)
# Fase 1: Construcción de datos funcionales paso a paso (versión en R)
#
# Hace en R base, sin librerías de datos funcionales, lo mismo que 03_construccion_funcional.ipynb:
#   1. lectura de las curvas preparadas y filtro del periodo 1981–2003
#   2. base B-spline cúbica por la recursión de Cox-de Boor programada a mano
#   3. coeficientes por mínimos cuadrados: c = (Phi'Phi)^-1 Phi' y  (y2cMap de clase)
#   4. selección de K con ECM, VC y VCG por fuente, SOLO con curvas de entrenamiento (1981–1995)
#   5. verificación de S (= H): simétrica, idempotente y traza = K
#   6. sigma^2, covarianza residual y bandas de confianza punto a punto
#   7. matriz de Gram W = integral de phi phi' por cuadratura (para la fase 4)
#   8. exportación de evidencia R; la comparación Python/R corre en tests/validar_salidas_cruzadas.R
#
# Ejecutar de arriba abajo desde la carpeta raíz del repositorio (Source o Ctrl+Shift+Enter).
# Entradas (las genera el notebook 03, carpeta fda_fase1/):
#   curvas_estaciones_mm.csv.gz : STATION_ID, YEAR, días 1..365 (mm/día, NA = faltante)
#   curvas_satelite_mm.csv.gz   : anio, fila, col, 72 péntadas en su día central (mm/día; solo píxeles de tierra)
#   coef_estaciones.csv, resumen_seleccion_K.csv, gram_*.csv : resultados de Python para verificar
# Salidas: fda_fase1/figuras_R/ (PDF), fda_fase1/tablas_R/ (CSV) y los coeficientes que lee 04_eda_funcional.R
#####################################################################################

# setwd("ruta/al/repositorio")          # descomentar si RStudio no abre en la raíz del repositorio
if (!file.exists("fda_fase1/curvas_estaciones_mm.csv.gz"))
  stop("Ejecute desde la raíz del repositorio y corra antes el notebook 03 (genera fda_fase1/curvas_*_mm.csv.gz)")

ANIOS      <- 1981:2003      # periodo del proyecto (fe de erratas)
ENTREN_FIN <- 1995           # entrenamiento 1981–1995: K se elige solo con estos años
RANGO      <- c(1, 365)      # dominio común: día del año (29 de febrero eliminado)
ORDEN      <- 4L             # B-splines cúbicas
DIR        <- "fda_fase1"
DIR_FIG    <- file.path(DIR, "figuras_R"); DIR_TAB <- file.path(DIR, "tablas_R")
dir.create(DIR_FIG, showWarnings = FALSE); dir.create(DIR_TAB, showWarnings = FALSE)
# Notación: Y es N x n (curvas en filas); coef es N x K. mapa = y2cMap = C de clase; H = S.

#------------------------------------------------------------------------------------
# Funciones auxiliares
#------------------------------------------------------------------------------------
guardar_csv <- function(x, nombre) write.csv(x, file.path(DIR_TAB, nombre), row.names = FALSE, na = "")

# Guarda la figura en PDF y, en una sesión interactiva, la dibuja también en la ventana Plots
figura <- function(nombre, dibujar, width = 11, height = 4) {
  pdf(file.path(DIR_FIG, nombre), width = width, height = height); dibujar(); dev.off()
  if (interactive()) dibujar()
  invisible(NULL)
}

mostrar_tabla <- function(x, titulo, filas = 6L) {
  cat("\n", titulo, " (", nrow(x), " filas)\n", sep = ""); print(utils::head(x, filas), row.names = FALSE); invisible(x)
}

mapa_calor <- function(tiempo, z, titulo, limites = range(z)) {
  colores <- hcl.colors(64, "Blue-Red 3")
  layout(matrix(1:2, 1), widths = c(5, 1)); par(mar = c(4, 4, 3, 1))
  image(tiempo, tiempo, z, zlim = limites, col = colores, main = titulo, xlab = "Día s", ylab = "Día t")
  par(mar = c(4, 0, 3, 4)); escala <- seq(limites[1], limites[2], length.out = 65)
  image(c(0, 1), escala, matrix(rep(escala, each = 2), nrow = 2), col = colores, zlim = limites, axes = FALSE, xlab = "", ylab = "")
  axis(4, las = 1, cex.axis = .8); box(); layout(1)
}

leer_matriz <- function(archivo, n_id) {
  datos <- read.csv(gzfile(archivo), check.names = FALSE)
  Y <- as.matrix(datos[, -seq_len(n_id), drop = FALSE]); storage.mode(Y) <- "double"
  claves <- datos[, seq_len(n_id), drop = FALSE]
  if (anyDuplicated(claves)) stop("Claves duplicadas en ", archivo)
  list(id = claves, Y = Y, tiempos = names(datos)[-seq_len(n_id)])
}

# El año se toma por nombre de columna, nunca por posición
anio_de <- function(id) as.integer(if ("YEAR" %in% names(id)) id$YEAR else id$anio)

filtrar_periodo <- function(datos) {
  conservar <- anio_de(datos$id) %in% ANIOS
  datos$id <- datos$id[conservar, , drop = FALSE]; datos$Y <- datos$Y[conservar, , drop = FALSE]
  datos
}

#------------------------------------------------------------------------------------
# Base B-spline paso a paso (Cox-de Boor) — script de clase "Construcción de datos funcionales"
#------------------------------------------------------------------------------------
nodos <- function(K, orden = ORDEN, rango = RANGO) {
  stopifnot(K == as.integer(K), K >= orden)
  interiores <- seq(rango[1], rango[2], length.out = K - orden + 2L)
  c(rep(rango[1], orden), interiores[-c(1L, length(interiores))], rep(rango[2], orden))
}
base_en_nodos <- function(tiempo, tau, orden) {
  B <- outer(tiempo, seq_len(length(tau) - 1L), function(x, j) (x >= tau[j]) & (x < tau[j + 1L])) * 1   # orden 1: indicadoras
  B[tiempo == max(tau), max(which(diff(tau) > 0))] <- 1                                              # incluye el extremo derecho
  if (orden > 1L) for (m in 2:orden) {                                                               # sube el orden
    siguiente <- matrix(0, length(tiempo), length(tau) - m)
    for (j in seq_len(ncol(siguiente))) {
      d1 <- tau[j + m - 1L] - tau[j]; d2 <- tau[j + m] - tau[j + 1L]
      if (d1 > 0) siguiente[, j] <- siguiente[, j] + (tiempo - tau[j]) / d1 * B[, j]
      if (d2 > 0) siguiente[, j] <- siguiente[, j] + (tau[j + m] - tiempo) / d2 * B[, j + 1L]
    }
    B <- siguiente
  }
  B
}
base_bspline <- function(tiempo, K) base_en_nodos(tiempo, nodos(K), ORDEN)

# Mínimos cuadrados sobre los días observados de una curva (sirve también con faltantes)
ajuste_mco <- function(y, Phi) {
  observado <- !is.na(y); P <- Phi[observado, , drop = FALSE]
  if (nrow(P) <= ncol(P) || qr(P)$rank != ncol(P)) stop("Curva no estimable con este K")
  mapa <- solve(crossprod(P), t(P))                     # (Phi'Phi)^-1 Phi'  (y2cMap)
  coef <- as.vector(mapa %*% y[observado]); H <- P %*% mapa
  residuo <- y[observado] - as.vector(P %*% coef)
  list(coef = coef, H = H, mapa = mapa, observado = observado,
       sigma2 = sum(residuo^2) / (nrow(P) - ncol(P)))
}
control_H <- function(H, K)
  c(simetria = max(abs(H - t(H))), idempotencia = max(abs(H %*% H - H)), traza = abs(sum(diag(H)) - K))

# ECM, VC (atajo con diag(H)) y VCG de clase, para curvas completas.
# La SSE por día se acumula con G = Y'Y: da lo mismo que formar todos los residuos y es mucho más rápido.
criterios <- function(Y, tiempo, Ks) {
  G <- crossprod(Y); n <- ncol(Y); N <- nrow(Y)
  do.call(rbind, lapply(Ks, function(K) {
    P <- base_bspline(tiempo, K); mapa <- solve(crossprod(P), t(P)); h <- rowSums(P * t(mapa))   # h = diag(H)
    T <- G %*% t(mapa)
    sse_t <- pmax(0, diag(G) - 2 * rowSums(P * T) + rowSums((P %*% (mapa %*% T)) * P))       # diag((I-H) G (I-H)')
    data.frame(K = K, ECM = sum(sse_t) / (N * n), VC = sum(sse_t / (1 - h)^2) / (N * n),
               VCG = n * sum(sse_t) / (n - K)^2)
  }))
}

#------------------------------------------------------------------------------------
# 1. Datos: leer y filtrar el periodo ANTES de elegir K
#------------------------------------------------------------------------------------
est <- filtrar_periodo(leer_matriz(file.path(DIR, "curvas_estaciones_mm.csv.gz"), 2L))
sat <- filtrar_periodo(leer_matriz(file.path(DIR, "curvas_satelite_mm.csv.gz"), 3L))
t_est <- 1:365
t_sat <- as.numeric(sub("^t", "", sat$tiempos))
stopifnot(identical(as.integer(est$tiempos), 1:365), length(t_sat) == 72L, all(diff(t_sat) > 0))
completas_est <- rowSums(is.na(est$Y)) == 0L
completas_sat <- rowSums(is.na(sat$Y)) == 0L
cat("Estaciones:", nrow(est$Y), "curvas estación-año (completas:", sum(completas_est), ") | Satélite (tierra):",
    nrow(sat$Y), "curvas píxel-año (completas:", sum(completas_sat), ")\n")

cobertura <- do.call(rbind, lapply(ANIOS, function(a) {
  k <- anio_de(est$id) == a
  data.frame(anio = a, conjunto = ifelse(a <= ENTREN_FIN, "entrenamiento", ifelse(a <= 2000, "validación", "prueba")),
             curvas_estaciones = sum(k), completas = sum(k & completas_est), curvas_satelite = sum(anio_de(sat$id) == a))
}))
guardar_csv(cobertura, "cobertura_por_anio.csv")
mostrar_tabla(cobertura, "Curvas por año y conjunto de la partición", 23L)

# Mirar los datos discretos antes de suavizarlos
figura("observados.pdf", function() {
  par(mfrow = c(1, 2))
  i <- which(completas_est)[1]
  plot(t_est, est$Y[i, ], type = "h", col = "grey50", xlab = "Día del año", ylab = "mm/día",
       main = paste("Estación", est$id$STATION_ID[i], est$id$YEAR[i], "(diario)"))
  plot(t_sat, sat$Y[1, ], type = "b", pch = 16, cex = .6, col = "darkorange", xlab = "Día del año", ylab = "mm/día",
       main = paste("Satélite", sat$id$anio[1], "píxel", sat$id$fila[1], sat$id$col[1], "(péntada en su día central)"))
})

#------------------------------------------------------------------------------------
# 2. Controles algebraicos de la base manual (sin usar otra implementación)
#------------------------------------------------------------------------------------
for (K in c(6, 11, 17, 40)) {
  Phi_control <- base_bspline(t_est, K)
  suma_filas <- rowSums(Phi_control)
  stopifnot(all(is.finite(Phi_control)), min(Phi_control) >= -1e-12,
            max(abs(suma_filas - 1)) < 1e-10)
  cat(sprintf("K = %2d: Phi manual finita/no negativa; suma de filas entre %.12f y %.12f\n",
              K, min(suma_filas), max(suma_filas)))
}

#------------------------------------------------------------------------------------
# 3. Mínimos cuadrados paso a paso y atajo de la VC
#------------------------------------------------------------------------------------
i_ej <- which(completas_est)[1]; y1 <- est$Y[i_ej, ]
P17 <- base_bspline(t_est, 17); a17 <- ajuste_mco(y1, P17)
print(control_H(a17$H, 17))
vc_atajo <- mean(((y1 - a17$H %*% y1) / (1 - diag(a17$H)))^2)
pred <- sapply(1:365, function(j) P17[j, ] %*% solve(crossprod(P17[-j, ]), crossprod(P17[-j, ], y1[-j])))
cat(sprintf("VC atajo = %.8f | VC dejar-uno-fuera (365 ajustes) = %.8f\n", vc_atajo, mean((y1 - pred)^2)))

#------------------------------------------------------------------------------------
# 4. Selección de K: ECM, VC y VCG solo con curvas de entrenamiento (sin fuga de información)
#------------------------------------------------------------------------------------
ent_est <- completas_est & anio_de(est$id) <= ENTREN_FIN
ent_sat <- completas_sat & anio_de(sat$id) <= ENTREN_FIN
cat("Curvas de entrenamiento para elegir K — estaciones completas:", sum(ent_est), "| satélite:", sum(ent_sat), "\n")
sel_est <- criterios(est$Y[ent_est, ], t_est, c(4:40, seq(45, 100, 5)))
sel_sat <- criterios(sat$Y[ent_sat, ], t_sat, 4:40)
K_EST <- sel_est$K[which.min(sel_est$VCG)]
K_SAT <- sel_sat$K[which.min(sel_sat$VCG)]
resumen_K <- data.frame(fuente = c("estaciones", "satelite"),
  K_VC = c(sel_est$K[which.min(sel_est$VC)], sel_sat$K[which.min(sel_sat$VC)]), K_VCG = c(K_EST, K_SAT))
guardar_csv(sel_est, "seleccion_K_estaciones.csv"); guardar_csv(sel_sat, "seleccion_K_satelite.csv")
mostrar_tabla(resumen_K, "K elegido (se usa el de la VCG)")
cat("VC de estaciones a menos de 1 % de su mínimo para K entre",
    paste(range(sel_est$K[sel_est$VC <= 1.01 * min(sel_est$VC)]), collapse = " y "), "\n")

figura("seleccion_K.pdf", function() {
  par(mfrow = c(2, 3))
  for (d in list(list(sel_est, "Estaciones (diario)", ""), list(sel_sat, "Satélite (péntada)", "y"))) {
    for (crit in c("ECM", "VC", "VCG")) {
      plot(d[[1]]$K, d[[1]][[crit]], type = "b", pch = 16, cex = .6, col = 4, xlab = "Número de bases K", ylab = crit,
           log = if (crit != "ECM") d[[3]] else "", main = paste(d[[2]], "-", crit))
      if (crit != "ECM") {
        j <- which.min(d[[1]][[crit]]); abline(v = d[[1]]$K[j], col = 2, lty = 2)
        legend("topright", paste("Mínimo: K =", d[[1]]$K[j]), col = 2, lty = 2, bty = "n")
      }
    }
  }
}, height = 7)

#------------------------------------------------------------------------------------
# 5. Construcción de todas las curvas con el K elegido y verificación de S
#------------------------------------------------------------------------------------
construir_fuente <- function(datos, tiempo, K, nombre) {
  Y <- datos$Y; completas <- rowSums(is.na(Y)) == 0L
  Phi <- base_bspline(tiempo, K)
  C <- matrix(NA_real_, nrow(Y), K); sigma2 <- rep(NA_real_, nrow(Y))
  mapa <- solve(crossprod(Phi), t(Phi)); H <- Phi %*% mapa            # curvas completas: una sola H
  C[completas, ] <- Y[completas, , drop = FALSE] %*% t(mapa)
  Res <- Y[completas, , drop = FALSE] - C[completas, , drop = FALSE] %*% t(Phi)
  sigma2[completas] <- rowSums(Res^2) / (length(tiempo) - K)
  controles <- control_H(H, K)
  for (i in which(!completas)) {                                      # curvas con faltantes: una H por curva
    a <- ajuste_mco(Y[i, ], Phi); C[i, ] <- a$coef; sigma2[i] <- a$sigma2
    controles <- pmax(controles, control_H(a$H, K))
  }
  stopifnot(all(controles < 1e-8), all(is.finite(C)))
  colnames(C) <- paste0("c", seq_len(K))
  mostrar_tabla(data.frame(control = names(controles), error_max = controles), paste(nombre, ": controles de S (peor caso)"))
  list(id = datos$id, coef = C, K = K, Phi = Phi, mapa = mapa, H = H, sigma2 = sigma2, Res = Res,
       completas = completas, n_observado = rowSums(!is.na(Y)))
}
e <- construir_fuente(est, t_est, K_EST, "Estaciones")
s <- construir_fuente(sat, t_sat, K_SAT, "Satélite")

figura("bases.pdf", function() {
  par(mfrow = c(1, 2))
  for (K in c(K_EST, K_SAT)) {
    matplot(t_est, base_bspline(t_est, K), type = "l", lty = 1, xlab = "Día del año", ylab = "phi_k(t)",
            main = paste("Base B-spline cúbica, K =", K))
    abline(v = unique(nodos(K)), lty = 3, col = "grey60")
  }
})

# Coeficientes de R: los usa 04_eda_funcional.R (el EDA en R no depende de los coeficientes de Python)
write.csv(cbind(e$id, e$coef), file.path(DIR, "coef_estaciones_R.csv"), row.names = FALSE)
write.csv(cbind(s$id, s$coef), gzfile(file.path(DIR, "coef_satelite_R.csv.gz")), row.names = FALSE)
guardar_csv(cbind(e$id, sigma2 = e$sigma2, n_observado = e$n_observado), "sigma2_estaciones.csv")

# Curvas evaluadas en los 365 días y valores negativos (se recortan en 0 en la fase 4)
F_est <- e$coef %*% t(e$Phi); F_sat <- s$coef %*% t(base_bspline(t_est, K_SAT))
negativos <- data.frame(fuente = c("estaciones", "satelite"),
  pct_puntos_negativos = 100 * c(mean(F_est < 0), mean(F_sat < 0)),
  pct_curvas_con_negativos = 100 * c(mean(apply(F_est < 0, 1, any)), mean(apply(F_sat < 0, 1, any))),
  minimo_mm_dia = c(min(F_est), min(F_sat)))
guardar_csv(negativos, "valores_negativos.csv"); mostrar_tabla(negativos, "Valores negativos de las curvas (mm/día)")

#------------------------------------------------------------------------------------
# 6. Covarianza residual y bandas de confianza punto a punto
#------------------------------------------------------------------------------------
# Momento residual con las curvas completas: N^-1 sum e_i e_i'
Sigma_res_est <- crossprod(e$Res) / nrow(e$Res); Sigma_res_sat <- crossprod(s$Res) / nrow(s$Res)
Cor_res <- cov2cor(Sigma_res_est); Cor_res_sat <- cov2cor(Sigma_res_sat)
cat(sprintf("|Correlación residual| media fuera de la diagonal — estaciones %.3f (lag 1 día %.3f) | satélite %.3f\n",
            mean(abs(Cor_res[upper.tri(Cor_res)])), mean(diag(Cor_res[-1, -365])), mean(abs(Cor_res_sat[upper.tri(Cor_res_sat)]))))
v_t <- diag(Sigma_res_est) / (1 - diag(e$H)); forma <- v_t / mean(v_t)   # forma estacional de la varianza del error
cat(sprintf("Varianza del error (estaciones): mínima el día %d, máxima el día %d\n", which.min(v_t), which.max(v_t)))
guardar_csv(data.frame(dia = t_est, varianza_error = v_t), "varianza_error_estaciones.csv")

# EE punto a punto: sqrt(diag(Phi C Sigma C' Phi')), con C restringida a los días observados
ee_banda <- function(Phi_eval, mapa, Sigma) { M <- Phi_eval %*% mapa; sqrt(rowSums((M %*% Sigma) * M)) }

i <- i_ej; y <- est$Y[i, ]; curva <- F_est[i, ]
ee_iid <- ee_banda(e$Phi, e$mapa, e$sigma2[i] * diag(365))
ee_het <- ee_banda(e$Phi, e$mapa, diag(e$sigma2[i] * forma))
# Curva con el hueco consecutivo más largo: su banda usa la C de los días observados
racha <- function(m) { r <- rle(m); max(c(0, r$lengths[r$values])) }
inc <- which(!e$completas); ig <- inc[which.max(sapply(inc, function(k) racha(is.na(est$Y[k, ]))))]
ag <- ajuste_mco(est$Y[ig, ], e$Phi)
ee_g <- ee_banda(e$Phi, ag$mapa, ag$sigma2 * diag(sum(ag$observado)))
ee_g_completa <- ee_banda(e$Phi, e$mapa, ag$sigma2 * diag(365))
cat(sprintf("Curva con hueco de %d días: EE con hueco / EE sin hueco = %.2f en los días faltantes\n",
            sum(!ag$observado), mean((ee_g / ee_g_completa)[!ag$observado])))
j <- 20000L; ee_s <- ee_banda(base_bspline(t_est, K_SAT), s$mapa, s$sigma2[j] * diag(72))
guardar_csv(data.frame(est$id[rep(i, 365), ], dia = t_est, curva = curva, ee_iid = ee_iid, ee_het = ee_het,
                       inferior = curva - 2 * ee_iid, superior = curva + 2 * ee_iid), "banda_ejemplo_estacion.csv")

figura("bandas.pdf", function() {
  par(mfrow = c(1, 3))
  plot(t_est, y, pch = 16, cex = .3, col = "grey50", ylim = c(-5, 40), xlab = "Día del año", ylab = "mm/día",
       main = paste("Estación", est$id$STATION_ID[i], est$id$YEAR[i]))
  lines(t_est, curva, lwd = 2, col = 4)
  lines(t_est, curva + 2 * ee_iid, lty = 2, col = 4); lines(t_est, curva - 2 * ee_iid, lty = 2, col = 4)
  lines(t_est, curva + 2 * ee_het, lty = 3); lines(t_est, curva - 2 * ee_het, lty = 3)
  legend("topright", c("Curva", "±2 EE (sigma^2 I)", "±2 EE (sigma^2(t))"), lty = 1:3, col = c(4, 4, 1), bty = "n", cex = .8)
  cs <- F_sat[j, ]
  plot(t_sat, sat$Y[j, ], pch = 16, cex = .6, col = "grey50", ylim = range(c(sat$Y[j, ], cs - 2 * ee_s, cs + 2 * ee_s)),
       xlab = "Día del año", ylab = "mm/día", main = paste("Satélite", sat$id$anio[j], "píxel", sat$id$fila[j], sat$id$col[j]))
  lines(t_est, cs, lwd = 2, col = "darkorange"); lines(t_est, cs + 2 * ee_s, lty = 2); lines(t_est, cs - 2 * ee_s, lty = 2)
  cg <- F_est[ig, ]
  plot(t_est, est$Y[ig, ], pch = 16, cex = .3, col = "grey50", ylim = c(-5, 40), xlab = "Día del año", ylab = "mm/día",
       main = paste("Estación", est$id$STATION_ID[ig], est$id$YEAR[ig], "-", sum(!ag$observado), "días faltantes"))
  lines(t_est, cg, lwd = 2, col = 4); lines(t_est, cg + 2 * ee_g, lty = 2, col = 4); lines(t_est, cg - 2 * ee_g, lty = 2, col = 4)
  rug(t_est[!ag$observado], col = 2, lwd = 2)
})
figura("covarianza_residual.pdf", function() {
  par(mfrow = c(1, 2))
  image(t_est, t_est, Cor_res, zlim = c(-.3, .3), col = hcl.colors(30, "RdBu", rev = TRUE),
        main = "Correlación residual (estaciones)", xlab = "Día", ylab = "Día")
  plot(t_est, v_t, type = "l", xlab = "Día del año", ylab = "(mm/día)^2", main = "Varianza del error por día (estaciones)")
})
figura("momento_residual_satelite.pdf", function() mapa_calor(t_sat, Sigma_res_sat, "Momento residual: satélite"),
       width = 7, height = 5)

# Compromiso sesgo-varianza: ancho medio de la banda sigma^2 I vs K (estaciones completas)
Yc <- est$Y[e$completas, ]
anchos <- sapply(c(6, 8, 11, 14, 17, 20, 25, 30, 40, 60), function(K) {
  P <- base_bspline(t_est, K); H <- P %*% solve(crossprod(P), t(P))
  s2 <- rowSums((Yc - Yc %*% t(H))^2) / (365 - K); 4 * mean(sqrt(mean(diag(H)) * s2))
})
figura("ancho_banda_vs_K.pdf", function() {
  plot(c(6, 8, 11, 14, 17, 20, 25, 30, 40, 60), anchos, type = "b", pch = 16, xlab = "K", ylab = "mm/día",
       main = "Ancho medio de la banda vs K (estaciones)"); abline(v = K_EST, col = 2, lty = 2)
}, width = 6, height = 4)

#------------------------------------------------------------------------------------
# 7. Matriz de Gram W = integral de phi(t) phi(t)' dt (codificador y pérdida de la fase 4)
#------------------------------------------------------------------------------------
# Gauss-Legendre de 4 puntos en cada intervalo entre nodos: phi_j phi_k es un polinomio de grado 6
# por tramo y 4 puntos integran exactamente hasta grado 7. Nodos y pesos de Gauss-Legendre (n = 4):
GL_x <- c(-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053)
GL_w <- c(0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454)
gram <- function(K) {
  tau <- unique(nodos(K)); a <- head(tau, -1); b <- tail(tau, -1)
  tq <- as.vector(outer(GL_x, seq_along(a), function(x, k) (b[k] - a[k]) / 2 * x + (a[k] + b[k]) / 2))
  wq <- as.vector(outer(GL_w, seq_along(a), function(w, k) (b[k] - a[k]) / 2 * w))
  P <- base_bspline(tq, K); crossprod(P, wq * P)
}
W_est <- gram(K_EST); W_sat <- gram(K_SAT)
for (d in list(list("estaciones", W_est), list("satelite", W_sat))) {
  W <- d[[2]]
  controles_W <- c(simetria = max(abs(W - t(W))), autovalor_minimo = min(eigen(W, symmetric = TRUE)$values),
                   integral_constante = abs(sum(W) - 364))
  stopifnot(controles_W["simetria"] < 1e-10, controles_W["autovalor_minimo"] > 0,
            controles_W["integral_constante"] < 1e-8)
  cat(sprintf("W %s: simétrica %.1e | autovalor mínimo %.3f | 1'W1 = %.8f (debe ser 364)\n",
              d[[1]], controles_W["simetria"], controles_W["autovalor_minimo"], sum(W)))
}
write.table(W_est, file.path(DIR, "gram_estaciones_R.csv"), sep = ",", row.names = FALSE, col.names = FALSE)
write.table(W_sat, file.path(DIR, "gram_satelite_R.csv"), sep = ",", row.names = FALSE, col.names = FALSE)
t_fino <- seq(1, 365, by = .01); x_fino <- as.vector(base_bspline(t_fino, K_EST) %*% e$coef[i_ej, ])
cat(sprintf("c'Wc = %.4f | integral de x(t)^2 por trapecio fino = %.4f\n", drop(e$coef[i_ej, ] %*% W_est %*% e$coef[i_ej, ]),
            sum(diff(t_fino) * (head(x_fino, -1)^2 + tail(x_fino, -1)^2) / 2)))

cat("Construcción manual en R terminada. Para comparar coeficientes, K y W con Python, ejecutar:",
    "Rscript tests/validar_salidas_cruzadas.R\n")
cat("Figuras en", DIR_FIG, "y tablas en", DIR_TAB, "\n")

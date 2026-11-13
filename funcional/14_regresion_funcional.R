# =============================================================================
# 14. REGRESIÓN FUNCIONAL FUNCIÓN A FUNCIÓN (SATÉLITE -> ESTACIÓN)
# =============================================================================
# Modelo: y_i(t) = mu(t) + alfa(t) + integral beta(s,t) x_i(s) ds
#                  + gamma(t) altitud_i + delta(t) ONI_i + e_i(t)
# sobre la raíz de la lluvia (curvas recortadas en 0 y luego con raíz, D-23).
# beta(s,t) = sum_m beta_m(s) psi_m(t): psi_m son las funciones propias del
# FPCA de las estaciones (solo entrenamiento) y beta_m(s) vive en la base
# B-spline del satélite con penalización de rugosidad lambda * integral beta_m''^2.
# Es el modelo de rango bajo al que equivale la cabeza lineal de la red
# funcional. lambda se elige con validación (1996-2000) y se evalúa en prueba
# (2001-2003) frente al satélite sin corregir. Los años defectuosos del
# satélite (1985, 1988, 2000) no entran.
# Entradas: salidas_R_funcional/construccion.rds, comparacion.rds
# Salidas:  salidas_R_funcional/regresion.rds y salidas_R_funcional/para_python/
# =============================================================================

source("funcional/00_funciones.R")
titulo("14. Regresión funcional")
set.seed(20261113)

o <- readRDS(file.path(SALIDA_F, "construccion.rds"))
pares <- readRDS(file.path(SALIDA_F, "comparacion.rds"))$pares
K_e <- o$K_est; K_s <- o$K_sat
Phi_e <- base_bspline(T_DIA, K_e); Phi_s <- base_bspline(T_DIA, K_s)
W_e <- o$W_est; W_s <- o$W_sat
R_s <- gram_derivada(K_s, 2L)
ref_d2 <- splines::splineDesign(nodos(K_s), T_DIA, ord = ORDEN, derivs = rep(2, 365), outer.ok = TRUE)
cat("Segunda derivada de la base propia vs splines::splineDesign:", format(max(abs(derivada_base(T_DIA, K_s, 2L) - ref_d2)), digits = 2), "\n")

titulo("14.1 Escala raíz")
a_raiz <- function(C, Phi) {
  F <- C %*% t(Phi); G <- sqrt(pmax(F, 0)); Cr <- G %*% t(solve(crossprod(Phi), t(Phi)))
  list(coef = Cr, pct_recortado = 100 * mean(F < 0), rmse_reproyeccion = sqrt(mean((Cr %*% t(Phi) - G)^2)))
}
re <- a_raiz(o$coef_est[pares$fila_est, ], Phi_e); rs <- a_raiz(o$coef_sat[pares$fila_sat, ], Phi_s)
print(data.table(fuente = c("estaciones", "satélite"), pct_dias_recortados_en_0 = round(c(re$pct_recortado, rs$pct_recortado), 2),
                 rmse_reproyeccion_raiz = signif(c(re$rmse_reproyeccion, rs$rmse_reproyeccion), 3)))
Y_mm <- pmax(o$coef_est[pares$fila_est, ] %*% t(Phi_e), 0)
X_mm <- pmax(o$coef_sat[pares$fila_sat, ] %*% t(Phi_s), 0)

pares[, usar := !defectuoso & !is.na(altitud)]
cat("Pares fuera del modelo: años defectuosos", pares[defectuoso == TRUE, .N], "| sin altitud oficial (Univalle)", pares[!defectuoso & is.na(altitud), .N], "\n")
print(pares[usar == TRUE, .(pares = .N, estaciones = uniqueN(codigo), anios = paste(range(anio), collapse = "-")), by = particion])
ent <- pares$usar & pares$particion == "entrenamiento"
val <- pares$usar & pares$particion == "validación"
pru <- pares$usar & pares$particion == "prueba"

titulo("14.2 FPCA de las estaciones (solo entrenamiento)")
ev <- eigen(W_e, symmetric = TRUE); W12 <- ev$vectors %*% diag(sqrt(ev$values)) %*% t(ev$vectors); Wm12 <- ev$vectors %*% diag(1 / sqrt(ev$values)) %*% t(ev$vectors)
mu_c <- colMeans(re$coef[ent, ])
Cc <- sweep(re$coef[ent, ], 2, mu_c)
eg <- eigen(W12 %*% (crossprod(Cc) / (sum(ent) - 1)) %*% W12, symmetric = TRUE)
psi_c <- Wm12 %*% eg$vectors
var_fpca <- data.table(m = seq_len(K_e), autovalor = eg$values, pct = 100 * eg$values / sum(eg$values))[, acumulado := cumsum(pct)]
print(tabla_f(var_fpca, "14_fpca_varianza"))
stopifnot(max(abs(t(psi_c) %*% W_e %*% psi_c - diag(K_e))) < 1e-10)
fd_r <- fda::fd(t(re$coef[ent, ]), fda::create.bspline.basis(RANGO, nbasis = K_e, norder = ORDEN))
pca_fda <- fda::pca.fd(fd_r, nharm = 4, centerfns = TRUE)
psi_fda <- fda::eval.fd(T_DIA, pca_fda$harmonics)
cat(sprintf("Contra fda::pca.fd (4 primeras): autovalores x N/(N-1) a %.1e relativo, funciones propias a %.1e (salvo signo); fda integra W numéricamente\n",
            max(abs(pca_fda$values[1:4] * sum(ent) / (sum(ent) - 1) - eg$values[1:4]) / eg$values[1:4]),
            max(abs(abs(psi_fda) - abs(Phi_e %*% psi_c[, 1:4])))))
S_todos <- sweep(re$coef, 2, mu_c) %*% W_e %*% psi_c

titulo("14.3 Ajuste y elección de M y lambda con validación")
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
tabla_f(sel_lambda, "14_seleccion_M_lambda")
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
          "14_seleccion_M_lambda", ancho = 9, alto = 5)

B <- ajustar_B(LAMBDA, M); B_sin_sat <- ajustar_B(0, M, con_sat = FALSE)
gl_df <- function(Z, P, l) sum(diag(solve(crossprod(Z) + l * P, crossprod(Z))))
cat(sprintf("Grados de libertad efectivos por componente: %.1f de %d parámetros\n", gl_df(diseno(which(ent)), penal(), LAMBDA), ncol(diseno(1))))

titulo("14.4 Superficie beta(s,t) y verificación de la equivalencia")
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

grafico_base("14_fpca", function() {
  par(mfrow = c(1, 2))
  plot(T_DIA, drop(Phi_e %*% mu_c), type = "l", lwd = 2, xlab = "Día del año", ylab = "raíz(mm/día)", main = "Media de entrenamiento (escala raíz)")
  matplot(T_DIA, psi_t, type = "l", lty = 1, lwd = 2, col = 1:M, xlab = "Día del año", ylab = expression(psi[m](t)),
          main = paste("Funciones propias:", paste(sprintf("%.0f %%", var_fpca$pct[1:M]), collapse = ", ")))
  abline(h = 0, lty = 3)
}, ancho = 12, alto = 4.5)
grafico_base("14_beta_superficie", function() {
  par(mfrow = c(1, 2))
  sem <- seq(1, 365, by = 7)
  persp(sem, sem, BETA[sem, sem], theta = 35, phi = 30, col = "lightblue", border = "grey50", xlab = "s (satélite)", ylab = "t (estación)",
        zlab = "beta(s,t)", main = expression(hat(beta)(s, t)), ticktype = "simple")
  lim <- max(abs(BETA)) * c(-1, 1)
  image(T_DIA, T_DIA, BETA, zlim = lim, col = hcl.colors(64, "Blue-Red 3"), xlab = "s: día del satélite", ylab = "t: día de la estación", main = "Curvas de nivel")
  contour(T_DIA, T_DIA, BETA, add = TRUE, nlevels = 10, drawlabels = FALSE, col = "grey30"); abline(0, 1, lty = 2)
}, ancho = 13, alto = 5.5)
grafico_base("14_efectos", function() {
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

titulo("14.5 Evaluación en prueba (2001-2003)")
kge <- function(P, Yo) {
  r <- sapply(seq_len(nrow(P)), function(k) cor(P[k, ], Yo[k, ]))
  1 - sqrt((r - 1)^2 + (apply(P, 1, sd) / apply(Yo, 1, sd) - 1)^2 + (rowMeans(P) / rowMeans(Yo) - 1)^2)
}
pico <- function(M, ventana) apply(M[, ventana, drop = FALSE], 1, which.max) + min(ventana) - 1
modelos <- list(`Satélite sin corregir` = X_mm[i_pru, ], `Regresión funcional` = a_mm(pred_puntajes),
                `Solo altitud y ONI` = a_mm(predecir_raiz(B_sin_sat, i_pru, FALSE)))
recorte_pred <- 100 * mean(pred_puntajes < 0)
metricas_curva <- rbindlist(lapply(names(modelos), function(nm) {
  P <- modelos[[nm]]; Yo <- Y_mm[i_pru, ]
  data.table(modelo = nm, codigo = pares$codigo[i_pru], anio = pares$anio[i_pru], fase = pares$fase[i_pru],
             rmse = rmse_l2(P, i_pru), sesgo_anual = (rowMeans(P) - rowMeans(Yo)) * 365, kge = kge(P, Yo),
             err_pico_1 = abs(pico(P, 60:200) - pico(Yo, 60:200)), err_pico_2 = abs(pico(P, 230:350) - pico(Yo, 230:350)))
}))
boot_media <- function(v, g, B = 1000) {
  m <- tapply(v, g, mean); q <- quantile(replicate(B, mean(sample(m, replace = TRUE))), c(.025, .975))
  sprintf("%.2f [%.2f, %.2f]", mean(m), q[1], q[2])
}
metricas <- metricas_curva[, .(curvas = .N, estaciones = uniqueN(codigo), RMSE_L2 = boot_media(rmse, codigo), KGE = boot_media(kge, codigo),
                               sesgo_anual_mm = boot_media(sesgo_anual, codigo), error_pico_1_dias = boot_media(err_pico_1, codigo),
                               error_pico_2_dias = boot_media(err_pico_2, codigo)), by = modelo]
print(tabla_f(metricas, "14_metricas_prueba"))
cat(sprintf("Reconstrucción en mm no negativa: %.2f %% de los valores predichos en la escala raíz se recortan en 0\n", recorte_pred))
por_fase <- metricas_curva[, .(curvas = .N, RMSE_L2 = round(mean(rmse), 2), KGE = round(median(kge), 2), sesgo_anual_mm = round(mean(sesgo_anual))), by = .(modelo, anio, fase)]
print(tabla_f(por_fase[order(anio, modelo)], "14_metricas_prueba_por_anio_ENSO"))
mejora <- dcast(metricas_curva[modelo != "Solo altitud y ONI"], codigo + anio ~ modelo, value.var = "rmse")
cat(sprintf("La regresión baja el RMSE frente al satélite sin corregir en %.0f %% de las curvas de prueba\n",
            100 * mean(mejora$`Regresión funcional` < mejora$`Satélite sin corregir`)))

R2_t <- function(P, Yo) 1 - colSums((Yo - P)^2) / colSums(sweep(Yo, 2, colMeans(Yo))^2)
grafico_base("14_prueba", function() {
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
  boxplot(rmse ~ modelo, metricas_curva, ylab = "RMSE L2 (mm/día)", xlab = "", main = "Error por curva (prueba)", cex.axis = .8)
  boxplot(sesgo_anual ~ modelo, metricas_curva, ylab = "mm/año", xlab = "", main = "Sesgo anual (prueba)", cex.axis = .8); abline(h = 0, lty = 3)
}, ancho = 15, alto = 8)

titulo("14.6 Exportación para la red funcional (Python)")
dir_py <- file.path(SALIDA_F, "para_python"); dir.create(dir_py, showWarnings = FALSE)
fwrite(pares[, .(codigo, anio, particion, defectuoso, usar, celda, altitud, ONI_anual, fase, vertiente)], file.path(dir_py, "pares.csv"))
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

writeLines(capture.output(sessionInfo()), file.path(SALIDA_F, "sessionInfo.txt"))

# QUÉ DICEN LOS DATOS
#   Escala: al recortar en 0 solo cambia el 2 % de los días de las estaciones y
#   casi ningún día del satélite; la raíz se reproyecta en la misma base con un
#   error de 0,09 (estaciones) y 0,03 (satélite) en raíz de mm/día.
#   Datos del modelo: 2 863 pares sin años defectuosos y con altitud oficial
#   (1 949 de entrenamiento, 542 de validación, 372 de prueba); sale Univalle, la
#   única incluida sin altitud oficial (D-35), con 17 pares.
#   FPCA (solo entrenamiento): la primera componente explica 61 % (el nivel de
#   lluvia de cada sitio), la segunda 13 % y la tercera 8 %. Con K = 8 las
#   últimas componentes tienen autovalores casi iguales (6,4; 6,0; 5,8), es decir,
#   ruido. M y lambda se eligen juntos con validación: el error deja de bajar
#   desde M = 5 (91 % de la varianza), con lambda = 1e8 en la escala de
#   integral beta''^2 (10,9 grados de libertad efectivos de 14).
#   Equivalencia: la predicción con la integral de beta(s,t) por cuadratura y la
#   predicción por puntajes FPCA difieren en 5e-15. Es la identidad
#   beta(s,t) = sum_m beta_m(s) psi_m(t) que debe reproducir la cabeza lineal de
#   la red.
#   Interpretación: la integral de beta sobre s va de 0,57 a 1,06. El modelo
#   traslada el nivel del satélite casi completo en los picos y lo amortigua en
#   marzo y diciembre, cuando el satélite arranca tarde las temporadas. El ONI
#   resta lluvia todo el año y sobre todo a mitad de año (-0,26 raíz(mm/día) por
#   unidad el 28-jun): El Niño profundiza la sequía de junio-julio. La altitud
#   pesa poco una vez que se conoce el píxel.
#   Prueba (2001-2003, 372 curvas de 136 estaciones, IC 95 % por bootstrap de
#   estaciones): la regresión baja el RMSE L2 de 2,44 a 2,25 mm/día, sube el KGE
#   de 0,38 a 0,45 y acerca el pico de la primera temporada de 23 a 18 días. Mejora
#   en 59 % de las curvas, pero los intervalos se solapan: la ganancia es
#   moderada. Sin el satélite (solo altitud y ONI) el error es mayor que el del
#   satélite crudo (2,61), así que el satélite sí aporta. En 2002 (Niño) el
#   modelo corrige de más (-392 mm/año): el efecto del ONI aprendido en
#   entrenamiento no se cumplió ese año. La reconstrucción en mm es no negativa
#   y solo 0,05 % de los valores se recorta en 0.
#   Esta regresión es la referencia lineal. La red funcional (codificador con W,
#   decodificador FPCA con estas psi_m y cabezas lineal, KAN y MLP) sigue en el
#   archivo de Python con los insumos de para_python/.

saveRDS(list(M = M, lambda = LAMBDA, B = B, mu_c = mu_c, psi_c = psi_c, var_fpca = var_fpca, BETA = BETA, metricas = metricas,
             metricas_curva = metricas_curva, dif_equivalencia = dif_equiv, alt_m = alt_m, alt_s = alt_s),
        file.path(SALIDA_F, "regresion.rds"), compress = "xz")

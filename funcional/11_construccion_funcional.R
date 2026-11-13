# =============================================================================
# 11. CONSTRUCCIÓN DE DATOS FUNCIONALES (ESTACIONES Y SATÉLITE)
# =============================================================================
# Base B-spline cúbica (Cox-de Boor) y mínimos cuadrados programados sin
# librerías; K por ECM, VC y VCG solo con entrenamiento (1981-1995, sin los años
# defectuosos del satélite); verificación de S; covarianza residual; bandas
# punto a punto; matriz de Gram W. Ambas fuentes en el dominio 1-365.
# Entradas: salidas_R_datos/lluvia_valle_limpia.rds, salidas_R_funcional/satelite.rds
# Salida:   salidas_R_funcional/construccion.rds
# =============================================================================

source("funcional/00_funciones.R")
titulo("11. Construcción de datos funcionales")

lluvia <- readRDS("salidas_R_datos/lluvia_valle_limpia.rds")
sat <- readRDS(file.path(SALIDA_F, "satelite.rds"))

d <- lluvia$diario
ancho <- dcast(d, codigo + anio ~ dia, value.var = "P")
id_est <- ancho[, .(codigo, anio)]
id_est[, `:=`(particion = PARTICION(anio), curva_id = paste(codigo, anio, sep = "_"))]
Y <- as.matrix(ancho[, -(1:2)]); storage.mode(Y) <- "double"
stopifnot(ncol(Y) == 365, identical(colnames(Y), as.character(1:365)))
X <- sat$X; id_sat <- sat$meta; T_SAT <- sat$t_sat

completas_est <- rowSums(is.na(Y)) == 0
print(id_est[, .(curvas = .N, completas = sum(completas_est[.I]), estaciones = uniqueN(codigo)), by = particion])
cat("Estaciones:", nrow(Y), "curvas x 365 días, NA:", sum(is.na(Y)), "| satélite:", nrow(X), "curvas x 72 péntadas\n")

grafico_base("11_datos_discretos", function() {
  par(mfrow = c(1, 2))
  i <- which(completas_est & id_est$codigo == "2622330101")[1]
  plot(T_DIA, Y[i, ], type = "h", col = "grey40", xlab = "Día del año", ylab = "mm/día",
       main = paste("Estación", id_est$codigo[i], id_est$anio[i], "(diario)"))
  j <- which(id_sat$celda == sat$emparejamiento[codigo == "2622330101", celda] & id_sat$anio == id_est$anio[i])
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
ent_sat <- id_sat$particion == "entrenamiento" & !id_sat$defectuoso
sel_est <- criterios_K(Y[ent_est, ], T_DIA, c(4:40, seq(45, 100, 5)))
sel_sat <- criterios_K(X[ent_sat, ], T_SAT, 4:40)
sel <- rbind(sel_est[, fuente := "estaciones"], sel_sat[, fuente := "satélite"])
tabla_f(sel, "11_seleccion_K")
resumen_K <- sel[, .(K_VC = K[which.min(VC)], K_VCG = K[which.min(VCG)],
                     K_VC_1pct = paste(range(K[VC <= 1.01 * min(VC)]), collapse = "-")), by = fuente]
resumen_K[, curvas_entrenamiento := c(sum(ent_est), sum(ent_sat))]
print(tabla_f(resumen_K, "11_K_elegido"))
K_EST <- resumen_K[fuente == "estaciones", K_VCG]; K_SAT <- resumen_K[fuente == "satélite", K_VCG]

sel_largo <- melt(sel, id.vars = c("fuente", "K"), variable.name = "criterio")
marcas_K <- sel_largo[criterio != "ECM", .SD[which.min(value)], by = .(fuente, criterio)]
grafico_f(ggplot(sel_largo[K <= 60], aes(K, value)) + geom_line(colour = "steelblue") + geom_point(size = .8, colour = "steelblue") +
            geom_vline(data = marcas_K, aes(xintercept = K), colour = "red", linetype = 2) +
            geom_text(data = marcas_K, aes(K, Inf, label = paste("K =", K)), vjust = 1.5, hjust = -.1, colour = "red", size = 3) +
            facet_wrap(fuente ~ criterio, scales = "free_y", ncol = 3) + scale_y_log10() +
            labs(x = "Número de bases K", y = "Criterio (log)", title = "Selección de K con curvas de entrenamiento"),
          "11_seleccion_K", ancho = 11, alto = 6)

e <- construir(Y, T_DIA, K_EST)
s <- construir(X, T_SAT, K_SAT)
controles <- rbind(data.table(fuente = "estaciones", K = K_EST, t(e$controles)), data.table(fuente = "satélite", K = K_SAT, t(s$controles)))
print(tabla_f(controles, "11_controles_S"))
stopifnot(all(controles[, .(simetria, idempotencia, traza_menos_K)] < 1e-8))

Phi_e <- base_bspline(T_DIA, K_EST); Phi_s365 <- base_bspline(T_DIA, K_SAT)
F_est <- e$coef %*% t(Phi_e); F_sat <- s$coef %*% t(Phi_s365)
negativos <- data.table(fuente = c("estaciones", "satélite"),
                        pct_puntos_negativos = round(100 * c(mean(F_est < 0), mean(F_sat < 0)), 2),
                        pct_curvas_con_negativos = round(100 * c(mean(apply(F_est < 0, 1, any)), mean(apply(F_sat < 0, 1, any))), 1),
                        minimo_mm_dia = round(c(min(F_est), min(F_sat)), 2))
print(tabla_f(negativos, "11_valores_negativos"))

grafico_base("11_bases", function() {
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
tabla_f(data.table(dia = T_DIA, varianza_error = v_t), "11_varianza_error_estaciones")

grafico_base("11_covarianza_residual", function() {
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
js <- which(id_sat$celda == sat$emparejamiento[codigo == id_est$codigo[i_ej], celda] & id_sat$anio == id_est$anio[i_ej])
ee_s <- ee_banda(Phi_s365, s$mapa, s$sigma2[js] * diag(72))
cat(sprintf("Curva %s con hueco de %d días: el EE en los días faltantes es %.2f veces el de la misma curva sin hueco\n",
            id_est$curva_id[ig], sum(!ag$obs), mean((ee_g / ee_g_completa)[!ag$obs])))
cat(sprintf("EE medio de la banda: estación %.2f mm/día (sigma %.1f) | su píxel %.2f mm/día (sigma %.1f)\n",
            mean(ee_iid), sqrt(e$sigma2[i_ej]), mean(ee_s), sqrt(s$sigma2[js])))
tabla_f(data.table(curva_id = id_est$curva_id[i_ej], dia = T_DIA, curva = F_est[i_ej, ], ee_iid, ee_het,
                   inferior = F_est[i_ej, ] - 2 * ee_iid, superior = F_est[i_ej, ] + 2 * ee_iid), "11_banda_ejemplo_estacion")

grafico_base("11_bandas_punto_a_punto", function() {
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

grafico_base("11_mismo_dominio", function() {
  par(mfrow = c(1, 1))
  plot(T_DIA, F_est[i_ej, ], type = "l", lwd = 2, col = 4, ylim = range(c(F_est[i_ej, ], F_sat[js, ], 0)),
       xlab = "Día del año (1-365)", ylab = "mm/día", main = paste("Estación", id_est$curva_id[i_ej], "y su píxel en el mismo dominio"))
  lines(T_DIA, F_sat[js, ], lwd = 2, col = "darkorange")
  points(T_SAT, X[js, ], pch = 16, cex = .6, col = "darkorange")
  legend("topright", c(paste("Estación, K =", K_EST), paste("Satélite, K =", K_SAT)), col = c(4, "darkorange"), lwd = 2, bty = "n")
}, ancho = 9, alto = 4.5)

# QUÉ DICEN LOS DATOS
#   3 326 curvas estación-año (176 estaciones; 2 278 de entrenamiento) y 15 663
#   curvas píxel-año. La base propia coincide con splines::splineDesign (3e-16) y
#   el atajo de la VC con diag(S) da lo mismo que dejar cada día fuera.
#   K se elige solo con entrenamiento (1 899 curvas completas de estación; 8 853
#   de satélite sin 1985 y 1988): VC y VCG coinciden, K = 8 en estaciones y K = 11
#   en satélite. En estaciones la VC es plana entre K = 5 y 17 (a menos de 1 % del
#   mínimo): la lluvia diaria es tan ruidosa que más bases solo persiguen
#   aguaceros; en el satélite el mínimo es nítido, porque con 72 péntadas un K
#   grande deja pocos grados de libertad y la VC se dispara.
#   S es simétrica, idempotente y tr(S) = K hasta 1e-14 en el peor de los ajustes
#   (también en las curvas con días faltantes, que usan su propia S).
#   El error de las estaciones no es homogéneo: su varianza va de 26 (mm/día)^2 en
#   agosto a 264 en octubre, la temporada húmeda; la correlación residual entre
#   días es casi nula (0,04; 0,08 a un día), así que la banda con sigma^2 I es
#   razonable. En el satélite la correlación residual es mayor (0,14): la péntada
#   ya promedia cinco días. Las bandas son punto a punto: cada día tiene 95 % de
#   cobertura por separado, no la curva completa. Un hueco de 15 días ensancha la
#   banda solo 15 % en esos días.
#   Valores negativos: 47 % de las curvas de estación baja de 0 en algún día, casi
#   siempre poco (mediana del mínimo -0,07 mm/día; 8 curvas bajo -5), en la
#   temporada seca de julio-agosto y en los bordes 1-ene y 31-dic (efecto de borde
#   de la base no periódica). En el satélite casi no ocurre (1 %). Para el
#   modelado la curva se recorta en 0 antes de la raíz (D-23).
#   W integra exactamente: 1'W1 = 364 y c'Wc coincide con el trapecio fino.

saveRDS(list(id_est = id_est, id_sat = id_sat, coef_est = e$coef, coef_sat = s$coef, K_est = K_EST, K_sat = K_SAT,
             sigma2_est = e$sigma2, sigma2_sat = s$sigma2, n_obs_est = e$n_obs, W_est = W_est, W_sat = W_sat,
             t_sat = T_SAT, seleccion_K = sel, controles_S = controles, Sigma_res_est = Sigma_res_est, Sigma_res_sat = Sigma_res_sat,
             emparejamiento = sat$emparejamiento, estaciones = lluvia$estaciones),
        file.path(SALIDA_F, "construccion.rds"), compress = "xz")

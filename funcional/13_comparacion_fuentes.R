# =============================================================================
# 13. COMPARACIÓN FUNCIONAL POR FUENTE (ESTACIÓN VS. SU PÍXEL)
# =============================================================================
# Pares estación-año / píxel-año con el emparejamiento de la sección 10.
# Curva de diferencia d(t) = estación - satélite con banda d_barra(t) ± 2EE(t);
# EE por estación (las curvas de una misma estación no son independientes).
# Correlación cruzada funcional y su diagonal, covarianzas de cada fuente,
# comparación por vertiente, altitud y fase ENSO. Los años defectuosos del
# satélite (1985, 1988, 2000) se excluyen del análisis principal y se reportan
# aparte.
# Entrada: salidas_R_funcional/construccion.rds
# Salida:  salidas_R_funcional/comparacion.rds
# =============================================================================

source("funcional/00_funciones.R")
titulo("13. Comparación funcional por fuente")

o <- readRDS(file.path(SALIDA_F, "construccion.rds"))
oni <- fread("datos_web/oni_anual_fase.csv")
Phi_e <- base_bspline(T_DIA, o$K_est); Phi_s <- base_bspline(T_DIA, o$K_sat)

pares <- o$emparejamiento[, .(codigo, nombre, vertiente, altitud, celda, dist_km, estaciones_en_el_pixel)][o$id_est[, .(fila_est = .I, codigo, anio, particion, curva_id)], on = "codigo"]
pares <- o$id_sat[, .(fila_sat = .I, anio, celda, defectuoso)][pares, on = .(anio, celda)]
pares <- oni[, .(anio, ONI_anual, fase)][pares, on = "anio"]
pares[, banda_altitud := cut(altitud, c(-Inf, 1200, 2000, Inf), labels = c("< 1200 m", "1200-2000 m", "> 2000 m"))]
stopifnot(!anyNA(pares$fila_sat), nrow(pares) == nrow(o$id_est))
cat("Pares estación-año / píxel-año:", nrow(pares), "en", uniqueN(pares$codigo), "estaciones y", uniqueN(pares$celda), "píxeles |",
    "en años defectuosos del satélite:", pares[defectuoso == TRUE, .N], "\n")
print(pares[, .(pares = .N, estaciones = uniqueN(codigo)), by = .(particion, defectuoso)][order(particion, defectuoso)])
tabla_f(pares[, .(codigo, nombre, anio, particion, celda, dist_km, estaciones_en_el_pixel, vertiente, altitud, fase, defectuoso)], "13_pares_estacion_pixel")

Y <- o$coef_est[pares$fila_est, ] %*% t(Phi_e)
X <- o$coef_sat[pares$fila_sat, ] %*% t(Phi_s)
D <- Y - X
ok <- !pares$defectuoso

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
print(tabla_f(int_d, "13_intervalos_diferencia"))
cat(sprintf("d(t) medio: %.2f mm/día (%.0f mm/año); días en que la banda excluye el 0: %d de 365\n",
            mean(b$d), mean(b$d) * 365, sum(int_d$dias)))
b_def <- banda_diferencia(D[!ok, ], pares$codigo[!ok])
cat(sprintf("En los años defectuosos: d(t) medio %.2f mm/día (%.0f mm/año)\n", mean(b_def$d), mean(b_def$d) * 365))
tabla_f(data.table(dia = T_DIA, d = b$d, ee = b$ee, inferior = b$d - 2 * b$ee, superior = b$d + 2 * b$ee, d_anios_defectuosos = b_def$d), "13_diferencia_por_dia")

grafico_base("13_diferencia", function() {
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
  legend("bottomleft", c("Años válidos (banda por estación)", "1985, 1988 y 2000", "Días en que la banda excluye el 0"),
         col = c(1, "firebrick", "firebrick"), lty = c(1, 2, 1), lwd = c(2, 1, 1), bty = "n", cex = .8)
}, ancho = 13, alto = 4.8)

titulo("13.2 Correlación cruzada y covarianzas")
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
tabla_f(data.table(dia = T_DIA, r_diagonal = diag_r, r_diagonal_anomalias = R_anom, var_estacion = diag(Cov_YY), var_satelite = diag(Cov_XX)), "13_diagonal_y_varianzas")

cov_comp <- data.table(fuente = c("estación", "satélite"), varianza_total = c(sum(diag(Cov_YY)), sum(diag(Cov_XX))),
                       DE_media = c(mean(sqrt(diag(Cov_YY))), mean(sqrt(diag(Cov_XX)))),
                       cor_abr_oct = c(cov2cor(Cov_YY)[105, 290], cov2cor(Cov_XX)[105, 290]),
                       cor_abr_jul = c(cov2cor(Cov_YY)[105, 200], cov2cor(Cov_XX)[105, 200]),
                       var_explicada_1a_componente = c(eigen(Cov_YY, symmetric = TRUE, only.values = TRUE)$values[1] / sum(diag(Cov_YY)),
                                                       eigen(Cov_XX, symmetric = TRUE, only.values = TRUE)$values[1] / sum(diag(Cov_XX))))
print(tabla_f(cov_comp, "13_covarianzas_comparadas"))
cat(sprintf("Cociente de varianzas satélite / estación por día: mediana %.2f (rango %.2f a %.2f)\n",
            median(diag(Cov_XX) / diag(Cov_YY)), min(diag(Cov_XX) / diag(Cov_YY)), max(diag(Cov_XX) / diag(Cov_YY))))

grafico_base("13_correlacion_cruzada", function() {
  par(mfrow = c(1, 1)); mapa_calor(T_DIA, R_YX, "Correlación cruzada r(s,t): satélite en s, estación en t", c(-1, 1), xlab = "Día t (estación)", ylab = "Día s (satélite)")
}, ancho = 7, alto = 5.5)
grafico_base("13_diagonal_correlacion", function() {
  plot(T_DIA, diag_r, type = "l", lwd = 2, ylim = c(0, 1), xlab = "Día del año", ylab = "r(t, t)", main = "Diagonal de la correlación cruzada")
  lines(T_DIA, R_anom, lwd = 2, col = "firebrick")
  legend("bottomleft", c("Curvas", "Anomalías (sin la climatología de cada estación)"), col = c(1, "firebrick"), lwd = 2, bty = "n")
}, ancho = 9, alto = 4.5)
lim_cov <- range(c(Cov_YY, Cov_XX))
grafico_base("13_covarianzas", function() {
  layout(matrix(1:3, 1)); par(mar = c(4, 4, 3, 1))
  image(T_DIA, T_DIA, Cov_YY, zlim = lim_cov, col = hcl.colors(64, "YlOrRd", rev = TRUE), main = "Covarianza: estaciones", xlab = "Día s", ylab = "Día t")
  image(T_DIA, T_DIA, Cov_XX, zlim = lim_cov, col = hcl.colors(64, "YlOrRd", rev = TRUE), main = "Covarianza: su píxel", xlab = "Día s", ylab = "Día t")
  lim_d <- max(abs(Cov_XX - Cov_YY)) * c(-1, 1)
  image(T_DIA, T_DIA, Cov_XX - Cov_YY, zlim = lim_d, col = hcl.colors(64, "Blue-Red 3"), main = "Satélite - estación", xlab = "Día s", ylab = "Día t")
  layout(1)
}, ancho = 15, alto = 4.8)

titulo("13.3 Por vertiente, altitud y fase ENSO")
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
print(tabla_f(grupos, "13_comparacion_por_grupo"))

pares[, `:=`(P_est = rowMeans(Y) * 365, P_sat = rowMeans(X) * 365)]
por_est <- pares[ok, .(P_est = mean(P_est), P_sat = mean(P_sat), altitud = altitud[1], vertiente = vertiente[1], nombre = nombre[1]), by = codigo]
por_est[, cociente := P_sat / P_est]
ajuste_alt <- lm(log(cociente) ~ altitud, data = por_est[vertiente != "Pacífico"])
cat(sprintf("Estaciones del Cauca: log(satélite/estación) cambia %.3f por cada 1 000 m de altitud (p = %.3g); r = %.2f\n",
            1000 * coef(ajuste_alt)[2], summary(ajuste_alt)$coefficients[2, 4], cor(por_est[vertiente != "Pacífico" & !is.na(altitud), .(altitud, log(cociente))])[1, 2]))
cat(sprintf("Lluvia anual estación vs su píxel (promedio por estación): r = %.2f; en el Cauca r = %.2f\n",
            cor(por_est$P_est, por_est$P_sat), por_est[vertiente != "Pacífico", cor(P_est, P_sat)]))
tabla_f(por_est[order(cociente)], "13_sesgo_por_estacion")
print(rbind(head(por_est[order(cociente)], 4), tail(por_est[order(cociente)], 4)))

grafico_f(ggplot(por_est, aes(P_est, P_sat, colour = vertiente)) + geom_abline(linetype = 2) + geom_point() +
            scale_x_log10() + scale_y_log10() + labs(x = "Estación (mm/año, log)", y = "Su píxel (mm/año, log)", colour = NULL,
                                                   title = "Lluvia anual media: estación frente a su píxel (1981-2003, sin años defectuosos)"),
          "13_anual_estacion_vs_pixel", ancho = 8, alto = 6)
grafico_f(ggplot(por_est[!is.na(altitud)], aes(altitud, cociente, colour = vertiente)) + geom_hline(yintercept = 1, linetype = 2) +
            geom_point() + geom_smooth(method = "lm", formula = y ~ x, se = FALSE) + scale_y_log10() +
            labs(x = "Altitud de la estación (m)", y = "Satélite / estación (log)", colour = NULL, title = "Sesgo relativo del satélite según la altitud"),
          "13_sesgo_altitud", ancho = 8, alto = 5)

curvas_d <- rbindlist(lapply(c("vertiente", "banda_altitud", "fase"), function(var) {
  g <- as.character(pares[[var]]); k <- ok & !is.na(g)
  rbindlist(lapply(sort(unique(g[k])), function(v) { i <- k & g == v; bb <- banda_diferencia(D[i, , drop = FALSE], pares$codigo[i])
    data.table(agrupacion = var, grupo = v, dia = T_DIA, d = bb$d, inf = bb$d - 2 * bb$ee, sup = bb$d + 2 * bb$ee) }))
}))
grafico_f(ggplot(curvas_d, aes(dia, d, colour = grupo, fill = grupo)) + geom_hline(yintercept = 0, linetype = 3) +
            geom_ribbon(aes(ymin = inf, ymax = sup), alpha = .15, colour = NA) + geom_line(linewidth = .8) +
            facet_wrap(~agrupacion, ncol = 3, labeller = as_labeller(c(vertiente = "Vertiente", banda_altitud = "Altitud", fase = "Fase ENSO"))) + labs(x = "Día del año", y = "d(t), mm/día", colour = NULL, fill = NULL,
                                                      title = "Diferencia estación - satélite por grupo (±2 EE por estación)"),
          "13_diferencia_por_grupo", ancho = 14, alto = 4.8)

# QUÉ DICEN LOS DATOS
#   Emparejamiento (sección 10): 3 326 pares estación-año / píxel-año, 176
#   estaciones en 148 píxeles; 446 pares caen en 1985, 1988 o 2000 y se analizan
#   aparte. El EE se calcula entre estaciones (promedio de cada estación primero):
#   es 3,4 veces el que saldría tratando cada curva como independiente, porque los
#   años de una misma estación repiten su sesgo.
#   d(t): en promedio el satélite da 83 mm/año más que la estación, pero el signo
#   cambia con la temporada. La estación supera al satélite del 10-feb al 24-mar
#   (inicio de la primera temporada: el satélite la "arranca" tarde); el satélite
#   supera a la estación del 12-abr al 13-jun, del 3 al 20-ago y del 22-oct al
#   22-nov, es decir, exagera los picos de las dos temporadas y la sequía de
#   agosto. Es el mismo retraso de 2-3 semanas que mostraron las derivadas. La
#   banda excluye el 0 en 176 de 365 días. Los extremos del 1-ene y del 31-dic son
#   efecto de borde de las bases. En 1985, 1988 y 2000 la diferencia se invierte
#   (+336 mm/año): el satélite pierde lluvia en todo el año, lo que confirma
#   marcarlos como defectuosos.
#   Correlación cruzada: la diagonal r(t,t) vale 0,63 en promedio (0,80 a mediados
#   de julio, 0,33 en el borde de enero), y la correlación máxima está en el mismo
#   día (desfase mediano de -1 día). Casi toda esa correlación viene de que el
#   satélite reconoce qué sitios son húmedos y cuáles secos: con anomalías
#   (restando la climatología de cada estación) la diagonal baja a 0,26, con un
#   máximo de 0,54 en febrero y casi 0 en noviembre. El satélite sigue mal la
#   variación de un año a otro en la segunda temporada.
#   Covarianzas: la forma es casi igual (correlación abril-octubre 0,73 en ambas;
#   la primera componente explica 75 % y 77 %), pero la varianza del satélite es
#   solo una cuarta parte de la de las estaciones (mediana 0,25 por día): el píxel
#   de ~5,5 km y la péntada promedian los aguaceros locales.
#   Subregión: en el valle del Cauca el satélite sobreestima 10 % (+157 mm/año) y la
#   banda excluye el 0 en 316 días. En la vertiente Pacífica el sesgo depende del
#   sitio. El satélite subestima a la mitad las estaciones más lluviosas: Bahía
#   Málaga y San Francisco en la costa de Buenaventura, y El Tabano en el bajo
#   Calima (cociente 0,46 a 0,51). En cambio, duplica o triplica la lluvia de los
#   sitios secos a sotavento: Loboguerrero en el cañón del Dagua (2,8), Puerto
#   Nuevo y Moravito en el norte (2,3 y 2,0) y La Despensa en La Unión (2,4). Al
#   otro extremo está La Patagonia, una ladera de Buga a 1 950 m (0,41). El
#   satélite suaviza hacia la media regional: sube lo seco y baja lo muy húmedo.
#   Entre estaciones, la lluvia
#   anual del píxel explica la de la estación con r = 0,79, pero dentro del valle
#   del Cauca solo con r = 0,42.
#   Altitud: en el Cauca el cociente satélite/estación baja 14 % por cada 1 000 m
#   (p = 0,002): el satélite exagera en el fondo del valle y se queda corto en las
#   laderas altas, donde la lluvia orográfica no se ve a 5 km de resolución.
#   ENSO: en años Niña las estaciones reciben 2 003 mm y el satélite 1 895
#   (subestima); en años Niño 1 580 frente a 1 829 (sobreestima 16 %). El satélite
#   amortigua la señal ENSO aun sin los años defectuosos. Esto justifica dar al
#   modelo la altitud y el ONI junto con la curva satelital, como pide el
#   enunciado.

saveRDS(list(pares = pares, d = b, d_defectuosos = b_def, intervalos = int_d, R_YX = R_YX, Cov_YY = Cov_YY, Cov_XX = Cov_XX,
             grupos = grupos, sesgo_estacion = por_est),
        file.path(SALIDA_F, "comparacion.rds"), compress = "xz")

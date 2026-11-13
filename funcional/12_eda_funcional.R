# =============================================================================
# 12. ANÁLISIS EXPLORATORIO DE DATOS FUNCIONALES (CADA FUENTE)
# =============================================================================
# Media paso a paso (verificada con fda), varianza, desviación estándar,
# covarianza y correlación, derivadas y temporadas, profundidad FM y MBD,
# boxplot funcional con atípicas por estación y año, y contexto climatológico
# (vertiente, altitud y ENSO). Se usan todas las estaciones incluidas y todos
# los años 1981-2003; el satélite se analiza en los píxeles del Valle.
# Entrada: salidas_R_funcional/construccion.rds
# Salida:  salidas_R_funcional/eda_funcional.rds
# =============================================================================

source("funcional/00_funciones.R")
titulo("12. EDA funcional")

o <- readRDS(file.path(SALIDA_F, "construccion.rds"))
oni <- fread("datos_web/oni_anual_fase.csv")
est <- o$estaciones

Phi_e <- base_bspline(T_DIA, o$K_est); Phi_s <- base_bspline(T_DIA, o$K_sat)
valle_s <- o$id_sat$en_valle
C_e <- o$coef_est; C_s <- o$coef_sat[valle_s, ]
id_e <- est[, .(codigo, nombre, vertiente, altitud)][o$id_est, on = "codigo"]
id_s <- o$id_sat[valle_s]
id_s[, vertiente := fifelse(lon < -76.8, "Pacífico", "Cauca (interandina)")]
Ye <- C_e %*% t(Phi_e); Xs <- C_s %*% t(Phi_s)
cat("Estaciones:", nrow(Ye), "curvas,", uniqueN(id_e$codigo), "estaciones, años", paste(range(id_e$anio), collapse = "-"),
    "| satélite (Valle):", nrow(Xs), "curvas,", uniqueN(id_s$celda), "píxeles\n")
print(id_e[, .N, by = vertiente])

titulo("12.1 Media funcional")
media_e <- drop(Phi_e %*% colMeans(C_e)); media_s <- drop(Phi_s %*% colMeans(C_s))
fd_e <- fda::fd(t(C_e), fda::create.bspline.basis(RANGO, nbasis = o$K_est, norder = ORDEN))
fd_s <- fda::fd(t(C_s), fda::create.bspline.basis(RANGO, nbasis = o$K_sat, norder = ORDEN))
verif_media <- data.table(fuente = c("estaciones", "satélite"),
                          vs_curvas_evaluadas = c(max(abs(media_e - colMeans(Ye))), max(abs(media_s - colMeans(Xs)))),
                          vs_fda_mean.fd = c(max(abs(media_e - fda::eval.fd(T_DIA, fda::mean.fd(fd_e)))),
                                             max(abs(media_s - fda::eval.fd(T_DIA, fda::mean.fd(fd_s))))))
print(tabla_f(verif_media, "12_verificacion_media"))
stopifnot(verif_media$vs_fda_mean.fd < 1e-10)
cat(sprintf("Lluvia media: estaciones %.0f mm/año | satélite en el Valle %.0f mm/año\n", mean(media_e) * 365, mean(media_s) * 365))

titulo("12.2 Varianza, desviación estándar, covarianza y correlación")
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
print(tabla_f(resumen_disp, "12_dispersion_resumen"))

grafico_base("12_media_y_dispersion", function() {
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
  grafico_base(paste0("12_covarianza_", f[[2]]), function() mapa_calor(T_DIA, f[[1]]$cov, paste("Covarianza funcional:", f[[2]]), paleta = "YlOrRd"), ancho = 7, alto = 5.5)
  grafico_base(paste0("12_correlacion_", f[[2]]), function() mapa_calor(T_DIA, f[[1]]$cor, paste("Correlación funcional:", f[[2]]), c(-1, 1)), ancho = 7, alto = 5.5)
}

titulo("12.3 Derivadas y temporadas de lluvia")
DPhi_e <- deriv_base(T_DIA, o$K_est); DPhi_s <- deriv_base(T_DIA, o$K_sat)
ref_d <- splines::splineDesign(nodos(o$K_est), T_DIA, ord = ORDEN, derivs = rep(1, 365), outer.ok = TRUE)
cat("Derivada de la base propia vs splines::splineDesign:", format(max(abs(DPhi_e - ref_d)), digits = 2), "\n")
dmu_e <- drop(DPhi_e %*% colMeans(C_e)); dmu_s <- drop(DPhi_s %*% colMeans(C_s))
cruces <- function(media, dmu, fuente) {
  k <- which(diff(sign(dmu)) != 0)
  data.table(fuente, dia = k + 1, fecha = fecha_dia(k + 1), tipo = fifelse(dmu[k] > 0, "máximo", "mínimo"), mm_dia = round(media[k + 1], 2))
}
extremos_media <- rbind(cruces(media_e, dmu_e, "estaciones"), cruces(media_s, dmu_s, "satélite"))
print(tabla_f(extremos_media, "12_extremos_de_la_media"))

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
print(tabla_f(dcast(resumen_tem, momento ~ fuente + vertiente, value.var = "fecha_mediana"), "12_temporadas_fechas"))
tabla_f(resumen_tem, "12_temporadas_detalle")

grafico_base("12_derivada_media", function() {
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
          "12_temporadas_por_curva", ancho = 11, alto = 5)

titulo("12.4 Profundidad, mediana y media recortada")
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
grafico_base("12_centralidad", function() {
  par(mfrow = c(1, 2))
  for (f in list(list(media_e, ce_e, "Estaciones"), list(media_s, ce_s, "Satélite (Valle)"))) {
    matplot(T_DIA, cbind(f[[1]], f[[2]]$mediana_FM, f[[2]]$mediana_MBD, f[[2]]$recortada_25), type = "l", lty = c(1, 2, 3, 1),
            col = c(4, 2, 1, 3), lwd = 2, xlab = "Día del año", ylab = "mm/día", main = paste("Centralidad:", f[[3]]))
    legend("topright", c("Media", "Mediana FM", "Mediana MBD", "Media recortada 25 % (FM)"), col = c(4, 2, 1, 3), lty = c(1, 2, 3, 1), lwd = 2, bty = "n", cex = .75)
  }
}, ancho = 12, alto = 4.5)

titulo("12.5 Boxplot funcional y curvas atípicas")
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
tabla_f(atip_e[order(-P_anual_mm)], "12_atipicas_estaciones_todas")
tabla_f(atip_v, "12_atipicas_estaciones_por_vertiente")
cat(sprintf("Atípicas, todas las estaciones juntas: %d de %d curvas (%d estaciones; %d del Pacífico)\n",
            nrow(atip_e), nrow(Ye), uniqueN(atip_e$codigo), atip_e[vertiente == "Pacífico", .N]))
print(head(atip_e[, .(curvas = .N, anios = paste(sort(anio), collapse = " "), P_anual_mediana = median(P_anual_mm)), by = .(codigo, nombre, vertiente)][order(-curvas)], 12))
cat("Atípicas dentro de cada vertiente:\n"); print(atip_v[, .(curvas = .N, estaciones = uniqueN(codigo)), by = vertiente])
print(atip_v[vertiente != "Pacífico", .N, by = .(codigo, nombre, anio)])
atip_s <- id_s[bx_s$atip]
cat(sprintf("Atípicas del satélite: %d de %d curvas píxel-año (%.1f %%)\n", nrow(atip_s), nrow(Xs), 100 * mean(bx_s$atip)))
print(atip_s[, .N, by = .(vertiente)]); print(head(atip_s[, .N, by = anio][order(-N)], 6))
tabla_f(atip_s[, .(anio, celda, fila, col, lon, lat, vertiente, n_tope)], "12_atipicas_satelite")

anios_e <- sort(unique(id_e$anio))
media_anual <- function(Y, g, a) t(sapply(a, function(x) colMeans(Y[g == x, , drop = FALSE])))
Ya_e <- media_anual(Ye[id_e$vertiente != "Pacífico", ], id_e$anio[id_e$vertiente != "Pacífico"], anios_e)
Ya_s <- media_anual(Xs, id_s$anio, ANIOS)
pa_e <- profundidades(Ya_e)$MBD; pa_s <- profundidades(Ya_s)$MBD
bxa_e <- boxplot_funcional(Ya_e, pa_e); bxa_s <- boxplot_funcional(Ya_s, pa_s)
prof_anual <- oni[data.table(anio = anios_e, MBD_estaciones_Cauca = round(pa_e, 3), atipico_estaciones = bxa_e$atip,
                             MBD_satelite = round(pa_s, 3), atipico_satelite = bxa_s$atip,
                             P_estaciones_Cauca = round(rowMeans(Ya_e) * 365), P_satelite = round(rowMeans(Ya_s) * 365)), on = "anio"]
print(tabla_f(prof_anual[order(MBD_estaciones_Cauca)], "12_profundidad_anual"))
cat(sprintf("Correlación (Spearman) entre la profundidad anual de estaciones y de satélite: %.2f\n", cor(pa_e, pa_s, method = "spearman")))

grafico_base("12_boxplot_funcional", function() {
  par(mfrow = c(2, 3))
  dibujar_boxplot(Ye, bx_e, paste("Estaciones, todas:", sum(bx_e$atip), "atípicas"))
  for (v in names(bx_v)) dibujar_boxplot(Ye[grupos_v[[v]], ], bx_v[[v]], paste("Estaciones,", v, ":", sum(bx_v[[v]]$atip), "atípicas"))
  dibujar_boxplot(Xs, bx_s, paste("Satélite (Valle):", sum(bx_s$atip), "atípicas"))
  dibujar_boxplot(Ya_e, bxa_e, "Curvas medias anuales: estaciones del Cauca")
  dibujar_boxplot(Ya_s, bxa_s, "Curvas medias anuales: satélite")
}, ancho = 15, alto = 9)

titulo("12.6 Contexto: vertiente, altitud y ENSO")
id_e[, banda_altitud := cut(altitud, c(-Inf, 1200, 2000, Inf), labels = c("< 1200 m", "1200-2000 m", "> 2000 m"))]
id_e <- oni[, .(anio, fase)][id_e, on = "anio"]; id_s <- oni[, .(anio, fase)][id_s, on = "anio"]
grupo_media <- function(Y, g) { g <- as.character(g); k <- !is.na(g); sapply(split(seq_len(nrow(Y))[k], g[k]), function(i) colMeans(Y[i, , drop = FALSE])) }
m_vert_e <- grupo_media(Ye, id_e$vertiente); m_vert_s <- grupo_media(Xs, id_s$vertiente)
m_alt <- grupo_media(Ye[id_e$vertiente != "Pacífico", ], id_e$banda_altitud[id_e$vertiente != "Pacífico"])
m_enso_e <- grupo_media(Ye[id_e$vertiente != "Pacífico", ], id_e$fase[id_e$vertiente != "Pacífico"])
m_enso_s <- grupo_media(Xs[id_s$vertiente != "Pacífico", ], id_s$fase[id_s$vertiente != "Pacífico"])
k_sd <- id_s$vertiente != "Pacífico" & !id_s$defectuoso
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
cat("Años defectuosos del satélite y su fase ENSO:\n"); print(oni[anio %in% ANIOS_SAT_DEFECTUOSOS])
contexto[, P_anual := round(P_anual)]
print(tabla_f(contexto, "12_contexto_vertiente_altitud_ENSO"))
picos_enso <- rbindlist(lapply(colnames(m_enso_e), function(f) data.table(fase = f, pico_1 = fecha_dia(59 + which.max(m_enso_e[60:200, f])),
  mm_pico_1 = round(max(m_enso_e[60:200, f]), 2), pico_2 = fecha_dia(229 + which.max(m_enso_e[230:365, f])), mm_pico_2 = round(max(m_enso_e[230:365, f]), 2),
  minimo_jul_ago = round(min(m_enso_e[160:260, f]), 2))))
print(tabla_f(picos_enso, "12_picos_por_fase_ENSO"))

grafico_base("12_contexto", function() {
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
#   Representatividad: se usan las 176 estaciones incluidas y sus 3 326 curvas de
#   1981-2003, sin submuestrear (la sección 8 muestra que las incluidas tienen la
#   misma altitud mediana y casi la misma lluvia que las 220 con datos), y los 676
#   píxeles del Valle en los 23 años.
#   Media: Phi c_barra coincide con fda::mean.fd (5e-15). Ambas fuentes muestran
#   el régimen bimodal del valle interandino: estaciones con máximos el 22-abr
#   (6,3 mm/día) y el 28-oct (7,2), mínimos el 25-ene y el 1-ago (2,2); el
#   satélite pone los máximos 2-3 semanas después (13-may y 5-nov) y con casi el
#   doble de magnitud, porque su media pesa por área y la vertiente Pacífica (un
#   tercio de los píxeles) llueve 4 600 mm/año. Por vertiente, el satélite supera
#   a las estaciones en ambas: 1 856 frente a 1 554 mm en el Cauca y 4 645 frente
#   a 3 099 en el Pacífico.
#   Dispersión: la DE de las estaciones es máxima en la segunda temporada (5,2 mm/día
#   el 18-oct, sin contar el borde del 1-ene) y mínima el 27-feb; el CV medio es
#   0,96, es decir, la variación entre curvas es tan grande como la media, sobre
#   todo por las diferencias entre estaciones (la costa frente al valle). La del
#   satélite es mayor en noviembre (7,4) y más concentrada en las temporadas.
#   Covarianza y correlación: la correlación es positiva en todo el año (mínimo
#   0,39 en estaciones y 0,43 en satélite): una estación húmeda lo es todo el año,
#   así que el nivel de cada sitio domina la covarianza. Dentro de esa base, las
#   dos temporadas se parecen entre sí (r = 0,74 entre abril y octubre) más que
#   cada una con la sequía de julio (0,63).
#   Derivadas y temporadas: en el valle del Cauca la primera temporada empieza a
#   finales de febrero, tiene su pico a mediados de abril y termina a mediados de
#   junio; la segunda empieza a mediados de septiembre y tiene su pico el 1-nov.
#   En 62 % de las curvas de estación la lluvia sigue bajando hasta el borde del
#   año y el fin de la segunda temporada no se resuelve antes del 16-dic. El
#   satélite fecha los picos 2-3 semanas tarde y con menos dispersión entre curvas
#   (péntadas más suaves).
#   Profundidad: FM y MBD ordenan casi igual (Spearman 0,99 y 1,00). La curva más
#   profunda es Betania 2002. Boxplot funcional (mismas atípicas que fda::fbplot):
#   con todas las estaciones juntas hay 62 atípicas, 61 de cinco estaciones de la
#   costa Pacífica (Bahía Málaga, Yurumanguí, San Francisco, Triana y El Tabano,
#   7 300 a 12 600 mm/año): no son errores sino otro clima. Por eso se separa por
#   vertiente: en el Cauca solo quedan Tenjo 1982 (4 482 mm, con el máximo en el
#   borde del 1-ene, candidata a revisión de datos) y Las Brisas 1981. Las 330
#   atípicas del satélite son todas del Pacífico y 229 tocan el tope de 200 mm. Ningún
#   año es atípico en las curvas medias anuales; los menos profundos de las
#   estaciones son 1988, 1984 y 1999, tres años Niña.
#   Orografía: en el Cauca la franja de 1 200-2 000 m es la más lluviosa
#   (1 718 mm), por encima de > 2 000 m (1 616) y del fondo del valle (1 332).
#   ENSO: en las estaciones del Cauca La Niña aporta +11 % (1 718 mm) y El Niño
#   -9 % (1 402) frente a neutro (1 545); la Niña se nota en la segunda temporada
#   (7,1 frente a 6,3 mm/día) y en una sequía de mitad de año más húmeda (2,4 frente
#   a 1,4). El satélite NO muestra esa señal (Niña 1 770 < neutro 1 892) porque
#   sus tres años defectuosos (1985, 1988 y 2000) son años Niña; sin ellos
#   recupera el orden esperado (Niña 1 949 > neutro 1 892 > Niño 1 875), aunque
#   más débil que en las estaciones.

saveRDS(list(media_e = media_e, media_s = media_s, me = me, ms = ms, pr_e = pr_e, pr_s = pr_s, atip_e = atip_e, atip_s = atip_s,
             temporadas = resumen_tem, prof_anual = prof_anual, contexto = contexto),
        file.path(SALIDA_F, "eda_funcional.rds"), compress = "xz")

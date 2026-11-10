#####################################################################################
# Proyecto final Analítica de Datos 2026-2S — lluvia del Valle del Cauca (1981–2003)
# Fase 2: Análisis exploratorio de datos funcionales (versión en R)
#
# Hace en R lo mismo que 04_eda_funcional.ipynb, a partir de los coeficientes que construyó
# 03_construccion_funcional.R (no usa los coeficientes de Python):
#   1. media funcional paso a paso (c_barra' phi)
#   2. varianza, desviación estándar, CV, covarianza Phi Sigma_c Phi' y correlación
#   3. derivada de la base programada a mano y fechas de las temporadas (picos, inicio y fin)
#   4. profundidad Fraiman-Muniz y MBD (fórmula de rangos), mediana y media recortada
#   5. boxplot funcional (Sun y Genton): todas las estaciones, estaciones del Cauca, satélite y curvas medias anuales
#   6. contexto: cuenca y fase ENSO (ONI)
#   7. comparación numérica con los resultados de Python
# La librería fda solo se usa, si está instalada, para VERIFICAR la media y el boxplot; nunca para construir.
#
# Ejecutar de arriba abajo desde la carpeta raíz del repositorio, después de 03_construccion_funcional.R.
# Salidas: fda_fase2/figuras_R/ (PDF) y fda_fase2/tablas_R/ (CSV)
#####################################################################################

rm(list = ls())
# setwd("ruta/al/repositorio")          # descomentar si RStudio no abre en la raíz del repositorio
if (!file.exists("fda_fase1/coef_satelite_R.csv.gz")) stop("Ejecute antes 03_construccion_funcional.R")

RANGO <- c(1, 365); ORDEN <- 4L; ANIOS <- 1981:2003; t <- 1:365
DIR_FIG <- "fda_fase2/figuras_R"; DIR_TAB <- "fda_fase2/tablas_R"
dir.create(DIR_FIG, showWarnings = FALSE, recursive = TRUE); dir.create(DIR_TAB, showWarnings = FALSE)

#------------------------------------------------------------------------------------
# Funciones auxiliares
#------------------------------------------------------------------------------------
guardar_csv <- function(x, nombre) write.csv(x, file.path(DIR_TAB, nombre), row.names = FALSE, na = "")
figura <- function(nombre, dibujar, width = 11, height = 4) {
  pdf(file.path(DIR_FIG, nombre), width = width, height = height); dibujar(); dev.off()
  if (interactive()) dibujar()
  invisible(NULL)
}
mostrar_tabla <- function(x, titulo, filas = 6L) {
  cat("\n", titulo, " (", nrow(x), " filas)\n", sep = ""); print(utils::head(x, filas), row.names = FALSE); invisible(x)
}
mapa_calor <- function(z, titulo, limites = range(z)) {
  colores <- hcl.colors(64, "Blue-Red 3")
  layout(matrix(1:2, 1), widths = c(5, 1)); par(mar = c(4, 4, 3, 1))
  image(t, t, z, zlim = limites, col = colores, main = titulo, xlab = "Día s", ylab = "Día t")
  par(mar = c(4, 0, 3, 4)); escala <- seq(limites[1], limites[2], length.out = 65)
  image(c(0, 1), escala, matrix(rep(escala, each = 2), nrow = 2), col = colores, zlim = limites, axes = FALSE, xlab = "", ylab = "")
  axis(4, las = 1, cex.axis = .8); box(); layout(1)
}

# Base B-spline y su derivada (Cox-de Boor programado a mano, igual que en la fase 1)
nodos <- function(K, orden = ORDEN) {
  interiores <- seq(RANGO[1], RANGO[2], length.out = K - orden + 2L)
  c(rep(RANGO[1], orden), interiores[-c(1L, length(interiores))], rep(RANGO[2], orden))
}
base_en_nodos <- function(tiempo, tau, orden) {
  B <- outer(tiempo, seq_len(length(tau) - 1L), function(x, j) (x >= tau[j]) & (x < tau[j + 1L])) * 1
  B[tiempo == max(tau), max(which(diff(tau) > 0))] <- 1
  if (orden > 1L) for (m in 2:orden) {
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
# Derivada: phi_j'(t) = (m-1) [B_{j,m-1}(t)/(tau_{j+m-1}-tau_j) - B_{j+1,m-1}(t)/(tau_{j+m}-tau_{j+1})]
deriv_base <- function(tiempo, K) {
  tau <- nodos(K); B <- base_en_nodos(tiempo, tau, ORDEN - 1L); D <- matrix(0, length(tiempo), K)
  for (j in seq_len(K)) {
    a <- tau[j + ORDEN - 1L] - tau[j]; b <- tau[j + ORDEN] - tau[j + 1L]
    if (a > 0) D[, j] <- D[, j] + (ORDEN - 1L) * B[, j] / a
    if (b > 0) D[, j] <- D[, j] - (ORDEN - 1L) * B[, j + 1L] / b
  }
  D
}

# Profundidad Fraiman-Muniz y MBD (bandas de 2 curvas, inclusivas) por la fórmula de rangos.
# En cada día: L curvas por debajo y U por encima; las bandas que NO contienen a la curva son C(L,2) + C(U,2).
profundidades <- function(Y) {
  N <- nrow(Y); fm <- mbd <- numeric(N)
  for (j in seq_len(ncol(Y))) {
    r <- rank(Y[, j], ties.method = "max")
    L <- rank(Y[, j], ties.method = "min") - 1; U <- N - r
    fm  <- fm + 1 - abs(.5 - r / N)
    mbd <- mbd + 1 - (L * (L - 1) + U * (U - 1)) / (N * (N - 1))
  }
  list(FM = fm / ncol(Y), MBD = mbd / ncol(Y))
}
# Boxplot funcional: región central = envolvente del 50 % más profundo (MBD); atípica si sale de los límites 1,5
boxplot_funcional <- function(Y, factor = 1.5, prof = profundidades(Y)$MBD) {
  centro <- order(-prof)[seq_len(ceiling(nrow(Y) / 2))]
  lo <- apply(Y[centro, , drop = FALSE], 2, min); hi <- apply(Y[centro, , drop = FALSE], 2, max)
  inf <- lo - factor * (hi - lo); sup <- hi + factor * (hi - lo)
  atip <- apply(Y, 1, function(y) any(y < inf | y > sup))
  buenas <- Y[!atip, , drop = FALSE]
  list(prof = prof, lo = lo, hi = hi, inf = inf, sup = sup, atip = atip, mediana = Y[which.max(prof), ],
       bigote_inf = apply(buenas, 2, min), bigote_sup = apply(buenas, 2, max))
}
dibujar_boxplot <- function(Y, bx, titulo) {
  plot(t, bx$mediana, type = "n", ylim = range(c(bx$inf, bx$sup)), xlab = "Día del año", ylab = "mm/día", main = titulo)
  polygon(c(t, rev(t)), c(bx$lo, rev(bx$hi)), col = "#b39ddb", border = NA)
  lines(t, bx$inf, lty = 2); lines(t, bx$sup, lty = 2)
  lines(t, bx$bigote_inf, col = 4); lines(t, bx$bigote_sup, col = 4)
  if (any(bx$atip)) matlines(t, t(Y[bx$atip, , drop = FALSE]), col = "#cc333355", lty = 1)
  lines(t, bx$mediana, lwd = 2)
  legend("topleft", c("Región central 50 %", "Mediana (MBD)", "Límites 1,5", "Bigotes", "Atípicas"),
         col = c("#b39ddb", 1, 1, 4, "#cc3333"), lty = c(NA, 1, 2, 1, 1), pch = c(15, NA, NA, NA, NA), cex = .7, bty = "n")
}

#------------------------------------------------------------------------------------
# Datos funcionales: coeficientes de la fase 1 (R)
#------------------------------------------------------------------------------------
ce <- read.csv("fda_fase1/coef_estaciones_R.csv"); cs <- read.csv(gzfile("fda_fase1/coef_satelite_R.csv.gz"))
C_est <- as.matrix(ce[, -(1:2)]); C_sat <- as.matrix(cs[, -(1:3)])
K_EST <- ncol(C_est); K_SAT <- ncol(C_sat)
spec <- paste(readLines("fda_fase1/base_especificacion.json", warn = FALSE, encoding = "UTF-8"), collapse = "")
stopifnot(K_EST == as.integer(sub('.*"K_estaciones": *([0-9]+).*', "\\1", spec, useBytes = TRUE)),
          K_SAT == as.integer(sub('.*"K_satelite": *([0-9]+).*', "\\1", spec, useBytes = TRUE)))
# Cuenca de cada estación: los dos primeros dígitos del código IDEAM son la zona hidrográfica (26 = Cauca)
cuenca <- ifelse(substr(ce$STATION_ID, 1, 2) == "26", "Cauca", "Pacífico")
# Píxeles: al oeste de -76,8° vertiente Pacífica; al este, valle del Cauca y cordilleras
sat_lon <- -77.55 + .05 * (cs$col + .5)
region_sat <- ifelse(sat_lon < -76.8, "Pacífico", "Cauca")
Phi_e <- base_bspline(t, K_EST); Phi_s <- base_bspline(t, K_SAT)
Ye <- C_est %*% t(Phi_e); Xs <- C_sat %*% t(Phi_s)            # curvas evaluadas en los 365 días
cat("K estaciones =", K_EST, "| K satélite =", K_SAT, "| curvas: estaciones", nrow(Ye), "(Cauca", sum(cuenca == "Cauca"),
    ") | satélite en tierra", nrow(Xs), "\n")

#------------------------------------------------------------------------------------
# 1. Media funcional paso a paso
#------------------------------------------------------------------------------------
media_e <- as.vector(Phi_e %*% colMeans(C_est)); media_s <- as.vector(Phi_s %*% colMeans(C_sat))
cat("Media paso a paso vs media de las curvas evaluadas:", format(max(abs(media_e - colMeans(Ye))), digits = 2), "\n")
estado_fda <- "fda no instalado (la verificación con librería se hace en Python con scikit-fda)"
if (requireNamespace("fda", quietly = TRUE)) {
  fd_e <- fda::fd(t(C_est), fda::create.bspline.basis(RANGO, nbasis = K_EST, norder = ORDEN))
  err <- max(abs(media_e - fda::eval.fd(t, fda::mean.fd(fd_e)))); stopifnot(err < 1e-9)
  estado_fda <- paste("media igual a fda::mean.fd, error", format(err, digits = 2))
}
cat("Verificación con librería:", estado_fda, "\n")

#------------------------------------------------------------------------------------
# 2. Varianza, DE, CV, covarianza y correlación
#------------------------------------------------------------------------------------
momentos <- function(C, Phi) {
  Y <- C %*% t(Phi); mu <- as.vector(Phi %*% colMeans(C))
  varianza <- colSums(sweep(Y, 2, mu)^2) / (nrow(Y) - 1)
  Cc <- sweep(C, 2, colMeans(C)); covar <- Phi %*% (crossprod(Cc) / (nrow(C) - 1)) %*% t(Phi)   # Phi Sigma_c Phi'
  stopifnot(max(abs(diag(covar) - varianza)) < 1e-8)
  list(mu = mu, var = varianza, de = sqrt(varianza), cv = sqrt(varianza) / mu, cov = covar, cor = cov2cor(covar))
}
me <- momentos(C_est, Phi_e); ms <- momentos(C_sat, Phi_s)
cat("Cov paso a paso vs cov() de R (estaciones):", format(max(abs(me$cov - cov(Ye))), digits = 2), "\n")
cat(sprintf("Correlación 15-abr / 17-oct (días 105 y 290): estaciones %.2f, satélite %.2f\n", me$cor[105, 290], ms$cor[105, 290]))
cat(sprintf("Día de mayor DE: estaciones %d, satélite %d | CV el 1-ago: estaciones %.2f, satélite %.2f\n",
            which.max(me$de), which.max(ms$de), me$cv[213], ms$cv[213]))

figura("media_de_covarianza.pdf", function() {
  par(mfrow = c(2, 3)); set.seed(1)
  for (d in list(list(Ye, me, "Estaciones", 4), list(Xs, ms, "Satélite (tierra)", "darkorange"))) {
    matplot(t, t(d[[1]][sample(nrow(d[[1]]), 200), ]), type = "l", lty = 1, col = "grey80", xlab = "Día del año",
            ylab = "mm/día", main = paste(d[[3]], ": media ± 2 DE"))
    lines(t, d[[2]]$mu, lwd = 3, col = d[[4]])
    lines(t, d[[2]]$mu + 2 * d[[2]]$de, lty = 2, col = 2, lwd = 2); lines(t, d[[2]]$mu - 2 * d[[2]]$de, lty = 2, col = 2, lwd = 2)
    image(t, t, d[[2]]$cov, col = hcl.colors(30), main = paste("Covarianza:", d[[3]]), xlab = "s", ylab = "t")
    image(t, t, d[[2]]$cor, zlim = c(-1, 1), col = hcl.colors(30, "Blue-Red 3"), main = paste("Correlación:", d[[3]]), xlab = "s", ylab = "t")
  }
}, width = 12, height = 7)
figura("dispersion.pdf", function() {
  par(mfrow = c(1, 3))
  for (v in c("var", "de", "cv")) {
    plot(t, me[[v]], type = "l", lwd = 2, col = 4, ylim = range(c(me[[v]], ms[[v]])), xlab = "Día del año",
         ylab = switch(v, var = "(mm/día)^2", de = "mm/día", cv = "DE / media"),
         main = switch(v, var = "Varianza", de = "Desviación estándar", cv = "Coeficiente de variación"))
    lines(t, ms[[v]], lwd = 2, col = "darkorange")
    legend("topright", c("Estaciones", "Satélite"), col = c(4, "darkorange"), lwd = 2, bty = "n")
  }
})
for (d in list(list(me, "estaciones"), list(ms, "satelite"))) {
  figura(paste0("covarianza_", d[[2]], ".pdf"), function() mapa_calor(d[[1]]$cov, paste("Covarianza:", d[[2]])), width = 7, height = 5)
  figura(paste0("correlacion_", d[[2]], ".pdf"), function() mapa_calor(d[[1]]$cor, paste("Correlación:", d[[2]]), c(-1, 1)),
         width = 7, height = 5)
}

#------------------------------------------------------------------------------------
# 3. Derivadas y fechas de las temporadas
#------------------------------------------------------------------------------------
DPhi_e <- deriv_base(t, K_EST); DPhi_s <- deriv_base(t, K_SAT)
ref <- splines::splineDesign(nodos(K_EST), t, ord = ORDEN, derivs = rep(1, length(t)), outer.ok = TRUE)
cat("Derivada de la base propia vs splines::splineDesign:", format(max(abs(DPhi_e - ref)), digits = 2), "\n")
DYe <- C_est %*% t(DPhi_e); DXs <- C_sat %*% t(DPhi_s)
dmu_e <- as.vector(DPhi_e %*% colMeans(C_est)); dmu_s <- as.vector(DPhi_s %*% colMeans(C_sat))

# Pico 1: máximo entre los días 60 y 200; pico 2: máximo desde el día 230. Inicio: máxima pendiente antes del pico;
# fin: mínima pendiente después del pico. Si el pico 2 cae en el día 365, su fin no está definido (efecto de borde).
temporadas <- function(Y, D) {
  t(sapply(seq_len(nrow(Y)), function(i) {
    y <- Y[i, ]; d <- D[i, ]
    p1 <- which(t >= 60 & t <= 200)[which.max(y[t >= 60 & t <= 200])]
    p2 <- which(t >= 230)[which.max(y[t >= 230])]
    c(inicio_1 = which.max(d[1:(p1 - 1)]), pico_1 = p1,
      fin_1 = p1 + which.min(d[(p1 + 1):230]), inicio_2 = 183 + which.max(d[184:(p2 - 1)]),
      pico_2 = p2, fin_2 = if (p2 < 365) p2 + which.min(d[(p2 + 1):365]) else NA)
  }))
}
tem_e <- temporadas(Ye, DYe)
fechas <- function(x) format(as.Date(round(x) - 1, origin = "2001-01-01"), "%d-%b")
mediana_tem <- apply(tem_e, 2, median, na.rm = TRUE)
mostrar_tabla(data.frame(fecha = names(mediana_tem), dia_mediano = mediana_tem, fecha_mediana = fechas(mediana_tem)),
              "Estaciones: fechas medianas de las temporadas")
figura("derivada_media.pdf", function() {
  par(mfrow = c(1, 2))
  for (d in list(list(media_e, dmu_e, "Estaciones", 4), list(media_s, dmu_s, "Satélite", "darkorange"))) {
    plot(t, d[[2]], type = "l", lwd = 2, col = d[[4]], xlab = "Día del año", ylab = "mm/día por día",
         main = paste(d[[3]], ": derivada de la curva media")); abline(h = 0, lty = 2, col = "grey40")
  }
})

#------------------------------------------------------------------------------------
# 4. Profundidad, mediana y media recortada
#------------------------------------------------------------------------------------
# Verificación de la fórmula de rangos contra la definición por pares (40 curvas)
Yp <- Ye[1:40, ]; pares <- combn(40, 2)
bruta <- sapply(1:40, function(i) mean(apply(pares, 2, function(p)
  mean(pmin(Yp[p[1], ], Yp[p[2], ]) <= Yp[i, ] & Yp[i, ] <= pmax(Yp[p[1], ], Yp[p[2], ])))))
cat("MBD por rangos vs fuerza bruta (todas las bandas de 2 curvas):", format(max(abs(profundidades(Yp)$MBD - bruta)), digits = 2), "\n")

pr_e <- profundidades(Ye); pr_s <- profundidades(Xs)
centralidad <- function(Y, pr) {
  incluidas <- which(pr$FM >= quantile(pr$FM, .25, names = FALSE))        # 75 % más profundo según FM
  list(mediana_FM = Y[which.max(pr$FM), ], mediana_MBD = Y[which.max(pr$MBD), ],
       recortada_FM_25 = colMeans(Y[incluidas, , drop = FALSE]))
}
ce_e <- centralidad(Ye, pr_e); ce_s <- centralidad(Xs, pr_s)
resumen_e <- data.frame(dia = t, media = media_e, de = me$de, cv = me$cv, derivada_media = dmu_e,
                        mediana_FM = ce_e$mediana_FM, mediana_MBD = ce_e$mediana_MBD, media_recortada_FM_25 = ce_e$recortada_FM_25)
resumen_s <- data.frame(dia = t, media = media_s, de = ms$de, cv = ms$cv, derivada_media = dmu_s,
                        mediana_FM = ce_s$mediana_FM, mediana_MBD = ce_s$mediana_MBD, media_recortada_FM_25 = ce_s$recortada_FM_25)
guardar_csv(resumen_e, "resumen_diario_estaciones.csv"); guardar_csv(resumen_s, "resumen_diario_satelite.csv")
cat("Curva más profunda (MBD) — estaciones:", paste(ce$STATION_ID[which.max(pr_e$MBD)], ce$YEAR[which.max(pr_e$MBD)]),
    "| satélite: año", cs$anio[which.max(pr_s$MBD)], "píxel", cs$fila[which.max(pr_s$MBD)], cs$col[which.max(pr_s$MBD)], "\n")
figura("centralidad.pdf", function() {
  par(mfrow = c(1, 2))
  for (d in list(list(resumen_e, "Estaciones"), list(resumen_s, "Satélite (tierra)"))) {
    matplot(t, d[[1]][, c("media", "mediana_FM", "mediana_MBD", "media_recortada_FM_25")], type = "l", lty = c(1, 2, 3, 1),
            col = c(4, 2, 1, 3), lwd = 2, xlab = "Día del año", ylab = "mm/día", main = paste("Centralidad:", d[[2]]))
    legend("topright", c("Media", "Mediana FM", "Mediana MBD", "Media recortada FM 25 %"), col = c(4, 2, 1, 3),
           lty = c(1, 2, 3, 1), lwd = 2, bty = "n", cex = .75)
  }
})

#------------------------------------------------------------------------------------
# 5. Boxplot funcional y curvas atípicas
#------------------------------------------------------------------------------------
# Con todas las estaciones juntas, las atípicas son casi todas del Pacífico (llueve mucho más que en el valle):
# por eso también se hace el boxplot solo con las estaciones del Cauca, que es la comparación útil.
bx_e  <- boxplot_funcional(Ye, prof = pr_e$MBD)
bx_ec <- boxplot_funcional(Ye[cuenca == "Cauca", ])
bx_s  <- boxplot_funcional(Xs, prof = pr_s$MBD)
cat(sprintf("Atípicas — todas las estaciones: %d de %d (%d del Pacífico, en %d estaciones) | Cauca: %d de %d | satélite: %d de %d\n",
            sum(bx_e$atip), nrow(Ye), sum(bx_e$atip & cuenca == "Pacífico"), length(unique(ce$STATION_ID[bx_e$atip])),
            sum(bx_ec$atip), sum(cuenca == "Cauca"), sum(bx_s$atip), nrow(Xs)))
if (requireNamespace("fda", quietly = TRUE)) {
  fb <- fda::fbplot(t(Ye), t, depth = pr_e$MBD, plot = FALSE)
  cat("Atípicas iguales a fda::fbplot:", identical(sort(as.integer(fb$outpoint)), which(bx_e$atip)), "\n")
}
atipicas <- data.frame(ce[bx_e$atip, 1:2], cuenca = cuenca[bx_e$atip], MBD = pr_e$MBD[bx_e$atip],
                       media_anual_mm = round(rowMeans(Ye[bx_e$atip, , drop = FALSE]) * 365))
at_cauca <- ce[cuenca == "Cauca", 1:2][bx_ec$atip, ]
guardar_csv(atipicas, "atipicas_todas_las_estaciones.csv"); guardar_csv(at_cauca, "atipicas_estaciones_cauca.csv")
mostrar_tabla(at_cauca, "Curvas atípicas de las estaciones del Cauca", 10L)
guardar_csv(data.frame(ce[, 1:2], cuenca = cuenca, FM = pr_e$FM, MBD = pr_e$MBD, atipica = bx_e$atip, tem_e),
            "profundidad_y_fechas_estaciones.csv")

# Años atípicos: boxplot de las curvas medias anuales (estaciones del Cauca y satélite)
anios_e <- sort(unique(ce$YEAR))
Ye_anual <- t(sapply(anios_e, function(a) colMeans(Ye[cuenca == "Cauca" & ce$YEAR == a, , drop = FALSE])))
Xs_anual <- t(sapply(ANIOS, function(a) colMeans(Xs[cs$anio == a, , drop = FALSE])))
bx_ya <- boxplot_funcional(Ye_anual); bx_sa <- boxplot_funcional(Xs_anual)
oni <- read.csv("eda_fase0/oni_anual_fase.csv", fileEncoding = "UTF-8")
prof_anual <- merge(data.frame(anio = anios_e, MBD_estaciones = round(bx_ya$prof, 3), atipico_estaciones = bx_ya$atip),
                    data.frame(anio = ANIOS, MBD_satelite = round(bx_sa$prof, 3), atipico_satelite = bx_sa$atip))
prof_anual <- merge(prof_anual, oni, by = "anio")
guardar_csv(prof_anual, "profundidad_anual.csv")
mostrar_tabla(prof_anual[order(prof_anual$MBD_estaciones), ], "Años menos profundos (estaciones del Cauca)")

figura("boxplot_funcional.pdf", function() {
  par(mfrow = c(2, 2))
  dibujar_boxplot(Ye, bx_e, paste("Todas las estaciones:", sum(bx_e$atip), "atípicas"))
  dibujar_boxplot(Ye[cuenca == "Cauca", ], bx_ec, paste("Estaciones del Cauca:", sum(bx_ec$atip), "atípicas"))
  dibujar_boxplot(Xs, bx_s, paste("Satélite (tierra):", sum(bx_s$atip), "atípicas"))
  dibujar_boxplot(Ye_anual, bx_ya, "Curvas medias anuales (estaciones del Cauca)")
}, width = 12, height = 8)

#------------------------------------------------------------------------------------
# 6. Contexto: cuenca y fase ENSO (descriptivo, no causal)
#------------------------------------------------------------------------------------
fase_e <- oni$fase[match(ce$YEAR, oni$anio)]; fase_s <- oni$fase[match(cs$anio, oni$anio)]
stopifnot(!anyNA(fase_e), !anyNA(fase_s))
medias_grupo <- function(Y, g) sapply(sort(unique(g)), function(k) colMeans(Y[g == k, , drop = FALSE]))
enso <- data.frame(fase = sort(unique(fase_e)),
  estaciones_Cauca_mm_anio = round(tapply(rowMeans(Ye[cuenca == "Cauca", ]) * 365, fase_e[cuenca == "Cauca"], mean)),
  satelite_Cauca_mm_anio = round(tapply(rowMeans(Xs[region_sat == "Cauca", ]) * 365, fase_s[region_sat == "Cauca"], mean)))
guardar_csv(enso, "lluvia_anual_por_fase_ENSO.csv"); mostrar_tabla(enso, "Lluvia media anual por fase ENSO (vertiente del Cauca)")
figura("contexto_cuenca_ENSO.pdf", function() {
  par(mfrow = c(1, 3))
  for (d in list(list(medias_grupo(Ye, cuenca), "Estaciones por cuenca"), list(medias_grupo(Xs, region_sat), "Satélite por vertiente"))) {
    matplot(t, d[[1]], type = "l", lty = 1, lwd = 2, col = c(4, "darkgreen"), xlab = "Día del año", ylab = "mm/día", main = d[[2]])
    legend("topright", colnames(d[[1]]), col = c(4, "darkgreen"), lwd = 2, bty = "n")
  }
  m <- medias_grupo(Ye[cuenca == "Cauca", ], fase_e[cuenca == "Cauca"])
  matplot(t, m, type = "l", lty = 1, lwd = 2, col = c(2, 1, 4), xlab = "Día del año", ylab = "mm/día",
          main = "Estaciones del Cauca por fase ENSO (ONI anual)")
  legend("topright", colnames(m), col = c(2, 1, 4), lwd = 2, bty = "n")
})

#------------------------------------------------------------------------------------
# 7. Verificación cruzada con Python
#------------------------------------------------------------------------------------
py <- read.csv("fda_fase2/profundidad_y_fechas_estaciones.csv")
py_tem <- as.matrix(py[, colnames(tem_e)])
cat("FM R vs Python:", format(max(abs(pr_e$FM - py$FM)), digits = 2), "| MBD R vs Python:", format(max(abs(pr_e$MBD - py$MBD)), digits = 2), "\n")
cat("Fechas de temporadas iguales a Python:", round(mean(tem_e == py_tem, na.rm = TRUE), 4),
    "| mismos faltantes:", identical(is.na(tem_e), is.na(py_tem)), "\n")
py_at <- read.csv("fda_fase2/atipicas_todas_las_estaciones.csv"); py_atc <- read.csv("fda_fase2/atipicas_estaciones_cauca.csv")
cat("Mismas atípicas que Python — todas las estaciones:", setequal(paste(atipicas$STATION_ID, atipicas$YEAR), paste(py_at$STATION_ID, py_at$YEAR)),
    "| Cauca:", setequal(paste(at_cauca$STATION_ID, at_cauca$YEAR), paste(py_atc$STATION_ID, py_atc$YEAR)), "\n")

diagnostico <- data.frame(fuente = c("estaciones", "satelite_tierra"), curvas = c(nrow(Ye), nrow(Xs)), K = c(K_EST, K_SAT),
  atipicas = c(sum(bx_e$atip), sum(bx_s$atip)), minimo_mm_dia = c(min(Ye), min(Xs)),
  pct_puntos_negativos = 100 * c(mean(Ye < 0), mean(Xs < 0)), dia_mayor_DE = c(which.max(me$de), which.max(ms$de)),
  cor_dias_105_290 = c(me$cor[105, 290], ms$cor[105, 290]), lluvia_media_mm_anio = round(c(mean(media_e), mean(media_s)) * 365))
guardar_csv(diagnostico, "diagnosticos_EDA.csv"); mostrar_tabla(diagnostico, "Diagnóstico final del EDA funcional")
cat("EDA funcional en R terminado. Figuras en", DIR_FIG, "y tablas en", DIR_TAB, "\n")

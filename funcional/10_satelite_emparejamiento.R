# =============================================================================
# 10. SATÉLITE: LECTURA, mm/día, RECORTE AL VALLE, AÑOS DEFECTUOSOS Y EMPAREJAMIENTO
# =============================================================================
# Entradas: datos_curso/ImagenesSatelitePentadal/ppt{AAAA}{PP}.tif (mm por péntada)
#           salidas_R_datos/lluvia_valle_limpia.rds
#           datos_web/limites/deptos_valle_y_vecinos_DANE_2020.gpkg
# Salida:   salidas_R_funcional/satelite.rds
# =============================================================================

source("funcional/00_funciones.R")
suppressPackageStartupMessages({ library(terra); library(sf) })
sf_use_s2(FALSE)
titulo("10. Satélite y emparejamiento estación-píxel")

lluvia <- readRDS("salidas_R_datos/lluvia_valle_limpia.rds")
est <- lluvia$estaciones

archivos <- data.table(ruta = list.files("datos_curso/ImagenesSatelitePentadal", "^ppt\\d{6}\\.tif$", full.names = TRUE))
archivos[, `:=`(anio = as.integer(substr(basename(ruta), 4, 7)), pentada = as.integer(substr(basename(ruta), 8, 9)))]
archivos <- archivos[anio %in% ANIOS][order(anio, pentada)]
stopifnot(nrow(archivos) == length(ANIOS) * 72, archivos[, .N, by = anio][, all(N == 72)])

pila <- rast(archivos$ruta)
V <- values(pila)
cat("Malla:", nrow(pila), "filas x", ncol(pila), "columnas =", ncell(pila), "píxeles de", res(pila)[1], "grados |",
    "péntadas:", ncol(V), "| faltantes:", sum(is.na(V)), "| valor máximo:", max(V), "mm/péntada\n")

archivos[, `:=`(mes = (pentada - 1L) %/% 6L + 1L, k = (pentada - 1L) %% 6L)]
archivos[, dias_mes := as.integer(format(as.Date(sprintf("%d-%02d-01", anio + (mes == 12L), mes %% 12L + 1L)) - 1L, "%d"))]
archivos[, dias := fifelse(k < 5L, 5L, dias_mes - 25L)]
stopifnot(archivos[k == 5 & mes == 2, all(dias == fifelse(anio %% 4 == 0, 4L, 3L))], archivos[k == 5 & mes == 1, all(dias == 6L)])
dias_mes_365 <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
ini_p <- 5 * (0:71 %% 6) + 1
fin_p <- ifelse(0:71 %% 6 < 5, ini_p + 4, dias_mes_365[0:71 %/% 6 + 1])
T_SAT <- c(0, cumsum(dias_mes_365))[0:71 %/% 6 + 1] + (ini_p + fin_p) / 2
tope <- V >= 199.99
Vd <- sweep(V, 2, archivos$dias, "/")

xy <- as.data.table(xyFromCell(pila, seq_len(ncell(pila))))
xy[, celda := seq_len(.N)][, `:=`(fila = rowFromCell(pila, celda), col = colFromCell(pila, celda))]
deptos <- st_read("datos_web/limites/deptos_valle_y_vecinos_DANE_2020.gpkg", quiet = TRUE)
valle_pol <- deptos[enc2utf8(deptos$shapeName) == "Valle del Cauca", ]
xy[, en_valle := lengths(suppressMessages(st_intersects(st_as_sf(xy, coords = c("x", "y"), crs = 4326), valle_pol))) > 0]
cat("Píxeles con centro dentro del polígono DANE del Valle:", xy[en_valle == TRUE, .N], "de", nrow(xy), "\n")

est[, celda := cellFromXY(pila, cbind(lon, lat))]
stopifnot(!anyNA(est$celda))
est <- xy[, .(celda, fila, col, lon_px = x, lat_px = y, px_en_valle = en_valle)][est, on = "celda"]
est[, dist_km := 12742 * asin(sqrt(sin((lat_px - lat) * pi / 360)^2 + cos(lat * pi / 180) * cos(lat_px * pi / 180) * sin((lon_px - lon) * pi / 360)^2))]
est[, comparten := .N, by = celda]
emparejamiento <- est[, .(codigo, nombre, municipio, vertiente, lat, lon, altitud, nivel_c, revisar,
                          celda, fila, col, lon_px, lat_px, px_en_valle, dist_km = round(dist_km, 2), estaciones_en_el_pixel = comparten)]
tabla_f(emparejamiento[order(celda, codigo)], "10_emparejamiento_estacion_pixel")
compartidos <- emparejamiento[estaciones_en_el_pixel > 1, .(estaciones = .N, codigos = paste(codigo, collapse = " / "),
                                                            nombres = paste(nombre, collapse = " / ")), by = .(celda, fila, col)]
tabla_f(compartidos, "10_pixeles_compartidos")
cat("Estaciones:", nrow(emparejamiento), "| píxeles distintos:", uniqueN(emparejamiento$celda),
    "| píxeles con más de una estación:", nrow(compartidos), "(", sum(compartidos$estaciones), "estaciones ) |",
    "distancia al centro del píxel: mediana", round(median(emparejamiento$dist_km), 2), "km, máx", round(max(emparejamiento$dist_km), 2), "km |",
    "píxel con centro fuera del polígono:", emparejamiento[px_en_valle == FALSE, .N], "\n")
print(compartidos[order(-estaciones)][1:5])

celdas_usar <- sort(unique(c(xy[en_valle == TRUE, celda], est$celda)))
n_an <- length(ANIOS)
X <- matrix(NA_real_, n_an * length(celdas_usar), 72)
n_tope <- integer(nrow(X))
for (a in seq_len(n_an)) {
  cols <- which(archivos$anio == ANIOS[a]); filas <- (a - 1) * length(celdas_usar) + seq_along(celdas_usar)
  X[filas, ] <- Vd[celdas_usar, cols]; n_tope[filas] <- rowSums(tope[celdas_usar, cols])
}
meta_sat <- data.table(anio = rep(ANIOS, each = length(celdas_usar)), celda = rep(celdas_usar, n_an))
meta_sat <- xy[, .(celda, fila, col, lon = x, lat = y, en_valle)][meta_sat, on = "celda"]
meta_sat[, `:=`(n_tope = n_tope, particion = PARTICION(anio), defectuoso = anio %in% ANIOS_SAT_DEFECTUOSOS,
                con_estacion = celda %in% est$celda, curva_id = paste(anio, celda, sep = "_"))]
stopifnot(!anyNA(X))
cat("Curvas satelitales:", nrow(X), "=", n_an, "años x", length(celdas_usar), "píxeles (Valle + píxeles de estación fuera del polígono) |",
    "con alguna péntada en el tope de 200 mm:", sum(n_tope > 0), sprintf("(%.1f %%)", 100 * mean(n_tope > 0)), "\n")

d <- lluvia$diario
anual_est <- d[, .(P_est = mean(P, na.rm = TRUE) * 365), by = .(codigo, anio)]
anual_est <- est[, .(codigo, celda)][anual_est, on = "codigo"]
anual_px <- meta_sat[, .(anio, celda, P_sat = rowMeans(X) * 365)]
anual_par <- anual_px[anual_est, on = .(anio, celda)]
anual_par[, cociente := P_sat / P_est]
por_anio <- anual_par[, .(estaciones = .N, P_est_mediana = round(median(P_est)), P_sat_mediana = round(median(P_sat)),
                          cociente_mediano = round(median(cociente), 3)), by = anio][order(anio)]
valle_anual <- meta_sat[en_valle == TRUE, .(P_sat_valle = round(mean(rowMeans(X[.I, , drop = FALSE])) * 365)), by = anio]
por_anio <- valle_anual[por_anio, on = "anio"]
por_anio[, `:=`(z_cociente = round((cociente_mediano - median(cociente_mediano)) / mad(cociente_mediano), 2),
                defectuoso = anio %in% ANIOS_SAT_DEFECTUOSOS, particion = PARTICION(anio))]
tabla_f(por_anio, "10_satelite_vs_estaciones_por_anio")
print(por_anio)

grafico_f(ggplot(por_anio, aes(anio, cociente_mediano, fill = defectuoso)) + geom_col() +
            geom_hline(yintercept = median(por_anio$cociente_mediano), linetype = 2) +
            scale_fill_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey60"), labels = c(`TRUE` = "1985, 1988, 2000", `FALSE` = "resto"), name = NULL) +
            labs(x = NULL, y = "Mediana de satélite / estación (lluvia anual)", title = "Lluvia anual del píxel frente a la de su estación"),
          "10_cociente_satelite_estacion")

media_px <- meta_sat[, .(P_media = mean(rowMeans(X[.I, , drop = FALSE])) * 365, n_tope = sum(n_tope)), by = .(celda, lon, lat, en_valle)]
grafico_f(ggplot() + geom_tile(data = media_px, aes(lon, lat, fill = P_media, alpha = en_valle)) +
            geom_sf(data = valle_pol, fill = NA, colour = "black", linewidth = .4) +
            geom_point(data = est, aes(lon, lat), size = .9, colour = "white") +
            scale_fill_viridis_c(trans = "log10", name = "mm/año") + scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = .35), guide = "none") +
            labs(title = "Lluvia media anual del satélite 1981–2003 y estaciones incluidas", x = NULL, y = NULL),
          "10_mapa_satelite_estaciones", ancho = 7, alto = 7)

grafico_f(ggplot(emparejamiento, aes(dist_km)) + geom_histogram(binwidth = .25, boundary = 0, fill = "steelblue") +
            labs(x = "Distancia de la estación al centro de su píxel (km)", y = "Estaciones", title = "Emparejamiento: celda que contiene la estación"),
          "10_distancia_estacion_pixel", ancho = 7, alto = 4)

# QUÉ DICEN LOS DATOS
#   La malla es un rectángulo de 39 x 37 píxeles de 0,05° (~5,5 km) sin faltantes;
#   676 píxeles tienen el centro en el polígono DANE del Valle y se suman 5 píxeles
#   de estaciones cercanas al borde: 681 píxeles x 23 años = 15 663 curvas. Los
#   valores vienen en mm por péntada con tope en 200 mm (7,7 % de las curvas lo
#   tocan alguna vez, todas al oeste de -76,8°, en la vertiente Pacífica); al dividir por los días
#   reales de cada péntada (5, o 3 a 6 en la última del mes) quedan en mm/día.
#   Emparejamiento: cada estación toma la celda que la contiene. Las 176 caen en
#   148 píxeles; 23 píxeles tienen 2 a 4 estaciones (51 estaciones, 16 de ellas en
#   Cali), que comparten la misma curva satelital. El centro
#   del píxel queda a 2,2 km en mediana y nunca a más de 3,9 km.
#   Años defectuosos: 1985, 1988 y 2000 son los tres años más secos del satélite
#   en el Valle (2 219, 2 302 y 2 374 mm frente a ~2 730 de mediana). Frente a su
#   estación, 1988 y 2000 están entre los cocientes más bajos (0,88 y 0,92), pero
#   1985 no (1,03) porque las estaciones también fueron secas ese año; 1984 es tan
#   bajo como 1988 sin estar marcado. Se mantiene D-16: se marcan y se excluyen del
#   entrenamiento de los modelos; en el EDA se conservan. En casi todos los años el
#   satélite da MÁS lluvia que su estación (cociente mediano 1,16).

saveRDS(list(X = X, meta = meta_sat, t_sat = T_SAT, dias_pentada = archivos[, .(anio, pentada, dias)],
             emparejamiento = emparejamiento, por_anio = por_anio, celdas = xy),
        file.path(SALIDA_F, "satelite.rds"), compress = "xz")

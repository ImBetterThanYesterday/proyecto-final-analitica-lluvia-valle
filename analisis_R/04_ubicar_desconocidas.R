# =============================================================================
# 04_ubicar_desconocidas.R
# ¿Algunas de las 65 estaciones sin ubicación son del Valle?
# Tres pistas independientes:
#   1. Datos útiles: ¿pasarían la regla de cobertura (D-32)? Si no, no importa.
#   2. Código: los códigos CVC de 10 dígitos codifican la cuenca (los primeros
#      dígitos = zona/subzona/cuenca). Se mira el departamento de las estaciones
#      conocidas con el mismo prefijo.
#   3. Lluvia: la correlación de anomalías mensuales cae con la distancia
#      (D-06). Las estaciones conocidas más correlacionadas indican dónde está;
#      su centroide ponderado da una ubicación aproximada.
# Salida: analisis_R/salidas/desconocidas_ubicacion_inferida.csv
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(sf) })
sf_use_s2(FALSE)

est <- fread("analisis_R/salidas/estaciones_metadatos.csv", colClasses = list(character = "codigo"))
d <- rbind(fread("P_1980_1990.csv"), fread("P_1991_2003.csv"))[YEAR %between% c(1981, 2003)]
d[, codigo := as.character(STATION_ID)]
d <- d[!(MONTH == 2 & DAY == 29)]

# --- 1. Calidad (regla D-32 propuesta: ≥90 % de días, hueco máximo ≤15) ------
d[, fecha := as.IDate(sprintf("%d-%02d-%02d", YEAR, MONTH, DAY))]
setkey(d, codigo, fecha)
hueco_max <- function(f, anio) {
  ini <- as.IDate(sprintf("%d-01-01", anio)); fin <- as.IDate(sprintf("%d-12-31", anio))
  diffs <- diff(c(ini - 1L, f, fin + 1L)) - 1L
  max(diffs)
}
cal <- d[, .(n = .N, hmax = hueco_max(fecha, YEAR[1])), by = .(codigo, anio = YEAR)]
cal[, valido := n >= 0.9 * 365 & hmax <= 15]
util <- cal[, .(anios_validos = sum(valido), validos_entren = sum(valido & anio <= 1995)), by = codigo]
est <- util[est, on = "codigo"]
est[is.na(anios_validos), `:=`(anios_validos = 0L, validos_entren = 0L)]
est[, pasa := validos_entren >= 7]

# --- 2. Prefijo del código ----------------------------------------------------
conocidas <- est[departamento != "Desconocido"]
prefijo <- function(cod, k) {
  ref <- conocidas[substr(codigo, 1, k) == substr(cod, 1, k)]
  if (!nrow(ref)) return(c(NA, 0, NA))
  t <- ref[, .N, by = departamento][order(-N)]
  c(paste(sprintf("%s:%d", t$departamento, t$N), collapse = " "), nrow(ref),
    paste(unique(ref$municipio)[1:min(4, uniqueN(ref$municipio))], collapse = ", "))
}

# --- 3. Huella de lluvia -----------------------------------------------------
mm <- d[, .(P = mean(VALUE), n = .N), by = .(codigo, YEAR, MONTH)][n >= 25]
mm[, anom := log1p(P) - mean(log1p(P)), by = .(codigo, MONTH)]
W <- dcast(mm, YEAR + MONTH ~ codigo, value.var = "anom")
M <- as.matrix(W[, -(1:2)])
C <- cor(M, use = "pairwise.complete.obs")
N <- crossprod(!is.na(M))

desc <- est[departamento == "Desconocido"]
con_xy <- conocidas[!is.na(lat) & codigo %in% colnames(C)]
res <- rbindlist(lapply(desc$codigo, function(cod) {
  p8 <- prefijo(cod, 7); p6 <- prefijo(cod, 6); p4 <- prefijo(cod, 4)
  out <- data.table(codigo = cod, prefijo7 = p8[1], n7 = as.integer(p8[2]), municipios7 = p8[3],
                    prefijo6 = p6[1], prefijo4 = p4[1])
  if (!cod %in% colnames(C)) return(out)
  ok <- con_xy$codigo[N[cod, con_xy$codigo] >= 48]
  if (length(ok) < 5) return(out)
  r <- sort(C[cod, ok], decreasing = TRUE)[1:5]
  top <- con_xy[match(names(r), codigo)]
  w <- pmax(r, 0)^4
  # ¿Gemela? Lluvia DIARIA comparada con la conocida más parecida
  dd <- merge(d[codigo == cod, .(fecha, a = VALUE)], d[codigo == top$codigo[1], .(fecha, b = VALUE)], by = "fecha")
  out[, `:=`(mas_parecida = sprintf("%s %s (%s)", top$codigo[1], top$nombre[1], top$municipio[1]),
             r_diaria = round(cor(dd$a, dd$b), 2), pct_dias_iguales = round(100 * mean(dd$a == dd$b), 1),
             dias_comunes = nrow(dd))]
  out[, `:=`(r_max = round(r[1], 2), top5 = paste(sprintf("%s (%s, %.2f)", top$nombre, top$municipio, r), collapse = "; "),
             top5_valle = sum(top$departamento == "Valle del Cauca"),
             lat_est = sum(w * top$lat) / sum(w), lon_est = sum(w * top$lon) / sum(w),
             dispersion_km = round(max(sqrt(((top$lat - top$lat[1]) * 111)^2 + ((top$lon - top$lon[1]) * 111)^2)), 1))]
  out
}), fill = TRUE)

# Departamento del centroide estimado
deptos <- st_read("datos_web/limites/deptos_valle_y_vecinos_DANE_2020.gpkg", quiet = TRUE)
xy <- res[!is.na(lat_est)]
i <- st_intersects(st_as_sf(xy, coords = c("lon_est", "lat_est"), crs = 4326), deptos)
xy[, depto_centroide := vapply(i, function(k) if (length(k)) enc2utf8(deptos$shapeName[k[1]]) else NA_character_, "")]
res <- merge(res, xy[, .(codigo, depto_centroide)], by = "codigo", all.x = TRUE)
res <- merge(desc[, .(codigo, anios_con_dato, anios_validos, validos_entren, pasa)], res, by = "codigo")

# Veredicto
res[, pista_codigo := fcase(grepl("^Valle del Cauca:\\d+$", prefijo7), "Valle",
                            grepl("^Cauca:\\d+$", prefijo7), "Cauca",
                            is.na(prefijo7) & grepl("^Valle del Cauca:\\d+$", prefijo6), "Valle (6 díg.)",
                            is.na(prefijo7) & grepl("^Cauca:\\d+$", prefijo6), "Cauca (6 díg.)",
                            default = "Mixta/sin referencia")]
res[, pista_lluvia := fcase(is.na(r_max), "Sin datos suficientes",
                            r_max < 0.6, "Débil (r<0.6)",
                            top5_valle >= 4, "Valle", top5_valle <= 1, "Cauca", default = "Mixta")]
res[, veredicto := fcase(
  startsWith(pista_codigo, "Valle") & pista_lluvia == "Valle", "Probable Valle (código + lluvia)",
  startsWith(pista_codigo, "Cauca") & pista_lluvia == "Cauca", "Probable Cauca (código + lluvia)",
  default = "Indeterminada")]
setorder(res, -pasa, veredicto, codigo)
fwrite(res, "analisis_R/salidas/desconocidas_ubicacion_inferida.csv")

cat("Desconocidas:", nrow(res), "| pasarían la regla de calidad (≥7 años válidos en entrenamiento):", sum(res$pasa), "\n\n")
print(res[, .N, by = .(pasa, veredicto)][order(-pasa, veredicto)])
cat("\nLas que pasarían calidad:\n")
print(res[pasa == TRUE, .(codigo, ent = validos_entren, prefijo7, r_max, top5_valle, depto_centroide,
                          disp_km = dispersion_km, veredicto)])
cat("\nPosibles gemelas (lluvia diaria casi igual a una conocida):\n")
print(res[r_diaria >= 0.8, .(codigo, pasa, mas_parecida, r_mensual = r_max, r_diaria, pct_dias_iguales, dias_comunes)])

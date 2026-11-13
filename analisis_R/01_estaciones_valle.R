# =============================================================================
# 01_estaciones_valle.R
# Población del EDA: las 232 estaciones del Valle de D-34/D-35
# (datos_web/estaciones_valle_finales.csv), con su cobertura cruda en 1981–2003.
# También deja el departamento de los 344 códigos (inventario) para contexto.
#
# Comprobaciones:
#   1. Cada coordenada cae dentro del polígono DANE del Valle (MGN 2020).
#   2. Distancia al borde con otro departamento (EPSG:9377).
#
# Salida: analisis_R/salidas/estaciones_metadatos.csv  (una fila por código con
#         datos en 1981–2003; `departamento` = "Valle del Cauca" solo para las 232)
# Ejecutar desde la raíz:  Rscript analisis_R/01_estaciones_valle.R
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(readxl); library(sf) })
sf_use_s2(FALSE)
dir.create("analisis_R/salidas", showWarnings = FALSE, recursive = TRUE)

# --- Cobertura cruda 1981–2003 ------------------------------------------------
crudo <- rbind(fread("P_1980_1990.csv"), fread("P_1991_2003.csv"))[YEAR %between% c(1981, 2003)]
crudo[, codigo := as.character(STATION_ID)]
cob <- crudo[, .(dias_con_dato = .N, anio_ini = min(YEAR), anio_fin = max(YEAR),
                 anios_con_dato = uniqueN(YEAR)), by = codigo]

# --- Inventario de los 344 códigos (D-31/D-32) --------------------------------
inv <- as.data.table(read_excel("datos_web/estaciones_valle_del_cauca.xlsx", sheet = 2))
inv <- inv[, .(codigo = format(inv[[1]], scientific = FALSE, trim = TRUE),
               depto_inventario = Departamento, nivel_inventario = `Nivel de evidencia`)]

# --- Las 232 finales (D-34/D-35) ----------------------------------------------
fin <- fread("datos_web/estaciones_valle_finales.csv", encoding = "UTF-8",
             colClasses = list(character = "codigo"))
fin <- fin[, .(codigo, nombre, municipio, lat, lon, altitud = alt, nivel_c, fuente,
               fuente_alt, srtm, revisar)]

m <- merge(cob, inv, by = "codigo", all.x = TRUE)
m <- merge(m, fin, by = "codigo", all = TRUE)   # all = TRUE: también finales sin datos 1981–2003
m[is.na(dias_con_dato), `:=`(dias_con_dato = 0L, anios_con_dato = 0L)]
m[, departamento := fifelse(!is.na(nivel_c), "Valle del Cauca",
                    fcase(depto_inventario == "Valle", "Valle (excluida, D-34)",
                          depto_inventario == "Cauca", "Cauca",
                          default = "Desconocido"))]

# --- Prueba geográfica --------------------------------------------------------
deptos <- st_read("datos_web/limites/deptos_valle_y_vecinos_DANE_2020.gpkg", quiet = TRUE)
deptos$shapeName <- enc2utf8(deptos$shapeName)
xy <- m[!is.na(lat)]
pts <- st_as_sf(xy, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
idx <- st_intersects(pts, deptos)
xy[, depto_geo := vapply(idx, function(i) if (length(i)) deptos$shapeName[i[1]] else "Fuera de los vecinos", "")]
Dd <- units::drop_units(st_distance(st_transform(pts, 9377), st_transform(deptos, 9377)))
propio <- match(xy$depto_geo, deptos$shapeName)
xy[, km_borde := round(vapply(seq_len(nrow(Dd)), function(i) min(Dd[i, -propio[i]]), 0) / 1000, 2)]
m <- merge(m, xy[, .(codigo, depto_geo, km_borde)], by = "codigo", all.x = TRUE)

setorder(m, departamento, codigo)
fwrite(m, "analisis_R/salidas/estaciones_metadatos.csv")

v <- m[departamento == "Valle del Cauca"]
cat("Códigos con datos 1981–2003:", m[dias_con_dato > 0, .N], "\n")
print(m[, .(codigos = .N, con_datos_81_03 = sum(dias_con_dato > 0), dias = sum(dias_con_dato)), by = departamento])
cat("\nLas 232 del Valle: polígono DANE\n"); print(v[, .N, by = depto_geo])
cat("\nSin datos en 1981–2003:", v[dias_con_dato == 0, .N], "\n")
cat("A ≤2 km de otro departamento:", v[km_borde <= 2, .N], "\n")
print(v[, .N, by = nivel_c])

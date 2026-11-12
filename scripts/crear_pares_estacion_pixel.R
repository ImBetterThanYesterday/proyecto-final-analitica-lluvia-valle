# Crosswalk local entre coordenadas CVC, matriz de curvas Jacobo y píxeles satelitales.
# No modifica los datos fuente. Los códigos sin coordenada quedan como no emparejados.

GRID <- list(
  xmin = -77.55, xmax = -75.70,
  ymin = 3.05, ymax = 5.00,
  resx = 0.05, resy = 0.05,
  nrow = 39L, ncol = 37L
)

haversine_km <- function(lon1, lat1, lon2, lat2) {
  rad <- pi / 180
  dlat <- (lat2 - lat1) * rad
  dlon <- (lon2 - lon1) * rad
  a <- sin(dlat / 2)^2 + cos(lat1 * rad) * cos(lat2 * rad) * sin(dlon / 2)^2
  6371.0088 * 2 * atan2(sqrt(pmin(1, a)), sqrt(pmax(0, 1 - a)))
}

mapear_a_pixel <- function(lon, lat, grid = GRID) {
  if (length(lon) != length(lat)) stop("lon y lat deben tener la misma longitud")
  en_malla <- is.finite(lon) & is.finite(lat) &
    lon >= grid$xmin & lon < grid$xmax & lat > grid$ymin & lat <= grid$ymax
  fila <- col <- rep(NA_integer_, length(lon))
  fila[en_malla] <- as.integer(floor((grid$ymax - lat[en_malla]) / grid$resy))
  col[en_malla] <- as.integer(floor((lon[en_malla] - grid$xmin) / grid$resx))
  lon_centro <- lat_centro <- distancia_km <- rep(NA_real_, length(lon))
  lon_centro[en_malla] <- grid$xmin + (col[en_malla] + 0.5) * grid$resx
  lat_centro[en_malla] <- grid$ymax - (fila[en_malla] + 0.5) * grid$resy
  distancia_km[en_malla] <- haversine_km(lon[en_malla], lat[en_malla],
                                         lon_centro[en_malla], lat_centro[en_malla])
  data.frame(en_malla, fila, col, lon_centro, lat_centro, distancia_km)
}

construir_pares <- function(estacion_anio, crosswalk, sat_keys) {
  stopifnot(all(c("STATION_ID", "YEAR") %in% names(estacion_anio)))
  stopifnot(all(c("STATION_ID", "tiene_cvc", "en_malla", "en_tierra", "fila", "col") %in% names(crosswalk)))
  stopifnot(all(c("YEAR", "fila", "col") %in% names(sat_keys)))
  if (anyDuplicated(estacion_anio[c("STATION_ID", "YEAR")])) stop("Claves estación-año duplicadas")
  if (anyDuplicated(sat_keys[c("YEAR", "fila", "col")])) stop("Claves satélite año-fila-col duplicadas")

  pares <- merge(estacion_anio, crosswalk, by = "STATION_ID", all.x = TRUE, sort = FALSE)
  claves_sat <- paste(sat_keys$YEAR, sat_keys$fila, sat_keys$col, sep = "_")
  clave_par <- paste(pares$YEAR, pares$fila, pares$col, sep = "_")
  pares$curva_satelital <- !is.na(pares$fila) & clave_par %in% claves_sat
  pares$estado <- ifelse(!pares$tiene_cvc | is.na(pares$tiene_cvc), "sin_coordenadas",
    ifelse(!pares$en_malla, "fuera_malla",
      ifelse(!pares$en_tierra, "pixel_oceano",
        ifelse(!pares$curva_satelital, "sin_curva_satelital", "pareada"))))
  pares
}

validar_grid_tif <- function(path, grid = GRID) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Se requiere terra para validar el georreferenciado del TIFF antes de emparejar")
  }
  if (!file.exists(path)) stop("Falta TIFF de referencia: ", path)
  r <- terra::rast(path)
  e <- terra::ext(r); rr <- terra::res(r)
  bounds <- c(terra::xmin(e), terra::xmax(e), terra::ymin(e), terra::ymax(e))
  esperado <- c(grid$xmin, grid$xmax, grid$ymin, grid$ymax)
  if (nrow(r) != grid$nrow || ncol(r) != grid$ncol || any(abs(bounds - esperado) > 1e-4) ||
      any(abs(rr - c(grid$resx, grid$resy)) > 1e-8) || !terra::is.lonlat(r)) {
    stop("La geometría del TIFF no coincide con la malla documentada; no se crean pares")
  }
  invisible(list(nrow = nrow(r), ncol = ncol(r), bounds = bounds, resolution = rr, crs = terra::crs(r)))
}

ejecutar_crosswalk <- function(root = ".") {
  files <- c(
    estaciones = file.path(root, "fda_fase1", "curvas_estaciones_mm.csv.gz"),
    satelite = file.path(root, "fda_fase1", "curvas_satelite_mm.csv.gz"),
    cvc = file.path(root, "CVC_Wheather_Stations_Long_Lat.csv"),
    catalogo = file.path(root, "CatalogoEstacionesElevacion.xlsx"),
    tiff = file.path(root, "ImagenesSatelitePentadal", "ppt198101.tif")
  )
  faltantes <- files[!file.exists(files)]
  if (length(faltantes)) stop("Faltan entradas locales: ", paste(names(faltantes), collapse = ", "))
  grid_meta <- validar_grid_tif(files[["tiff"]])

  cvc <- read.csv(files[["cvc"]], fileEncoding = "latin1", colClasses = c(CODIGO = "character"),
                  check.names = FALSE, stringsAsFactors = FALSE)
  required_cvc <- c("CODIGO", "STATION_NAME", "BASIN", "STATUS", "Lon", "lat")
  if (!all(required_cvc %in% names(cvc))) stop("El CSV CVC no trae las columnas requeridas: ",
                                                 paste(setdiff(required_cvc, names(cvc)), collapse = ", "))
  if (anyDuplicated(cvc$CODIGO)) stop("El CSV CVC tiene códigos duplicados")
  cvc$Lon <- as.numeric(cvc$Lon); cvc$lat <- as.numeric(cvc$lat)
  dta_path <- file.path(root, "CVC_Wheather_Stations_Long_Lat.dta")
  cvc_dta_audit <- list(run = FALSE, id_overlap = NA_integer_, lon_error_max = NA_real_, lat_error_max = NA_real_)
  if (file.exists(dta_path) && requireNamespace("haven", quietly = TRUE)) {
    dta <- haven::read_dta(dta_path)
    dta$station_id <- as.character(dta$station_id)
    dta_match <- match(cvc$CODIGO, dta$station_id)
    if (anyNA(dta_match)) stop("Hay códigos CVC del CSV ausentes en el DTA")
    lon_err <- max(abs(cvc$Lon - dta$lon[dta_match]), na.rm = TRUE)
    lat_err <- max(abs(cvc$lat - dta$lat[dta_match]), na.rm = TRUE)
    if (lon_err > 1e-4 || lat_err > 1e-4) stop("Las coordenadas CVC CSV y DTA no concuerdan")
    cvc_dta_audit <- list(run = TRUE, id_overlap = length(unique(dta$station_id)),
                          lon_error_max = lon_err, lat_error_max = lat_err)
  }

  if (!requireNamespace("readxl", quietly = TRUE)) stop("Se requiere readxl para leer el catálogo local de elevación")
  alt <- as.data.frame(readxl::read_excel(files[["catalogo"]], sheet = 1,
    col_types = c("text", "text", "text", "text", "numeric")), stringsAsFactors = FALSE)
  if (!all(c("codigo", "nombre", "municipio", "alt") %in% names(alt))) stop("Esquema inesperado en CatalogoEstacionesElevacion.xlsx")
  if (anyDuplicated(alt$codigo)) stop("El catálogo de elevación tiene códigos duplicados")

  ids <- read.csv(gzfile(files[["estaciones"]]),
                  colClasses = c("character", "integer", rep("NULL", 365)), check.names = FALSE)
  station_years <- unique(ids[c("STATION_ID", "YEAR")])
  station_ids <- unique(ids$STATION_ID)
  sat_keys <- read.csv(gzfile(files[["satelite"]]),
    colClasses = c("integer", "integer", "integer", rep("NULL", 72)), check.names = FALSE)
  sat_keys <- unique(sat_keys[c("anio", "fila", "col")])
  names(sat_keys)[names(sat_keys) == "anio"] <- "YEAR"
  if (anyDuplicated(station_years)) stop("Claves estación-año duplicadas en las curvas")

  idx <- match(station_ids, cvc$CODIGO)
  lon <- cvc$Lon[idx]; lat <- cvc$lat[idx]
  geo <- mapear_a_pixel(lon, lat)
  alt_idx <- match(station_ids, alt$codigo)
  map_base <- data.frame(
    STATION_ID = station_ids,
    tiene_cvc = !is.na(idx),
    STATION_NAME = cvc$STATION_NAME[idx],
    STATUS = cvc$STATUS[idx],
    BASIN = cvc$BASIN[idx],
    lon = lon, lat = lat,
    altitud_m = alt$alt[alt_idx],
    nombre_catalogo = alt$nombre[alt_idx],
    municipio_catalogo = alt$municipio[alt_idx],
    geo,
    stringsAsFactors = FALSE
  )
  land_key <- paste(sat_keys$fila, sat_keys$col, sep = "_")
  map_base$en_tierra <- ifelse(map_base$en_malla,
                               paste(map_base$fila, map_base$col, sep = "_") %in% land_key,
                               NA)
  map_base$estado_geo <- ifelse(!map_base$tiene_cvc, "sin_coordenadas",
    ifelse(!map_base$en_malla, "fuera_malla",
      ifelse(!map_base$en_tierra, "pixel_oceano", "coordenada_pareada")))
  pares <- construir_pares(station_years, map_base, sat_keys)

  out <- file.path(root, "fda_fase2", "tablas_R")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write.csv(map_base, file.path(out, "crosswalk_estacion_pixel.csv"), row.names = FALSE, na = "")
  write.csv(pares, file.path(out, "pares_estacion_pixel.csv"), row.names = FALSE, na = "")
  resumen <- data.frame(
    indicador = c("estaciones_en_curvas", "estaciones_con_CVC", "coordenadas_dentro_malla",
                  "pixeles_de_tierra", "estaciones_con_altitud", "estaciones_en_malla_con_altitud",
                  "pares_estacion_anio", "pares_estacion_pixel_validos", "pares_validos_con_altitud",
                  "pares_fuera_malla", "pares_sin_coordenada", "pares_pixel_oceano",
                  "codigos_CVC_no_presentes_en_curvas", "distancia_pixel_mediana_km", "distancia_pixel_maxima_km",
                  "cvc_dta_comprobacion_ejecutada", "cvc_dta_error_max_lon", "cvc_dta_error_max_lat"),
    valor = c(nrow(map_base), sum(map_base$tiene_cvc), sum(map_base$en_malla, na.rm = TRUE),
              sum(map_base$en_tierra, na.rm = TRUE), sum(is.finite(map_base$altitud_m)),
              sum(map_base$en_malla & is.finite(map_base$altitud_m), na.rm = TRUE), nrow(station_years),
              sum(pares$estado == "pareada"), sum(pares$estado == "pareada" & is.finite(pares$altitud_m)),
              sum(pares$estado == "fuera_malla"), sum(pares$estado == "sin_coordenadas"),
              sum(pares$estado == "pixel_oceano"), sum(!cvc$CODIGO %in% station_ids),
              median(map_base$distancia_km[map_base$en_tierra], na.rm = TRUE),
              max(map_base$distancia_km[map_base$en_tierra], na.rm = TRUE), as.integer(cvc_dta_audit$run),
              cvc_dta_audit$lon_error_max, cvc_dta_audit$lat_error_max)
  )
  cobertura_anual <- do.call(rbind, lapply(sort(unique(pares$YEAR)), function(y) {
    py <- pares[pares$YEAR == y, , drop = FALSE]
    data.frame(YEAR = y, estacion_anio = nrow(py), pares = sum(py$estado == "pareada"),
               sin_coordenada = sum(py$estado == "sin_coordenadas"), fuera_malla = sum(py$estado == "fuera_malla"),
               pixel_oceano = sum(py$estado == "pixel_oceano"), sin_curva_satelital = sum(py$estado == "sin_curva_satelital"))
  }))
  write.csv(resumen, file.path(out, "resumen_pares_estacion_pixel.csv"), row.names = FALSE)
  write.csv(cobertura_anual, file.path(out, "cobertura_pares_estacion_pixel.csv"), row.names = FALSE)
  cat("Raster verificado:", grid_meta$nrow, "x", grid_meta$ncol,
      "WGS84, resolución", grid_meta$resolution[1], "grados\n")
  print(resumen, row.names = FALSE)
  invisible(list(crosswalk = map_base, pares = pares, resumen = resumen))
}

if (sys.nframe() == 0L) ejecutar_crosswalk()

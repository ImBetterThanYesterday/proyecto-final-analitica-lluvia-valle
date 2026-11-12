root <- normalizePath(".")
source(file.path(root, "scripts", "crear_pares_estacion_pixel.R"), local = TRUE)

grid_test <- list(xmin = 0, xmax = 2, ymin = 0, ymax = 2, resx = 1, resy = 1, nrow = 2, ncol = 2)
mapped <- mapear_a_pixel(lon = c(0.5, 1.5, 2.1), lat = c(1.5, 0.5, 0.5), grid = grid_test)
stopifnot(identical(mapped$en_malla, c(TRUE, TRUE, FALSE)))
stopifnot(identical(mapped$fila[1:2], c(0L, 1L)), identical(mapped$col[1:2], c(0L, 1L)))
stopifnot(abs(mapped$distancia_km[1]) < 1e-12, is.na(mapped$fila[3]), is.na(mapped$col[3]))
stopifnot(abs(haversine_km(0, 0, 0, 0)) < 1e-12)

crosswalk <- data.frame(
  STATION_ID = c("A", "B", "C"),
  tiene_cvc = c(TRUE, TRUE, TRUE),
  en_malla = c(TRUE, TRUE, FALSE),
  en_tierra = c(TRUE, FALSE, NA),
  fila = c(0L, 0L, NA_integer_),
  col = c(1L, 2L, NA_integer_),
  distancia_km = c(0.2, 0.1, NA_real_),
  stringsAsFactors = FALSE
)
station_years <- data.frame(
  STATION_ID = c("A", "A", "B", "C", "D"),
  YEAR = c(1981L, 1982L, 1981L, 1981L, 1981L),
  stringsAsFactors = FALSE
)
satellite_keys <- data.frame(YEAR = 1981L, fila = 0L, col = 1L)
pairs <- construir_pares(station_years, crosswalk, satellite_keys)
stopifnot(nrow(pairs) == nrow(station_years))
stopifnot(pairs$estado[pairs$STATION_ID == "A" & pairs$YEAR == 1981L] == "pareada")
stopifnot(pairs$estado[pairs$STATION_ID == "A" & pairs$YEAR == 1982L] == "sin_curva_satelital")
stopifnot(pairs$estado[pairs$STATION_ID == "B"] == "pixel_oceano")
stopifnot(pairs$estado[pairs$STATION_ID == "C"] == "fuera_malla")
stopifnot(pairs$estado[pairs$STATION_ID == "D"] == "sin_coordenadas")
cat("OK: mapeo de celdas, distancia y estados de pares estación-píxel\n")

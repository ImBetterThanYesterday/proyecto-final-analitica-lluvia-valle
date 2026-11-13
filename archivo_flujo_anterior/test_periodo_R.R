# Ventana normativa y cobertura: casos que min/max no detectaría.
local({
  e <- new.env()
  source("soporte_pruebas_R.R", local = TRUE)
  cargar_funciones_prueba(e)
  rechaza <- function(expr) inherits(tryCatch({force(expr); NULL}, error = identity), "error")
  d <- list(id = data.frame(STATION_ID = c("001", "001", "002", "002"),
                            YEAR = c(1980, 1981, 2003, 2004)),
            Y = matrix(c(1, 2, NA, 4, 5, 6, 7, 8), 4), tiempos = c("1", "2"))
  f <- e$filtrar_periodo(d, "estaciones")
  stopifnot(identical(f$id, d$id[2:3, ]), identical(f$Y, d$Y[2:3, , drop = FALSE]))
  malo <- d; malo$id$YEAR[2] <- 1981.5
  stopifnot(rechaza(e$filtrar_periodo(malo, "estaciones")))
  malo <- d; names(malo$id)[2] <- "otra"
  stopifnot(rechaza(e$filtrar_periodo(malo, "estaciones")))
  malo <- d; malo$id$YEAR <- 2004:2007
  stopifnot(rechaza(e$filtrar_periodo(malo, "estaciones")))
  malo <- d; malo$id <- rbind(d$id[2, ], d$id[2:4, ])
  stopifnot(rechaza(e$filtrar_periodo(malo, "estaciones")))
  c <- e$cobertura_periodo(f, "estaciones")
  stopifnot(nrow(c$anual) == 23, c$anual$curvas[c$anual$anio == 1982] == 0,
            nrow(c$ausentes) == 44, sum(c$anual$valores_faltantes) == 1,
            c$anual$curvas[c$anual$anio == 1981] == 1)
  sat <- list(id = data.frame(anio = c("2003", "2004"), fila = 0, col = 1),
              Y = matrix(1:4, 2), tiempos = c("t1", "t2"))
  stopifnot(nrow(e$filtrar_periodo(sat, "satelite")$id) == 1)
  oni <- data.frame(anio = 1981:2004, fase = "Neutro")
  stopifnot(nrow(e$validar_oni(oni)) == 23,
            rechaza(e$validar_oni(oni[-2, ])),
            rechaza(e$validar_oni(rbind(oni, oni[1, ]))))
  carpeta <- tempfile("periodo_"); dir.create(carpeta)
  on.exit(unlink(carpeta, recursive = TRUE), add = TRUE)
  e$SALIDA <- carpeta
  saveRDS(list(estaciones = list(), satelite = list()), file.path(carpeta, "construccion.rds"))
  stopifnot(rechaza(e$eda())) # RDS histórico sin ventana declarada.
  saveRDS(list(periodo = 1981:2012), file.path(carpeta, "construccion.rds"))
  stopifnot(rechaza(e$eda()))
  cat("OK: filtrado de ambas fuentes, claves explícitas, rechazos, huecos anuales/sitio-año/valores y ONI completo\n")
})

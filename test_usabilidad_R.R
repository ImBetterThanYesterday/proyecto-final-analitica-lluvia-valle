# Regresión de importación exclusiva de tests y de dispositivos con datos sintéticos.
local({
  entorno <- new.env()
  source("soporte_pruebas_R.R", local = TRUE)
  cargar_funciones_prueba(entorno)
  stopifnot(is.function(entorno$nodos), !exists("est", entorno, inherits = FALSE),
            !exists("clase", entorno, inherits = FALSE), !exists("principal", entorno, inherits = FALSE))
  carpeta <- tempfile("usabilidad_")
  dir.create(carpeta)
  anterior <- getwd()
  on.exit({setwd(anterior); unlink(carpeta, recursive = TRUE)}, add = TRUE)
  setwd(carpeta)
  stopifnot(!dir.exists("resultados_R"), is.null(dev.list()))
  entorno$SALIDA <- carpeta
  entorno$exportar_grafico("solo.pdf", function() plot(1:3), mostrar = FALSE)
  stopifnot(file.info("solo.pdf")$size > 0, is.null(dev.list()))
  # PNG simula un dispositivo visible ya abierto, sin necesitar RStudio.
  png("pantalla.png")
  pantalla <- dev.cur()
  llamadas <- 0L
  dibujar <- function() {llamadas <<- llamadas + 1L; plot(1:3)}
  entorno$exportar_grafico("doble.pdf", dibujar, mostrar = TRUE)
  stopifnot(dev.cur() == pantalla, llamadas == 2L, file.info("doble.pdf")$size > 0)
  dev.off()
  stopifnot(file.info("pantalla.png")$size > 0)
  # Rutinas reales con datos sintéticos: misma matemática y ruta gráfica,
  # sin recalcular las decenas de miles de curvas del proyecto.
  set.seed(19)
  tiempo <- seq(1, 365, length.out = 21)
  datos <- list(id = data.frame(STATION_ID = rep(c("2601", "2701"), 4), YEAR = 1981:1988),
                Y = matrix(runif(8 * 21, 1, 20), 8))
  png("pantalla_fases.png", width = 1100, height = 700)
  pantalla <- dev.cur()
  opciones <- options(lluvia.mostrar = TRUE)
  on.exit(options(opciones), add = TRUE)
  objeto <- entorno$construir_fuente(datos, tiempo, 4:6, "estaciones")
  diagnostico <- entorno$eda_fuente(objeto, "estaciones",
    data.frame(anio = 1981:1988, fase = rep(c("neutra", "calida"), 4)))
  stopifnot(dev.cur() == pantalla, diagnostico$n == 8,
            all(file.info(c("seleccion_K_estaciones.pdf", "curva_banda_residuos_estaciones.pdf",
                            "EDA_estaciones.pdf", "contexto_estaciones.pdf"))$size > 0),
            nrow(read.csv("resumen_estaciones.csv")) == 365)
  dev.off()
  stopifnot(file.info("pantalla_fases.png")$size > 0)
  options(opciones)
  fallo <- try(entorno$exportar_grafico("fallo.pdf", function() stop("deliberado"), mostrar = FALSE), silent = TRUE)
  stopifnot(inherits(fallo, "try-error"), is.null(dev.list()))
  # Rama nueva: exportación simultánea PDF + PNG + dispositivo activo.
  options(lluvia.png = TRUE)
  png("pantalla_tres.png")
  pantalla <- dev.cur(); llamadas <- 0L
  entorno$exportar_grafico("tres.pdf", dibujar, mostrar = TRUE)
  stopifnot(llamadas == 3L, dev.cur() == pantalla,
            all(file.info(c("tres.pdf", "tres.png"))$size > 0))
  dev.off()
  llamadas <- 0L
  fallo_png <- try(entorno$exportar_grafico("fallo_png.pdf", function() {
    llamadas <<- llamadas + 1L
    if (llamadas == 2L) stop("fallo deliberado en PNG")
    plot(1:3)
  }, mostrar = FALSE), silent = TRUE)
  stopifnot(inherits(fallo_png, "try-error"), is.null(dev.list()))
  options(lluvia.png = NULL)
  write.csv(data.frame(x = 1:3), "tabla.csv", row.names = FALSE)
  stopifnot(identical(entorno$ver_tabla("tabla.csv", abrir = FALSE)$x, 1:3))
  cat("OK: arnés aislado, PDF y dispositivo visible, construcción/EDA sintéticos reales, limpieza por error y tablas\n")
})

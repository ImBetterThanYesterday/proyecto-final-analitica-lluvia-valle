# Verificación integral de los archivos reales, salida independiente del historial Python.
source("test_fda_R.R")
source("test_periodo_R.R")
entradas <- c("fda_fase1/curvas_estaciones_mm.csv.gz", "fda_fase1/curvas_satelite_mm.csv.gz",
              "fda_fase2/coef_satelite_tierra.csv.gz", "eda_fase0/oni_anual_fase.csv",
              "Proyecto_Final_Analitica.docx", "Proyecto_Final_Analitica_FeErratas.docx",
              "GUIA_RECONSTRUCCION_EDA.docx")
antes <- tools::md5sum(entradas)
# Source normal en sesión limpia de objetos del análisis, sin precargar funciones.
sesion <- new.env(parent = globalenv())
source("proyecto_funcional.R", local = sesion)
stopifnot(exists("est", sesion, inherits = FALSE), exists("diagnosticos", sesion, inherits = FALSE),
          sesion$e$K == 6L, sesion$s$K == 11L,
          sesion$eda_est$diagnostico$atipicas == 170L, sesion$eda_sat$diagnostico$atipicas == 0L)
# Los controles siguientes usan los auxiliares ya importados por la prueba matemática.
list2env(as.list(entorno), envir = environment())
stopifnot(identical(antes, tools::md5sum(entradas)))
obj <- readRDS("resultados_R/construccion.rds")
stopifnot(obj$estaciones$K == 6, obj$satelite$K == 11,
          identical(obj$periodo, 1981:2003),
          nrow(obj$estaciones$coef) == 4426, nrow(obj$satelite$coef) == 33189,
          length(obj$satelite$tierra) == 31119)
controles <- list()
anotar <- function(nombre, error, tol) {
  stopifnot(is.finite(error), error <= tol)
  controles[[length(controles) + 1L]] <<- data.frame(control = nombre, error = error, tolerancia = tol, estado = "OK")
}
for (fuente in c("estaciones", "satelite")) {
  o <- obj[[fuente]]
  d <- leer_matriz(paste0("fda_fase1/curvas_", fuente, "_mm.csv.gz"), if (fuente == "estaciones") 2L else 3L)
  # Referencia temporal independiente del helper productivo.
  anio <- as.integer(d$id[[if (fuente == "estaciones") "YEAR" else "anio"]])
  incluir <- anio >= 1981 & anio <= 2003
  d$id <- d$id[incluir, , drop = FALSE]; d$Y <- d$Y[incluir, , drop = FALSE]
  stopifnot(setequal(anio[incluir], 1981:2003))
  stopifnot(identical(o$id, d$id), identical(o$observados, !is.na(d$Y)))
  P <- base_bspline(o$tiempo, o$K)
  # Resolver de nuevo cada curva con QR, independiente de las ecuaciones normales.
  error <- 0
  for (i in seq_len(nrow(d$Y))) {
    m <- !is.na(d$Y[i, ])
    ref <- qr.solve(P[m, , drop = FALSE], d$Y[i, m])
    error <- max(error, max(abs(ref - o$coef[i, ])))
  }
  anotar(paste("Coeficientes todas las curvas contra QR", fuente), error, 1e-9)
  h <- read.csv(paste0("resultados_R/controles_H_", fuente, ".csv"))
  anotar(paste("H todas las máscaras", fuente), max(h$error_max), 1e-8)
  if (fuente == "satelite") {o$coef <- o$coef[o$tierra, , drop = FALSE]; o$id <- o$id[o$tierra, , drop = FALSE]}
  nombre <- if (fuente == "satelite") "satelite_tierra" else fuente
  medias_grupo <- read.csv(paste0("resultados_R/media_grupos_descriptivos_", nombre, ".csv"))
  stopifnot(sum(medias_grupo$n_curvas[medias_grupo$dia == 1]) == nrow(o$coef))
  anual <- read.csv(paste0("resultados_R/cobertura_anual_", nombre, ".csv"))
  ausentes <- read.csv(paste0("resultados_R/sitio_anio_ausentes_", nombre, ".csv"))
  stopifnot(identical(anual$anio, 1981:2003), all(anual$curvas > 0),
            sum(anual$curvas) == nrow(o$coef),
            sum(anual$sitios_sin_curva) == nrow(ausentes),
            sum(anual$valores_faltantes) == if (fuente == "estaciones") 2871 else 0,
            nrow(ausentes) == if (fuente == "estaciones") 1048 else 0)
  oni <- read.csv(paste0("resultados_R/media_ONI_anual_", nombre, ".csv"))
  stopifnot(sum(oni$n_curvas[oni$dia == 1]) == nrow(o$coef))
  anotar(paste("Periodo cobertura completa y grupos ONI", fuente), 0, 0)
  verif_lib <- read.csv(paste0("resultados_R/verificacion_media_libreria_", nombre, ".csv"))
  stopifnot(verif_lib$estado == "OK")
  F <- o$coef %*% t(base_bspline(1:365, o$K))
  resumen <- read.csv(paste0("resultados_R/resumen_", nombre, ".csv"))
  anotar(paste("Media exportada", fuente), max(abs(resumen$media - colMeans(F))), 1e-9)
  anotar(paste("Varianza exportada", fuente), max(abs(resumen$varianza - apply(F, 2, var))), 1e-8)
  esperado <- boxplot_funcional(F)$atip
  tabla <- read.csv(paste0("resultados_R/profundidades_", nombre, ".csv"))
  stopifnot(identical(tabla$atipica, esperado), nrow(tabla) == nrow(F))
  stopifnot(sum(tabla$atipica) == if (fuente == "estaciones") 170 else 0)
  anotar(paste("Atípicas exportadas", fuente), sum(tabla$atipica != esperado), 0)
}
tierra <- leer_matriz("fda_fase2/coef_satelite_tierra.csv.gz", 3L)$id
tierra <- tierra[as.integer(tierra$anio) %in% 1981:2003, , drop = FALSE]
# Los índices de fila de los dos CSV difieren; nombres, tipos, valores y orden de claves sí deben ser idénticos.
stopifnot(identical(as.list(obj$satelite$id[obj$satelite$tierra, ]), as.list(tierra)))
stopifnot(identical(read.csv("resultados_R/ONI_periodo.csv")$anio, 1981:2003))
guardar_csv(do.call(rbind, controles), "verificacion_flujo.csv")
guardar_csv(data.frame(archivo = names(antes), md5 = as.vector(antes)), "integridad_entradas.csv")
cat("OK: flujo R completo,", length(controles), "controles de resultados y entradas intactas.\n")
cat("Consultar verificacion_media_libreria_*.csv: un PENDIENTE no se convierte en OK.\n")

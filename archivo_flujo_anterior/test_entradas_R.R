# Procesos R limpios: ejecución directa, Source normal y entradas por fase.
local({
  raiz <- normalizePath(".")
  entradas <- c("proyecto_funcional.R", "03_construccion_funcional.R", "04_eda_funcional.R")
  rutas <- file.path(raiz, entradas)
  resultados <- list.files("resultados_R", pattern = "\\.csv$", full.names = TRUE)
  resultados <- resultados[!grepl("verificacion_bloques", resultados)]
  resultados <- normalizePath(resultados)
  antes <- tools::md5sum(resultados)
  carpeta <- tempfile("entrada_limpia_")
  dir.create(carpeta)
  anterior <- getwd()
  on.exit({setwd(anterior); unlink(carpeta, recursive = TRUE)}, add = TRUE)
  setwd(carpeta)
  correr <- function(args, nombre) {
    log <- file.path(raiz, "resultados_R", paste0("prueba_", nombre, ".log"))
    estado <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", args), stdout = log, stderr = log)
    stopifnot(estado == 0L, !dir.exists("resultados_R"), identical(antes, tools::md5sum(resultados)))
    cat("OK:", nombre, "en proceso limpio y desde otra carpeta; CSV idénticos\n")
  }
  correr(shQuote(rutas[1]), "rscript_directo")
  codigo <- paste0("source(", encodeString(rutas[1], quote = '"'), ");",
    "stopifnot(exists('e'),exists('s'),exists('diagnosticos'),e$K==6L,s$K==11L,",
    "eda_est$diagnostico$atipicas==170L,eda_sat$diagnostico$atipicas==0L,is.null(dev.list()));",
    "cat('OK: Source ejecutó y dejó los objetos visibles en la sesión\\n')")
  correr(c("-e", shQuote(codigo)), "source_normal")
  correr(shQuote(rutas[2]), "entrada_03")
  correr(shQuote(rutas[3]), "entrada_04")
})

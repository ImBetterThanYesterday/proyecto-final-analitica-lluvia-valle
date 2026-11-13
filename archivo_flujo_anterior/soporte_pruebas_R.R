# Importación exclusiva del arnés: el alumno ejecuta el script normal completo.
cargar_funciones_prueba <- function(entorno) {
  lineas <- readLines("proyecto_funcional.R", warn = FALSE)
  fin <- match("# FIN DE PREPARACIÓN ----", lineas)
  stopifnot(!is.na(fin))
  conexion <- textConnection(lineas[seq_len(fin)])
  on.exit(close(conexion))
  source(conexion, local = entorno)
  invisible(entorno)
}

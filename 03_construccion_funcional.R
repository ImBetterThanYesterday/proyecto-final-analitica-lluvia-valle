# Entrada opcional: ejecuta solo preparación y construcción del script único.
local({
  archivo <- "/Users/mariavalencia/Downloads/lluvia-valle-git/proyecto_funcional.R"
  lineas <- readLines(archivo, warn = FALSE)
  fin <- match("# EDA FUNCIONAL ----", lineas) - 1L
  stopifnot(!is.na(fin))
  conexion <- textConnection(lineas[seq_len(fin)])
  on.exit(close(conexion))
  source(conexion, local = parent.env(environment()))
})

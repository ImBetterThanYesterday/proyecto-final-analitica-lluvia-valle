# Entrada opcional: ejecuta preparación y EDA usando la construcción R guardada.
local({
  archivo <- "/Users/mariavalencia/Downloads/proyecto-final-analitica-lluvia-valle-main/proyecto_funcional.R"
  lineas <- readLines(archivo, warn = FALSE)
  fin <- match("# FIN DE PREPARACIÓN ----", lineas)
  inicio <- match("# EDA FUNCIONAL ----", lineas)
  stopifnot(!anyNA(c(fin, inicio)), fin < inicio)
  conexion <- textConnection(c(lineas[seq_len(fin)], lineas[inicio:length(lineas)]))
  on.exit(close(conexion))
  source(conexion, local = parent.env(environment()))
})

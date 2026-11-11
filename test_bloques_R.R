# Reproduce selecciones normales consecutivas del editor, desde la primera línea.
# Un source por sección de texto; no usa API especial del script ni precarga sesión.
local({
  lineas <- readLines("proyecto_funcional.R", warn = FALSE)
  inicios <- c(1L, grep("^# [0-9]+\\. ", lineas))
  finales <- c(inicios[-1L] - 1L, length(lineas))
  sesion <- new.env(parent = globalenv())
  registro <- vector("list", length(inicios))
  for (i in seq_along(inicios)) {
    inicio <- Sys.time()
    conexion <- textConnection(lineas[inicios[i]:finales[i]])
    texto <- tryCatch(capture.output(source(conexion, local = sesion)), finally = close(conexion))
    registro[[i]] <- data.frame(bloque = i - 1L, lineas_consola = length(texto),
      segundos = as.numeric(difftime(Sys.time(), inicio, units = "secs")),
      objetos_visibles = length(ls(sesion)), salida = paste(texto, collapse = " | "))
    if (i > 1L) stopifnot(length(texto) > 0L)
    if (i == 4L) stopifnot(exists("seleccion_est", sesion), !exists("e", sesion, inherits = FALSE))
    if (i == 8L) stopifnot(exists("mu", sesion$eda_est), !exists("pr", sesion$eda_est, inherits = FALSE))
  }
  write.csv(do.call(rbind, registro), file.path(sesion$SALIDA, "verificacion_bloques.csv"), row.names = FALSE)
  stopifnot(sesion$e$K == 6L, sesion$s$K == 11L,
            sesion$eda_est$diagnostico$atipicas == 170L,
            sesion$eda_sat$diagnostico$atipicas == 0L)
  imagenes <- list.files(sesion$SALIDA, pattern = "\\.png$", full.names = TRUE)
  stopifnot(length(imagenes) >= 25L, all(file.info(imagenes)$size > 5000), is.null(dev.list()))
  cat("OK: script secuencial desde primera línea, salidas inmediatas, selección anterior al ajuste, momentos anteriores a profundidad, PNG y resultados esperados\n")
})

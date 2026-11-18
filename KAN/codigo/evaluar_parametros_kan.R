#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
raiz_default <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
raiz <- if (length(args)) normalizePath(args[1], mustWork = TRUE) else raiz_default
entrada <- file.path(raiz, "salidas_R_funcional", "para_python")
salida <- file.path(
  raiz, "KAN", "salidas", "analisis_parametros_kan"
)

pares <- read.csv(
  file.path(entrada, "pares.csv"),
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)
coef_obs <- as.matrix(read.csv(
  file.path(entrada, "coef_estacion_mm.csv"),
  check.names = FALSE
))
if (nrow(pares) != nrow(coef_obs)) stop("Pares y coeficientes observados no concuerdan")
filas_validacion <- which(pares$particion == "validación")
if (length(filas_validacion) != 672L) stop("Se esperaban exactamente 672 filas de validación")

orden <- 4L
K <- ncol(coef_obs)
interiores <- seq(1, 365, length.out = K - orden + 2L)
tau <- c(rep(1, orden), interiores[-c(1L, length(interiores))], rep(365, orden))
base_en_nodos <- function(t, tau, orden) {
  B <- outer(t, seq_len(length(tau) - 1L), function(x, j) (x >= tau[j]) & (x < tau[j + 1L])) * 1
  B[t == max(tau), max(which(diff(tau) > 0))] <- 1
  if (orden > 1L) for (m in 2:orden) {
    S <- matrix(0, length(t), length(tau) - m)
    for (j in seq_len(ncol(S))) {
      d1 <- tau[j + m - 1L] - tau[j]
      d2 <- tau[j + m] - tau[j + 1L]
      if (d1 > 0) S[, j] <- S[, j] + (t - tau[j]) / d1 * B[, j]
      if (d2 > 0) S[, j] <- S[, j] + (tau[j + m] - t) / d2 * B[, j + 1L]
    }
    B <- S
  }
  B
}
Phi <- base_en_nodos(1:365, tau, orden)
observada <- pmax(coef_obs[filas_validacion, , drop = FALSE] %*% t(Phi), 0)

ruta_pred <- file.path(salida, "predicciones_validacion_raiz.csv")
if (!file.exists(ruta_pred)) stop(paste("Falta el archivo de predicciones:", ruta_pred))
pred <- read.csv(ruta_pred, stringsAsFactors = FALSE, check.names = FALSE)
nombres_c <- paste0("c", seq_len(K))
requeridas <- c(
  "fila", "codigo", "anio", "particion", "parametro", "valor", "semilla", nombres_c
)
if (!all(requeridas %in% names(pred))) stop("Faltan columnas requeridas en las predicciones")

pred$fila <- as.integer(pred$fila)
if (any(pred$fila < 1L | pred$fila > nrow(pares))) stop("Hay filas fuera del rango de pares.csv")
if (any(pares$particion[pred$fila] != "validación")) {
  stop("Se detectó una predicción fuera de validación")
}
if (any(pred$particion != "validación")) stop("La partición declarada no es validación")
if (any(as.character(pred$codigo) != as.character(pares$codigo[pred$fila])) ||
    any(pred$anio != pares$anio[pred$fila])) {
  stop("Los identificadores fila/código/año no concuerdan con pares.csv")
}

grupos <- split(
  pred,
  interaction(pred$parametro, pred$valor, pred$semilla, drop = TRUE, lex.order = TRUE)
)
por_semilla <- do.call(rbind, lapply(grupos, function(g) {
  if (nrow(g) != 672L || length(unique(g$fila)) != 672L) {
    stop("Cada escenario y semilla debe contener las 672 filas de validación")
  }
  raiz_pred <- as.matrix(g[, nombres_c, drop = FALSE]) %*% t(Phi)
  predicha <- pmax(raiz_pred, 0)^2
  fila_validacion <- match(g$fila, filas_validacion)
  if (anyNA(fila_validacion)) stop("Se encontró una fila que no pertenece a validación")
  obs <- observada[fila_validacion, , drop = FALSE]
  l2_curva <- sqrt(rowSums((predicha - obs)^2) / rowSums(obs^2))
  por_estacion <- aggregate(l2_curva, list(codigo = as.character(pares$codigo[g$fila])), mean)
  data.frame(
    parametro = g$parametro[1],
    valor = g$valor[1],
    semilla = g$semilla[1],
    curvas = nrow(g),
    estaciones = nrow(por_estacion),
    L2_relativo_media_estacion = mean(por_estacion$x),
    stringsAsFactors = FALSE
  )
}))

resumen <- do.call(rbind, lapply(
  split(por_semilla, interaction(por_semilla$parametro, por_semilla$valor, drop = TRUE)),
  function(g) data.frame(
    parametro = g$parametro[1],
    valor = g$valor[1],
    semillas = nrow(g),
    L2_relativo_media = mean(g$L2_relativo_media_estacion),
    L2_relativo_DE = sd(g$L2_relativo_media_estacion),
    stringsAsFactors = FALSE
  )
))
mejores <- do.call(rbind, lapply(split(resumen, resumen$parametro), function(g) {
  g[which.min(g$L2_relativo_media), , drop = FALSE]
}))

dir.create(salida, recursive = TRUE, showWarnings = FALSE)
write.csv(por_semilla, file.path(salida, "L2_relativo_validacion_por_semilla_R.csv"), row.names = FALSE)
write.csv(resumen, file.path(salida, "L2_relativo_validacion_resumen_R.csv"), row.names = FALSE)
write.csv(mejores, file.path(salida, "mejor_valor_por_parametro_R.csv"), row.names = FALSE)
print(mejores, row.names = FALSE)
cat("\nSolo se evaluaron las 672 filas de validación. No se usaron predicciones de prueba.\n")

# Estadística de diferencias pareadas funcionales con remuestreo por estación.
# Se usa desde 04_eda_funcional.R y se prueba con datos sintéticos.
bootstrap_cluster_difference <- function(D, station_id, B = 1000L, seed = 20261111L) {
  if (!is.matrix(D) || nrow(D) != length(station_id) || nrow(D) < 2L || ncol(D) < 2L)
    stop("D debe ser una matriz pares×tiempo con una estación por fila")
  if (any(!is.finite(D))) stop("La diferencia funcional pareada contiene valores no finitos")
  station_id <- as.character(station_id)
  ids <- unique(station_id)
  if (length(ids) < 2L) stop("Se necesitan al menos dos estaciones para bootstrap por conglomerado")
  por_estacion <- t(vapply(ids, function(id) colMeans(D[station_id == id, , drop = FALSE]),
                           numeric(ncol(D))))
  media <- colMeans(por_estacion)
  ee <- apply(por_estacion, 2, stats::sd) / sqrt(nrow(por_estacion))
  set.seed(seed)
  boot <- matrix(NA_real_, B, ncol(D))
  for (b in seq_len(B)) {
    muestra <- sample.int(nrow(por_estacion), nrow(por_estacion), replace = TRUE)
    boot[b, ] <- colMeans(por_estacion[muestra, , drop = FALSE])
  }
  qs <- apply(boot, 2, stats::quantile, probs = c(.025, .975), names = FALSE)
  inferior <- qs[1, ]; superior <- qs[2, ]
  list(ids = ids, por_estacion = por_estacion, media = media, ee = ee,
       inferior_95 = inferior, superior_95 = superior,
       zero_excluded = inferior > 0 | superior < 0, bootstrap = boot)
}

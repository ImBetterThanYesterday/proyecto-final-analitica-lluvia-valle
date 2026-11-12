# Cross-language conformance check. Run after the Python notebook 03 and R script 03.
local({
  dir <- "fda_fase1"
  required <- file.path(dir, c(
    "coef_estaciones.csv", "coef_satelite_python.csv.gz", "resumen_seleccion_K.csv",
    "coef_estaciones_R.csv", "coef_satelite_R.csv.gz", "gram_estaciones.csv", "gram_satelite.csv",
    "gram_estaciones_R.csv", "gram_satelite_R.csv"
  ))
  missing <- required[!file.exists(required)]
  if (length(missing)) stop("Ejecute Python 03 y R 03 antes de comparar; faltan: ", paste(missing, collapse = ", "))

  py_est <- read.csv(file.path(dir, "coef_estaciones.csv"), check.names = FALSE)
  r_est <- read.csv(file.path(dir, "coef_estaciones_R.csv"), check.names = FALSE)
  py_sat <- read.csv(gzfile(file.path(dir, "coef_satelite_python.csv.gz")), check.names = FALSE)
  r_sat <- read.csv(gzfile(file.path(dir, "coef_satelite_R.csv.gz")), check.names = FALSE)
  py_k <- read.csv(file.path(dir, "resumen_seleccion_K.csv"), check.names = FALSE)

  compare_table <- function(r, py, n_id, label) {
    if (nrow(r) != nrow(py)) stop(label, ": distinta cantidad de filas R/Python")
    if (!isTRUE(all.equal(r[, seq_len(n_id), drop = FALSE], py[, seq_len(n_id), drop = FALSE],
                          check.attributes = FALSE))) stop(label, ": las claves no coinciden en orden/valor")
    cr <- as.matrix(r[, -seq_len(n_id), drop = FALSE]); cp <- as.matrix(py[, -seq_len(n_id), drop = FALSE])
    if (!identical(dim(cr), dim(cp))) stop(label, ": dimensiones de coeficientes distintas")
    max(abs(cr - cp))
  }
  error_est <- compare_table(r_est, py_est, 2L, "estaciones")
  error_sat <- compare_table(r_sat, py_sat, 3L, "satélite")

  compare_gram <- function(r_path, py_path, label) {
    r <- as.matrix(read.csv(r_path, header = FALSE, check.names = FALSE))
    py <- as.matrix(read.csv(py_path, header = FALSE, check.names = FALSE))
    if (!identical(dim(r), dim(py))) stop(label, ": dimensiones W distintas")
    max(abs(r - py))
  }
  error_w_est <- compare_gram(file.path(dir, "gram_estaciones_R.csv"), file.path(dir, "gram_estaciones.csv"), "W estaciones")
  error_w_sat <- compare_gram(file.path(dir, "gram_satelite_R.csv"), file.path(dir, "gram_satelite.csv"), "W satélite")

  k_vcg <- as.integer(py_k[["K por VCG"]][1:2])
  k_r <- as.integer(c(ncol(r_est) - 2L, ncol(r_sat) - 3L))
  # El notebook ajusta el satélite con los valores float32 del TIFF; R, con curvas_satelite_mm.csv.gz
  # (5 decimales). Cota exacta de ese redondeo: max_k sum_j |y2cMap_kj| * 5e-6. La igualdad con la misma
  # entrada (< 1e-8) se prueba en tests/validar_salidas_cruzadas.py.
  curvas <- read.csv(gzfile(file.path(dir, "curvas_satelite_mm.csv.gz")), check.names = FALSE, nrows = 1)
  t_sat <- as.numeric(sub("^t", "", names(curvas)[-(1:3)]))
  K_sat <- ncol(r_sat) - 3L; orden <- 4L; rango <- c(1, 365)
  interiores <- seq(rango[1], rango[2], length.out = K_sat - orden + 2L)
  tau <- c(rep(rango[1], orden), interiores[-c(1L, length(interiores))], rep(rango[2], orden))
  B <- outer(t_sat, seq_len(length(tau) - 1L), function(x, j) (x >= tau[j]) & (x < tau[j + 1L])) * 1
  B[t_sat == max(tau), max(which(diff(tau) > 0))] <- 1
  for (m in 2:orden) {
    S <- matrix(0, length(t_sat), length(tau) - m)
    for (j in seq_len(ncol(S))) {
      d1 <- tau[j + m - 1L] - tau[j]; d2 <- tau[j + m] - tau[j + 1L]
      if (d1 > 0) S[, j] <- S[, j] + (t_sat - tau[j]) / d1 * B[, j]
      if (d2 > 0) S[, j] <- S[, j] + (tau[j + m] - t_sat) / d2 * B[, j + 1L]
    }
    B <- S
  }
  mapa <- solve(crossprod(B), t(B))
  cota_sat <- max(rowSums(abs(mapa))) * 5e-6
  stopifnot(identical(k_r, k_vcg), error_est < 1e-8, error_sat <= cota_sat,
            error_w_est < 1e-8, error_w_sat < 1e-8)
  results <- data.frame(
    control = c("filas_estaciones", "filas_satelite", "error_coef_estaciones", "error_coef_satelite",
                "cota_redondeo_satelite", "K_VCG_estaciones", "K_VCG_satelite", "error_W_estaciones", "error_W_satelite"),
    valor = c(nrow(r_est), nrow(r_sat), error_est, error_sat, cota_sat, k_r, error_w_est, error_w_sat)
  )
  dir.create(file.path("fda_fase1", "tablas_R"), showWarnings = FALSE, recursive = TRUE)
  write.csv(results, file.path("fda_fase1", "tablas_R", "validacion_cruzada_R_Python.csv"), row.names = FALSE)
  print(results, row.names = FALSE)
  cat("OK: claves, K, coeficientes de ambas fuentes y matrices W coinciden dentro de tolerancia.\n")
})

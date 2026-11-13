root <- normalizePath(".")
source(file.path(root, "scripts", "analizar_pares_estacion_pixel.R"), local = TRUE)

D <- rbind(
  c(1, 2, 3),
  c(3, 4, 5),
  c(7, 8, 9),
  c(9, 10, 11)
)
station_id <- c("A", "A", "B", "B")
cluster_means <- rbind(colMeans(D[station_id == "A", , drop = FALSE]),
                       colMeans(D[station_id == "B", , drop = FALSE]))
a <- bootstrap_cluster_difference(D, station_id, B = 50L, seed = 7L)
stopifnot(nrow(a$por_estacion) == 2L,
          max(abs(a$media - colMeans(cluster_means))) < 1e-12,
          nrow(a$bootstrap) == 50L,
          all(a$inferior_95 <= a$superior_95),
          identical(a$zero_excluded, a$inferior_95 > 0 | a$superior_95 < 0))
cat("OK: media pareada por estación y bootstrap por conglomerados reproducible\n")

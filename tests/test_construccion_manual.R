# Verificaciones independientes del constructor manual. Los paquetes se usan
# exclusivamente aquí como referencia de prueba, nunca para producir las curvas.
source("03_construccion_funcional.R", local = .GlobalEnv)

for (K in c(6L, 8L, 11L, 17L, 40L)) {
  Phi_manual <- base_bspline(t_est, K)
  Phi_ref <- splines::splineDesign(nodos(K), t_est, ord = ORDEN, outer.ok = TRUE)
  stopifnot(max(abs(Phi_manual - Phi_ref)) < 1e-12,
            min(Phi_manual) >= -1e-12,
            max(abs(rowSums(Phi_manual) - 1)) < 1e-10)
}

stopifnot(all(anio_de(est$id[ent_est, , drop = FALSE]) <= ENTREN_FIN),
          all(anio_de(sat$id[ent_sat, , drop = FALSE]) <= ENTREN_FIN))

for (d in list(e, s)) {
  stopifnot(max(abs(d$H - t(d$H))) < 1e-10,
            max(abs(d$H %*% d$H - d$H)) < 1e-10,
            abs(sum(diag(d$H)) - d$K) < 1e-8)
}

i <- which(completas_est)[1]
qr_coef <- qr.coef(qr(base_bspline(t_est, K_EST)), est$Y[i, ])
stopifnot(max(abs(e$coef[i, ] - qr_coef)) < 1e-8)

cat("OK: base manual vs splineDesign, entrenamiento-only K, H/S y MCO vs QR.\n")

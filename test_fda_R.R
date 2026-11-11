# Pruebas matemáticas independientes; no requiere paquetes de pruebas.
entorno <- new.env()
source("soporte_pruebas_R.R")
cargar_funciones_prueba(entorno)
Y <- matrix(c(0, 0, 1, 2), ncol = 1)
pares <- combn(nrow(Y), 2)
esperada <- sapply(seq_len(nrow(Y)), function(i) mean(apply(pares, 2, function(p)
  mean(pmin(Y[p[1], ], Y[p[2], ]) <= Y[i, ] & Y[i, ] <= pmax(Y[p[1], ], Y[p[2], ])))))
stopifnot(max(abs(entorno$profundidades(Y)$MBD - esperada)) < 1e-12)
stopifnot(all(entorno$profundidades(matrix(0, 4, 5))$MBD == 1))
rechaza <- function(expr) inherits(tryCatch({force(expr); NULL}, error = identity), "error")
stopifnot(rechaza(entorno$profundidades(matrix(0, 1, 5))))
stopifnot(rechaza(entorno$profundidades(matrix(c(NA, 1, 2, 3), 2))))
cat("OK: MBD inclusiva por pares, empates y rejilla de un punto\n")
# depth.FM de clase usa FM1 con scale=TRUE; el PDF muestra la escala sin transformar.
V <- rbind(c(0, 0, 2), c(0, 1, 1), c(1, 2, 0), c(3, 3, 3))
F <- sapply(seq_len(ncol(V)), function(j) ecdf(V[, j])(V[, j]))
fm_pdf <- rowMeans(1 - abs(.5 - F))
stopifnot(max(abs(entorno$profundidades(V)$FM - (2 * fm_pdf - 1))) < 1e-12)
stopifnot(max(abs(entorno$profundidades(V)$FM_pdf - fm_pdf)) < 1e-12)
# El recorte por cuantil conserva todos los empates del punto de corte, como depth.FM.
pr <- c(.1, .2, .2, .2, .4, .5, .6, .7)
stopifnot(identical(entorno$centralidad_FM(matrix(1:24, 8), pr)$incluidas, 2:8))
# Igualdad exacta con un límite no degenerado: fbplot usa <= / >=.
W <- matrix(rep(c(-1.5, 0, 1, 2.5), 3), nrow = 4)
bx <- entorno$boxplot_funcional(W, c(0, 1, .9, 0))
stopifnot(identical(bx$atip, c(TRUE, FALSE, FALSE, TRUE)),
          all(bx$bigote_inferior == 0), all(bx$bigote_superior == 1))
cat("OK: FM contra ECDF, escalas PDF/script, recorte con empates y límites/bigotes del boxplot\n")
with(entorno, {
  tt <- seq(1, 365, length.out = 41)
  P <- base_bspline(tt, 6)
  stopifnot(max(abs(rowSums(P) - 1)) < 1e-12, min(P) >= 0)
  # Un polinomio cúbico está en el espacio de B-splines cúbicas.
  x <- (tt - 1) / 364; y <- 2 + 3 * x - x^2 + x^3
  a <- ajuste_mco(y, P)
  stopifnot(max(abs(P %*% a$coef - y)) < 1e-11, all(control_H(a$H, 6) < 1e-11))
  derivada <- as.vector(deriv_base(tt, 6) %*% a$coef)
  stopifnot(max(abs(derivada - (3 - 2 * x + 3 * x^2) / 364)) < 1e-11)
  yy <- y; yy[c(10, 11, 20)] <- NA
  b <- ajuste_mco(yy, P)
  stopifnot(max(abs(P %*% b$coef - y)) < 1e-11,
            all(control_H(b$H, 6) < 1e-11), sum(b$observado) == 38)
  stopifnot(rechaza(ajuste_mco(rep(NA_real_, length(tt)), P)),
            rechaza(ajuste_mco(y, cbind(P, P[, 1]))))
  set.seed(2026)
  Z <- matrix(rnorm(5 * length(tt)), 5)
  tabla <- criterios(Z, tt, 4:8)
  for (k in 4:8) {
    B <- base_bspline(tt, k); H <- B %*% solve(crossprod(B), t(B))
    E <- Z - Z %*% t(H); n <- length(tt)
    esperada <- c(mean(E^2), mean(sweep(E, 2, 1 - diag(H), "/")^2), n * sum(E^2) / (n - k)^2)
    stopifnot(max(abs(as.numeric(tabla[tabla$K == k, -1]) - esperada)) < 1e-10)
  }
  # VC verdaderamente independiente: ajustar ocultando cada punto.
  B <- base_bspline(tt, 6); v <- Z[1, ]; a <- ajuste_mco(v, B)
  loo <- sapply(seq_along(v), function(j) B[j, ] %*% qr.solve(B[-j, ], v[-j]))
  atajo <- mean(((v - B %*% a$coef) / (1 - diag(a$H)))^2)
  stopifnot(abs(atajo - mean((v - loo)^2)) < 1e-10)
  stopifnot(max(abs(a$ee^2 - a$sigma2 * diag(a$H %*% t(a$H)))) < 1e-10)
})
cat("OK: base, cúbico conocido, derivada analítica incluyendo bordes, MCO con huecos, rango, H, ECM/VC/VCG y bandas\n")
.libPaths(c(".Rlib", .libPaths()))
if (!requireNamespace("fda", quietly = TRUE)) stop("Falta fda para los contrastes docentes de EDA")
with(entorno, {
  # Contrastes de EDA independientes de los ejemplos polinómicos anteriores.
  set.seed(731)
  C <- matrix(rnorm(8 * 6), 8)
  tt <- seq(1, 365, length.out = 31)
  P <- base_bspline(tt, 6); Y <- C %*% t(P)
  fd <- fda::fd(t(C), fda::create.bspline.basis(RANGO, 6, ORDEN))
  Sc <- cov(C)
  stopifnot(max(abs(P %*% Sc %*% t(P) - fda::eval.bifd(tt, tt, fda::var.fd(fd)))) < 1e-10,
            max(abs(deriv_base(tt, 6) %*% colMeans(C) -
                      fda::eval.fd(tt, fda::mean.fd(fd), Lfdobj = 1))) < 1e-10)
  # deriv.fd reproyecta sobre la base cúbica original: no es el contraste exacto.
  error_reproyeccion <- max(abs(fda::eval.fd(tt, fda::deriv.fd(fda::mean.fd(fd), 1)) -
                                fda::eval.fd(tt, fda::mean.fd(fd), Lfdobj = 1)))
  cat("Diagnóstico docente: deriv.fd reproyectada vs derivada exacta =", error_reproyeccion, "\n")
  ref <- fda::fbplot(t(Y), tt, plot = FALSE)
  pr <- profundidades(Y)
  stopifnot(max(abs(pr$MBD - ref$depth)) < 1e-12,
            identical(as.integer(which(boxplot_funcional(Y)$atip)), as.integer(ref$outpoint)))
  ref <- fda::fbplot(t(W), depth = c(0, 1, .9, 0), plot = FALSE)
  stopifnot(identical(as.integer(which(boxplot_funcional(W, c(0, 1, .9, 0))$atip)),
                      as.integer(ref$outpoint)))
  # Diferencia conocida: rangos medios de fda frente a definición inclusiva con empates.
  V <- matrix(rep(c(0, 0, 1, 2), 3), 4)
  ref <- fda::fbplot(t(V), plot = FALSE)
  stopifnot(max(abs(profundidades(V)$MBD - ref$depth)) > .01)
})
cat("OK: covarianza y derivada contra fda; MBD sin empates y boxplot contra fbplot; discrepancia de empates documentada\n")

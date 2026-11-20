# =============================================================================
# REGRESIÓN FUNCIONAL CONCURRENTE (SATÉLITE -> ESTACIÓN), ESTILO DE CLASE
# Proyecto final de Analítica de Datos (UAO 2026-2S) — lluvia en el Valle del Cauca 1981-2003
# Rama: maria-regresion-lineal
# =============================================================================
# Modelo (ejercicio 3 de Script_Regresion_Funcional.R de clase):
#
#     y_i(t) = b0(t) + b1(t) x_i(t) + e_i(t)                      (modelo A, el de clase)
#     y_i(t) = b0(t) + b1(t) x_i(t) + g(t) altitud_i + d(t) ONI_i + e_i(t)   (modelo B)
#
# y_i(t): curva anual de la estación; x_i(t): curva anual del satélite en el píxel de la
# estación; las dos en escala raíz (raíz de mm/día), igual que la regresión del equipo.
# A diferencia de la regresión función a función de la sección 5 de lluvia_valle.R
# (β(s,t): todo el año del satélite predice cada día de la estación), aquí cada día t de
# la estación se explica SOLO con el día t del satélite. Eso permite leer b1(t) y R²(t):
# en qué época del año el satélite sirve más y en cuál menos.
#
# Reglas del proyecto que se respetan:
#   - Partición 1981-1995 (entrenamiento) / 1996-2000 (validación) / 2001-2003 (prueba).
#   - Los parámetros (K de las funciones b, lambda) se eligen SOLO con validación.
#   - Todos los años entran (D-47). Reconstrucción final en mm no negativa: max(ŷ, 0)².
#   - Las mismas curvas y las mismas métricas que 5.5 de lluvia_valle.R.
#   - Paso a paso en R base; fda::fRegress solo como contraste (herramienta de clase).
#
# Entradas (las deja lluvia_valle.R): salidas_R_funcional/para_python/ y tablas/
# Salidas: salidas_R_concurrente/ (tablas y figuras). No modifica nada de lo existente.
# Uso, desde la raíz del repositorio:   LC_ALL=en_US.UTF-8 Rscript regresion_concurrente.R
# =============================================================================

for (loc in c("en_US.UTF-8", "es_ES.UTF-8", "C.UTF-8")) if (nzchar(suppressWarnings(Sys.setlocale("LC_ALL", loc)))) break   # R necesita UTF-8 para leer los acentos
set.seed(20261114)

# ---- 0. CONFIGURACIÓN Y UTILIDADES ------------------------------------------
# Conjunto de datos (D-61, D-63): principal (≥10 años) o r6_omm (sensibilidad de R6), con la variable CONJUNTO
CONJUNTO_C <- Sys.getenv("CONJUNTO", "principal")
SUFIJO_C <- if (CONJUNTO_C == "principal") "" else paste0("_", CONJUNTO_C)
DIR_IN  <- file.path(paste0("salidas_R_funcional", SUFIJO_C), "para_python")
DIR_OUT <- paste0("salidas_R_concurrente", SUFIJO_C)
dir.create(file.path(DIR_OUT, "tablas"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DIR_OUT, "figuras"), recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(file.path(DIR_IN, "pares.csv")))

RANGO <- c(1, 365); ORDEN <- 4L; T_DIA <- 1:365
N_BOOT <- 2000L
VENT_PICO_1 <- 60:200; VENT_PICO_2 <- 230:350      # mismas ventanas de los picos que 5.5
MESES <- data.frame(mes = c("ene", "feb", "mar", "abr", "may", "jun", "jul", "ago", "sep", "oct", "nov", "dic"),
                    ini = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335))

titulo <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
tabla_f <- function(df, nombre) {
  write.csv(df, file.path(DIR_OUT, "tablas", paste0(nombre, ".csv")), row.names = FALSE, fileEncoding = "UTF-8")
  invisible(df)
}
figura <- function(nombre, f, ancho = 11, alto = 5) {
  png(file.path(DIR_OUT, "figuras", paste0(nombre, ".png")), width = ancho, height = alto, units = "in", res = 150)
  on.exit(dev.off()); f()
}

# Base B-spline cúbica de K funciones con nodos equiespaciados en [1, 365] (como la del equipo)
nodos <- function(K, orden = ORDEN) c(rep(RANGO[1], orden), seq(RANGO[1], RANGO[2], length.out = K - orden + 2)[-c(1, K - orden + 2)], rep(RANGO[2], orden))
base_bspline <- function(t, K, orden = ORDEN) splines::splineDesign(nodos(K, orden), t, ord = orden, outer.ok = TRUE)
# Rugosidad R = ∫ φ''(t) φ''(t)ᵀ dt, por trapecios sobre una malla fina (φ'' es lineal a trozos)
gram_rugosidad <- function(K) {
  tt <- seq(RANGO[1], RANGO[2], by = 0.05)
  D2 <- splines::splineDesign(nodos(K), tt, ord = ORDEN, derivs = rep(2, length(tt)), outer.ok = TRUE)
  w <- rep(0.05, length(tt)); w[c(1, length(tt))] <- 0.025
  crossprod(D2, D2 * w)
}
W_T <- rep(1, 365); W_T[c(1, 365)] <- 0.5            # pesos del trapecio sobre los 365 días: ∫ ≈ Σ w_t f(t)

fecha_dia <- function(d) format(as.Date(d - 1, origin = "2001-01-01"), "%d-%b")

# ---- 1. DATOS: LOS MISMOS PARES Y CURVAS DEL EQUIPO --------------------------
titulo("1. Datos")
pares <- read.csv(file.path(DIR_IN, "pares.csv"), encoding = "UTF-8", stringsAsFactors = FALSE)
C_e_r <- as.matrix(read.csv(file.path(DIR_IN, "coef_estacion_raiz.csv")))   # coef. B-spline de la estación, escala raíz
C_s_r <- as.matrix(read.csv(file.path(DIR_IN, "coef_satelite_raiz.csv")))   # idem satélite
C_e_m <- as.matrix(read.csv(file.path(DIR_IN, "coef_estacion_mm.csv")))     # escala mm/día (para evaluar)
C_s_m <- as.matrix(read.csv(file.path(DIR_IN, "coef_satelite_mm.csv")))
stopifnot(nrow(pares) == nrow(C_e_r), nrow(pares) == nrow(C_s_r), nrow(C_e_m) == nrow(pares), nrow(C_s_m) == nrow(pares))
K_e <- ncol(C_e_r); K_s <- ncol(C_s_r)
Phi_e <- base_bspline(T_DIA, K_e); Phi_s <- base_bspline(T_DIA, K_s)

Y_r <- C_e_r %*% t(Phi_e)             # n × 365, escala raíz
X_r <- C_s_r %*% t(Phi_s)
Y_mm <- pmax(C_e_m %*% t(Phi_e), 0)   # igual que 5.1 de lluvia_valle.R
X_mm <- pmax(C_s_m %*% t(Phi_s), 0)

ent <- pares$usar & pares$particion == "entrenamiento"
val <- pares$usar & grepl("^valid", pares$particion)
pru <- pares$usar & pares$particion == "prueba"
stopifnot(sum(ent) + sum(val) + sum(pru) == sum(pares$usar), all(pares$anio[ent] <= 1995), all(pares$anio[pru] >= 2001))   # partición por años
cat(sprintf("Curvas: entrenamiento %d, validación %d, prueba %d | estaciones %d | K estación = %d, K satélite = %d\n",
            sum(ent), sum(val), sum(pru), length(unique(pares$codigo)), K_e, K_s))

alt_m <- mean(pares$altitud[ent]); alt_s <- sd(pares$altitud[ent])       # estandarización solo con entrenamiento
ALT <- (pares$altitud - alt_m) / alt_s; ONI <- pares$ONI_anual
cat(sprintf("Altitud de entrenamiento: media %.1f, DE %.1f (la del equipo: 1394.2 y 560.0)\n", alt_m, alt_s))

# ---- 2. ESTIMACIÓN PASO A PASO (MÍNIMOS CUADRADOS PENALIZADOS) ---------------
# Cada coeficiente funcional b_a(t) = Σ_k θ_ak φ_k(t) en una base B-spline de K_b funciones.
# Se minimiza  Σ_i ∫ (y_i - Σ_a c_ia(t) b_a(t))² dt + lambda Σ_a ∫ b_a''(t)² dt,
# con c_i0 = 1, c_i1(t) = x_i(t) y, en el modelo B, c_i2 = altitud_i y c_i3 = ONI_i.
# Ecuaciones normales por bloques: G_ab = Φᵀ diag(w · Σ_i c_ia c_ib) Φ ;  r_a = Φᵀ (w · Σ_i c_ia y_i).
covariables <- function(i, modelo) {          # lista de matrices n × 365 (una por coeficiente)
  n <- length(i); L <- list(const = matrix(1, n, 365), sat = X_r[i, , drop = FALSE])
  if (modelo == "B") { L$altitud <- matrix(ALT[i], n, 365); L$ONI <- matrix(ONI[i], n, 365) }
  L
}
preparar <- function(i, modelo, K_b) {        # todo lo que no depende de lambda
  Cv <- covariables(i, modelo); p <- length(Cv); Phi_b <- base_bspline(T_DIA, K_b)
  s_ab <- function(a, b) colSums(Cv[[a]] * Cv[[b]])
  G <- matrix(0, p * K_b, p * K_b); r <- numeric(p * K_b)
  for (a in 1:p) {
    r[(a - 1) * K_b + 1:K_b] <- crossprod(Phi_b, W_T * colSums(Cv[[a]] * Y_r[i, , drop = FALSE]))
    for (b in a:p) {
      blq <- crossprod(Phi_b, Phi_b * (W_T * s_ab(a, b)))
      G[(a - 1) * K_b + 1:K_b, (b - 1) * K_b + 1:K_b] <- blq
      if (b != a) G[(b - 1) * K_b + 1:K_b, (a - 1) * K_b + 1:K_b] <- t(blq)
    }
  }
  list(G = G, r = r, p = p, K_b = K_b, Phi_b = Phi_b, R = kronecker(diag(p), gram_rugosidad(K_b)), nombres = names(Cv))
}
resolver <- function(pre, lambda) solve(pre$G + lambda * pre$R, pre$r)       # θ̂
coef_t <- function(pre, theta) {                                             # b_a(t): 365 × p
  B <- sapply(1:pre$p, function(a) pre$Phi_b %*% theta[(a - 1) * pre$K_b + 1:pre$K_b]); colnames(B) <- pre$nombres; B
}
predecir_raiz <- function(B, i, modelo) {                                    # ŷ_i(t) en escala raíz
  Cv <- covariables(i, modelo); P <- 0
  for (a in seq_along(Cv)) P <- P + sweep(Cv[[a]], 2, B[, a], "*")
  P
}
a_mm <- function(R) pmax(R, 0)^2
rmse_l2 <- function(P, i) sqrt(rowMeans((P - Y_mm[i, , drop = FALSE])^2))    # RMSE L² en mm/día por curva

# ---- 3. ELECCIÓN DE K_b Y lambda CON VALIDACIÓN (1996-2000) ------------------
titulo("3. Selección de K_b y lambda con validación")
LAMBDAS <- c(0, 10^seq(-2, 12, by = 0.5))
KBS <- c(8L, 12L, 16L)
i_ent <- which(ent); i_val <- which(val); i_pru <- which(pru)
seleccion <- do.call(rbind, lapply(c("A", "B"), function(mod) do.call(rbind, lapply(KBS, function(K_b) {
  pre <- preparar(i_ent, mod, K_b)
  do.call(rbind, lapply(LAMBDAS, function(l) {
    B <- coef_t(pre, resolver(pre, l))
    data.frame(modelo = mod, K_b = K_b, lambda = l,
               rmse_entrenamiento = mean(rmse_l2(a_mm(predecir_raiz(B, i_ent, mod)), i_ent)),
               rmse_validacion    = mean(rmse_l2(a_mm(predecir_raiz(B, i_val, mod)), i_val)))
  }))
}))))
tabla_f(seleccion, "C1_seleccion_Kb_lambda")
# Regla (la misma idea que M en 5.3): el menor K_b cuyo mejor RMSE de validación queda a menos de 0,5 % del mínimo del modelo.
elegir <- function(mod) {
  s <- seleccion[seleccion$modelo == mod, ]
  mejor <- do.call(rbind, lapply(KBS, function(k) { z <- s[s$K_b == k, ]; z[which.min(z$rmse_validacion), ] }))
  K_sel <- min(mejor$K_b[mejor$rmse_validacion <= min(mejor$rmse_validacion) * 1.005])
  list(K_b = K_sel, lambda = mejor$lambda[mejor$K_b == K_sel], mejor = mejor)
}
eleccion <- list(A = elegir("A"), B = elegir("B"))
for (mod in c("A", "B")) {
  cat(sprintf("Modelo %s: K_b = %d, lambda = %s (RMSE validación %.3f mm/día)\n", mod, eleccion[[mod]]$K_b,
              format(eleccion[[mod]]$lambda, scientific = TRUE), with(eleccion[[mod]]$mejor, rmse_validacion[K_b == eleccion[[mod]]$K_b])))
  print(eleccion[[mod]]$mejor, row.names = FALSE)
}
figura("C1_seleccion_lambda", function() {
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  for (mod in c("A", "B")) {
    s <- seleccion[seleccion$modelo == mod & seleccion$lambda > 0, ]
    plot(NA, xlim = range(log10(s$lambda)), ylim = range(s$rmse_validacion), xlab = expression(log[10](lambda)),
         ylab = "RMSE L2 medio en validación (mm/día)", main = paste("Modelo", mod))
    for (j in seq_along(KBS)) lines(log10(s$lambda[s$K_b == KBS[j]]), s$rmse_validacion[s$K_b == KBS[j]], col = j + 1, lwd = 2)
    abline(v = log10(eleccion[[mod]]$lambda), lty = 2); legend("topleft", paste("K_b =", KBS), col = seq_along(KBS) + 1, lwd = 2, bty = "n")
  }
})

# ---- 4. AJUSTE FINAL (SOLO ENTRENAMIENTO) Y COEFICIENTES FUNCIONALES ---------
titulo("4. Coeficientes funcionales b(t) con bandas punto a punto")
ajuste <- lapply(c(A = "A", B = "B"), function(mod) {
  pre <- preparar(i_ent, mod, eleccion[[mod]]$K_b); th <- resolver(pre, eleccion[[mod]]$lambda)
  list(pre = pre, theta = th, B = coef_t(pre, th), lambda = eleccion[[mod]]$lambda)
})
# Bandas ±2 EE punto a punto (estimador sándwich, como el SigmaE de clase): la covarianza de los
# residuos entre días, Σ_E (365 × 365), se estima con los residuos de entrenamiento.
# Cov(θ̂) = H⁻¹ [Σ_i Z_iᵀ W Σ_E W Z_i] H⁻¹, con H = G + lambda R. Son bandas punto a punto, no simultáneas.
bandas <- function(mod) {
  a <- ajuste[[mod]]; pre <- a$pre; Cv <- covariables(i_ent, mod); p <- pre$p; K_b <- pre$K_b
  Res <- Y_r[i_ent, ] - predecir_raiz(a$B, i_ent, mod)
  SigE <- crossprod(sweep(Res, 2, colMeans(Res))) / (length(i_ent) - 1)
  Hinv <- solve(pre$G + a$lambda * pre$R); WSW <- (W_T %o% W_T) * SigE
  M <- matrix(0, p * K_b, p * K_b)
  for (k in seq_along(i_ent)) {
    Z <- do.call(cbind, lapply(1:p, function(a_) Cv[[a_]][k, ] * pre$Phi_b))   # 365 × (p K_b)
    M <- M + crossprod(Z, WSW %*% Z)
  }
  V <- Hinv %*% M %*% Hinv
  sapply(1:p, function(a_) { idx <- (a_ - 1) * K_b + 1:K_b; sqrt(rowSums((pre$Phi_b %*% V[idx, idx]) * pre$Phi_b)) })
}
EE <- lapply(c(A = "A", B = "B"), bandas)
ajuste$A$EE <- EE$A; ajuste$B$EE <- EE$B
coef_tab <- do.call(rbind, lapply(c("A", "B"), function(mod) do.call(rbind, lapply(seq_len(ajuste[[mod]]$pre$p), function(a) {
  data.frame(modelo = mod, coeficiente = ajuste[[mod]]$pre$nombres[a], dia = T_DIA, valor = ajuste[[mod]]$B[, a], EE = ajuste[[mod]]$EE[, a])
}))))
tabla_f(coef_tab, "C2_coeficientes_funcionales")
figura("C2_coeficientes_modeloA", function() {
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  for (a in 1:2) {
    b <- ajuste$A$B[, a]; e <- ajuste$A$EE[, a]
    matplot(T_DIA, cbind(b, b + 2 * e, b - 2 * e), type = "l", lty = c(1, 2, 2), lwd = c(2, 1, 1), col = c("black", "red", "red"),
            xlab = "Día del año", ylab = "", main = c("Intercepto b0(t)", "Efecto del satélite b1(t)")[a]); abline(h = c(0, if (a == 2) 1), lty = 3)
  }
})
figura("C2_coeficientes_modeloB", function() {
  par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))
  for (a in 1:4) {
    b <- ajuste$B$B[, a]; e <- ajuste$B$EE[, a]
    matplot(T_DIA, cbind(b, b + 2 * e, b - 2 * e), type = "l", lty = c(1, 2, 2), lwd = c(2, 1, 1), col = c("black", "red", "red"),
            xlab = "Día del año", ylab = "", main = c("b0(t)", "b1(t): satélite", "g(t): altitud (por DE)", "d(t): ONI (por unidad)")[a]); abline(h = 0, lty = 3)
  }
})

# ---- 4b. BANDAS POR BOOTSTRAP DE ESTACIONES ----------------------------------
titulo("4b. Bandas de b(t) por bootstrap de estaciones")
# Las bandas sándwich de arriba usan Σ_E entre días, pero suman las curvas como si fueran
# independientes. Los años de una misma estación repiten su sesgo: en la comparación por fuente
# (sección 4.1 de lluvia_valle.R) el EE por estación es 3,7 veces el ingenuo. Aquí se remuestrean
# ESTACIONES de entrenamiento con todos sus años y se reajusta el modelo A con el mismo K_b y lambda.
B_BOOT_B <- 300L
filas_est_ent <- split(i_ent, pares$codigo[i_ent])
set.seed(20261121)
boot_b <- replicate(B_BOOT_B, {
  f <- unlist(filas_est_ent[sample(length(filas_est_ent), replace = TRUE)], use.names = FALSE)
  pre <- preparar(f, "A", eleccion$A$K_b); coef_t(pre, resolver(pre, eleccion$A$lambda))
})
b_lo <- apply(boot_b, c(1, 2), quantile, .025); b_hi <- apply(boot_b, c(1, 2), quantile, .975)
comp_bandas <- data.frame(coeficiente = c("b0(t)", "b1(t)"),
                          ancho_medio_sandwich = colMeans(4 * ajuste$A$EE), ancho_medio_bootstrap_estaciones = colMeans(b_hi - b_lo),
                          cociente_anchos = colMeans(b_hi - b_lo) / colMeans(4 * ajuste$A$EE),
                          dias_excluye_0_sandwich = colSums(abs(ajuste$A$B) > 2 * ajuste$A$EE), dias_excluye_0_bootstrap = colSums(b_lo > 0 | b_hi < 0),
                          dias_excluye_1_sandwich = c(NA, sum(abs(ajuste$A$B[, 2] - 1) > 2 * ajuste$A$EE[, 2])),
                          dias_excluye_1_bootstrap = c(NA, sum(b_lo[, 2] > 1 | b_hi[, 2] < 1)))
comp_bandas[, -1] <- round(comp_bandas[, -1], 3)
print(tabla_f(comp_bandas, "C2_bandas_sandwich_vs_bootstrap"))

# ---- 5. ¿EN QUÉ ÉPOCA DEL AÑO EL SATÉLITE EXPLICA MEJOR A LA ESTACIÓN? --------
titulo("5. R²(t) y época del año")
# R²(t) = 1 - Σ_i (y_i(t) - ŷ_i(t))² / Σ_i (y_i(t) - ȳ(t))², con ȳ(t) la media de entrenamiento (referencia fija).
R2_t <- function(P, i) 1 - colSums((Y_r[i, , drop = FALSE] - P)^2) / colSums(sweep(Y_r[i, , drop = FALSE], 2, colMeans(Y_r[i_ent, ]))^2)
R2_tab <- data.frame(dia = T_DIA,
                     R2_A_prueba = R2_t(predecir_raiz(ajuste$A$B, i_pru, "A"), i_pru), R2_B_prueba = R2_t(predecir_raiz(ajuste$B$B, i_pru, "B"), i_pru),
                     R2_A_validacion = R2_t(predecir_raiz(ajuste$A$B, i_val, "A"), i_val),
                     R2_satelite_prueba = R2_t(X_r[i_pru, ], i_pru))     # el satélite sin corregir, tal cual
tabla_f(R2_tab, "C3_R2_por_dia")
mes_de <- cut(T_DIA, c(MESES$ini, 366), labels = MESES$mes, right = FALSE)
por_mes <- data.frame(mes = MESES$mes, R2_A_prueba = tapply(R2_tab$R2_A_prueba, mes_de, mean), R2_satelite_prueba = tapply(R2_tab$R2_satelite_prueba, mes_de, mean),
                      b1_medio = tapply(ajuste$A$B[, 2], mes_de, mean), b1_sobre_EE = tapply(abs(ajuste$A$B[, 2]) / ajuste$A$EE[, 2], mes_de, mean))
por_mes[, -1] <- round(por_mes[, -1], 3); rownames(por_mes) <- NULL
print(tabla_f(por_mes, "C3_R2_y_b1_por_mes"))
cat(sprintf("R²(t) del modelo A en prueba: máximo %.2f el %s, mínimo %.2f el %s | promedio anual %.2f\n",
            max(R2_tab$R2_A_prueba), fecha_dia(which.max(R2_tab$R2_A_prueba)), min(R2_tab$R2_A_prueba), fecha_dia(which.min(R2_tab$R2_A_prueba)), mean(R2_tab$R2_A_prueba)))
figura("C3_R2_por_dia", function() {
  par(mar = c(4.5, 4.5, 3, 1))
  matplot(T_DIA, cbind(R2_tab$R2_A_prueba, R2_tab$R2_B_prueba, R2_tab$R2_satelite_prueba), type = "l", lty = c(1, 1, 2), lwd = 2, col = c("black", "blue", "grey50"),
          ylim = c(min(0, as.matrix(R2_tab[, 2:5])), 1), xlab = "Día del año", ylab = expression(R^2*(t)), main = "Varianza explicada por día, años de prueba (2001-2003)")
  abline(h = 0, lty = 3); legend("bottomleft", c("Modelo A (clase)", "Modelo B (+ altitud y ONI)", "Satélite sin corregir"), col = c("black", "blue", "grey50"), lty = c(1, 1, 2), lwd = 2, bty = "n")
})

figura("C6_concurrente_estilo_clase", function() {
  # La misma figura que la clase arma para el modelo concurrente (Script_Regresion_Funcional.R, marcha):
  # x_i(t), y_i(t), R²(t), b0(t) y b1(t) con bandas, y la correlación de los residuos entre días.
  par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3, 1)); set.seed(3); m <- sample(i_ent, 150)
  matplot(T_DIA, t(X_r[m, ]), type = "l", lty = 1, col = "grey60", xlab = "Día del año", ylab = "raíz(mm/día)", main = "Satélite: x_i(t)")
  matplot(T_DIA, t(Y_r[m, ]), type = "l", lty = 1, col = "grey60", xlab = "Día del año", ylab = "raíz(mm/día)", main = "Estación: y_i(t)")
  plot(T_DIA, R2_tab$R2_A_prueba, type = "l", lwd = 2, ylim = c(min(0, R2_tab$R2_A_prueba), 1), xlab = "Día del año", ylab = expression(R^2*(t)),
       main = "Varianza explicada por el satélite (prueba)"); abline(h = 0, lty = 3)
  for (a in 1:2) {
    b <- ajuste$A$B[, a]; e <- ajuste$A$EE[, a]
    plot(T_DIA, b, type = "n", ylim = range(c(b_lo[, a], b_hi[, a], b - 2 * e, b + 2 * e, 0)), xlab = "Día del año", ylab = "",
         main = c("Intercepto: b0(t)", "Efecto del satélite: b1(t)")[a])
    polygon(c(T_DIA, rev(T_DIA)), c(b_lo[, a], rev(b_hi[, a])), col = adjustcolor("steelblue", .3), border = NA)
    lines(T_DIA, b, lwd = 2); lines(T_DIA, b + 2 * e, lty = 2, col = "red"); lines(T_DIA, b - 2 * e, lty = 2, col = "red")
    abline(h = c(0, if (a == 2) 1), lty = 3)
    if (a == 2) legend("bottomleft", c("Bootstrap de estaciones 95 %", "Sándwich ±2 EE"), fill = c(adjustcolor("steelblue", .3), NA),
                       border = NA, lty = c(NA, 2), col = c(NA, "red"), bty = "n", cex = .8)
  }
  Res_A <- Y_r[i_ent, ] - predecir_raiz(ajuste$A$B, i_ent, "A")
  contour(T_DIA, T_DIA, cor(Res_A), nlevels = 10, xlab = "Día s", ylab = "Día t", main = "Correlación de los residuos")
}, ancho = 14, alto = 8.5)

# ---- 5b. EL MODELO CONCURRENTE AÑO POR AÑO (DESCRIPTIVO) -----------------------
# ¿Cambia con el ENSO cuánto aporta el satélite? Se ajusta y_i(t) = b0(t) + b1(t) x_i(t)
# por separado en cada año, con todas las estaciones de ese año, y se mira b1(t) con su
# banda (modelo concurrente, Ramsay y Silverman 2005, cap. 14). Se hace con los
# 23 años: mínimos cuadrados en cada día t (con curvas suaves, b(t) sale suave) y banda
# ±2 EE punto a punto. Es DESCRIPTIVO: usa todos los años, incluidos validación y prueba,
# y no se usa para predecir ni para escoger nada. Las bandas suponen estaciones independientes
# dentro del año; con la correlación espacial de los residuos (traza-variograma, 5.6b del
# funcional) son algo optimistas.
titulo("5b. Modelo concurrente año por año")
anios_c <- sort(unique(pares$anio[pares$usar]))
por_anio_c <- lapply(anios_c, function(a) {
  i <- which(pares$usar & pares$anio == a); n <- length(i)
  x <- X_r[i, , drop = FALSE]; y <- Y_r[i, , drop = FALSE]
  xc <- sweep(x, 2, colMeans(x)); yc <- sweep(y, 2, colMeans(y))
  b1 <- colSums(xc * yc) / colSums(xc^2); b0 <- colMeans(y) - b1 * colMeans(x)
  res <- y - (rep(1, n) %o% b0 + sweep(x, 2, b1, "*"))
  ee <- sqrt(colSums(res^2) / (n - 2) / colSums(xc^2))
  list(anio = a, n = n, b1 = b1, ee = ee, R2 = 1 - colSums(res^2) / colSums(yc^2))
})
b1_anio <- do.call(rbind, lapply(por_anio_c, function(p) data.frame(anio = p$anio, estaciones = p$n, b1_medio = round(mean(p$b1), 3),
  dias_b1_excluye_0 = sum(p$b1 - 2 * p$ee > 0 | p$b1 + 2 * p$ee < 0), dias_b1_excluye_1 = sum(p$b1 - 2 * p$ee > 1 | p$b1 + 2 * p$ee < 1),
  R2_medio = round(mean(p$R2), 3))))
b1_anio <- merge(unique(pares[, c("anio", "fase", "ONI_anual")]), b1_anio, by = "anio")
tabla_f(b1_anio, "C7_concurrente_por_anio")
print(b1_anio)
ct_b1 <- cor.test(b1_anio$b1_medio, b1_anio$ONI_anual); ct_R2 <- cor.test(b1_anio$R2_medio, b1_anio$ONI_anual)
cat(sprintf("Con el ONI (23 años): b1 medio r = %.2f (p = %.3f) | R² medio r = %.2f (p = %.3f)\n",
            ct_b1$estimate, ct_b1$p.value, ct_R2$estimate, ct_R2$p.value))
print(aggregate(cbind(b1_medio, R2_medio) ~ fase, data = b1_anio, FUN = function(v) round(mean(v), 3)))
lim_b <- c(-0.5, 2)
figura("C7_b1_por_anio", function() {
  par(mfrow = c(6, 4), mar = c(3, 3.4, 2, .6), mgp = c(1.9, .5, 0), oma = c(0, 0, 3, 0))
  for (p in por_anio_c) {
    f <- b1_anio$fase[b1_anio$anio == p$anio]
    plot(NA, xlim = RANGO, ylim = lim_b, xaxt = "n", xlab = "Día del año", ylab = expression(b[1](t)),
         main = sprintf("%d · %s (R² medio %.2f)", p$anio, f, mean(p$R2)), cex.main = .95, cex.axis = .75, cex.lab = .8)
    axis(1, at = c(1, 91, 182, 274, 365), labels = c("ene", "abr", "jul", "oct", "dic"), cex.axis = .75)
    polygon(c(T_DIA, rev(T_DIA)), pmin(pmax(c(p$b1 + 2 * p$ee, rev(p$b1 - 2 * p$ee)), lim_b[1]), lim_b[2]), col = "grey85", border = NA)
    abline(h = 0, lty = 3); abline(h = 1, lty = 2, col = "grey40"); lines(T_DIA, p$b1, lwd = 1.5)
  }
  par(mar = c(3.5, 3.4, 2, .6))
  plot(b1_anio$ONI_anual, b1_anio$b1_medio, pch = 19, xlab = "ONI del año", ylab = "b1 medio del año", cex.axis = .75, cex.lab = .8,
       main = sprintf("b1 medio vs ONI: r = %.2f", ct_b1$estimate), cex.main = .95)
  abline(lm(b1_medio ~ ONI_anual, data = b1_anio), lty = 2)
  mtext("Coeficiente del satélite b1(t) del modelo concurrente ajustado en cada año (banda: ±2 EE punto a punto; punteada: 0; discontinua: 1)",
        outer = TRUE, cex = .95, font = 2)
}, ancho = 11, alto = 15)

# ---- 6. EVALUACIÓN EN PRUEBA CON LAS MÉTRICAS DE LA RÚBRICA 4.4 --------------
titulo("6. Evaluación en prueba (2001-2003)")
kge <- function(P, Yo) {
  r <- sapply(seq_len(nrow(P)), function(k) cor(P[k, ], Yo[k, ]))
  1 - sqrt((r - 1)^2 + (apply(P, 1, sd) / apply(Yo, 1, sd) - 1)^2 + (rowMeans(P) / rowMeans(Yo) - 1)^2)
}
pico <- function(M, ventana) apply(M[, ventana, drop = FALSE], 1, which.max) + min(ventana) - 1
metricas_de <- function(P, i, modelo) {
  Yo <- Y_mm[i, , drop = FALSE]
  data.frame(modelo = modelo, codigo = pares$codigo[i], anio = pares$anio[i], fase = pares$fase[i],
             l2_rel = sqrt(rowSums((P - Yo)^2) / rowSums(Yo^2)), rmse = sqrt(rowMeans((P - Yo)^2)), kge = kge(P, Yo),
             sesgo_pct = 100 * (rowMeans(P) - rowMeans(Yo)) / rowMeans(Yo),
             err_pico_1 = abs(pico(P, VENT_PICO_1) - pico(Yo, VENT_PICO_1)), err_pico_2 = abs(pico(P, VENT_PICO_2) - pico(Yo, VENT_PICO_2)))
}
# IC del 95 % por bootstrap de ESTACIONES (las curvas de una estación no son independientes)
boot_media <- function(v, g, dec = 2) {
  m <- tapply(v, g, mean); q <- quantile(replicate(N_BOOT, mean(sample(m, replace = TRUE))), c(.025, .975))
  sprintf(paste0("%.", dec, "f [%.", dec, "f, %.", dec, "f]"), mean(m), q[1], q[2])
}
resumir <- function(mc) data.frame(curvas = nrow(mc), estaciones = length(unique(mc$codigo)), L2_relativo = boot_media(mc$l2_rel, mc$codigo),
                                   RMSE_L2_mm_dia = boot_media(mc$rmse, mc$codigo), KGE = boot_media(mc$kge, mc$codigo),
                                   sesgo_anual_pct = boot_media(mc$sesgo_pct, mc$codigo, 1), error_pico_1_dias = boot_media(mc$err_pico_1, mc$codigo, 1),
                                   error_pico_2_dias = boot_media(mc$err_pico_2, mc$codigo, 1))
P_A <- predecir_raiz(ajuste$A$B, i_pru, "A"); P_B <- predecir_raiz(ajuste$B$B, i_pru, "B")
recorte <- c(A = 100 * mean(P_A < 0), B = 100 * mean(P_B < 0))
modelos <- list(`Satélite sin corregir` = X_mm[i_pru, ], `Concurrente A (clase)` = a_mm(P_A), `Concurrente B (+ altitud y ONI)` = a_mm(P_B))
mc <- do.call(rbind, lapply(names(modelos), function(nm) metricas_de(modelos[[nm]], i_pru, nm)))
metricas <- do.call(rbind, lapply(split(mc, factor(mc$modelo, levels = names(modelos))), function(d) cbind(modelo = d$modelo[1], resumir(d))))
rownames(metricas) <- NULL
print(tabla_f(metricas, "C4_metricas_prueba"))
cat(sprintf("Reconstrucción en mm no negativa: se recorta en 0 el %.2f %% (A) y el %.2f %% (B) de los valores predichos en escala raíz\n", recorte["A"], recorte["B"]))

# Frente a la regresión función a función del equipo (sección 5.5 de lluvia_valle.R): se lee su tabla, no se recalcula
ref_f <- file.path(paste0("salidas_R_funcional", SUFIJO_C), "tablas", "F5_metricas_prueba.csv")
if (file.exists(ref_f)) {
  ref <- read.csv(ref_f, encoding = "UTF-8", check.names = FALSE)
  cat("\nTabla de referencia del equipo (F5_metricas_prueba.csv):\n"); print(ref[, c("modelo", "L2_relativo", "KGE", "sesgo_anual_pct")], row.names = FALSE)
}

# Diferencia de cada métrica por estación contra el satélite sin corregir, con IC 95 % por
# bootstrap de estaciones (Davison y Hinkley 1997, sec. 3.8). Sin Wilcoxon: el IC ya decide.
por_est <- aggregate(cbind(l2_rel, rmse, kge, sesgo_pct, err_pico_1, err_pico_2) ~ modelo + codigo, data = mc, FUN = mean)
comparar <- function(mod, v, mejor_si_baja = TRUE, absoluto = FALSE) {
  a <- por_est[por_est$modelo == mod, ]; s <- por_est[por_est$modelo == "Satélite sin corregir", ]
  a <- a[order(a$codigo), ]; s <- s[order(s$codigo), ]; stopifnot(identical(a$codigo, s$codigo))
  x <- a[[v]]; z <- s[[v]]; if (absoluto) { x <- abs(x); z <- abs(z) }
  d <- x - z; q <- quantile(replicate(N_BOOT, mean(sample(d, replace = TRUE))), c(.025, .975))
  data.frame(modelo = mod, metrica = if (absoluto) paste0("|", v, "|") else v, satelite = mean(z), concurrente = mean(x), diferencia = mean(d),
             ic_inf = q[1], ic_sup = q[2], estaciones_que_mejoran = sum(if (mejor_si_baja) d < 0 else d > 0), ic_excluye_0 = q[1] > 0 | q[2] < 0,
             mejora = if (mejor_si_baja) mean(d) < 0 else mean(d) > 0)
}
pareada <- do.call(rbind, lapply(c("Concurrente A (clase)", "Concurrente B (+ altitud y ONI)"), function(mod)
  rbind(comparar(mod, "l2_rel"), comparar(mod, "rmse"), comparar(mod, "kge", FALSE), comparar(mod, "sesgo_pct", absoluto = TRUE),
        comparar(mod, "err_pico_1"), comparar(mod, "err_pico_2"))))
pareada[, c("satelite", "concurrente", "diferencia", "ic_inf", "ic_sup")] <- round(pareada[, c("satelite", "concurrente", "diferencia", "ic_inf", "ic_sup")], 3)
rownames(pareada) <- NULL
cat("\nFrente al satélite sin corregir, por estación (", length(unique(por_est$codigo)), " estaciones):\n", sep = "")
print(tabla_f(pareada, "C4_concurrente_vs_satelite_por_estacion"))

# Por fase ENSO. La prueba (2001-2003) no tiene años Niña; se muestra también la validación, que sí los tiene,
# solo como referencia (la validación eligió K_b y lambda, no ajustó los coeficientes).
mc_val <- do.call(rbind, list(metricas_de(X_mm[i_val, ], i_val, "Satélite sin corregir"),
                              metricas_de(a_mm(predecir_raiz(ajuste$A$B, i_val, "A")), i_val, "Concurrente A (clase)"),
                              metricas_de(a_mm(predecir_raiz(ajuste$B$B, i_val, "B")), i_val, "Concurrente B (+ altitud y ONI)")))
mc_val$particion <- "validación (referencia)"; mc$particion <- "prueba"
mc_todo <- rbind(mc[, names(mc_val)], mc_val)
grupos <- unique(mc_todo[, c("particion", "fase", "modelo")])
por_fase <- do.call(rbind, lapply(seq_len(nrow(grupos)), function(g) {
  d <- mc_todo[mc_todo$particion == grupos$particion[g] & mc_todo$fase == grupos$fase[g] & mc_todo$modelo == grupos$modelo[g], ]
  cbind(grupos[g, ], anios = paste(sort(unique(d$anio)), collapse = " "), resumir(d))
}))
por_fase <- por_fase[order(por_fase$particion, por_fase$fase, por_fase$modelo), ]; rownames(por_fase) <- NULL
tabla_f(por_fase, "C4_metricas_por_fase_ENSO")
cat("Métricas por fase ENSO guardadas en C4_metricas_por_fase_ENSO.csv\n")

# ---- 7. CONTRASTE CON fda::fRegress (HERRAMIENTA DE CLASE) -------------------
titulo("7. Contraste con fda::fRegress")
# Mismo modelo A, mismas curvas, mismo K_b y mismo lambda. fRegress usa la penalización int2Lfd(2) = ∫ b''².
# Las diferencias esperadas son pequeñas: aquí la integral se aproxima por trapecios sobre los 365 días.
# Esta sección es solo verificación: si fda no está instalado, el script sigue.
if (requireNamespace("fda", quietly = TRUE)) {
  suppressPackageStartupMessages(library(fda))
  K_b <- eleccion$A$K_b; lam <- eleccion$A$lambda
  base_e <- create.bspline.basis(RANGO, nbasis = K_e, norder = ORDEN); base_s <- create.bspline.basis(RANGO, nbasis = K_s, norder = ORDEN)
  base_b <- create.bspline.basis(RANGO, nbasis = K_b, norder = ORDEN)
  yfd <- fd(t(C_e_r[i_ent, ]), base_e); xfd <- fd(t(C_s_r[i_ent, ]), base_s)
  xlist <- list(const = rep(1, length(i_ent)), sat = xfd)
  blist <- list(const = fdPar(base_b, int2Lfd(2), lam), sat = fdPar(base_b, int2Lfd(2), lam))
  fit_fda <- fRegress(yfd, xlist, blist)
  b0_f <- drop(eval.fd(T_DIA, fit_fda$betaestlist$const$fd)); b1_f <- drop(eval.fd(T_DIA, fit_fda$betaestlist$sat$fd))
  P_f <- sweep(X_r[i_pru, ], 2, b1_f, "*"); P_f <- sweep(P_f, 2, b0_f, "+")
  # Error de ajuste en entrenamiento (la función objetivo sin penalizar) de cada solución, con los mismos pesos
  sse <- function(b0, b1) { P <- sweep(sweep(X_r[i_ent, ], 2, b1, "*"), 2, b0, "+"); sum(sweep((Y_r[i_ent, ] - P)^2, 2, W_T, "*")) }
  centro <- 11:355                                          # sin los 10 días de cada extremo, donde b(t) está menos determinada
  cmp <- data.frame(cantidad = c("SSE de entrenamiento, solución propia", "SSE de entrenamiento, fRegress", "diferencia relativa del SSE (%)",
                                 "RMS de la diferencia de predicción en prueba (escala raíz)", "máx. |b0 propio - b0 fRegress|, centro del año", "máx. |b1 propio - b1 fRegress|, centro del año",
                                 "RMSE L2 en prueba, mm/día (propio)", "RMSE L2 en prueba, mm/día (fRegress)"),
                    valor = c(sse(ajuste$A$B[, 1], ajuste$A$B[, 2]), sse(b0_f, b1_f), 100 * (sse(b0_f, b1_f) / sse(ajuste$A$B[, 1], ajuste$A$B[, 2]) - 1),
                              sqrt(mean((P_A - P_f)^2)), max(abs(ajuste$A$B[centro, 1] - b0_f[centro])), max(abs(ajuste$A$B[centro, 2] - b1_f[centro])),
                              mean(rmse_l2(a_mm(P_A), i_pru)), mean(rmse_l2(a_mm(P_f), i_pru))))
  print(tabla_f(cmp, "C5_contraste_fRegress"))
  cat("Nota: con poca penalización y muchos parámetros fRegress integra con tolerancia numérica (inprod) y queda algo por encima del mínimo exacto de la solución propia;\n",
      "      con penalización fuerte coinciden. La solución propia se verifica contra mínimos cuadrados directos en la sección 8.\n")
} else {
  cat("fda no está instalado: se omite el contraste (instalar con install.packages('fda')).\n")
}

# ---- 8. PRUEBAS ---------------------------------------------------------------
titulo("8. Pruebas")
pruebas <- character(0)
chk <- function(nombre, cond) { pruebas <<- c(pruebas, sprintf("%s %s", if (isTRUE(cond)) "OK  " else "FALLA", nombre)); stopifnot(isTRUE(cond)) }
chk("las filas de pares, coeficientes y curvas están alineadas", nrow(pares) == nrow(Y_r) && nrow(Y_r) == nrow(X_r))
chk("el b1(t) por año coincide con lm() en un año y un día (diferencia < 1e-10), y hay un ajuste por cada uno de los 23 años", {
  p <- por_anio_c[[which(anios_c == 1999)]]; i <- which(pares$usar & pares$anio == 1999)
  cf <- summary(lm(Y_r[i, 105] ~ X_r[i, 105]))$coefficients
  length(por_anio_c) == 23 && abs(cf[2, 1] - p$b1[105]) < 1e-10 && abs(cf[2, 2] - p$ee[105]) < 1e-10 })
chk("la partición temporal no se cruza (entrenamiento < validación < prueba)",
    max(pares$anio[ent]) < min(pares$anio[val]) && max(pares$anio[val]) < min(pares$anio[pru]))
chk("el ajuste usa solo filas de entrenamiento (la matriz normal no cambia si se alteran validación y prueba)", {
  X_bak <- X_r; Y_bak <- Y_r; X_r[c(i_val, i_pru), ] <- 0; Y_r[c(i_val, i_pru), ] <- 0
  th2 <- resolver(preparar(i_ent, "A", eleccion$A$K_b), eleccion$A$lambda); X_r <- X_bak; Y_r <- Y_bak
  max(abs(th2 - ajuste$A$theta)) < 1e-10 })
chk("las ecuaciones normales dan lo mismo que mínimos cuadrados directos sobre todas las filas (curva, día), diferencia < 1e-8", {
  K_t <- 12L; pre_t <- preparar(i_ent, "A", K_t); n_t <- length(i_ent); Pb <- base_bspline(T_DIA, K_t); sw <- sqrt(W_T[rep(T_DIA, n_t)])
  Xd <- cbind(Pb[rep(T_DIA, n_t), ] * sw, (Pb[rep(T_DIA, n_t), ] * as.vector(t(X_r[i_ent, ]))) * sw)
  max(abs(qr.coef(qr(Xd), as.vector(t(Y_r[i_ent, ])) * sw) - resolver(pre_t, 0))) < 1e-8 })
chk("la rugosidad del equipo y la propia coinciden para una recta (R b = 0)", {
  Pb <- base_bspline(T_DIA, 8L); th_recta <- solve(crossprod(Pb), crossprod(Pb, 3 + 0.01 * T_DIA)); max(abs(gram_rugosidad(8L) %*% th_recta)) < 1e-6 })
chk("con lambda enorme b1(t) tiende a una recta", {
  pre <- preparar(i_ent, "A", 12L); b1 <- coef_t(pre, resolver(pre, 1e14))[, 2]; max(abs(diff(b1, differences = 2))) < 1e-6 })
chk("la reconstrucción en mm es no negativa", all(modelos[[2]] >= 0) && all(modelos[[3]] >= 0))
chk("bandas de b(t) por bootstrap de estaciones con al menos 300 remuestreos", B_BOOT_B >= 300)
chk("el modelo A baja el error L² frente al satélite sin corregir en prueba (IC 95 % por estación excluye el 0)",
    with(pareada[pareada$modelo == "Concurrente A (clase)" & pareada$metrica == "l2_rel", ], mejora && ic_excluye_0))
cat(paste(pruebas, collapse = "\n"), "\n")
cat(sprintf("\nPruebas: %d de %d\n", sum(grepl("^OK", pruebas)), length(pruebas)))
cat("Listo. Tablas en", file.path(DIR_OUT, "tablas"), "y figuras en", file.path(DIR_OUT, "figuras"), "\n")

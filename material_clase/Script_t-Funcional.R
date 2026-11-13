#Script t-funcional


#Limpiar área de trabajo
rm(list = ls())


## Librerias
library(data.table)
library(openxlsx)
library(readxl)
library(tidyverse)
library(ggplot2)
library(fda)


#Directorio de trabajo (ajustar a la carpeta donde están los datos)
#setwd("D:/X_ENME-UT/Ponencia_Calidad_Aire")

#Semilla para reproducir las permutaciones aleatorias
set.seed(2026)

#Estaciones a comparar (nombres de columna del archivo) y etiquetas para las figuras
Estacion1 <- "Univalle"
Estacion2 <- "Pance"
Nombre1   <- "Univalle"
Nombre2   <- "Pance"



#Cargar datos

#Cargar datos ("ND" y "NA" = dato no disponible)
ICA_2010_2018 <- read_excel("Indice_de_calidad_del_aire__ICA__de_Santiago_de_Cali.xlsx",
                            col_types = c("date", "numeric", "numeric", "numeric", "numeric",
                                          "numeric", "numeric", "numeric", "numeric", "numeric"),
                            na = c("ND", "NA", ""))

#Convertir a data.table
ICA_2010_2018 <- as.data.table(ICA_2010_2018)


#Recodificar la variable fecha
#"tiempo" es el día del año
ICA_2010_2018[, año    := format(Fecha, "%Y")]
ICA_2010_2018[, mes    := format(Fecha, "%m")]
ICA_2010_2018[, dia    := as.numeric(format(Fecha, "%d"))]
ICA_2010_2018[, tiempo := as.numeric(format(Fecha, "%j"))]

head(ICA_2010_2018)



#Análisis exploratorio de datos

ICA_2010_2018 %>%
  filter(!is.na(.data[[Estacion1]])) %>%
  group_by(año) %>%
  summarize(Prom = mean(.data[[Estacion1]], na.rm = TRUE),
            Min = min(.data[[Estacion1]], na.rm = TRUE),
            Max = max(.data[[Estacion1]], na.rm = TRUE),
            CV = sd(.data[[Estacion1]], na.rm = TRUE)/mean(.data[[Estacion1]], na.rm = TRUE)*100)


ICA_2010_2018 %>%
  filter(!is.na(.data[[Estacion2]])) %>%
  group_by(año) %>%
  summarize(Prom = mean(.data[[Estacion2]], na.rm = TRUE),
            Min = min(.data[[Estacion2]], na.rm = TRUE),
            Max = max(.data[[Estacion2]], na.rm = TRUE),
            CV = sd(.data[[Estacion2]], na.rm = TRUE)/mean(.data[[Estacion2]], na.rm = TRUE)*100)


#curva anual con huecos muy largos no se puede reconstruir de forma confiable.
#Se calcula, por año, el número de días con dato y la racha más larga de días faltantes.
racha_max <- function(x){
  r <- rle(is.na(x))
  if(any(r$values)) max(r$lengths[r$values]) else 0
}

Disponibilidad <- ICA_2010_2018 %>%
  group_by(año) %>%
  summarize(n_1     = sum(!is.na(.data[[Estacion1]])),
            Hueco_1 = racha_max(.data[[Estacion1]]),
            n_2     = sum(!is.na(.data[[Estacion2]])),
            Hueco_2 = racha_max(.data[[Estacion2]]))

Disponibilidad


#Criterio de inclusión: se usan los años cuyo hueco más largo no supera 45 días
hueco_max <- 45

años_1 <- Disponibilidad$año[Disponibilidad$Hueco_1 <= hueco_max]
años_2 <- Disponibilidad$año[Disponibilidad$Hueco_2 <= hueco_max]

años_1
años_2

n1 <- length(años_1)
n2 <- length(años_2)
n  <- n1 + n2


#Diagrama de cajas mensual del ICA por estación (años seleccionados)
Datos_largo <- rbind(
  data.frame(Estacion = Nombre1, mes = ICA_2010_2018$mes, año = ICA_2010_2018$año, ICA = ICA_2010_2018[[Estacion1]]),
  data.frame(Estacion = Nombre2, mes = ICA_2010_2018$mes, año = ICA_2010_2018$año, ICA = ICA_2010_2018[[Estacion2]])) %>%
  filter((Estacion == Nombre1 & año %in% años_1) | (Estacion == Nombre2 & año %in% años_2))

pdf("Fig1_Boxplot_Mensual.pdf", width = 9, height = 4.5)
ggplot(Datos_largo, aes(x = mes, y = ICA, fill = Estacion)) +
  geom_boxplot(outlier.size = 0.6, na.rm = TRUE) +
  geom_hline(yintercept = 50, colour = "red", linetype = 2) +
  labs(x = "Mes", y = "ICA diario", fill = "Estación") +
  theme_bw()
dev.off()



#Construcción de datos funcionales

#Base B-spline cúbica (orden 4) con nodos semanales sobre [1, 366]
rangval <- c(1, 366)
nodos   <- seq(1, 366, length.out = 53)
BaseICA <- create.bspline.basis(rangval, norder = 4, breaks = nodos)
BaseICA$nbasis

#Como cada año tiene días faltantes distintos, cada curva se ajusta por separado
#con los días observados. La penalización sobre la segunda derivada (lambda)
#estabiliza la curva en los huecos.
Suavizar_ICA <- function(estacion, años, lambda){
  fdParICA <- fdPar(BaseICA, int2Lfd(2), lambda)
  Coefs <- NULL
  GCV   <- 0
  for(a in años){
    Datos_a <- ICA_2010_2018[año == a & !is.na(get(estacion))]
    Ajuste  <- smooth.basis(Datos_a$tiempo, Datos_a[[estacion]], fdParICA)
    Coefs   <- cbind(Coefs, Ajuste$fd$coefs)
    GCV     <- GCV + sum(Ajuste$gcv)
  }
  colnames(Coefs) <- años
  list(fd = fd(Coefs, BaseICA, list("Día del año", "Año", "ICA")), gcv = GCV)
}


#Selección de lambda por validación cruzada generalizada (VCG), común a ambas estaciones
log_lambda <- seq(0, 4, by = 0.5)
VCG <- sapply(log_lambda, function(l) Suavizar_ICA(Estacion1, años_1, 10^l)$gcv +
                                      Suavizar_ICA(Estacion2, años_2, 10^l)$gcv)

pdf("Fig2_VCG_lambda.pdf", width = 6, height = 4.5)
plot(log_lambda, VCG, type = "b", pch = 19, xlab = expression(log[10](lambda)), ylab = "VCG total",
     main = "Selección del parámetro de suavizado")
dev.off()

lambda_opt <- 10^log_lambda[which.min(VCG)]
lambda_opt

ICA1_fd <- Suavizar_ICA(Estacion1, años_1, lambda_opt)$fd
ICA2_fd <- Suavizar_ICA(Estacion2, años_2, lambda_opt)$fd


#Malla de evaluación: días 1 a 365 (todas las curvas comparten este dominio)
t_eval <- 1:365
M      <- length(t_eval)

X1 <- eval.fd(t_eval, ICA1_fd)         #matriz 365 x n1
X2 <- eval.fd(t_eval, ICA2_fd)         #matriz 365 x n2
colnames(X1) <- años_1
colnames(X2) <- años_2


#Datos observados y datos funcionales
pdf("Fig3_Datos_Funcionales.pdf", width = 10, height = 7)
par(mfrow = c(2, 2))
matplot(t_eval, sapply(años_1, function(a) ICA_2010_2018[año == a][[Estacion1]][t_eval]), type = "l", lty = 1,
        col = "grey60", xlab = "Día del año", ylab = "ICA", main = paste("Datos observados -", Nombre1), ylim = c(0, 110))
matplot(t_eval, X1, type = "l", lty = 1, col = "grey40", xlab = "Día del año", ylab = "ICA",
        main = paste("Datos funcionales -", Nombre1), ylim = c(0, 110))
matplot(t_eval, sapply(años_2, function(a) ICA_2010_2018[año == a][[Estacion2]][t_eval]), type = "l", lty = 1,
        col = "grey60", xlab = "Día del año", ylab = "ICA", main = paste("Datos observados -", Nombre2), ylim = c(0, 110))
matplot(t_eval, X2, type = "l", lty = 1, col = "grey40", xlab = "Día del año", ylab = "ICA",
        main = paste("Datos funcionales -", Nombre2), ylim = c(0, 110))
dev.off()



#Promedio y DE funcional
Media1 <- rowMeans(X1)
Media2 <- rowMeans(X2)
DE1    <- apply(X1, 1, sd)
DE2    <- apply(X2, 1, sd)

#La media también se obtiene con fda: coinciden
max(abs(Media1 - eval.fd(t_eval, mean.fd(ICA1_fd))))


pdf("Fig4_Media_DE.pdf", width = 10, height = 4.5)
par(mfrow = c(1, 2))
matplot(t_eval, cbind(Media1, Media2), type = "l", lty = 1, lwd = 2, col = c("red", "blue"),
        xlab = "Día del año", ylab = "ICA", main = "Media funcional",
        ylim = c(min(Media1, Media2), max(Media1, Media2)*1.15))
legend("top", c(Nombre1, Nombre2), col = c("red", "blue"), lwd = 2, bty = "n", horiz = TRUE)
matplot(t_eval, cbind(DE1, DE2), type = "l", lty = 1, lwd = 2, col = c("red", "blue"),
        xlab = "Día del año", ylab = "ICA", main = "Desviación estándar funcional",
        ylim = c(0, max(DE1, DE2)*1.2))
legend("top", c(Nombre1, Nombre2), col = c("red", "blue"), lwd = 2, bty = "n", horiz = TRUE)
dev.off()



#Prueba t funcional

#Función que calcula, en cada día, el estadístico t con varianza combinada
#y el estadístico de tipo Welch T(t) (el que usa tperm.fd)
Estadisticos_t <- function(A, B){
  nA <- ncol(A); nB <- ncol(B)
  mA <- rowMeans(A); mB <- rowMeans(B)
  vA <- rowSums((A - mA)^2)/(nA - 1)
  vB <- rowSums((B - mB)^2)/(nB - 1)
  sp2 <- ((nA - 1)*vA + (nB - 1)*vB)/(nA + nB - 2)
  list(t_comb  = (mA - mB)/sqrt(sp2*(1/nA + 1/nB)),
       T_welch = abs(mA - mB)/sqrt(vA/nA + vB/nB))
}

Est_obs <- Estadisticos_t(X1, X2)

#Bajo normalidad, t(t) ~ t_{n-2} en cada día fijo
p_t_puntual <- 2*pt(abs(Est_obs$t_comb), df = n - 2, lower.tail = FALSE)

pdf("Fig5_t_Puntual.pdf", width = 10, height = 4.5)
par(mfrow = c(1, 2))
plot(t_eval, Est_obs$t_comb, type = "l", lwd = 2, xlab = "Día del año", ylab = "t(t)",
     main = "Estadístico t puntual (varianza combinada)")
abline(h = c(-1, 1)*qt(0.975, n - 2), col = "red", lty = 2)
plot(t_eval, p_t_puntual, type = "l", lwd = 2, xlab = "Día del año", ylab = "Valor p",
     main = "Valor p puntual (t de Student)", ylim = c(0, 1))
abline(h = 0.05, col = "red", lty = 2)
dev.off()

#Proporción de días con valor p < 0.05 (sin corrección por multiplicidad)
mean(p_t_puntual < 0.05)*100



#Prueba t funcional con permutaciones

#Con n = n1 + n2 curvas hay choose(n, n1) asignaciones posibles de las etiquetas.
#Si ese número es manejable se recorren todas (prueba exacta); si no, se toman B al azar.
choose(n, n1)

X_todas <- cbind(X1, X2)

if(choose(n, n1) <= 20000){
  Grupos1 <- combn(n, n1)                                      #todas las asignaciones
} else {
  Grupos1 <- replicate(5000, sort(sample(n, n1)))              #permutaciones aleatorias
}
B <- ncol(Grupos1)
B

#En la prueba exacta, el menor valor p posible es 1/B (2/B si n1 = n2, porque intercambiar
#los grupos completos produce el mismo |T(t)|)

#Matriz M x B de estadísticos bajo H0
T_nulos <- matrix(NA, M, B)
for(b in 1:B){
  g1 <- Grupos1[, b]
  T_nulos[, b] <- Estadisticos_t(X_todas[, g1], X_todas[, -g1])$T_welch
}

#Máximo de cada permutación
T_max_nulos <- apply(T_nulos, 2, max)

#Valores observados
T_obs     <- Est_obs$T_welch
T_max_obs <- max(T_obs)

#Curva crítica puntual y valor crítico máximo al 5 %
alfa       <- 0.05
c_puntual  <- apply(T_nulos, 1, quantile, probs = 1 - alfa)
c_maximo   <- quantile(T_max_nulos, 1 - alfa)

#Valores p puntual y global
#(B incluye la asignación observada, por eso no se suma 1 en este caso exacto)
p_puntual <- rowMeans(T_nulos >= T_obs - 1e-12)
p_global  <- mean(T_max_nulos >= T_max_obs - 1e-12)

c_maximo
p_global


#Días con diferencia significativa
Dias_a_fecha <- function(dias) format(as.Date(dias - 1, origin = "2019-01-01"), "%d-%b")

Intervalos <- function(indicador){
  r <- rle(indicador)
  fin <- cumsum(r$lengths); ini <- fin - r$lengths + 1
  data.frame(Desde = Dias_a_fecha(t_eval[ini[r$values]]),
             Hasta = Dias_a_fecha(t_eval[fin[r$values]]),
             Dias  = r$lengths[r$values])
}

#Diferencias puntuales (sin corrección por multiplicidad)
Intervalos(p_puntual < alfa)

#Diferencias con control del error en todo el año (criterio del máximo)
Intervalos(T_obs > c_maximo)


pdf("Fig6_Permutacion_tFuncional.pdf", width = 10, height = 4.5)
par(mfrow = c(1, 2))
plot(t_eval, p_puntual, type = "l", lwd = 2, col = "blue", xlab = "Día del año", ylab = "Valor p",
     main = "Valor p puntual (permutación)", ylim = c(0, 1))
abline(h = alfa, col = "red", lty = 2)
plot(t_eval, T_obs, type = "l", lwd = 2, col = "red", xlab = "Día del año", ylab = "T(t)",
     main = "Prueba t funcional por permutación", ylim = c(0, max(T_obs, c_maximo)*1.6))
lines(t_eval, c_puntual, col = "blue", lty = 3, lwd = 2)
abline(h = c_maximo, col = "black", lty = 2)
legend("topleft", c("Estadístico observado", "Valor crítico puntual (0.05)", "Valor crítico máximo"),
       col = c("red", "blue", "black"), lty = c(1, 3, 2), lwd = c(2, 2, 1), bty = "n", cex = 0.75)
dev.off()


#Comparación con la función tperm.fd() de fda (permutaciones aleatorias)
#tperm.fd usa el mismo estadístico T(t) y el mismo criterio del máximo
tperm_ICA <- tperm.fd(ICA1_fd, ICA2_fd, nperm = 1000, q = alfa,
                      argvals = t_eval, plotres = FALSE)

c(p_global_manual = p_global, p_global_tperm = tperm_ICA$pval)
c(c_max_manual = unname(c_maximo), c_max_tperm = unname(tperm_ICA$qval))



#Pruebas globales: normal L2 y tipo F

#Pesos de la regla del trapecio para aproximar las integrales sobre la malla
w <- rep(1, M); w[c(1, M)] <- 0.5

#Covarianza combinada en la malla (N = n - 2 grados de libertad)
R1 <- X1 - Media1
R2 <- X2 - Media2
N_gl  <- n - 2
G_comb <- (R1 %*% t(R1) + R2 %*% t(R2))/N_gl         #matriz 365 x 365

#Trazas: tr(gamma) y tr(gamma^{⊗2})
DG        <- G_comb*w                                 #equivale a diag(w) %*% G
tr_gamma  <- sum(diag(DG))
tr_gamma2 <- sum(DG*t(DG))                            #tr((D G)^2)

#Estimadores con corrección de sesgo (Zhang, 2013)
tr_gamma2_c  <- N_gl^2/((N_gl - 1)*(N_gl + 2))*(tr_gamma2 - tr_gamma^2/N_gl)
tr2_gamma_c  <- N_gl*(N_gl + 1)/((N_gl - 1)*(N_gl + 2))*(tr_gamma^2 - 2*tr_gamma2/(N_gl + 1))

beta_hat <- tr_gamma2_c/tr_gamma
d_hat    <- tr2_gamma_c/tr_gamma2_c

#Sin corrección, para comparar
d_sust <- tr_gamma^2/tr_gamma2

c(beta = beta_hat, d = d_hat, d_sustitucion = d_sust)


#Estadístico L2: T_n = (n1 n2 / n) ||X1_bar - X2_bar||^2
T_L2 <- n1*n2/n*sum(w*(Media1 - Media2)^2)
p_L2 <- pchisq(T_L2/beta_hat, df = d_hat, lower.tail = FALSE)

#Estadístico de tipo F: F_n = T_n / tr(gamma_hat) ~ F(d, (n-2) d)
F_n <- T_L2/tr_gamma
p_F <- pf(F_n, d_hat, N_gl*d_hat, lower.tail = FALSE)

c(T_L2 = T_L2, p_L2 = p_L2, F_n = F_n, p_F = p_F)


#Calibración de F_n por permutación (con las mismas asignaciones de la Sección 6)
F_perm <- numeric(B)
for(b in 1:B){
  g1 <- Grupos1[, b]
  A  <- X_todas[, g1]; Bm <- X_todas[, -g1]
  mA <- rowMeans(A);   mB <- rowMeans(Bm)
  SSE_b <- sum(w*(rowSums((A - mA)^2) + rowSums((Bm - mB)^2)))
  F_perm[b] <- (n1*n2/n*sum(w*(mA - mB)^2))/(SSE_b/N_gl)
}
p_F_perm <- mean(F_perm >= F_n - 1e-12)
p_F_perm


#Valores propios de la covarianza combinada: explican por qué d es pequeño
Desc  <- eigen(sqrt(w)*t(sqrt(w)*G_comb), symmetric = TRUE)    #D^{1/2} G D^{1/2}
lambda_r <- pmax(Desc$values, 0)
Var_acum <- cumsum(lambda_r)/sum(lambda_r)
round(Var_acum[1:6], 3)

pdf("Fig7_Valores_Propios_F.pdf", width = 10, height = 4.5)
par(mfrow = c(1, 2))
barplot(lambda_r[1:10]/sum(lambda_r)*100, names.arg = 1:10, col = "grey70",
        xlab = "Componente", ylab = "% de varianza", main = "Valores propios de la covarianza combinada")
hist(F_perm, breaks = 40, freq = FALSE, col = "grey85", border = "white",
     xlab = expression(F[n]), main = "Distribución de permutación de F")
curve(df(x, d_hat, N_gl*d_hat), add = TRUE, col = "blue", lwd = 2)
abline(v = F_n, col = "red", lwd = 2)
legend("topright", c(paste0("F(", round(d_hat, 2), ", ", round(N_gl*d_hat, 2), ")"), "F observado"),
       col = c("blue", "red"), lwd = 2, bty = "n", cex = 0.8)
dev.off()



#Prueba de hotelling sobre PCA funcional

#Funciones propias (normalizadas en L2) y número de componentes K (80 % de varianza)
phi_r <- Desc$vectors/sqrt(w)
K     <- which(Var_acum >= 0.80)[1]
K

#Proyección de la diferencia de medias y estadístico T^2_K
delta_r <- colSums(w*(Media1 - Media2)*phi_r[, 1:K, drop = FALSE])
T2_K    <- n1*n2/n*sum(delta_r^2/lambda_r[1:K])

#Aproximación chi-cuadrado (asintótica) y F de Hotelling (muestras pequeñas)
p_T2_chi <- pchisq(T2_K, df = K, lower.tail = FALSE)
p_T2_F   <- pf((n - K - 1)/((n - 2)*K)*T2_K, K, n - K - 1, lower.tail = FALSE)

c(T2_K = T2_K, p_chi2 = p_T2_chi, p_F = p_T2_F)



#prueba pareada (es para datos relacionados)
#mis años en ambas estaciones

#Se usan los años comunes a ambas estacione
#comparar el mismo año elimina la variación

#interanual compartida por ambas (por ejemplo, la asociada a El Niño 2015-2016).
años_comunes <- intersect(años_1, años_2)
años_comunes
n_p <- length(años_comunes)

D <- X1[, años_comunes] - X2[, años_comunes]       #curvas de diferencias Estacion1 - Estacion2
D_media <- rowMeans(D)
D_de    <- apply(D, 1, sd)

#Estadístico t pareado en cada día
t_par <- D_media/(D_de/sqrt(n_p))

#Prueba exacta por cambio de signo: bajo H0, D_j y -D_j son igualmente probables.
#Hay 2^n_p combinaciones de signos.
Signos <- as.matrix(expand.grid(rep(list(c(-1, 1)), n_p)))
B_p    <- nrow(Signos)
B_p

#Cambiar todos los signos produce el mismo |t(t)|, así que el menor valor p posible es 2/2^n_p.
#Con pocos años, la prueba pareada exacta puede no alcanzar el nivel del 5 %.
2/B_p

T_par_nulos <- matrix(NA, M, B_p)
for(b in 1:B_p){
  Db <- sweep(D, 2, Signos[b, ], "*")
  T_par_nulos[, b] <- abs(rowMeans(Db))/(apply(Db, 1, sd)/sqrt(n_p))
}

T_par_obs     <- abs(t_par)
c_puntual_par <- apply(T_par_nulos, 1, quantile, probs = 1 - alfa)
c_maximo_par  <- quantile(apply(T_par_nulos, 2, max), 1 - alfa)
p_puntual_par <- rowMeans(T_par_nulos >= T_par_obs - 1e-12)
p_global_par  <- mean(apply(T_par_nulos, 2, max) >= max(T_par_obs) - 1e-12)

p_global_par
Intervalos(p_puntual_par < alfa)
Intervalos(T_par_obs > c_maximo_par)


#Prueba F pareada: F_n = n ||D_bar||^2 / tr(gamma_D) ~ F(d, (n-1) d)
G_D   <- (D - D_media) %*% t(D - D_media)/(n_p - 1)
DG_D  <- G_D*w
trD   <- sum(diag(DG_D)); trD2 <- sum(DG_D*t(DG_D))
N_D   <- n_p - 1
trD2_c <- N_D^2/((N_D - 1)*(N_D + 2))*(trD2 - trD^2/N_D)
tr2D_c <- N_D*(N_D + 1)/((N_D - 1)*(N_D + 2))*(trD^2 - 2*trD2/(N_D + 1))
d_D    <- tr2D_c/trD2_c

F_par <- n_p*sum(w*D_media^2)/trD
p_F_par <- pf(F_par, d_D, N_D*d_D, lower.tail = FALSE)
c(F_pareada = F_par, d = d_D, p = p_F_par)


pdf("Fig8_Prueba_Pareada.pdf", width = 10, height = 4.5)
par(mfrow = c(1, 2))
matplot(t_eval, D, type = "l", lty = 1, col = "grey60", xlab = "Día del año", ylab = paste(Nombre1, "-", Nombre2),
        main = "Diferencias pareadas por año")
lines(t_eval, D_media, lwd = 2, col = "red")
lines(t_eval, D_media + 2*D_de/sqrt(n_p), lty = 2, col = "red")
lines(t_eval, D_media - 2*D_de/sqrt(n_p), lty = 2, col = "red")
abline(h = 0, lty = 3)
plot(t_eval, T_par_obs, type = "l", lwd = 2, col = "red", xlab = "Día del año", ylab = "|t(t)|",
     main = "Prueba pareada por cambio de signo", ylim = c(0, max(T_par_obs, c_maximo_par)*1.6))
lines(t_eval, c_puntual_par, col = "blue", lty = 3, lwd = 2)
abline(h = c_maximo_par, col = "black", lty = 2)
legend("topleft", c("Estadístico observado", "Valor crítico puntual (0.05)", "Valor crítico máximo"),
       col = c("red", "blue", "black"), lty = c(1, 3, 2), lwd = c(2, 2, 1), bty = "n", cex = 0.75)
dev.off()


#Actividad de estudio independiente

#1. Compare ERA-Obrero con La Flora (Estacion1 <- "ERA", Estacion2 <- "Flora"). Esas estaciones
#   aportan más años. ¿Hay diferencias? Compare las conclusiones de la prueba para muestras
#   independientes y de la prueba pareada, y explique el papel de la variación interanual.
#2. Repita el análisis con hueco_max = 30 y hueco_max = 60. ¿Cambian los años incluidos y las conclusiones?
#3. Cambie lambda a 10 veces y a 1/10 del valor elegido por la VCG. ¿Cómo cambian d_hat y la región
#   significativa? Relacione el resultado con el número efectivo de grados de libertad.
#4. Reemplace T(t) por el estadístico con varianza combinada en la Sección 6. ¿Se modifican los resultados?
#5. ¿Cuál es el menor valor p que puede alcanzar la prueba pareada por cambio de signo con 5 años?
#   ¿Cuántos años se necesitarían para poder rechazar H0 al 5 %?
#6. Construya curvas mensuales (día 1 a 31, como en data_ERA.xlsx) y compare la temporada seca
#   (enero, febrero, julio y agosto) con la temporada de lluvias (abril, mayo, octubre y noviembre) en ERA-Obrero.

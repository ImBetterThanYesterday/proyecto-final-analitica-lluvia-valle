#Script ANOVA funcional

#Limpiar área de trabajo
rm(list = ls())


## Librerias
library(data.table)
library(openxlsx)
library(readxl)
library(tidyverse)
library(ggplot2)
library(fda)
library(fdANOVA)


#Directorio de trabajo (ajustar a la carpeta donde están los datos)
#setwd("D:/X_ENME-UT/Ponencia_Calidad_Aire")

#Semilla para reproducir las permutaciones y simulaciones
set.seed(2026)

#Estaciones candidatas (nombres de columna del archivo)
Estaciones_cand <- c("ERA", "Flora", "Pance", "Univalle", "Base_Aerea",
                     "Compartir", "La_Ermita", "Cañaveralejo", "Transitoria")



#Cargar datos

#Cargar datos ("ND" y "NA" = dato no disponible)
ICA_2010_2018 <- read_excel("Indice_de_calidad_del_aire__ICA__de_Santiago_de_Cali.xlsx",
                            col_types = c("date", "numeric", "numeric", "numeric", "numeric",
                                          "numeric", "numeric", "numeric", "numeric", "numeric"),
                            na = c("ND", "NA", ""))

#Convertir a data.table
ICA_2010_2018 <- as.data.table(ICA_2010_2018)

#Un ICA igual a 0 no es físicamente plausible (falla del equipo)
Columnas_ICA <- c("ERA", "Flora", "Pance", "Univalle", "Base_Aerea",
                  "Compartir", "La_Ermita", "Cañaveralejo", "Transitoria")
ICA_2010_2018[, (Columnas_ICA) := lapply(.SD, function(v) replace(v, v == 0, NA)), .SDcols = Columnas_ICA]


#Recodificar la variable fecha
ICA_2010_2018[, año    := format(Fecha, "%Y")]
ICA_2010_2018[, mes    := format(Fecha, "%m")]
ICA_2010_2018[, dia    := as.numeric(format(Fecha, "%d"))]
ICA_2010_2018[, tiempo := as.numeric(format(Fecha, "%j"))]

head(ICA_2010_2018)



#preprocesamiento

#Racha más larga de días faltantes en un vector
racha_max <- function(x){
  r <- rle(is.na(x))
  if(any(r$values)) max(r$lengths[r$values]) else 0
}

#Días sin dato al inicio o al final del año. En esos extremos la curva no se interpola
#sino que se extrapola, lo que es mucho menos confiable que rellenar un hueco interior.
borde_max <- function(x){
  obs <- which(!is.na(x))
  if(length(obs) == 0) return(length(x))
  max(obs[1] - 1, length(x) - obs[length(obs)])
}

#Datos con estructura long  por estación y año
Disponibilidad <- ICA_2010_2018 %>%
  pivot_longer(cols = all_of(Estaciones_cand), names_to = "Estacion", values_to = "ICA") %>%
  group_by(Estacion, año) %>%
  summarize(n_dias = sum(!is.na(ICA)),
            Hueco  = racha_max(ICA),
            Borde  = borde_max(ICA), .groups = "drop")

#Hueco más largo por estación y año
Disponibilidad %>%
  dplyr::select(Estacion, año, Hueco) %>%   #dplyr:: porque fda carga MASS, que también tiene select()
  pivot_wider(names_from = año, values_from = Hueco)


#Criterios:
#un año entra si su hueco más largo no supera 45 días y si le faltan como máximo
#7 días al inicio o al final del año;
#una estación entra si tiene al menos 4 años válidos.
hueco_max <- 45
borde_lim <- 7
min_años  <- 4

Seleccion <- Disponibilidad %>% filter(Hueco <= hueco_max, Borde <= borde_lim)

Estaciones <- Seleccion %>%
  count(Estacion) %>%
  filter(n >= min_años) %>%
  pull(Estacion)

Estaciones <- Estaciones_cand[Estaciones_cand %in% Estaciones]     #conservar el orden original

Seleccion <- Seleccion %>% filter(Estacion %in% Estaciones)

#Años incluidos por estación
Seleccion %>% group_by(Estacion) %>% summarize(Años = paste(año, collapse = ", "), n_i = n())

k <- length(Estaciones)
k


#Estadísticas descriptivas de los registros diarios incluidos
Datos_largo <- ICA_2010_2018 %>%
  pivot_longer(cols = all_of(Estaciones), names_to = "Estacion", values_to = "ICA") %>%
  semi_join(Seleccion, by = c("Estacion", "año"))

Datos_largo %>%
  group_by(Estacion) %>%
  summarize(Prom = mean(ICA, na.rm = TRUE),
            Min  = min(ICA, na.rm = TRUE),
            Max  = max(ICA, na.rm = TRUE),
            CV   = sd(ICA, na.rm = TRUE)/mean(ICA, na.rm = TRUE)*100)


pdf("Fig1_Boxplot_Estaciones.pdf", width = 9, height = 4.5)
ggplot(Datos_largo, aes(x = factor(Estacion, levels = Estaciones), y = ICA, fill = Estacion)) +
  geom_boxplot(outlier.size = 0.6, na.rm = TRUE, show.legend = FALSE) +
  geom_hline(yintercept = 50, colour = "red", linetype = 2) +
  labs(x = "Estación", y = "ICA diario") +
  theme_bw()
dev.off()



#Construcción de datos funcionales

#Base B-spline cúbica (orden 4) con nodos semanales sobre [1, 366]
rangval <- c(1, 366)
nodos   <- seq(1, 366, length.out = 53)
BaseICA <- create.bspline.basis(rangval, norder = 4, breaks = nodos)

#Cada curva se ajusta con los días observados de su año; la penalización
#sobre la segunda derivada (lambda) estabiliza la curva en los huecos.
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
  list(fd = fd(Coefs, BaseICA), gcv = GCV)
}

Años_estacion <- lapply(Estaciones, function(e) Seleccion$año[Seleccion$Estacion == e])
names(Años_estacion) <- Estaciones


#Selección de lambda por VCG, común a todas las estaciones
log_lambda <- seq(0, 4, by = 0.5)
VCG <- sapply(log_lambda, function(l)
  sum(sapply(Estaciones, function(e) Suavizar_ICA(e, Años_estacion[[e]], 10^l)$gcv)))

pdf("Fig2_VCG_lambda.pdf", width = 6, height = 4.5)
plot(log_lambda, VCG, type = "b", pch = 19, xlab = expression(log[10](lambda)), ylab = "VCG total",
     main = "Selección del parámetro de suavizado")
dev.off()

lambda_opt <- 10^log_lambda[which.min(VCG)]
lambda_opt


#Malla de evaluación común: días 1 a 365
t_eval <- 1:365
M      <- length(t_eval)

#Matriz de curvas (365 x n) y etiqueta de grupo de cada curva
X     <- NULL
Grupo <- NULL
for(e in Estaciones){
  fd_e  <- Suavizar_ICA(e, Años_estacion[[e]], lambda_opt)$fd
  X_e   <- eval.fd(t_eval, fd_e)
  colnames(X_e) <- paste(e, Años_estacion[[e]], sep = "_")
  X     <- cbind(X, X_e)
  Grupo <- c(Grupo, rep(e, ncol(X_e)))
}
Grupo <- factor(Grupo, levels = Estaciones)

n_i <- table(Grupo)
n   <- sum(n_i)
n_i
n


#Curvas por estación con su media
pdf("Fig3_Curvas_por_Estacion.pdf", width = 11, height = 6.5)
par(mfrow = c(2, ceiling(k/2)))
for(e in Estaciones){
  matplot(t_eval, X[, Grupo == e], type = "l", lty = 1, col = "grey60", ylim = c(0, 110),
          xlab = "Día del año", ylab = "ICA", main = paste0(e, " (n = ", n_i[e], ")"))
  lines(t_eval, rowMeans(X[, Grupo == e]), lwd = 2, col = "red")
}
dev.off()



#Promedio y vrianza funcional por grupos

#Matriz indicadora de grupos (n x k) y medias de grupo (365 x k)
Z_grupo <- model.matrix(~ Grupo - 1)
Medias  <- X %*% Z_grupo %*% diag(1/as.numeric(n_i))
colnames(Medias) <- Estaciones
Media_global <- rowMeans(X)

#Varianza funcional dentro de cada grupo (para revisar el supuesto de homocedasticidad)
Varianzas <- sapply(Estaciones, function(e) apply(X[, Grupo == e], 1, var))

colores <- c("red", "blue", "darkgreen", "orange", "purple", "brown")[1:k]

pdf("Fig4_Medias_Varianzas.pdf", width = 11, height = 4.5)
par(mfrow = c(1, 2))
matplot(t_eval, Medias, type = "l", lty = 1, lwd = 2, col = colores, ylim = c(15, 95),
        xlab = "Día del año", ylab = "ICA", main = "Medias funcionales por estación")
lines(t_eval, Media_global, lwd = 2, lty = 2)
legend("top", c(Estaciones, "Global"), col = c(colores, "black"), lty = c(rep(1, k), 2), lwd = 2,
       bty = "n", ncol = 4, cex = 0.75)
matplot(t_eval, sqrt(Varianzas), type = "l", lty = 1, lwd = 2, col = colores,
        ylim = c(0, max(sqrt(Varianzas))*1.3),
        xlab = "Día del año", ylab = "ICA", main = "Desviación estándar funcional por estación")
legend("top", Estaciones, col = colores, lwd = 2, bty = "n", ncol = 3, cex = 0.75)
dev.off()



#Descomposicion de la varianza y F funcional

#Sumas de cuadrados puntuales
SSR_t <- as.numeric((Medias - Media_global)^2 %*% as.numeric(n_i))
SSE_t <- rowSums((X - Medias[, as.integer(Grupo)])^2)
SST_t <- rowSums((X - Media_global)^2)

#Comprobación: SST(t) = SSR(t) + SSE(t)
max(abs(SST_t - (SSR_t + SSE_t)))

#Estadístico F puntual: F(t) = [SSR(t)/(k-1)] / [SSE(t)/(n-k)]
F_t <- (SSR_t/(k - 1))/(SSE_t/(n - k))

#Bajo normalidad y H0, F(t) ~ F_{k-1, n-k} en cada día fijo
p_F_t <- pf(F_t, k - 1, n - k, lower.tail = FALSE)

#Proporción de días con valor p < 0.05 (sin corrección por multiplicidad)
mean(p_F_t < 0.05)

pdf("Fig5_F_Puntual.pdf", width = 11, height = 4.5)
par(mfrow = c(1, 2))
matplot(t_eval, cbind(SSR_t, SSE_t)/1000, type = "l", lty = 1, lwd = 2, col = c("red", "blue"),
        xlab = "Día del año", ylab = "Suma de cuadrados (miles)", main = "SSR(t) y SSE(t)")
legend("topright", c("SSR(t): entre estaciones", "SSE(t): dentro de estaciones"),
       col = c("red", "blue"), lwd = 2, bty = "n", cex = 0.8)
plot(t_eval, F_t, type = "l", lwd = 2, xlab = "Día del año", ylab = "F(t)",
     main = "Prueba F puntual")
abline(h = qf(0.95, k - 1, n - k), col = "red", lty = 2)
dev.off()



#Pruebas globales: normal L2 y tipo F

#Pesos de la regla del trapecio para aproximar las integrales sobre la malla
w <- rep(1, M); w[c(1, M)] <- 0.5

#Sumas de cuadrados integradas
SSR_n <- sum(w*SSR_t)
SSE_n <- sum(w*SSE_t)

#Covarianza combinada (N = n - k grados de libertad) y sus trazas
N_gl   <- n - k
Resid  <- X - Medias[, as.integer(Grupo)]
G_comb <- Resid %*% t(Resid)/N_gl

DG        <- G_comb*w                                #diag(w) %*% G
tr_gamma  <- sum(diag(DG))
tr_gamma2 <- sum(DG*t(DG))

#Comprobación: SSE_n = (n - k) tr(gamma_hat)
c(SSE_n, N_gl*tr_gamma)

#Estimadores con corrección de sesgo (Zhang, 2013)
tr_gamma2_c <- N_gl^2/((N_gl - 1)*(N_gl + 2))*(tr_gamma2 - tr_gamma^2/N_gl)
tr2_gamma_c <- N_gl*(N_gl + 1)/((N_gl - 1)*(N_gl + 2))*(tr_gamma^2 - 2*tr_gamma2/(N_gl + 1))

beta_hat <- tr_gamma2_c/tr_gamma
d_hat    <- tr2_gamma_c/tr_gamma2_c
c(beta = beta_hat, d = d_hat)


#Prueba L2: SSR_n ~ beta chi2_{(k-1) d}
p_L2 <- pchisq(SSR_n/beta_hat, df = (k - 1)*d_hat, lower.tail = FALSE)

#Prueba de tipo F: F_n ~ F_{(k-1) d, (n-k) d}
F_n <- (SSR_n/(k - 1))/(SSE_n/(n - k))
p_F <- pf(F_n, (k - 1)*d_hat, (n - k)*d_hat, lower.tail = FALSE)

c(SSR_n = SSR_n, p_L2 = p_L2, F_n = F_n, p_F = p_F)


#Tabla de ANOVA funcional
Tabla_ANOVA <- data.frame(
  Fuente         = c("Entre estaciones", "Dentro de estaciones", "Total"),
  Suma_Cuadrados = c(SSR_n, SSE_n, SSR_n + SSE_n),
  gl_aprox       = c((k - 1)*d_hat, (n - k)*d_hat, (n - 1)*d_hat),
  Cuadrado_Medio = c(SSR_n/(k - 1), SSE_n/(n - k), NA),
  F              = c(F_n, NA, NA),
  Valor_p        = c(p_F, NA, NA))
Tabla_ANOVA



#Prueba F puntual globalizada 

#Estadístico: T_GPF = integral de F(t) sobre el dominio
T_GPF <- sum(w*F_t)

#La distribución nula depende de la función de correlación rho(s,t)
s_t   <- sqrt(diag(G_comb))
Rho   <- G_comb/(s_t %o% s_t)
DR    <- Rho*w
tr_rho  <- sum(diag(DR))                 #igual a la longitud del dominio (b - a)
tr_rho2 <- sum(DR*t(DR))

beta_rho <- tr_rho2/((k - 1)*tr_rho)
d_rho    <- (k - 1)*tr_rho^2/tr_rho2

p_GPF <- pchisq(T_GPF/beta_rho, df = d_rho, lower.tail = FALSE)
c(T_GPF = T_GPF, beta = beta_rho, d = d_rho, p_GPF = p_GPF)



#Ajuste por permutación:F_n, GPF Y F_max

#Bajo H0 las etiquetas de estación son intercambiables. SST(t) no cambia al permutar,
#de modo que basta recalcular SSR(t) en cada permutación.
B <- 5000

F_t_nulos   <- matrix(NA, M, B)
SSR_n_nulos <- numeric(B)
for(b in 1:B){
  Zb      <- Z_grupo[sample(n), ]                       #etiquetas permutadas
  Med_b   <- X %*% Zb %*% diag(1/colSums(Zb))
  SSR_b   <- as.numeric((Med_b - Media_global)^2 %*% colSums(Zb))
  F_t_nulos[, b] <- (SSR_b/(k - 1))/((SST_t - SSR_b)/(n - k))
  SSR_n_nulos[b] <- sum(w*SSR_b)
}

#Estadísticos globales bajo H0 a partir de las mismas permutaciones
F_n_nulos   <- (SSR_n_nulos/(k - 1))/((sum(w*SST_t) - SSR_n_nulos)/(n - k))
GPF_nulos   <- colSums(w*F_t_nulos)
Fmax_nulos  <- apply(F_t_nulos, 2, max)

#Valores p por permutación (se suma 1 por la asignación observada)
p_perm <- function(nulos, obs) (1 + sum(nulos >= obs - 1e-12))/(B + 1)

F_max   <- max(F_t)
p_F_perm    <- p_perm(F_n_nulos,  F_n)
p_GPF_perm  <- p_perm(GPF_nulos,  T_GPF)
p_Fmax_perm <- p_perm(Fmax_nulos, F_max)

c(p_F_perm = p_F_perm, p_GPF_perm = p_GPF_perm, p_Fmax_perm = p_Fmax_perm)


#Curva crítica puntual y valor crítico máximo al 5 %
alfa      <- 0.05
c_puntual <- apply(F_t_nulos, 1, quantile, probs = 1 - alfa)
c_maximo  <- quantile(Fmax_nulos, 1 - alfa)
p_puntual <- (1 + rowSums(F_t_nulos >= F_t - 1e-12))/(B + 1)


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
Intervalos(F_t > c_maximo)


pdf("Fig6_Permutacion_ANOVA.pdf", width = 11, height = 4.5)
par(mfrow = c(1, 2))
plot(t_eval, F_t, type = "l", lwd = 2, col = "red", xlab = "Día del año", ylab = "F(t)",
     main = "Prueba F puntual por permutación", ylim = c(0, max(F_t, c_maximo)*1.6))
lines(t_eval, c_puntual, col = "blue", lty = 3, lwd = 2)
abline(h = c_maximo, col = "black", lty = 2)
legend("topleft", c("Estadístico observado", "Valor crítico puntual (0.05)", "Valor crítico máximo"),
       col = c("red", "blue", "black"), lty = c(1, 3, 2), lwd = c(2, 2, 1), bty = "n", cex = 0.75)
hist(F_n_nulos, breaks = 50, freq = FALSE, col = "grey85", border = "white",
     xlim = range(0, F_n_nulos, F_n*1.1),
     xlab = expression(F[n]), main = "Distribución de permutación de F")
curve(df(x, (k - 1)*d_hat, (n - k)*d_hat), add = TRUE, col = "blue", lwd = 2)
abline(v = F_n, col = "red", lwd = 2)
legend("topright", c("Aproximación F", "F observado"), col = c("blue", "red"), lwd = 2, bty = "n", cex = 0.8)
dev.off()



#Pruebas de: cuevs, febreros y fraiman

#Estadístico: V_n = suma sobre pares i < l de n_i ||X_i_bar - X_l_bar||^2
pares <- combn(k, 2)
V_n <- sum(apply(pares, 2, function(p)
  n_i[p[1]]*sum(w*(Medias[, p[1]] - Medias[, p[2]])^2)))
V_n

#Bajo H0, V_n se aproxima por V = suma ||Z_i - sqrt(n_i/n_l) Z_l||^2, con Z_i ~ GP(0, gamma_i).
#Cada gamma_i se estima por separado, así que la prueba admite covarianzas distintas.
#Los procesos gaussianos se simulan con la descomposición espectral de cada covarianza.
B_cff <- 2000

Raiz_cov <- lapply(Estaciones, function(e){
  R_e <- X[, Grupo == e] - Medias[, e]
  E_e <- eigen(R_e %*% t(R_e)/(n_i[e] - 1), symmetric = TRUE)
  r   <- n_i[e] - 1                                    #rango de la covarianza muestral
  E_e$vectors[, 1:r] %*% diag(sqrt(pmax(E_e$values[1:r], 0)), r)
})

V_nulos <- numeric(B_cff)
for(b in 1:B_cff){
  Zs <- sapply(Raiz_cov, function(A) A %*% rnorm(ncol(A)))
  V_nulos[b] <- sum(apply(pares, 2, function(p)
    sum(w*(Zs[, p[1]] - sqrt(n_i[p[1]]/n_i[p[2]])*Zs[, p[2]])^2)))
}
p_CFF <- (1 + sum(V_nulos >= V_n))/(B_cff + 1)
p_CFF



#Usando librerias: fdANOVA

#fanova.tests() recibe la matriz de curvas discretizadas (puntos x curvas) y las etiquetas.
#Los estadísticos están en otra escala (no usan los pesos de cuadratura), pero los valores p
#son comparables con los obtenidos a mano.
Pruebas_fdANOVA <- fanova.tests(x = X, group.label = as.character(Grupo),
                                test = c("L2B", "FB", "GPF", "Fmaxb", "CS"),
                                params = list(paramFmaxb = 2000, paramCS = 2000))

Comparacion <- data.frame(
  Prueba    = c("L2 (sesgo corregido)", "Tipo F (sesgo corregido)", "GPF", "F_max", "CFF heterocedástica"),
  p_manual  = c(p_L2, p_F, p_GPF, p_Fmax_perm, p_CFF),
  p_fdANOVA = c(Pruebas_fdANOVA$L2B$pvalueL2B, Pruebas_fdANOVA$FB$pvalueFB,
                Pruebas_fdANOVA$GPF$pvalueGPF, Pruebas_fdANOVA$Fmaxb$pvalueFmaxb,
                Pruebas_fdANOVA$CS$pvalueCS))
Comparacion



#Comparaciones múltiples

#Para cada par de estaciones se aplica la prueba de tipo F de dos muestras con la
#covarianza combinada de TODAS las estaciones (N = n - k) y se ajustan los valores p por Holm.
Pares <- data.frame(t(combn(Estaciones, 2)))
names(Pares) <- c("Estacion_A", "Estacion_B")

Pares$F <- apply(Pares, 1, function(p){
  nA <- n_i[p[1]]; nB <- n_i[p[2]]
  (nA*nB/(nA + nB))*sum(w*(Medias[, p[1]] - Medias[, p[2]])^2)/tr_gamma
})
Pares$p     <- pf(Pares$F, d_hat, N_gl*d_hat, lower.tail = FALSE)
Pares$p_Holm <- p.adjust(Pares$p, method = "holm")
Pares$F      <- round(Pares$F, 3)
Pares$p      <- signif(Pares$p, 3)
Pares$p_Holm <- signif(Pares$p_Holm, 3)
Pares[order(Pares$p), ]


#Matriz de valores p ajustados
Matriz_p <- matrix(NA, k, k, dimnames = list(Estaciones, Estaciones))
for(r in 1:nrow(Pares)){
  Matriz_p[Pares$Estacion_A[r], Pares$Estacion_B[r]] <- Pares$p_Holm[r]
  Matriz_p[Pares$Estacion_B[r], Pares$Estacion_A[r]] <- Pares$p_Holm[r]
}

pdf("Fig7_Comparaciones_Pares.pdf", width = 7, height = 6)
Datos_p <- as.data.frame(as.table(Matriz_p)) %>%
  filter(!is.na(Freq)) %>%
  mutate(Decision = ifelse(Freq < 0.05, "p < 0.05", "p >= 0.05"))

ggplot(Datos_p, aes(Var1, Var2, fill = Decision)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = format(round(Freq, 3), nsmall = 3)), size = 3.5) +
  scale_fill_manual(values = c("p < 0.05" = "cyan", "p >= 0.05" = "grey85"), name = "Holm") +
  labs(x = "", y = "", title = "Valores p ajustados por Holm (prueba F de dos muestras)") +
  theme_minimal()
dev.off()



#Actividad de estudio independiente

#1. Los años incluidos no son los mismos en todas las estaciones: ERA y La Flora incluyen años
#   anteriores a 2014, mientras que las demás empiezan en 2014 o 2015. Repita el ANOVA solo con 2016-2018, años comunes a todas.
#   ¿Cambian las conclusiones? ¿Qué parte de la diferencia podría deberse al periodo y no a la estación?
#2. Compare las desviaciones estándar funcionales de la Figura 4. ¿Es razonable el supuesto de
#   homocedasticidad? ¿Qué prueba es más confiable en este caso?
#3. Con k = 2 estaciones, verifique numéricamente que SSR_n coincide con el estadístico T_n de la
#   prueba t funcional de la práctica anterior.
#4. Agrupe las estaciones por zona de la ciudad (por ejemplo, norte, centro y sur) y realice un ANOVA
#   funcional con esos grupos.
#5. Aumente min_años a 5. ¿Qué estaciones salen del análisis y cómo cambian d_hat y los valores p?

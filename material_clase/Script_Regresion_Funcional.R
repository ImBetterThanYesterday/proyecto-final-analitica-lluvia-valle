#Script - Regresión funcional

#Limpiar área de trabajo
rm(list = ls())

## Librerias
library(fda)

#Carpeta de figuras
dir.create("Figuras_Regresion", showWarnings = FALSE)
figdir <- "Figuras_Regresion"

set.seed(2026)



#Datos de temperatura

#Base de Fourier con 65 funciones sobre [0, 365] y operador de aceleración armónica,
#L = (2*pi/365)^2 D + D^3, que anula las funciones seno y coseno de periodo anual
daybasis65   <- create.fourier.basis(c(0, 365), 65)
harmaccelLfd <- vec2Lfd(c(0, (2*pi/365)^2, 0), c(0, 365))

#Curvas de temperatura (35 estaciones)
tempSmooth <- smooth.basis(day.5, CanadianWeather$dailyAv[, , "Temperature.C"], daybasis65)
tempfd     <- tempSmooth$fd

#Respuesta escalar: log10 de la precipitación anual total
annualprec <- log10(apply(CanadianWeather$dailyAv[, , "Precipitation.mm"], 2, sum))

#Región climática de cada estación
region  <- CanadianWeather$region
regiones <- c("Atlantic", "Pacific", "Continental", "Arctic")
col_reg  <- c(Atlantic = "red", Pacific = "blue", Continental = "darkgreen", Arctic = "purple")

pdf(file.path(figdir, "Fig_Datos_Canada.pdf"), width = 10, height = 4.5)
par(mfrow = c(1, 2))
plot(tempfd, col = col_reg[region], lty = 1, xlab = "Día", ylab = "Temperatura (°C)",
     main = "Curvas de temperatura")
legend("bottom", regiones, col = col_reg[regiones], lwd = 2, bty = "n", ncol = 2, cex = 0.8)
temp_media <- apply(CanadianWeather$dailyAv[, , "Temperature.C"], 2, mean)
plot(temp_media, annualprec, pch = 19, col = col_reg[region], xlab = "Temperatura media anual (°C)",
     ylab = expression(log[10]~"precipitación anual"), main = "Respuesta escalar")
dev.off()



#Modelo escalar-funcional

#Covariables: una constante (intercepto) y la curva de temperatura
templist <- list(const = rep(1, 35), temp = tempfd)
conbasis <- create.constant.basis(c(0, 365))


#beta(t) en una base de dimensión baja (Fourier con 5 funciones, sin penalización)
betabasis5 <- create.fourier.basis(c(0, 365), 5)
betalist5  <- list(const = conbasis, temp = betabasis5)
fit5       <- fRegress(annualprec, templist, betalist5)

yhat5 <- as.numeric(fit5$yhatfdobj)
R2_5  <- 1 - sum((annualprec - yhat5)^2)/sum((annualprec - mean(annualprec))^2)
R2_5


#base amplia (Fourier con 35 funciones) + penalización de rugosidad.
#lambda se elige por validación cruzada ordinaria (OCV), que fRegress entrega directamente.
betabasis35 <- create.fourier.basis(c(0, 365), 35)
loglam      <- seq(5, 15, by = 0.5)
OCV         <- sapply(loglam, function(l){
  bl <- list(const = conbasis, temp = fdPar(betabasis35, harmaccelLfd, 10^l))
  fRegress(annualprec, templist, bl)$OCV
})

loglam_opt <- loglam[which.min(OCV)]
loglam_opt

betalist  <- list(const = conbasis, temp = fdPar(betabasis35, harmaccelLfd, 10^loglam_opt))
fitPen    <- fRegress(annualprec, templist, betalist)

yhatPen <- as.numeric(fitPen$yhatfdobj)
resPen  <- annualprec - yhatPen
R2_Pen  <- 1 - sum(resPen^2)/sum((annualprec - mean(annualprec))^2)
df_Pen  <- fitPen$df
c(R2 = R2_Pen, df = df_Pen)

#Errores estándar puntuales (estimador sándwich): sigma2 = SSE / (n - df).
#En fda >= 6, fRegress() los calcula cuando recibe SigmaE y la matriz y2cMap
#(para respuesta escalar basta la identidad, que mapea y en sí misma).
sigma2  <- sum(resPen^2)/(35 - df_Pen)
SigmaE  <- sigma2*diag(35)
fitPen  <- fRegress(annualprec, templist, betalist, y2cMap = diag(35), SigmaE = SigmaE)

betaPen   <- fitPen$betaestlist$temp$fd
betaSE    <- fitPen$betastderrlist[[2]]


#Regresión sobre componentes principales funcionales (fPCR)
K       <- 4
#Las funciones propias se suavizan con la misma penalización de aceleración armónica
tempPCA <- pca.fd(tempfd, nharm = K, harmfdPar = fdPar(daybasis65, harmaccelLfd, 1e5))
round(cumsum(tempPCA$varprop), 3)

Scores  <- tempPCA$scores
fitPCR  <- lm(annualprec ~ Scores)
summary(fitPCR)$r.squared

#beta(t) = suma_j b_j xi_j(t)
betaPCR <- tempPCA$harmonics[1]*coef(fitPCR)[2]
for(j in 2:K) betaPCR <- betaPCR + tempPCA$harmonics[j]*coef(fitPCR)[j + 1]

#Banda de 2 errores estándar para beta(t) en fPCR
Vb       <- vcov(fitPCR)[-1, -1]
Xi_t     <- eval.fd(day.5, tempPCA$harmonics)
se_PCR   <- sqrt(rowSums((Xi_t %*% Vb)*Xi_t))


#rueba de permutación para H0: beta(t) = 0
FpermEsc <- Fperm.fd(annualprec, templist, betalist, nperm = 500, plotres = FALSE)
FpermEsc$pval


#Figuras del modelo escalar-funcional
pdf(file.path(figdir, "Fig_EscFun_OCV.pdf"), width = 10, height = 4.5)
par(mfrow = c(1, 2))
plot(loglam, OCV, type = "b", pch = 19, xlab = expression(log[10](lambda)), ylab = "OCV",
     main = "Selección de lambda por validación cruzada")
abline(v = loglam_opt, lty = 2, col = "red")
plot(yhatPen, annualprec, pch = 19, col = col_reg[region], xlab = "Valor ajustado", ylab = "Valor observado",
     main = bquote("Observado vs. ajustado (" * R^2 == .(round(R2_Pen, 3)) * ")"))
abline(0, 1, lty = 2)
dev.off()

pdf(file.path(figdir, "Fig_EscFun_Beta.pdf"), width = 10, height = 4.5)
par(mfrow = c(1, 2))
b_pen  <- eval.fd(day.5, betaPen)
se_pen <- eval.fd(day.5, betaSE)
matplot(day.5, cbind(b_pen, b_pen + 2*se_pen, b_pen - 2*se_pen), type = "l", lty = c(1, 2, 2),
        col = c("black", "red", "red"), lwd = c(2, 1, 1), xlab = "Día", ylab = expression(beta(t)),
        main = "Penalización de rugosidad")
abline(h = 0, lty = 3)
b_pcr <- eval.fd(day.5, betaPCR)
b_5   <- eval.fd(day.5, fit5$betaestlist$temp$fd)
matplot(day.5, cbind(b_pcr, b_pcr + 2*se_PCR, b_pcr - 2*se_PCR, b_5), type = "l", lty = c(1, 2, 2, 1),
        col = c("black", "red", "red", "blue"), lwd = c(2, 1, 1, 2), xlab = "Día", ylab = expression(beta(t)),
        main = paste0("fPCR (K = ", K, ") y base de Fourier con 5 funciones"))
abline(h = 0, lty = 3)
legend("topleft", c("fPCR", "Fourier, K = 5"), col = c("black", "blue"), lwd = 2, bty = "n", cex = 0.8)
dev.off()



#Modelo funcional-escalar

#ANOVA funcional de la temperatura por región, escrito como regresión.
#Restricción de identificación: la suma de los efectos regionales es 0 en todo t.
#Se impone agregando una observación ficticia con respuesta 0 (Ramsay y Silverman, 2005, cap. 13).
p <- length(regiones) + 1
regionList <- vector("list", p)
names(regionList) <- c("Canada", regiones)
regionList[[1]] <- c(rep(1, 35), 0)
for(j in 2:p) regionList[[j]] <- c(as.numeric(region == regiones[j - 1]), 1)

coef36    <- cbind(tempfd$coefs, matrix(0, 65, 1))
temp36fd  <- fd(coef36, daybasis65, tempfd$fdnames)

betabasis11 <- create.fourier.basis(c(0, 365), 11)
betaListReg <- vector("list", p)
for(j in 1:p) betaListReg[[j]] <- fdPar(betabasis11)

fitReg <- fRegress(temp36fd, regionList, betaListReg)

#Errores estándar puntuales a partir de la covarianza de los residuos
yhatReg  <- eval.fd(day.5, fitReg$yhatfdobj)[, 1:35]
tempmat  <- eval.fd(day.5, tempfd)
residReg <- tempmat - yhatReg
SigmaReg <- var(t(residReg))

stderrReg <- fRegress(temp36fd, regionList, betaListReg,
                      y2cMap = tempSmooth$y2cMap, SigmaE = SigmaReg)

#Proporción de varianza explicada en cada día: R2(t)
R2_t <- 1 - rowSums(residReg^2)/rowSums((tempmat - rowMeans(tempmat))^2)
range(R2_t)

#Prueba de permutación para H0: todos los efectos regionales son nulos
FpermReg <- Fperm.fd(temp36fd, regionList, betaListReg, nperm = 200, plotres = FALSE)
FpermReg$pval

pdf(file.path(figdir, "Fig_FunEsc_Region.pdf"), width = 11, height = 6.5)
par(mfrow = c(2, 3))
for(j in 1:p){
  bj  <- eval.fd(day.5, fitReg$betaestlist[[j]]$fd)
  sej <- eval.fd(day.5, stderrReg$betastderrlist[[j]])
  matplot(day.5, cbind(bj, bj + 2*sej, bj - 2*sej), type = "l", lty = c(1, 2, 2), lwd = c(2, 1, 1),
          col = c("black", "red", "red"), xlab = "Día", ylab = "°C",
          main = ifelse(j == 1, "Media de Canadá: beta_0(t)", paste("Efecto", names(regionList)[j])))
  abline(h = 0, lty = 3)
}
plot(day.5, R2_t, type = "l", lwd = 2, ylim = c(0, 1), xlab = "Día", ylab = expression(R^2*(t)),
     main = "Varianza explicada por región")
dev.off()



#Modelo concurrente

#Ángulos de la cadera y la rodilla de 39 niños durante un ciclo de marcha (20 instantes)
gaittime  <- as.numeric(dimnames(gait)[[1]])*20
gaitrange <- c(0, 20)
gaitbasis <- create.fourier.basis(gaitrange, nbasis = 21)
harmaccelLfd20 <- vec2Lfd(c(0, (2*pi/20)^2, 0), rangeval = gaitrange)

gaitSmooth <- smooth.basis(gaittime, gait, fdPar(gaitbasis, harmaccelLfd20, 1e-2))
gaitfd     <- gaitSmooth$fd
hipfd  <- gaitfd[, 1]
kneefd <- gaitfd[, 2]

xfdlistGait <- list(const = rep(1, 39), hip = hipfd)
betafdPar   <- fdPar(gaitbasis, harmaccelLfd20)
betalistGait <- list(const = betafdPar, hip = betafdPar)

fitGait <- fRegress(kneefd, xfdlistGait, betalistGait)

#Errores estándar puntuales
gaitfine  <- seq(0, 20, length.out = 101)
kneemat   <- eval.fd(gaittime, kneefd)
kneehat   <- eval.fd(gaittime, fitGait$yhatfdobj)
residGait <- kneemat - kneehat
SigmaGait <- var(t(residGait))
y2cMapGait <- smooth.basis(gaittime, gait[, , 2], fdPar(gaitbasis, harmaccelLfd20, 1e-2))$y2cMap
stderrGait <- fRegress(kneefd, xfdlistGait, betalistGait, y2cMap = y2cMapGait, SigmaE = SigmaGait)

#R2(t)
kmat  <- eval.fd(gaitfine, kneefd)
khat  <- eval.fd(gaitfine, fitGait$yhatfdobj)
R2_gait <- 1 - rowSums((kmat - khat)^2)/rowSums((kmat - rowMeans(kmat))^2)

pdf(file.path(figdir, "Fig_Concurrente_Gait.pdf"), width = 11, height = 6.5)
par(mfrow = c(2, 3))
plot(hipfd, col = "grey50", lty = 1, xlab = "Tiempo (1/20 del ciclo)", ylab = "Grados", main = "Cadera: x_i(t)")
plot(kneefd, col = "grey50", lty = 1, xlab = "Tiempo (1/20 del ciclo)", ylab = "Grados", main = "Rodilla: y_i(t)")
plot(gaitfine, R2_gait, type = "l", lwd = 2, ylim = c(0, 1), xlab = "Tiempo (1/20 del ciclo)",
     ylab = expression(R^2*(t)), main = "Varianza explicada por la cadera")
for(j in 1:2){
  bj  <- eval.fd(gaitfine, fitGait$betaestlist[[j]]$fd)
  sej <- eval.fd(gaitfine, stderrGait$betastderrlist[[j]])
  matplot(gaitfine, cbind(bj, bj + 2*sej, bj - 2*sej), type = "l", lty = c(1, 2, 2), lwd = c(2, 1, 1),
          col = c("black", "red", "red"), xlab = "Tiempo (1/20 del ciclo)", ylab = "",
          main = ifelse(j == 1, "Intercepto: beta_0(t)", "Efecto de la cadera: beta_1(t)"))
  abline(h = 0, lty = 3)
}
#Correlación de los residuos entre instantes del ciclo (Hooker, 2010)
resfine <- kmat - khat
contour(gaitfine, gaitfine, cor(t(resfine)), xlab = "Tiempo (s)", ylab = "Tiempo (t)",
        main = "Correlación de los residuos", nlevels = 10)
dev.off()

range(R2_gait)



#Modelo funcional-funcional

#Respuesta funcional: log10 de la precipitación diaria (se suma 0.05 para evitar log(0))
precmat   <- CanadianWeather$dailyAv[, , "Precipitation.mm"]
logprec   <- log10(precmat + 0.05)
lnprecfd  <- smooth.basis(day.5, logprec, fdPar(daybasis65, harmaccelLfd, 1e5))$fd

#Bases para beta_0(t) y beta_1(s,t)
betabasis_s <- create.fourier.basis(c(0, 365), 11)
betabasis_t <- create.fourier.basis(c(0, 365), 11)
beta0Par    <- fdPar(create.fourier.basis(c(0, 365), 11), harmaccelLfd, 0)

Ajuste_linmod <- function(lam, idx = 1:35){
  beta1bifd <- bifd(matrix(0, 11, 11), betabasis_s, betabasis_t)
  beta1Par  <- bifdPar(beta1bifd, harmaccelLfd, harmaccelLfd, lam, lam)
  linmod(tempfd[idx], lnprecfd[idx], list(beta0Par, beta1Par))
}

#Elección de lambda (común en s y t) por validación cruzada dejando una curva fuera
loglamFF <- seq(3, 9, by = 1)
CV_FF <- sapply(loglamFF, function(l){
  e <- 0
  for(i in 1:35){
    fit_i <- Ajuste_linmod(10^l, idx = (1:35)[-i])
    b0  <- eval.fd(day.5, fit_i$beta0estfd)
    B1  <- eval.bifd(day.5, day.5, fit_i$beta1estbifd)
    xi  <- eval.fd(day.5, tempfd[i])
    yi_hat <- b0 + t(B1) %*% xi           #integral sobre s aproximada con paso de 1 día
    yi     <- eval.fd(day.5, lnprecfd[i])
    e <- e + sum((yi - yi_hat)^2)
  }
  e
})
loglamFF_opt <- loglamFF[which.min(CV_FF)]
cbind(loglamFF, CV_FF)

fitFF <- Ajuste_linmod(10^loglamFF_opt)

#R2(t)
ymat    <- eval.fd(day.5, lnprecfd)
yhatFF  <- eval.fd(day.5, fitFF$yhatfdobj)
R2_FF_t <- 1 - rowSums((ymat - yhatFF)^2)/rowSums((ymat - rowMeans(ymat))^2)
mean(R2_FF_t)

semanas <- seq(0, 365, length.out = 53)
B1mat   <- eval.bifd(semanas, semanas, fitFF$beta1estbifd)

pdf(file.path(figdir, "Fig_FunFun_Superficie.pdf"), width = 11, height = 5)
par(mfrow = c(1, 2))
persp(semanas, semanas, B1mat, theta = 30, phi = 30, col = "lightblue", border = "grey40",
      xlab = "s (temperatura)", ylab = "t (precipitación)", zlab = "beta_1(s,t)",
      main = expression(hat(beta)[1](s, t)), ticktype = "simple")
contour(semanas, semanas, B1mat, xlab = "s: día de la temperatura", ylab = "t: día de la precipitación",
        main = "Curvas de nivel", nlevels = 12)
abline(0, 1, lty = 2, col = "red")
dev.off()

pdf(file.path(figdir, "Fig_FunFun_Ajuste.pdf"), width = 11, height = 4.5)
par(mfrow = c(1, 3))
for(est in c("Montreal", "Vancouver", "Resolute")){
  i <- which(CanadianWeather$place == est)
  matplot(day.5, cbind(logprec[, i], ymat[, i], yhatFF[, i]), type = c("p", "l", "l"), pch = 20, cex = 0.4,
          col = c("grey60", "black", "red"), lty = c(1, 1, 2), lwd = c(1, 2, 2),
          xlab = "Día", ylab = expression(log[10]~"precipitación"), main = est)
}
legend("bottomright", c("Datos", "Curva suavizada", "Ajuste del modelo"), col = c("grey60", "black", "red"),
       pch = c(20, NA, NA), lty = c(NA, 1, 2), lwd = 2, bty = "n", cex = 0.8)
dev.off()



#Ejercicios de estudio independiente

#1. En el modelo escalar-funcional, cambie el operador de aceleración armónica por int2Lfd(2).
#   ¿Cambian la forma de beta(t) y el lambda óptimo?
#2. Repita la fPCR con K = 2, 3, ..., 8 y grafique el R2 y el error de validación cruzada en función de K.
#3. Con los datos de precipitación del curso, ajuste el modelo concurrente Y_est(t) = b0(t) + b1(t) X_sat(t).
#   ¿En qué parte del año el satélite explica mejor las mediciones en tierra?
#4. Con las curvas anuales del ICA, ajuste el modelo funcional-escalar ICA(t) ~ estación y compare
#   los efectos estimados con los resultados del ANOVA funcional.
#5. En el modelo funcional-funcional, use lambdas distintos en s y en t. ¿Qué dirección requiere más suavizado?

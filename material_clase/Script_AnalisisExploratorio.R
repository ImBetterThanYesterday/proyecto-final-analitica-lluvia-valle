#Script: Análisis exploratorio de datos funcionales
#Fecha: Octubre 19 de 2026

#Limpiar área de trabajo
rm(list=ls())

#Instalar paquetes
install.packages(c("fda", "fda.usc"))

#Cargar librerías
library(fda)
library(fda.usc)

#Directorio de trabajo

#Nota cambiar el directorio de trabajo
setwd("C:/Users/Johann/Dropbox/Curso_AnaliticaDatos/AnalisisExploratorio")

#Verificación
getwd()
dir()

#Carpeta para guardar las figuras
dir.create("Figuras_Sesion3", showWarnings = FALSE)


#Cargar datos
SO2_Trihorario <- read.csv("SO2_Trihorario_19Mayo2003_ERA.csv")
X_sat          <- read.csv("X_sat.csv")   #Precipitación acumulada de satélite
Y_est          <- read.csv("Y_est.csv")   #Precipitación acumulada de estaciones en tierra

head(SO2_Trihorario)
tail(SO2_Trihorario)
dim(X_sat)
dim(Y_est)


#Poner nombres a las filas
row.names(SO2_Trihorario) <- SO2_Trihorario[,1]
row.names(X_sat) <- paste("Sitio", 1:nrow(X_sat), sep = "")
row.names(Y_est) <- paste("Sitio", 1:nrow(Y_est), sep = "")

#Quitar la primera columna
SO2_Trihorario <- SO2_Trihorario[,-1]


#Transponer los conjuntos de datos (tiempos en las filas y curvas en las columnas)
SO2_Trihorario_t <- t(as.matrix(SO2_Trihorario))   #24 x 8
X_sat_t          <- t(as.matrix(X_sat))            #72 x 10
Y_est_t          <- t(as.matrix(Y_est))            #72 x 10

TiempoEval_SO2  = seq(1, dim(SO2_Trihorario_t)[1])
rangval_SO2     = c(1, dim(SO2_Trihorario_t)[1])

TiempoEval_Prec = seq(1, dim(X_sat_t)[1])
rangval_Prec    = c(1, dim(X_sat_t)[1])

#n = 24: Número de puntos muestreados del SO2
#N = 8:  Número de días
#n = 72: Número de pentadas de precipitación
#N = 10: Número de sitios (la fila i de X_sat y de Y_est es el mismo sitio)


#Gráficar los datos observados
windows()
matplot(SO2_Trihorario_t, type = "o", pch = 1,
        xlab = "Tiempo (horas)",
        ylab = expression(SO[2]), main = "Datos observados")

windows()
par(mfrow = c(1,2))
matplot(X_sat_t, type = "o", pch = 1, cex = 0.4, ylim = c(0, 2500),
        xlab = "Tiempo (pentadas)", ylab = "Precipitación (mm)", main = "Satélite")
matplot(Y_est_t, type = "o", pch = 1, cex = 0.4, ylim = c(0, 2500),
        xlab = "Tiempo (pentadas)", ylab = "Precipitación (mm)", main = "Estaciones")

pdf("Figuras_Sesion3/Fig_Precip_Observados.pdf", width = 10, height = 4.5)
par(mfrow = c(1,2), mar = c(4.2,4.2,2,1))
matplot(X_sat_t, type = "o", pch = 1, cex = 0.4, ylim = c(0, 2500),
        xlab = "Tiempo (pentadas)", ylab = "Precipitación (mm)", main = "Satélite")
matplot(Y_est_t, type = "o", pch = 1, cex = 0.4, ylim = c(0, 2500),
        xlab = "Tiempo (pentadas)", ylab = "Precipitación (mm)", main = "Estaciones")
dev.off()






#Construcción de datos funcionales usando las funciones de la libreria "fda"

#SO2: B-splines de orden 2 con K = 14 (valor elegido por la VCG)
BasesSO2   = create.bspline.basis(rangval_SO2, nbasis = 14, norder = 2)
fd_nombres_SO2 = list("Tiempo (horas)", "Día", "SO2")

SO2_Trihorario_smoothList <- smooth.basis(TiempoEval_SO2,
                                          SO2_Trihorario_t,
                                          BasesSO2,
                                          fdnames = fd_nombres_SO2)

#Precipitación: B-splines cúbicos (orden 4) con K = 20 y sin penalización
BasesPrec  = create.bspline.basis(rangval_Prec, nbasis = 20, norder = 4)
fd_nombres_Prec = list("Tiempo (pentadas)", "Sitio", "Precipitación (mm)")

X_sat_smoothList <- smooth.basis(TiempoEval_Prec,
                                 X_sat_t,
                                 BasesPrec,
                                 fdnames = fd_nombres_Prec)

Y_est_smoothList <- smooth.basis(TiempoEval_Prec,
                                 Y_est_t,
                                 BasesPrec,
                                 fdnames = fd_nombres_Prec)

names(SO2_Trihorario_smoothList)

#Extraer los datos funcionales
SO2_Trihorario_fd <- SO2_Trihorario_smoothList$fd
X_sat_fd          <- X_sat_smoothList$fd
Y_est_fd          <- Y_est_smoothList$fd

#Validación cruzada generalizada de cada ajuste
sum(SO2_Trihorario_smoothList$gcv)
sum(X_sat_smoothList$gcv)
sum(Y_est_smoothList$gcv)

#Malla fina para evaluar las funciones
tSO2  <- seq(rangval_SO2[1],  rangval_SO2[2],  length.out = 101)
tPrec <- seq(rangval_Prec[1], rangval_Prec[2], length.out = 144)

SO2_fvalores <- eval.fd(tSO2,  SO2_Trihorario_fd)   #101 x 8
X_fvalores   <- eval.fd(tPrec, X_sat_fd)            #144 x 10
Y_fvalores   <- eval.fd(tPrec, Y_est_fd)            #144 x 10

windows()
par(mfrow = c(1,2))
plot(X_sat_fd, ylim = c(0, 2500), xlab = "Tiempo (pentadas)", ylab = "Precipitación (mm)")
title("Satélite")
plot(Y_est_fd, ylim = c(0, 2500), xlab = "Tiempo (pentadas)", ylab = "Precipitación (mm)")
title("Estaciones")

pdf("Figuras_Sesion3/Fig_Precip_Suavizados.pdf", width = 10, height = 4.5)
par(mfrow = c(1,2), mar = c(4.2,4.2,2,1))
matplot(tPrec, X_fvalores, type = "l", lty = 1, ylim = c(0, 2500),
        xlab = "Tiempo (pentadas)", ylab = "Precipitación (mm)", main = "Satélite")
matplot(tPrec, Y_fvalores, type = "l", lty = 1, ylim = c(0, 2500),
        xlab = "Tiempo (pentadas)", ylab = "Precipitación (mm)", main = "Estaciones")
dev.off()






#Media funcional paso a paso: x_barra(t) = c_barra' phi(t)
CoefSO2     = SO2_Trihorario_fd$coefs        #Matriz de coeficientes de 14x8
CoefMedia   = rowMeans(CoefSO2)              #Vector de coeficientes medios
bvals       = eval.basis(tSO2, BasesSO2)     #Matriz phi
dim(bvals)

MediaSO2_paso = bvals%*%CoefMedia

#Media funcional con la función mean.fd()
MediaSO2_fd = mean.fd(SO2_Trihorario_fd)
MediaSO2    = eval.fd(tSO2, MediaSO2_fd)

#Comprobación: la diferencia debe ser cero
max(abs(MediaSO2_paso - MediaSO2))

MediaX_fd = mean.fd(X_sat_fd)
MediaY_fd = mean.fd(Y_est_fd)
MediaX    = eval.fd(tPrec, MediaX_fd)
MediaY    = eval.fd(tPrec, MediaY_fd)






#Varianza y desviación estándar funcional
#Var(t) = (N-1)^{-1} sum_i [x_i(t) - x_barra(t)]^2

N_SO2     = ncol(SO2_fvalores)
VarSO2    = apply((SO2_fvalores - as.vector(MediaSO2))^2, 1, sum)/(N_SO2 - 1)
DESO2     = sqrt(VarSO2)

#Con la función std.fd() se obtiene una versión suavizada de la desviación estándar
DESO2_fd  = std.fd(SO2_Trihorario_fd)
summary(as.vector(eval.fd(tSO2, DESO2_fd)) - DESO2)

DEX = sqrt(apply((X_fvalores - as.vector(MediaX))^2, 1, sum)/(ncol(X_fvalores) - 1))
DEY = sqrt(apply((Y_fvalores - as.vector(MediaY))^2, 1, sum)/(ncol(Y_fvalores) - 1))

#Coeficiente de variación funcional
CVX = DEX/MediaX
CVY = DEY/MediaY

#Gráfico de la media +- 2 desviaciones estándar para el SO2
windows()
matplot(tSO2, SO2_fvalores, type = "l", lty = 1, col = "grey70",
        xlab = "Tiempo (horas)", ylab = expression(SO[2]), ylim = c(-2, 20))
lines(tSO2, MediaSO2, lwd = 3, col = "blue")
lines(tSO2, MediaSO2 + 2*DESO2, lwd = 2, lty = 2, col = "red")
lines(tSO2, MediaSO2 - 2*DESO2, lwd = 2, lty = 2, col = "red")

pdf("Figuras_Sesion3/Fig_SO2_MediaDE.pdf", width = 10, height = 4.5)
par(mfrow = c(1,2), mar = c(4.2,4.2,2,1))
matplot(tSO2, SO2_fvalores, type = "l", lty = 1, col = "grey65",
        xlab = "Tiempo (horas)", ylab = expression(SO[2]), ylim = c(-2, 20),
        main = "Media y media \u00B1 2 DE")
lines(tSO2, MediaSO2, lwd = 3, col = "blue")
lines(tSO2, MediaSO2 + 2*DESO2, lwd = 2, lty = 2, col = "red")
lines(tSO2, MediaSO2 - 2*DESO2, lwd = 2, lty = 2, col = "red")
legend("topright", c("Curvas", "Media", "Media \u00B1 2 DE"),
       col = c("grey65", "blue", "red"), lty = c(1,1,2), lwd = c(1,3,2), bty = "n")
plot(tSO2, DESO2, type = "l", lwd = 3, col = "darkgreen",
     xlab = "Tiempo (horas)", ylab = "Desviación estándar", main = "Desviación estándar funcional")
dev.off()

#Gráfico de la media y la desviación estándar de la precipitación
windows()
par(mfrow = c(1,2))
plot(tPrec, MediaY, type = "l", lwd = 3, col = "slategray", ylim = c(0, 2500),
     xlab = "Tiempo (pentadas)", ylab = "Precipitación (mm)", main = "Media funcional")
lines(tPrec, MediaX, lwd = 3, lty = 2, col = "darkturquoise")
plot(tPrec, DEY, type = "l", lwd = 3, col = "slategray", ylim = c(0, 700),
     xlab = "Tiempo (pentadas)", ylab = "Precipitación (mm)", main = "Desviación estándar funcional")
lines(tPrec, DEX, lwd = 3, lty = 2, col = "darkturquoise")

pdf("Figuras_Sesion3/Fig_Precip_Media.pdf", width = 10, height = 4.5)
par(mfrow = c(1,2), mar = c(4.2,4.2,2,1))
plot(tPrec, MediaY, type = "l", lwd = 3, col = "slategray", ylim = c(0, 2500),
     xlab = "Tiempo (pentadas)", ylab = "Precipitación (mm)", main = "Media funcional")
lines(tPrec, MediaX, lwd = 3, lty = 2, col = "darkturquoise")
legend("topleft", c("Estaciones", "Satélite"), col = c("slategray", "darkturquoise"),
       lty = c(1,2), lwd = 3, bty = "n")
plot(tPrec, DEY, type = "l", lwd = 3, col = "slategray", ylim = c(0, 700),
     xlab = "Tiempo (pentadas)", ylab = "Precipitación (mm)", main = "Desviación estándar funcional")
lines(tPrec, DEX, lwd = 3, lty = 2, col = "darkturquoise")
legend("topleft", c("Estaciones", "Satélite"), col = c("slategray", "darkturquoise"),
       lty = c(1,2), lwd = 3, bty = "n")
dev.off()

#Valores al final del año (pentada 72)
c(MediaX = MediaX[length(tPrec)], MediaY = MediaY[length(tPrec)],
  DEX = DEX[length(tPrec)],       DEY = DEY[length(tPrec)])






#Covarianza funcional paso a paso: Cov(s,t) = phi(s)' Sigma_c phi(t)
CoefCentrados = CoefSO2 - CoefMedia
SigmaC        = CoefCentrados%*%t(CoefCentrados)/(N_SO2 - 1)   #Matriz de 14x14
CovSO2_paso   = bvals%*%SigmaC%*%t(bvals)

#Covarianza funcional con la función var.fd()
CovSO2_bifd = var.fd(SO2_Trihorario_fd)
CovSO2      = eval.bifd(tSO2, tSO2, CovSO2_bifd)

max(abs(CovSO2_paso - CovSO2))

#La diagonal de la covarianza es la varianza funcional
max(abs(diag(CovSO2) - VarSO2))

#Correlación funcional
CorSO2 = cor.fd(tSO2, SO2_Trihorario_fd)

windows()
par(mfrow = c(1,2))
persp(tSO2, tSO2, CovSO2, theta = -45, phi = 25, r = 3, expand = 0.6,
      xlab = "s (horas)", ylab = "t (horas)", zlab = "Covarianza", col = "grey80",
      border = NA, shade = 0.5, ticktype = "detailed")
contour(tSO2, tSO2, CorSO2, xlab = "s (horas)", ylab = "t (horas)", main = "Correlación")

pdf("Figuras_Sesion3/Fig_SO2_CovCor.pdf", width = 10, height = 4.8)
par(mfrow = c(1,2), mar = c(4.2,4.2,2,1))
persp(tSO2, tSO2, CovSO2, theta = -45, phi = 25, r = 3, expand = 0.6,
      xlab = "s (horas)", ylab = "t (horas)", zlab = "Covarianza", col = "grey85",
      border = NA, shade = 0.5, ticktype = "detailed", cex.axis = 0.7, main = "Covarianza")
image(tSO2, tSO2, CorSO2, zlim = c(-1, 1), col = hcl.colors(40, "Blue-Red 3"),
      xlab = "s (horas)", ylab = "t (horas)", main = "Correlación")
contour(tSO2, tSO2, CorSO2, levels = seq(-1, 1, 0.5), labcex = 0.7, add = TRUE)
dev.off()

#Covarianza de la precipitación de estaciones y de satélite
CovY = eval.bifd(tPrec, tPrec, var.fd(Y_est_fd))
CovX = eval.bifd(tPrec, tPrec, var.fd(X_sat_fd))

windows()
par(mfrow = c(1,2))
persp(tPrec, tPrec, CovY, theta = -45, phi = 25, r = 3, expand = 0.6,
      xlab = "s", ylab = "t", zlab = "Covarianza", col = "grey80", border = NA, shade = 0.5)
title("Estaciones")
persp(tPrec, tPrec, CovX, theta = -45, phi = 25, r = 3, expand = 0.6,
      xlab = "s", ylab = "t", zlab = "Covarianza", col = "grey80", border = NA, shade = 0.5)
title("Satélite")

pdf("Figuras_Sesion3/Fig_Precip_Cov.pdf", width = 10, height = 4.8)
par(mfrow = c(1,2), mar = c(4.2,4.2,2,1))
contour(tPrec, tPrec, CovY/1000, nlevels = 10, labcex = 0.7,
        xlab = "s (pentadas)", ylab = "t (pentadas)", main = expression("Estaciones (mm"^2*" x 1000)"))
contour(tPrec, tPrec, CovX/1000, nlevels = 10, labcex = 0.7,
        xlab = "s (pentadas)", ylab = "t (pentadas)", main = expression("Satélite (mm"^2*" x 1000)"))
dev.off()






#Covarianza y correlación cruzada: satélite (X) Vs estaciones (Y)
#Cov_XY(s,t) = (N-1)^{-1} sum_i [x_i(s) - x_barra(s)][y_i(t) - y_barra(t)]

CovCruzada = eval.bifd(tPrec, tPrec, var.fd(X_sat_fd, Y_est_fd))
CorCruzada = cor.fd(tPrec, X_sat_fd, tPrec, Y_est_fd)

#Comprobación paso a paso
X_centrados = X_fvalores - as.vector(MediaX)
Y_centrados = Y_fvalores - as.vector(MediaY)
CovCruzada_paso = X_centrados%*%t(Y_centrados)/(ncol(X_fvalores) - 1)
max(abs(CovCruzada_paso - CovCruzada))

#Correlación en la diagonal: r_XY(t,t)
CorDiag = diag(CorCruzada)
summary(CorDiag)

windows()
persp(tPrec, tPrec, CorCruzada, theta = -45, phi = 25, r = 3, expand = 0.6,
      xlab = "s satélite", ylab = "t estaciones", zlab = "Correlación", col = "grey80",
      border = NA, shade = 0.5)

pdf("Figuras_Sesion3/Fig_Precip_CorCruzada.pdf", width = 10, height = 4.8)
par(mfrow = c(1,2), mar = c(4.2,4.2,2,1))
image(tPrec, tPrec, CorCruzada, zlim = c(-1, 1), col = hcl.colors(40, "Blue-Red 3"),
      xlab = "s: satélite (pentadas)", ylab = "t: estaciones (pentadas)",
      main = "Correlación cruzada")
contour(tPrec, tPrec, CorCruzada, levels = seq(-1, 1, 0.2), labcex = 0.7, add = TRUE)
abline(0, 1, lty = 3, col = "red")
plot(tPrec, CorDiag, type = "l", lwd = 3, col = "darkblue", ylim = c(-1, 1),
     xlab = "Tiempo (pentadas)", ylab = expression(r[XY](t,t)), main = "Diagonal de la correlación cruzada")
abline(h = 0, lty = 3)
dev.off()






#Derivadas de los datos funcionales: D x(t) = c' D phi(t)
#En la precipitación acumulada la primera derivada es la intensidad de lluvia (mm/pentada)

#Exploración de acuerdo al orden de la base: orden 2 Vs orden 4
BasesPrec_o2 = create.bspline.basis(rangval_Prec, nbasis = 20, norder = 2)
Y_est_fd_o2  = smooth.basis(TiempoEval_Prec, Y_est_t, BasesPrec_o2)$fd

DMediaY_o2 = eval.fd(tPrec, deriv.fd(mean.fd(Y_est_fd_o2), 1))
DMediaY_o4 = eval.fd(tPrec, deriv.fd(MediaY_fd, 1))

#La derivada de la media es igual a la media de las derivadas
DY_fd = deriv.fd(Y_est_fd, 1)
DX_fd = deriv.fd(X_sat_fd, 1)
max(abs(eval.fd(tPrec, mean.fd(DY_fd)) - DMediaY_o4))

DMediaX = eval.fd(tPrec, mean.fd(DX_fd))
DMediaY = DMediaY_o4

windows()
par(mfrow = c(1,2))
plot(tPrec, DMediaY_o2, type = "l", lwd = 2, xlab = "Tiempo (pentadas)", ylab = "mm/pentada",
     main = "Orden 2")
plot(tPrec, DMediaY_o4, type = "l", lwd = 2, xlab = "Tiempo (pentadas)", ylab = "mm/pentada",
     main = "Orden 4")

pdf("Figuras_Sesion3/Fig_Precip_Derivadas.pdf", width = 10, height = 4.5)
par(mfrow = c(1,2), mar = c(4.2,4.2,2,1))
plot(tPrec, DMediaY_o2, type = "l", lwd = 2, col = "grey40",
     xlab = "Tiempo (pentadas)", ylab = "mm/pentada",
     main = "Derivada de la media (estaciones)", ylim = range(c(DMediaY_o2, DMediaY_o4, 0)))
lines(tPrec, DMediaY_o4, lwd = 3, col = "slategray")
abline(h = 0, lty = 3)
legend("topleft", c("B-splines orden 2", "B-splines orden 4"), col = c("grey40", "slategray"),
       lwd = c(2,3), bty = "n")
plot(tPrec, DMediaY, type = "l", lwd = 3, col = "slategray", ylim = range(c(DMediaY, DMediaX, 0)),
     xlab = "Tiempo (pentadas)", ylab = "mm/pentada", main = "Intensidad media de la lluvia")
lines(tPrec, DMediaX, lwd = 3, lty = 2, col = "darkturquoise")
abline(h = 0, lty = 3)
legend("topleft", c("Estaciones", "Satélite"), col = c("slategray", "darkturquoise"),
       lty = c(1,2), lwd = 3, bty = "n")
dev.off()

#Pentadas de máxima intensidad media
tPrec[which.max(DMediaY)]
tPrec[which.max(DMediaX)]

#Proporción de la malla con derivada negativa (incoherente con una curva acumulada)
mean(eval.fd(tPrec, DY_fd) < 0)
mean(eval.fd(tPrec, DX_fd) < 0)






#Diferencia satélite - estaciones: d_i(t) = x_i(t) - y_i(t)
Dif_fd       = X_sat_fd - Y_est_fd
Dif_fvalores = eval.fd(tPrec, Dif_fd)
MediaDif     = rowMeans(Dif_fvalores)
EEDif        = apply(Dif_fvalores, 1, sd)/sqrt(ncol(Dif_fvalores))   #Error estándar de la media

windows()
matplot(tPrec, Dif_fvalores, type = "l", lty = 1, col = "grey60",
        xlab = "Tiempo (pentadas)", ylab = "Satélite - estación (mm)")
lines(tPrec, MediaDif, lwd = 3, col = "darkred")
lines(tPrec, MediaDif + 2*EEDif, lty = 2, col = "darkred")
lines(tPrec, MediaDif - 2*EEDif, lty = 2, col = "darkred")
abline(h = 0, lty = 3)

pdf("Figuras_Sesion3/Fig_Precip_Diferencia.pdf", width = 10, height = 4.5)
par(mfrow = c(1,2), mar = c(4.2,4.2,2,1))
matplot(tPrec, Dif_fvalores, type = "l", lty = 1, col = "grey60",
        xlab = "Tiempo (pentadas)", ylab = "Satélite - estación (mm)", main = "Diferencias por sitio")
lines(tPrec, MediaDif, lwd = 3, col = "darkred")
abline(h = 0, lty = 3)
plot(tPrec, MediaDif, type = "l", lwd = 3, col = "darkred",
     ylim = range(c(MediaDif - 2*EEDif, MediaDif + 2*EEDif)),
     xlab = "Tiempo (pentadas)", ylab = "mm", main = "Diferencia media \u00B1 2 EE")
lines(tPrec, MediaDif + 2*EEDif, lty = 2, col = "darkred")
lines(tPrec, MediaDif - 2*EEDif, lty = 2, col = "darkred")
abline(h = 0, lty = 3)
dev.off()

MediaDif[length(tPrec)]
EEDif[length(tPrec)]

#Pentadas donde la banda no contiene el cero
range(tPrec[MediaDif - 2*EEDif > 0])






#Profundidad funcional usando las funciones de la libreria "fda.usc"
SO2_fdata = fdata(t(SO2_fvalores), tSO2, rangval_SO2)
Y_fdata   = fdata(t(Y_fvalores), tPrec, rangval_Prec)

#Profundidad de Fraiman-Muniz y profundidad modal
ProfFM_SO2   = depth.FM(SO2_fdata, trim = 0.25)
ProfMode_SO2 = depth.mode(SO2_fdata, trim = 0.25)

names(ProfFM_SO2)
round(ProfFM_SO2$dep, 3)
round(ProfMode_SO2$dep, 3)

#Orden de las curvas de la más profunda a la menos profunda
order(ProfFM_SO2$dep, decreasing = TRUE)

#Mediana funcional (curva más profunda) y media recortada al 25 %
MedianaSO2  = ProfFM_SO2$median$data
MediaRecSO2 = ProfFM_SO2$mtrim$data

#Día de la mediana y días incluidos en la media recortada
colnames(SO2_fvalores)[ProfFM_SO2$lmed]
colnames(SO2_fvalores)[ProfFM_SO2$ltrim]

windows()
matplot(tSO2, SO2_fvalores, type = "l", lty = 1, col = "grey70", xlab = "Tiempo (horas)",
        ylab = expression(SO[2]))
lines(tSO2, MediaSO2, lwd = 3, col = "blue")
lines(tSO2, MedianaSO2, lwd = 3, col = "red")
lines(tSO2, MediaRecSO2, lwd = 3, lty = 2, col = "darkgreen")

pdf("Figuras_Sesion3/Fig_SO2_Profundidad.pdf", width = 10, height = 4.5)
par(mfrow = c(1,2), mar = c(4.2,4.2,2,1))
grises = gray(seq(0.15, 0.85, length.out = N_SO2))
orden  = order(ProfFM_SO2$dep, decreasing = TRUE)
matplot(tSO2, SO2_fvalores[, orden], type = "l", lty = 1, col = grises, lwd = 2,
        xlab = "Tiempo (horas)", ylab = expression(SO[2]), main = "Curvas ordenadas por profundidad FM")
legend("topright", c("más profunda", "menos profunda"), col = c(grises[1], grises[N_SO2]),
       lwd = 2, bty = "n")
matplot(tSO2, SO2_fvalores, type = "l", lty = 1, col = "grey80",
        xlab = "Tiempo (horas)", ylab = expression(SO[2]), main = "Medidas de centralidad")
lines(tSO2, MediaSO2, lwd = 3, col = "blue")
lines(tSO2, MedianaSO2, lwd = 3, col = "red")
lines(tSO2, MediaRecSO2, lwd = 3, lty = 2, col = "darkgreen")
legend("topright", c("Media", "Mediana (FM)", "Media recortada 25 %"),
       col = c("blue", "red", "darkgreen"), lty = c(1,1,2), lwd = 3, bty = "n")
dev.off()


#Boxplot funcional con la profundidad de bandas modificada (MBD)
windows()
par(mfrow = c(1,2))
BoxSO2 = fbplot(SO2_fvalores, tSO2, method = "MBD", xlab = "Tiempo (horas)", ylab = "SO2",
                xlim = rangval_SO2, ylim = c(0, 20))
BoxY   = fbplot(Y_fvalores, tPrec, method = "MBD", xlab = "Tiempo (pentadas)",
                ylab = "Precipitación (mm)", xlim = rangval_Prec, ylim = c(0, 2500))

pdf("Figuras_Sesion3/Fig_Boxplot.pdf", width = 10, height = 4.5)
par(mfrow = c(1,2), mar = c(4.2,4.2,2,1))
BoxSO2 = fbplot(SO2_fvalores, tSO2, method = "MBD", xlab = "Tiempo (horas)", ylab = "SO2",
                xlim = rangval_SO2, ylim = c(0, 20))
title("SO2 trihorario")
BoxY   = fbplot(Y_fvalores, tPrec, method = "MBD", xlab = "Tiempo (pentadas)",
                ylab = "Precipitación (mm)", xlim = rangval_Prec, ylim = c(0, 2500))
title("Precipitación de estaciones")
dev.off()

#Profundidad MBD y curvas atípicas
round(BoxSO2$depth, 3)
BoxSO2$outpoint
round(BoxY$depth, 3)
BoxY$outpoint

#Profundidad de las curvas de satélite
round(fbplot(X_fvalores, tPrec, method = "MBD", plot = FALSE)$depth, 3)






#Tabla de valores de media, DE, CV e intensidad máxima de la precipitación
cbind.data.frame("Fuente"         = c("Satélite", "Estaciones"),
                 "Media_final"    = round(c(MediaX[length(tPrec)], MediaY[length(tPrec)]), 1),
                 "DE_final"       = round(c(DEX[length(tPrec)],    DEY[length(tPrec)]), 1),
                 "CV_final"       = round(c(CVX[length(tPrec)],    CVY[length(tPrec)]), 3),
                 "Max_intensidad" = round(c(max(DMediaX), max(DMediaY)), 1))

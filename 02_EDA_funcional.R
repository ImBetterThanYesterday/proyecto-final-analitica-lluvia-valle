# =============================================================================
# 02_EDA_FUNCIONAL.R
# =============================================================================
#
# Proyecto:
# Lluvia Valle del Cauca 1981–2003
#
# Objetivo:
# Exploración funcional de las curvas construidas mediante
# B-splines cúbicas.
#
# Entradas:
#
# salidas_R_datos/construccion_funcional.rds
# salidas_R_datos/lluvia_valle_limpia.rds
#
# Salidas:
#
# medias funcionales
# covarianzas
# correlaciones
# profundidades
# boxplot funcional
# ENSO
# grupos altitudinales
# vertientes
#
# =============================================================================
CARPETA_PROYECTO <-
  "/Users/mariavalencia/Downloads/proyecto-final-analitica-lluvia-valle-gus-estaciones-valle"

setwd(CARPETA_PROYECTO)

SALIDA <- file.path(
  CARPETA_PROYECTO,
  "salidas_R_datos"
)

dir.create(
  file.path(
    SALIDA,
    "graficos_EDA"
  ),
  recursive = TRUE,
  showWarnings = FALSE
)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
obj <- readRDS(
  file.path(
    SALIDA,
    "construccion_funcional.rds"
  )
)

lluvia <- readRDS(
  file.path(
    SALIDA,
    "lluvia_valle_limpia.rds"
  )
)
Phi <- obj$Phi

C <- obj$coeficientes

Y <- obj$Y_hat

tiempo <- obj$tiempo

meta <- lluvia$estaciones

id <- obj$id

K <- obj$K

cat(
  "\nCurvas:",
  nrow(Y),
  "\nDias:",
  ncol(Y),
  "\nK:",
  K,
  "\n"
)
mu <- colMeans(
  Y
)

varianza <- apply(
  Y,
  2,
  var
)

de <- sqrt(
  varianza
)

cv <- de/mu
png(
  file.path(
    SALIDA,
    "graficos_EDA",
    "media_funcional.png"
  ),
  width = 1200,
  height = 800
)

plot(
  tiempo,
  mu,
  type = "l",
  lwd = 3,
  col = "blue",
  xlab = "Dia",
  ylab = "mm/dia",
  main = "Media funcional"
)

dev.off()
png(
  file.path(
    SALIDA,
    "graficos_EDA",
    "varianza_funcional.png"
  ),
  width = 1200,
  height = 800
)

plot(
  tiempo,
  varianza,
  type = "l",
  lwd = 3,
  xlab = "Dia",
  ylab = "Varianza",
  main = "Varianza funcional"
)

dev.off()
Cc <- sweep(
  C,
  2,
  colMeans(C)
)

covar <- Phi %*%
  (crossprod(Cc)/(nrow(C)-1)) %*%
  t(Phi)

write.csv(
  covar,
  file.path(
    SALIDA,
    "covarianza_funcional.csv"
  )
)
correlacion <- cov2cor(
  covar
)

write.csv(
  correlacion,
  file.path(
    SALIDA,
    "correlacion_funcional.csv"
  )
)
png(
  file.path(
    SALIDA,
    "graficos_EDA",
    "covarianza_funcional.png"
  ),
  width = 1000,
  height = 900
)

image(
  tiempo,
  tiempo,
  covar,
  main = "Covarianza funcional"
)

dev.off()


# =============================================================================
# 01_DATOS_Y_CONSTRUCCION_FUNCIONAL.R
# =============================================================================
#
# Proyecto:
# Lluvia Valle del Cauca 1981–2003
#
# Objetivo:
# Construcción funcional de las curvas anuales de precipitación
# mediante B-splines cúbicas implementadas con la recursión de
# Cox-de Boor y estimación por mínimos cuadrados.
#
# Entradas:
#   salidas_R_datos/lluvia_valle_limpia.rds
#
# Salidas:
#   construccion_funcional.rds
#   Phi_estaciones.csv
#   criterios_K_estaciones.csv
#   controles_S.csv
#   cov_residual_estaciones.csv
#   banda_puntual_estaciones.csv
#   graficos_funcionales/
#
# Metodología:
#   - Selección de K mediante ECM, VC y VCG.
#   - K seleccionado exclusivamente con entrenamiento (1981–1995).
#   - Valor final: mínimo de VCG.
#   - Construcción sin librerías FDA.
#   - Verificación numérica de la matriz S.
#
# =============================================================================

CARPETA_PROYECTO <-
  "/Users/mariavalencia/Downloads/proyecto-final-analitica-lluvia-valle-gus-estaciones-valle"

setwd(CARPETA_PROYECTO)

lluvia <- readRDS(
  file.path(
    CARPETA_PROYECTO,
    "salidas_R_datos",
    "lluvia_valle_limpia.rds"
  )
)

names(lluvia)

# =============================================================================
# CONFIGURACION GENERAL
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

theme_set(theme_minimal())

options(stringsAsFactors = FALSE)

SALIDA <- file.path(
  CARPETA_PROYECTO,
  "salidas_R_datos"
)

dir.create(
  file.path(
    SALIDA,
    "graficos_funcionales"
  ),
  recursive = TRUE,
  showWarnings = FALSE
)

cat("\n----------------------------------\n")
cat("CONSTRUCCION FUNCIONAL\n")
cat("----------------------------------\n")

# =============================================================================
# DATOS FUNCIONALES
# =============================================================================

d <- as.data.table(
  lluvia$diario
)

meta <- as.data.table(
  lluvia$estaciones
)

cat(
  "\nEstaciones:",
  uniqueN(d$codigo)
)

cat(
  "\nCurvas estacion-anio:",
  uniqueN(
    paste(
      d$codigo,
      d$anio
    )
  )
)

cat(
  "\nObservaciones:",
  nrow(d),
  "\n"
)

# =============================================================================
# MATRIZ FUNCIONAL
# =============================================================================

d[
  ,
  curva_id :=
    paste(
      codigo,
      anio,
      sep = "_"
    )
]

matriz_dt <- dcast(
  d,
  curva_id +
    codigo +
    anio +
    particion ~ dia,
  value.var = "P"
)

id_curvas <- matriz_dt[
  ,
  .(
    curva_id,
    codigo,
    anio,
    particion
  )
]

Y <- as.matrix(
  matriz_dt[
    ,
    !c(
      "curva_id",
      "codigo",
      "anio",
      "particion"
    ),
    with = FALSE
  ]
)

storage.mode(Y) <- "double"

tiempo <- 1:365

cat(
  "\nDimension de matriz:",
  dim(Y)[1],
  "x",
  dim(Y)[2]
)

cat(
  "\nNA totales:",
  sum(is.na(Y)),
  "\n"
)

# =============================================================================
# PARTICIONES
# =============================================================================

idx_ent <-
  id_curvas$particion ==
  "entrenamiento"

idx_val <-
  id_curvas$particion ==
  "validación"

idx_pru <-
  id_curvas$particion ==
  "prueba"

Y_ent <- Y[idx_ent,]

Y_val <- Y[idx_val,]

Y_pru <- Y[idx_pru,]

cat(
  "\nEntrenamiento:",
  nrow(Y_ent)
)

cat(
  "\nValidacion:",
  nrow(Y_val)
)

cat(
  "\nPrueba:",
  nrow(Y_pru),
  "\n"
)

# =============================================================================
# B-SPLINES PROPIAS
# =============================================================================

RANGO <- c(1,365)

ORDEN <- 4L

nodos <- function(
    K,
    orden = ORDEN,
    rango = RANGO
){
  
  interiores <- seq(
    rango[1],
    rango[2],
    length.out =
      K - orden + 2L
  )
  
  c(
    rep(
      rango[1],
      orden
    ),
    interiores[
      -c(
        1,
        length(interiores)
      )
    ],
    rep(
      rango[2],
      orden
    )
  )
  
}

base_en_nodos <- function(
    tiempo,
    tau,
    orden
){
  
  B <- outer(
    tiempo,
    seq_len(
      length(tau)-1L
    ),
    function(x,j)
      (x >= tau[j]) &
      (x < tau[j+1L])
  ) * 1
  
  B[
    tiempo == max(tau),
    max(which(diff(tau)>0))
  ] <- 1
  
  if(orden > 1){
    
    for(m in 2:orden){
      
      siguiente <- matrix(
        0,
        length(tiempo),
        length(tau)-m
      )
      
      for(j in seq_len(
        ncol(
          siguiente
        )
      )){
        
        d1 <-
          tau[j+m-1L] -
          tau[j]
        
        d2 <-
          tau[j+m] -
          tau[j+1L]
        
        if(d1 > 0){
          
          siguiente[,j] <-
            siguiente[,j] +
            (tiempo-tau[j])/
            d1*
            B[,j]
          
        }
        
        if(d2 > 0){
          
          siguiente[,j] <-
            siguiente[,j] +
            (tau[j+m]-tiempo)/
            d2*
            B[,j+1L]
          
        }
        
      }
      
      B <- siguiente
      
    }
    
  }
  
  B
  
}

base_bspline <- function(
    tiempo,
    K
){
  
  base_en_nodos(
    tiempo,
    nodos(K),
    ORDEN
  )
  
}

# =============================================================================
# MINIMOS CUADRADOS
# =============================================================================

ajuste_mco <- function(
    y,
    Phi
){
  
  observado <- !is.na(y)
  
  P <- Phi[
    observado,
    ,
    drop = FALSE
  ]
  
  if(
    nrow(P) <= ncol(P)
  ){
    return(NULL)
  }
  
  if(
    qr(P)$rank !=
    ncol(P)
  ){
    return(NULL)
  }
  
  mapa <- solve(
    crossprod(P),
    t(P)
  )
  
  coeficientes <-
    as.vector(
      mapa %*%
        y[observado]
    )
  
  list(
    coef =
      coeficientes,
    mapa =
      mapa
  )
  
}

# =============================================================================
# CRITERIOS DE SELECCION DE K
# =============================================================================
#
# K se selecciona automaticamente.
#
# Se comparan:
#
# ECM
# VC
# VCG
#
# utilizando EXCLUSIVAMENTE
# entrenamiento 1981-1995.
#
# El K final se elige
# con el minimo de VCG.
#
# =============================================================================

criterios <- function(
    Y,
    tiempo,
    Ks
){
  
  completas <- rowSums(is.na(Y)) == 0
  
  Yc <- Y[
    completas,
    ,
    drop = FALSE
  ]
  
  G <- crossprod(Yc)
  
  n <- ncol(Yc)
  
  N <- nrow(Yc)
  
  salida <- list()
  
  for(K in Ks){
    
    Phi <- base_bspline(
      tiempo,
      K
    )
    
    if(qr(Phi)$rank != K){
      next
    }
    
    mapa <- solve(
      crossprod(Phi),
      t(Phi)
    )
    
    h <- rowSums(
      Phi * t(mapa)
    )
    
    T0 <- G %*% t(mapa)
    
    sse <- pmax(
      0,
      diag(G) -
        2 * rowSums(Phi * T0) +
        rowSums(
          (
            Phi %*%
              (mapa %*% T0)
          ) * Phi
        )
    )
    
    salida[[length(salida)+1]] <-
      data.frame(
        K = K,
        ECM = sum(sse)/(N*n),
        VC = sum(sse/(1-h)^2)/(N*n),
        VCG = n*sum(sse)/(n-K)^2
      )
    
  }
  
  rbindlist(salida)
  
}

Ks <- 4:40

seleccion_K <- criterios(
  Y_ent,
  tiempo,
  Ks
)

# =============================================================================
# GRAFICAS DE SELECCION DE K
# =============================================================================

K_opt <- seleccion_K$K[
  which.min(
    seleccion_K$VCG
  )
]

cat(
  "\nK seleccionado por VCG:",
  K_opt,
  "\n"
)

for(crit in c(
  "ECM",
  "VC",
  "VCG"
)){
  
  png(
    file.path(
      SALIDA,
      "graficos_funcionales",
      paste0(
        "criterio_",
        crit,
        "_estaciones.png"
      )
    ),
    width = 1200,
    height = 800
  )
  
  plot(
    seleccion_K$K,
    seleccion_K[[crit]],
    type = "b",
    pch = 16,
    col = "steelblue",
    xlab = "K",
    ylab = crit,
    main = paste(
      crit,
      "- Seleccion de K"
    )
  )
  
  abline(
    v = K_opt,
    col = "red",
    lty = 2,
    lwd = 2
  )
  
  dev.off()
  
}

# =============================================================================
# MATRIZ BASE PHI
# =============================================================================

Phi <- base_bspline(
  tiempo,
  K_opt
)

write.csv(
  data.frame(
    dia = tiempo,
    Phi
  ),
  file.path(
    SALIDA,
    "Phi_estaciones.csv"
  ),
  row.names = FALSE
)

cat(
  "\nDimension Phi:",
  dim(Phi)[1],
  "x",
  dim(Phi)[2],
  "\n"
)

# =============================================================================
# MATRIZ DE SUAVIZADO S
# =============================================================================

mapa <- solve(
  crossprod(Phi),
  t(Phi)
)

S <- Phi %*% mapa

control_S <- data.frame(
  
  simetria =
    max(
      abs(
        S - t(S)
      )
    ),
  
  idempotencia =
    max(
      abs(
        S %*% S - S
      )
    ),
  
  traza =
    sum(
      diag(S)
    ),
  
  K =
    K_opt
  
)

print(control_S)

write.csv(
  control_S,
  file.path(
    SALIDA,
    "controles_S.csv"
  ),
  row.names = FALSE
)

# =============================================================================
# COEFICIENTES FUNCIONALES
# =============================================================================

coeficientes <- matrix(
  NA_real_,
  nrow(Y),
  K_opt
)

for(i in seq_len(
  nrow(Y)
)){
  
  ajuste <- ajuste_mco(
    Y[i,],
    Phi
  )
  
  if(
    is.null(ajuste)
  ){
    next
  }
  
  coeficientes[i,] <-
    ajuste$coef
  
}

colnames(
  coeficientes
) <- paste0(
  "c",
  seq_len(K_opt)
)

cat(
  "\nCoeficientes estimados:",
  nrow(coeficientes),
  "curvas\n"
)

# =============================================================================
# RECONSTRUCCION DE CURVAS
# =============================================================================

Y_hat <- coeficientes %*% t(Phi)

cat(
  "\nDimension Y_hat:",
  dim(Y_hat)[1],
  "x",
  dim(Y_hat)[2],
  "\n"
)

# =============================================================================
# RESIDUOS
# =============================================================================

residuos <- Y

for(i in seq_len(nrow(Y))){
  
  obs <- !is.na(Y[i,])
  
  residuos[
    i,
    obs
  ] <-
    Y[i,obs] -
    Y_hat[i,obs]
  
}
# =============================================================================
# COVARIANZA RESIDUAL
# =============================================================================

completas <- rowSums(
  is.na(Y)
) == 0

Res_comp <- residuos[
  completas,
  ,
  drop = FALSE
]

Sigma_res <- crossprod(
  Res_comp
) /
  nrow(
    Res_comp
  )

write.csv(
  Sigma_res,
  file.path(
    SALIDA,
    "cov_residual_estaciones.csv"
  )
)

cat(
  "\nDimension covarianza residual:",
  dim(Sigma_res)[1],
  "x",
  dim(Sigma_res)[2],
  "\n"
)
# =============================================================================
# BANDAS PUNTO A PUNTO
# =============================================================================

ejemplo <- which(
  rowSums(is.na(Y)) == 0
)[1]

ajuste_ejemplo <-
  ajuste_mco(
    Y[ejemplo,],
    Phi
  )

sigma2 <- mean(
  residuos[
    ejemplo,
  ]^2,
  na.rm = TRUE
)

ee <- sqrt(
  
  sigma2 *
    
    rowSums(
      
      (
        Phi %*%
          ajuste_ejemplo$mapa
      )^2
      
    )
  
)

banda <- data.frame(
  
  dia = tiempo,
  
  curva = Y_hat[
    ejemplo,
  ],
  
  ee = ee,
  
  inferior =
    Y_hat[
      ejemplo,
    ] -
    2*ee,
  
  superior =
    Y_hat[
      ejemplo,
    ] +
    2*ee
  
)

write.csv(
  banda,
  file.path(
    SALIDA,
    "banda_puntual_estaciones.csv"
  ),
  row.names = FALSE
)
png(
  file.path(
    SALIDA,
    "graficos_funcionales",
    "banda_puntual_estaciones.png"
  ),
  width = 1200,
  height = 800
)

plot(
  tiempo,
  Y[ejemplo,],
  pch = 16,
  cex = .4,
  main = "Banda puntual",
  xlab = "Dia",
  ylab = "mm/dia"
)

lines(
  tiempo,
  banda$curva,
  col = "blue",
  lwd = 2
)

lines(
  tiempo,
  banda$inferior,
  lty = 2
)

lines(
  tiempo,
  banda$superior,
  lty = 2
)

dev.off()

# =============================================================================
# CONSTRUCCION FINAL
# =============================================================================

obj_construccion <- list(
  
  id = id_curvas,
  
  meta = meta,
  
  K = K_opt,
  
  tiempo = tiempo,
  
  Phi = Phi,
  
  S = S,
  
  control_S = control_S,
  
  coeficientes = coeficientes,
  
  Y_hat = Y_hat,
  
  seleccion_K = seleccion_K,
  
  Sigma_res = Sigma_res
  
)

saveRDS(
  
  obj_construccion,
  
  file.path(
    SALIDA,
    "construccion_funcional.rds"
  )
  
)

cat(
  "\n======================================\n"
)

cat(
  "CONSTRUCCION FUNCIONAL TERMINADA\n"
)

cat(
  "======================================\n"
)

cat(
  "\nK final:",
  K_opt
)

cat(
  "\nCurvas:",
  nrow(Y)
)

cat(
  "\nArchivo:",
  file.path(
    SALIDA,
    "construccion_funcional.rds"
  ),
  "\n"
)

#verificacioneas
names(
  readRDS(
    file.path(
      SALIDA,
      "construccion_funcional.rds"
    )
  )
)
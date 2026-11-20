# =============================================================================
# LLUVIA DEL VALLE DEL CAUCA 1981–2003 · ENTENDIMIENTO, DEPURACIÓN Y DATASETS FINALES
# =============================================================================
# Un solo archivo, de arriba abajo, desde los datos crudos hasta el dataset limpio.
#
#   0. Configuración
#   1. Población: las 232 estaciones del Valle (D-34/D-35) y prueba geográfica
#   2. Estructura de los archivos crudos
#   3. Disponibilidad: rejilla estación × día, cobertura y huecos
#   4. Valores: distribución, extremos y contraste con vecinas
#   5. Patrones sospechosos: totales mensuales, rachas de ceros, copias
#   6. Climatología y coherencia espacial
#   7. Experimento de huecos: ¿cuánta falta de datos aguanta una curva?
#   8. Depuración (reglas R0 a R7, D-36 y D-61: R7 rellena por kriging) y dataset final
#   9. Satélite: lectura, estructura y fuente (D-43)
#  10. Satélite: recorte al Valle, calendario de péntadas y mm/día
#  11. Satélite: EDA
#  12. Satélite frente a estaciones: emparejamiento, sesgo y años sospechosos
#  13. Datos preparados para la construcción funcional (matrices y pares)
#  14. Pruebas de los datos preparados
#  15. Resumen
#
# Ejecutar desde la raíz del repositorio:  Rscript lluvia_valle_datos.R
# Conjunto de sensibilidad de R6 (regla de la OMM, ≥19 de 23 años, D-63):  CONJUNTO=r6_omm Rscript lluvia_valle_datos.R
# (o abrirlo en RStudio con la raíz del repositorio como directorio de trabajo
#  y pulsar Source). Tarda unos minutos, sobre todo la sección 7.
#
# Entradas
#   datos_curso/P_1980_1990.csv, P_1991_2003.csv  lluvia diaria (STATION_ID, YEAR, MONTH, DAY, VALUE)
#   datos_web/estaciones_valle_finales.csv        las 232 estaciones del Valle (D-34/D-35)
#   datos_web/estaciones_valle_del_cauca.xlsx     inventario de los 344 códigos (D-31/D-32)
#   datos_web/limites/deptos_valle_y_vecinos_DANE_2020.gpkg   límites DANE (MGN 2020)
#   datos_web/oni_anual_fase.csv                  ONI anual (contexto ENSO)
#   datos_curso/ImagenesSatelitePentadal/ppt*.tif satélite por péntadas (2 328 GeoTIFF)
#   datos_web/chirp_v2_valle_muestra_7anios.bin   CHIRP v2.0 del CHC, para verificar la fuente (D-43)
#   datos_web/chirps_v2_valle_1981_2013.bin       CHIRPS v2.0 del CHC, para verificar la fuente
# Salidas (salidas_R_datos/)
#   lluvia_valle_limpia.rds   dataset final de estaciones (ver sección 8)
#   satelite_valle.rds        satélite recortado, en mm/día, y pares estación–píxel (sección 12)
#   datos_preparados.rds      matrices listas para la construcción funcional (sección 13)
#   tablas/*.csv              tablas de cada sección
#   graficos/*.png            gráficos de cada sección
#
# Librerías: data.table, ggplot2, readxl, sf y terra (EDA). La base B-spline del
# experimento de huecos se programa sin librerías (recursión de Cox–de Boor),
# como pide el curso para la construcción funcional.
# Decisiones y su justificación: DECISIONES.md (D-01 a D-61).
# =============================================================================


# ---- 0. CONFIGURACIÓN -------------------------------------------------------
# Locale: sin UTF-8 (LANG vacío en la Mac del grupo) las tildes del código no coinciden
# con las de los datos. Se fija aquí, antes de cualquier cadena con tilde.
if (!isTRUE(l10n_info()$`UTF-8`)) {
  for (loc in c("en_US.UTF-8", "es_CO.UTF-8", "es_ES.UTF-8", "C.UTF-8"))
    if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break
  if (!isTRUE(l10n_info()$`UTF-8`)) stop("No hay un locale UTF-8 disponible; corra con LC_ALL=en_US.UTF-8 Rscript ...")
}
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(readxl)
  library(sf)
})
sf_use_s2(FALSE)
theme_set(theme_minimal(base_size = 11))
set.seed(20261112)

if (!file.exists("datos_curso/P_1991_2003.csv")) stop("Ejecute desde la raíz del repositorio y copie los datos del curso en datos_curso/")

ANIOS       <- 1981:2003                    # fe de erratas
PARTICION   <- function(anio) fcase(anio <= 1995, "entrenamiento", anio <= 2000, "validación", default = "prueba")
RADIO_VEC   <- 25                           # km: vecinas para cualquier contraste (H9)
# R5 y R6 siguen la guía de la OMM para datos faltantes (WMO-No. 1203, 2017, sección 4.4;
# fundamentos en docs/FUNDAMENTOS_REGLAS.md). Los días faltantes de un año válido se
# ESTIMAN con las vecinas (R7, kriging); la OMM (4.4.1 c y 4.6) acepta estimar días solo
# si en el mes no son 11 o más, ni 5 o más seguidos.
MES_FALT_MAX   <- 10                        # R5: faltantes máximos en un mes (OMM: 11 o más lo anulan)
MES_SEGUID_MAX <- 4                         # R5: faltantes seguidos máximos en un mes (OMM: 5 o más lo anulan)
# R6 tiene dos conjuntos (D-61, D-63). Se elige con la variable de entorno CONJUNTO:
#   principal (por defecto): ≥10 años válidos en 1981–2003. Es el mínimo de la OMM para
#       cualquier promedio climatológico ("period averages ... for any period of at least ten
#       years", Reglamento Técnico, citado en WMO-No. 1203, 2017, sección 3) y el de las normales
#       "representativas" de NOAA 1991–2020 (≥10 años con los huecos rellenados, como hace R7).
#       Nuestra unidad es la curva estación-año, no la normal de 30 años de una estación.
#   r6_omm (sensibilidad):   años válidos en ≥80 % de 1981–2003 = 19 de 23 (OMM 2017, 4.4.2,
#       la regla para normales estándar).
#   CONJUNTO=r6_omm Rscript lluvia_valle_datos.R   → escribe en salidas_R_datos_r6_omm/
CONJUNTO <- Sys.getenv("CONJUNTO", "principal")
if (!CONJUNTO %in% c("principal", "r6_omm")) stop("CONJUNTO debe ser 'principal' o 'r6_omm'")
SUFIJO   <- if (CONJUNTO == "principal") "" else paste0("_", CONJUNTO)
PCT_ANIOS_MIN  <- 0.80                      # R6 de sensibilidad: fracción mínima de años válidos (OMM 4.4.2)
ANIOS_MIN      <- if (CONJUNTO == "principal") 10L else ceiling(PCT_ANIOS_MIN * length(ANIOS))   # 10 o 19 de 23
R6_TEXTO       <- if (CONJUNTO == "principal") "≥10 años válidos (mínimo OMM para promedios; NOAA)" else "≥80 % de 1981–2003 (19 de 23, OMM, sensibilidad)"
# R0: máximo físico = récord mundial de lluvia en 24 h reconocido por la OMM: 1 825 mm,
# Foc-Foc (La Reunión), 7-8 ene 1966, ciclón Denise (WMO Weather and Climate Extremes
# Archive, https://wmo.asu.edu). Un dato diario por encima no puede ser lluvia real.
P_MAX_FISICO <- 1825

SALIDA <- paste0("salidas_R_datos", SUFIJO)
cat("Conjunto:", CONJUNTO, "| R6:", R6_TEXTO, "| salidas en", SALIDA, "\n")
dir.create(file.path(SALIDA, "tablas"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(SALIDA, "graficos"), showWarnings = FALSE, recursive = TRUE)
tabla   <- function(x, nombre) { fwrite(x, file.path(SALIDA, "tablas", paste0(nombre, ".csv"))); invisible(x) }
grafico <- function(p, nombre, ancho = 9, alto = 5, dpi = 110) {
  ggsave(file.path(SALIDA, "graficos", paste0(nombre, ".png")), p, width = ancho, height = alto, dpi = dpi, bg = "white")
  if (interactive()) print(p)
  invisible(p)
}
titulo <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
# Gráficas con R base en blanco y negro, con un panel por año (Ramsay y Silverman 2005, cap. 1)
grafico_r <- function(nombre, dibujar, ancho = 11, alto = 15) {
  png(file.path(SALIDA, "graficos", paste0(nombre, ".png")), width = ancho, height = alto, units = "in", res = 110)
  on.exit(dev.off()); dibujar(); invisible(NULL)
}

dist_km <- function(lat1, lon1, lat2, lon2) {   # distancia de gran círculo
  r <- pi / 180
  a <- sin((lat2 - lat1) * r / 2)^2 + cos(lat1 * r) * cos(lat2 * r) * sin((lon2 - lon1) * r / 2)^2
  12742 * asin(sqrt(a))
}
hueco_mas_largo <- function(obs) { r <- rle(!obs); if (any(r$values)) max(r$lengths[r$values]) else 0L }
MES_DIA <- month(as.IDate("2001-01-01") + 0:364)   # mes de cada día del dominio 1–365
# R5 (OMM 2017, 4.4.1): el año vale si en TODOS sus meses faltan ≤10 días y no hay 5 o más seguidos
mes_valido_omm <- function(o, falt = MES_FALT_MAX, seguid = MES_SEGUID_MAX) sum(!o) <= falt && hueco_mas_largo(o) <= seguid
anio_valido_omm <- function(obs, ...) all(vapply(split(obs, MES_DIA), mes_valido_omm, logical(1), ...))


# ---- 1. POBLACIÓN: LAS 232 ESTACIONES DEL VALLE -----------------------------
titulo("1. Población")
# Los 344 códigos de los tres archivos de lluvia se ubicaron con fuentes
# oficiales (D-31/D-32). Quedaron 232 del Valle con coordenada: 226 con ubicación
# oficial y 6 identificadas por cruce con el CNE del IDEAM (D-34). La altitud
# solo sale de fuentes oficiales (D-35).

crudo <- rbind(fread("datos_curso/P_1980_1990.csv"), fread("datos_curso/P_1991_2003.csv"))
crudo[, codigo := as.character(STATION_ID)][, STATION_ID := NULL]
crudo <- crudo[YEAR %in% ANIOS]
cobertura_cruda <- crudo[, .(dias_con_dato = .N, anio_ini = min(YEAR), anio_fin = max(YEAR),
                             anios_con_dato = uniqueN(YEAR)), by = codigo]

inventario <- as.data.table(read_excel("datos_web/estaciones_valle_del_cauca.xlsx", sheet = 2))
inventario <- inventario[, .(codigo = format(inventario[[1]], scientific = FALSE, trim = TRUE),
                             depto_inventario = Departamento, nivel_inventario = `Nivel de evidencia`)]

finales <- fread("datos_web/estaciones_valle_finales.csv", encoding = "UTF-8", colClasses = list(character = "codigo"))
finales <- finales[, .(codigo, nombre, municipio, lat, lon, altitud = alt, nivel_c, fuente, fuente_alt, srtm, revisar)]
# Univalle (2622320201) no tiene altitud en el Excel del curso ni en el CNE. Se
# toma la de la estación del IDEAM en el mismo campus, "UNIVERSIDAD DEL VALLE - AUT"
# (CNE 0026055120, 996 m, a 1,6 km). El SRTM en su coordenada da 977 m: 19 m de
# diferencia, dentro del error típico del SRTM en Cali (mediana 15 m). D-46.
finales[codigo == "2622320201", `:=`(altitud = 996, fuente_alt = "CNE IDEAM: estación 0026055120 del mismo campus (D-46)")]
# Altitud imposible → terreno (D-59). En la cuenca del Cauca (código 26) ningún terreno de
# estación baja de PISO_CAUCA, el mínimo del SRTM en las estaciones de esa cuenca (904 m). Con la
# misma tolerancia de 300 m del control de coordenadas (D-35), una altitud oficial por debajo de
# PISO_CAUCA − 300 no puede ser real: es un error de la fuente, como un valor de lluvia imposible
# (R0). Solo cae El Alcázar (Zarzal): el CNE del IDEAM (0026100400) le da 17 m, el terreno en su
# coordenada está a 977 m y El Medio Hda., la estación oficial más cercana (3 km), a 952 m.
PISO_CAUCA <- finales[substr(codigo, 1, 2) == "26", min(srtm, na.rm = TRUE)]
alt_imposible <- finales[substr(codigo, 1, 2) == "26" & !is.na(altitud) & altitud < PISO_CAUCA - 300,
                         .(codigo, nombre, altitud_oficial = altitud, srtm)]
cat("Piso de la cuenca del Cauca (mínimo SRTM en sus estaciones):", PISO_CAUCA, "m. Altitudes imposibles:\n"); print(alt_imposible)
finales[alt_imposible, on = "codigo", `:=`(altitud = i.srtm,
        fuente_alt = sprintf("SRTM en la coordenada: la oficial (%g m) es imposible en la cuenca del Cauca (D-59)", i.altitud_oficial))]

est <- merge(merge(cobertura_cruda, inventario, by = "codigo", all.x = TRUE), finales, by = "codigo", all = TRUE)
est[is.na(dias_con_dato), `:=`(dias_con_dato = 0L, anios_con_dato = 0L)]
est[, grupo := fifelse(!is.na(nivel_c), "Valle del Cauca",
               fcase(depto_inventario == "Valle", "Valle sin coordenada (D-34)",
                     depto_inventario == "Cauca", "Cauca", default = "Sin ubicar"))]

# Prueba geográfica: ¿cada coordenada cae dentro del polígono DANE del Valle?
deptos <- st_read("datos_web/limites/deptos_valle_y_vecinos_DANE_2020.gpkg", quiet = TRUE)
deptos$shapeName <- enc2utf8(deptos$shapeName)
xy <- est[!is.na(lat)]
pts <- st_as_sf(xy, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
dentro <- suppressMessages(st_intersects(pts, deptos))
xy[, depto_geo := vapply(dentro, function(i) if (length(i)) deptos$shapeName[i[1]] else "Fuera", "")]
# Distancia al departamento distinto del propio más cercano (borde terrestre,
# no la costa), en MAGNA-SIRGAS origen nacional (EPSG:9377)
Dd <- units::drop_units(st_distance(st_transform(pts, 9377), st_transform(deptos, 9377)))
propio <- match(xy$depto_geo, deptos$shapeName)
xy[, km_borde := round(vapply(seq_len(nrow(Dd)), function(i) min(Dd[i, -propio[i]]), 0) / 1000, 2)]
est <- merge(est, xy[, .(codigo, depto_geo, km_borde)], by = "codigo", all.x = TRUE)
tabla(est, "1_estaciones_metadatos")

print(est[, .(codigos = .N, con_datos_81_03 = sum(dias_con_dato > 0), dias = sum(dias_con_dato)), by = grupo])
print(est[grupo == "Valle del Cauca", .N, by = depto_geo])
stopifnot(est[grupo == "Valle del Cauca", all(depto_geo == "Valle del Cauca")])

valle <- est[grupo == "Valle del Cauca" & dias_con_dato > 0]
valle[, vertiente := fifelse(substr(codigo, 1, 2) == "26", "Cauca (interandina)", "Pacífico")]
# Nombre de cada zona en las gráficas. Las dos zonas están DENTRO del Valle del Cauca:
# "Cauca" es el río, no el departamento. El valor interno no cambia porque lo usan las
# partes funcionales (comparan con "Pacífico").
ETQ_VERT <- c("Cauca (interandina)" = "Zona andina (valle del río Cauca)", "Pacífico" = "Zona pacífica (hacia el océano)")
COL_VERT <- setNames(c("#B8860B", "#1F78B4"), ETQ_VERT)
cat("Valle con datos en 1981–2003:", nrow(valle), "de 232 | sin altitud:", valle[is.na(altitud), .N],
    "| coordenada a revisar:", valle[revisar != "", .N], "| a ≤2 km de otro departamento:", valle[km_borde <= 2, .N], "\n")
print(valle[, .N, by = .(nivel_c, vertiente)])

# Mapa de las dos zonas. El fondo pinta cada punto del Valle con la zona de la estación con
# datos más cercana: la parte funcional verifica esa regla (acierta la zona del 97,7 % de las
# estaciones dejando cada una fuera; ningún corte en longitud pasa de 89,2 %). Las ciudades
# son solo referencia, con su coordenada aproximada.
poli_valle <- deptos[deptos$shapeName == "Valle del Cauca", ]
rej <- st_make_grid(poli_valle, cellsize = 0.02, what = "centers")
rej <- st_coordinates(rej[lengths(st_intersects(rej, poli_valle)) > 0])
cerca <- apply(rej, 1, function(p) which.min(((valle$lon - p[1]) * cos(p[2] * pi / 180))^2 + (valle$lat - p[2])^2))
fondo <- data.table(lon = rej[, 1], lat = rej[, 2], zona = ETQ_VERT[valle$vertiente[cerca]])
ciudades <- data.table(ciudad = c("Buenaventura", "Cali", "Palmira", "Buga", "Tuluá", "Cartago"),
                       lon = c(-77.03, -76.53, -76.30, -76.30, -76.20, -75.91),
                       lat = c(3.88, 3.45, 3.54, 3.90, 4.08, 4.75))
print(valle[, .(estaciones_con_datos = .N), by = .(zona = ETQ_VERT[vertiente])])
grafico(ggplot() +
          geom_tile(data = fondo, aes(lon, lat, fill = zona), width = .02, height = .02, alpha = .35) +
          geom_sf(data = poli_valle, inherit.aes = FALSE, fill = NA, colour = "grey30") +
          geom_point(data = est[grupo == "Valle del Cauca"], aes(lon, lat, shape = dias_con_dato > 0), size = 1.5) +
          geom_point(data = ciudades, aes(lon, lat), shape = 23, size = 3, fill = "red", colour = "black") +
          geom_text(data = ciudades, aes(lon, lat, label = ciudad), size = 3.6, fontface = "bold", nudge_y = .08) +
          scale_fill_manual(values = COL_VERT, name = NULL) +
          scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = NULL,
                             labels = c(`TRUE` = "estación con datos en 1981–2003 (220)", `FALSE` = "estación sin datos en el periodo (12)")) +
          guides(fill = guide_legend(override.aes = list(alpha = .6))) +
          labs(title = "Las dos zonas del Valle del Cauca y sus 232 estaciones",
               subtitle = "Ambas zonas están dentro del departamento; las separa la Cordillera Occidental.\nRombos rojos: ciudades de referencia.",
               x = NULL, y = NULL) +
          theme(legend.position = "bottom", legend.box = "vertical"),
        "1_mapa_estaciones", ancho = 7, alto = 8.5)

# Matriz de distancias entre las estaciones del Valle con datos
DK <- outer(seq_len(nrow(valle)), seq_len(nrow(valle)),
            function(i, j) dist_km(valle$lat[i], valle$lon[i], valle$lat[j], valle$lon[j]))
dimnames(DK) <- list(valle$codigo, valle$codigo)
vecinas <- function(cod) setdiff(names(which(DK[cod, ] <= RADIO_VEC)), cod)

# QUÉ DICEN LOS DATOS (H1)
#   Las cifras y su interpretación están en EXPLICACION_DATOS.md, sección 1 (corrida del
#   19-nov-2026, con las 98 estaciones de la regla de la OMM); las tablas de esta
#   sección las imprimen las líneas de arriba.


# ---- 2. ESTRUCTURA DE LOS ARCHIVOS CRUDOS ------------------------------------
titulo("2. Estructura del archivo")
d <- crudo[codigo %in% valle$codigo, .(codigo, anio = YEAR, mes = MONTH, dia_mes = DAY, P = as.numeric(VALUE))]
crudo_r0 <- crudo[, .(codigo, P = as.numeric(VALUE))]   # solo para el conteo de R0
rm(crudo); invisible(gc())
d[, fecha := as.IDate(sprintf("%d-%02d-%02d", anio, mes, dia_mes))]
estructura <- data.table(
  filas = nrow(d), estaciones = uniqueN(d$codigo),
  fechas_invalidas = sum(is.na(d$fecha)),
  duplicados = sum(duplicated(d[, .(codigo, fecha)])),
  valores_NA = sum(is.na(d$P)), negativos = sum(d$P < 0, na.rm = TRUE),
  no_enteros = sum(d$P != round(d$P)), dias_29feb = d[mes == 2 & dia_mes == 29, .N],
  maximo = max(d$P))
print(tabla(estructura, "2_estructura"))

# R0. Valor físicamente imposible → faltante: P < 0 o P > P_MAX_FISICO (récord OMM 24 h).
# La regla se aplica antes de todo lo demás y se cuenta aunque dé 0, para que quede
# verificado en el código y no supuesto. Como referencia se cuenta también en el archivo
# crudo completo (todas las estaciones, 1981-2003), donde sí aparecen valores imposibles.
imposible <- function(P) !is.na(P) & (P < 0 | P > P_MAX_FISICO)
m_r0 <- d[imposible(P), .(codigo, fecha, P, regla = "R0 valor imposible",
                          detalle = fifelse(P < 0, "negativo", sprintf("> %d mm (récord OMM 24 h)", P_MAX_FISICO)))]
r0_resumen <- data.table(
  archivo = c("crudo completo (todas las estaciones)", "estaciones del Valle"),
  negativos = c(sum(crudo_r0$P < 0), m_r0[P < 0, .N]),
  sobre_maximo_fisico = c(sum(crudo_r0$P > P_MAX_FISICO), m_r0[P > P_MAX_FISICO, .N]),
  estaciones_afectadas = c(uniqueN(crudo_r0[imposible(P), codigo]), uniqueN(m_r0$codigo)),
  valor_maximo = c(max(crudo_r0$P), d[, max(P)]))
print(tabla(r0_resumen, "2_R0_valores_imposibles"))
rm(crudo_r0)
d <- d[!imposible(P)]
traza_r0 <- d[, .(estaciones = uniqueN(codigo), pares_estacion_anio = uniqueN(paste(codigo, anio)), dias_con_dato = .N)]

# Los decimales vienen de las estaciones climatológicas (tipo 02/03 en los
# dígitos 7-8 del código), que leen a 0,1 mm. No son error.
decimales <- d[P > 0, .(pct_decimales = round(100 * mean(P != round(P)), 1), dias_lluvia = .N),
               by = .(tipo = substr(codigo, 7, 8))]
print(tabla(decimales, "2_decimales_por_tipo"))

# Redondeo de los observadores: cuántos días con lluvia (en mm enteros) terminan en cada
# dígito. El 5 y el 0 sobresalen frente a sus vecinos. No cambia la curva anual, así que
# queda solo como tabla de respaldo.
digitos <- d[P > 0 & P == round(P), .(dias = .N), by = .(ultimo_digito = P %% 10)][order(ultimo_digito)]
print(tabla(digitos, "2_ultimo_digito"))

# QUÉ DICEN LOS DATOS (H2)
#   Las cifras y su interpretación están en EXPLICACION_DATOS.md, sección 2 (corrida del
#   19-nov-2026, con las 98 estaciones de la regla de la OMM); las tablas de esta
#   sección las imprimen las líneas de arriba.


# ---- 3. DISPONIBILIDAD: REJILLA, COBERTURA Y HUECOS --------------------------
titulo("3. Disponibilidad")
# Dominio común 1–365: se quita el 29 de febrero (D-14, regla R1)
d <- d[!(mes == 2 & dia_mes == 29)]
d[, dia := yday(fecha) - (anio %% 4 == 0 & mes > 2)]
stopifnot(d[, all(dia %between% c(1, 365))])
setkey(d, codigo, anio, dia)

# Rejilla completa estación × año × día de todo el periodo
g <- CJ(codigo = valle$codigo, anio = ANIOS, dia = 1:365)
g <- d[, .(codigo, anio, dia, P)][g, on = .(codigo, anio, dia)]
op <- g[!is.na(P), .(ini = min(anio * 1000 + dia), fin = max(anio * 1000 + dia)), by = codigo]
g <- op[g, on = "codigo"]
g[, t := anio * 1000 + dia]
g[, estado := fifelse(!is.na(P), "observado", fifelse(t < ini | t > fin, "fuera de operación", "hueco interno"))]
g[, c("ini", "fin", "t") := NULL]
print(g[, .N, by = estado][, pct := round(100 * N / sum(N), 1)][])

cob <- g[, .(obs = sum(estado == "observado"), hueco_int = sum(estado == "hueco interno"),
             hueco_max = hueco_mas_largo(estado == "observado")), by = .(codigo, anio)]
cob[, pct := 100 * obs / 365]
cob <- valle[, .(codigo, nombre, vertiente)][cob, on = "codigo"]
tabla(cob, "3_cobertura_estacion_anio")

# Cuatro niveles en vez de una escala continua: así se ve de inmediato qué años sirven.
NIV_COB <- c("completo (100 %)", "casi completo (90–99 %)", "incompleto (menos de 90 %)", "sin datos: no operaba")
cob[, nivel := factor(fcase(pct == 0, NIV_COB[4], pct < 90, NIV_COB[3], pct < 100, NIV_COB[2], default = NIV_COB[1]), levels = NIV_COB)]
print(cob[, .N, by = nivel][, pct := round(100 * N / sum(N), 1)][order(nivel)])
EJE_PARTICION <- scale_x_continuous(breaks = c(1981, 1985, 1990, 1995, 2000, 2003), expand = c(0, 0),
                                    sec.axis = dup_axis(breaks = c(1988, 1998, 2002), name = NULL,
                                                        labels = c("ENTRENAMIENTO\n1981–1995", "VALIDACIÓN\n1996–2000", "PRUEBA\n2001–2003")))
orden <- cob[, .(m = mean(pct)), by = codigo][order(m), codigo]
grafico(ggplot(cob, aes(anio, factor(codigo, levels = orden), fill = nivel)) + geom_tile() +
          facet_grid(ETQ_VERT[vertiente] ~ ., scales = "free_y", space = "free_y", labeller = label_wrap_gen(16)) +
          scale_fill_manual(values = setNames(c("#1A9850", "#A6D96A", "#FDAE61", "grey85"), NIV_COB), name = "Días con dato\nen ese año") +
          geom_vline(xintercept = c(1995.5, 2000.5), colour = "black", linetype = 2) + EJE_PARTICION +
          labs(x = NULL, y = "Cada fila es una estación (arriba, las que tienen más datos)",
               title = "¿Qué estaciones tienen datos y en qué años? (las 220 del Valle con datos)",
               subtitle = "Cada casilla es un año de una estación. El nombre de cada estación incluida está en la figura siguiente.") +
          theme(axis.text.y = element_blank(), panel.grid = element_blank(), strip.text.y = element_text(angle = 0, face = "bold")),
        "3_cobertura_heatmap", alto = 8.5)
grafico(ggplot(cob[pct > 0], aes(pct)) + geom_histogram(binwidth = 2) +
          labs(x = "% de días con dato en el año", y = "Estación-años", title = "La cobertura es bimodal"),
        "3_cobertura_histograma")
umbrales <- c(50, 70, 80, 85, 90, 95, 99, 100)
print(cob[pct > 0, lapply(setNames(umbrales, paste0(">=", umbrales)), function(u) sum(pct >= u))])

grafico(ggplot(melt(cob[, .(`con algún dato` = sum(pct > 0), `≥90 % de días` = sum(pct >= 90)), by = anio],
                    id.vars = "anio"), aes(anio, value, colour = variable)) +
          geom_line() + geom_point() + labs(y = "Estaciones", x = NULL, colour = NULL, title = "Estaciones activas por año"),
        "3_estaciones_por_anio")

# Huecos internos (dentro del periodo de operación de cada estación)
rachas <- function(x) { r <- rle(x); data.table(valor = r$values, largo = r$lengths) }
hue <- g[estado != "fuera de operación", rachas(estado == "hueco interno"), by = codigo][valor == TRUE]
print(hue[, .(huecos = .N, dias = sum(largo), p50 = median(largo), p90 = quantile(largo, .9),
              max = max(largo), pct_de_1_dia = round(100 * mean(largo == 1), 1))])
grafico(ggplot(hue, aes(largo)) + geom_histogram(bins = 40) + scale_x_log10() +
          labs(x = "Largo del hueco (días, log)", y = "Huecos", title = "Largo de los huecos internos"),
        "3_huecos_largo")

# Años válidos de entrenamiento con el criterio provisional (solo cobertura)
e90 <- cob[anio <= 1995 & pct >= 90, .N, by = codigo]
cat("Estaciones con ≥7 / ≥10 años de entrenamiento con ≥90 % de días:", e90[N >= 7, .N], "/", e90[N >= 10, .N], "\n")

# QUÉ DICEN LOS DATOS (H3, H8)
#   Las cifras y su interpretación están en EXPLICACION_DATOS.md, sección 3 (corrida del
#   19-nov-2026, con las 98 estaciones de la regla de la OMM); las tablas de esta
#   sección las imprimen las líneas de arriba.


# ---- 4. VALORES: DISTRIBUCIÓN, EXTREMOS Y VECINAS ----------------------------
titulo("4. Valores")
print(d[, .(dias = .N, pct_ceros = round(100 * mean(P == 0), 1), p50_lluvia = median(P[P > 0]),
            p90 = quantile(P[P > 0], .9), p99 = quantile(P[P > 0], .99), p999 = quantile(P[P > 0], .999), max = max(P))])

# El p99 de cada estación es un umbral de las reglas R2 y R3: se calcula SOLO con
# entrenamiento (1981–1995) para no usar validación ni prueba (D-28). Si la
# estación tiene menos de 100 días con lluvia en entrenamiento, se usa todo su registro.
p99 <- d[, .(p99_todo = quantile(P[P > 0], .99), max = max(P), pct_ceros = 100 * mean(P == 0)), by = codigo]
p99 <- d[P > 0 & anio <= 1995, .(p99 = quantile(P, .99), dias_lluvia_ent = .N), by = codigo][p99, on = "codigo"]
p99[, p99 := fifelse(is.na(p99) | dias_lluvia_ent < 100, p99_todo, p99)]
cat("Estaciones con p99 de entrenamiento:", p99[!is.na(dias_lluvia_ent) & dias_lluvia_ent >= 100, .N], "de", nrow(p99),
    "| diferencia mediana con el p99 de todo el registro:", round(p99[, median(abs(p99 - p99_todo))], 1), "mm\n")
p99 <- valle[, .(codigo, nombre, vertiente, altitud)][p99, on = "codigo"]
p99[, razon_max_p99 := max / p99]
tabla(p99, "4_valores_por_estacion")
print(p99[order(-razon_max_p99)][1:10])

# ¿Los extremos son tormentas o errores? Si llovió de verdad, alguna vecina a
# ≤25 km debería registrar algo grande ese día, el anterior o el siguiente.
d <- p99[, .(codigo, p99)][d, on = "codigo"]
setkey(d, codigo, fecha)
extremos <- d[P > 3 * p99 & P >= 150, .(codigo, fecha, P, p99)]
extremos$vecinas <- mapply(function(cod, f) uniqueN(d[codigo %in% vecinas(cod) & fecha %between% c(f - 1L, f + 1L), codigo]),
                           extremos$codigo, extremos$fecha)
extremos$vecina_max <- mapply(function(cod, f) {
  v <- d[codigo %in% vecinas(cod) & fecha %between% c(f - 1L, f + 1L), P]
  if (length(v)) max(v) else NA_real_
}, extremos$codigo, extremos$fecha)
extremos <- valle[, .(codigo, nombre)][extremos, on = "codigo"]
print(tabla(extremos, "4_extremos_vs_vecinas"))

grafico(ggplot(valle[, .(codigo, vertiente)][d[P > 0], on = "codigo"], aes(P, colour = ETQ_VERT[vertiente])) +
          stat_ecdf() + scale_x_log10() + scale_colour_manual(values = COL_VERT) +
          labs(x = "Lluvia diaria (mm, días con lluvia, log)", y = "Fracción de días con esa lluvia o menos", colour = NULL,
               title = "Distribución de la lluvia diaria por zona"),
        "4_distribucion_lluvia")

# QUÉ DICEN LOS DATOS (H4)
#   Las cifras y su interpretación están en EXPLICACION_DATOS.md, sección 4 (corrida del
#   19-nov-2026, con las 98 estaciones de la regla de la OMM); las tablas de esta
#   sección las imprimen las líneas de arriba.


# ---- 5. PATRONES SOSPECHOSOS -------------------------------------------------
titulo("5. Patrones sospechosos")

# 5a. Totales mensuales anotados en un solo día: el mes tiene un solo día con
# lluvia y ese valor supera el p99 de la estación.
dias_mes <- c(31L, 28L, 31L, 30L, 31L, 30L, 31L, 31L, 30L, 31L, 30L, 31L)
mesd <- d[, .(dias_obs = .N, dias_lluvia = sum(P > 0), max = max(P), p99 = p99[1]), by = .(codigo, anio, mes)]
meses_sosp <- mesd[dias_lluvia == 1 & max > p99]
meses_sosp$p_lluvia_vecinas <- mapply(function(cod, a, m) {
  v <- d[codigo %in% vecinas(cod) & anio == a & mes == m, P]
  if (length(v)) mean(v > 0) else NA_real_
}, meses_sosp$codigo, meses_sosp$anio, meses_sosp$mes)
meses_sosp[, mes_vacio := dias_obs < dias_mes[mes] / 2]
meses_sosp <- valle[, .(codigo, nombre)][meses_sosp, on = "codigo"]
print(tabla(meses_sosp, "5a_meses_un_dia_con_lluvia"))

# 5b. Rachas de ceros: ¿sequía real o pluviómetro sin leer anotado como 0?
# Se contrasta con las vecinas a ≤25 km esos mismos días. (Contrastar con toda
# la vertiente no sirve: la Pacífica mezcla la costa, >7 000 mm/año, con el
# cañón seco del Dagua, ~800 mm/año.)
L <- dcast(g, anio + dia ~ codigo, value.var = "P")
clave_dia <- L$anio * 1000 + L$dia
L <- as.matrix(L[, -(1:2)])[, valle$codigo] > 0
gz <- g[estado != "fuera de operación"]
gz[, es0 := !is.na(P) & P == 0][, id := rleid(es0), by = codigo]
rz <- gz[es0 == TRUE, .(ini = min(anio * 1000 + dia), fin = max(anio * 1000 + dia), largo = .N), by = .(codigo, id)][largo >= 15]
contraste <- function(cod, a, b) {
  j <- vecinas(cod)
  if (!length(j)) return(c(NA_real_, 0))
  c(mean(L[clave_dia >= a & clave_dia <= b, j, drop = FALSE], na.rm = TRUE), length(j))
}
cc <- t(mapply(contraste, rz$codigo, rz$ini, rz$fin))
rz[, `:=`(p_lluvia_vecinas = cc[, 1], n_vecinas = cc[, 2])]
a_fecha <- function(k) as.IDate(sprintf("%d-01-01", k %/% 1000)) + (k %% 1000) - 1L + ((k %/% 1000) %% 4 == 0 & k %% 1000 > 59)
rz[, `:=`(ini = a_fecha(ini), fin = a_fecha(fin))]
rz <- valle[, .(codigo, nombre, vertiente)][rz, on = "codigo"]
tabla(rz, "5b_rachas_ceros_15d")
print(rz[, .(rachas = .N, de_30_o_mas = sum(largo >= 30), sin_vecinas = sum(is.na(p_lluvia_vecinas)))])
grafico(ggplot(rz[!is.na(p_lluvia_vecinas)], aes(largo, p_lluvia_vecinas, colour = ETQ_VERT[vertiente])) +
          geom_point(alpha = .5) + scale_x_log10() + scale_colour_manual(values = COL_VERT) +
          labs(x = "Largo de la racha de ceros (días, log)", y = "Fracción de días-vecina con lluvia (≤25 km)", colour = NULL,
               title = "Rachas ≥15 días de ceros: ¿las vecinas también estaban secas?"),
        "5b_rachas_ceros")

# Prueba de cada racha de ≥30 días en cero (R4). Si la estación de verdad estuvo seca, sus
# vecinas también debieron estar secas o casi: es la regla de periodos secos de GSDR-QC
# (Lewis et al. 2021, Environ. Model. Softw. 144:105169), que pide que VARIAS vecinas la
# contradigan. Para cada una de las 3 vecinas más cercanas a ≤25 km con datos en ≥80 % de
# la racha se mide:
#   k = días en que llovió la vecina durante la racha;
#   q = P(llueve en la estación | llueve en la vecina), con los días de entrenamiento
#       (1981–1995, sin fuga) del mes central de la racha y sus dos meses vecinos: q de la
#       temporada, como la ventana estacional de Durre et al. (2010); exige ≥30 días;
#   p = (1 - q)^k: probabilidad de que la estación marque 0 todos esos días.
# Una vecina CONTRADICE la racha si p < 0,05 / número de rachas (Bonferroni). La racha pasa
# a faltante si la contradicen al menos 2 vecinas (mayoría): una sola vecina en desacuerdo
# puede ser un aguacero local, y una sequía real (El Niño) deja secas también a las vecinas.
# Supone días independientes; como la lluvia viene en rachas, p exagera lo raro, y por eso
# el umbral es exigente. Antes R4 decidía con la vecina más cercana y q de todo el año
# (D-36); la comparación queda en 5b_rachas_ceros_versiones.csv.
R4_LARGO <- 30L; R4_VECINAS <- 3L; R4_MAYORIA <- 2L; R4_MIN_DIAS_Q <- 30L
ll_ent <- d[anio <= 1995 & !is.na(P), .(codigo, fecha, mes = month(fecha), ll = P > 0)]; setkey(ll_ent, codigo, fecha)
d_obs <- d[!is.na(P), .(codigo, fecha, P)]; setkey(d_obs, codigo, fecha)
r4 <- rz[largo >= R4_LARGO][, i := .I]
R4_UMBRAL <- 0.05 / nrow(r4)
r4_vec <- rbindlist(lapply(seq_len(nrow(r4)), function(i) {
  cod <- r4$codigo[i]; desde <- r4$ini[i]; hasta <- r4$fin[i]   # copias locales: dentro de d_obs[...] "codigo" sería su columna
  mc <- month(desde + (hasta - desde) %/% 2L); meses <- (mc + -2:0) %% 12L + 1L   # mes central y sus dos vecinos
  out <- list()
  for (v in names(sort(DK[cod, vecinas(cod)]))) {
    x <- d_obs[.(v)][fecha %between% c(desde, hasta)]
    if (nrow(x) < .8 * r4$largo[i]) next
    a <- merge(ll_ent[.(cod)], ll_ent[.(v)], by = "fecha")[ll.y == TRUE]
    q_t <- if (a[mes.x %in% meses, .N] >= R4_MIN_DIAS_Q) a[mes.x %in% meses, mean(ll.x)] else NA_real_
    out[[length(out) + 1L]] <- data.table(i, orden = length(out) + 1L, vecina = v, km = round(DK[cod, v], 1), dias_vecina = nrow(x),
                                          dias_lluvia_vecina = sum(x$P > 0), q = q_t, q_anual = a[, mean(ll.x)])
    if (length(out) == R4_VECINAS) break
  }
  rbindlist(out)
}))
r4_vec[, `:=`(p_azar = (1 - q)^dias_lluvia_vecina, nombre_vecina = valle[match(vecina, codigo), nombre])]
r4_vec[, contradice := !is.na(p_azar) & p_azar < R4_UMBRAL]
r4 <- r4_vec[, .(n_vecinas_r4 = .N, vecinas_contradicen = sum(contradice),
                 detalle_vecinas = paste(sprintf("%s (%.0f km): llovió %d días, q %.2f", nombre_vecina, km, dias_lluvia_vecina, q), collapse = "; ")),
             by = i][r4, on = "i"]
r4[is.na(n_vecinas_r4), `:=`(n_vecinas_r4 = 0L, vecinas_contradicen = 0L)]
# la vecina más cercana evaluada, para las tablas y la gráfica
r4 <- r4_vec[orden == 1L, .(i, vecina, nombre_vecina, km, dias_vecina, dias_lluvia_vecina, q, q_anual, p_azar)][r4, on = "i"]
r4[, falso := vecinas_contradicen >= R4_MAYORIA]
tabla(r4[order(-vecinas_contradicen, p_azar)], "5b_rachas_ceros_prueba"); tabla(r4_vec, "5b_rachas_ceros_vecinas")
cat(sprintf("Rachas de ≥%d días en 0: %d (%d con menos de 2 vecinas con datos). Umbral 0,05/%d = %.1e → %d rachas falsas (%d días) en %d estaciones\n",
            R4_LARGO, nrow(r4), r4[n_vecinas_r4 < R4_MAYORIA, .N], nrow(r4), R4_UMBRAL, r4[falso == TRUE, .N], r4[falso == TRUE, sum(largo)],
            r4[falso == TRUE, uniqueN(codigo)]))
print(r4[falso == TRUE][order(ini), .(nombre, ini, fin, largo, vecinas_contradicen, detalle_vecinas)])
# Versiones de la regla (evidencia de la elección, D-61)
r4_versiones <- rbind(
  r4_vec[orden == 1L, .(version = "Vecina más cercana, q de todo el año (R4 anterior, D-36)", rachas = sum((1 - q_anual)^dias_lluvia_vecina < R4_UMBRAL, na.rm = TRUE))],
  r4_vec[orden == 1L, .(version = "Vecina más cercana, q de la temporada", rachas = sum(contradice))],
  r4[, .(version = "Al menos 2 de las 3 vecinas, q de la temporada (R4 vigente)", rachas = sum(falso))],
  r4[, .(version = "Las 3 vecinas, q de la temporada", rachas = sum(vecinas_contradicen == R4_VECINAS))])
print(tabla(r4_versiones, "5b_rachas_ceros_versiones"))
grafico(ggplot(r4_vec[!is.na(p_azar)], aes(dias_lluvia_vecina, pmax(p_azar, 1e-12), colour = i %in% r4[falso == TRUE, i])) + geom_point(alpha = .6) +
          geom_hline(yintercept = R4_UMBRAL, linetype = 2) + scale_y_log10() +
          scale_colour_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey50"), labels = c(`TRUE` = "Racha falsa → faltante", `FALSE` = "Se mantiene"), name = NULL) +
          labs(x = "Días con lluvia en la vecina durante la racha", y = "p = (1 - q)^k (log)",
               title = "R4: rachas de ≥30 días en 0 frente a sus 3 vecinas más cercanas",
               subtitle = sprintf("Un punto por racha y vecina. Bajo la línea (0,05 / %d) la vecina contradice la racha; es falsa si la contradicen ≥2", nrow(r4))),
        "5b_rachas_ceros_prueba")

# Inventario de escenarios de R4. Una racha larga en 0 puede ser una sequía real (temporada
# seca, El Niño) o un pluviómetro sin leer. Se cruzan la temporada (seca o lluviosa según el
# ciclo de la PROPIA estación en entrenamiento: mes por debajo o por encima de su media
# mensual) con la fase ENSO del año en que empieza la racha, y se ve qué fracción vuelve
# faltante R4 en cada escenario. R4 no usa la fase ENSO: decide con las vecinas, que en una
# sequía real también están secas (regla de periodos secos de GSDR-QC, Lewis et al. 2021;
# ninguna de las fuentes revisadas condiciona este control a la fase ENSO).
clim_mes <- d[anio <= 1995 & !is.na(P), .(m = mean(P)), by = .(codigo, mes = month(fecha))][, seco := m < mean(m), by = codigo]
inv_r4 <- copy(r4)[, `:=`(mes = month(ini + (fin - ini) %/% 2L), anio = year(ini))]
inv_r4 <- clim_mes[, .(codigo, mes, seco)][inv_r4, on = .(codigo, mes)]
inv_r4 <- fread("datos_web/oni_anual_fase.csv")[, .(anio, fase)][inv_r4, on = "anio"]
inv_r4 <- inv_r4[, .(rachas = .N, a_faltante_R4 = sum(falso), pct_a_faltante = round(100 * mean(falso), 1), largo_mediano = as.double(median(largo)),
                     dias_lluvia_vecina_mediana = as.double(median(dias_lluvia_vecina, na.rm = TRUE))),
                 by = .(temporada = fcase(seco, "seca", !seco, "lluviosa", default = "sin clima de entrenamiento"), fase)][order(temporada, fase)]
print(tabla(inv_r4, "5b_rachas_ceros_escenarios"))

# 5c. Años idénticos entre estaciones (copias)
huella <- d[, .(n = .N, h = paste(P, collapse = ",")), by = .(codigo, anio)][n >= 300]
copias <- huella[, .(estaciones = paste(codigo, collapse = " / "), k = .N), by = .(anio, h)][k > 1, .(anio, estaciones, k)]
cat("Años idénticos entre estaciones:", nrow(copias), "\n")

# 5d. Valores positivos repetidos varios días seguidos
dr <- d[order(codigo, fecha)]
dr[, id := rleid(P, fecha - shift(fecha, fill = fecha[1] - 1L) != 1), by = codigo]
repetidos <- dr[P > 0, .(valor = P[1], ini = min(fecha), largo = .N), by = .(codigo, id)][largo >= 4]
print(repetidos[order(-largo)][1:5])
rm(dr)

# QUÉ DICEN LOS DATOS (H5–H7)
#   Las cifras y su interpretación están en EXPLICACION_DATOS.md, sección 5 (corrida del
#   19-nov-2026, con las 98 estaciones de la regla de la OMM); las tablas de esta
#   sección las imprimen las líneas de arriba.


# ---- 6. CLIMATOLOGÍA Y COHERENCIA ESPACIAL -----------------------------------
titulo("6. Climatología y coherencia espacial")
buenos <- cob[pct >= 90, .(codigo, anio)]   # solo para describir, no es regla
anual <- d[buenos, on = .(codigo, anio), nomatch = 0][, .(P_anual = mean(P) * 365), by = .(codigo, anio)]
clim <- anual[, .(P_media = mean(P_anual), cv = sd(P_anual) / mean(P_anual), n = .N), by = codigo]
clim <- valle[, .(codigo, nombre, vertiente, lon, lat, altitud)][clim, on = "codigo"]
tabla(clim, "6_climatologia_estacion")
cat("Lluvia media anual: de", round(min(clim$P_media)), "a", round(max(clim$P_media)), "mm\n")

grafico(ggplot(clim, aes(lon, lat, colour = P_media)) + geom_point(size = 2.5) + coord_equal() +
          scale_colour_viridis_c(trans = "log10", name = "mm/año") + labs(title = "Lluvia media anual por estación"),
        "6_mapa_lluvia_anual", ancho = 7, alto = 7)
grafico(ggplot(clim[!is.na(altitud)], aes(altitud, P_media, colour = ETQ_VERT[vertiente])) + geom_point() +
          scale_y_log10() + geom_smooth(se = FALSE, method = "loess", formula = y ~ x) + scale_colour_manual(values = COL_VERT) +
          labs(x = "Altitud de la estación (m)", y = "Lluvia media anual (mm, escala log)", colour = NULL,
               title = "La lluvia anual depende de la altitud, y distinto en cada zona",
               subtitle = "Cada punto es una estación; la línea es la tendencia de su zona") +
          theme(legend.position = "bottom"),
        "6_lluvia_altitud")

men <- d[buenos, on = .(codigo, anio), nomatch = 0][, .(P = mean(P)), by = .(codigo, mes)]
men <- valle[, .(codigo, vertiente)][men, on = "codigo"]
grafico(ggplot(men, aes(factor(mes), P)) + geom_boxplot(outlier.size = .5, fill = "grey92") +
          facet_wrap(~vertiente, scales = "free_y", labeller = as_labeller(ETQ_VERT)) +
          scale_x_discrete(labels = c("E", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) +
          labs(x = "Mes", y = "Lluvia media diaria (mm/día)",
               title = "Ciclo anual por zona: dos temporadas de lluvia en la zona andina",
               subtitle = "Cada caja resume las estaciones de la zona en ese mes. Ojo: el eje de la zona pacífica llega a más del doble."),
        "6_ciclo_anual")

oni <- fread("datos_web/oni_anual_fase.csv")
reg <- anual[valle[, .(codigo, vertiente)], on = "codigo", nomatch = 0]
reg[, anom := P_anual / clim$P_media[match(codigo, clim$codigo)] - 1]
reg <- reg[, .(anom = median(anom)), by = .(vertiente, anio)]
reg <- oni[reg, on = "anio"]
print(reg[, .(r_ONI_anomalia = round(cor(ONI_anual, anom), 2)), by = vertiente])
grafico(ggplot(reg, aes(anio, 100 * anom, fill = fase)) + geom_col() + facet_wrap(~vertiente, ncol = 1, labeller = as_labeller(ETQ_VERT)) +
          labs(x = NULL, y = "% respecto a la media de la estación (mediana)", title = "Anomalía anual y fase ENSO"),
        "6_anomalia_anual_ENSO", alto = 6)

# Coherencia espacial: correlación de anomalías mensuales de log(1+P) vs. distancia
mm <- d[, .(P = mean(P), n = .N), by = .(codigo, anio, mes)][n >= 25]
mm[, anom := log1p(P) - mean(log1p(P)), by = .(codigo, mes)]
W <- dcast(mm, anio + mes ~ codigo, value.var = "anom")
M <- as.matrix(W[, -(1:2)])
C <- cor(M, use = "pairwise.complete.obs"); N <- crossprod(!is.na(M))
D <- DK[colnames(M), colnames(M)]
ut <- upper.tri(C) & N >= 60
vert <- valle$vertiente[match(colnames(M), valle$codigo)]
pares <- data.table(r = C[ut], km = D[ut], misma_vertiente = outer(vert, vert, "==")[ut])
print(pares[, .(r_mediana = round(median(r), 2), pares = .N), by = .(rango_km = cut(km, c(0, 10, 20, 40, 80, Inf)))][order(rango_km)])
grafico(ggplot(pares, aes(km, r, colour = misma_vertiente)) + geom_point(alpha = .15, size = .6) +
          scale_colour_discrete(labels = c(`FALSE` = "pares de zonas distintas", `TRUE` = "pares de la misma zona"), name = NULL) +
          geom_smooth(se = FALSE, method = "gam", formula = y ~ s(x, bs = "cs")) +
          labs(x = "Distancia (km)", y = "r de anomalías mensuales", title = "Las estaciones cercanas se parecen más"),
        "6_coherencia_distancia")
coherencia <- rbindlist(lapply(seq_len(ncol(M)), function(i) {
  j <- which(D[i, ] <= RADIO_VEC & N[i, ] >= 60 & seq_len(ncol(M)) != i)
  if (!length(j)) return(NULL)
  data.table(codigo = colnames(M)[i], vecinas = length(j), r_max = max(C[i, j]))
}))
coherencia <- valle[, .(codigo, nombre)][coherencia, on = "codigo"]
tabla(coherencia, "6_coherencia_vecinas")
print(coherencia[order(r_max)][1:5])

# QUÉ DICEN LOS DATOS (H9)
#   Las cifras y su interpretación están en EXPLICACION_DATOS.md, sección 6 (corrida del
#   19-nov-2026, con las 98 estaciones de la regla de la OMM); las tablas de esta
#   sección las imprimen las líneas de arriba.
rm(L, gz, M, W, C, N, D, pares, huella); invisible(gc())


# ---- 7. EXPERIMENTO DE HUECOS ------------------------------------------------
titulo("7. Experimento de huecos")
# Pregunta: la curva anual se ajusta solo con los días observados, así que un año
# incompleto igual produce una curva. ¿Desde cuándo deja de ser confiable?
# Diseño: años COMPLETOS de entrenamiento = verdad. Se les quitan días de tres
# maneras (A hueco continuo, B días sueltos, C patrón real de faltantes de otro
# año), se reajusta y se compara. Error = RMSE / variabilidad interanual de la
# estación (0,25 = el hueco mueve la curva un cuarto de lo que varía entre años).

# Base B-spline cúbica propia (Cox–de Boor), nodos equiespaciados en 1–365
nodos <- function(K, orden = 4L) {
  interiores <- seq(1, 365, length.out = K - orden + 2L)
  c(rep(1, orden), interiores[-c(1L, length(interiores))], rep(365, orden))
}
base_bspline <- function(t, K, orden = 4L) {
  tau <- nodos(K, orden)
  B <- outer(t, seq_len(length(tau) - 1L), function(x, j) (x >= tau[j]) & (x < tau[j + 1L])) * 1
  B[t == max(tau), max(which(diff(tau) > 0))] <- 1
  for (m in 2:orden) {
    S <- matrix(0, length(t), length(tau) - m)
    for (j in seq_len(ncol(S))) {
      d1 <- tau[j + m - 1L] - tau[j]; d2 <- tau[j + m] - tau[j + 1L]
      if (d1 > 0) S[, j] <- S[, j] + (t - tau[j]) / d1 * B[, j]
      if (d2 > 0) S[, j] <- S[, j] + (tau[j + m] - t) / d2 * B[, j + 1L]
    }
    B <- S
  }
  B
}
stopifnot(all(abs(rowSums(base_bspline(1:365, 8)) - 1) < 1e-12))   # partición de la unidad

# Mínimos cuadrados con los días observados; NA si la base no es identificable
ajustar <- function(y, obs, Phi) {
  P <- Phi[obs, , drop = FALSE]; A <- crossprod(P)
  if (rcond(A) < 1e-10) return(rep(NA_real_, 365))
  drop(Phi %*% solve(A, crossprod(P, y[obs])))
}

d[, clave := paste(codigo, anio, sep = "_")]
Ytodo <- dcast(d, clave ~ dia, value.var = "P")
claves <- Ytodo$clave; Ytodo <- as.matrix(Ytodo[, -1]); rownames(Ytodo) <- claves
info <- data.table(clave = claves, codigo = sub("_.*", "", claves), anio = as.integer(sub(".*_", "", claves)))
info[, n_obs := rowSums(!is.na(Ytodo))]
completos <- info[n_obs == 365 & anio <= 1995]
Ycomp <- Ytodo[completos$clave, ]
cat("Estación-años completos de entrenamiento:", nrow(completos), "de", uniqueN(completos$codigo), "estaciones\n")

Ks <- c(8, 12)   # K vigente y una sensibilidad
ref <- lapply(Ks, function(K) {
  Phi <- base_bspline(1:365, K)
  X <- t(apply(Ycomp, 1, ajustar, obs = rep(TRUE, 365), Phi = Phi))
  sd <- data.table(codigo = completos$codigo, i = seq_len(nrow(X)))[
    , .(sd = if (.N >= 3) sqrt(mean(sweep(X[i, , drop = FALSE], 2, colMeans(X[i, , drop = FALSE]))^2)) else NA_real_), by = codigo]
  list(Phi = Phi, X = X, sd = sd)
}); names(ref) <- Ks
err_rel <- function(i, obs, K) {
  r <- ref[[as.character(K)]]
  rmse <- sqrt(mean((ajustar(Ycomp[i, ], obs, r$Phi) - r$X[i, ])^2))
  c(rmse = rmse, rel = rmse / r$sd[codigo == completos$codigo[i], sd])
}
usar <- which(completos$codigo %in% ref[["8"]]$sd[!is.na(sd), codigo])   # estaciones con ≥3 años completos

# A. Un hueco continuo de largo L
Ls <- c(1, 3, 5, 7, 10, 15, 20, 30, 45, 60, 75, 90)
expA <- rbindlist(lapply(Ks, function(K) rbindlist(lapply(Ls, function(Lh) rbindlist(lapply(usar, function(i) {
  ini <- sample.int(365 - Lh + 1, 1); obs <- rep(TRUE, 365); obs[ini:(ini + Lh - 1)] <- FALSE
  e <- err_rel(i, obs, K); data.table(K, L = Lh, ini, rmse = e[["rmse"]], rel = e[["rel"]])
}))))))
tabla(expA, "7A_hueco_continuo")
print(dcast(expA[, .(p50 = round(median(rel, na.rm = TRUE), 2)), by = .(K, L)], L ~ K, value.var = "p50"))
grafico(ggplot(expA[!is.na(rel)], aes(factor(L), rel)) + geom_boxplot(outlier.size = .3) +
          facet_wrap(~paste("K =", K)) + geom_hline(yintercept = .25, colour = "red", linetype = 2) + scale_y_log10() +
          labs(x = "Largo del hueco (días)", y = "Error / variabilidad interanual (log)",
               title = "Un solo hueco continuo: el error crece con el largo", subtitle = "Línea roja: 0,25"),
        "7A_hueco_continuo")

# B. Días sueltos al azar (fracción f)
fs <- c(.01, .02, .05, .10, .15, .20, .30, .40, .50)
expB <- rbindlist(lapply(Ks, function(K) rbindlist(lapply(fs, function(f) rbindlist(lapply(usar, function(i) {
  obs <- rep(TRUE, 365); obs[sample.int(365, round(f * 365))] <- FALSE
  e <- err_rel(i, obs, K); data.table(K, f, rmse = e[["rmse"]], rel = e[["rel"]])
}))))))
tabla(expB, "7B_dias_sueltos")
print(dcast(expB[, .(p50 = round(median(rel, na.rm = TRUE), 2)), by = .(K, pct_faltante = 100 * f)], pct_faltante ~ K, value.var = "p50"))

# C. Patrones reales: cada año incompleto (≥40 % de días) presta su máscara a 5 años completos
inc <- info[n_obs < 365 & n_obs >= 0.4 * 365]
mascaras <- !is.na(Ytodo[inc$clave, , drop = FALSE])
inc[, `:=`(pct = 100 * n_obs / 365, hueco_max = apply(mascaras, 1, hueco_mas_largo))]
expC <- rbindlist(lapply(Ks, function(K) rbindlist(lapply(seq_len(nrow(inc)), function(m) rbindlist(lapply(sample(usar, 5), function(i) {
  e <- err_rel(i, mascaras[m, ], K)
  data.table(K, mascara = inc$clave[m], pct = inc$pct[m], hueco_max = inc$hueco_max[m], rmse = e[["rmse"]], rel = e[["rel"]])
}))))))
tabla(expC, "7C_mascaras_reales")
c8 <- expC[K == 8 & !is.na(rel)]
cat("R² del log-error explicado por el hueco más largo:", round(summary(lm(log(rel) ~ log1p(hueco_max), c8))$r.squared, 2),
    "| por el % de cobertura:", round(summary(lm(log(rel) ~ pct, c8))$r.squared, 2), "\n")
grafico(ggplot(expC[!is.na(rel)], aes(hueco_max, rel, colour = pct)) + geom_point(alpha = .3, size = .8) +
          facet_wrap(~paste("K =", K)) + scale_x_log10() + scale_y_log10() +
          geom_hline(yintercept = .25, colour = "red", linetype = 2) + scale_colour_viridis_c(name = "% días") +
          labs(x = "Hueco más largo (días, log)", y = "Error relativo (log)",
               title = "Máscaras reales: el hueco más largo explica el error mejor que el % de cobertura"),
        "7C_mascaras_reales")

# ¿Cuánto mueve la curva cada regla de la literatura? Se comparan la regla de la OMM de
# 2017 (R5), la más estricta de 1989 (5 faltantes y 3 seguidos por mes, WMO 1989), la regla
# propia que usábamos antes (≥90 % de días y hueco ≤15) y no exigir nada (≥40 % de días).
# Para cada regla: cuántos años de entrenamiento conserva y qué error dejan las máscaras
# reales que la cumplen (experimento C, K = 8).
reglas_def <- list(
  "OMM 2017 (R5): mes con ≤10 faltantes y <5 seguidos" = function(o) anio_valido_omm(o),
  "OMM 1989: mes con ≤5 faltantes y ≤3 seguidos" = function(o) anio_valido_omm(o, falt = 5, seguid = 3),
  "Regla propia anterior: ≥90 % de días y hueco ≤15" = function(o) sum(o) >= .9 * 365 && hueco_mas_largo(o) <= 15,
  "Sin regla: cualquier año con ≥40 % de días" = function(o) sum(o) >= .4 * 365)
obs_ent <- !is.na(Ytodo[info[anio <= 1995, clave], , drop = FALSE])
masc_ok <- sapply(reglas_def, function(f) apply(mascaras, 1, f))   # máscara × regla
rownames(masc_ok) <- inc$clave
reglas_cmp <- rbindlist(lapply(names(reglas_def), function(nm) {
  e <- c8[masc_ok[mascara, nm], rel]
  data.table(regla = nm, estacion_anios_ent = sum(apply(obs_ent, 1, reglas_def[[nm]])), mascaras_que_cumplen = length(e),
             mediana_err = round(median(e), 2), p90_err = round(quantile(e, .9), 2), p99_err = round(quantile(e, .99), 2))
}))
print(tabla(reglas_cmp, "7_reglas_comparadas"))
reglas_cmp[, vj := c(-.9, 1.9, -.9, -.9)]   # la de 1989 va debajo para no encimarse con la de 2017
grafico(ggplot(reglas_cmp, aes(estacion_anios_ent, p99_err)) + geom_point(size = 2.5) +
          geom_text(aes(label = regla, vjust = vj), hjust = 0, nudge_x = 6, size = 3) + geom_hline(yintercept = .25, linetype = 2) +
          scale_y_log10() + scale_x_continuous(expand = expansion(mult = c(.05, .55))) +
          labs(x = "Estación-años de entrenamiento que conserva la regla", y = "p99 del error relativo (log)",
               title = "Reglas de la literatura: datos que conservan y error que dejan en la curva",
               subtitle = "Error = cuánto mueve el hueco la curva, sobre la variación entre años de la estación. Línea: 0,25"),
        "7_reglas_comparadas", ancho = 10, alto = 4.5)

# QUÉ DICEN LOS DATOS
#   Las cifras y su interpretación están en EXPLICACION_DATOS.md, sección 7 (corrida del
#   19-nov-2026, con las 98 estaciones de la regla de la OMM); las tablas de esta
#   sección las imprimen las líneas de arriba.
rm(Ytodo, Ycomp, mascaras, ref, expA, expB, expC, c8, obs_ent, masc_ok); invisible(gc())


# ---- 8. DEPURACIÓN (D-36) Y DATASET FINAL -------------------------------------
titulo("8. Depuración y dataset final")
# R1. Dominio 1–365 sin 29-feb (hecho en la sección 3).
# R2. Extremo sin respaldo → faltante: P > 3·p99, P ≥ 150 mm y ninguna vecina
#     a ≤25 km supera P/3 ese día, el anterior o el siguiente.
m_ext <- extremos[!is.na(vecina_max) & vecina_max < P / 3,
                  .(codigo, fecha, P, regla = "R2 extremo sin respaldo",
                    detalle = sprintf("p99 %.0f; vecina máx %.0f mm", p99, vecina_max))]

# R3. Total mensual en un día → mes entero faltante si el mes tiene menos de la
#     mitad de sus días observados o si las vecinas llovieron en ≥40 % de sus
#     días (entonces los ceros del resto del mes también son falsos).
meses_sosp[, marcar := mes_vacio | (!is.na(p_lluvia_vecinas) & p_lluvia_vecinas >= 0.4)]
m_mes <- d[meses_sosp[marcar == TRUE, .(codigo, anio, mes, dias_obs, p_lluvia_vecinas, mes_vacio)],
           on = .(codigo, anio, mes), nomatch = 0][
  , .(codigo, fecha, P, regla = "R3 total mensual en un día",
      detalle = sprintf("%d días observados; vecinas llovieron %.0f %% de los días%s",
                        dias_obs, 100 * p_lluvia_vecinas, fifelse(mes_vacio, "; mes casi vacío", "")))]

# R4. Racha de ≥30 días en 0 que al menos 2 de sus 3 vecinas más cercanas contradicen →
#     faltante (prueba de la sección 5b, umbral 0,05 / número de rachas).
m_cero <- r4[falso == TRUE][, d[codigo == .BY$codigo & fecha %between% c(ini, fin), .(fecha, P)], by = .(codigo, ini, fin, vecinas_contradicen, detalle_vecinas)][
  , .(codigo, fecha, P, regla = "R4 racha de ceros sin respaldo",
      detalle = sprintf("%d de 3 vecinas la contradicen: %s", vecinas_contradicen, detalle_vecinas))]
stopifnot(all(m_cero$P == 0), nrow(m_cero) == r4[falso == TRUE, sum(largo)])

marcas <- valle[, .(codigo, nombre)][rbind(m_r0, m_ext, m_mes, m_cero), on = "codigo"]
print(marcas[, .(dias = .N, estaciones = uniqueN(codigo)), by = regla])

# Aplicar: un faltante es una fila ausente / NA en la rejilla (R7 rellena después de R5)
limpio <- d[!marcas, on = .(codigo, fecha), .(codigo, anio, dia, fecha, P)]
anios_op <- unique(limpio[, .(codigo, anio)])
diario_todo <- limpio[anios_op[, .(dia = 1:365), by = .(codigo, anio)], on = .(codigo, anio, dia)]
diario_todo[is.na(fecha), fecha := as.IDate(sprintf("%d-01-01", anio)) + dia - 1L + (anio %% 4 == 0 & dia > 59)]
diario_todo[, particion := PARTICION(anio)]
setkey(diario_todo, codigo, anio, dia)
stopifnot(!anyDuplicated(diario_todo[, .(codigo, anio, dia)]), diario_todo[, all(month(fecha) * 100 + mday(fecha) != 229)])

# R5 (OMM 2017, 4.4.1). Año válido: en cada uno de sus 12 meses faltan ≤10 días y no hay
# 5 o más días faltantes seguidos. Un valor que el control de calidad (R0–R4) volvió
# faltante cuenta como faltante (OMM 2017, nota 2 de la sección 4.4).
anios <- diario_todo[, .(n_obs = sum(!is.na(P)), hueco_max = hueco_mas_largo(!is.na(P)), valido = anio_valido_omm(!is.na(P))),
                     by = .(codigo, anio, particion)]

# R7. Relleno de los días faltantes de los años válidos con kriging ordinario de las vecinas.
#   - Qué: solo los faltantes de años que pasaron R5 (≤10 por mes y <5 seguidos). Los años
#     que no pasan R5, las rachas de R4, los meses de R3 y los años sin operar NO se rellenan.
#   - Por qué con vecinas: para lluvia diaria, "values from other days are rarely useful"
#     (OMM 2017, 4.6); la OMM pone primero la interpolación espacial, solo para huecos de
#     ≤10 días por mes y sin sesgo (estima igual días secos y con lluvia). En el Valle la
#     lluvia del MISMO día en una vecina a ≤5 km se correlaciona 0,54 con la estación; la de
#     la propia estación el día anterior, 0,12 (EXPLICACION_DATOS.md, sección 8).
#   - Cómo: kriging ordinario (Cressie 1993), como la tesis del profesor y Terán-Chaves et al.
#     (2021, datos diarios del IDEAM). Cada vecina entra escalada por la razón entre la
#     lluvia media de la estación y la suya en ese mes (razón normal), porque vecinas a
#     ≤25 km llueven en promedio 1,3 veces una lo de la otra y el kriging ordinario supone
#     media constante. Hasta 10 vecinas a ≤50 km con dato ese día; negativo → 0.
#   - Sin fuga: el semivariograma y las medias mensuales salen de 1981–1995. El modelo del
#     semivariograma se elige por validación cruzada (menor ECM), como en la tesis.
#   - Nunca con el satélite: es lo que el modelo usa para predecir (sería circular).
#   - Un día sin ninguna vecina con dato queda faltante y lo cubre la curva.
titulo("8. R7: relleno por kriging")
R7_RADIO <- 50; R7_MAX_VEC <- 10L; R7_MIN_MEDIA <- 0.05   # mm/día: evita dividir por cero en meses secos
R7_EST <- valle$codigo; R7_T <- length(ANIOS) * 365L; R7_ENT <- seq_len(sum(ANIOS <= 1995) * 365L)
PM <- matrix(NA_real_, R7_T, length(R7_EST), dimnames = list(NULL, R7_EST))
PM[cbind((limpio$anio - ANIOS[1]) * 365L + limpio$dia, match(limpio$codigo, R7_EST))] <- limpio$P
R7_MES <- rep(MES_DIA, length(ANIOS))
# Lluvia media de cada estación en cada mes (mm/día), con entrenamiento; si en ese mes tiene
# <30 días de entrenamiento se usa todo su registro (como el p99 de R2 y R3)
med_ent <- sapply(R7_EST, function(s) tapply(PM[R7_ENT, s], R7_MES[R7_ENT], function(v) if (sum(!is.na(v)) >= 30) mean(v, na.rm = TRUE) else NA_real_))
med_todo <- sapply(R7_EST, function(s) tapply(PM[, s], R7_MES, mean, na.rm = TRUE))
R7_MEDIA <- ifelse(is.na(med_ent), med_todo, med_ent)   # 12 meses × estaciones
cat("Medias mensuales: de entrenamiento", sum(!is.na(med_ent)), "| de todo el registro", sum(is.na(med_ent) & !is.na(med_todo)), "\n")

# ¿Por qué con vecinas y no con los días de antes y después? Correlación de la lluvia del
# mismo día entre estaciones, por distancia, frente a la de la estación con su día anterior
pares_r7 <- which(upper.tri(DK[R7_EST, R7_EST]) & DK[R7_EST, R7_EST] <= 50, arr.ind = TRUE)
cor_r7 <- rbindlist(lapply(seq_len(nrow(pares_r7)), function(k) {
  a <- PM[R7_ENT, pares_r7[k, 1]]; b <- PM[R7_ENT, pares_r7[k, 2]]; ok <- !is.na(a) & !is.na(b)
  if (sum(ok) >= 365) data.table(km = DK[R7_EST, R7_EST][pares_r7[k, , drop = FALSE]], r = cor(a[ok], b[ok]),
                                 razon_medias = max(mean(a[ok]), mean(b[ok])) / max(min(mean(a[ok]), mean(b[ok])), R7_MIN_MEDIA)) }))
auto_r7 <- vapply(R7_EST, function(s) { a <- PM[R7_ENT, s]; x <- a[-length(a)]; y <- a[-1]; ok <- !is.na(x) & !is.na(y)
  if (sum(ok) >= 365) cor(x[ok], y[ok]) else NA_real_ }, 0)
por_que_r7 <- rbind(
  cor_r7[, .(tramo = cut(km, c(0, 5, 10, 15, 25, 50)), r = r)][, .(pares = .N, r_mediana = round(median(r), 2)), by = tramo][
    order(tramo), .(comparacion = paste("Vecina a", tramo, "km, mismo día"), pares, r_mediana)],
  data.table(comparacion = "La misma estación con su día anterior", pares = sum(!is.na(auto_r7)), r_mediana = round(median(auto_r7, na.rm = TRUE), 2)))
print(tabla(por_que_r7, "8_R7_por_que_vecinas"))
cat("Vecinas a ≤25 km: una llueve en promedio", round(cor_r7[km <= 25, median(razon_medias)], 2),
    "veces lo que la otra (mediana) → cada vecina se escala por la razón de medias mensuales\n")

# Semivariograma empírico de la lluvia diaria (entrenamiento), pares a ≤100 km en tramos de 5 km
DKv <- DK[R7_EST, R7_EST]; tramos <- seq(0, 100, by = 5); acum <- matrix(0, length(tramos) - 1L, 2L)
for (i in seq_along(R7_EST)) {
  j <- which(DKv[i, ] > 0 & DKv[i, ] <= 100); if (!length(j)) next
  Z <- (PM[R7_ENT, j, drop = FALSE] - PM[R7_ENT, i])^2
  n_c <- colSums(!is.na(Z)); m_c <- colSums(Z, na.rm = TRUE) / pmax(n_c, 1)
  b <- findInterval(DKv[i, j], tramos, left.open = TRUE); ok <- n_c >= 365 & b >= 1 & b < length(tramos)
  for (k in which(ok)) acum[b[k], ] <- acum[b[k], ] + c(m_c[k] / 2, 1)
}
vg_emp <- data.table(km = (tramos[-1] + tramos[-length(tramos)]) / 2, gamma = acum[, 1] / acum[, 2], pares = acum[, 2] / 2)[pares > 0]
print(tabla(vg_emp, "8_R7_semivariograma_empirico")[1:8])
# Modelos teóricos (pepita c0, meseta parcial c1, rango a en km); γ(0) = 0
VG <- list(
  exponencial = function(h, p) ifelse(h == 0, 0, p[1] + p[2] * (1 - exp(-h / p[3]))),
  esferico    = function(h, p) ifelse(h == 0, 0, p[1] + p[2] * ifelse(h < p[3], 1.5 * h / p[3] - 0.5 * (h / p[3])^3, 1)),
  gaussiano   = function(h, p) ifelse(h == 0, 0, p[1] + p[2] * (1 - exp(-(h / p[3])^2))),
  matern_1.5  = function(h, p) ifelse(h == 0, 0, p[1] + p[2] * (1 - (1 + sqrt(3) * h / p[3]) * exp(-sqrt(3) * h / p[3]))))
ajusta_vg <- function(f) {   # mínimos cuadrados ponderados por el número de pares; rango entre 1 y 300 km
  par_de <- function(z) c(exp(z[1]), exp(z[2]), 1 + 299 * plogis(z[3]))
  obj <- function(z) sum(vg_emp$pares * (vg_emp$gamma - f(vg_emp$km, par_de(z)))^2)
  par_de(optim(c(log(min(vg_emp$gamma)), log(diff(range(vg_emp$gamma)) + 1e-3), 0), obj, control = list(maxit = 5000))$par)
}
vg_par <- lapply(VG, ajusta_vg)
# Pesos del kriging ordinario para estimar en x0 con las vecinas v (devuelve pesos y varianza)
pesos_ok <- function(s, v, f, p) {
  n <- length(v); G <- f(DKv[v, v, drop = FALSE], p); diag(G) <- 0
  A <- rbind(cbind(G, 1), c(rep(1, n), 0)); b <- c(f(DKv[s, v], p), 1)
  sol <- tryCatch(solve(A, b), error = function(e) c(rep(1 / n, n), 0))
  list(w = sol[seq_len(n)], var = sum(sol[seq_len(n)] * b[seq_len(n)]) + sol[n + 1L])
}
# Una estación sin NINGÚN año válido (R5) no sirve de vecina: su cobertura no se puede comprobar y
# puede esconder acumulados. Ejemplo: La Perdición (Calima) anota unos 52 días al año, casi ninguno
# en 0 y con mediana de 130 mm: la lluvia de la semana en un solo día (Viney y Bates 2004; IDEAM:
# "solo se aceptan hasta cuatro (4) días acumulados"). Son 7 estaciones (D-64).
R7_CON_VALIDOS <- intersect(R7_EST, anios[valido == TRUE, unique(codigo)])
cat("R7: estaciones sin ningún año válido que no se usan como vecinas:", paste(valle[codigo %in% setdiff(R7_EST, R7_CON_VALIDOS), nombre], collapse = ", "), "\n")
vecinas_r7 <- lapply(R7_EST, function(s) { v <- setdiff(names(sort(DKv[s, ])), s); v[DKv[s, v] <= R7_RADIO & v %in% R7_CON_VALIDOS] }); names(vecinas_r7) <- R7_EST
estima_ok <- function(s, t, f = VG[[vg_mejor]], p = vg_par[[vg_mejor]], razon = TRUE) {
  v <- vecinas_r7[[s]]; v <- head(v[!is.na(PM[t, v])], R7_MAX_VEC)
  if (!length(v)) return(c(NA_real_, NA_real_, 0))
  m <- R7_MES[t]; z <- PM[t, v] * if (razon) R7_MEDIA[m, s] / pmax(R7_MEDIA[m, v], R7_MIN_MEDIA) else 1
  k <- pesos_ok(s, v, f, p); c(max(0, sum(k$w * z)), sqrt(max(k$var, 0)), length(v))
}
# Validación cruzada: se quita un dato observado de entrenamiento y se estima con las vecinas
set.seed(20261120)
cv_obs <- which(!is.na(PM[R7_ENT, ]), arr.ind = TRUE)
cv_obs <- cv_obs[sample(nrow(cv_obs), 4000), , drop = FALSE]
vg_cv <- rbindlist(lapply(names(VG), function(nm) {
  e <- vapply(seq_len(nrow(cv_obs)), function(r) estima_ok(R7_EST[cv_obs[r, 2]], cv_obs[r, 1], VG[[nm]], vg_par[[nm]])[1], 0)
  data.table(modelo = nm, pepita = vg_par[[nm]][1], meseta_parcial = vg_par[[nm]][2], rango_km = vg_par[[nm]][3],
             ECM_validacion_cruzada = mean((e - PM[cbind(cv_obs[, 1], cv_obs[, 2])])^2, na.rm = TRUE), datos = sum(!is.na(e)))
}))
vg_mejor <- vg_cv[which.min(ECM_validacion_cruzada), modelo]
print(tabla(vg_cv, "8_R7_semivariograma_modelos")); cat("Modelo elegido por validación cruzada:", vg_mejor, "\n")
vg_curva <- CJ(modelo = names(VG), km = seq(0, 100, by = 1))[, gamma := VG[[modelo[1]]](km, vg_par[[modelo[1]]]), by = modelo]
grafico(ggplot(vg_emp, aes(km, gamma)) + geom_point(aes(size = pares), colour = "grey30") +
          geom_line(data = vg_curva, aes(colour = modelo == vg_mejor, group = modelo)) +
          scale_colour_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey70"), labels = c(`TRUE` = paste("elegido:", vg_mejor), `FALSE` = "otros modelos"), name = NULL) +
          labs(x = "Distancia entre estaciones (km)", y = "Semivarianza (mm²/día²)", size = "Pares",
               title = "R7: semivariograma de la lluvia diaria (1981–1995)",
               subtitle = "Aun a 2,5 km la semivarianza es alta: la lluvia de un día varía mucho entre vecinas"),
        "8_R7_semivariograma")

# Evidencia de la elección (D-61): años completos de entrenamiento a los que se les quitan días
# con el patrón real de una curva que pasa R5 o con el máximo que R5 acepta (10 por mes en
# bloques de 4 + 4 + 2); se compara la curva rellenada con la del año completo. Error = RMSE
# entre curvas (K = 8) / variación entre años de las curvas de la estación (como en la sección 7).
R7_PHI <- base_bspline(1:365, 8)
falt_r5 <- diario_todo[anios[valido == TRUE, .(codigo, anio)], on = .(codigo, anio)][, .(o = list(!is.na(P))), by = .(codigo, anio)]
masc_r7 <- falt_r5[vapply(o, function(z) !all(z), TRUE), o]
masc_extrema <- function() { o <- rep(TRUE, 365)
  for (m in 1:12) { k <- which(MES_DIA == m); a <- sample(seq_len(length(k) - 12L), 1); o[k[c(a:(a + 3), (a + 5):(a + 8), (a + 10):(a + 11))]] <- FALSE }; o }
comp_r7 <- rbindlist(lapply(R7_EST, function(s) { if (!length(vecinas_r7[[s]])) return(NULL)
  a <- which(vapply(seq_len(sum(ANIOS <= 1995)), function(y) !anyNA(PM[(y - 1L) * 365L + 1:365, s]), TRUE))
  if (length(a) >= 3) data.table(s, y = a) }))
sd_r7 <- comp_r7[, .(sd = { X <- vapply(y, function(z) ajustar(PM[(z - 1L) * 365L + 1:365, s[1]], rep(TRUE, 365), R7_PHI), numeric(365))
                         sqrt(mean((X - rowMeans(X))^2)) }), by = s]
idw3 <- function(s, t) { v <- vecinas_r7[[s]]; v <- head(v[DKv[s, v] <= RADIO_VEC & !is.na(PM[t, v])], 3)
  if (!length(v)) return(NA_real_); m <- R7_MES[t]; w <- 1 / DKv[s, v]^2
  sum(w * PM[t, v] * R7_MEDIA[m, s] / pmax(R7_MEDIA[m, v], R7_MIN_MEDIA)) / sum(w) }
exp_r7 <- rbindlist(lapply(c("real", "maximo_R5"), function(patron) rbindlist(lapply(seq_len(nrow(comp_r7)), function(r) {
  s <- comp_r7$s[r]; idx <- (comp_r7$y[r] - 1L) * 365L + 1:365; y <- PM[idx, s]
  o <- if (patron == "real") masc_r7[[sample(length(masc_r7), 1)]] else masc_extrema()
  f <- which(!o); ref <- ajustar(y, rep(TRUE, 365), R7_PHI); sdv <- sd_r7$sd[sd_r7$s == s]
  err <- function(z) { oo <- o | !is.na(z); z[!oo] <- 0; sqrt(mean((ajustar(z, oo, R7_PHI) - ref)^2)) / sdv }
  rellena <- function(fun) { z <- y; z[f] <- vapply(idx[f], fun, 0); z }
  data.table(patron, faltan = length(f),
             solo_curva = err(replace(y, f, NA)),
             idw_3_vecinas = err(rellena(function(t) idw3(s, t))),
             kriging_sin_razon = err(rellena(function(t) estima_ok(s, t, razon = FALSE)[1])),
             kriging_con_razon = err(rellena(function(t) estima_ok(s, t)[1])))
}))))
res_r7 <- melt(exp_r7, id.vars = c("patron", "faltan"), variable.name = "metodo", value.name = "error")[
  , .(casos = .N, faltan_mediana = as.double(median(faltan)), error_mediana = round(median(error), 3), error_p90 = round(quantile(error, .9), 3)), by = .(patron, metodo)]
res_r7[, pct_casos_mejor_que_idw := NA_real_]   # en qué % de los casos el método deja menos error que IDW
for (pt in unique(exp_r7$patron)) for (mt in setdiff(levels(res_r7$metodo), "idw_3_vecinas"))
  res_r7[patron == pt & metodo == mt, pct_casos_mejor_que_idw := round(100 * exp_r7[patron == pt, mean(get(mt) < idw_3_vecinas)])]
print(tabla(res_r7, "8_R7_metodos_relleno"))
stopifnot(res_r7[patron == "maximo_R5" & metodo == "kriging_con_razon", error_mediana] < res_r7[patron == "maximo_R5" & metodo == "solo_curva", error_mediana])

# ---- R7: el flujo de kriging de la clase y la validación de cada decisión (D-64) ----
# El proyecto 1 del curso decide el kriging en pasos: modelo de la media → Moran de lo que queda
# → semivariograma (¿pepita pura?) → kriging solo si baja el error ≥1 % frente al modelo sin él.
# Aquí la "media" es la lluvia media de la estación en el mes (la razón de medias quita el nivel).
# (1) ¿Hay dependencia espacial? Moran de la anomalía diaria (P / media del mes), 4 vecinas más
#     cercanas y permutación, en días de entrenamiento al azar.
set.seed(20261122)
moran_perm <- function(q, D, B = 199L) { n <- length(q); nn <- t(apply(D, 1, function(r) order(r)[2:5])); q <- q - mean(q)
  f <- function(u) sum(u * rowMeans(matrix(u[nn], n, 4))) / sum(u^2); I <- f(q); c(I = I, p = (1 + sum(replicate(B, f(sample(q))) >= I)) / (B + 1)) }
R7_ZM <- PM / pmax(R7_MEDIA[cbind(R7_MES, rep(seq_along(R7_EST), each = nrow(PM)))], R7_MIN_MEDIA); dim(R7_ZM) <- dim(PM); dimnames(R7_ZM) <- dimnames(PM)
dias_moran <- sample(R7_ENT[rowSums(!is.na(PM[R7_ENT, ])) >= 60], 300)
moran_r7 <- rbindlist(lapply(dias_moran, function(t) { ok <- which(!is.na(R7_ZM[t, ])); z <- R7_ZM[t, ok]
  m <- if (sd(z) > 0) moran_perm(z, DKv[ok, ok]) else c(I = NA, p = 1)
  data.table(t, estaciones = length(ok), I = m[["I"]], p = m[["p"]], pct_con_lluvia = mean(PM[t, ok] > 0)) }))
moran_r7[, tipo_de_dia := cut(pct_con_lluvia, c(-.01, .2, .6, 1), labels = c("seco (<20 % de estaciones con lluvia)", "intermedio", "húmedo (>60 %)"))]
moran_res <- rbind(moran_r7[, .(tipo_de_dia = "todos", dias = .N, I_mediano = round(median(I, na.rm = TRUE), 2), pct_dias_p_menor_005 = round(100 * mean(p < .05)))],
                   moran_r7[, .(dias = .N, I_mediano = round(median(I, na.rm = TRUE), 2), pct_dias_p_menor_005 = round(100 * mean(p < .05))), by = tipo_de_dia][order(tipo_de_dia)])
print(tabla(moran_res, "8_R7_moran"))
# (2) ¿Es pepita pura? En la clase, meseta parcial / meseta total < 0,05 → no hay kriging
cat(sprintf("Semivariograma elegido (%s): meseta parcial / meseta total = %.2f (pepita pura si < 0,05)\n",
            vg_mejor, vg_par[[vg_mejor]][2] / sum(vg_par[[vg_mejor]][1:2])))
# (3) Validación de cada decisión, sin evaluar con los datos de ajuste: las medias mensuales y los
#     semivariogramas se ajustan con los años IMPARES de 1981–1995 y se quitan días (el máximo que
#     acepta R5) a años completos PARES. Un cambio se adopta solo si su IC 95 % (bootstrap de
#     estaciones) excluye el 0 y baja el error al menos un 1 % (la regla del proyecto 1, MEJORA_MIN).
R7_ANIO <- rep(ANIOS, each = 365L); R7_AJ <- R7_ENT[R7_ANIO[R7_ENT] %% 2L == 1L]
R7_FASE <- fread("datos_web/oni_anual_fase.csv")[match(R7_ANIO, anio), fase]
R7_ESTAC <- c("DEF", "DEF", "MAM", "MAM", "MAM", "JJA", "JJA", "JJA", "SON", "SON", "SON", "DEF")[R7_MES]
R7_ZONA <- setNames(ifelse(valle$vertiente %in% "Pacífico", "Pacífico", "Andina"), valle$codigo)[R7_EST]
med_imp <- sapply(R7_EST, function(s) tapply(PM[R7_AJ, s], R7_MES[R7_AJ], function(v) if (sum(!is.na(v)) >= 20) mean(v, na.rm = TRUE) else NA_real_))
med_imp <- ifelse(is.na(med_imp), R7_MEDIA, med_imp)
emp_vg <- function(M, filas, zona = NULL) { acc <- matrix(0, length(tramos) - 1L, 2L)
  for (i in seq_along(R7_EST)) { if (!is.null(zona) && R7_ZONA[i] != zona) next
    j <- which(DKv[i, ] > 0 & DKv[i, ] <= 100 & (is.null(zona) | R7_ZONA == if (is.null(zona)) "" else zona)); if (!length(j)) next
    q <- (M[filas, j, drop = FALSE] - M[filas, i])^2; n_c <- colSums(!is.na(q)); m_c <- colSums(q, na.rm = TRUE) / pmax(n_c, 1)
    b <- findInterval(DKv[i, j], tramos, left.open = TRUE); ok <- n_c >= 60 & b >= 1 & b < length(tramos)
    for (k in which(ok)) acc[b[k], ] <- acc[b[k], ] + c(m_c[k] / 2, 1) }
  data.table(km = (tramos[-1] + tramos[-length(tramos)]) / 2, gamma = acc[, 1] / acc[, 2], pares = acc[, 2])[pares > 0] }
ajusta_exp <- function(e) { ob <- function(w) { p <- c(exp(w[1]), exp(w[2]), 1 + 299 * plogis(w[3])); sum(e$pares * (e$gamma - VG$exponencial(e$km, p))^2) }
  best <- NULL; for (c0 in c(.3, .6, .9) * min(e$gamma)) for (a in c(10, 40, 150)) {
    o <- optim(c(log(c0), log(max(e$gamma) - c0 + 1e-6), qlogis((a - 1) / 299)), ob, control = list(maxit = 3000)); if (is.null(best) || o$value < best$value) best <- o }
  c(exp(best$par[1:2]), 1 + 299 * plogis(best$par[3])) }
VV <- list(global = ajusta_exp(emp_vg(PM, R7_AJ)))
for (e in unique(R7_ESTAC)) VV[[paste0("temporada_", e)]] <- ajusta_exp(emp_vg(PM, R7_AJ[R7_ESTAC[R7_AJ] == e]))
for (f in c("Niño", "Niña", "Neutro")) VV[[paste0("ENSO_", f)]] <- ajusta_exp(emp_vg(PM, R7_AJ[R7_FASE[R7_AJ] == f]))
for (z in c("Andina", "Pacífico")) VV[[paste0("zona_", z)]] <- ajusta_exp(emp_vg(PM, R7_AJ, z))
# kriging simple sobre la anomalía (media conocida = 1): necesita la covarianza de la anomalía
R7_ZI <- PM / pmax(med_imp[cbind(R7_MES, rep(seq_along(R7_EST), each = nrow(PM)))], R7_MIN_MEDIA); dim(R7_ZI) <- dim(PM); dimnames(R7_ZI) <- dimnames(PM)
VV$anomalia <- ajusta_exp(emp_vg(R7_ZI, R7_AJ))
tabla(rbindlist(lapply(names(VV), function(n) data.table(semivariograma = n, pepita = VV[[n]][1], meseta_parcial = VV[[n]][2], rango_km = VV[[n]][3]))), "8_R7_semivariogramas_por_grupo")
est_val <- function(metodo, s, t, media_anual) {
  m <- R7_MES[t]; tend <- med_imp[m, s]
  if (metodo == "tendencia (media del mes)") return(tend)
  if (metodo == "tesis del profesor (media anual de la estación)") return(media_anual)
  p <- VV$global; k <- R7_MAX_VEC; corr <- FALSE
  if (metodo == "semivariograma por temporada") p <- VV[[paste0("temporada_", R7_ESTAC[t])]]
  if (metodo == "semivariograma por fase ENSO") p <- VV[[paste0("ENSO_", R7_FASE[t])]]
  if (metodo == "semivariograma por zona") p <- VV[[paste0("zona_", R7_ZONA[[s]])]]
  if (metodo == "5 vecinas") k <- 5L
  if (metodo == "20 vecinas") k <- 20L
  if (metodo == "pesos negativos a 0 (Deutsch 1996)") corr <- TRUE
  v <- vecinas_r7[[s]]; v <- head(v[!is.na(PM[t, v])], k); if (!length(v)) return(tend)
  zz <- PM[t, v] / pmax(med_imp[m, v], R7_MIN_MEDIA)
  if (metodo == "kriging simple") { pa <- VV$anomalia; C <- function(h) ifelse(h == 0, pa[1] + pa[2], pa[2] * exp(-h / pa[3]))
    w <- tryCatch(solve(C(DKv[v, v, drop = FALSE]), C(DKv[s, v])), error = function(e) rep(0, length(v))); return(max(0, tend * (1 + sum(w * (zz - 1))))) }
  w <- pesos_ok(s, v, VG$exponencial, p)$w; if (corr) { w <- pmax(w, 0); w <- w / sum(w) }
  max(0, tend * sum(w * zz)) }
METODOS_R7 <- c("R7 (un semivariograma, 10 vecinas)", "tendencia (media del mes)", "tesis del profesor (media anual de la estación)",
                "semivariograma por temporada", "semivariograma por fase ENSO", "semivariograma por zona", "5 vecinas", "20 vecinas",
                "pesos negativos a 0 (Deutsch 1996)", "kriging simple")
casos_val <- comp_r7[y %% 2L == 0L][sample(.N, min(.N, 400))]
val_r7 <- rbindlist(lapply(seq_len(nrow(casos_val)), function(r) {
  s <- casos_val$s[r]; idx <- (casos_val$y[r] - 1L) * 365L + 1:365; obs <- PM[idx, s]; o <- masc_extrema(); f <- which(!o)
  ref <- ajustar(obs, rep(TRUE, 365), R7_PHI); sdv <- sd_r7$sd[sd_r7$s == s]; ma <- mean(obs[o])
  rbindlist(lapply(METODOS_R7, function(mt) { z <- obs; z[f] <- vapply(idx[f], function(t) est_val(mt, s, t, ma), 0)
    data.table(s, anio = ANIOS[casos_val$y[r]], fase = R7_FASE[idx[1]], zona = R7_ZONA[[s]], metodo = mt,
               curva = sqrt(mean((ajustar(z, rep(TRUE, 365), R7_PHI) - ref)^2)) / sdv, dia = sqrt(mean((z[f] - obs[f])^2))) })) }))
val_r7 <- val_r7[val_r7[metodo == METODOS_R7[1], .(s, anio, c0 = curva)], on = .(s, anio)][, dif := curva - c0]
ic_est <- function(dd) { u <- unique(dd$s); q <- replicate(1000, { k <- sample(u, replace = TRUE); dd[data.table(s = k), on = "s", allow.cartesian = TRUE][, mean(dif)] }); quantile(q, c(.025, .975)) }
dec_r7 <- val_r7[, { ic <- ic_est(.SD); .(casos = .N, error_curva_media = round(mean(curva), 4), rmse_dia_mm = round(mean(dia), 2), dif_vs_R7 = round(mean(dif), 4),
                 ic_inf = round(ic[1], 4), ic_sup = round(ic[2], 4), cambio_pct = round(100 * mean(dif) / mean(c0), 1)) }, by = metodo]
dec_r7[, decision := fifelse(metodo == METODOS_R7[1], "elegido",
                     fifelse(ic_sup < 0 & -cambio_pct >= 100 * 0.01, "mejora ≥1 %: habría que adoptarlo",
                     fifelse(ic_inf > 0, "peor que R7", "no mejora ≥1 % con IC que excluya el 0")))]
print(tabla(dec_r7, "8_R7_validacion_decisiones"))
print(tabla(val_r7[metodo %in% METODOS_R7[c(1, 2, 5)], .(casos = .N, error_curva_media = round(mean(curva), 4)), by = .(fase, metodo)][order(fase, metodo)], "8_R7_por_fase_ENSO"))
stopifnot(dec_r7[decision == "mejora ≥1 %: habría que adoptarlo", .N] == 0L)
# (4) ¿Un semivariograma por día, como el proyecto 1 por corte? Dejar una estación fuera en días de
#     entrenamiento: si ese día Moran no es significativo, o su semivariograma es pepita pura, se usa la
#     tendencia; si no, kriging con el semivariograma de ESE día (exp, esf o gauss por MCP, menor SSE).
fam3 <- VG[c("exponencial", "esferico", "gaussiano")]
ajusta_dia <- function(z, ok) { D <- DKv[ok, ok]; ij <- which(upper.tri(D) & D <= 100, arr.ind = TRUE); h <- D[ij]; g <- (z[ok][ij[, 1]] - z[ok][ij[, 2]])^2 / 2
  e <- data.table(h, g, b = cut(h, tramos))[, .(km = mean(h), gamma = mean(g), pares = .N), by = b][pares >= 10]
  if (nrow(e) < 4) return(NULL); best <- NULL
  for (nm in names(fam3)) { ob <- function(w) { p <- c(exp(w[1]), exp(w[2]), 1 + 299 * plogis(w[3])); sum(e$pares * (e$gamma - fam3[[nm]](e$km, p))^2) }
    o <- optim(c(log(max(min(e$gamma), 1e-3)), log(max(diff(range(e$gamma)), 1e-3)), 0), ob); if (is.null(best) || o$value < best$sse) best <- list(nm = nm, p = c(exp(o$par[1:2]), 1 + 299 * plogis(o$par[3])), sse = o$value) }
  best }
dia_r7 <- rbindlist(lapply(dias_moran, function(t) { m <- R7_MES[t]; ok <- which(!is.na(R7_ZM[t, ])); z <- R7_ZM[t, ]
  mo <- moran_r7$p[moran_r7$t == t]; fit <- if (mo < 0.05) ajusta_dia(z, ok) else NULL
  usa <- !is.null(fit) && fit$p[2] / sum(fit$p[1:2]) >= 0.05
  val_t <- anios[valido == TRUE & anio == R7_ANIO[t], codigo]; obj <- ok[R7_EST[ok] %in% val_t]
  rbindlist(lapply(obj[sample.int(length(obj), min(30L, length(obj)))], function(si) { s <- R7_EST[si]; v <- setdiff(R7_EST[ok][order(DKv[s, R7_EST[ok]])], s); v <- v[v %in% R7_CON_VALIDOS]; v <- head(v[DKv[s, v] <= R7_RADIO], R7_MAX_VEC)
    if (!length(v)) return(NULL); tend <- R7_MEDIA[m, s]
    kg <- max(0, tend * sum(pesos_ok(s, v, VG[[vg_mejor]], vg_par[[vg_mejor]])$w * z[v]))
    kd <- if (usa) max(0, tend * sum(pesos_ok(s, v, fam3[[fit$nm]], fit$p)$w * z[v])) else tend
    data.table(t, obs = PM[t, s], tendencia = tend, R7 = kg, por_dia = kd, usa_kriging_ese_dia = usa) })) }))
rm2 <- function(a, b) sqrt(mean((a - b)^2))
ud <- unique(dia_r7$t); bs_d <- replicate(1000, { k <- sample(ud, replace = TRUE); w <- dia_r7[data.table(t = k), on = "t", allow.cartesian = TRUE]; rm2(w$por_dia, w$obs) - rm2(w$R7, w$obs) })
por_dia_r7 <- data.table(validacion = "dejar una estación fuera", dias = length(ud), predicciones = nrow(dia_r7), pct_dias_con_kriging_por_dia = round(100 * dia_r7[!duplicated(t), mean(usa_kriging_ese_dia)]),
  rmse_tendencia = round(rm2(dia_r7$tendencia, dia_r7$obs), 2), rmse_R7 = round(rm2(dia_r7$R7, dia_r7$obs), 2), rmse_por_dia = round(rm2(dia_r7$por_dia, dia_r7$obs), 2),
  dif_por_dia_menos_R7 = round(rm2(dia_r7$por_dia, dia_r7$obs) - rm2(dia_r7$R7, dia_r7$obs), 3), ic_inf = round(quantile(bs_d, .025), 3), ic_sup = round(quantile(bs_d, .975), 3))
print(tabla(por_dia_r7, "8_R7_un_semivariograma_vs_por_dia"))
rm(R7_ZM, R7_ZI); invisible(gc())

# Aplicar: se rellenan los faltantes de los años válidos de las estaciones del Valle
diario_todo[, imputado := FALSE]
obj_r7 <- diario_todo[is.na(P)][anios[valido == TRUE, .(codigo, anio)], on = .(codigo, anio), nomatch = 0][, .(codigo, anio, dia)]
obj_r7[, t := (anio - ANIOS[1]) * 365L + dia]
est_r7 <- t(mapply(estima_ok, obj_r7$codigo, obj_r7$t))
obj_r7[, `:=`(P = est_r7[, 1], sd_kriging = est_r7[, 2], vecinas = as.integer(est_r7[, 3]))]
diario_todo[obj_r7[!is.na(P)], on = .(codigo, anio, dia), `:=`(P = i.P, imputado = TRUE)]
imputados <- valle[, .(codigo, nombre)][obj_r7, on = "codigo"]
tabla(imputados, "8_R7_dias_rellenados")
cat(sprintf("R7: %d días faltantes en años válidos; %d rellenados por kriging (mediana %d vecinas), %d sin vecina con dato (quedan para la curva)\n",
            nrow(obj_r7), obj_r7[!is.na(P), .N], as.integer(obj_r7[!is.na(P), median(vecinas)]), obj_r7[is.na(P), .N]))
rm(PM, Z, est_r7); invisible(gc())

# R6 (D-63). Estación incluida: ≥10 años válidos en 1981–2003 (mínimo de la OMM para cualquier
# promedio climatológico; normales "representativas" de NOAA). Un año en que la estación no
# operaba cuenta como no válido; no se pide que los años sean seguidos (WMO 2007: no reduce la
# incertidumbre). El conjunto de sensibilidad (CONJUNTO=r6_omm) pide el 80 % de la OMM para
# normales estándar: 19 de 23 (ver la configuración).
por_est <- anios[, .(validos = sum(valido), validos_ent = sum(valido & particion == "entrenamiento"),
                     validos_val = sum(valido & particion == "validación"),
                     validos_pru = sum(valido & particion == "prueba")), by = codigo]
estaciones <- merge(valle, por_est, by = "codigo")[validos >= ANIOS_MIN]
r6_umbrales <- rbindlist(lapply(c(19L, 15L, 10L), function(u) {
  e <- por_est[validos >= u, codigo]
  data.table(anios_validos_min = u, estaciones = length(e), pacificas = valle[codigo %in% e, sum(vertiente == "Pacífico")],
             curvas = anios[valido == TRUE & codigo %in% e, .N], usado_en_este_conjunto = u == ANIOS_MIN) }))
print(tabla(r6_umbrales, "8_R6_umbrales"))
diario <- diario_todo[anios[valido == TRUE, .(codigo, anio)], on = .(codigo, anio), nomatch = 0][codigo %in% estaciones$codigo]

# Trazabilidad (rúbrica 4.1): estaciones y pares estación-año antes y después de cada paso
cuenta <- function(x) x[, .(estaciones = uniqueN(codigo), pares_estacion_anio = uniqueN(paste(codigo, anio)), dias_con_dato = .N)]
d_r2 <- d[!m_ext, on = .(codigo, fecha)]; d_r3 <- d_r2[!m_mes, on = .(codigo, fecha)]; d_r4 <- d_r3[!m_cero, on = .(codigo, fecha)]
trazabilidad <- rbind(
  data.table(paso = "0. Crudo: códigos con datos en 1981–2003", estaciones = nrow(cobertura_cruda),
             pares_estacion_anio = sum(cobertura_cruda$anios_con_dato), dias_con_dato = sum(cobertura_cruda$dias_con_dato)),
  data.table(paso = "1. Solo el Valle, con coordenada (D-34)", estaciones = nrow(valle),
             pares_estacion_anio = sum(valle$anios_con_dato), dias_con_dato = sum(valle$dias_con_dato)),
  cbind(paso = "1b. R0: negativos y valores sobre el récord OMM (1 825 mm) a faltante", traza_r0),
  cbind(paso = "2. R1: sin 29 de febrero", cuenta(d)),
  cbind(paso = "3. R2: extremos sin respaldo a faltante", cuenta(d_r2)),
  cbind(paso = "4. R3: totales mensuales en un día a faltante", cuenta(d_r3)),
  cbind(paso = "5. R4: rachas de ceros que ≥2 de 3 vecinas contradicen a faltante", cuenta(d_r4)),
  cbind(paso = "6. R5 (OMM): años con todos sus meses con ≤10 faltantes y <5 seguidos", cuenta(diario_todo[!is.na(P) & !imputado][anios[valido == TRUE, .(codigo, anio)], on = .(codigo, anio), nomatch = 0])),
  cbind(paso = "6b. R7: días faltantes de esos años rellenados por kriging (suman)", cuenta(diario_todo[!is.na(P)][anios[valido == TRUE, .(codigo, anio)], on = .(codigo, anio), nomatch = 0])),
  cbind(paso = paste("7. R6: estaciones con", R6_TEXTO), cuenta(diario[!is.na(P)])))
trazabilidad[, `:=`(estaciones_que_salen = c(NA, -diff(estaciones)), pares_que_salen = c(NA, -diff(pares_estacion_anio)),
                    dias_que_salen = c(NA, -diff(dias_con_dato)))]
print(tabla(trazabilidad, "8_trazabilidad"))
rm(d_r2, d_r3, d_r4)

reglas <- c(
  R0 = "Valor imposible → faltante: P < 0 o P > 1 825 mm (récord mundial OMM en 24 h)",
  R1 = "Dominio 1–365: se quita el 29 de febrero",
  R2 = "Extremo > 3·p99 y ≥150 mm sin vecina (≤25 km, ±1 día) que supere P/3 → faltante",
  R3 = "Mes con un solo día de lluvia > p99 y (mes con <50 % de días o vecinas con lluvia ≥40 % de los días) → mes faltante",
  R4 = "Racha ≥30 días en 0 → faltante si ≥2 de sus 3 vecinas más cercanas (≤25 km) la contradicen: (1 − q)^k < 0,05/rachas, q de la temporada en entrenamiento",
  R5 = "Año válido (OMM 2017, 4.4.1): en cada mes faltan ≤10 días y no hay 5 o más seguidos",
  R6 = paste("Estación: años válidos", R6_TEXTO),
  R7 = sprintf("Faltantes de los años válidos rellenados por kriging ordinario (%s) con ≤10 vecinas a ≤50 km, escaladas por la razón de medias mensuales; sin satélite", vg_mejor))

saveRDS(list(diario = diario, diario_todo = diario_todo, anios = anios, estaciones = estaciones,
             marcas = marcas, imputados = imputados, meses_sospechosos = meses_sosp, reglas = reglas, conjunto = CONJUNTO),
        file.path(SALIDA, "lluvia_valle_limpia.rds"), compress = "xz")
tabla(marcas, "8_marcas_limpieza"); tabla(anios, "8_anios_validez"); tabla(estaciones, "8_estaciones_incluidas")

# Representatividad: ¿las incluidas se parecen a todas las del Valle?
repres <- rbind(valle[, .(grupo = "220 con datos", estaciones = .N, altitud_mediana = median(altitud, na.rm = TRUE),
                          pct_pacifico = round(100 * mean(vertiente == "Pacífico")))],
                estaciones[, .(grupo = "incluidas", estaciones = .N, altitud_mediana = median(altitud, na.rm = TRUE),
                               pct_pacifico = round(100 * mean(vertiente == "Pacífico")))])
repres$P_anual_mediana <- c(median(clim$P_media), median(clim[codigo %in% estaciones$codigo, P_media]))
print(tabla(repres, "8_representatividad"))
grafico(ggplot(valle, aes(lon, lat, colour = codigo %in% estaciones$codigo)) + geom_point(size = 2) + coord_equal() +
          scale_colour_manual(values = c(`TRUE` = "steelblue", `FALSE` = "red"), labels = c(`TRUE` = "incluida", `FALSE` = "excluida (pocos años válidos)"), name = NULL) +
          labs(title = "Estaciones incluidas en el dataset final", x = NULL, y = NULL),
        "8_mapa_incluidas", ancho = 7, alto = 7)

# Qué curvas entran al análisis: las estaciones incluidas × 1981–2003, con el nombre de cada
# una. Cada casilla es una curva estación-año: el color dice si es válida (R5) y a qué
# conjunto va. Se parte en dos imágenes para que los nombres se puedan leer.
ESTADO_CURVA <- c("válida: entrenamiento", "válida: validación", "válida: prueba",
                  "no válida: algún mes con ≥11 faltantes o ≥5 seguidos (OMM)", "sin datos ese año")
ver <- anios[, .(codigo, anio, n_obs, valido)][CJ(codigo = estaciones$codigo, anio = ANIOS), on = .(codigo, anio)]
ver[, estado := factor(fcase(valido %in% TRUE, paste("válida:", PARTICION(anio)), !is.na(n_obs) & n_obs > 0, ESTADO_CURVA[4],
                             default = ESTADO_CURVA[5]), levels = ESTADO_CURVA)]
ver <- estaciones[, .(codigo, nombre, lat, zona = ETQ_VERT[vertiente])][ver, on = "codigo"]
ver[, etq := sprintf("%s · %s", tools::toTitleCase(tolower(nombre)), codigo)]
ver[, zona := factor(zona, levels = rev(ETQ_VERT))]                       # zona pacífica arriba
orden_v <- unique(ver[order(zona == ETQ_VERT[1], -lat), etq])          # zona pacífica primero, de norte a sur
n_val <- ver[, .(n = sum(valido %in% TRUE)), by = .(etq, zona)]
stopifnot(sum(n_val$n) == anios[codigo %in% estaciones$codigo, sum(valido)], !anyDuplicated(orden_v))
print(ver[, .N, by = estado][order(estado)])
partes <- split(orden_v, ceiling(seq_along(orden_v) / ceiling(length(orden_v) / 2)))
for (k in seq_along(partes)) {
  vk <- ver[etq %in% partes[[k]]]; vk[, etq := factor(etq, levels = rev(partes[[k]]))]
  grafico(ggplot(vk, aes(anio, etq, fill = estado)) + geom_tile(colour = "white", linewidth = .2) +
            geom_text(data = n_val[etq %in% partes[[k]]], aes(x = 2004.4, y = etq, label = n), inherit.aes = FALSE, size = 2.3) +
            facet_grid(zona ~ ., scales = "free_y", space = "free_y", labeller = label_wrap_gen(14)) +
            scale_fill_manual(values = setNames(c("#2166AC", "#F4A582", "#B2182B", "grey70", "grey95"), ESTADO_CURVA),
                              name = NULL, drop = FALSE) +
            geom_vline(xintercept = c(1995.5, 2000.5), colour = "black", linetype = 2) + EJE_PARTICION +
            coord_cartesian(xlim = c(1980.5, 2005), clip = "off") +
            labs(x = NULL, y = NULL,
                 title = sprintf("Curvas que entran al análisis (parte %d de %d)", k, length(partes)),
                 subtitle = "Cada fila es una estación incluida (nombre · código), de norte a sur.\nNúmero a la derecha: curvas válidas que aporta la estación.") +
            theme(axis.text.y = element_text(size = 6.2), panel.grid = element_blank(), legend.position = "bottom",
                  plot.title = element_text(size = 11), plot.subtitle = element_text(size = 8.5), legend.text = element_text(size = 8),
                  axis.text.x.top = element_text(size = 6.5),
                  strip.text.y = element_text(angle = 0, face = "bold")) +
            guides(fill = guide_legend(nrow = 3)),
          sprintf("3_curvas_validas_parte%d", k), ancho = 6.2, alto = 8.8, dpi = 200)
}


# ---- 9. SATÉLITE: LECTURA, ESTRUCTURA Y FUENTE --------------------------------
titulo("9. Satélite: lectura, estructura y fuente")
# 2 328 GeoTIFF (ppt + año + péntada 01–72), uno por péntada de 1981 a abril de
# 2013. Cada mes tiene 6 péntadas: días 1–5, 6–10, …, 21–25 y 26–fin de mes.
if (!requireNamespace("terra", quietly = TRUE)) stop("Falta el paquete terra (lectura de los GeoTIFF)")
arch_sat <- sort(list.files("datos_curso/ImagenesSatelitePentadal", "^ppt[0-9]{6}\\.tif$", full.names = TRUE))
inv_sat <- data.table(archivo = arch_sat, anio = as.integer(substr(basename(arch_sat), 4, 7)),
                      pentada = as.integer(substr(basename(arch_sat), 8, 9)))
print(inv_sat[, .N, by = anio][, .(anios = .N, con_72_pentadas = sum(N == 72),
                                   incompletos = paste(anio[N != 72], N[N != 72], sep = ": ", collapse = "; "))])
# Se lee por lotes de 200: abrir los 2 328 a la vez supera el límite de archivos
# abiertos del sistema (1 024 en macOS) y terra se detiene. gc() cierra cada lote.
lee_lote <- function(a) { v <- terra::values(terra::rast(a), mat = TRUE); invisible(gc()); v }
S_todo <- do.call(cbind, lapply(split(inv_sat$archivo, ceiling(seq_len(nrow(inv_sat)) / 200)), lee_lote))  # píxel × péntada
r_sat <- terra::rast(inv_sat$archivo[1])             # solo la geometría (todas comparten la rejilla)
stopifnot(ncol(S_todo) == nrow(inv_sat), identical(colnames(S_todo), sub("\\.tif$", "", basename(inv_sat$archivo))))
print(r_sat)
storage.mode(S_todo) <- "double"
estructura_sat <- data.table(
  archivos = ncol(S_todo), pixeles = nrow(S_todo), filas = terra::nrow(r_sat), columnas = terra::ncol(r_sat),
  resolucion_grados = terra::res(r_sat)[1], sin_dato = sum(is.na(S_todo)), minimo = min(S_todo), maximo = max(S_todo),
  no_enteros = sum(S_todo != round(S_todo)), ceros = sum(S_todo == 0),
  en_el_tope_200 = sum(S_todo == 200), pct_tope = round(100 * mean(S_todo == 200), 2))
print(tabla(estructura_sat, "9_estructura_satelite"))

grafico(ggplot(data.table(v = as.vector(S_todo))[, .N, by = v], aes(v, N)) + geom_col(width = 1) + scale_y_log10() +
          labs(x = "mm por péntada", y = "Valores (log)", title = "Satélite: valores enteros de 0 a 200, con tope en 200",
               subtitle = "La barra de 200 acumula todo lo que superaba el tope"),
        "9_satelite_valores")

# Píxeles: centro, y si caen dentro del polígono DANE del Valle (D-03)
pix <- data.table(pixel = seq_len(nrow(S_todo)), terra::xyFromCell(r_sat, seq_len(nrow(S_todo))))
setnames(pix, c("x", "y"), c("lon", "lat"))
poly_valle <- deptos[deptos$shapeName == "Valle del Cauca", ]
pix[, en_valle := lengths(suppressMessages(st_intersects(st_as_sf(pix, coords = c("lon", "lat"), crs = 4326), poly_valle))) > 0]
# También se conservan los píxeles que contienen una estación incluida aunque su
# centro quede fuera del polígono (estaciones cerca del borde del departamento)
pix[, con_estacion := pixel %in% terra::cellFromXY(r_sat, as.matrix(estaciones[, .(lon, lat)]))]
pix[, usar := en_valle | con_estacion]
cat("Píxeles de la ventana:", nrow(pix), "| con centro dentro del Valle:", pix[en_valle == TRUE, .N],
    "| con estación y centro fuera:", pix[con_estacion & !en_valle, .N], "| se usan:", pix[usar == TRUE, .N], "\n")

# Fuente (D-43): se compara con CHIRP v2.0 (solo infrarrojo) y con CHIRPS v2.0
# (con estaciones), leídos del servidor del CHC para la misma ventana.
leer_ref <- function(f_bin, f_claves) {
  cl <- readLines(f_claves, warn = FALSE); cl <- cl[nzchar(cl)]
  M <- matrix(readBin(f_bin, "numeric", n = length(cl) * nrow(S_todo), size = 4, endian = "little"), nrow = nrow(S_todo))
  M[M < -100] <- NA                                  # sin dato (mar en CHIRPS)
  p <- do.call(rbind, strsplit(cl, ".", fixed = TRUE))
  colnames(M) <- sprintf("%s%02d", p[, 1], (as.integer(p[, 2]) - 1L) * 6L + as.integer(p[, 3]))
  M
}
colnames(S_todo) <- sprintf("%d%02d", inv_sat$anio, inv_sat$pentada)
chirp  <- leer_ref("datos_web/chirp_v2_valle_muestra_7anios.bin", "datos_web/chirp_v2_valle_muestra_7anios_claves.txt")
chirps <- leer_ref("datos_web/chirps_v2_valle_1981_2013.bin", "datos_web/chirps_v2_valle_1981_2013_claves.txt")
kk <- Reduce(intersect, list(colnames(S_todo), colnames(chirp), colnames(chirps)))
fuente <- rbindlist(lapply(split(kk, substr(kk, 1, 4)), function(k) {
  e <- S_todo[, k]; a <- chirp[, k]; b <- chirps[, k]; ok <- !is.na(b)
  data.table(anio = substr(k[1], 1, 4), pentadas = length(k),
             r_con_CHIRP = round(cor(as.vector(e), as.vector(a)), 3), r_con_CHIRPS = round(cor(e[ok], b[ok]), 3),
             mm_entregado = round(sum(e) / nrow(e)), mm_CHIRP = round(sum(a) / nrow(a)),
             mm_entregado_tierra = round(sum(e[ok]) / sum(ok[, 1])), mm_CHIRPS_tierra = round(sum(b[ok]) / sum(ok[, 1])))
}))
print(tabla(fuente, "9_fuente_satelite"))
cat("r del mapa medio: con CHIRP", round(cor(rowMeans(S_todo[, kk]), rowMeans(chirp[, kk])), 3),
    "| con CHIRPS", round(cor(rowMeans(S_todo[, kk]), rowMeans(chirps[, kk]), use = "complete.obs"), 3),
    "| píxeles sin dato en CHIRPS (mar):", sum(is.na(chirps[, 1])), "| tope de CHIRP:", round(max(chirp)), "\n")

# QUÉ DICEN LOS DATOS (D-43)
#   Las cifras y su interpretación están en EXPLICACION_DATOS.md, sección 9 (corrida del
#   19-nov-2026, con las 98 estaciones de la regla de la OMM); las tablas de esta
#   sección las imprimen las líneas de arriba.
rm(chirp, chirps, r_sat); invisible(gc())


# ---- 10. SATÉLITE: RECORTE AL VALLE, CALENDARIO Y UNIDADES ---------------------
titulo("10. Satélite: recorte, calendario y unidades")
# Calendario de péntadas en el dominio 1–365 (sin 29-feb). La 6.ª péntada de cada
# mes va del 26 al fin de mes: 3 días en febrero (4 en bisiesto), 5 o 6 en el resto.
dias_mes_nb <- c(31L, 28L, 31L, 30L, 31L, 30L, 31L, 31L, 30L, 31L, 30L, 31L)
cal <- CJ(anio = ANIOS, pentada = 1:72)
cal[, `:=`(mes = (pentada - 1L) %/% 6L + 1L, p = (pentada - 1L) %% 6L + 1L)]
cal[, dias_dominio := fifelse(p < 6L, 5L, dias_mes_nb[mes] - 25L)]
cal[, dias_reales := dias_dominio + (p == 6L & mes == 2L & anio %% 4L == 0L)]
cal[, dia_ini := c(0L, cumsum(dias_mes_nb))[mes] + (p - 1L) * 5L + 1L]
cal[, dia_centro := dia_ini + (dias_dominio - 1) / 2]
stopifnot(cal[, sum(dias_dominio), by = anio][, all(V1 == 365)])
print(cal[anio == 1981, .N, by = .(dias = dias_dominio)][order(dias)])

# Recorte: periodo 1981–2003 y píxeles dentro del Valle. Unidades: mm/día (D-15),
# dividiendo cada péntada por sus días reales.
cols <- sprintf("%d%02d", cal$anio, cal$pentada)
stopifnot(all(cols %in% colnames(S_todo)))
S <- S_todo[pix[usar == TRUE, pixel], cols]
sat <- data.table(pixel = rep(pix[usar == TRUE, pixel], times = ncol(S)),
                  anio = rep(cal$anio, each = nrow(S)), pentada = rep(cal$pentada, each = nrow(S)),
                  mm_pentada = as.vector(S))
sat <- cal[, .(anio, pentada, mes, dia_centro, dias_reales)][sat, on = .(anio, pentada)]
sat[, mm_dia := mm_pentada / dias_reales]
cat("Satélite recortado:", uniqueN(sat$pixel), "píxeles ×", nrow(cal), "péntadas =", nrow(sat), "valores\n")
# Sin dividir, la 6.ª péntada de los meses de 31 días parece más lluviosa
print(sat[, .(mm_pentada = round(mean(mm_pentada), 1), mm_dia = round(mean(mm_dia), 2)), by = .(dias = dias_reales)][order(dias)])


# ---- 11. SATÉLITE: EDA ----------------------------------------------------------
titulo("11. Satélite: EDA")
sat_anual <- sat[, .(mm = sum(mm_pentada), pct_tope = 100 * mean(mm_pentada == 200)), by = .(pixel, anio)]
sat_pix <- sat_anual[, .(mm_anual = mean(mm), cv = sd(mm) / mean(mm), pct_tope = mean(pct_tope)), by = pixel][pix, on = "pixel", nomatch = 0]
tabla(sat_pix, "11_satelite_pixel")
cat("Lluvia anual del satélite por píxel: de", round(min(sat_pix$mm_anual)), "a", round(max(sat_pix$mm_anual)),
    "mm | estaciones incluidas: de", round(min(clim[codigo %in% estaciones$codigo, P_media])), "a",
    round(max(clim[codigo %in% estaciones$codigo, P_media])), "mm\n")

grafico(ggplot(sat_pix, aes(lon, lat, fill = mm_anual)) + geom_tile() +
          geom_sf(data = poly_valle, inherit.aes = FALSE, fill = NA, colour = "black") +
          geom_point(data = estaciones, aes(lon, lat), inherit.aes = FALSE, size = .6, colour = "white") +
          scale_fill_viridis_c(name = "mm/año") + labs(title = "Satélite: lluvia media anual por píxel (1981–2003)",
                                                       subtitle = "Puntos blancos: estaciones incluidas", x = NULL, y = NULL),
        "11_satelite_mapa_anual", ancho = 7, alto = 7)
grafico(ggplot(sat_pix, aes(lon, lat, fill = pct_tope)) + geom_tile() +
          geom_sf(data = poly_valle, inherit.aes = FALSE, fill = NA, colour = "black") +
          scale_fill_viridis_c(name = "% de péntadas\nen 200 mm", option = "magma") +
          labs(title = "¿Dónde pega el tope de 200 mm?", x = NULL, y = NULL),
        "11_satelite_mapa_tope", ancho = 7, alto = 7)

# Ciclo anual (mm/día) por péntada, con la dispersión entre píxeles
ciclo_sat <- sat[, .(mm_dia = mean(mm_dia)), by = .(pixel, pentada, dia_centro)]
grafico(ggplot(ciclo_sat[, .(p10 = quantile(mm_dia, .1), p50 = median(mm_dia), p90 = quantile(mm_dia, .9)), by = dia_centro],
               aes(dia_centro)) + geom_ribbon(aes(ymin = p10, ymax = p90), alpha = .25) + geom_line(aes(y = p50)) +
          labs(x = "Día del año", y = "mm/día", title = "Satélite: ciclo anual medio",
               subtitle = "Línea: mediana entre píxeles; banda: p10–p90"),
        "11_satelite_ciclo")

# Serie anual regional y ENSO
sat_reg <- sat_anual[, .(mm = mean(mm), pct_tope = mean(pct_tope)), by = anio]
sat_reg[, anom_pct := round(100 * (mm / mean(mm) - 1), 1)]
sat_reg <- oni[sat_reg, on = "anio"]
print(tabla(sat_reg, "11_satelite_anual")[, .(anio, mm = round(mm), anom_pct, pct_tope = round(pct_tope, 2), ONI_anual, fase)])
cat("r entre ONI y la anomalía anual del satélite:", round(sat_reg[, cor(ONI_anual, anom_pct)], 2), "\n")
grafico(ggplot(sat_reg, aes(anio, anom_pct, fill = fase)) + geom_col() +
          labs(x = NULL, y = "% respecto a la media", title = "Satélite: anomalía anual regional y fase ENSO"),
        "11_satelite_anomalia_anual")


# ---- 12. SATÉLITE FRENTE A ESTACIONES --------------------------------------------
titulo("12. Satélite frente a estaciones")
# Emparejamiento (P-7): cada estación incluida con el píxel que contiene su coordenada
ext_sat <- c(min(pix$lon) - 0.025, max(pix$lon) + 0.025, min(pix$lat) - 0.025, max(pix$lat) + 0.025)
col_px <- floor((estaciones$lon - ext_sat[1]) / 0.05) + 1L
fil_px <- floor((ext_sat[4] - estaciones$lat) / 0.05) + 1L
ncol_sat <- uniqueN(pix$lon); nfil_sat <- uniqueN(pix$lat)
dentro <- col_px >= 1 & col_px <= ncol_sat & fil_px >= 1 & fil_px <= nfil_sat
pares <- data.table(codigo = estaciones$codigo, pixel = ifelse(dentro, (fil_px - 1L) * ncol_sat + col_px, NA_integer_))
pares <- pix[, .(pixel, lon_px = lon, lat_px = lat, en_valle)][pares, on = "pixel"]
pares <- estaciones[, .(codigo, nombre, lon, lat, altitud, vertiente, nivel_c, revisar)][pares, on = "codigo"]
pares[, km_al_centro := round(dist_km(lat, lon, lat_px, lon_px), 2)]
stopifnot(pares[!is.na(pixel), all(abs(lon - lon_px) <= 0.0251 & abs(lat - lat_px) <= 0.0251)])
cat("Estaciones incluidas:", nrow(pares), "| fuera de la ventana del satélite:", pares[is.na(pixel), .N],
    "| en píxel con centro fuera del polígono:", pares[!is.na(pixel) & en_valle == FALSE, .N],
    "| píxeles distintos usados:", pares[!is.na(pixel), uniqueN(pixel)],
    "| píxeles con más de una estación:", pares[!is.na(pixel), .N, by = pixel][N > 1, .N], "\n")
if (pares[is.na(pixel), .N]) print(pares[is.na(pixel), .(codigo, nombre, lon, lat)])
tabla(pares, "12_pares_estacion_pixel")

# Estaciones a péntadas: solo péntadas con todos sus días observados
est_d <- diario[!is.na(P) & !imputado, .(codigo, anio, mes = month(fecha), dm = mday(fecha), P)]   # solo días observados
est_d[, pentada := (mes - 1L) * 6L + pmin((dm - 1L) %/% 5L, 5L) + 1L]
est_pen <- est_d[, .(mm_est = sum(P), n = .N), by = .(codigo, anio, pentada)]
est_pen <- cal[, .(anio, pentada, dias_dominio, dias_reales, dia_centro)][est_pen, on = .(anio, pentada)][n == dias_dominio]
# Píxel de cada estación (también si el centro del píxel queda fuera del polígono)
S_par <- S_todo[, cols]
cmp <- pares[!is.na(pixel), .(codigo, pixel)][est_pen, on = "codigo", nomatch = 0]
cmp[, mm_sat := S_par[cbind(pixel, match(sprintf("%d%02d", anio, pentada), cols))]]
cmp[, particion := PARTICION(anio)]
cat("Pares estación-péntada comparables:", nrow(cmp), "de", uniqueN(cmp$codigo), "estaciones\n")

# ¿Qué tan bien sigue el satélite a la estación, péntada a péntada?
print(cmp[, .(pares = .N, r = round(cor(mm_sat, mm_est), 2), mm_sat = round(mean(mm_sat), 1), mm_est = round(mean(mm_est), 1),
              sesgo_pct = round(100 * (sum(mm_sat) / sum(mm_est) - 1)), pct_est_cero = round(100 * mean(mm_est == 0), 1),
              pct_sat_cero = round(100 * mean(mm_sat == 0), 2), pct_est_sobre_200 = round(100 * mean(mm_est > 200), 2)), by = particion])
grafico(ggplot(cmp, aes(mm_est, mm_sat)) + geom_bin2d(bins = 60) + geom_abline(colour = "red") +
          scale_fill_viridis_c(trans = "log10", name = "péntadas") + coord_cartesian(xlim = c(0, 400)) +
          labs(x = "Estación (mm por péntada)", y = "Satélite (mm por péntada)",
               title = "Satélite frente a estación, péntada a péntada", subtitle = "Línea roja: igualdad"),
        "12_satelite_vs_estacion")

# El tope de 200 mm por péntada (P-15; el profesor pidió armar la serie y dibujarla).
# ¿Qué registra la estación cuando su píxel está en el tope? Si el tope fuera inocuo,
# la estación también estaría cerca de 200; si recorta, la estación pasaría de 200.
tope_est <- cmp[, .(pentadas = .N, sat_en_tope = sum(mm_sat == 200), est_sobre_200 = sum(mm_est > 200),
                    est_mediana_en_tope = if (any(mm_sat == 200)) median(mm_est[mm_sat == 200]) else NA_real_), by = codigo]
tope_est <- pares[, .(codigo, nombre, vertiente)][tope_est, on = "codigo"][order(-sat_en_tope)]
tabla(tope_est[sat_en_tope > 0], "12_tope_200_por_estacion")
print(tope_est[sat_en_tope > 0])
cat("Péntadas estación-satélite con el satélite en 200:", cmp[mm_sat == 200, .N], "de", nrow(cmp),
    "| en ellas la estación registra: mediana", cmp[mm_sat == 200, median(mm_est)], "mm, más de 200 en",
    cmp[mm_sat == 200, round(100 * mean(mm_est > 200))], "% |",
    "péntadas de estación sobre 200 mm:", cmp[mm_est > 200, .N], "(satélite en el tope en",
    cmp[mm_est > 200, round(100 * mean(mm_sat == 200))], "% de ellas)\n")
# La serie se dibuja con TODOS los años en que la estación midió (también los no válidos) y
# el satélite completo. Así un vacío en la línea de la estación es un vacío real del registro,
# y los años que R5 descarta quedan sombreados con su causa.
# Se dibuja Bahía Málaga, la serie por la que preguntó el profesor (P-15): la estación de la
# costa con más péntadas del satélite en el tope. Con R6 (OMM) ya no entra al modelo (empezó a
# medir en 1984), así que su píxel se ubica aquí con la misma regla de las incluidas.
ej_tope <- "5300000203"; st_tope <- valle[codigo == ej_tope]
px_tope <- (floor((ext_sat[4] - st_tope$lat) / 0.05)) * ncol_sat + floor((st_tope$lon - ext_sat[1]) / 0.05) + 1L
pen_de <- function(x) x[, pentada := (month(fecha) - 1L) * 6L + pmin((mday(fecha) - 1L) %/% 5L, 5L) + 1L][]
tope_dia <- pen_de(diario_todo[codigo == ej_tope & !is.na(P), .(anio, fecha, P)])
tope_est_p <- cal[, .(anio, pentada, dias_dominio, dia_centro)][tope_dia[, .(mm = sum(P), n = .N), by = .(anio, pentada)], on = .(anio, pentada)][n == dias_dominio]
tope_sat_p <- cal[, .(anio, pentada, dia_centro, mm = S_todo[px_tope, cols])]
con_tramo <- function(x) x[, `:=`(t = anio + (dia_centro - 1) / 365, tramo = cumsum(c(1L, diff((anio - 1981L) * 72L + pentada) != 1L)))][]
tope_est_p <- con_tramo(tope_est_p[order(anio, pentada)]); tope_sat_p <- con_tramo(tope_sat_p)
tope_anios <- anios[codigo == ej_tope][CJ(anio = ANIOS), on = "anio"]
tope_anios[, causa := fcase(is.na(n_obs), "sin registro", !valido, sprintf("%d días", n_obs), default = NA_character_)]
sombra <- tope_anios[!is.na(causa)]
# Años seguidos sin registro se juntan en un solo rótulo para que no se encimen
sombra[, grupo := rleid(causa, cumsum(c(1L, diff(anio) != 1L)))]
sombra_r <- sombra[, .(ini = min(anio), fin = max(anio) + 1L, causa = if (causa[1] == "sin registro" && .N > 1) "sin registro" else causa[1]), by = grupo]
sombra_r[, causa := sub(" ", "\n", causa)]
cat("Serie de", st_tope$nombre, "(incluida:", ej_tope %in% estaciones$codigo, "): años sin registro o no válidos:", paste(sombra$anio, sombra$causa, sep = " = ", collapse = "; "), "\n")
grafico(ggplot() +
          geom_rect(data = sombra_r, aes(xmin = ini, xmax = fin, ymin = -Inf, ymax = Inf), fill = "grey88") +
          geom_text(data = sombra_r, aes(x = (ini + fin) / 2, y = Inf, label = causa), vjust = 1.2, size = 2.5, lineheight = .85, colour = "grey30") +
          geom_line(data = tope_sat_p, aes(t, mm, colour = "Satélite (su píxel)", group = tramo), linewidth = .3) +
          geom_line(data = tope_est_p, aes(t, mm, colour = "Estación", group = tramo), linewidth = .3) +
          geom_hline(yintercept = 200, linetype = 2) +
          geom_point(data = tope_sat_p[mm == 200], aes(t, mm), colour = "red", size = 1) +
          scale_colour_manual(values = c("Estación" = "black", "Satélite (su píxel)" = "#1F78B4")) +
          labs(x = NULL, y = "mm por péntada", colour = NULL,
               title = sprintf("Serie péntada a péntada de %s (%s) y su píxel, 1981–2003", st_tope$nombre, ETQ_VERT[st_tope$vertiente]),
               subtitle = "Gris: años sin registro o que R5 no acepta (días medidos). Punteada: tope de 200 mm; puntos rojos: péntadas en el tope") +
          theme(legend.position = "bottom"),
        "12_serie_tope_200", ancho = 11, alto = 4.8)

# La serie real que se modela: lluvia DIARIA de una estación, año por año, con la lluvia
# diaria media de su píxel en cada péntada. Es el dato antes de volverlo curva; se escoge la
# estación incluida con más años válidos (empate: menor código).
ej_serie <- estaciones[codigo %in% pares[!is.na(pixel), codigo]][order(-validos, codigo)]$codigo[1]
px_serie <- pares[codigo == ej_serie, pixel]
serie_dia <- diario_todo[codigo == ej_serie, .(anio, dia, P)]
serie_sat <- cal[, .(anio, dia_centro, mm_dia = S_todo[px_serie, cols] / dias_reales)]
fase_de <- setNames(oni$fase, oni$anio)
lim_serie <- c(0, max(serie_dia$P, na.rm = TRUE))
grafico_r("6_serie_diaria_por_anio", function() {
  par(mfrow = c(6, 4), mar = c(3, 3.2, 2, .6), mgp = c(1.9, .5, 0), oma = c(0, 0, 3, 0))
  for (a in ANIOS) {
    x <- serie_dia[anio == a]; v <- anios[codigo == ej_serie & anio == a, valido]
    etq <- if (!nrow(x)) "sin registro" else if (!isTRUE(v)) "no válido (R5)" else ""
    plot(NA, xlim = c(1, 365), ylim = lim_serie, xaxt = "n", xlab = "Día del año", ylab = "mm/día",
         main = sprintf("%d · %s%s", a, fase_de[as.character(a)], if (nzchar(etq)) paste(" ·", etq) else ""), cex.axis = .75, cex.lab = .8)
    axis(1, at = c(1, 91, 182, 274, 365), labels = c("ene", "abr", "jul", "oct", "dic"), cex.axis = .75)
    if (nrow(x)) { lines(x$dia, x$P, type = "h", col = "grey45"); f <- x[is.na(P), dia]; if (length(f)) rug(f, col = "black", lwd = 1) }
    s <- serie_sat[anio == a]; lines(s$dia_centro, s$mm_dia, type = "s", lwd = 1.6)
  }
  plot.new(); legend("center", c("Estación: lluvia de cada día", "Satélite: lluvia diaria media\nde la péntada en su píxel", "Día sin dato (marca abajo)"),
                     col = c("grey45", "black", "black"), lty = c(1, 1, NA), lwd = c(1, 1.6, NA), pch = c(NA, NA, 124), bty = "n", cex = 1)
  mtext(sprintf("Lluvia diaria observada de %s (%s) y de su píxel, 1981–2003", estaciones[codigo == ej_serie, nombre], ETQ_VERT[estaciones[codigo == ej_serie, vertiente]]),
        outer = TRUE, cex = 1.1, font = 2)
})
cat("Serie diaria de ejemplo:", estaciones[codigo == ej_serie, paste(nombre, codigo)], "| años válidos:", estaciones[codigo == ej_serie, validos], "\n")

# Por estación: correlación y sesgo; ¿de qué dependen?
por_est_sat <- cmp[, .(pentadas = .N, r = cor(mm_sat, mm_est), sesgo = sum(mm_sat) / sum(mm_est),
                       mm_anual_est = mean(mm_est / dias_reales) * 365, mm_anual_sat = mean(mm_sat / dias_reales) * 365), by = codigo]
por_est_sat <- pares[, .(codigo, nombre, altitud, vertiente, revisar, nivel_c)][por_est_sat, on = "codigo"]
tabla(por_est_sat, "12_satelite_por_estacion")
print(por_est_sat[, .(estaciones = .N, r_mediana = round(median(r), 2), r_p10 = round(quantile(r, .1), 2),
                      sesgo_mediano = round(median(sesgo), 2), sobreestima = sum(sesgo > 1.25), subestima = sum(sesgo < 0.8)), by = vertiente])
cat("r entre el sesgo y la lluvia anual de la estación (log):", round(por_est_sat[, cor(log(sesgo), log(mm_anual_est))], 2), "\n")
print(por_est_sat[order(r)][1:8, .(codigo, nombre, vertiente, r = round(r, 2), sesgo = round(sesgo, 2), mm_anual_est = round(mm_anual_est), revisar)])
grafico(ggplot(por_est_sat, aes(mm_anual_est, mm_anual_sat, colour = ETQ_VERT[vertiente])) + geom_point() + geom_abline() +
          scale_x_log10() + scale_y_log10() + scale_colour_manual(values = COL_VERT) + theme(legend.position = "bottom") +
          labs(x = "Lluvia anual medida en la estación (mm/año, log)", y = "Lluvia anual del satélite en su píxel (mm/año, log)", colour = NULL,
               subtitle = "Cada punto es una estación. Sobre la diagonal, el satélite da más lluvia que la estación; debajo, menos.",
               title = "El satélite aplana: sobreestima donde llueve poco y subestima donde llueve mucho"),
        "12_sesgo_por_estacion")
grafico(ggplot(por_est_sat, aes(r, fill = ETQ_VERT[vertiente])) + geom_histogram(bins = 30) + scale_fill_manual(values = COL_VERT) +
          labs(x = "r entre satélite y estación (péntadas)", y = "Estaciones", fill = NULL, title = "Correlación por estación"),
        "12_correlacion_por_estacion")

# Ciclo anual: satélite y estación en los mismos pares
ciclo_cmp <- cmp[, .(`Estación` = mean(mm_est / dias_reales), `Satélite` = mean(mm_sat / dias_reales)), by = .(dia_centro)]
grafico(ggplot(melt(ciclo_cmp, id.vars = "dia_centro"), aes(dia_centro, value, colour = variable)) + geom_line() +
          labs(x = "Día del año", y = "mm/día", colour = NULL, title = "Ciclo anual medio en los pares estación–píxel"),
        "12_ciclo_satelite_estacion")

# Por año: ¿hay años en que el satélite se aparta de las estaciones? (D-16)
# Anomalía de cada par respecto a su propia media, con las mismas péntadas en las dos fuentes.
an_par <- cmp[, .(sat = mean(mm_sat), est = mean(mm_est), n = .N), by = .(codigo, anio)][n >= 60]
an_par[, `:=`(a_sat = sat / mean(sat) - 1, a_est = est / mean(est) - 1), by = codigo]
an_cmp <- an_par[, .(estaciones = .N, anom_sat_pct = round(100 * median(a_sat), 1), anom_est_pct = round(100 * median(a_est), 1),
                     r_espacial = round(cor(sat, est), 2)), by = anio][order(anio)]
an_cmp[, dif := anom_sat_pct - anom_est_pct]
an_cmp[, z := round((dif - median(dif)) / mad(dif), 1)]
an_cmp[, se_aparta := abs(z) > 2.5]
an_cmp <- oni[an_cmp, on = "anio"]
an_cmp[, particion := PARTICION(anio)]
print(tabla(an_cmp, "12_satelite_vs_estaciones_anual"))
cat("r entre las anomalías anuales del satélite y de las estaciones:", round(an_cmp[, cor(anom_sat_pct, anom_est_pct)], 2),
    "| desviación estándar: satélite", round(sd(an_cmp$anom_sat_pct), 1), "%, estaciones", round(sd(an_cmp$anom_est_pct), 1), "%\n")
# ¿De qué depende la diferencia? El satélite se queda corto en La Niña y se pasa en El Niño
cat("r con el ONI: estaciones", round(an_cmp[, cor(ONI_anual, anom_est_pct)], 2), "| satélite", round(an_cmp[, cor(ONI_anual, anom_sat_pct)], 2),
    "| diferencia satélite − estaciones", round(an_cmp[, cor(ONI_anual, dif)], 2), "\n")
print(an_cmp[order(dif), .(anio, fase, ONI_anual, anom_sat_pct, anom_est_pct, dif, z)][c(1:5, (.N - 4):.N)])
grafico(ggplot(an_cmp, aes(ONI_anual, dif)) + geom_text(aes(label = anio, colour = fase), size = 3.2) +
          geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey50", linewidth = .5) +
          labs(x = "ONI anual", y = "Anomalía del satélite − anomalía de las estaciones (puntos %)",
               title = "El satélite se queda corto en La Niña y se pasa en El Niño"),
        "12_diferencia_vs_ONI")

# Por año: el total anual de la estación y el de su píxel en los mismos pares (solo
# estación-años con las 72 péntadas completas, para sumar exactamente las mismas péntadas)
# y la correlación estación–satélite de cada año. Son resúmenes escalares del año; la lluvia
# se analiza en mm/día, nunca acumulada.
completos <- cmp[, .N, by = .(codigo, anio)][N == 72, .(codigo, anio)]
tot <- cmp[completos, on = .(codigo, anio)][, .(est = sum(mm_est), sat = sum(mm_sat)), by = .(codigo, anio)]
acu_anio <- tot[, .(estaciones = .N, P_est = mean(est), P_sat = mean(sat)), by = anio]
acu_anio <- oni[acu_anio, on = "anio"][order(anio)]
acu_anio[, sat_sobre_est_pct := round(100 * P_sat / P_est)]
r_anio <- oni[cmp[, .(r = cor(mm_sat, mm_est), pares = .N), by = anio], on = "anio"][order(anio)]
acu_anio <- r_anio[, .(anio, r_pentadas = round(r, 2))][acu_anio, on = "anio"]
print(tabla(acu_anio[, .(anio, fase, ONI_anual, estaciones, P_est = round(P_est), P_sat = round(P_sat), sat_sobre_est_pct, r_pentadas)], "12_por_anio"))
cat(sprintf("Total anual del satélite / estación: %d a %d %% | r con el ONI: %.2f | r estación–satélite por año: %.2f a %.2f (mediana %.2f)\n",
            min(acu_anio$sat_sobre_est_pct), max(acu_anio$sat_sobre_est_pct), cor(acu_anio$sat_sobre_est_pct, acu_anio$ONI_anual),
            min(r_anio$r), max(r_anio$r), median(r_anio$r)))
grafico_r("12_correlacion_por_anio", function() {
  par(mar = c(4, 6, 3, 1)); tonos <- c("Niña" = "grey20", "Neutro" = "grey60", "Niño" = "grey90")
  b <- barplot(r_anio$r, names.arg = r_anio$anio, col = tonos[r_anio$fase], ylim = c(0, 1), las = 2, cex.names = .8,
               ylab = "Correlación estación–satélite\n(péntadas del año)", main = "Correlación entre estación y satélite, año por año")
  abline(h = median(r_anio$r), lty = 3); text(b, r_anio$r + .03, sprintf("%.2f", r_anio$r), cex = .65)
  legend("topright", names(tonos), fill = tonos, bty = "n", horiz = TRUE, title = "Fase ENSO")
}, ancho = 10, alto = 4.5)

# ¿De dónde sale lo que el satélite sí acierta? Se separan tres componentes
cmp[, mes := (pentada - 1L) %/% 6L + 1L]
comp_mes <- cmp[, .(sat = mean(mm_sat), est = mean(mm_est)), by = .(codigo, anio, mes)]
comp_mes[, `:=`(a_sat = sat - mean(sat), a_est = est - mean(est)), by = .(codigo, mes)]
comp_ciclo <- cmp[, .(sat = mean(mm_sat), est = mean(mm_est)), by = .(codigo, pentada)][, .(r = cor(sat, est)), by = codigo]
comp_esp <- cmp[, .(sat = mean(mm_sat), est = mean(mm_est)), by = codigo]
aciertos <- data.table(
  componente = c("Dónde llueve más (media por estación)", "Cuándo llueve en el año (ciclo medio, mediana entre estaciones)",
                 "Meses más o menos lluviosos de lo normal (anomalías)", "Años más o menos lluviosos de lo normal (regional)"),
  r = round(c(cor(comp_esp$sat, comp_esp$est), median(comp_ciclo$r), comp_mes[, cor(a_sat, a_est)],
              an_cmp[, cor(anom_sat_pct, anom_est_pct)]), 2))
print(tabla(aciertos, "12_que_acierta_el_satelite"))
grafico(ggplot(melt(an_cmp[, .(anio, `Satélite` = anom_sat_pct, Estaciones = anom_est_pct)], id.vars = "anio"),
               aes(anio, value, colour = variable)) + geom_line() + geom_point() +
          geom_vline(xintercept = an_cmp[se_aparta == TRUE, anio], linetype = 3) +
          geom_vline(xintercept = c(1995.5, 2000.5), colour = "red", linetype = 2) +
          labs(x = NULL, y = "Anomalía anual (%, mediana de los pares)", colour = NULL,
               title = "Años en que el satélite se aparta de las estaciones",
               subtitle = "Punteado: años que más se apartan; rojo: límites de entrenamiento, validación y prueba"),
        "12_anomalia_anual_satelite_estaciones")
# ¿En qué meses falla el satélite en los años que más se apartan?
if (an_cmp[se_aparta == TRUE, .N]) {
  mes_sosp <- cmp[anio %in% an_cmp[se_aparta == TRUE, anio], .(sat = sum(mm_sat), est = sum(mm_est)), by = .(anio, mes = (pentada - 1L) %/% 6L + 1L)]
  mes_ref <- cmp[!anio %in% an_cmp[se_aparta == TRUE, anio], .(razon_normal = sum(mm_sat) / sum(mm_est)), by = .(mes = (pentada - 1L) %/% 6L + 1L)]
  mes_sosp <- mes_ref[mes_sosp, on = "mes"][, razon := round(sat / est / razon_normal, 2)]
  cat("Satélite/estación por mes en los años que más se apartan, relativo al resto (1 = igual que siempre):\n")
  print(dcast(mes_sosp, anio ~ mes, value.var = "razon"))
}

# QUÉ DICEN LOS DATOS (secciones 10–12)
#   Las cifras y su interpretación están en EXPLICACION_DATOS.md, sección 12 (corrida del
#   19-nov-2026, con las 98 estaciones de la regla de la OMM); las tablas de esta
#   sección las imprimen las líneas de arriba.

saveRDS(list(pentadal = sat, pixeles = pix, calendario = cal, pares = pares, pares_pentada = cmp,
             anual = an_cmp, por_estacion = por_est_sat, que_acierta = aciertos),
        file.path(SALIDA, "satelite_valle.rds"), compress = "xz")
rm(S_todo, S_par, est_d); invisible(gc())


# ---- 13. DATOS PREPARADOS PARA LA CONSTRUCCIÓN FUNCIONAL ------------------------
titulo("13. Datos preparados")
# Todo lo que necesita la construcción funcional, ya en forma de matriz:
#   estaciones  Y: estación-año × día (1–365), mm/día, NA = faltante
#   satélite    X: píxel-año × péntada (1–72), mm/día, sin faltantes; t = día central
#   pares       fila de Y ↔ fila de X (el píxel de la estación, el mismo año)

# Estaciones
id_est <- unique(diario[, .(codigo, anio, particion)])
setkey(id_est, codigo, anio)
Yw <- dcast(diario, codigo + anio ~ dia, value.var = "P")
stopifnot(identical(Yw$codigo, id_est$codigo), identical(Yw$anio, id_est$anio), identical(names(Yw)[-(1:2)], as.character(1:365)))
Y <- as.matrix(Yw[, -(1:2)]); dimnames(Y) <- list(paste(id_est$codigo, id_est$anio, sep = "_"), 1:365); rm(Yw)
Yimp <- as.matrix(dcast(diario, codigo + anio ~ dia, value.var = "imputado")[, -(1:2)]); dimnames(Yimp) <- dimnames(Y)   # TRUE = rellenado por R7
Yobs <- !is.na(Y) & !Yimp                                                                                                  # TRUE = observado
id_est <- anios[, .(codigo, anio, n_obs, hueco_max)][id_est, on = .(codigo, anio)]
id_est[, dias_rellenados := rowSums(Yimp)]
id_est[, fila := .I]

# Satélite
id_sat <- unique(sat[, .(pixel, anio)]); setkey(id_sat, pixel, anio)
Xw <- dcast(sat, pixel + anio ~ pentada, value.var = "mm_dia")
stopifnot(identical(Xw$pixel, id_sat$pixel), identical(Xw$anio, id_sat$anio))
X <- as.matrix(Xw[, -(1:2)]); dimnames(X) <- list(paste(id_sat$pixel, id_sat$anio, sep = "_"), 1:72); rm(Xw)
id_sat[, `:=`(particion = PARTICION(anio), fila = .I)]
t_sat <- cal[anio == ANIOS[1], dia_centro]
stopifnot(cal[, uniqueN(dia_centro), by = pentada][, all(V1 == 1)])   # el mismo día central todos los años

# Pares estación-año ↔ píxel-año
pares_fun <- pares[, .(codigo, pixel)][id_est[, .(codigo, anio, particion, fila_est = fila)], on = "codigo"]
pares_fun <- id_sat[, .(pixel, anio, fila_sat = fila)][pares_fun, on = .(pixel, anio)]
stopifnot(!anyNA(pares_fun$fila_sat))
# ONI del año (NOAA CPC, D-19): el enunciado lo agrega al vector z de la red, junto con la altitud
# (que está en meta_estaciones)
pares_fun <- oni[, .(anio, ONI = ONI_anual, fase_enso = fase)][pares_fun, on = "anio"]
setcolorder(pares_fun, c("codigo", "anio", "particion", "pixel", "fila_est", "fila_sat", "ONI", "fase_enso")); setorder(pares_fun, fila_est)

# Metadatos de las estaciones, con las marcas para el análisis de sensibilidad
meta_est <- por_est_sat[, .(codigo, r_satelite = round(r, 3), sesgo_satelite = round(sesgo, 3))][
  pares[, .(codigo, nombre, lon, lat, altitud, vertiente, pixel, km_al_centro,
            cruce_CNE = nivel_c == "Cruce con CNE", coordenada_a_revisar = revisar != "")], on = "codigo"]
meta_est <- por_est[, .(codigo, validos_ent, validos_val, validos_pru)][meta_est, on = "codigo"]
meta_est <- clim[, .(codigo, mm_anual = round(P_media))][meta_est, on = "codigo"]

preparados <- list(
  estaciones = list(Y = Y, rellenado = Yimp, id = id_est, t = 1:365,
                    unidades = "mm/día; NA = día faltante; rellenado = TRUE si el día lo estimó el kriging (R7)"),
  conjunto = CONJUNTO,
  satelite   = list(X = X, id = id_sat, t = t_sat, dias = cal[anio == ANIOS[1], dias_dominio], unidades = "mm/día; t = día central de la péntada"),
  pares = pares_fun, meta_estaciones = meta_est, pixeles = pix[usar == TRUE],
  particion = c(entrenamiento = "1981–1995", validacion = "1996–2000", prueba = "2001–2003"),
  reglas = reglas)
saveRDS(preparados, file.path(SALIDA, "datos_preparados.rds"), compress = "xz")
tabla(meta_est, "13_meta_estaciones"); tabla(id_est, "13_id_estacion_anio")

print(rbind(
  id_est[, .(fuente = "Estaciones (curvas estación-año)", curvas = .N, sitios = uniqueN(codigo)), by = particion],
  id_sat[, .(fuente = "Satélite (curvas píxel-año)", curvas = .N, sitios = uniqueN(pixel)), by = particion],
  pares_fun[, .(fuente = "Pares estación-año ↔ píxel-año", curvas = .N, sitios = uniqueN(codigo)), by = particion]))
cat("Y:", nrow(Y), "×", ncol(Y), "| faltantes:", sum(is.na(Y)), "(", round(100 * mean(is.na(Y)), 2), "% )",
    "| rellenados por kriging:", sum(Yimp), "(", round(100 * mean(Yimp), 3), "% )",
    "\nX:", nrow(X), "×", ncol(X), "| faltantes:", sum(is.na(X)), "\n")
cat("Marcas de sensibilidad: cruce CNE", meta_est[cruce_CNE == TRUE, .N], "| coordenada a revisar", meta_est[coordenada_a_revisar == TRUE, .N],
    "| sin altitud", meta_est[is.na(altitud), .N], "\n")

# Vista de las curvas crudas que entran a la construcción funcional
set.seed(7); ej <- sample(id_est[particion == "entrenamiento", fila], 4)
ej_dt <- rbindlist(lapply(ej, function(i) {
  j <- pares_fun[fila_est == i, fila_sat]
  rbind(data.table(curva = rownames(Y)[i], fuente = "Estación (diaria)", t = 1:365, v = Y[i, ]),
        data.table(curva = rownames(Y)[i], fuente = "Satélite (péntada)", t = t_sat, v = X[j, ]))
}))
grafico(ggplot(ej_dt, aes(t, v, colour = fuente)) + geom_line(alpha = .8) + facet_wrap(~curva, scales = "free_y") +
          labs(x = "Día del año", y = "mm/día", colour = NULL, title = "Ejemplos de pares estación-año y su píxel",
               subtitle = "Así entran los datos a la construcción funcional"),
        "13_ejemplos_pares", alto = 6)


# ---- 14. PRUEBAS -----------------------------------------------------------------
titulo("14. Pruebas")
# Comprobaciones independientes de los datos preparados. Si alguna falla, el script se detiene.
pruebas <- list()
prueba <- function(nombre, ok, detalle = "") {
  pruebas[[length(pruebas) + 1L]] <<- data.table(prueba = nombre, resultado = ifelse(isTRUE(ok), "OK", "FALLA"), detalle = as.character(detalle))
}

# -- Población
prueba("Las 232 del Valle caen dentro del polígono DANE", est[grupo == "Valle del Cauca", .N == 232 && all(depto_geo == "Valle del Cauca")])
prueba("Ninguna estación final es del Cauca ni está sin ubicar", all(id_est$codigo %in% est[grupo == "Valle del Cauca", codigo]))
prueba("Las incluidas tienen coordenada", !anyNA(meta_est[, .(lon, lat)]))
prueba("Las incluidas tienen altitud (entra al vector z)", !anyNA(meta_est$altitud), paste("Univalle:", meta_est[codigo == "2622320201", altitud], "m"))

# -- Periodo y partición
prueba("Solo años 1981–2003", all(id_est$anio %in% ANIOS) && all(id_sat$anio %in% ANIOS))
prueba("Partición correcta por año", id_est[, all((particion == "entrenamiento") == (anio <= 1995) & (particion == "prueba") == (anio >= 2001))])
prueba(paste("R6: cada estación tiene años válidos", R6_TEXTO), id_est[, .N, by = codigo][, min(N) >= ANIOS_MIN],
       paste("mínimo:", id_est[, .N, by = codigo][, min(N)], "de", length(ANIOS)))

# -- Estaciones: reglas recalculadas desde la matriz (R5 se mide con los días OBSERVADOS)
n_obs_Y <- rowSums(Yobs); hueco_Y <- apply(Yobs, 1, hueco_mas_largo)
falt_mes_Y <- apply(Yobs, 1, function(o) max(tapply(!o, MES_DIA, sum)))
seg_mes_Y <- apply(Yobs, 1, function(o) max(vapply(split(o, MES_DIA), hueco_mas_largo, 0L)))
prueba("R5 (OMM): ningún mes de una curva tiene 11 o más días faltantes", all(falt_mes_Y <= MES_FALT_MAX), paste("máximo:", max(falt_mes_Y)))
prueba("R5 (OMM): ningún mes de una curva tiene 5 o más días faltantes seguidos", all(seg_mes_Y <= MES_SEGUID_MAX), paste("máximo:", max(seg_mes_Y)))
prueba("R7: solo se rellenan días de años válidos y nunca con valores negativos",
       imputados[!is.na(P), all(P >= 0)] && nrow(imputados[!is.na(P)][!anios[valido == TRUE], on = .(codigo, anio)]) == 0 && all(id_est$n_obs + id_est$dias_rellenados == rowSums(!is.na(Y))),
       sprintf("%d días rellenados (%.3f %% de la matriz)", sum(Yimp), 100 * mean(Yimp)))
prueba("R7: el semivariograma y las medias mensuales usan solo entrenamiento", max(R7_ENT) == sum(ANIOS <= 1995) * 365L,
       paste("modelo:", vg_mejor))
prueba("Sin duplicados estación-año", !anyDuplicated(id_est[, .(codigo, anio)]))
prueba("R0: ningún valor negativo ni sobre el récord OMM de 24 h", min(Y, na.rm = TRUE) >= 0 && max(Y, na.rm = TRUE) <= P_MAX_FISICO,
       sprintf("en el Valle R0 quitó %d valores; en el crudo completo habría quitado %d", nrow(m_r0), sum(r0_resumen[1, .(negativos, sobre_maximo_fisico)])))
prueba("Valores entre 0 y el máximo crudo", min(Y, na.rm = TRUE) >= 0 && max(Y, na.rm = TRUE) <= estructura$maximo, paste("máximo:", max(Y, na.rm = TRUE)))

# -- Estaciones: contra los archivos crudos, leídos de nuevo
crudo2 <- rbind(fread("datos_curso/P_1980_1990.csv"), fread("datos_curso/P_1991_2003.csv"))
crudo2 <- crudo2[YEAR %in% ANIOS & !(MONTH == 2 & DAY == 29), .(codigo = as.character(STATION_ID), anio = YEAR, mes = MONTH, dm = DAY, P_crudo = as.numeric(VALUE))]
crudo2[, dia := c(0L, cumsum(dias_mes_nb))[mes] + dm]
crudo2 <- id_est[, .(codigo, anio, fila)][crudo2, on = .(codigo, anio), nomatch = 0]
crudo2[, P_final := Y[cbind(fila, dia)]]
m2 <- marcas[, .(codigo, anio = year(fecha), mes = month(fecha), dm = mday(fecha), marcado = TRUE)]
crudo2 <- m2[crudo2, on = .(codigo, anio, mes, dm)]
crudo2[, rellenado := Yimp[cbind(fila, dia)]]
prueba("Cada valor final es igual al crudo (salvo los marcados)", crudo2[is.na(marcado), all(!is.na(P_final) & P_final == P_crudo & !rellenado)],
       paste(crudo2[is.na(marcado), .N], "valores comparados"))
prueba("Los días marcados quedaron como faltante (o los rellenó R7)", crudo2[marcado == TRUE, all(is.na(P_final) | rellenado)],
       paste(crudo2[marcado == TRUE, .N], "marcados dentro de años válidos;", crudo2[marcado == TRUE & rellenado, .N], "rellenados"))
prueba("Todo día observado de la matriz es un dato del crudo no marcado", sum(Yobs) == crudo2[is.na(marcado), .N],
       paste("observados:", sum(Yobs), "| rellenados:", sum(Yimp)))
rm(crudo2)

prueba("La trazabilidad termina en los datos preparados", trazabilidad[.N, estaciones == uniqueN(id_est$codigo) && pares_estacion_anio == nrow(Y) && dias_con_dato == sum(!is.na(Y))])

# -- Umbrales sin fuga
p99_chk <- d[P > 0 & anio <= 1995, .(p = quantile(P, .99), n = .N), by = codigo][n >= 100]
prueba("El p99 de las reglas usa solo entrenamiento", p99[p99_chk, on = "codigo", all(abs(p99 - p) < 1e-9)], paste(nrow(p99_chk), "estaciones"))

# -- Satélite
prueba("Satélite sin faltantes y no negativo", !anyNA(X) && min(X) >= 0)
prueba("72 péntadas que suman 365 días", length(t_sat) == 72 && sum(preparados$satelite$dias) == 365 && all(diff(t_sat) > 0))
prueba("mm/día × días reales = mm de la péntada (enteros)", sat[, all(abs(mm_dia * dias_reales - round(mm_dia * dias_reales)) < 1e-9)])
# contra los TIFF, leídos de nuevo
set.seed(11); chk <- inv_sat[anio %in% ANIOS][sample(.N, 12)]
ok_tif <- vapply(seq_len(nrow(chk)), function(i) {
  v <- terra::values(terra::rast(chk$archivo[i]))[, 1]
  dr <- cal[anio == chk$anio[i] & pentada == chk$pentada[i], dias_reales]
  f <- id_sat[anio == chk$anio[i]]
  all(abs(X[f$fila, chk$pentada[i]] * dr - v[f$pixel]) < 1e-9)
}, TRUE)
prueba("La matriz del satélite coincide con 12 TIFF releídos", all(ok_tif))
prueba("Todos los píxeles del Valle están", all(pix[en_valle == TRUE, pixel] %in% id_sat$pixel), paste(uniqueN(id_sat$pixel), "píxeles"))

# -- Pares
prueba("Cada estación cae dentro de su píxel", meta_est[pix, on = "pixel", nomatch = 0][, all(abs(lon - i.lon) <= 0.0251 & abs(lat - i.lat) <= 0.0251)],
       paste("distancia máxima al centro:", max(meta_est$km_al_centro), "km"))
prueba("Cada par une el mismo año y el píxel de la estación",
       all(id_sat$anio[pares_fun$fila_sat] == id_est$anio[pares_fun$fila_est]) &&
       all(id_sat$pixel[pares_fun$fila_sat] == meta_est[match(id_est$codigo[pares_fun$fila_est], codigo), pixel]))
prueba("Ninguna altitud es imposible en la cuenca del Cauca (D-59)",
       valle[substr(codigo, 1, 2) == "26" & !is.na(altitud), all(altitud >= PISO_CAUCA - 300)],
       paste(nrow(alt_imposible), "altitud(es) oficial(es) reemplazada(s) por el terreno; piso", PISO_CAUCA, "m"))
prueba("Hay un par por cada curva de estación", nrow(pares_fun) == nrow(Y) && !anyDuplicated(pares_fun$fila_est))
prueba("Cada par tiene el ONI de su año (entra al vector z)", !anyNA(pares_fun$ONI) && pares_fun[oni, on = "anio", nomatch = 0][, all(ONI == ONI_anual)], paste(uniqueN(pares_fun$anio), "años"))

# -- Agregación a péntadas
comp <- which(n_obs_Y == 365)[1:50]   # años con los 365 días observados
ok_pen <- vapply(comp, function(i) {
  tot <- est_pen[codigo == id_est$codigo[i] & anio == id_est$anio[i], sum(mm_est)]
  abs(tot - sum(Y[i, ])) < 1e-6
}, TRUE)
prueba("En años completos, la suma por péntadas = la suma diaria", all(ok_pen), "50 curvas")

pruebas <- rbindlist(pruebas)
print(tabla(pruebas, "14_pruebas"), right = FALSE)
if (any(pruebas$resultado != "OK")) stop("Hay pruebas que fallan: revisar la tabla anterior")
cat("\nTodas las pruebas pasan (", nrow(pruebas), ").\n")


# ---- 15. RESUMEN --------------------------------------------------------------
titulo("15. Resumen")
cat("Días pasados a faltante:", nrow(marcas), "\n")
cat("Años válidos:", anios[valido == TRUE, .N], "de", nrow(anios), "\n")
print(anios[, .(validos = sum(valido), total = .N), by = particion])
cat("Estaciones con datos:", nrow(valle), "| incluidas:", nrow(estaciones),
    "(CNE:", estaciones[nivel_c == "Cruce con CNE", .N], "| coordenada a revisar:", estaciones[revisar != "", .N], ")\n")
cat("Con ≥1 año válido en validación:", estaciones[validos_val >= 1, .N], "| en prueba:", estaciones[validos_pru >= 1, .N],
    "| en ambos:", estaciones[validos_val >= 1 & validos_pru >= 1, .N], "\n")
print(diario[, .(estacion_anios = uniqueN(paste(codigo, anio)), dias = .N, faltantes = sum(is.na(P))), by = particion])
cat("\nDataset final:", file.path(SALIDA, "lluvia_valle_limpia.rds"), "\n",
    "  $diario       incluidas × años válidos, día 1–365, P en mm (NA = faltante)\n",
    "  $diario_todo  las", nrow(valle), "estaciones y todos sus años, ya depurados\n",
    "  $anios        validez de cada estación-año  |  $estaciones  metadatos de las incluidas\n",
    "  $marcas       días pasados a faltante y por qué  |  $reglas  texto de las reglas\n")
cat("\nSatélite:", file.path(SALIDA, "satelite_valle.rds"), "\n",
    "  $pentadal       ", uniqueN(sat$pixel), "píxeles del Valle ×", nrow(cal), "péntadas, en mm por péntada y mm/día\n",
    "  $pixeles, $calendario   centro de cada píxel y días de cada péntada\n",
    "  $pares          estación ↔ píxel (", pares[!is.na(pixel), .N], "estaciones en", pares[!is.na(pixel), uniqueN(pixel)], "píxeles )\n",
    "  $pares_pentada  estación y satélite en las mismas péntadas  |  $anual, $por_estacion, $que_acierta\n")
cat("\nDatos preparados:", file.path(SALIDA, "datos_preparados.rds"), "\n",
    "  $estaciones$Y   ", nrow(Y), "curvas estación-año × 365 días (mm/día, NA = faltante); $estaciones$id las identifica\n",
    "  $satelite$X     ", nrow(X), "curvas píxel-año × 72 péntadas (mm/día); $satelite$t = día central\n",
    "  $pares          fila de Y ↔ fila de X  |  $meta_estaciones, $pixeles, $particion, $reglas\n",
    "  Pruebas:", nrow(pruebas), "de", nrow(pruebas), "pasan\n")
writeLines(capture.output(sessionInfo()), file.path(SALIDA, "sessionInfo.txt"))

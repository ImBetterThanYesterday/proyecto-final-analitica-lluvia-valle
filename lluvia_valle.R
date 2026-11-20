# =============================================================================
# LLUVIA EN EL VALLE DEL CAUCA 1981-2003 — ARCHIVO ÚNICO DE R
# Proyecto final de Analítica de Datos (UAO, 2026-2S)
# =============================================================================
# Este archivo es la unión, SIN CAMBIOS en el código, de:
#   PARTE A = lluvia_valle_datos.R       (datos crudos → datos preparados)
#   PARTE B = lluvia_valle_funcional.R   (todo lo funcional y la evaluación de la red)
# Solo se agregaron este encabezado, la separación entre partes y las marcas
# "RÚBRICA x" antes de cada sección, para ubicar cada criterio.
#
# Mapa de la rúbrica:
#   1    Construcción de datos funcionales ........ Parte B, sección 2
#   2    EDA funcional ............................ Parte B, sección 3 (y Parte A, sección 11)
#   3    Comparación funcional por fuente ......... Parte B, secciones 1 y 4 (y Parte A, sección 12)
#   4.1  Preprocesamiento y selección de datos .... Parte A, secciones 1-8, 13 y 14
#   4.2  Red neuronal funcional ................... Parte B, sección 5 (regresión = cabeza lineal) y lluvia_valle_red.py
#   4.3  Ablación KAN, MLP y cabeza lineal ........ lluvia_valle_red.py (entrenamiento) y Parte B, sección 7
#   4.4  Evaluación ............................... Parte B, secciones 5.5 y 7
#   4.5  Interpretación climatológica ............. Parte B, sección 8
#
# Orden de ejecución (desde la raíz del repositorio):
#   Rscript lluvia_valle.R          # partes A y B; deja salidas_R_funcional/para_python/
#   python lluvia_valle_red.py      # entrena la red y exporta sus predicciones
#   Rscript lluvia_valle.R          # otra vez: la Parte B evalúa (7) e interpreta (8) la red
# La Parte A necesita datos_curso/ (ver README). Pruebas: 27 (A) + 13 (B) + 5 (red).
# =============================================================================

# █████ PARTE A — DATOS (antes lluvia_valle_datos.R) █████

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
#   8. Depuración (reglas R0 a R7, D-36) y dataset final
#   9. Satélite: lectura, estructura y fuente (D-43)
#  10. Satélite: recorte al Valle, calendario de péntadas y mm/día
#  11. Satélite: EDA
#  12. Satélite frente a estaciones: emparejamiento, sesgo y años sospechosos
#  13. Datos preparados para la construcción funcional (matrices y pares)
#  14. Pruebas de los datos preparados
#  15. Resumen
#
# Ejecutar desde la raíz del repositorio:  Rscript lluvia_valle_datos.R
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
# Decisiones y su justificación: DECISIONES.md (D-01 a D-49).
# =============================================================================


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ┃ RÚBRICA — — configuración de la parte A (datos)
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
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
COB_MIN     <- 0.90                         # R5: fracción mínima de días en el año
HUECO_MAX   <- 15                           # R5: hueco máximo permitido (días)
ANIOS_MIN   <- 7                            # R6: años válidos mínimos en entrenamiento
# R0: máximo físico = récord mundial de lluvia en 24 h reconocido por la OMM: 1 825 mm,
# Foc-Foc (La Reunión), 7-8 ene 1966, ciclón Denise (WMO Weather and Climate Extremes
# Archive, https://wmo.asu.edu). Un dato diario por encima no puede ser lluvia real.
P_MAX_FISICO <- 1825

SALIDA <- "salidas_R_datos"
dir.create(file.path(SALIDA, "tablas"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(SALIDA, "graficos"), showWarnings = FALSE, recursive = TRUE)
tabla   <- function(x, nombre) { fwrite(x, file.path(SALIDA, "tablas", paste0(nombre, ".csv"))); invisible(x) }
grafico <- function(p, nombre, ancho = 9, alto = 5, dpi = 110) {
  ggsave(file.path(SALIDA, "graficos", paste0(nombre, ".png")), p, width = ancho, height = alto, dpi = dpi, bg = "white")
  if (interactive()) print(p)
  invisible(p)
}
titulo <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
# Gráficas con R base, en el estilo de la tesis del profesor (paneles por año, blanco y negro)
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


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ┃ RÚBRICA 4.1 — Preprocesamiento: población de estaciones del Valle
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
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
#   323 códigos con datos en 1981–2003. De las 232 del Valle, 220 tienen datos en
#   el periodo; las 12 restantes solo tienen datos fuera. Las 232 caen dentro del
#   polígono DANE del Valle (13 a ≤2 km de otro departamento).


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ┃ RÚBRICA 4.1 — Preprocesamiento: estructura de los crudos y R0 (valores imposibles)
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
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
#   Un faltante es una FILA AUSENTE: no hay NA, negativos ni duplicados, y no
#   aparece el valor imposible (859 375 mm) de otra estación. R0 lo comprueba con
#   un umbral citado (récord OMM, 1 825 mm en 24 h): en el Valle no quita nada. Hay 832 registros
#   del 29-feb y 5 144 valores con decimales (estaciones climatológicas).


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ┃ RÚBRICA 4.1 — Preprocesamiento: días faltantes, cobertura y huecos
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
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
#   El 81 % de los estación-años con datos tiene ≥99 % de días. Los huecos
#   internos son el 2,7 % de la rejilla y son cortos (mediana de 3 días, p90 de
#   31, 28 % de un solo día). Con ≥90 %, 180 estaciones tienen ≥7 años en
#   entrenamiento. El umbral definitivo sale del experimento de la sección 7.


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ┃ RÚBRICA 4.1 — Preprocesamiento: errores de digitación y extremos (R2)
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
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
#   60,6 % de los días son ceros; p99 ≈ 76 mm; máximo 600 mm. Solo 5 días superan
#   3·p99 y 150 mm, y en todos las vecinas registran como máximo 78 mm (Pardo 600,
#   Heraclio Uribe 503, La Cascada 421, Costa Rica 283, Venecia 283) → regla R2.


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ┃ RÚBRICA 4.1 — Preprocesamiento: totales mensuales en un día (R3) y ceros falsos (R4)
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
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

# Prueba de cada racha de ≥30 días en cero (R4). Si la estación de verdad estuvo
# seca, su vecina más cercana también debió estar seca o casi. Se mide:
#   k = días en que llovió la vecina más cercana con datos (≥80 % de la racha);
#   q = P(llueve en la estación | llueve en la vecina), con los días de
#       entrenamiento (1981–1995) en que ambas tienen dato (sin fuga);
#   p = (1 - q)^k: probabilidad de que la estación marque 0 todos esos días.
# Se hacen muchas pruebas a la vez, así que el umbral es 0,05 / número de rachas
# (Bonferroni). Supone días independientes; como la lluvia viene en rachas, p
# exagera lo raro, y por eso el umbral es exigente.
R4_LARGO <- 30L
ll_ent <- d[anio <= 1995 & !is.na(P), .(codigo, fecha, ll = P > 0)]; setkey(ll_ent, codigo, fecha)
d_obs <- d[!is.na(P), .(codigo, fecha, P)]; setkey(d_obs, codigo, fecha)
r4 <- rz[largo >= R4_LARGO]
r4[, c("vecina", "km", "dias_vecina", "dias_lluvia_vecina", "q") := rbindlist(lapply(seq_len(.N), function(i) {
  cod <- codigo[i]; desde <- ini[i]; hasta <- fin[i]   # copias locales: dentro de d_obs[...] "codigo" sería su columna
  for (v in names(sort(DK[cod, vecinas(cod)]))) {
    x <- d_obs[.(v)][fecha %between% c(desde, hasta)]
    if (nrow(x) >= .8 * largo[i]) {
      a <- merge(ll_ent[.(cod)], ll_ent[.(v)], by = "fecha")
      return(data.table(v, round(DK[cod, v], 1), nrow(x), sum(x$P > 0), a[ll.y == TRUE, mean(ll.x)]))
    }
  }
  data.table(NA_character_, NA_real_, NA_integer_, NA_integer_, NA_real_)
}))]
R4_UMBRAL <- 0.05 / nrow(r4)
r4[, `:=`(p_azar = (1 - q)^dias_lluvia_vecina, nombre_vecina = valle[match(vecina, codigo), nombre])]
r4[, falso := !is.na(p_azar) & p_azar < R4_UMBRAL]
tabla(r4[order(p_azar)], "5b_rachas_ceros_prueba")
cat(sprintf("Rachas de ≥%d días en 0: %d (%d sin vecina con datos). Umbral 0,05/%d = %.1e → %d rachas falsas en %d estaciones\n",
            R4_LARGO, nrow(r4), r4[is.na(p_azar), .N], nrow(r4), R4_UMBRAL, r4[falso == TRUE, .N], r4[falso == TRUE, uniqueN(codigo)]))
print(r4[order(p_azar)][1:12, .(nombre, ini, fin, largo, vecina = nombre_vecina, km, dias_lluvia_vecina, q = round(q, 2), p_azar = signif(p_azar, 2), falso)])
grafico(ggplot(r4[!is.na(p_azar)], aes(dias_lluvia_vecina, pmax(p_azar, 1e-12), colour = falso)) + geom_point(alpha = .6) +
          geom_hline(yintercept = R4_UMBRAL, linetype = 2) + scale_y_log10() +
          scale_colour_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey50"), labels = c(`TRUE` = "Falsa → faltante", `FALSE` = "Se mantiene"), name = NULL) +
          labs(x = "Días con lluvia en la vecina más cercana durante la racha", y = "p = (1 - q)^k (log)",
               title = "R4: rachas de ≥30 días en 0 que la vecina no respalda", subtitle = sprintf("Línea: umbral 0,05 / %d rachas", nrow(r4))),
        "5b_rachas_ceros_prueba")

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
#   19 meses con un solo día de lluvia > p99; en varios solo hay un día observado
#   (se anotó el total mensual) → R3. De 917 rachas de ceros ≥30 días, 30 (en
#   27 estaciones) no tienen respaldo: su vecina más cercana llovió normal y la
#   probabilidad de que la estación marcara 0 todos esos días es < 0,05/917 → R4.
#   Ejemplo: Bella Vista 1982, 49 días en 0 mientras La Ceja (8,7 km) llovió 23
#   días (p = 7e-10). La mayoría de las 917 son sequías reales: la vecina también
#   estuvo seca (El Niño de 1997, temporadas secas). Loboguerrero oct–nov 1995,
#   antes revisada a mano, no pasa el umbral (p = 2e-4): Dagua coincide con ella
#   solo la mitad de los días de lluvia, y 12 días no alcanzan para descartar
#   una sequía real en un enclave seco. Se mantiene como dato.
#   No hay años copiados; lo más repetido son 11 días seguidos de 1 mm (llovizna).


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ┃ RÚBRICA 4.1 — Preprocesamiento: coherencia con las vecinas
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
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
#   La lluvia anual va de ~700 a ~12 000 mm. La vertiente Pacífica mezcla la costa
#   muy húmeda con el cañón seco del Dagua. La correlación con las vecinas cae con
#   la distancia: por eso los contrastes usan vecinas por distancia, no por vertiente.
rm(L, gz, M, W, C, N, D, pares, huella); invisible(gc())


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ┃ RÚBRICA 4.1 — Preprocesamiento: criterio de cobertura justificado (R5, R6)
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
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

# Reglas candidatas: cobertura ≥ X % y hueco máximo ≤ G días
tr <- merge(info[anio <= 1995], inc[, .(clave, hueco_max)], by = "clave", all.x = TRUE)
tr[is.na(hueco_max), hueco_max := fifelse(n_obs == 365, 0L, 999L)]
tr[, pct := 100 * n_obs / 365]
reglas_cand <- CJ(X = c(70, 80, 85, 90, 95), G = c(7, 15, 20, 30, 45, 60, 365))[, {
  acept <- tr[pct >= X & hueco_max <= G]
  e <- c8[pct >= X & hueco_max <= G, rel]
  .(estacion_anios = nrow(acept), estaciones_7anios = acept[, .N, by = codigo][N >= 7, .N],
    p90_err = round(quantile(e, .9), 2), p99_err = round(quantile(e, .99), 2))
}, by = .(X, G)]
print(tabla(reglas_cand, "7_reglas_candidatas")[X == 90])
grafico(ggplot(reglas_cand, aes(estacion_anios, p99_err, colour = factor(X), label = ifelse(G == 365, "sin límite", G))) +
          geom_point() + geom_line() + geom_text(nudge_y = .03, size = 3, show.legend = FALSE) +
          geom_hline(yintercept = .25, colour = "red", linetype = 2) + scale_y_log10() +
          labs(x = "Estación-años de entrenamiento que se conservan", y = "p99 del error relativo (log)",
               colour = "Cobertura ≥ %", title = "Datos que se conservan vs. error que se admite",
               subtitle = "Etiqueta = hueco máximo permitido (días)"),
        "7_reglas_candidatas")

# QUÉ DICEN LOS DATOS
#   Un hueco de 15 días mueve la curva 0,12; de 30 días, 0,23; de 60 días, 0,48
#   (K = 8). Con K = 12 los huecos largos pesan más. Los días sueltos casi no
#   importan: lo que daña la curva es el HUECO MÁS LARGO (R² 0,58 frente a 0,46
#   del % de cobertura). Con ≥90 % y hueco ≤15 días el p99 del error es 0,64;
#   con ≤30 días sube a 0,95 y solo se ganan 2 estaciones → regla R5.
rm(Ytodo, Ycomp, mascaras, ref, expA, expB, expC, c8, tr); invisible(gc())


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ┃ RÚBRICA 4.1 — Preprocesamiento: depuración y trazabilidad (estaciones y pares en cada paso)
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
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

# R4. Racha de ≥30 días en 0 que la vecina más cercana no respalda → faltante
#     (prueba de la sección 5b, umbral 0,05 / número de rachas).
m_cero <- r4[falso == TRUE][, d[codigo == .BY$codigo & fecha %between% c(ini, fin), .(fecha, P)], by = .(codigo, ini, fin, nombre_vecina, dias_lluvia_vecina, p_azar)][
  , .(codigo, fecha, P, regla = "R4 racha de ceros sin respaldo",
      detalle = sprintf("%s: llovió %d días; p = %.1e", nombre_vecina, dias_lluvia_vecina, p_azar))]
stopifnot(all(m_cero$P == 0), nrow(m_cero) == r4[falso == TRUE, sum(largo)])

marcas <- valle[, .(codigo, nombre)][rbind(m_r0, m_ext, m_mes, m_cero), on = "codigo"]
print(marcas[, .(dias = .N, estaciones = uniqueN(codigo)), by = regla])

# Aplicar: un faltante es una fila ausente / NA en la rejilla (R7: no se imputa)
limpio <- d[!marcas, on = .(codigo, fecha), .(codigo, anio, dia, fecha, P)]
anios_op <- unique(limpio[, .(codigo, anio)])
diario_todo <- limpio[anios_op[, .(dia = 1:365), by = .(codigo, anio)], on = .(codigo, anio, dia)]
diario_todo[is.na(fecha), fecha := as.IDate(sprintf("%d-01-01", anio)) + dia - 1L + (anio %% 4 == 0 & dia > 59)]
diario_todo[, particion := PARTICION(anio)]
setkey(diario_todo, codigo, anio, dia)
stopifnot(!anyDuplicated(diario_todo[, .(codigo, anio, dia)]), diario_todo[, all(month(fecha) * 100 + mday(fecha) != 229)])

# R5. Año válido: ≥90 % de los 365 días y ningún hueco > 15 días
anios <- diario_todo[, .(n_obs = sum(!is.na(P)), hueco_max = hueco_mas_largo(!is.na(P))), by = .(codigo, anio, particion)]
anios[, valido := n_obs >= COB_MIN * 365 & hueco_max <= HUECO_MAX]

# R6. Estación incluida: ≥7 años válidos en entrenamiento
por_est <- anios[, .(validos_ent = sum(valido & particion == "entrenamiento"),
                     validos_val = sum(valido & particion == "validación"),
                     validos_pru = sum(valido & particion == "prueba")), by = codigo]
estaciones <- merge(valle, por_est, by = "codigo")[validos_ent >= ANIOS_MIN]
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
  cbind(paso = "5. R4: rachas de ceros sin respaldo de la vecina a faltante", cuenta(d_r4)),
  cbind(paso = "6. R5: años con ≥90 % de días y hueco ≤15", cuenta(diario_todo[!is.na(P)][anios[valido == TRUE, .(codigo, anio)], on = .(codigo, anio), nomatch = 0])),
  cbind(paso = "7. R6: estaciones con ≥7 años válidos en entrenamiento", cuenta(diario[!is.na(P)])))
trazabilidad[, `:=`(estaciones_que_salen = c(NA, -diff(estaciones)), pares_que_salen = c(NA, -diff(pares_estacion_anio)),
                    dias_que_salen = c(NA, -diff(dias_con_dato)))]
print(tabla(trazabilidad, "8_trazabilidad"))
rm(d_r2, d_r3, d_r4)

reglas <- c(
  R0 = "Valor imposible → faltante: P < 0 o P > 1 825 mm (récord mundial OMM en 24 h)",
  R1 = "Dominio 1–365: se quita el 29 de febrero",
  R2 = "Extremo > 3·p99 y ≥150 mm sin vecina (≤25 km, ±1 día) que supere P/3 → faltante",
  R3 = "Mes con un solo día de lluvia > p99 y (mes con <50 % de días o vecinas con lluvia ≥40 % de los días) → mes faltante",
  R4 = "Racha ≥30 días en 0 → faltante si la vecina más cercana llovió k días y (1 − q)^k < 0,05/rachas (q de entrenamiento)",
  R5 = "Año válido: ≥90 % de 365 días y hueco máximo ≤15 días",
  R6 = "Estación: ≥7 años válidos en entrenamiento (1981–1995)",
  R7 = "Sin imputación")

saveRDS(list(diario = diario, diario_todo = diario_todo, anios = anios, estaciones = estaciones,
             marcas = marcas, meses_sospechosos = meses_sosp, reglas = reglas),
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
                  "no válida: menos de 90 % de días o hueco > 15", "sin datos ese año")
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


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ┃ RÚBRICA 4.1 / 3 — Satélite: lectura y fuente (insumo de la comparación por fuente)
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
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
#   Malla de 0,05° (unos 5,5 km), 39 × 37 píxeles, enteros de 0 a 200 mm por
#   péntada con tope en 200, sin faltantes, casi sin ceros y con datos sobre el mar.
#   Se parece a CHIRP (solo infrarrojo) y no a CHIRPS (con estaciones): los
#   archivos son del 6-may-2013, de una corrida de CHIRP anterior a la publicada.
rm(chirp, chirps, r_sat); invisible(gc())


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ┃ RÚBRICA 4.1 / 3 — Satélite: recorte al Valle y unidades
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
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


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ┃ RÚBRICA 2 / 3 — EDA del satélite
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
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


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ┃ RÚBRICA 3 — Satélite frente a estaciones (base del sesgo)
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
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
est_d <- diario[!is.na(P), .(codigo, anio, mes = month(fecha), dm = mday(fecha), P)]
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
ej_tope <- tope_est$codigo[1]
serie_tope <- cmp[codigo == ej_tope][order(anio, pentada)][, t := anio + (dia_centro - 1) / 365]
grafico(ggplot(melt(serie_tope[, .(t, anio, `Estación` = mm_est, `Satélite (su píxel)` = mm_sat)], id.vars = c("t", "anio")),
               aes(t, value, colour = variable, group = interaction(variable, anio))) + geom_line(linewidth = .3) +
          geom_hline(yintercept = 200, linetype = 2) +
          geom_point(data = serie_tope[mm_sat == 200], aes(t, mm_sat), inherit.aes = FALSE, colour = "red", size = 1) +
          labs(x = NULL, y = "mm por péntada", colour = NULL,
               title = paste0("Serie péntada a péntada de ", tope_est$nombre[1], " (", ETQ_VERT[tope_est$vertiente[1]], "), 1981–2003"),
               subtitle = "Línea punteada: tope de 200 mm del satélite; puntos rojos: péntadas en el tope"),
        "12_serie_tope_200", ancho = 11, alto = 4.5)

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

# Por año, como en la tesis del profesor (cap. 5: figuras 5.8 y 5.12, y la correlación
# estación–satélite de cada año). Solo estación-años con las 72 péntadas completas, para que
# el acumulado de la estación y el de su píxel sumen exactamente las mismas péntadas.
completos <- cmp[, .N, by = .(codigo, anio)][N == 72, .(codigo, anio)]
acu <- cmp[completos, on = .(codigo, anio)][order(codigo, anio, pentada)]
acu[, `:=`(acu_est = cumsum(mm_est), acu_sat = cumsum(mm_sat)), by = .(codigo, anio)]
acu_anio <- acu[pentada == 72, .(estaciones = .N, P_est = mean(acu_est), P_sat = mean(acu_sat)), by = anio]
acu_anio <- oni[acu_anio, on = "anio"][order(anio)]
acu_anio[, sat_sobre_est_pct := round(100 * P_sat / P_est)]
r_anio <- oni[cmp[, .(r = cor(mm_sat, mm_est), pares = .N), by = anio], on = "anio"][order(anio)]
acu_anio <- r_anio[, .(anio, r_pentadas = round(r, 2))][acu_anio, on = "anio"]
print(tabla(acu_anio[, .(anio, fase, ONI_anual, estaciones, P_est = round(P_est), P_sat = round(P_sat), sat_sobre_est_pct, r_pentadas)], "12_por_anio"))
cat(sprintf("Satélite / estación al final del año: %d a %d %% | r con el ONI: %.2f | r estación–satélite por año: %.2f a %.2f (mediana %.2f)\n",
            min(acu_anio$sat_sobre_est_pct), max(acu_anio$sat_sobre_est_pct), cor(acu_anio$sat_sobre_est_pct, acu_anio$ONI_anual),
            min(r_anio$r), max(r_anio$r), median(r_anio$r)))
lim_acu <- c(0, quantile(acu[pentada == 72, acu_est], .95))   # las pocas estaciones de la costa (hasta ~10 000 mm) se salen por arriba
grafico_r("12_acumulada_por_anio", function() {
  par(mfrow = c(6, 4), mar = c(3, 3.4, 2, .6), mgp = c(1.9, .5, 0), oma = c(0, 0, 3, 0))
  for (a in ANIOS) {
    d <- acu[anio == a]; f <- acu_anio[anio == a]
    plot(NA, xlim = c(1, 72), ylim = lim_acu, xlab = "Tiempo (péntada)", ylab = "Precipitación acumulada (mm)",
         main = sprintf("%d · %s", a, f$fase), cex.axis = .75, cex.lab = .8)
    for (k in split(d, d$codigo)) lines(k$pentada, k$acu_est, col = "grey80", lwd = .6)
    m <- d[, .(e = mean(acu_est), s = mean(acu_sat)), by = pentada][order(pentada)]
    lines(m$pentada, m$e, lwd = 2); lines(m$pentada, m$s, lwd = 2, lty = 2)
    text(1, lim_acu[2] * .95, sprintf("satélite = %d %% de la estación", f$sat_sobre_est_pct), adj = 0, cex = .8)
  }
  plot.new(); legend("center", c("Cada estación", "Promedio de las estaciones", "Promedio del satélite\nen sus píxeles"),
                     col = c("grey70", "black", "black"), lty = c(1, 1, 2), lwd = c(1, 2, 2), bty = "n", cex = 1)
  mtext("Precipitación acumulada por año: estaciones y satélite en los mismos pares (1981–2003)", outer = TRUE, cex = 1.1, font = 2)
})
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
#   - 676 píxeles con centro dentro del Valle. Las 176 estaciones caen en 148
#     píxeles distintos; 23 píxeles tienen más de una estación.
#   - El satélite aplana el espacio: va de ~1 100 a ~6 200 mm/año y las estaciones
#     de ~800 a ~12 100. Sobreestima donde llueve poco y subestima donde llueve mucho.
#   - Péntada a péntada sigue poco a la estación (r ≈ 0,5; mediana por estación 0,42)
#     y casi nunca marca cero (las estaciones, el 20–25 % de las péntadas).
#   - Acierta el clima medio (dónde y cuándo llueve, r ≈ 0,8) pero no las anomalías
#     (r ≈ 0,3 mensual y 0,25 anual). Su variación de un año a otro es la mitad.
#   - La diferencia con las estaciones sigue al ENSO (r = 0,74 con el ONI): se queda
#     corto en La Niña (1984, 1988, 1999, 2000) y se pasa en El Niño (1992). Los
#     "años defectuosos" de D-16 son en realidad este efecto, no archivos dañados.
#   - Santa Teresa (2612810102) tiene r = 0 con su píxel y la coordenada "a revisar".

saveRDS(list(pentadal = sat, pixeles = pix, calendario = cal, pares = pares, pares_pentada = cmp,
             anual = an_cmp, por_estacion = por_est_sat, que_acierta = aciertos),
        file.path(SALIDA, "satelite_valle.rds"), compress = "xz")
rm(S_todo, S_par, est_d); invisible(gc())


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ┃ RÚBRICA 4.1 / 4.4 — Datos preparados y partición 1981-1995 / 1996-2000 / 2001-2003 sin fuga
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
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
id_est <- anios[, .(codigo, anio, n_obs, hueco_max)][id_est, on = .(codigo, anio)]
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
  estaciones = list(Y = Y, id = id_est, t = 1:365, unidades = "mm/día; NA = día faltante"),
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


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ┃ RÚBRICA 4.1 — 27 pruebas automáticas de la parte A
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
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
prueba("Cada estación tiene ≥7 años válidos en entrenamiento", id_est[particion == "entrenamiento", .N, by = codigo][, .N == uniqueN(id_est$codigo) && min(N) >= ANIOS_MIN],
       paste("mínimo:", id_est[particion == "entrenamiento", .N, by = codigo][, min(N)]))

# -- Estaciones: reglas recalculadas desde la matriz
n_obs_Y <- rowSums(!is.na(Y)); hueco_Y <- apply(!is.na(Y), 1, hueco_mas_largo)
prueba("Cada curva tiene ≥90 % de días", all(n_obs_Y >= COB_MIN * 365), paste("mínimo:", min(n_obs_Y)))
prueba("Ninguna curva tiene un hueco > 15 días", all(hueco_Y <= HUECO_MAX), paste("máximo:", max(hueco_Y)))
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
prueba("Cada valor final es igual al crudo (salvo los marcados)", crudo2[is.na(marcado), all(!is.na(P_final) & P_final == P_crudo)],
       paste(crudo2[is.na(marcado), .N], "valores comparados"))
prueba("Los días marcados quedaron como faltante", crudo2[marcado == TRUE, all(is.na(P_final))], paste(crudo2[marcado == TRUE, .N], "marcados dentro de años válidos"))
prueba("Todo faltante de la matriz es un día ausente en el crudo o marcado", sum(!is.na(Y)) == crudo2[is.na(marcado), .N],
       paste("observados:", sum(!is.na(Y))))
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
comp <- which(n_obs_Y == 365)[1:50]
ok_pen <- vapply(comp, function(i) {
  tot <- est_pen[codigo == id_est$codigo[i] & anio == id_est$anio[i], sum(mm_est)]
  abs(tot - sum(Y[i, ])) < 1e-6
}, TRUE)
prueba("En años completos, la suma por péntadas = la suma diaria", all(ok_pen), "50 curvas")

pruebas <- rbindlist(pruebas)
print(tabla(pruebas, "14_pruebas"), right = FALSE)
if (any(pruebas$resultado != "OK")) stop("Hay pruebas que fallan: revisar la tabla anterior")
cat("\nTodas las pruebas pasan (", nrow(pruebas), ").\n")


# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ┃ RÚBRICA — — resumen de la parte A
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
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

# █████████████████████████████████████████████████████████████████████████████
# FIN DE LA PARTE A — INICIO DE LA PARTE B (lluvia_valle_funcional.R)
# La Parte B lee del disco lo que dejó la Parte A (salidas_R_datos/datos_preparados.rds).
# Se vacía la memoria de R (no los archivos) para que la Parte B corra igual que
# cuando era un archivo aparte y ningún objeto de la Parte A se mezcle con los suyos.
# █████████████████████████████████████████████████████████████████████████████
rm(list = ls()); invisible(gc())

# █████ PARTE B — FUNCIONAL (antes lluvia_valle_funcional.R) █████

# =============================================================================
# LLUVIA EN EL VALLE DEL CAUCA 1981-2003 — PARTE FUNCIONAL
# Proyecto final de Analítica de Datos (UAO, 2026-2S)
# =============================================================================
# Este archivo contiene TODA la parte funcional en R, en orden:
#   1. Datos preparados y emparejamiento estación-píxel
#   2. Construcción de datos funcionales            (rúbrica 1)
#   3. Análisis exploratorio de datos funcionales   (rúbrica 2)
#   4. Comparación funcional por fuente             (rúbrica 3)
#   5. Regresión funcional función a función        (rúbrica 4.2 y 4.4, referencia de la red)
#
# Entrada única: salidas_R_datos/datos_preparados.rds, que produce
#   lluvia_valle_datos.R (estaciones depuradas, satélite en mm/día, pares,
#   altitud y ONI). Aquí no se vuelve a leer ni a limpiar ningún dato crudo:
#   así las 27 pruebas del archivo de datos cubren todo lo que entra.
# Salidas: salidas_R_funcional/tablas, /graficos y /para_python (insumos de la red).
#
# Cómo se corre (desde la raíz del repositorio):
#   Rscript lluvia_valle_datos.R        # una vez, deja datos_preparados.rds
#   Rscript lluvia_valle_funcional.R    # este archivo, ~1 minuto (fija UTF-8 solo)
#
# Reglas del proyecto que se cumplen aquí:
#   - La base B-spline (Cox-de Boor), los mínimos cuadrados, S, W, la media, la
#     covarianza, las derivadas, las profundidades y el FPCA están programados
#     paso a paso. Las librerías fda y splines solo se usan para VERIFICAR.
#   - Toda elección (K, M, lambda, escalas) sale de entrenamiento (1981-1995) o
#     de validación (1996-2000). La prueba (2001-2003) solo se usa para evaluar.
#   - No se elimina ningún año (D-47, asesoría del 13-nov-2026).
#   - Cada decisión tiene su verificación numérica en el código: una tabla, una
#     gráfica o un stopifnot(). El bloque "QUÉ DICEN LOS DATOS" de cada sección
#     resume lo que salió.
# =============================================================================

# Locale: si R arranca sin UTF-8 (en la Mac del grupo LANG está vacío), las tildes del
# código ("validación") no coinciden con las de los datos y el script se detiene en la
# primera comparación de texto. Se fija UTF-8 aquí, antes de cualquier cadena con tilde.
if (!isTRUE(l10n_info()$`UTF-8`)) {
  for (loc in c("en_US.UTF-8", "es_CO.UTF-8", "es_ES.UTF-8", "C.UTF-8"))
    if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break
  if (!isTRUE(l10n_info()$`UTF-8`)) stop("No hay un locale UTF-8 disponible; corra con LC_ALL=en_US.UTF-8 Rscript ...")
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(fda)        # herramientas de clase: verificación (media, fbplot) y bases del contraste 5.5b
  library(fda.usc)    # herramientas de clase: fdata y fregre.basis.fr (contraste 5.5b)
})
theme_set(theme_minimal(base_size = 11))

ANIOS     <- 1981:2003
PARTICION <- function(anio) fcase(anio <= 1995, "entrenamiento", anio <= 2000, "validación", default = "prueba")
RANGO     <- c(1, 365)
ORDEN     <- 4L
T_DIA     <- 1:365
# Años en que el satélite se queda más corto (La Niña). NO se excluyen (D-47): solo se usan
# para el análisis de sensibilidad "con y sin esos años".
ANIOS_MARCADOS <- c(1985L, 1988L, 2000L)
COB_MIN_R5 <- 0.90
HUECO_MAX_R5 <- 15L

SALIDA_F <- "salidas_R_funcional"
dir.create(file.path(SALIDA_F, "tablas"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(SALIDA_F, "graficos"), showWarnings = FALSE, recursive = TRUE)
tabla_f   <- function(x, nombre) { fwrite(x, file.path(SALIDA_F, "tablas", paste0(nombre, ".csv"))); invisible(x) }
grafico_f <- function(p, nombre, ancho = 9, alto = 5) {
  ggsave(file.path(SALIDA_F, "graficos", paste0(nombre, ".png")), p, width = ancho, height = alto, dpi = 110, bg = "white")
  if (interactive()) print(p)
  invisible(p)
}
grafico_base <- function(nombre, dibujar, ancho = 11, alto = 5) {
  png(file.path(SALIDA_F, "graficos", paste0(nombre, ".png")), width = ancho, height = alto, units = "in", res = 110)
  dibujar(); dev.off()
  if (interactive()) dibujar()
  invisible(NULL)
}
titulo <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
fecha_dia <- function(d) format(as.Date(round(d) - 1, origin = "2001-01-01"), "%d-%b")

nodos <- function(K, orden = ORDEN) {
  interiores <- seq(RANGO[1], RANGO[2], length.out = K - orden + 2L)
  c(rep(RANGO[1], orden), interiores[-c(1L, length(interiores))], rep(RANGO[2], orden))
}
base_en_nodos <- function(t, tau, orden) {
  B <- outer(t, seq_len(length(tau) - 1L), function(x, j) (x >= tau[j]) & (x < tau[j + 1L])) * 1
  B[t == max(tau), max(which(diff(tau) > 0))] <- 1
  if (orden > 1L) for (m in 2:orden) {
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
base_bspline <- function(t, K, orden = ORDEN) base_en_nodos(t, nodos(K, orden), orden)
deriv_base <- function(t, K) {
  tau <- nodos(K); B <- base_en_nodos(t, tau, ORDEN - 1L); D <- matrix(0, length(t), K)
  for (j in seq_len(K)) {
    a <- tau[j + ORDEN - 1L] - tau[j]; b <- tau[j + ORDEN] - tau[j + 1L]
    if (a > 0) D[, j] <- D[, j] + (ORDEN - 1L) * B[, j] / a
    if (b > 0) D[, j] <- D[, j] - (ORDEN - 1L) * B[, j + 1L] / b
  }
  D
}

derivada_base <- function(t, K, n = 1L, tau = nodos(K), m = ORDEN) {
  if (n == 0L) return(base_en_nodos(t, tau, m))
  B <- derivada_base(t, K, n - 1L, tau, m - 1L); D <- matrix(0, length(t), length(tau) - m)
  for (j in seq_len(ncol(D))) {
    a <- tau[j + m - 1L] - tau[j]; b <- tau[j + m] - tau[j + 1L]
    if (a > 0) D[, j] <- D[, j] + (m - 1L) * B[, j] / a
    if (b > 0) D[, j] <- D[, j] - (m - 1L) * B[, j + 1L] / b
  }
  D
}
gram_derivada <- function(K, n = 2L) {
  tau <- unique(nodos(K)); a <- head(tau, -1); b <- tail(tau, -1)
  tq <- as.vector(outer(GL_x, seq_along(a), function(x, k) (b[k] - a[k]) / 2 * x + (a[k] + b[k]) / 2))
  wq <- as.vector(outer(GL_w, seq_along(a), function(w, k) (b[k] - a[k]) / 2 * w))
  P <- derivada_base(tq, K, n); crossprod(P, wq * P)
}

ajuste_mco <- function(y, Phi) {
  obs <- !is.na(y); P <- Phi[obs, , drop = FALSE]
  if (nrow(P) <= ncol(P) || qr(P)$rank != ncol(P)) stop("Curva no estimable con este K")
  mapa <- solve(crossprod(P), t(P))
  coef <- drop(mapa %*% y[obs]); H <- P %*% mapa
  res <- y[obs] - drop(P %*% coef)
  list(coef = coef, H = H, mapa = mapa, obs = obs, sigma2 = sum(res^2) / (nrow(P) - ncol(P)))
}
control_S <- function(H, K) c(simetria = max(abs(H - t(H))), idempotencia = max(abs(H %*% H - H)), traza_menos_K = abs(sum(diag(H)) - K))
# Escala raíz de una curva en mm: se evalúa, se recorta en 0, se saca la raíz y se reproyecta en la misma base.
a_raiz <- function(C, Phi) {
  F <- C %*% t(Phi); G <- sqrt(pmax(F, 0)); Cr <- G %*% t(solve(crossprod(Phi), t(Phi)))
  list(coef = Cr, pct_recortado = 100 * mean(F < 0), rmse_reproyeccion = sqrt(mean((Cr %*% t(Phi) - G)^2)))
}

criterios_K <- function(Y, t, Ks) {
  G <- crossprod(Y); n <- ncol(Y); N <- nrow(Y)
  rbindlist(lapply(Ks, function(K) {
    P <- base_bspline(t, K); mapa <- solve(crossprod(P), t(P)); h <- rowSums(P * t(mapa))
    T0 <- G %*% t(mapa)
    sse <- pmax(0, diag(G) - 2 * rowSums(P * T0) + rowSums((P %*% (mapa %*% T0)) * P))
    data.table(K = K, ECM = sum(sse) / (N * n), VC = sum(sse / (1 - h)^2) / (N * n), VCG = n * sum(sse) / (n - K)^2 / N)
  }))
}

construir <- function(Y, t, K) {
  Phi <- base_bspline(t, K); completas <- rowSums(is.na(Y)) == 0L
  C <- matrix(NA_real_, nrow(Y), K); sigma2 <- rep(NA_real_, nrow(Y))
  mapa <- solve(crossprod(Phi), t(Phi)); H <- Phi %*% mapa
  C[completas, ] <- Y[completas, , drop = FALSE] %*% t(mapa)
  Res <- Y[completas, , drop = FALSE] - C[completas, , drop = FALSE] %*% t(Phi)
  sigma2[completas] <- rowSums(Res^2) / (length(t) - K)
  controles <- control_S(H, K)
  for (i in which(!completas)) {
    a <- ajuste_mco(Y[i, ], Phi); C[i, ] <- a$coef; sigma2[i] <- a$sigma2
    controles <- pmax(controles, control_S(a$H, K))
  }
  colnames(C) <- paste0("c", seq_len(K))
  list(coef = C, K = K, Phi = Phi, mapa = mapa, H = H, sigma2 = sigma2, Res = Res, completas = completas,
       controles = controles, n_obs = rowSums(!is.na(Y)))
}

GL_x <- c(-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053)
GL_w <- c(0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454)
gram <- function(K) {
  tau <- unique(nodos(K)); a <- head(tau, -1); b <- tail(tau, -1)
  tq <- as.vector(outer(GL_x, seq_along(a), function(x, k) (b[k] - a[k]) / 2 * x + (a[k] + b[k]) / 2))
  wq <- as.vector(outer(GL_w, seq_along(a), function(w, k) (b[k] - a[k]) / 2 * w))
  P <- base_bspline(tq, K); crossprod(P, wq * P)
}

profundidades <- function(Y) {
  N <- nrow(Y); fm <- mbd <- numeric(N)
  for (j in seq_len(ncol(Y))) {
    r <- rank(Y[, j], ties.method = "max"); L <- rank(Y[, j], ties.method = "min") - 1; U <- N - r
    fm  <- fm + 1 - abs(.5 - r / N)
    mbd <- mbd + 1 - (L * (L - 1) + U * (U - 1)) / (N * (N - 1))
  }
  list(FM = fm / ncol(Y), MBD = mbd / ncol(Y))
}
boxplot_funcional <- function(Y, prof, factor = 1.5) {
  centro <- order(-prof, seq_along(prof))[seq_len(ceiling(nrow(Y) / 2))]
  lo <- apply(Y[centro, , drop = FALSE], 2, min); hi <- apply(Y[centro, , drop = FALSE], 2, max)
  inf <- lo - factor * (hi - lo); sup <- hi + factor * (hi - lo)
  atip <- apply(Y, 1, function(y) any(y < inf | y > sup))
  list(lo = lo, hi = hi, inf = inf, sup = sup, atip = atip, mediana = Y[which.max(prof), ],
       bigote_inf = apply(Y[!atip, , drop = FALSE], 2, min), bigote_sup = apply(Y[!atip, , drop = FALSE], 2, max))
}
dibujar_boxplot <- function(Y, bx, titulo, t = T_DIA) {
  plot(t, bx$mediana, type = "n", ylim = range(c(0, bx$bigote_inf, bx$bigote_sup, Y[bx$atip, ])), xlab = "Día del año", ylab = "mm/día", main = titulo)
  if (any(bx$atip)) matlines(t, t(Y[bx$atip, , drop = FALSE]), col = "#cc333330", lty = 1)
  polygon(c(t, rev(t)), c(bx$lo, rev(bx$hi)), col = "#b39ddbcc", border = NA)
  lines(t, bx$bigote_inf, col = 4); lines(t, bx$bigote_sup, col = 4)
  lines(t, bx$mediana, lwd = 2)
  legend("topleft", c("Región central 50 %", "Mediana (MBD)", "Bigotes", "Atípicas"),
         col = c("#b39ddb", 1, 4, "#cc3333"), lty = c(NA, 1, 1, 1), pch = c(15, NA, NA, NA), cex = .75, bty = "n")
}
mapa_calor <- function(t, z, titulo, limites = range(z, na.rm = TRUE), paleta = "Blue-Red 3", xlab = "Día s", ylab = "Día t") {
  colores <- hcl.colors(64, paleta)
  layout(matrix(1:2, 1), widths = c(5, 1)); par(mar = c(4, 4, 3, 1))
  image(t, t, z, zlim = limites, col = colores, main = titulo, xlab = xlab, ylab = ylab)
  par(mar = c(4, 0, 3, 4)); escala <- seq(limites[1], limites[2], length.out = 65)
  image(c(0, 1), escala, matrix(rep(escala, each = 2), nrow = 2), col = colores, zlim = limites, axes = FALSE, xlab = "", ylab = "")
  axis(4, las = 1, cex.axis = .8); box(); layout(1)
}

# =============================================================================
# 1. DATOS PREPARADOS Y EMPAREJAMIENTO ESTACIÓN-PÍXEL  — rúbrica 3 (emparejamiento)
# =============================================================================
# Lee datos_preparados.rds y deja las matrices con las que trabaja todo el
# archivo: Y (estación-año x 365 días, mm/día, NA = día faltante) y X (píxel-año
# x 72 péntadas, mm/día, cada péntada en su día central). Documenta el
# emparejamiento (regla, distancias y píxeles compartidos) y verifica con datos
# dos decisiones: la vertiente de los píxeles y el trato de 1985, 1988 y 2000.
# =============================================================================
# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ┃ RÚBRICA 3 — Emparejamiento estación–píxel documentado
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
titulo("1. Datos preparados y emparejamiento estación-píxel")

dp  <- readRDS("salidas_R_datos/datos_preparados.rds")
oni <- fread("datos_web/oni_anual_fase.csv")
Y <- dp$estaciones$Y; X <- dp$satelite$X; T_SAT <- dp$satelite$t
storage.mode(Y) <- "double"
Y_dia <- Y; X_pen <- X      # copias de los datos discretos (las secciones 4 y 5 reutilizan los nombres Y y X para curvas)
stopifnot(ncol(Y) == 365, ncol(X) == 72, !anyNA(X), all(Y >= 0, na.rm = TRUE), all(X >= 0),
          identical(dp$estaciones$id$fila, seq_len(nrow(Y))), identical(dp$satelite$id$fila, seq_len(nrow(X))),
          all(PARTICION(dp$estaciones$id$anio) == dp$estaciones$id$particion))

est <- copy(dp$meta_estaciones)
id_est <- dp$estaciones$id[, .(codigo, anio, particion, n_obs, hueco_max)]
id_est[, curva_id := paste(codigo, anio, sep = "_")]

# --- Emparejamiento: cada estación toma la celda de la malla que la contiene ---
px <- dp$pixeles[, .(celda = pixel, lon, lat, en_valle, con_estacion)]
emparejamiento <- px[, .(celda, lon_px = lon, lat_px = lat, px_en_valle = en_valle)][
  est[, .(codigo, nombre, vertiente, altitud, lat, lon, celda = pixel, dist_km = km_al_centro)], on = "celda"]
emparejamiento[, estaciones_en_el_pixel := .N, by = celda]
stopifnot(!anyNA(emparejamiento$lon_px), abs(emparejamiento$lon - emparejamiento$lon_px) <= 0.025 + 1e-6,
          abs(emparejamiento$lat - emparejamiento$lat_px) <= 0.025 + 1e-6)   # la estación cae DENTRO de su celda de 0,05°
tabla_f(emparejamiento[order(celda, codigo)], "F1_emparejamiento_estacion_pixel")
compartidos <- emparejamiento[estaciones_en_el_pixel > 1, .(estaciones = .N, codigos = paste(codigo, collapse = " / "),
                                                            nombres = paste(nombre, collapse = " / ")), by = celda]
tabla_f(compartidos, "F1_pixeles_compartidos")
cat("Estaciones:", nrow(emparejamiento), "| píxeles distintos:", uniqueN(emparejamiento$celda),
    "| píxeles con más de una estación:", nrow(compartidos), "(", sum(compartidos$estaciones), "estaciones ) |",
    "distancia al centro del píxel: mediana", round(median(emparejamiento$dist_km), 2), "km, máx", round(max(emparejamiento$dist_km), 2), "km |",
    "píxel con centro fuera del polígono del Valle:", emparejamiento[px_en_valle == FALSE, .N], "\n")
print(compartidos[order(-estaciones)][1:5])

# --- Vertiente de cada píxel: la de la estación más cercana (regla verificada) ---
# La vertiente de las estaciones viene del código CVC (cuenca hidrográfica). Los
# píxeles no traen ese dato. Se comparan dos reglas dejando una estación fuera
# cada vez y prediciendo su vertiente: (a) un corte en la longitud; (b) la
# vertiente de la estación más cercana.
haversine <- function(la1, lo1, la2, lo2) 12742 * asin(sqrt(sin((la2 - la1) * pi / 360)^2 + cos(la1 * pi / 180) * cos(la2 * pi / 180) * sin((lo2 - lo1) * pi / 360)^2))
D_ee <- outer(seq_len(nrow(est)), seq_len(nrow(est)), function(i, j) haversine(est$lat[i], est$lon[i], est$lat[j], est$lon[j])); diag(D_ee) <- Inf
cortes <- seq(-77.1, -76.4, by = .05)
acierto_corte <- sapply(cortes, function(k) mean((est$lon < k) == (est$vertiente == "Pacífico")))
acierto_vecina <- mean(est$vertiente[apply(D_ee, 1, which.min)] == est$vertiente)
regla_vertiente <- data.table(regla = c(paste0("corte en longitud ", cortes[which.max(acierto_corte)], "° (el mejor corte)"), "corte en longitud -76.8° (el que se usaba)", "estación más cercana (dejando una fuera)"),
                              pct_acierto = round(100 * c(max(acierto_corte), acierto_corte[abs(cortes + 76.8) < 1e-9], acierto_vecina), 1))
print(tabla_f(regla_vertiente, "F1_regla_vertiente_pixeles"))
D_pe <- outer(seq_len(nrow(px)), seq_len(nrow(est)), function(i, j) haversine(px$lat[i], px$lon[i], est$lat[j], est$lon[j]))
px[, `:=`(vertiente = est$vertiente[apply(D_pe, 1, which.min)], km_estacion_cercana = round(apply(D_pe, 1, min), 1))]
stopifnot(acierto_vecina > max(acierto_corte))
cat("Píxeles por vertiente (estación más cercana):\n"); print(px[, .(pixeles = .N, km_mediana_a_la_estacion = median(km_estacion_cercana)), by = vertiente])

id_sat <- px[dp$satelite$id[, .(celda = pixel, anio, particion)], on = "celda"]
stopifnot(identical(id_sat$celda, dp$satelite$id$pixel), dp$satelite$dias[12] == 3)
dias_reales <- matrix(dp$satelite$dias, nrow(X), 72, byrow = TRUE); dias_reales[id_sat$anio %% 4 == 0, 12] <- 4
id_sat[, `:=`(n_tope = rowSums(X * dias_reales >= 199.99), marcado = anio %in% ANIOS_MARCADOS, curva_id = paste(anio, celda, sep = "_"))]
cat("Curvas de estación:", nrow(Y), "(", uniqueN(id_est$codigo), "estaciones ) | curvas de satélite:", nrow(X), "=", uniqueN(id_sat$anio), "años x", uniqueN(id_sat$celda), "píxeles |",
    "con alguna péntada en el tope de 200 mm:", sum(id_sat$n_tope > 0), sprintf("(%.1f %%)", 100 * mean(id_sat$n_tope > 0)), "\n")
print(id_est[, .(curvas = .N, completas = sum(n_obs == 365), estaciones = uniqueN(codigo)), by = particion])

# --- 1985, 1988 y 2000: ¿archivos dañados o error que sigue al ENSO? ---
# Se compara, año por año, la lluvia anual de cada píxel con la de su estación
# (datos discretos, sin suavizar) y se cruza con el ONI.
anual <- dp$pares[, .(codigo, anio, particion, ONI, fase = fase_enso, P_est = rowMeans(Y[fila_est, ], na.rm = TRUE) * 365, P_sat = rowMeans(X[fila_sat, ]) * 365)]
por_anio <- anual[, .(pares = .N, P_est_mediana = round(median(P_est)), P_sat_mediana = round(median(P_sat)), cociente_mediano = round(median(P_sat / P_est), 3),
                      ONI = ONI[1], fase = fase[1], particion = particion[1]), by = anio][order(anio)]
por_anio[, marcado := anio %in% ANIOS_MARCADOS]
print(tabla_f(por_anio, "F1_satelite_vs_estaciones_por_anio"))
r_oni <- por_anio[, cor(ONI, cociente_mediano)]; r_oni_ent <- por_anio[particion == "entrenamiento", cor(ONI, cociente_mediano)]
cat(sprintf("Correlación entre el ONI del año y el cociente satélite/estación: r = %.2f con los 23 años (%.2f solo con entrenamiento)\n", r_oni, r_oni_ent))
cat("Cociente mediano por fase ENSO:\n"); print(por_anio[, .(anios = .N, cociente = round(mean(cociente_mediano), 3), anios_lista = paste(anio, collapse = " ")), by = fase])
cat("Años con el cociente más bajo:", paste(por_anio[order(cociente_mediano)][1:5, sprintf("%d (%.2f, %s)", anio, cociente_mediano, fase)], collapse = ", "), "\n")
grafico_f(ggplot(por_anio, aes(ONI, cociente_mediano)) + geom_hline(yintercept = 1, linetype = 3) +
            geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey50", linewidth = .5) +
            geom_text(aes(colour = fase, label = anio), size = 3.2, fontface = ifelse(por_anio$marcado, "bold", "plain")) +
            scale_colour_manual(values = c(Neutro = "grey30", `Niña` = "steelblue", `Niño` = "firebrick")) +
            labs(x = "ONI anual", y = "Mediana de satélite / estación (lluvia anual)", colour = NULL,
                 title = "El error del satélite sigue al ENSO", subtitle = sprintf("r = %.2f; en negrilla 1985, 1988 y 2000", r_oni)),
          "F1_cociente_satelite_estacion_vs_ONI", ancho = 8, alto = 5)

# QUÉ DICEN LOS DATOS
#   Datos: 3 299 curvas estación-año de 176 estaciones (2 257 de entrenamiento,
#   672 de validación y 370 de prueba; 2 812 sin ningún día faltante) y 15 663
#   curvas píxel-año (681 píxeles x 23 años). Todo viene de datos_preparados.rds.
#   Emparejamiento: cada estación toma la celda de 0,05° que la contiene (se
#   comprueba que la coordenada cae dentro). Las 176 estaciones caen en 148
#   píxeles; 23 píxeles tienen 2 a 4 estaciones (51 estaciones), que comparten la
#   misma curva satelital. El centro del píxel queda a 2,2 km en mediana y nunca
#   a más de 3,9 km. Cinco píxeles tienen el centro fuera del polígono del Valle,
#   pero su estación está dentro.
#   Vertiente de los píxeles: antes se usaba un corte en la longitud -76,8°. Al
#   comprobarlo con las estaciones, ese corte solo acierta el 88,6 % y deja como
#   "Cauca" a 20 de las 26 estaciones que la CVC pone en la vertiente Pacífica
#   (cuencas del Dagua, Calima y Garrapatas, sobre la cordillera Occidental).
#   Ningún corte en longitud pasa de 89,2 %. La regla "vertiente de la estación
#   más cercana" acierta el 97,7 % dejando una estación fuera cada vez, así que
#   se usa esa: una sola definición (la cuenca de la CVC) para las dos fuentes.
#   1985, 1988 y 2000: no son archivos dañados. El cociente satélite/estación de
#   cada año sigue al ONI (r = 0,68): en años Niña el satélite se queda corto
#   (cociente medio 0,99) y en años Niño se pasa (1,29); en neutros da 1,16. Los
#   cinco cocientes más bajos son de 1988, 1984, 2000, 1999 y 1996: cuatro Niña.
#   1984 y 1999 se apartan tanto como los tres años marcados y nunca se marcaron.
#   Es un error sistemático del satélite ligado al ENSO, justo lo que el modelo
#   debe aprender con el ONI. Por eso no se excluye ningún año (D-47).


# =============================================================================
# 2. CONSTRUCCIÓN DE DATOS FUNCIONALES (ESTACIONES Y SATÉLITE)  — rúbrica 1
# =============================================================================
# Base B-spline cúbica (Cox-de Boor) y mínimos cuadrados programados sin
# librerías; K por ECM, VC y VCG solo con entrenamiento (1981-1995, todos los
# años); verificación de S; covarianza residual; bandas punto a punto; matriz
# de Gram W. Ambas fuentes quedan en el mismo dominio 1-365.
# Días faltantes: cada curva se ajusta solo con sus días observados (no se imputa).
# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ┃ RÚBRICA 1 — Construcción de datos funcionales (paso a paso, sin librerías)
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
titulo("2. Construcción de datos funcionales")

completas_est <- rowSums(is.na(Y)) == 0
stopifnot(all(completas_est == (id_est$n_obs == 365L)))

grafico_base("F2_datos_discretos", function() {
  par(mfrow = c(1, 2))
  i <- which(completas_est & id_est$codigo == "2622330101")[1]
  plot(T_DIA, Y[i, ], type = "h", col = "grey40", xlab = "Día del año", ylab = "mm/día",
       main = paste("Estación", id_est$codigo[i], id_est$anio[i], "(diario)"))
  j <- which(id_sat$celda == emparejamiento[codigo == "2622330101", celda] & id_sat$anio == id_est$anio[i])
  plot(T_SAT, X[j, ], type = "b", pch = 16, cex = .6, col = "darkorange", xlab = "Día del año", ylab = "mm/día",
       main = paste("Su píxel, mismo año (péntada en su día central)"))
})

for (K in c(6, 8, 11, 20, 40)) {
  P <- base_bspline(T_DIA, K)
  stopifnot(all(is.finite(P)), min(P) >= -1e-12, max(abs(rowSums(P) - 1)) < 1e-10)
}
ref_splines <- splines::splineDesign(nodos(11), T_DIA, ord = ORDEN, outer.ok = TRUE)
cat("Base propia vs splines::splineDesign (K = 11), diferencia máxima:", format(max(abs(base_bspline(T_DIA, 11) - ref_splines)), digits = 2), "\n")

i_ej <- which(completas_est)[1]; y1 <- Y[i_ej, ]
P8 <- base_bspline(T_DIA, 8); a8 <- ajuste_mco(y1, P8)
vc_atajo <- mean(((y1 - a8$H %*% y1) / (1 - diag(a8$H)))^2)
vc_bruta <- mean(sapply(1:365, function(j) (y1[j] - P8[j, ] %*% solve(crossprod(P8[-j, ]), crossprod(P8[-j, ], y1[-j])))^2))
cat(sprintf("VC por el atajo de diag(S) = %.8f | VC dejando un día fuera (365 ajustes) = %.8f\n", vc_atajo, vc_bruta))

ent_est <- completas_est & id_est$particion == "entrenamiento"
ent_sat <- id_sat$particion == "entrenamiento"
sel_est <- criterios_K(Y[ent_est, ], T_DIA, c(4:40, seq(45, 100, 5)))
sel_sat <- criterios_K(X[ent_sat, ], T_SAT, 4:40)
sel <- rbind(sel_est[, fuente := "estaciones"], sel_sat[, fuente := "satélite"])
tabla_f(sel, "F2_seleccion_K")
resumen_K <- sel[, .(K_VC = K[which.min(VC)], K_VCG = K[which.min(VCG)],
                     K_VC_1pct = paste(range(K[VC <= 1.01 * min(VC)]), collapse = "-")), by = fuente]
resumen_K[, curvas_entrenamiento := c(sum(ent_est), sum(ent_sat))]
print(tabla_f(resumen_K, "F2_K_elegido"))
K_EST <- resumen_K[fuente == "estaciones", K_VCG]; K_SAT <- resumen_K[fuente == "satélite", K_VCG]

# Sensibilidad de K en las estaciones. La VC es casi plana alrededor del mínimo, así que
# se comprueba qué cambia con una base más flexible: las fechas de los picos (métrica de
# la rúbrica 4.4) y los valores negativos de la curva. Solo con curvas de entrenamiento.
sens_K <- rbindlist(lapply(sort(unique(c(K_EST, 12L, 16L, 24L))), function(K) {
  P <- base_bspline(T_DIA, K); Fk <- Y[ent_est, ] %*% t(solve(crossprod(P), t(P))) %*% t(P)
  p1 <- 59 + apply(Fk[, 60:200], 1, which.max); p2 <- 229 + apply(Fk[, 230:350], 1, which.max)
  data.table(K = K, VC_sobre_el_minimo_pct = round(100 * (sel_est$VC[sel_est$K == K] / min(sel_est$VC) - 1), 2),
             pico_1_mediana = fecha_dia(median(p1)), pico_1_RIC_dias = IQR(p1), pico_2_mediana = fecha_dia(median(p2)), pico_2_RIC_dias = IQR(p2),
             pct_picos_en_el_borde_de_la_ventana = round(100 * mean(p1 %in% c(60, 200) | p2 %in% c(230, 350)), 1),
             pct_puntos_negativos = round(100 * mean(Fk < 0), 2))
}))
print(tabla_f(sens_K, "F2_sensibilidad_K_estaciones"))

# P-20 (D-51): con K mayor la fecha del pico se acerca a una media móvil de 31 días de
# los datos crudos (prueba del stash "fase A"). ¿Es una mejora real o depende de la
# referencia? Se mide contra tres medias móviles (15, 31 y 61 días; ninguna depende de K)
# en las curvas completas de la vertiente interandina (régimen bimodal), solo de
# entrenamiento. Si el K que gana cambia con la referencia, la fecha del pico no sirve
# para elegir K y se queda el de la VCG.
vert_curva <- est$vertiente[match(id_est$codigo, est$codigo)]
cur_p20 <- ent_est & vert_curva != "Pacífico"
media_movil <- function(M, w) t(apply(M, 1, function(y) as.numeric(stats::filter(y, rep(1 / w, w), sides = 2))))
pico_en <- function(M, v) apply(M[, v, drop = FALSE], 1, which.max) + min(v) - 1
V1 <- 60:200; V2 <- 230:350; ANCHOS_REF <- c(15L, 31L, 61L)
picos_ref <- lapply(ANCHOS_REF, function(w) { R <- media_movil(Y[cur_p20, ], w); cbind(pico_en(R, V1), pico_en(R, V2)) })
K_picos <- rbindlist(lapply(c(K_EST, 10L, 12L, 16L, 20L), function(k) {
  P <- base_bspline(T_DIA, k); Fk <- Y[cur_p20, ] %*% t(P %*% solve(crossprod(P), t(P)))
  pk <- cbind(pico_en(Fk, V1), pico_en(Fk, V2))
  rbindlist(lapply(seq_along(ANCHOS_REF), function(r) data.table(K = k, referencia_dias = ANCHOS_REF[r],
    err_pico_1 = median(abs(pk[, 1] - picos_ref[[r]][, 1])), err_pico_2 = median(abs(pk[, 2] - picos_ref[[r]][, 2])))))
}))
K_picos[, err_total := err_pico_1 + err_pico_2]
gana_por_ref <- K_picos[, .(K_que_gana = K[which.min(err_total)]), by = referencia_dias]
dif_refs <- rbindlist(lapply(list(c(1, 2), c(1, 3), c(2, 3)), function(ab) data.table(
  referencias = paste(ANCHOS_REF[ab], collapse = " vs "),
  pico_1 = median(abs(picos_ref[[ab[1]]][, 1] - picos_ref[[ab[2]]][, 1])), pico_2 = median(abs(picos_ref[[ab[1]]][, 2] - picos_ref[[ab[2]]][, 2])))))
print(tabla_f(dcast(K_picos, K ~ referencia_dias, value.var = c("err_pico_1", "err_pico_2")), "F2_K_picos_vs_referencia"))
print(tabla_f(gana_por_ref, "F2_K_que_gana_por_referencia"))
print(tabla_f(dif_refs, "F2_diferencia_entre_referencias"))
cat(sprintf("Curvas: %d. Dispersión de la fecha del pico entre curvas (RIC, ref. 31 días): %d y %d días\n",
            sum(cur_p20), as.integer(IQR(picos_ref[[2]][, 1])), as.integer(IQR(picos_ref[[2]][, 2]))))
# La decisión depende de este resultado: si una misma K ganara con las tres referencias,
# el argumento de D-51 no valdría y el script se detiene para revisarlo.
stopifnot(uniqueN(gana_por_ref$K_que_gana) > 1)

sel_largo <- melt(sel, id.vars = c("fuente", "K"), variable.name = "criterio")
marcas_K <- sel_largo[criterio != "ECM", .SD[which.min(value)], by = .(fuente, criterio)]
grafico_f(ggplot(sel_largo[K <= 60], aes(K, value)) + geom_line(colour = "steelblue") + geom_point(size = .8, colour = "steelblue") +
            geom_vline(data = marcas_K, aes(xintercept = K), colour = "red", linetype = 2) +
            geom_text(data = marcas_K, aes(K, Inf, label = paste("K =", K)), vjust = 1.5, hjust = -.1, colour = "red", size = 3) +
            facet_wrap(fuente ~ criterio, scales = "free_y", ncol = 3) + scale_y_log10() +
            labs(x = "Número de bases K", y = "Criterio (log)", title = "Selección de K con curvas de entrenamiento"),
          "F2_seleccion_K", ancho = 11, alto = 6)

e <- construir(Y, T_DIA, K_EST)
s <- construir(X, T_SAT, K_SAT)
controles <- rbind(data.table(fuente = "estaciones", K = K_EST, t(e$controles)), data.table(fuente = "satélite", K = K_SAT, t(s$controles)))
print(tabla_f(controles, "F2_controles_S"))
stopifnot(all(controles[, .(simetria, idempotencia, traza_menos_K)] < 1e-8))

Phi_e <- base_bspline(T_DIA, K_EST); Phi_s365 <- base_bspline(T_DIA, K_SAT)
F_est <- e$coef %*% t(Phi_e); F_sat <- s$coef %*% t(Phi_s365)
negativos <- data.table(fuente = c("estaciones", "satélite"),
                        pct_puntos_negativos = round(100 * c(mean(F_est < 0), mean(F_sat < 0)), 2),
                        pct_curvas_con_negativos = round(100 * c(mean(apply(F_est < 0, 1, any)), mean(apply(F_sat < 0, 1, any))), 1),
                        minimo_mm_dia = round(c(min(F_est), min(F_sat)), 2))
print(tabla_f(negativos, "F2_valores_negativos"))

grafico_base("F2_bases", function() {
  par(mfrow = c(1, 2))
  for (K in c(K_EST, K_SAT)) {
    matplot(T_DIA, base_bspline(T_DIA, K), type = "l", lty = 1, xlab = "Día del año", ylab = expression(phi[k](t)),
            main = paste("Base B-spline cúbica, K =", K, if (K == K_EST) "(estaciones)" else "(satélite)"))
    abline(v = unique(nodos(K)), lty = 3, col = "grey60")
  }
})

Sigma_res_est <- crossprod(e$Res) / nrow(e$Res); Sigma_res_sat <- crossprod(s$Res) / nrow(s$Res)
Cor_res_est <- cov2cor(Sigma_res_est); Cor_res_sat <- cov2cor(Sigma_res_sat)
v_t <- diag(Sigma_res_est) / (1 - diag(e$H)); forma <- v_t / mean(v_t)
cat(sprintf("|Correlación residual| fuera de la diagonal: estaciones %.3f (rezago de 1 día %.3f) | satélite %.3f (rezago de 1 péntada %.3f)\n",
            mean(abs(Cor_res_est[upper.tri(Cor_res_est)])), mean(diag(Cor_res_est[-1, -365])),
            mean(abs(Cor_res_sat[upper.tri(Cor_res_sat)])), mean(diag(Cor_res_sat[-1, -72]))))
cat(sprintf("Varianza del error por día (estaciones): mínima %.1f el %s, máxima %.1f el %s\n",
            min(v_t), fecha_dia(which.min(v_t)), max(v_t), fecha_dia(which.max(v_t))))
tabla_f(data.table(dia = T_DIA, varianza_error = v_t), "F2_varianza_error_estaciones")

grafico_base("F2_covarianza_residual", function() {
  par(mfrow = c(1, 3))
  image(T_DIA, T_DIA, Cor_res_est, zlim = c(-.3, .3), col = hcl.colors(30, "RdBu", rev = TRUE),
        main = "Correlación residual: estaciones", xlab = "Día s", ylab = "Día t")
  image(T_SAT, T_SAT, Cor_res_sat, zlim = c(-1, 1), col = hcl.colors(30, "RdBu", rev = TRUE),
        main = "Correlación residual: satélite", xlab = "Día s", ylab = "Día t")
  plot(T_DIA, v_t, type = "l", xlab = "Día del año", ylab = "(mm/día)^2", main = "Varianza del error por día (estaciones)")
}, ancho = 14, alto = 4.8)

ee_banda <- function(Phi_eval, mapa, Sigma) { M <- Phi_eval %*% mapa; sqrt(rowSums((M %*% Sigma) * M)) }
ee_iid <- ee_banda(Phi_e, e$mapa, e$sigma2[i_ej] * diag(365))
ee_het <- ee_banda(Phi_e, e$mapa, diag(e$sigma2[i_ej] * forma))
racha <- function(m) { r <- rle(m); max(c(0, r$lengths[r$values])) }
inc <- which(!e$completas); ig <- inc[which.max(sapply(inc, function(k) racha(is.na(Y[k, ]))))]
ag <- ajuste_mco(Y[ig, ], Phi_e)
ee_g <- ee_banda(Phi_e, ag$mapa, ag$sigma2 * diag(sum(ag$obs)))
ee_g_completa <- ee_banda(Phi_e, e$mapa, ag$sigma2 * diag(365))
js <- which(id_sat$celda == emparejamiento[codigo == id_est$codigo[i_ej], celda] & id_sat$anio == id_est$anio[i_ej])
ee_s <- ee_banda(Phi_s365, s$mapa, s$sigma2[js] * diag(72))
cat(sprintf("Curva %s con hueco de %d días: el EE en los días faltantes es %.2f veces el de la misma curva sin hueco\n",
            id_est$curva_id[ig], sum(!ag$obs), mean((ee_g / ee_g_completa)[!ag$obs])))
cat(sprintf("EE medio de la banda: estación %.2f mm/día (sigma %.1f) | su píxel %.2f mm/día (sigma %.1f)\n",
            mean(ee_iid), sqrt(e$sigma2[i_ej]), mean(ee_s), sqrt(s$sigma2[js])))
tabla_f(data.table(curva_id = id_est$curva_id[i_ej], dia = T_DIA, curva = F_est[i_ej, ], ee_iid, ee_het,
                   inferior = F_est[i_ej, ] - 2 * ee_iid, superior = F_est[i_ej, ] + 2 * ee_iid), "F2_banda_ejemplo_estacion")

grafico_base("F2_bandas_punto_a_punto", function() {
  par(mfrow = c(1, 3))
  yl <- range(c(F_est[i_ej, ] - 2 * ee_het, Y[i_ej, ]), na.rm = TRUE)
  plot(T_DIA, Y[i_ej, ], pch = 16, cex = .4, col = "grey50", ylim = yl, xlab = "Día del año", ylab = "mm/día",
       main = paste("Estación", id_est$curva_id[i_ej]))
  lines(T_DIA, F_est[i_ej, ], lwd = 2, col = 4)
  lines(T_DIA, F_est[i_ej, ] + 2 * ee_iid, lty = 2, col = 4); lines(T_DIA, F_est[i_ej, ] - 2 * ee_iid, lty = 2, col = 4)
  lines(T_DIA, F_est[i_ej, ] + 2 * ee_het, lty = 3); lines(T_DIA, F_est[i_ej, ] - 2 * ee_het, lty = 3)
  legend("topright", c("Curva", "±2 EE, sigma^2 I", "±2 EE, sigma^2(t)"), lty = 1:3, col = c(4, 4, 1), bty = "n", cex = .8)
  cs <- F_sat[js, ]
  plot(T_SAT, X[js, ], pch = 16, cex = .7, col = "grey50", ylim = range(c(X[js, ], cs - 2 * ee_s, cs + 2 * ee_s)),
       xlab = "Día del año", ylab = "mm/día", main = paste("Su píxel, año", id_sat$anio[js]))
  lines(T_DIA, cs, lwd = 2, col = "darkorange"); lines(T_DIA, cs + 2 * ee_s, lty = 2); lines(T_DIA, cs - 2 * ee_s, lty = 2)
  plot(T_DIA, Y[ig, ], pch = 16, cex = .4, col = "grey50", ylim = c(-5, max(Y[ig, ], na.rm = TRUE)), xlab = "Día del año",
       ylab = "mm/día", main = paste("Estación", id_est$curva_id[ig], "-", sum(!ag$obs), "días faltantes"))
  lines(T_DIA, F_est[ig, ], lwd = 2, col = 4)
  lines(T_DIA, F_est[ig, ] + 2 * ee_g, lty = 2, col = 4); lines(T_DIA, F_est[ig, ] - 2 * ee_g, lty = 2, col = 4)
  rug(T_DIA[!ag$obs], col = 2, lwd = 2)
}, ancho = 15, alto = 4.8)

W_est <- gram(K_EST); W_sat <- gram(K_SAT)
for (w in list(list("estaciones", W_est), list("satélite", W_sat))) {
  W <- w[[2]]
  stopifnot(max(abs(W - t(W))) < 1e-10, min(eigen(W, symmetric = TRUE)$values) > 0, abs(sum(W) - 364) < 1e-8)
  cat(sprintf("W %s (%d x %d): simétrica, definida positiva (autovalor mínimo %.3f), 1'W1 = %.8f = longitud del dominio 364\n",
              w[[1]], nrow(W), ncol(W), min(eigen(W, symmetric = TRUE)$values), sum(W)))
}
t_fino <- seq(1, 365, by = .01); x_fino <- drop(base_bspline(t_fino, K_EST) %*% e$coef[i_ej, ])
cat(sprintf("Integral de x(t)^2: c'Wc = %.4f | trapecio con paso 0,01 = %.4f\n",
            drop(e$coef[i_ej, ] %*% W_est %*% e$coef[i_ej, ]), sum(diff(t_fino) * (head(x_fino, -1)^2 + tail(x_fino, -1)^2) / 2)))

grafico_base("F2_mismo_dominio", function() {
  par(mfrow = c(1, 1))
  plot(T_DIA, F_est[i_ej, ], type = "l", lwd = 2, col = 4, ylim = range(c(F_est[i_ej, ], F_sat[js, ], 0)),
       xlab = "Día del año (1-365)", ylab = "mm/día", main = paste("Estación", id_est$curva_id[i_ej], "y su píxel en el mismo dominio"))
  lines(T_DIA, F_sat[js, ], lwd = 2, col = "darkorange")
  points(T_SAT, X[js, ], pch = 16, cex = .6, col = "darkorange")
  legend("topright", c(paste("Estación, K =", K_EST), paste("Satélite, K =", K_SAT)), col = c(4, "darkorange"), lwd = 2, bty = "n")
}, ancho = 9, alto = 4.5)

# QUÉ DICEN LOS DATOS
#   La base propia coincide con splines::splineDesign (3e-16) y el atajo de la VC
#   con diag(S) da lo mismo que dejar cada día fuera (365 ajustes).
#   K se elige solo con entrenamiento (1 882 curvas completas de estación y
#   10 215 de satélite, con todos los años): VC y VCG coinciden, K = 8 en
#   estaciones y K = 11 en satélite. El ECM de ajuste baja siempre con K, por eso
#   no sirve para elegir. En el satélite el mínimo es nítido. En estaciones la VC
#   es casi plana entre K = 5 y 17 (a menos de 1 % del mínimo): la lluvia diaria
#   es tan ruidosa que más bases solo persiguen aguaceros.
#   Sensibilidad de K (estaciones): con K = 12, 16 o 24 la fecha mediana de los
#   picos se mueve a lo sumo 9 días (18 a 27 de abril; 24 a 29 de octubre), mucho
#   menos que la dispersión entre curvas (29 a 45 días), y los valores negativos
#   suben de 2,0 % a 5,4 %. Una base más flexible no cambia las conclusiones y
#   empeora la VC, así que se mantiene K = 8.
#   P-20 (D-51): frente a una media móvil de 31 días, K = 16 deja el pico a 6 y 6
#   días por curva y K = 8 a 12 y 9, lo que parecía pedir subir K. Pero el K que
#   gana depende de la referencia: con 15 y 31 días gana K = 20, con 61 días gana
#   K = 10 (y K = 16 queda peor que K = 8: 13 y 15 días). Las tres referencias
#   difieren entre sí 8 a 19 días en la fecha del pico, tanto como lo que se gana
#   subiendo K, y la fecha varía 32 a 36 días entre curvas. La fecha del pico no
#   sirve para elegir K: se queda el mínimo de la VCG, K = 8. Así tampoco hay que
#   rehacer el experimento de huecos de R5, hecho con K = 8.
#   S es simétrica, idempotente y tr(S) = K hasta 1e-14 en el peor de los ajustes
#   (también en las curvas con días faltantes, que usan su propia S).
#   El error de las estaciones no es homogéneo: su varianza va de 26 (mm/día)^2 en
#   agosto a 264 en octubre, la temporada húmeda. La correlación residual entre
#   días es casi nula (0,04; 0,08 a un día), así que la banda con sigma^2 I es
#   razonable. En el satélite es mayor (0,14): la péntada ya promedia cinco días.
#   Las bandas son punto a punto: cada día tiene 95 % de cobertura por separado,
#   no la curva completa. Un hueco de 15 días ensancha la banda solo 15 % en esos
#   días.
#   Valores negativos: 47 % de las curvas de estación baja de 0 en algún día (2 %
#   de los puntos), en la temporada seca de julio-agosto y en los bordes 1-ene y
#   31-dic (efecto de borde de la base no periódica). En el satélite casi no
#   ocurre (1 % de las curvas). Para modelar, la curva se recorta en 0 (D-23).
#   Por qué B-splines y no Fourier: una curva de un solo año no es periódica (el
#   31-dic no tiene por qué empalmar con el 1-ene del mismo año), el ejemplo de
#   clase con esta lluvia usa B-splines cúbicas y el enunciado las exige en la
#   red. El costo es el efecto de borde, que queda documentado.
#   W integra exactamente: 1'W1 = 364 y c'Wc coincide con el trapecio fino.

o <- list(id_est = id_est, id_sat = id_sat, coef_est = e$coef, coef_sat = s$coef, K_est = K_EST, K_sat = K_SAT,
             sigma2_est = e$sigma2, sigma2_sat = s$sigma2, n_obs_est = e$n_obs, W_est = W_est, W_sat = W_sat,
             t_sat = T_SAT, seleccion_K = sel, controles_S = controles, Sigma_res_est = Sigma_res_est, Sigma_res_sat = Sigma_res_sat,
             emparejamiento = emparejamiento, estaciones = est)

# =============================================================================
# 3. ANÁLISIS EXPLORATORIO DE DATOS FUNCIONALES (CADA FUENTE)  — rúbrica 2
# =============================================================================
# Media paso a paso (verificada con fda), varianza, desviación estándar,
# covarianza y correlación, derivadas y temporadas, profundidad FM y MBD,
# boxplot funcional con atípicas por estación y año, y contexto climatológico
# (vertiente, altitud y ENSO). Se usan todas las estaciones incluidas y todos
# los años 1981-2003; el satélite se analiza en los píxeles del Valle.
# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ┃ RÚBRICA 2 — Análisis exploratorio de datos funcionales
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
titulo("3. EDA funcional")

Phi_e <- base_bspline(T_DIA, o$K_est); Phi_s <- base_bspline(T_DIA, o$K_sat)
valle_s <- o$id_sat$en_valle
C_e <- o$coef_est; C_s <- o$coef_sat[valle_s, ]
id_e <- est[, .(codigo, nombre, vertiente, altitud)][o$id_est, on = "codigo"]
id_s <- o$id_sat[valle_s]
Ye <- C_e %*% t(Phi_e); Xs <- C_s %*% t(Phi_s)
cat("Estaciones:", nrow(Ye), "curvas,", uniqueN(id_e$codigo), "estaciones, años", paste(range(id_e$anio), collapse = "-"),
    "| satélite (Valle):", nrow(Xs), "curvas,", uniqueN(id_s$celda), "píxeles\n")
print(id_e[, .N, by = vertiente])

titulo("3.1 Media funcional")
media_e <- drop(Phi_e %*% colMeans(C_e)); media_s <- drop(Phi_s %*% colMeans(C_s))
fd_e <- fda::fd(t(C_e), fda::create.bspline.basis(RANGO, nbasis = o$K_est, norder = ORDEN))
fd_s <- fda::fd(t(C_s), fda::create.bspline.basis(RANGO, nbasis = o$K_sat, norder = ORDEN))
verif_media <- data.table(fuente = c("estaciones", "satélite"),
                          vs_curvas_evaluadas = c(max(abs(media_e - colMeans(Ye))), max(abs(media_s - colMeans(Xs)))),
                          vs_fda_mean.fd = c(max(abs(media_e - fda::eval.fd(T_DIA, fda::mean.fd(fd_e)))),
                                             max(abs(media_s - fda::eval.fd(T_DIA, fda::mean.fd(fd_s))))))
print(tabla_f(verif_media, "F3_verificacion_media"))
stopifnot(verif_media$vs_fda_mean.fd < 1e-10)
cat(sprintf("Lluvia media: estaciones %.0f mm/año | satélite en el Valle %.0f mm/año\n", mean(media_e) * 365, mean(media_s) * 365))

titulo("3.2 Varianza, desviación estándar, covarianza y correlación")
momentos <- function(C, Phi) {
  Cc <- sweep(C, 2, colMeans(C)); Sc <- crossprod(Cc) / (nrow(C) - 1)
  covar <- Phi %*% Sc %*% t(Phi)
  list(var = diag(covar), de = sqrt(diag(covar)), cov = covar, cor = cov2cor(covar))
}
me <- momentos(C_e, Phi_e); ms <- momentos(C_s, Phi_s)
cat("Covarianza Phi Sigma_c Phi' vs cov() de las curvas evaluadas (estaciones):", format(max(abs(me$cov - cov(Ye))), digits = 2), "\n")
resumen_disp <- data.table(fuente = c("estaciones", "satélite"),
                           dia_DE_max = c(which.max(me$de), which.max(ms$de)), DE_max = c(max(me$de), max(ms$de)),
                           dia_DE_min = c(which.min(me$de), which.min(ms$de)), DE_min = c(min(me$de), min(ms$de)),
                           CV_medio = c(mean(me$de / media_e), mean(ms$de / media_s)),
                           cor_abr_oct = c(me$cor[105, 290], ms$cor[105, 290]), cor_abr_jul = c(me$cor[105, 200], ms$cor[105, 200]),
                           cor_ene_dic = c(me$cor[15, 350], ms$cor[15, 350]))
resumen_disp[, `:=`(fecha_DE_max = fecha_dia(dia_DE_max), fecha_DE_min = fecha_dia(dia_DE_min))]
print(tabla_f(resumen_disp, "F3_dispersion_resumen"))

# ¿Qué mide la covarianza con todas las estaciones juntas? Se separa por vertiente y se
# repite con ANOMALÍAS (cada curva menos la curva media de su estación o de su píxel),
# que quitan las diferencias entre sitios y dejan la variación de un año a otro.
anom_coef <- function(C, g) { g <- as.character(g); m <- rowsum(C, g) / as.vector(table(g)[rownames(rowsum(C, g))]); C - m[g, ] }
espectro <- function(C, Phi, nombre) {
  W <- gram(ncol(C)); ev <- eigen(W, symmetric = TRUE); W12 <- ev$vectors %*% diag(sqrt(ev$values)) %*% t(ev$vectors)
  l <- eigen(W12 %*% cov(C) %*% W12, symmetric = TRUE, only.values = TRUE)$values; m <- momentos(C, Phi)
  data.table(grupo = nombre, curvas = nrow(C), DE_media_mm_dia = round(mean(m$de), 2), pct_varianza_1a_componente = round(100 * l[1] / sum(l), 1),
             pct_varianza_2a_componente = round(100 * l[2] / sum(l), 1), cor_abr_oct = round(m$cor[105, 290], 2), cor_abr_jul = round(m$cor[105, 200], 2), cor_minima = round(min(m$cor), 2))
}
Ca_e <- anom_coef(C_e, id_e$codigo); Ca_s <- anom_coef(C_s, id_s$celda)
cov_grupos <- rbind(espectro(C_e, Phi_e, "Estaciones: todas"), espectro(C_e[id_e$vertiente != "Pacífico", ], Phi_e, "Estaciones: Cauca"),
                    espectro(C_e[id_e$vertiente == "Pacífico", ], Phi_e, "Estaciones: Pacífico"), espectro(Ca_e, Phi_e, "Estaciones: anomalías por estación"),
                    espectro(C_s, Phi_s, "Satélite: todos los píxeles"), espectro(C_s[id_s$vertiente != "Pacífico", ], Phi_s, "Satélite: Cauca"),
                    espectro(C_s[id_s$vertiente == "Pacífico", ], Phi_s, "Satélite: Pacífico"), espectro(Ca_s, Phi_s, "Satélite: anomalías por píxel"))
print(tabla_f(cov_grupos, "F3_covarianza_por_vertiente_y_anomalias"))
ma_e <- momentos(Ca_e, Phi_e); ma_s <- momentos(Ca_s, Phi_s)
cat(sprintf("Parte de la varianza que es diferencia ENTRE sitios: estaciones %.0f %% | satélite %.0f %%\n",
            100 * (1 - sum(ma_e$var) / sum(me$var)), 100 * (1 - sum(ma_s$var) / sum(ms$var))))
grafico_base("F3_correlacion_anomalias", function() {
  layout(matrix(1:3, 1), widths = c(5, 5, 1)); par(mar = c(4, 4, 3, 1)); colores <- hcl.colors(64, "Blue-Red 3")
  image(T_DIA, T_DIA, ma_e$cor, zlim = c(-1, 1), col = colores, main = "Correlación de las anomalías: estaciones", xlab = "Día s", ylab = "Día t")
  image(T_DIA, T_DIA, ma_s$cor, zlim = c(-1, 1), col = colores, main = "Correlación de las anomalías: satélite", xlab = "Día s", ylab = "Día t")
  par(mar = c(4, 0, 3, 4)); escala <- seq(-1, 1, length.out = 65)
  image(c(0, 1), escala, matrix(rep(escala, each = 2), nrow = 2), col = colores, zlim = c(-1, 1), axes = FALSE, xlab = "", ylab = ""); axis(4, las = 1, cex.axis = .8); box(); layout(1)
}, ancho = 13, alto = 5.5)

grafico_base("F3_media_y_dispersion", function() {
  par(mfrow = c(1, 4)); set.seed(1)
  for (f in list(list(Ye, media_e, me, "Estaciones", 4), list(Xs, media_s, ms, "Satélite (Valle)", "darkorange"))) {
    matplot(T_DIA, t(f[[1]][sample(nrow(f[[1]]), 250), ]), type = "l", lty = 1, col = "grey85", ylim = c(-2, 22),
            xlab = "Día del año", ylab = "mm/día", main = paste(f[[4]], ": media ± 2 DE"))
    lines(T_DIA, f[[2]], lwd = 3, col = f[[5]])
    lines(T_DIA, f[[2]] + 2 * f[[3]]$de, lty = 2, lwd = 2, col = 2); lines(T_DIA, pmax(f[[2]] - 2 * f[[3]]$de, -2), lty = 2, lwd = 2, col = 2)
  }
  plot(T_DIA, me$de, type = "l", lwd = 2, col = 4, ylim = range(c(0, me$de, ms$de)), xlab = "Día del año", ylab = "mm/día",
       main = "Desviación estándar funcional")
  lines(T_DIA, ms$de, lwd = 2, col = "darkorange")
  legend("topleft", c("Estaciones", "Satélite"), col = c(4, "darkorange"), lwd = 2, bty = "n")
  # Coeficiente de variación funcional CV(t) = DE(t) / media(t) (clase, semanas 12-13): compara la
  # dispersión de fuentes con escalas distintas.
  plot(T_DIA, me$de / media_e, type = "l", lwd = 2, col = 4, ylim = range(c(0, me$de / media_e, ms$de / media_s)), xlab = "Día del año", ylab = "CV(t)",
       main = "Coeficiente de variación funcional")
  lines(T_DIA, ms$de / media_s, lwd = 2, col = "darkorange")
}, ancho = 19, alto = 4.8)
for (f in list(list(me, "estaciones"), list(ms, "satelite"))) {
  grafico_base(paste0("F3_covarianza_", f[[2]]), function() mapa_calor(T_DIA, f[[1]]$cov, paste("Covarianza funcional:", f[[2]]), paleta = "YlOrRd"), ancho = 7, alto = 5.5)
  grafico_base(paste0("F3_correlacion_", f[[2]]), function() mapa_calor(T_DIA, f[[1]]$cor, paste("Correlación funcional:", f[[2]]), c(-1, 1)), ancho = 7, alto = 5.5)
}

titulo("3.3 Derivadas y temporadas de lluvia")
DPhi_e <- deriv_base(T_DIA, o$K_est); DPhi_s <- deriv_base(T_DIA, o$K_sat)
ref_d <- splines::splineDesign(nodos(o$K_est), T_DIA, ord = ORDEN, derivs = rep(1, 365), outer.ok = TRUE)
cat("Derivada de la base propia vs splines::splineDesign:", format(max(abs(DPhi_e - ref_d)), digits = 2), "\n")
dmu_e <- drop(DPhi_e %*% colMeans(C_e)); dmu_s <- drop(DPhi_s %*% colMeans(C_s))
cruces <- function(media, dmu, fuente) {
  k <- which(diff(sign(dmu)) != 0)
  data.table(fuente, dia = k + 1, fecha = fecha_dia(k + 1), tipo = fifelse(dmu[k] > 0, "máximo", "mínimo"), mm_dia = round(media[k + 1], 2))
}
extremos_media <- rbind(cruces(media_e, dmu_e, "estaciones"), cruces(media_s, dmu_s, "satélite"))
print(tabla_f(extremos_media, "F3_extremos_de_la_media"))

temporadas <- function(Y, D) {
  t(sapply(seq_len(nrow(Y)), function(i) {
    y <- Y[i, ]; d <- D[i, ]
    p1 <- 59 + which.max(y[60:200]); p2 <- 229 + which.max(y[230:365])
    r <- c(inicio_1 = if (p1 > 21) 19 + which.max(d[20:(p1 - 1)]) else NA_integer_, pico_1 = p1, fin_1 = p1 + which.min(d[(p1 + 1):229]),
           inicio_2 = 200 + which.max(d[201:(p2 - 1)]), pico_2 = p2, fin_2 = if (p2 < 349) p2 + which.min(d[(p2 + 1):350]) else NA_integer_)
    r[c("inicio_1", "fin_2")][r[c("inicio_1", "fin_2")] %in% c(20, 350)] <- NA
    r
  }))
}
tem_e <- temporadas(Ye, C_e %*% t(DPhi_e)); tem_s <- temporadas(Xs, C_s %*% t(DPhi_s))
fechas_tem <- rbind(data.table(fuente = "estaciones", vertiente = id_e$vertiente, tem_e), data.table(fuente = "satélite", vertiente = id_s$vertiente, tem_s))
fechas_tem <- melt(fechas_tem, id.vars = c("fuente", "vertiente"), variable.name = "momento", value.name = "dia")
resumen_tem <- fechas_tem[, .(pct_no_identificable = round(100 * mean(is.na(dia)), 1), dia_mediano = median(dia, na.rm = TRUE), p25 = quantile(dia, .25, na.rm = TRUE), p75 = quantile(dia, .75, na.rm = TRUE)),
                          by = .(fuente, vertiente, momento)]
resumen_tem[, `:=`(fecha_mediana = fecha_dia(dia_mediano), rango_intercuartil = paste(fecha_dia(p25), "a", fecha_dia(p75)))]
print(tabla_f(dcast(resumen_tem, momento ~ fuente + vertiente, value.var = "fecha_mediana"), "F3_temporadas_fechas"))
tabla_f(resumen_tem, "F3_temporadas_detalle")

grafico_base("F3_derivada_media", function() {
  par(mfrow = c(1, 2))
  for (f in list(list(media_e, dmu_e, "Estaciones", 4), list(media_s, dmu_s, "Satélite (Valle)", "darkorange"))) {
    plot(T_DIA, f[[2]], type = "l", lwd = 2, col = f[[4]], xlab = "Día del año", ylab = "mm/día por día", main = paste(f[[3]], ": derivada de la media"))
    abline(h = 0, lty = 2, col = "grey40")
    k <- which(diff(sign(f[[2]])) != 0) + 1; abline(v = k, lty = 3, col = "grey60")
    axis(3, at = k, labels = fecha_dia(k), cex.axis = .7, las = 2, tick = FALSE, line = -.8)
  }
}, ancho = 12, alto = 4.5)
grafico_f(ggplot(fechas_tem[!is.na(dia)], aes(dia, momento, fill = fuente)) + geom_boxplot(outlier.size = .2, alpha = .7) +
            facet_wrap(~vertiente) + scale_x_continuous(breaks = c(1, 60, 121, 182, 244, 305, 365), labels = fecha_dia(c(1, 60, 121, 182, 244, 305, 365))) +
            labs(x = NULL, y = NULL, fill = NULL, title = "Inicio, pico y fin de cada temporada en cada curva (derivada)"),
          "F3_temporadas_por_curva", ancho = 11, alto = 5)

titulo("3.4 Profundidad, mediana y media recortada")
Yp <- Ye[1:40, ]; pares <- combn(40, 2)
mbd_bruta <- sapply(1:40, function(i) mean(apply(pares, 2, function(p) mean(pmin(Yp[p[1], ], Yp[p[2], ]) <= Yp[i, ] & Yp[i, ] <= pmax(Yp[p[1], ], Yp[p[2], ])))))
cat("MBD por rangos vs definición con todas las bandas de dos curvas:", format(max(abs(profundidades(Yp)$MBD - mbd_bruta)), digits = 2), "\n")
pr_e <- profundidades(Ye); pr_s <- profundidades(Xs)
centralidad <- function(Y, pr) {
  dentro <- pr$FM >= quantile(pr$FM, .25)
  list(mediana_FM = Y[which.max(pr$FM), ], mediana_MBD = Y[which.max(pr$MBD), ], recortada_25 = colMeans(Y[dentro, , drop = FALSE]))
}
ce_e <- centralidad(Ye, pr_e); ce_s <- centralidad(Xs, pr_s)
cat("Curva más profunda (MBD), estaciones:", id_e$curva_id[which.max(pr_e$MBD)], id_e$nombre[which.max(pr_e$MBD)],
    "| (FM):", id_e$curva_id[which.max(pr_e$FM)], "| satélite (MBD): año", id_s$anio[which.max(pr_s$MBD)], "píxel", id_s$celda[which.max(pr_s$MBD)], "\n")
cat(sprintf("Correlación de rangos FM vs MBD: estaciones %.2f | satélite %.2f\n",
            cor(pr_e$FM, pr_e$MBD, method = "spearman"), cor(pr_s$FM, pr_s$MBD, method = "spearman")))
grafico_base("F3_centralidad", function() {
  par(mfrow = c(1, 2))
  for (f in list(list(media_e, ce_e, "Estaciones"), list(media_s, ce_s, "Satélite (Valle)"))) {
    matplot(T_DIA, cbind(f[[1]], f[[2]]$mediana_FM, f[[2]]$mediana_MBD, f[[2]]$recortada_25), type = "l", lty = c(1, 2, 3, 1),
            col = c(4, 2, 1, 3), lwd = 2, xlab = "Día del año", ylab = "mm/día", main = paste("Centralidad:", f[[3]]))
    legend("topright", c("Media", "Mediana FM", "Mediana MBD", "Media recortada 25 % (FM)"), col = c(4, 2, 1, 3), lty = c(1, 2, 3, 1), lwd = 2, bty = "n", cex = .75)
  }
}, ancho = 12, alto = 4.5)

titulo("3.5 Boxplot funcional y curvas atípicas")
bx_e <- boxplot_funcional(Ye, pr_e$MBD); bx_s <- boxplot_funcional(Xs, pr_s$MBD)
grupos_v <- split(seq_len(nrow(Ye)), id_e$vertiente)
bx_v <- lapply(grupos_v, function(i) boxplot_funcional(Ye[i, ], profundidades(Ye[i, ])$MBD))
fb_e <- fda::fbplot(t(Ye), T_DIA, depth = pr_e$MBD, plot = FALSE)
fb_s <- fda::fbplot(t(Xs), T_DIA, depth = pr_s$MBD, plot = FALSE)
verif_fb <- c(estaciones = setequal(as.integer(fb_e$outpoint), which(bx_e$atip)), satelite = setequal(as.integer(fb_s$outpoint), which(bx_s$atip)))
cat("Atípicas iguales a fda::fbplot:", paste(names(verif_fb), verif_fb, collapse = " | "), "\n")
stopifnot(verif_fb)

atip_e <- id_e[bx_e$atip, .(codigo, nombre, anio, vertiente, altitud, particion)]
atip_e[, `:=`(MBD = round(pr_e$MBD[bx_e$atip], 3), P_anual_mm = round(rowMeans(Ye[bx_e$atip, , drop = FALSE]) * 365), dia_max = apply(Ye[bx_e$atip, , drop = FALSE], 1, which.max))]
atip_v <- rbindlist(lapply(names(grupos_v), function(v) id_e[grupos_v[[v]][bx_v[[v]]$atip], .(codigo, nombre, anio, vertiente)]))
tabla_f(atip_e[order(-P_anual_mm)], "F3_atipicas_estaciones_todas")
tabla_f(atip_v, "F3_atipicas_estaciones_por_vertiente")
cat(sprintf("Atípicas, todas las estaciones juntas: %d de %d curvas (%d estaciones; %d del Pacífico)\n",
            nrow(atip_e), nrow(Ye), uniqueN(atip_e$codigo), atip_e[vertiente == "Pacífico", .N]))
print(head(atip_e[, .(curvas = .N, anios = paste(sort(anio), collapse = " "), P_anual_mediana = median(P_anual_mm)), by = .(codigo, nombre, vertiente)][order(-curvas)], 12))
cat("Atípicas dentro de cada vertiente:\n"); print(atip_v[, .(curvas = .N, estaciones = uniqueN(codigo)), by = vertiente])
print(atip_v[vertiente != "Pacífico", .N, by = .(codigo, nombre, anio)])
atip_s <- id_s[bx_s$atip]
cat(sprintf("Atípicas del satélite: %d de %d curvas píxel-año (%.1f %%)\n", nrow(atip_s), nrow(Xs), 100 * mean(bx_s$atip)))
print(atip_s[, .N, by = .(vertiente)]); print(head(atip_s[, .N, by = anio][order(-N)], 6))
tabla_f(atip_s[, .(anio, celda, lon, lat, vertiente, n_tope)], "F3_atipicas_satelite")

anios_e <- sort(unique(id_e$anio))
media_anual <- function(Y, g, a) t(sapply(a, function(x) colMeans(Y[g == x, , drop = FALSE])))
Ya_e <- media_anual(Ye[id_e$vertiente != "Pacífico", ], id_e$anio[id_e$vertiente != "Pacífico"], anios_e)
Ya_s <- media_anual(Xs, id_s$anio, ANIOS)
pa_e <- profundidades(Ya_e)$MBD; pa_s <- profundidades(Ya_s)$MBD
bxa_e <- boxplot_funcional(Ya_e, pa_e); bxa_s <- boxplot_funcional(Ya_s, pa_s)
prof_anual <- oni[data.table(anio = anios_e, MBD_estaciones_Cauca = round(pa_e, 3), atipico_estaciones = bxa_e$atip,
                             MBD_satelite = round(pa_s, 3), atipico_satelite = bxa_s$atip,
                             P_estaciones_Cauca = round(rowMeans(Ya_e) * 365), P_satelite = round(rowMeans(Ya_s) * 365)), on = "anio"]
print(tabla_f(prof_anual[order(MBD_estaciones_Cauca)], "F3_profundidad_anual"))
cat(sprintf("Correlación (Spearman) entre la profundidad anual de estaciones y de satélite: %.2f\n", cor(pa_e, pa_s, method = "spearman")))

grafico_base("F3_boxplot_funcional", function() {
  par(mfrow = c(2, 3))
  dibujar_boxplot(Ye, bx_e, paste("Estaciones, todas:", sum(bx_e$atip), "atípicas"))
  for (v in names(bx_v)) dibujar_boxplot(Ye[grupos_v[[v]], ], bx_v[[v]], paste("Estaciones,", v, ":", sum(bx_v[[v]]$atip), "atípicas"))
  dibujar_boxplot(Xs, bx_s, paste("Satélite (Valle):", sum(bx_s$atip), "atípicas"))
  dibujar_boxplot(Ya_e, bxa_e, "Curvas medias anuales: estaciones del Cauca")
  dibujar_boxplot(Ya_s, bxa_s, "Curvas medias anuales: satélite")
}, ancho = 15, alto = 9)

titulo("3.6 Contexto: vertiente, altitud y ENSO")
id_e[, banda_altitud := cut(altitud, c(-Inf, 1200, 2000, Inf), labels = c("< 1200 m", "1200-2000 m", "> 2000 m"))]
id_e <- oni[, .(anio, fase)][id_e, on = "anio"]; id_s <- oni[, .(anio, fase)][id_s, on = "anio"]

# Curvas por año, como en la tesis del profesor (figuras 6.2 y 6.3): un panel por año con las
# curvas en gris, el promedio en negro y la desviación estándar punteada. Todos los paneles
# tienen el mismo eje, para que los años se comparen a simple vista; el título lleva la fase ENSO.
por_anio <- function(M, anios, fases, cada, titulo_fig) {
  lim <- c(0, quantile(M, .995))
  par(mfrow = c(6, 4), mar = c(3, 3.2, 2, .6), mgp = c(1.9, .5, 0), oma = c(0, 0, 3, 0))
  for (a in ANIOS) {
    k <- anios == a
    matplot(T_DIA, t(M[k, , drop = FALSE]), type = "l", lty = 1, col = "grey80", ylim = lim, xaxt = "n", xlab = "Día del año",
            ylab = "Precipitación (mm/día)", main = sprintf("%d · %s", a, fases[which(k)[1]]), cex.axis = .75, cex.lab = .8)
    axis(1, at = c(1, 91, 182, 274, 365), labels = c("ene", "abr", "jul", "oct", "dic"), cex.axis = .75)
    lines(T_DIA, colMeans(M[k, , drop = FALSE]), lwd = 2); lines(T_DIA, apply(M[k, , drop = FALSE], 2, sd), lwd = 1.5, lty = 2)
  }
  plot.new(); legend("center", c(cada, "Promedio", "Desviación estándar"), col = c("grey70", "black", "black"),
                     lty = c(1, 1, 2), lwd = c(1, 2, 1.5), bty = "n", cex = 1)
  mtext(titulo_fig, outer = TRUE, cex = 1.1, font = 2)
}
grafico_base("F3_por_anio_estaciones", function() por_anio(Ye, id_e$anio, id_e$fase, "Cada estación",
             "Datos funcionales de las estaciones: promedio y desviación estándar por año, 1981–2003"), ancho = 11, alto = 15)
grafico_base("F3_por_anio_satelite", function() por_anio(Xs, id_s$anio, id_s$fase, "Cada píxel del Valle",
             "Datos funcionales del satélite: promedio y desviación estándar por año, 1981–2003"), ancho = 11, alto = 15)
resumen_anio <- oni[data.table(anio = ANIOS, estaciones = sapply(ANIOS, function(a) sum(id_e$anio == a)),
                               P_estaciones = sapply(ANIOS, function(a) round(mean(Ye[id_e$anio == a, ]) * 365)),
                               P_satelite = sapply(ANIOS, function(a) round(mean(Xs[id_s$anio == a, ]) * 365)),
                               DE_estaciones = sapply(ANIOS, function(a) round(mean(apply(Ye[id_e$anio == a, , drop = FALSE], 2, sd)), 2)),
                               DE_satelite = sapply(ANIOS, function(a) round(mean(apply(Xs[id_s$anio == a, , drop = FALSE], 2, sd)), 2))), on = "anio"]
print(tabla_f(resumen_anio, "F3_por_anio_resumen"))
cat(sprintf("Lluvia anual media vs ONI (23 años): estaciones r = %.2f | satélite r = %.2f\n",
            cor(resumen_anio$P_estaciones, resumen_anio$ONI_anual), cor(resumen_anio$P_satelite, resumen_anio$ONI_anual)))
# Orden de los grupos fijo (niveles del factor u orden "radix", que no depende del locale):
# así las tablas salen con las mismas filas en cualquier máquina.
niveles_de <- function(g) if (is.factor(g)) levels(droplevels(g)) else sort(unique(g[!is.na(g)]), method = "radix")
grupo_media <- function(Y, g) {
  k <- !is.na(g); g <- factor(as.character(g), levels = niveles_de(g))
  sapply(split(seq_len(nrow(Y))[k], g[k], drop = TRUE), function(i) colMeans(Y[i, , drop = FALSE]))
}
m_vert_e <- grupo_media(Ye, id_e$vertiente); m_vert_s <- grupo_media(Xs, id_s$vertiente)
m_alt <- grupo_media(Ye[id_e$vertiente != "Pacífico", ], id_e$banda_altitud[id_e$vertiente != "Pacífico"])
m_enso_e <- grupo_media(Ye[id_e$vertiente != "Pacífico", ], id_e$fase[id_e$vertiente != "Pacífico"])
m_enso_s <- grupo_media(Xs[id_s$vertiente != "Pacífico", ], id_s$fase[id_s$vertiente != "Pacífico"])
k_sd <- id_s$vertiente != "Pacífico" & !id_s$marcado
m_enso_s_sd <- grupo_media(Xs[k_sd, ], id_s$fase[k_sd])
contexto <- rbind(
  data.table(fuente = "estaciones", agrupacion = "vertiente", grupo = colnames(m_vert_e), P_anual = colMeans(m_vert_e) * 365,
             curvas = as.vector(table(id_e$vertiente)[colnames(m_vert_e)])),
  data.table(fuente = "satélite", agrupacion = "vertiente", grupo = colnames(m_vert_s), P_anual = colMeans(m_vert_s) * 365,
             curvas = as.vector(table(id_s$vertiente)[colnames(m_vert_s)])),
  data.table(fuente = "estaciones (Cauca)", agrupacion = "altitud", grupo = colnames(m_alt), P_anual = colMeans(m_alt) * 365,
             curvas = as.vector(table(id_e[vertiente != "Pacífico", banda_altitud])[colnames(m_alt)])),
  data.table(fuente = "estaciones (Cauca)", agrupacion = "ENSO", grupo = colnames(m_enso_e), P_anual = colMeans(m_enso_e) * 365,
             curvas = as.vector(table(id_e[vertiente != "Pacífico", fase])[colnames(m_enso_e)])),
  data.table(fuente = "satélite (Cauca)", agrupacion = "ENSO", grupo = colnames(m_enso_s), P_anual = colMeans(m_enso_s) * 365,
             curvas = as.vector(table(id_s[vertiente != "Pacífico", fase])[colnames(m_enso_s)])),
  data.table(fuente = "satélite (Cauca) sin 1985, 1988 y 2000", agrupacion = "ENSO", grupo = colnames(m_enso_s_sd), P_anual = colMeans(m_enso_s_sd) * 365,
             curvas = as.vector(table(id_s[k_sd, fase])[colnames(m_enso_s_sd)])))
cat("Años marcados (1985, 1988 y 2000) y su fase ENSO:\n"); print(oni[anio %in% ANIOS_MARCADOS])
contexto[, P_anual := round(P_anual)]
print(tabla_f(contexto, "F3_contexto_vertiente_altitud_ENSO"))
picos_enso <- rbindlist(lapply(colnames(m_enso_e), function(f) data.table(fase = f, pico_1 = fecha_dia(59 + which.max(m_enso_e[60:200, f])),
  mm_pico_1 = round(max(m_enso_e[60:200, f]), 2), pico_2 = fecha_dia(229 + which.max(m_enso_e[230:365, f])), mm_pico_2 = round(max(m_enso_e[230:365, f]), 2),
  minimo_jul_ago = round(min(m_enso_e[160:260, f]), 2))))
print(tabla_f(picos_enso, "F3_picos_por_fase_ENSO"))

grafico_base("F3_contexto", function() {
  par(mfrow = c(1, 3))
  matplot(T_DIA, cbind(m_vert_e, m_vert_s), type = "l", lty = c(1, 1, 2, 2), lwd = 2, col = c(4, "darkgreen", 4, "darkgreen"),
          xlab = "Día del año", ylab = "mm/día", main = "Media por vertiente")
  legend("topleft", c(paste("Estaciones,", colnames(m_vert_e)), paste("Satélite,", colnames(m_vert_s))), col = c(4, "darkgreen", 4, "darkgreen"),
         lty = c(1, 1, 2, 2), lwd = 2, bty = "n", cex = .75)
  matplot(T_DIA, m_alt, type = "l", lty = 1, lwd = 2, col = c("goldenrod", "darkorange3", "brown"), xlab = "Día del año", ylab = "mm/día",
          main = "Estaciones del Cauca por altitud")
  legend("topleft", colnames(m_alt), col = c("goldenrod", "darkorange3", "brown"), lwd = 2, bty = "n")
  matplot(T_DIA, m_enso_e, type = "l", lty = 1, lwd = 2, col = c(1, 4, 2)[match(colnames(m_enso_e), c("Neutro", "Niña", "Niño"))],
          xlab = "Día del año", ylab = "mm/día", main = "Estaciones del Cauca por fase ENSO")
  legend("topleft", colnames(m_enso_e), col = c(1, 4, 2)[match(colnames(m_enso_e), c("Neutro", "Niña", "Niño"))], lwd = 2, bty = "n")
}, ancho = 15, alto = 4.8)

# QUÉ DICEN LOS DATOS
#   Representatividad: se usan las 176 estaciones incluidas y sus 3 299 curvas de
#   1981-2003, sin submuestrear (la sección 8 del archivo de datos muestra que
#   las incluidas representan a las 220 con datos), y los 676 píxeles del Valle
#   en los 23 años.
#   Media: Phi c_barra coincide con fda::mean.fd (5e-15). Ambas fuentes muestran
#   el régimen bimodal: estaciones con máximos el 22-abr (6,3 mm/día) y el 28-oct
#   (7,2), mínimos el 25-ene y el 1-ago (2,2). El satélite pone los máximos 1 a 3
#   semanas después (13-may y 5-nov) y con casi el doble de magnitud, porque su
#   media pesa por área y la vertiente Pacífica (43 % de los píxeles) llueve
#   4 000 mm/año. Por vertiente, el satélite supera a las estaciones en ambas:
#   1 735 frente a 1 551 mm en el Cauca y 3 998 frente a 3 105 en el Pacífico.
#   Dispersión: la DE de las estaciones es alta en la segunda temporada y mínima
#   el 27-feb; el máximo del 1-ene es efecto de borde. El CV medio es 0,96: la
#   variación entre curvas es tan grande como la media.
#   Qué mide la covarianza: con todas las estaciones juntas la correlación es
#   positiva todo el año (mínimo 0,39) y la primera componente explica 76 %. Eso
#   es sobre todo "el nivel de lluvia de cada sitio": el 74 % de la varianza de
#   las estaciones (88 % en el satélite) es diferencia ENTRE sitios. Dentro del
#   Cauca la primera componente baja a 50 % y la correlación abril-octubre a
#   0,54. Con anomalías (cada curva menos la media de su estación) queda la
#   variación de un año a otro: la primera componente explica solo 29 %, abril y
#   octubre casi no se relacionan (0,14) y aparecen correlaciones negativas
#   (-0,27). Es decir: un sitio húmedo lo es todo el año, pero un año con una
#   primera temporada lluviosa no anuncia una segunda temporada lluviosa.
#   Derivadas y temporadas: en el Cauca la primera temporada empieza a finales de
#   febrero, tiene su pico a mediados de abril y termina a mediados de junio; la
#   segunda empieza a mediados de septiembre y tiene su pico el 1-nov. En 62 % de
#   las curvas de estación la lluvia sigue bajando hasta el borde del año y el
#   fin de la segunda temporada no se puede fechar. El satélite fecha el inicio y
#   el pico de la primera temporada 2 a 3 semanas tarde.
#   Profundidad: FM y MBD ordenan casi igual (Spearman 0,99 y 1,00). La curva más
#   profunda es Betania 2002. Boxplot funcional (mismas atípicas que fda::fbplot):
#   con todas las estaciones juntas hay 63 atípicas, 62 de cinco estaciones de la
#   costa Pacífica (Bahía Málaga, Yurumanguí, San Francisco, Triana y El Tabano,
#   7 000 a 12 600 mm/año): no son errores sino otro clima. Por eso se separa por
#   vertiente: en el Cauca solo quedan Tenjo 1982 (4 482 mm, candidata a revisión)
#   y Las Brisas 1981. Las 330 atípicas del satélite son todas del Pacífico y 229
#   tocan el tope de 200 mm. Ningún año es atípico en las curvas medias anuales;
#   los menos profundos de las estaciones son 1988, 1984, 1983 y 1999 (tres Niña).
#   Orografía: en el Cauca la franja de 1 200-2 000 m es la más lluviosa
#   (1 715 mm), por encima de > 2 000 m (1 607) y del fondo del valle (1 333).
#   ENSO: en las estaciones del Cauca La Niña aporta +11 % (1 718 mm) y El Niño
#   -10 % (1 387) frente a neutro (1 544); la Niña se nota en la segunda temporada
#   (7,1 frente a 6,3 mm/día) y en una sequía de mitad de año más húmeda (2,4
#   frente a 1,4). El satélite NO muestra esa señal: da menos lluvia en años Niña
#   (1 664) que en neutros (1 767). Es el mismo error ligado al ENSO de la
#   sección 1, y es la razón para darle el ONI al modelo.

eda <- list(media_e = media_e, media_s = media_s, me = me, ms = ms, pr_e = pr_e, pr_s = pr_s, atip_e = atip_e, atip_s = atip_s,
             temporadas = resumen_tem, prof_anual = prof_anual, contexto = contexto)

# =============================================================================
# 4. COMPARACIÓN FUNCIONAL POR FUENTE (ESTACIÓN VS. SU PÍXEL)  — rúbrica 3
# =============================================================================
# Pares estación-año / píxel-año con el emparejamiento de la sección 1.
# Curva de diferencia d(t) = estación - satélite con banda d_barra(t) ± 2EE(t);
# el EE se calcula entre estaciones (las curvas de una misma estación no son
# independientes). Correlación cruzada funcional y su diagonal, covarianzas de
# cada fuente, comparación por vertiente, altitud y fase ENSO, regresión concurrente
# (en qué época del año el satélite explica a la estación) y ANOVA funcional del error.
# Entran todos los años; 1985, 1988 y 2000 se muestran además aparte, como descripción.
# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ┃ RÚBRICA 3 — Comparación funcional por fuente
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
titulo("4. Comparación funcional por fuente")

Phi_e <- base_bspline(T_DIA, o$K_est); Phi_s <- base_bspline(T_DIA, o$K_sat)

pares <- o$emparejamiento[, .(codigo, nombre, vertiente, altitud, celda, dist_km, estaciones_en_el_pixel)][o$id_est[, .(fila_est = .I, codigo, anio, particion, curva_id)], on = "codigo"]
pares <- o$id_sat[, .(fila_sat = .I, anio, celda, marcado)][pares, on = .(anio, celda)]
pares <- oni[, .(anio, ONI_anual, fase)][pares, on = "anio"]
pares[, banda_altitud := cut(altitud, c(-Inf, 1200, 2000, Inf), labels = c("< 1200 m", "1200-2000 m", "> 2000 m"))]
stopifnot(!anyNA(pares$fila_sat), nrow(pares) == nrow(o$id_est))
cat("Pares estación-año / píxel-año:", nrow(pares), "en", uniqueN(pares$codigo), "estaciones y", uniqueN(pares$celda), "píxeles |",
    "en 1985, 1988 y 2000:", pares[marcado == TRUE, .N], "\n")
print(pares[, .(pares = .N, estaciones = uniqueN(codigo)), by = .(particion, marcado)][order(particion, marcado)])
tabla_f(pares[, .(codigo, nombre, anio, particion, celda, dist_km, estaciones_en_el_pixel, vertiente, altitud, ONI_anual, fase, marcado)], "F4_pares_estacion_pixel")

Y <- o$coef_est[pares$fila_est, ] %*% t(Phi_e)
X <- o$coef_sat[pares$fila_sat, ] %*% t(Phi_s)
D <- Y - X
ok <- rep(TRUE, nrow(pares))   # todos los años entran (D-47)
marc <- pares$marcado

banda_diferencia <- function(D, grupo) {
  medias <- rowsum(D, grupo) / as.vector(table(grupo)[rownames(rowsum(D, grupo))])
  n <- nrow(medias)
  list(d = colMeans(medias), ee = apply(medias, 2, sd) / sqrt(n), n = n, ee_ingenuo = apply(D, 2, sd) / sqrt(nrow(D)))
}
intervalos <- function(inf, sup) {
  signo <- ifelse(inf > 0, "estación > satélite", ifelse(sup < 0, "satélite > estación", NA))
  r <- rle(signo); fin <- cumsum(r$lengths); ini <- fin - r$lengths + 1
  data.table(desde = fecha_dia(ini), hasta = fecha_dia(fin), dias = r$lengths, signo = r$values)[!is.na(signo)]
}

b <- banda_diferencia(D[ok, ], pares$codigo[ok])
cat(sprintf("EE por estación / EE tratando cada curva como independiente: %.2f en promedio\n", mean(b$ee / b$ee_ingenuo)))
int_d <- intervalos(b$d - 2 * b$ee, b$d + 2 * b$ee)
print(tabla_f(int_d, "F4_intervalos_diferencia"))
cat(sprintf("d(t) medio: %.2f mm/día (%.0f mm/año); días en que la banda excluye el 0: %d de 365\n",
            mean(b$d), mean(b$d) * 365, sum(int_d$dias)))
b_def <- banda_diferencia(D[marc, ], pares$codigo[marc]); b_resto <- banda_diferencia(D[!marc, ], pares$codigo[!marc])
cat(sprintf("En 1985, 1988 y 2000: d(t) medio %.2f mm/día (%.0f mm/año) | en el resto: %.2f mm/día (%.0f mm/año)\n", mean(b_def$d), mean(b_def$d) * 365, mean(b_resto$d), mean(b_resto$d) * 365))
# Contraste de la banda ±2 EE (rúbrica 3): (a) percentiles de un bootstrap de estaciones, que
# no supone normalidad de d̄(t), y (b) la misma banda con el PÍXEL como unidad, porque varias
# estaciones comparten píxel y su diferencia no es independiente.
B_BOOT_D <- 2000L
medias_est <- rowsum(D, pares$codigo) / as.vector(table(pares$codigo)[rownames(rowsum(D, pares$codigo))])
set.seed(20261123)
boot_d <- replicate(B_BOOT_D, colMeans(medias_est[sample(nrow(medias_est), replace = TRUE), ]))
d_lo <- apply(boot_d, 1, quantile, .025); d_hi <- apply(boot_d, 1, quantile, .975)
b_pix <- banda_diferencia(D, pares$celda)
excluye_0 <- function(inf, sup) (inf > 0) | (sup < 0)
bandas_d <- data.table(banda = c(sprintf("±2 EE entre estaciones (%d), rúbrica 3", b$n), sprintf("Bootstrap de estaciones, percentiles 95 %% (%d)", B_BOOT_D),
                                 sprintf("±2 EE entre píxeles (%d)", b_pix$n)),
                       ancho_medio_mm_dia = round(c(mean(4 * b$ee), mean(d_hi - d_lo), mean(4 * b_pix$ee)), 3),
                       dias_excluye_0 = c(sum(excluye_0(b$d - 2 * b$ee, b$d + 2 * b$ee)), sum(excluye_0(d_lo, d_hi)), sum(excluye_0(b_pix$d - 2 * b_pix$ee, b_pix$d + 2 * b_pix$ee))),
                       dias_mismo_veredicto = c(365L, sum(excluye_0(d_lo, d_hi) == excluye_0(b$d - 2 * b$ee, b$d + 2 * b$ee)),
                                                sum(excluye_0(b_pix$d - 2 * b_pix$ee, b_pix$d + 2 * b_pix$ee) == excluye_0(b$d - 2 * b$ee, b$d + 2 * b$ee))))
print(tabla_f(bandas_d, "F4_diferencia_bandas_contraste"))
print(tabla_f(intervalos(b_pix$d - 2 * b_pix$ee, b_pix$d + 2 * b_pix$ee), "F4_intervalos_diferencia_por_pixel"))
tabla_f(data.table(dia = T_DIA, d = b$d, ee = b$ee, inferior = b$d - 2 * b$ee, superior = b$d + 2 * b$ee, inferior_bootstrap = d_lo, superior_bootstrap = d_hi,
                   d_por_pixel = b_pix$d, inferior_pixel = b_pix$d - 2 * b_pix$ee, superior_pixel = b_pix$d + 2 * b_pix$ee,
                   d_1985_1988_2000 = b_def$d, d_resto = b_resto$d), "F4_diferencia_por_dia")

grafico_base("F4_diferencia", function() {
  par(mfrow = c(1, 2))
  plot(T_DIA, colMeans(Y[ok, ]), type = "l", lwd = 2, col = 4, ylim = range(c(0, colMeans(Y[ok, ]), colMeans(X[ok, ]))),
       xlab = "Día del año", ylab = "mm/día", main = "Media de los pares: estación y su píxel")
  lines(T_DIA, colMeans(X[ok, ]), lwd = 2, col = "darkorange")
  legend("bottomleft", c("Estaciones", "Satélite (su píxel)"), col = c(4, "darkorange"), lwd = 2, bty = "n")
  plot(T_DIA, b$d, type = "n", ylim = range(c(b$d - 2 * b$ee, b$d + 2 * b$ee, d_lo, d_hi, b_pix$d - 2 * b_pix$ee, b_pix$d + 2 * b_pix$ee, b_def$d, 0)), xlab = "Día del año", ylab = "mm/día",
       main = expression(paste(bar(d), "(t) = estación - satélite  ±2 EE(t)")))
  polygon(c(T_DIA, rev(T_DIA)), c(b$d - 2 * b$ee, rev(b$d + 2 * b$ee)), col = adjustcolor("steelblue", .3), border = NA)
  lines(T_DIA, b$d, lwd = 2); lines(T_DIA, b_def$d, lty = 2, col = "firebrick"); abline(h = 0, lty = 3)
  lines(T_DIA, d_lo, lty = 3, col = "navy"); lines(T_DIA, d_hi, lty = 3, col = "navy")
  lines(T_DIA, b_pix$d - 2 * b_pix$ee, lty = 4, col = "darkgreen"); lines(T_DIA, b_pix$d + 2 * b_pix$ee, lty = 4, col = "darkgreen")
  excl <- (b$d - 2 * b$ee > 0) | (b$d + 2 * b$ee < 0); rug(T_DIA[excl], col = "firebrick")
  legend("bottomleft", c("Todos los años (banda ±2 EE por estación)", "Bootstrap de estaciones 95 %", "±2 EE por píxel", "Solo 1985, 1988 y 2000", "Días en que la banda excluye el 0"),
         col = c(1, "navy", "darkgreen", "firebrick", "firebrick"), lty = c(1, 3, 4, 2, 1), lwd = c(2, 1, 1, 1, 1), bty = "n", cex = .75)
}, ancho = 13, alto = 4.8)

titulo("4.1b Serie diaria de un año y ciclo anual de cada estación frente a su píxel")
# Se muestra la lluvia puntual en mm/día, sin agregar en totales: el dato diario de la
# estación, el dato pentadal de su píxel y las dos curvas suavizadas de un año, y al lado el
# ciclo anual medio de la estación y de su píxel con los años individuales de fondo. Se toma
# una estación de cada zona (costa Pacífica bajo 500 m, fondo del valle y ladera alta del Cauca): la de más
# años; el año es el más reciente con los 365 días observados.
zonas_ej <- list(`Pacífico, costa` = pares$vertiente == "Pacífico" & pares$altitud < 500, `Cauca, fondo del valle` = pares$vertiente != "Pacífico" & pares$banda_altitud == "< 1200 m",
                 `Cauca, ladera alta` = pares$vertiente != "Pacífico" & pares$banda_altitud == "> 2000 m")
n_completos <- o$id_est$n_obs[pares$fila_est] == 365L
ejemplos <- rbindlist(lapply(names(zonas_ej), function(z) {
  k <- which(zonas_ej[[z]] & !is.na(zonas_ej[[z]]))
  cand <- data.table(pares[k, .(codigo, nombre, altitud, celda)], completo = n_completos[k])[
    , .(anios = .N, con_completo = any(completo)), by = .(codigo, nombre, altitud, celda)][con_completo == TRUE][order(-anios, codigo)][1]
  fil <- which(pares$codigo == cand$codigo); anio_ej <- max(pares$anio[fil][n_completos[fil]])
  cand[, .(zona = z, codigo, nombre, altitud, celda, anios, anio_ejemplo = anio_ej, fila = fil[pares$anio[fil] == anio_ej],
           P_estacion_mm_dia = round(mean(Y[fil, ]), 2), P_satelite_mm_dia = round(mean(X[fil, ]), 2),
           r_ciclo_anual = round(cor(colMeans(Y[fil, , drop = FALSE]), colMeans(X[fil, , drop = FALSE])), 2))]
}))
print(tabla_f(ejemplos[, !"fila"], "F4_ejemplos_estacion_pixel"))
grafico_base("F4_serie_diaria_y_ciclo_anual", function() {
  par(mfrow = c(nrow(ejemplos), 2), mar = c(4, 4, 3, 1))
  for (k in seq_len(nrow(ejemplos))) with(ejemplos[k], {
    fe <- pares$fila_est[fila]; fs <- pares$fila_sat[fila]; y_obs <- Y_dia[fe, ]; x_obs <- X_pen[fs, ]
    plot(T_DIA, y_obs, type = "h", col = "grey60", ylim = c(0, max(y_obs, x_obs, Y[fila, ], X[fila, ], na.rm = TRUE)), xlab = "Día del año", ylab = "mm/día",
         main = sprintf("%s (%s, %.0f m) · %d", nombre, zona, altitud, anio_ejemplo), cex.main = .95)
    points(T_SAT, x_obs, pch = 19, cex = .5, col = "darkorange")
    lines(T_DIA, Y[fila, ], lwd = 2, col = 4); lines(T_DIA, X[fila, ], lwd = 2, col = "darkorange")
    if (k == 1) legend("topright", c("Estación: dato diario", "Satélite: dato pentadal", "Curva de la estación", "Curva del satélite"),
                       col = c("grey60", "darkorange", 4, "darkorange"), lty = c(1, NA, 1, 1), pch = c(NA, 19, NA, NA), lwd = c(1, NA, 2, 2), bty = "n", cex = .75)
    fil <- which(pares$codigo == codigo)
    matplot(T_DIA, t(Y[fil, , drop = FALSE]), type = "l", lty = 1, col = adjustcolor(4, .15), ylim = c(0, max(Y[fil, ], X[fil, ])), xlab = "Día del año", ylab = "mm/día",
            main = sprintf("Ciclo anual medio, %d años (r = %.2f)", anios, r_ciclo_anual), cex.main = .95)
    matlines(T_DIA, t(X[fil, , drop = FALSE]), lty = 1, col = adjustcolor("darkorange", .15))
    lines(T_DIA, colMeans(Y[fil, , drop = FALSE]), lwd = 3, col = 4); lines(T_DIA, colMeans(X[fil, , drop = FALSE]), lwd = 3, col = "darkorange")
    if (k == 1) legend("topright", c("Estación (media de sus años)", "Su píxel (media de los mismos años)"), col = c(4, "darkorange"), lwd = 3, bty = "n", cex = .75)
  })
}, ancho = 13, alto = 3.6 * 3)

titulo("4.2 Correlación cruzada y covarianzas")
Ce <- o$coef_est[pares$fila_est[ok], ]; Cs <- o$coef_sat[pares$fila_sat[ok], ]
Cec <- sweep(Ce, 2, colMeans(Ce)); Csc <- sweep(Cs, 2, colMeans(Cs)); n <- nrow(Ce)
Cov_YY <- Phi_e %*% (crossprod(Cec) / (n - 1)) %*% t(Phi_e)
Cov_XX <- Phi_s %*% (crossprod(Csc) / (n - 1)) %*% t(Phi_s)
Cov_YX <- Phi_e %*% (crossprod(Cec, Csc) / (n - 1)) %*% t(Phi_s)
R_YX <- Cov_YX / sqrt(outer(diag(Cov_YY), diag(Cov_XX)))
cat("Covarianza cruzada por coeficientes vs cov() de las curvas evaluadas:", format(max(abs(Cov_YX - cov(Y[ok, ], X[ok, ]))), digits = 2), "\n")
diag_r <- diag(R_YX)
cat(sprintf("Diagonal r_XY(t,t): media %.2f, máxima %.2f el %s, mínima %.2f el %s\n",
            mean(diag_r), max(diag_r), fecha_dia(which.max(diag_r)), min(diag_r), fecha_dia(which.min(diag_r))))
mejor_s <- apply(R_YX, 1, which.max)
cat(sprintf("Para cada día t de la estación, el día s del satélite más correlacionado está a %+.1f días en mediana\n", median(mejor_s - T_DIA)))

anom <- function(M, g) M - (rowsum(M, g) / as.vector(table(g)[rownames(rowsum(M, g))]))[as.character(g), ]
Ya <- anom(Y[ok, ], pares$codigo[ok]); Xa <- anom(X[ok, ], pares$codigo[ok])
R_anom <- diag(cov(Ya, Xa)) / sqrt(apply(Ya, 2, var) * apply(Xa, 2, var))
cat(sprintf("Diagonal con anomalías (restando la climatología de cada estación): media %.2f, máxima %.2f el %s, mínima %.2f el %s\n",
            mean(R_anom), max(R_anom), fecha_dia(which.max(R_anom)), min(R_anom), fecha_dia(which.min(R_anom))))
tabla_f(data.table(dia = T_DIA, r_diagonal = diag_r, r_diagonal_anomalias = R_anom, var_estacion = diag(Cov_YY), var_satelite = diag(Cov_XX)), "F4_diagonal_y_varianzas")

cov_comp <- data.table(fuente = c("estación", "satélite"), varianza_total = c(sum(diag(Cov_YY)), sum(diag(Cov_XX))),
                       DE_media = c(mean(sqrt(diag(Cov_YY))), mean(sqrt(diag(Cov_XX)))),
                       cor_abr_oct = c(cov2cor(Cov_YY)[105, 290], cov2cor(Cov_XX)[105, 290]),
                       cor_abr_jul = c(cov2cor(Cov_YY)[105, 200], cov2cor(Cov_XX)[105, 200]),
                       var_explicada_1a_componente = c(eigen(Cov_YY, symmetric = TRUE, only.values = TRUE)$values[1] / sum(diag(Cov_YY)),
                                                       eigen(Cov_XX, symmetric = TRUE, only.values = TRUE)$values[1] / sum(diag(Cov_XX))))
print(tabla_f(cov_comp, "F4_covarianzas_comparadas"))
cat(sprintf("Cociente de varianzas satélite / estación por día: mediana %.2f (rango %.2f a %.2f)\n",
            median(diag(Cov_XX) / diag(Cov_YY)), min(diag(Cov_XX) / diag(Cov_YY)), max(diag(Cov_XX) / diag(Cov_YY))))

# ¿Las diferencias entre fuentes salen de usar K distinto (8 y 11)? Se repite la
# comparación con el satélite construido con el mismo K de las estaciones.
X_k <- construir(X_pen[pares$fila_sat, ], T_SAT, o$K_est)$coef %*% t(Phi_e); b_k <- banda_diferencia(Y - X_k, pares$codigo)
excluye <- function(bb) sum((bb$d - 2 * bb$ee > 0) | (bb$d + 2 * bb$ee < 0))
mismo_K <- data.table(K_estacion = o$K_est, K_satelite = c(o$K_sat, o$K_est),
                      d_medio_mm_anio = round(c(mean(b$d), mean(b_k$d)) * 365), dias_banda_sin_cero = c(excluye(b), excluye(b_k)),
                      dias_con_el_mismo_veredicto = c(365L, sum(((b$d - 2 * b$ee > 0) - (b$d + 2 * b$ee < 0)) == ((b_k$d - 2 * b_k$ee > 0) - (b_k$d + 2 * b_k$ee < 0)))),
                      r_diagonal_media = round(c(mean(diag_r), mean(diag(cor(Y, X_k)))), 3),
                      cociente_varianzas_mediana = round(c(median(diag(Cov_XX) / diag(Cov_YY)), median(apply(X_k, 2, var) / apply(Y, 2, var))), 3))
print(tabla_f(mismo_K, "F4_sensibilidad_mismo_K"))

grafico_base("F4_correlacion_cruzada", function() {
  par(mfrow = c(1, 1)); mapa_calor(T_DIA, R_YX, "Correlación cruzada r(s,t): satélite en s, estación en t", c(-1, 1), xlab = "Día t (estación)", ylab = "Día s (satélite)")
}, ancho = 7, alto = 5.5)
grafico_base("F4_diagonal_correlacion", function() {
  plot(T_DIA, diag_r, type = "l", lwd = 2, ylim = c(0, 1), xlab = "Día del año", ylab = "r(t, t)", main = "Diagonal de la correlación cruzada")
  lines(T_DIA, R_anom, lwd = 2, col = "firebrick")
  legend("bottomleft", c("Curvas", "Anomalías (sin la climatología de cada estación)"), col = c(1, "firebrick"), lwd = 2, bty = "n")
}, ancho = 9, alto = 4.5)
lim_cov <- range(c(Cov_YY, Cov_XX))
grafico_base("F4_covarianzas", function() {
  layout(matrix(1:3, 1)); par(mar = c(4, 4, 3, 1))
  image(T_DIA, T_DIA, Cov_YY, zlim = lim_cov, col = hcl.colors(64, "YlOrRd", rev = TRUE), main = "Covarianza: estaciones", xlab = "Día s", ylab = "Día t")
  image(T_DIA, T_DIA, Cov_XX, zlim = lim_cov, col = hcl.colors(64, "YlOrRd", rev = TRUE), main = "Covarianza: su píxel", xlab = "Día s", ylab = "Día t")
  lim_d <- max(abs(Cov_XX - Cov_YY)) * c(-1, 1)
  image(T_DIA, T_DIA, Cov_XX - Cov_YY, zlim = lim_d, col = hcl.colors(64, "Blue-Red 3"), main = "Satélite - estación", xlab = "Día s", ylab = "Día t")
  layout(1)
}, ancho = 15, alto = 4.8)

titulo("4.3 Por vertiente, altitud y fase ENSO")
por_grupo <- function(var) {
  g <- as.character(pares[[var]]); k <- ok & !is.na(g)
  rbindlist(lapply(intersect(niveles_de(pares[[var]]), unique(g[k])), function(v) {
    i <- k & g == v; bb <- banda_diferencia(D[i, , drop = FALSE], pares$codigo[i])
    data.table(agrupacion = var, grupo = v, pares = sum(i), estaciones = bb$n,
               P_estacion = mean(Y[i, ]) * 365, P_satelite = mean(X[i, ]) * 365, sesgo_mm_anio = mean(bb$d) * 365,
               cociente = mean(X[i, ]) / mean(Y[i, ]), dias_banda_sin_cero = sum((bb$d - 2 * bb$ee > 0) | (bb$d + 2 * bb$ee < 0)),
               r_diagonal_media = mean(diag(cor(Y[i, ], X[i, ]))))
  }))
}
grupos <- rbind(por_grupo("vertiente"), por_grupo("banda_altitud"), por_grupo("fase"))
grupos[, c("P_estacion", "P_satelite", "sesgo_mm_anio") := lapply(.SD, round), .SDcols = c("P_estacion", "P_satelite", "sesgo_mm_anio")]
grupos[, c("cociente", "r_diagonal_media") := lapply(.SD, round, 2), .SDcols = c("cociente", "r_diagonal_media")]
print(tabla_f(grupos, "F4_comparacion_por_grupo"))

pares[, `:=`(P_est = rowMeans(Y) * 365, P_sat = rowMeans(X) * 365)]
por_est <- pares[ok, .(P_est = mean(P_est), P_sat = mean(P_sat), altitud = altitud[1], vertiente = vertiente[1], nombre = nombre[1]), by = codigo]
por_est[, cociente := P_sat / P_est]
ajuste_alt <- lm(log(cociente) ~ altitud, data = por_est[vertiente != "Pacífico"])
cat(sprintf("Estaciones del Cauca: log(satélite/estación) cambia %.3f por cada 1 000 m de altitud (p = %.3g); r = %.2f\n",
            1000 * coef(ajuste_alt)[2], summary(ajuste_alt)$coefficients[2, 4], cor(por_est[vertiente != "Pacífico" & !is.na(altitud), .(altitud, log(cociente))])[1, 2]))
cat(sprintf("Lluvia anual estación vs su píxel (promedio por estación): r = %.2f; en el Cauca r = %.2f\n",
            cor(por_est$P_est, por_est$P_sat), por_est[vertiente != "Pacífico", cor(P_est, P_sat)]))
tabla_f(por_est[order(cociente)], "F4_sesgo_por_estacion")
print(rbind(head(por_est[order(cociente)], 4), tail(por_est[order(cociente)], 4)))

grafico_f(ggplot(por_est, aes(P_est, P_sat, colour = vertiente)) + geom_abline(linetype = 2) + geom_point() +
            scale_x_log10() + scale_y_log10() + labs(x = "Estación (mm/año, log)", y = "Su píxel (mm/año, log)", colour = NULL,
                                                   title = "Lluvia anual media: estación frente a su píxel (1981-2003)"),
          "F4_anual_estacion_vs_pixel", ancho = 8, alto = 6)
grafico_f(ggplot(por_est[!is.na(altitud)], aes(altitud, cociente, colour = vertiente)) + geom_hline(yintercept = 1, linetype = 2) +
            geom_point() + geom_smooth(method = "lm", formula = y ~ x, se = FALSE) + scale_y_log10() +
            labs(x = "Altitud de la estación (m)", y = "Satélite / estación (log)", colour = NULL, title = "Sesgo relativo del satélite según la altitud"),
          "F4_sesgo_altitud", ancho = 8, alto = 5)

curvas_d <- rbindlist(lapply(c("vertiente", "banda_altitud", "fase"), function(var) {
  g <- as.character(pares[[var]]); k <- ok & !is.na(g)
  rbindlist(lapply(intersect(niveles_de(pares[[var]]), unique(g[k])), function(v) { i <- k & g == v; bb <- banda_diferencia(D[i, , drop = FALSE], pares$codigo[i])
    data.table(agrupacion = var, grupo = v, dia = T_DIA, d = bb$d, inf = bb$d - 2 * bb$ee, sup = bb$d + 2 * bb$ee) }))
}))
grafico_f(ggplot(curvas_d, aes(dia, d, colour = grupo, fill = grupo)) + geom_hline(yintercept = 0, linetype = 3) +
            geom_ribbon(aes(ymin = inf, ymax = sup), alpha = .15, colour = NA) + geom_line(linewidth = .8) +
            facet_wrap(~agrupacion, ncol = 3, labeller = as_labeller(c(vertiente = "Vertiente", banda_altitud = "Altitud", fase = "Fase ENSO"))) + labs(x = "Día del año", y = "d(t), mm/día", colour = NULL, fill = NULL,
                                                      title = "Diferencia estación - satélite por grupo (±2 EE por estación)"),
          "F4_diferencia_por_grupo", ancho = 14, alto = 4.8)

titulo("4.4 Regresión concurrente: ¿en qué época del año el satélite explica a la estación? (clase, ejercicio 3)")
# Script_Regresion_Funcional.R (modelo concurrente y ejercicio 3 con los datos de precipitación):
#     y_i(t) = b0(t) + b1(t) x_i(t) + e_i(t)                                     (modelo A, el de clase)
#     y_i(t) = b0(t) + b1(t) x_i(t) + g(t) altitud_i + d(t) ONI_i + e_i(t)       (modelo B)
# y_i(t) es la curva de la estación y x_i(t) la de su píxel en el mismo año, las dos en raíz de
# mm/día (la escala de la sección 5). Cada día t de la estación se explica solo con el mismo
# día t del satélite: b1(t) dice cuánto de la señal del satélite pasa a la estación ese día y
# R²(t) en qué época del año el satélite la explica mejor. El satélite es la covariable y la
# estación la respuesta: no se supone que sean dos mediciones del mismo individuo. Como el
# satélite trae error de medida, b1(t) < 1 se lee como "conviene encoger la señal del satélite
# hacia la media", no como un error de escala; la igualdad de las medias la responde d(t) (4.1).
# Cada b_a(t) vive en una base B-spline de K_b funciones con penalización lambda ∫ b_a''(t)² dt.
# K_b y lambda se eligen con validación (1996-2000) y los coeficientes se ajustan solo con
# entrenamiento. Bandas: sándwich ±2 EE (la de clase, con Σ_E entre días pero curvas
# independientes) y bootstrap de estaciones y de píxeles (los años de una estación y las
# estaciones de un mismo píxel no son independientes). Se contrasta con fda::fRegress.
concurrente <- local({
  W_T <- c(.5, rep(1, 363), .5)
  MESES <- data.table(mes = c("ene", "feb", "mar", "abr", "may", "jun", "jul", "ago", "sep", "oct", "nov", "dic"),
                      ini = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335))
  Y_r <- a_raiz(o$coef_est[pares$fila_est, ], Phi_e)$coef %*% t(Phi_e)
  X_r <- a_raiz(o$coef_sat[pares$fila_sat, ], Phi_s)$coef %*% t(Phi_s)
  Y_mm <- pmax(Y, 0); X_mm <- pmax(X, 0)
  usar <- !is.na(pares$altitud)
  i_ent <- which(usar & pares$particion == "entrenamiento"); i_val <- which(usar & pares$particion == "validación"); i_pru <- which(usar & pares$particion == "prueba")
  alt_m <- mean(pares$altitud[i_ent]); alt_s <- sd(pares$altitud[i_ent])
  ALT <- (pares$altitud - alt_m) / alt_s; ONI <- pares$ONI_anual
  cat(sprintf("Curvas: entrenamiento %d, validación %d, prueba %d\n", length(i_ent), length(i_val), length(i_pru)))

  # Mínimos cuadrados penalizados por bloques: G_ab = Φᵀ diag(w Σ_i c_ia c_ib) Φ y r_a = Φᵀ (w Σ_i c_ia y_i),
  # con c_i0 = 1, c_i1(t) = x_i(t) y, en el modelo B, c_i2 = altitud_i y c_i3 = ONI_i.
  covariables <- function(i, modelo) {
    n <- length(i); L <- list(const = matrix(1, n, 365), sat = X_r[i, , drop = FALSE])
    if (modelo == "B") { L$altitud <- matrix(ALT[i], n, 365); L$ONI <- matrix(ONI[i], n, 365) }
    L
  }
  preparar <- function(i, modelo, K_b) {
    Cv <- covariables(i, modelo); p <- length(Cv); Phi_b <- base_bspline(T_DIA, K_b)
    G <- matrix(0, p * K_b, p * K_b); r <- numeric(p * K_b)
    for (a in 1:p) {
      r[(a - 1) * K_b + 1:K_b] <- crossprod(Phi_b, W_T * colSums(Cv[[a]] * Y_r[i, , drop = FALSE]))
      for (b in a:p) {
        blq <- crossprod(Phi_b, Phi_b * (W_T * colSums(Cv[[a]] * Cv[[b]])))
        G[(a - 1) * K_b + 1:K_b, (b - 1) * K_b + 1:K_b] <- blq
        if (b != a) G[(b - 1) * K_b + 1:K_b, (a - 1) * K_b + 1:K_b] <- t(blq)
      }
    }
    list(G = G, r = r, p = p, K_b = K_b, Phi_b = Phi_b, R = kronecker(diag(p), gram_derivada(K_b, 2L)), nombres = names(Cv))
  }
  resolver <- function(pre, lambda) solve(pre$G + lambda * pre$R, pre$r)
  coef_t <- function(pre, theta) { B <- sapply(1:pre$p, function(a) pre$Phi_b %*% theta[(a - 1) * pre$K_b + 1:pre$K_b]); colnames(B) <- pre$nombres; B }
  predecir <- function(B, i, modelo) { Cv <- covariables(i, modelo); P <- 0; for (a in seq_along(Cv)) P <- P + sweep(Cv[[a]], 2, B[, a], "*"); P }
  a_mm <- function(R) pmax(R, 0)^2
  rmse_l2 <- function(P, i) sqrt(rowMeans((P - Y_mm[i, , drop = FALSE])^2))

  LAMBDAS <- c(0, 10^seq(-2, 12, by = .5)); KBS <- c(8L, 12L, 16L)
  seleccion <- rbindlist(lapply(c("A", "B"), function(mod) rbindlist(lapply(KBS, function(K_b) {
    pre <- preparar(i_ent, mod, K_b)
    rbindlist(lapply(LAMBDAS, function(l) {
      Bl <- coef_t(pre, resolver(pre, l))
      data.table(modelo = mod, K_b = K_b, lambda = l, rmse_entrenamiento = mean(rmse_l2(a_mm(predecir(Bl, i_ent, mod)), i_ent)),
                 rmse_validacion = mean(rmse_l2(a_mm(predecir(Bl, i_val, mod)), i_val)))
    }))
  }))))
  tabla_f(seleccion, "F4_concurrente_seleccion_Kb_lambda")
  # Regla (la misma de M en 5.3): el menor K_b cuyo mejor RMSE de validación queda a menos de 0,5 % del mínimo.
  elegir <- function(mod) {
    mejor <- seleccion[modelo == mod, .SD[which.min(rmse_validacion)], by = K_b]
    K_sel <- mejor[rmse_validacion <= min(rmse_validacion) * 1.005, min(K_b)]
    list(K_b = K_sel, lambda = mejor[K_b == K_sel, lambda], mejor = mejor)
  }
  eleccion <- list(A = elegir("A"), B = elegir("B"))
  for (mod in c("A", "B")) cat(sprintf("Modelo %s: K_b = %d, lambda = %s (RMSE de validación %.3f mm/día)\n", mod, eleccion[[mod]]$K_b,
                                       format(eleccion[[mod]]$lambda, scientific = TRUE), eleccion[[mod]]$mejor[K_b == eleccion[[mod]]$K_b, rmse_validacion]))
  grafico_base("F4_concurrente_seleccion_lambda", function() {
    par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
    for (mod in c("A", "B")) {
      s <- seleccion[modelo == mod & lambda > 0]
      plot(NA, xlim = range(log10(s$lambda)), ylim = range(s$rmse_validacion), xlab = expression(log[10](lambda)),
           ylab = "RMSE L2 medio en validación (mm/día)", main = paste("Modelo concurrente", mod))
      for (j in seq_along(KBS)) lines(log10(s[K_b == KBS[j], lambda]), s[K_b == KBS[j], rmse_validacion], col = j + 1, lwd = 2)
      abline(v = log10(eleccion[[mod]]$lambda), lty = 2); legend("topleft", paste("K_b =", KBS), col = seq_along(KBS) + 1, lwd = 2, bty = "n")
    }
  })

  ajuste <- lapply(c(A = "A", B = "B"), function(mod) {
    pre <- preparar(i_ent, mod, eleccion[[mod]]$K_b); th <- resolver(pre, eleccion[[mod]]$lambda)
    list(pre = pre, theta = th, B = coef_t(pre, th), lambda = eleccion[[mod]]$lambda)
  })
  # Banda sándwich de clase: Cov(θ̂) = H⁻¹ [Σ_i Z_iᵀ W Σ_E W Z_i] H⁻¹, H = G + lambda R, Σ_E con los residuos de entrenamiento.
  ee_sandwich <- function(mod) {
    a <- ajuste[[mod]]; pre <- a$pre; Cv <- covariables(i_ent, mod); p <- pre$p; K_b <- pre$K_b
    Res <- Y_r[i_ent, ] - predecir(a$B, i_ent, mod); SigE <- crossprod(sweep(Res, 2, colMeans(Res))) / (length(i_ent) - 1)
    Hinv <- solve(pre$G + a$lambda * pre$R); WSW <- (W_T %o% W_T) * SigE; Mz <- matrix(0, p * K_b, p * K_b)
    for (k in seq_along(i_ent)) { Z <- do.call(cbind, lapply(1:p, function(a_) Cv[[a_]][k, ] * pre$Phi_b)); Mz <- Mz + crossprod(Z, WSW %*% Z) }
    V <- Hinv %*% Mz %*% Hinv
    sapply(1:p, function(a_) { idx <- (a_ - 1) * K_b + 1:K_b; sqrt(rowSums((pre$Phi_b %*% V[idx, idx]) * pre$Phi_b)) })
  }
  for (mod in c("A", "B")) ajuste[[mod]]$EE <- ee_sandwich(mod)
  tabla_f(rbindlist(lapply(c("A", "B"), function(mod) rbindlist(lapply(seq_len(ajuste[[mod]]$pre$p), function(a)
    data.table(modelo = mod, coeficiente = ajuste[[mod]]$pre$nombres[a], dia = T_DIA, valor = ajuste[[mod]]$B[, a], EE_sandwich = ajuste[[mod]]$EE[, a]))))),
    "F4_concurrente_coeficientes")

  # Bootstrap por bloques del modelo A con el mismo K_b y lambda: se remuestrean estaciones de
  # entrenamiento con todos sus años y, como sensibilidad, píxeles con todas sus estaciones y años.
  B_BOOT <- 300L
  remuestreo <- function(grupo, semilla) {
    bloques <- split(i_ent, grupo[i_ent]); set.seed(semilla)
    bb <- replicate(B_BOOT, { f <- unlist(bloques[sample(length(bloques), replace = TRUE)], use.names = FALSE)
      pre <- preparar(f, "A", eleccion$A$K_b); coef_t(pre, resolver(pre, eleccion$A$lambda)) })
    list(lo = apply(bb, c(1, 2), quantile, .025), hi = apply(bb, c(1, 2), quantile, .975), bloques = length(bloques))
  }
  bt_est <- remuestreo(pares$codigo, 20261121); bt_pix <- remuestreo(pares$celda, 20261122)
  bA <- ajuste$A$B; eA <- ajuste$A$EE
  bandas <- rbindlist(lapply(1:2, function(a) data.table(
    coeficiente = c("b0(t)", "b1(t)")[a], banda = c("Sándwich ±2 EE (clase)", sprintf("Bootstrap de estaciones (%d)", bt_est$bloques), sprintf("Bootstrap de píxeles (%d)", bt_pix$bloques)),
    ancho_medio = round(c(mean(4 * eA[, a]), mean(bt_est$hi[, a] - bt_est$lo[, a]), mean(bt_pix$hi[, a] - bt_pix$lo[, a])), 3),
    dias_excluye_0 = c(sum(abs(bA[, a]) > 2 * eA[, a]), sum(bt_est$lo[, a] > 0 | bt_est$hi[, a] < 0), sum(bt_pix$lo[, a] > 0 | bt_pix$hi[, a] < 0)),
    dias_excluye_1 = if (a == 2) c(sum(abs(bA[, 2] - 1) > 2 * eA[, 2]), sum(bt_est$lo[, 2] > 1 | bt_est$hi[, 2] < 1), sum(bt_pix$lo[, 2] > 1 | bt_pix$hi[, 2] < 1)) else NA_integer_)))
  print(tabla_f(bandas, "F4_concurrente_bandas"))
  b1_bajo_1 <- intervalos(bt_est$lo[, 2] - 1, bt_est$hi[, 2] - 1)
  b1_bajo_1[, signo := fifelse(signo == "estación > satélite", "b1(t) > 1", "b1(t) < 1")]
  print(tabla_f(b1_bajo_1, "F4_concurrente_intervalos_b1"))

  R2_t <- function(P, i) 1 - colSums((Y_r[i, , drop = FALSE] - P)^2) / colSums(sweep(Y_r[i, , drop = FALSE], 2, colMeans(Y_r[i_ent, ]))^2)
  R2_tab <- data.table(dia = T_DIA, R2_A_prueba = R2_t(predecir(bA, i_pru, "A"), i_pru), R2_B_prueba = R2_t(predecir(ajuste$B$B, i_pru, "B"), i_pru),
                       R2_A_validacion = R2_t(predecir(bA, i_val, "A"), i_val), R2_satelite_prueba = R2_t(X_r[i_pru, ], i_pru))
  tabla_f(R2_tab, "F4_concurrente_R2_por_dia")
  mes_de <- cut(T_DIA, c(MESES$ini, 366), labels = MESES$mes, right = FALSE)
  por_mes_de <- function(v, f = mean) as.vector(tapply(v, mes_de, f))
  por_mes <- data.table(mes = MESES$mes, R2_A_prueba = por_mes_de(R2_tab$R2_A_prueba), R2_A_validacion = por_mes_de(R2_tab$R2_A_validacion),
                        R2_satelite_prueba = por_mes_de(R2_tab$R2_satelite_prueba), b1_medio = por_mes_de(bA[, 2]),
                        dias_b1_bajo_1 = por_mes_de(bt_est$hi[, 2] < 1, sum), dias_b1_sobre_1 = por_mes_de(bt_est$lo[, 2] > 1, sum))
  por_mes[, c("R2_A_prueba", "R2_A_validacion", "R2_satelite_prueba", "b1_medio") := lapply(.SD, round, 3), .SDcols = c("R2_A_prueba", "R2_A_validacion", "R2_satelite_prueba", "b1_medio")]
  print(tabla_f(por_mes, "F4_concurrente_R2_y_b1_por_mes"))
  cat(sprintf("R²(t) del modelo A en prueba: máximo %.2f el %s, mínimo %.2f el %s, promedio %.2f (satélite sin corregir %.2f)\n",
              max(R2_tab$R2_A_prueba), fecha_dia(which.max(R2_tab$R2_A_prueba)), min(R2_tab$R2_A_prueba), fecha_dia(which.min(R2_tab$R2_A_prueba)),
              mean(R2_tab$R2_A_prueba), mean(R2_tab$R2_satelite_prueba)))

  Res_A <- Y_r[i_ent, ] - predecir(bA, i_ent, "A")
  grafico_base("F4_concurrente_estilo_clase", function() {
    par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3, 1)); set.seed(3); m <- sample(i_ent, 150)
    matplot(T_DIA, t(X_r[m, ]), type = "l", lty = 1, col = "grey60", xlab = "Día del año", ylab = "raíz(mm/día)", main = "Satélite: x_i(t)")
    matplot(T_DIA, t(Y_r[m, ]), type = "l", lty = 1, col = "grey60", xlab = "Día del año", ylab = "raíz(mm/día)", main = "Estación: y_i(t)")
    matplot(T_DIA, cbind(R2_tab$R2_A_prueba, R2_tab$R2_B_prueba, R2_tab$R2_satelite_prueba), type = "l", lty = c(1, 1, 2), lwd = 2, col = c(1, 4, "grey50"),
            ylim = c(min(0, unlist(R2_tab[, -1])), 1), xlab = "Día del año", ylab = expression(R^2 * (t)), main = "Varianza explicada por día (prueba 2001-2003)")
    abline(h = 0, lty = 3); legend("bottomleft", c("Modelo A (clase)", "Modelo B (+ altitud y ONI)", "Satélite sin corregir"), col = c(1, 4, "grey50"), lty = c(1, 1, 2), lwd = 2, bty = "n", cex = .8)
    for (a in 1:2) {
      plot(T_DIA, bA[, a], type = "n", ylim = range(c(bt_est$lo[, a], bt_est$hi[, a], bt_pix$lo[, a], bt_pix$hi[, a], bA[, a] - 2 * eA[, a], bA[, a] + 2 * eA[, a], 0)),
           xlab = "Día del año", ylab = "", main = c("Intercepto: b0(t)", "Efecto del satélite: b1(t)")[a])
      polygon(c(T_DIA, rev(T_DIA)), c(bt_est$lo[, a], rev(bt_est$hi[, a])), col = adjustcolor("steelblue", .3), border = NA)
      lines(T_DIA, bt_pix$lo[, a], lty = 3, col = "darkgreen"); lines(T_DIA, bt_pix$hi[, a], lty = 3, col = "darkgreen")
      lines(T_DIA, bA[, a], lwd = 2); lines(T_DIA, bA[, a] + 2 * eA[, a], lty = 2, col = "red"); lines(T_DIA, bA[, a] - 2 * eA[, a], lty = 2, col = "red")
      abline(h = c(0, if (a == 2) 1), lty = 3)
      if (a == 2) legend("bottomleft", c("Bootstrap de estaciones 95 %", "Bootstrap de píxeles 95 %", "Sándwich ±2 EE (clase)"), fill = c(adjustcolor("steelblue", .3), NA, NA),
                         border = NA, lty = c(NA, 3, 2), col = c(NA, "darkgreen", "red"), bty = "n", cex = .8)
    }
    contour(T_DIA, T_DIA, cor(Res_A), nlevels = 10, xlab = "Día s", ylab = "Día t", main = "Correlación de los residuos (entrenamiento)")
  }, ancho = 14, alto = 8.5)
  # Diagnóstico de linealidad en la escala raíz: residuos frente al satélite en cada temporada.
  temp_dia <- cut(T_DIA, c(0, 59, 151, 243, 304, 365), labels = c("dic-feb (seca)", "mar-may (1.a lluviosa)", "jun-ago (seca)", "sep-oct (transición)", "nov-dic (2.a lluviosa)"))
  diag_lin <- rbindlist(lapply(levels(temp_dia), function(tp) {
    d <- which(temp_dia == tp); x <- as.vector(X_r[i_ent, d]); r <- as.vector(Res_A[, d])
    q <- cut(x, unique(quantile(x, 0:5 / 5)), include.lowest = TRUE)
    data.table(temporada = factor(tp, levels(temp_dia)), quintil_satelite = seq_len(nlevels(q)), x_medio = round(as.vector(tapply(x, q, mean)), 2), residuo_medio = round(as.vector(tapply(r, q, mean)), 3),
               residuo_DE = round(as.vector(tapply(r, q, sd)), 3))
  }))
  tabla_f(diag_lin, "F4_concurrente_linealidad")
  print(dcast(diag_lin, temporada ~ quintil_satelite, value.var = "residuo_medio"))
  cat(sprintf("DE de los residuos en entrenamiento: %.2f raíz(mm/día)
", sd(Res_A)))
  grafico_base("F4_concurrente_linealidad", function() {
    par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3, 1)); cols <- c("grey30", "steelblue", "darkorange", "darkgreen", "firebrick")
    plot(NA, xlim = range(diag_lin$x_medio), ylim = range(c(diag_lin$residuo_medio, -.1, .1)), xlab = "Satélite x(t), raíz(mm/día) (media de cada quintil)",
         ylab = "Residuo medio, raíz(mm/día)", main = "Linealidad: residuo medio por quintil del satélite y temporada (entrenamiento)")
    for (k in seq_along(levels(temp_dia))) with(diag_lin[temporada == levels(temp_dia)[k]], lines(x_medio, residuo_medio, type = "b", pch = 19, col = cols[k], lwd = 2))
    abline(h = 0, lty = 3); legend("topright", levels(temp_dia), col = cols, lwd = 2, pch = 19, bty = "n", cex = .8)
  }, ancho = 9, alto = 5)

  # El modelo A año por año (tesis del profesor, figura 9.2): mínimos cuadrados en cada día con las
  # estaciones de ese año. Es descriptivo (usa los 23 años) y no se usa para elegir ni predecir nada.
  anios_c <- sort(unique(pares$anio[usar]))
  por_anio <- lapply(anios_c, function(an) {
    i <- which(usar & pares$anio == an); n <- length(i); x <- X_r[i, , drop = FALSE]; y <- Y_r[i, , drop = FALSE]
    xc <- sweep(x, 2, colMeans(x)); yc <- sweep(y, 2, colMeans(y))
    b1 <- colSums(xc * yc) / colSums(xc^2); b0 <- colMeans(y) - b1 * colMeans(x); res <- y - (rep(1, n) %o% b0 + sweep(x, 2, b1, "*"))
    list(anio = an, n = n, b1 = b1, ee = sqrt(colSums(res^2) / (n - 2) / colSums(xc^2)), R2 = 1 - colSums(res^2) / colSums(yc^2), i = i)
  })
  b1_anio <- rbindlist(lapply(por_anio, function(p) data.table(anio = p$anio, estaciones = p$n, b1_medio = round(mean(p$b1), 3),
    dias_b1_excluye_0 = sum(p$b1 - 2 * p$ee > 0 | p$b1 + 2 * p$ee < 0), dias_b1_excluye_1 = sum(p$b1 - 2 * p$ee > 1 | p$b1 + 2 * p$ee < 1), R2_medio = round(mean(p$R2), 3))))
  b1_anio <- unique(oni[, .(anio, fase, ONI_anual)])[b1_anio, on = "anio"]
  print(tabla_f(b1_anio, "F4_concurrente_por_anio"))
  ct_b1 <- cor.test(b1_anio$b1_medio, b1_anio$ONI_anual); ct_R2 <- cor.test(b1_anio$R2_medio, b1_anio$ONI_anual)
  cat(sprintf("Con el ONI (23 años): b1 medio r = %.2f (p = %.3f) | R² medio r = %.2f (p = %.3f)\n", ct_b1$estimate, ct_b1$p.value, ct_R2$estimate, ct_R2$p.value))
  print(b1_anio[, .(b1_medio = round(mean(b1_medio), 3), R2_medio = round(mean(R2_medio), 3), anios = .N), by = fase])
  lim_b <- c(-.5, 2)
  grafico_base("F4_concurrente_b1_por_anio", function() {
    par(mfrow = c(6, 4), mar = c(3, 3.4, 2, .6), mgp = c(1.9, .5, 0), oma = c(0, 0, 3, 0))
    for (p in por_anio) {
      plot(NA, xlim = RANGO, ylim = lim_b, xaxt = "n", xlab = "Día del año", ylab = expression(b[1](t)),
           main = sprintf("%d · %s (R² medio %.2f)", p$anio, b1_anio[anio == p$anio, fase], mean(p$R2)), cex.main = .95, cex.axis = .75, cex.lab = .8)
      axis(1, at = c(1, 91, 182, 274, 365), labels = c("ene", "abr", "jul", "oct", "dic"), cex.axis = .75)
      polygon(c(T_DIA, rev(T_DIA)), pmin(pmax(c(p$b1 + 2 * p$ee, rev(p$b1 - 2 * p$ee)), lim_b[1]), lim_b[2]), col = "grey85", border = NA)
      abline(h = 0, lty = 3); abline(h = 1, lty = 2, col = "grey40"); lines(T_DIA, p$b1, lwd = 1.5)
    }
    par(mar = c(3.5, 3.4, 2, .6))
    plot(b1_anio$ONI_anual, b1_anio$b1_medio, pch = 19, xlab = "ONI del año", ylab = "b1 medio del año", cex.axis = .75, cex.lab = .8,
         main = sprintf("b1 medio vs ONI: r = %.2f", ct_b1$estimate), cex.main = .95)
    abline(lm(b1_medio ~ ONI_anual, data = b1_anio), lty = 2)
    mtext("Coeficiente del satélite b1(t) del modelo concurrente ajustado en cada año (banda: ±2 EE punto a punto; punteada: 0; discontinua: 1)", outer = TRUE, cex = .95, font = 2)
  }, ancho = 11, alto = 15)

  # Contraste con fda::fRegress (herramienta de clase): mismo modelo A, curvas, K_b y lambda.
  base_b <- create.bspline.basis(RANGO, nbasis = eleccion$A$K_b, norder = ORDEN)
  fit_fda <- fRegress(fd(t(a_raiz(o$coef_est[pares$fila_est[i_ent], ], Phi_e)$coef), create.bspline.basis(RANGO, nbasis = o$K_est, norder = ORDEN)),
                      list(const = rep(1, length(i_ent)), sat = fd(t(a_raiz(o$coef_sat[pares$fila_sat[i_ent], ], Phi_s)$coef), create.bspline.basis(RANGO, nbasis = o$K_sat, norder = ORDEN))),
                      list(const = fdPar(base_b, int2Lfd(2), eleccion$A$lambda), sat = fdPar(base_b, int2Lfd(2), eleccion$A$lambda)))
  b0_f <- drop(eval.fd(T_DIA, fit_fda$betaestlist$const$fd)); b1_f <- drop(eval.fd(T_DIA, fit_fda$betaestlist$sat$fd))
  sse <- function(b0, b1) sum(sweep((Y_r[i_ent, ] - sweep(sweep(X_r[i_ent, ], 2, b1, "*"), 2, b0, "+"))^2, 2, W_T, "*"))
  P_A <- predecir(bA, i_pru, "A"); P_f <- sweep(sweep(X_r[i_pru, ], 2, b1_f, "*"), 2, b0_f, "+"); centro <- 11:355
  contraste <- data.table(cantidad = c("SSE de entrenamiento, solución propia", "SSE de entrenamiento, fRegress", "Diferencia relativa del SSE (%)",
                                       "RMS de la diferencia de predicción en prueba (raíz)", "Máx. |b0 propio - b0 fRegress| (centro del año)", "Máx. |b1 propio - b1 fRegress| (centro del año)",
                                       "RMSE L2 en prueba, mm/día (propio)", "RMSE L2 en prueba, mm/día (fRegress)"),
                          valor = signif(c(sse(bA[, 1], bA[, 2]), sse(b0_f, b1_f), 100 * (sse(b0_f, b1_f) / sse(bA[, 1], bA[, 2]) - 1), sqrt(mean((P_A - P_f)^2)),
                                           max(abs(bA[centro, 1] - b0_f[centro])), max(abs(bA[centro, 2] - b1_f[centro])),
                                           mean(rmse_l2(a_mm(P_A), i_pru)), mean(rmse_l2(a_mm(P_f), i_pru))), 4))
  print(tabla_f(contraste, "F4_concurrente_contraste_fRegress"))

  # Comprobaciones paso a paso (van a la sección 6).
  p99 <- por_anio[[which(anios_c == 1999)]]; cf <- summary(lm(Y_r[p99$i, 105] ~ X_r[p99$i, 105]))$coefficients
  X_bak <- X_r; Y_bak <- Y_r; X_r[c(i_val, i_pru), ] <- 0; Y_r[c(i_val, i_pru), ] <- 0
  th_nf <- resolver(preparar(i_ent, "A", eleccion$A$K_b), eleccion$A$lambda); X_r <- X_bak; Y_r <- Y_bak
  K_t <- 12L; n_t <- length(i_ent); Pb <- base_bspline(T_DIA, K_t); sw <- sqrt(W_T[rep(T_DIA, n_t)])
  Xd <- cbind(Pb[rep(T_DIA, n_t), ] * sw, (Pb[rep(T_DIA, n_t), ] * as.vector(t(X_r[i_ent, ]))) * sw)
  dif_mco <- max(abs(qr.coef(qr(Xd), as.vector(t(Y_r[i_ent, ])) * sw) - resolver(preparar(i_ent, "A", K_t), 0))); rm(Xd)
  th_recta <- solve(crossprod(Pb), crossprod(Pb, 3 + .01 * T_DIA)); b1_rigido <- coef_t(preparar(i_ent, "A", K_t), resolver(preparar(i_ent, "A", K_t), 1e14))[, 2]
  P_mm <- list(A_val = a_mm(predecir(bA, i_val, "A")), A_pru = a_mm(P_A), B_val = a_mm(predecir(ajuste$B$B, i_val, "B")), B_pru = a_mm(predecir(ajuste$B$B, i_pru, "B")))
  pruebas <- c(lm_por_anio = length(por_anio) == 23L && abs(cf[2, 1] - p99$b1[105]) < 1e-10 && abs(cf[2, 2] - p99$ee[105]) < 1e-10,
               sin_fuga = max(abs(th_nf - ajuste$A$theta)) < 1e-10, ecuaciones_normales = dif_mco < 1e-8,
               rugosidad_recta = max(abs(gram_derivada(K_t, 2L) %*% th_recta)) < 1e-6, lambda_enorme = max(abs(diff(b1_rigido, differences = 2))) < 1e-6,
               no_negativa = all(vapply(P_mm, function(P) all(P >= 0), TRUE)), bootstrap = B_BOOT >= 300L)
  cat("Comprobaciones de la regresión concurrente:", paste(names(pruebas), ifelse(pruebas, "OK", "FALLA"), collapse = " | "), "\n")
  list(eleccion = lapply(eleccion, `[`, c("K_b", "lambda")), B_A = bA, EE_A = eA, B_B = ajuste$B$B, bandas = bandas, intervalos_b1 = b1_bajo_1,
       R2 = R2_tab, por_mes = por_mes, por_anio = b1_anio, contraste_fRegress = contraste, i_val = i_val, i_pru = i_pru, P_mm = P_mm,
       recorte_pct = 100 * c(A = mean(P_A < 0), B = mean(predecir(ajuste$B$B, i_pru, "B") < 0)), pruebas = pruebas)
})

titulo("4.5 ANOVA funcional del error del satélite (clase, semana 14)")
# ¿Cambia la diferencia estación - satélite según la fase ENSO, la vertiente o la altitud?
# Pruebas de la clase: F puntual F_n(t), tipo F global (Zhang), GPF = ∫ F_n(t) dt y
# F_max = sup F_n(t), calibradas por permutación de las etiquetas de grupo.
# Unidades independientes: para el ENSO, los 23 AÑOS (curva media de d de cada año; la
# fase es del año, no de la estación); para vertiente y altitud, las ESTACIONES: D_j(t) es la
# diferencia media de los años de la estación j. Las trazas del operador de covarianza usan los
# estimadores insesgados de Zhang (2013), como en la clase.
B_PERM <- 5000L
w_trap <- c(.5, rep(1, 363), .5)                         # regla del trapecio en la malla diaria
trazas_zhang <- function(G, N, w = w_trap) {             # G: covarianza evaluada en la malla; N: sus grados de libertad
  DG <- G * w; t1 <- sum(diag(DG)); t2 <- sum(DG * t(DG))  # tr(D G) y tr((D G)^2)
  t2_ins <- N^2 / ((N - 1) * (N + 2)) * (t2 - t1^2 / N)
  t1c_ins <- N * (N + 1) / ((N - 1) * (N + 2)) * (t1^2 - 2 * t2 / (N + 1))
  list(tr = t1, beta = t2_ins / t1, d = t1c_ins / t2_ins)
}
D_est <- rowsum(D, pares$codigo) / as.vector(table(pares$codigo)[rownames(rowsum(D, pares$codigo))])
vert_est <- pares[, .(vertiente = vertiente[1], banda_altitud = banda_altitud[1]), by = codigo][match(rownames(D_est), codigo)]
anova_funcional <- function(Xg, g, B_perm = B_PERM) {
  g <- factor(g); k <- nlevels(g); n <- nrow(Xg)
  estad <- function(gg) {
    nk <- tabulate(gg, k); med <- rowsum(Xg, gg, reorder = TRUE) / nk; mg <- colMeans(Xg)
    SSR <- colSums(nk * sweep(med, 2, mg)^2); SSE <- colSums((Xg - med[gg, , drop = FALSE])^2)
    list(Ft = (SSR / (k - 1)) / (SSE / (n - k)), SSRn = sum(w_trap * SSR), SSEn = sum(w_trap * SSE), med = med)
  }
  o_ <- estad(as.integer(g)); E <- Xg - o_$med[as.integer(g), , drop = FALSE]; z <- trazas_zhang(crossprod(E) / (n - k), n - k)
  Fn <- (o_$SSRn / (k - 1)) / (o_$SSEn / (n - k))
  perm <- replicate(B_perm, { e <- estad(sample(as.integer(g))); c(GPF = sum(w_trap * e$Ft), Fmax = max(e$Ft), Fn = (e$SSRn / (k - 1)) / (e$SSEn / (n - k))) })
  Ft_perm <- replicate(min(B_perm, 2000L), estad(sample(as.integer(g)))$Ft)
  list(k = k, n = n, niveles = levels(g), Ft = o_$Ft, med = o_$med, Fn = Fn, d = z$d, p_F = pf(Fn, (k - 1) * z$d, (n - k) * z$d, lower.tail = FALSE),
       GPF = sum(w_trap * o_$Ft), Fmax = max(o_$Ft), p_GPF = (1 + sum(perm["GPF", ] >= sum(w_trap * o_$Ft))) / (B_perm + 1),
       p_Fmax = (1 + sum(perm["Fmax", ] >= max(o_$Ft))) / (B_perm + 1), p_F_perm = (1 + sum(perm["Fn", ] >= Fn)) / (B_perm + 1),
       c_alfa = apply(Ft_perm, 1, quantile, .95), c_max = unname(quantile(perm["Fmax", ], .95)), SSEn = o_$SSEn, tr = z$tr)
}
set.seed(20261118)
D_anio <- rowsum(D, pares$anio) / as.vector(table(pares$anio)[rownames(rowsum(D, pares$anio))])
fase_anio <- oni[match(as.integer(rownames(D_anio)), anio), fase]
cauca_est <- vert_est$vertiente != "Pacífico"
av <- list(`Fase ENSO (23 años)` = anova_funcional(D_anio, fase_anio),
           `Vertiente (176 estaciones)` = anova_funcional(D_est, vert_est$vertiente),
           `Altitud en el Cauca (150 estaciones)` = anova_funcional(D_est[cauca_est, ], vert_est$banda_altitud[cauca_est]))
tabla_anova <- rbindlist(lapply(names(av), function(nm) with(av[[nm]], data.table(
  factor = nm, grupos = paste(niveles, collapse = " / "), k = k, n = n, F_n = round(Fn, 2), d = round(d, 2), p_F = signif(p_F, 2),
  p_F_permutacion = signif(p_F_perm, 2), GPF = round(GPF), p_GPF = signif(p_GPF, 2), F_max = round(Fmax, 2), p_F_max = signif(p_Fmax, 2),
  dias_F_puntual = sum(Ft > c_alfa), dias_criterio_maximo = sum(Ft > c_max)))))
print(tabla_f(tabla_anova, "F4_anova_funcional"))
# Comprobación de la clase (ejercicio 1): SSE_n = (n - k) tr(γ̂); y con k = 2 la F puntual es t² de dos muestras.
chequeo_anova <- c(SSE = max(sapply(av, function(a) abs(a$SSEn - (a$n - a$k) * a$tr) / a$SSEn)),
                   t2 = { a <- av[[2]]; gv <- factor(vert_est$vertiente); x1 <- D_est[gv == levels(gv)[1], 100]; x2 <- D_est[gv == levels(gv)[2], 100]
                          abs(a$Ft[100] - t.test(x1, x2, var.equal = TRUE)$statistic^2) / a$Ft[100] })
cat(sprintf("Comprobaciones: SSE_n = (n-k) tr(γ̂) a %.1e relativo; con k = 2, F(t) = t² a %.1e relativo\n", chequeo_anova[1], chequeo_anova[2]))
grafico_base("F4_anova_funcional", function() {
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
  cols <- c("grey30", "steelblue", "firebrick", "darkgreen")
  for (nm in names(av)) with(av[[nm]], {
    matplot(T_DIA, t(med), type = "l", lty = 1, lwd = 2, col = cols[seq_len(k)], xlab = "Día del año", ylab = "d(t), mm/día", main = nm)
    abline(h = 0, lty = 3); legend("bottomleft", niveles, col = cols[seq_len(k)], lwd = 2, bty = "n", cex = .8)
  })
  for (nm in names(av)) with(av[[nm]], {
    plot(T_DIA, Ft, type = "l", lwd = 2, ylim = c(0, max(Ft, c_max) * 1.05), xlab = "Día del año", ylab = "F(t)",
         main = sprintf("F puntual: GPF p = %.3f, Fmax p = %.3f", p_GPF, p_Fmax))
    lines(T_DIA, c_alfa, lty = 2, col = "steelblue", lwd = 2); abline(h = c_max, lty = 2, col = "firebrick", lwd = 2)
  })
}, ancho = 15, alto = 8)

# QUÉ DICEN LOS DATOS
#   Pares: 3 299 estación-año / píxel-año, 176 estaciones en 148 píxeles, todos
#   los años. El EE se calcula entre estaciones (promedio de cada estación
#   primero): es 3,7 veces el que saldría tratando cada curva como independiente,
#   porque los años de una misma estación repiten su sesgo.
#   d(t): en el promedio anual casi no hay diferencia (el satélite da 22 mm/año
#   más), pero el signo cambia con la temporada. La estación supera al satélite
#   del 8-feb al 28-mar (inicio de la primera temporada: el satélite la arranca
#   tarde); el satélite supera a la estación del 15-abr al 9-jun y del 28-oct al
#   18-nov, es decir, exagera los picos de las dos temporadas. Es el mismo
#   retraso que mostraron las derivadas. La banda excluye el 0 en 147 de 365
#   días. Los tramos del 1 al 12 de enero y del 24 al 31 de diciembre son efecto
#   de borde de las bases. En 1985, 1988 y 2000 el satélite pierde 337 mm/año
#   frente a la estación; en el resto de años da 83 mm/año de más.
#   Contraste de la banda: con percentiles de un bootstrap de estaciones la banda
#   excluye el 0 en 152 días (356 con el mismo veredicto que ±2 EE): la
#   aproximación normal no cambia nada. Con el píxel como unidad (148) la banda
#   se ensancha 14 % y excluye el 0 en 137 días; siguen los mismos tramos (el
#   satélite arranca tarde la primera temporada, 7-feb a 30-mar, y exagera su pico,
#   20-abr a 5-jun), pero el de la segunda temporada se reduce a 5-12 de noviembre.
#   Las estaciones que comparten píxel no cambian la lectura de d(t).
#   Serie diaria y ciclo anual (4.1b): en Yurumanguí (costa Pacífica, 21 m) la
#   estación promedia 22,4 mm/día y su píxel 14,2: el satélite ve la forma
#   bimodal pero se queda en dos tercios del nivel. En Puerto Molina (fondo del
#   valle, 939 m) pasa lo contrario, 3,9 frente a 5,4 mm/día: el satélite
#   exagera los dos picos. En La Gitana (2 783 m) el píxel da 4,1 frente a 5,1 y
#   aplana el pico de noviembre. La forma del ciclo anual coincide en los tres
#   (r = 0,83 a 0,86): el error del satélite es sobre todo de nivel y de
#   amplitud, y cambia de signo según el sitio, como resume d(t).
#   Correlación cruzada: la diagonal r(t,t) vale 0,62 en promedio (0,78 a
#   mediados de julio) y la correlación máxima está en el mismo día (desfase
#   mediano de -1 día). Casi toda viene de que el satélite reconoce qué sitios
#   son húmedos y cuáles secos: con anomalías (restando la climatología de cada
#   estación) baja a 0,24, con un máximo de 0,47 en febrero y 0,07 en noviembre.
#   El satélite sigue mal la variación de un año a otro, sobre todo en la segunda
#   temporada.
#   Covarianzas: la forma es casi igual (correlación abril-octubre 0,74 y 0,72;
#   la primera componente explica 75 % y 77 %), pero la varianza del satélite es
#   una cuarta parte de la de las estaciones (mediana 0,24 por día): el píxel de
#   ~5,5 km y la péntada promedian los aguaceros locales.
#   Mismo K: al construir el satélite también con K = 8, la diferencia media
#   (-20 frente a -22 mm/año), la diagonal (0,64 frente a 0,62) y el cociente de
#   varianzas (0,25 frente a 0,24) casi no cambian: el diagnóstico global no sale
#   de usar K distinto. El veredicto día a día sí es sensible: coincide en 281 de
#   365 días (77 %) y la banda excluye el 0 en 159 días en vez de 147. Por eso las
#   fechas exactas de los días significativos se reportan con esa salvedad.
#   Subregión: en el Cauca el satélite sobreestima 8 % (+109 mm/año) y la banda
#   excluye el 0 en 263 días. En la vertiente Pacífica el sesgo depende del
#   sitio. El satélite subestima a la mitad las estaciones más lluviosas (El
#   Tabano, San Francisco y Bahía Málaga, cociente 0,43 a 0,50) y duplica la
#   lluvia de los sitios secos a sotavento (Loboguerrero 2,6; La Despensa 2,3;
#   Puerto Nuevo 2,3; Moravito 1,9). El satélite suaviza hacia la media regional:
#   sube lo seco y baja lo muy húmedo. Entre estaciones, la lluvia anual del píxel
#   explica la de la estación con r = 0,79, pero dentro del Cauca solo con 0,42.
#   Altitud: en el Cauca el cociente satélite/estación baja 15 % por cada 1 000 m
#   (p = 0,001): el satélite exagera en el fondo del valle y se queda corto en
#   las laderas altas, donde la lluvia orográfica no se ve a 5 km de resolución.
#   ENSO: en años Niña las estaciones reciben 1 926 mm y el satélite 1 731
#   (subestima 10 %); en años Niño 1 570 frente a 1 833 (sobreestima 17 %). El
#   satélite amortigua la señal ENSO. Esto justifica dar al modelo la altitud y
#   el ONI junto con la curva satelital, como pide el enunciado.
#   Regresión concurrente (4.4, ejercicio 3 de clase): la validación elige el
#   modelo A con K_b = 12 y lambda = 0 (el B con K_b = 8 y lambda = 1e8). En
#   prueba el satélite explica mejor a la estación de mayo a octubre (R²(t) de
#   0,44 a 0,59; máximo 0,62 el 24-jun) y peor en enero-marzo (0,17 a 0,19) y en
#   diciembre (0,27); el satélite sin corregir da R² medio 0,29 y el modelo 0,37.
#   Con la banda por bootstrap de estaciones, b1(t) queda bajo 1 del 10-feb al
#   29-abr y del 1-nov al 1-dic: en las dos temporadas lluviosas conviene encoger
#   la señal del satélite (b1 de 0,68 en marzo y 0,73 en noviembre). De junio a
#   septiembre b1 vale cerca de 1 y la banda lo contiene; ningún día b1 > 1.
#   Bandas: la sándwich de clase es 4 veces más angosta que la del bootstrap de
#   estaciones y marcaría b1 distinto de 1 en 262 días en vez de 119. El
#   bootstrap de píxeles (148) da casi lo mismo (112 días): las estaciones que
#   comparten píxel no cambian la conclusión. Los residuos están muy
#   correlacionados entre días cercanos (0,7 a 0,9) y poco a varios meses (0,2 a
#   0,5); por eso las bandas se leen punto a punto.
#   Linealidad: el residuo medio por quintil del satélite se aparta de 0 hasta
#   0,16 raíz(mm/día), con una DE de los residuos de 0,63. Donde el satélite marca
#   poco, la estación recibe algo más de lo que da la recta, y en el cuarto
#   quintil algo menos. Es una no linealidad moderada que el modelo lineal no
#   captura y que justifica probar una cabeza no lineal (KAN).
#   Año por año, el b1 medio baja cuando sube el ONI (r = -0,44, p = 0,036): en
#   años Niño el satélite exagera más las variaciones entre sitios (b1 medio 0,91
#   frente a 1,02 en Niña). fRegress llega a la misma solución (SSE 0,4 % mayor;
#   RMSE de prueba 2,145 frente a 2,140 mm/día).

comparacion <- list(pares = pares, d = b, d_marcados = b_def, intervalos = int_d, R_YX = R_YX, Cov_YY = Cov_YY, Cov_XX = Cov_XX,
             grupos = grupos, sesgo_estacion = por_est, concurrente = concurrente[setdiff(names(concurrente), "P_mm")], anova = tabla_anova)

# =============================================================================
# 5. REGRESIÓN FUNCIONAL FUNCIÓN A FUNCIÓN (SATÉLITE -> ESTACIÓN)  — rúbrica 4.2 y 4.4
# =============================================================================
# Modelo: y_i(t) = mu(t) + alfa(t) + integral beta(s,t) x_i(s) ds
#                  + gamma(t) altitud_i + delta(t) ONI_i + e_i(t)
# sobre la raíz de la lluvia (curvas en mm recortadas en 0 y luego con raíz, D-23).
# beta(s,t) = sum_m beta_m(s) psi_m(t): psi_m son las funciones propias del
# FPCA de las estaciones (solo entrenamiento) y beta_m(s) vive en la base
# B-spline del satélite con penalización de rugosidad lambda * integral beta_m''^2.
# Es el modelo de rango bajo al que equivale la cabeza lineal de la red
# funcional. M y lambda se eligen con validación (1996-2000) y se evalúa en
# prueba (2001-2003) frente al satélite sin corregir. Entran todos los años.
# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ┃ RÚBRICA 4.2 / 4.4 — Regresión funcional (= cabeza lineal de la red) y su evaluación
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
titulo("5. Regresión funcional")
set.seed(20261113)

K_e <- o$K_est; K_s <- o$K_sat
Phi_e <- base_bspline(T_DIA, K_e); Phi_s <- base_bspline(T_DIA, K_s)
W_e <- o$W_est; W_s <- o$W_sat
R_s <- gram_derivada(K_s, 2L)
ref_d2 <- splines::splineDesign(nodos(K_s), T_DIA, ord = ORDEN, derivs = rep(2, 365), outer.ok = TRUE)
cat("Segunda derivada de la base propia vs splines::splineDesign:", format(max(abs(derivada_base(T_DIA, K_s, 2L) - ref_d2)), digits = 2), "\n")

titulo("5.1 Escala raíz")
re <- a_raiz(o$coef_est[pares$fila_est, ], Phi_e); rs <- a_raiz(o$coef_sat[pares$fila_sat, ], Phi_s)
print(data.table(fuente = c("estaciones", "satélite"), pct_dias_recortados_en_0 = round(c(re$pct_recortado, rs$pct_recortado), 2),
                 rmse_reproyeccion_raiz = signif(c(re$rmse_reproyeccion, rs$rmse_reproyeccion), 3)))

# ¿Suavizar en mm y luego sacar la raíz, o sacar la raíz al dato diario y luego suavizar?
# Se comprueba cuál conserva el acumulado anual (métrica de la rúbrica 4.4), con las
# curvas completas de entrenamiento y el mismo K.
i_cmp <- which(id_est$n_obs == 365L & id_est$particion == "entrenamiento")
mapa_e <- solve(crossprod(Phi_e), t(Phi_e)); obs_mm <- rowMeans(Y_dia[i_cmp, ])
camino_A <- pmax(Y_dia[i_cmp, ] %*% t(mapa_e) %*% t(Phi_e), 0)
camino_B <- (sqrt(Y_dia[i_cmp, ]) %*% t(mapa_e) %*% t(Phi_e))^2
orden_raiz <- data.table(camino = c("A. Suavizar en mm, recortar en 0 y luego raíz (el que se usa)", "B. Raíz al dato diario, suavizar y elevar al cuadrado"),
                         curvas = length(i_cmp), cociente_acumulado_anual_mediano = round(c(median(rowMeans(camino_A) / obs_mm), median(rowMeans(camino_B) / obs_mm)), 3),
                         cociente_p05 = round(c(quantile(rowMeans(camino_A) / obs_mm, .05), quantile(rowMeans(camino_B) / obs_mm, .05)), 3),
                         cociente_p95 = round(c(quantile(rowMeans(camino_A) / obs_mm, .95), quantile(rowMeans(camino_B) / obs_mm, .95)), 3))
print(tabla_f(orden_raiz, "F5_orden_raiz_y_suavizado"))
stopifnot(abs(orden_raiz$cociente_acumulado_anual_mediano[1] - 1) < .02)
Y_mm <- pmax(o$coef_est[pares$fila_est, ] %*% t(Phi_e), 0)
X_mm <- pmax(o$coef_sat[pares$fila_sat, ] %*% t(Phi_s), 0)

pares[, usar := !is.na(altitud)]   # todos los años entran (D-47)
cat("Pares fuera del modelo por no tener altitud:", pares[is.na(altitud), .N], "\n")
print(pares[usar == TRUE, .(pares = .N, estaciones = uniqueN(codigo), anios = paste(range(anio), collapse = "-")), by = particion])
ent <- pares$usar & pares$particion == "entrenamiento"
val <- pares$usar & pares$particion == "validación"
pru <- pares$usar & pares$particion == "prueba"

titulo("5.2 FPCA de las estaciones (solo entrenamiento)")
ev <- eigen(W_e, symmetric = TRUE); W12 <- ev$vectors %*% diag(sqrt(ev$values)) %*% t(ev$vectors); Wm12 <- ev$vectors %*% diag(1 / sqrt(ev$values)) %*% t(ev$vectors)
# El FPCA es una función de los coeficientes y de las filas que lo ajustan: así la prueba
# de no fuga (sección 6) lo puede volver a correr con los datos de fuera alterados.
fpca_ent <- function(coef, filas) {
  mu <- colMeans(coef[filas, ]); Cc <- sweep(coef[filas, ], 2, mu)
  eg <- eigen(W12 %*% (crossprod(Cc) / (sum(filas) - 1)) %*% W12, symmetric = TRUE)
  list(mu = mu, values = eg$values, psi = Wm12 %*% eg$vectors)
}
eg <- fpca_ent(re$coef, ent); mu_c <- eg$mu; psi_c <- eg$psi; psi_todas <- psi_c
var_fpca <- data.table(m = seq_len(K_e), autovalor = eg$values, pct = 100 * eg$values / sum(eg$values))[, acumulado := cumsum(pct)]
print(tabla_f(var_fpca, "F5_fpca_varianza"))
stopifnot(max(abs(t(psi_c) %*% W_e %*% psi_c - diag(K_e))) < 1e-10)
fd_r <- fda::fd(t(re$coef[ent, ]), fda::create.bspline.basis(RANGO, nbasis = K_e, norder = ORDEN))
pca_fda <- fda::pca.fd(fd_r, nharm = 4, centerfns = TRUE)
psi_fda <- fda::eval.fd(T_DIA, pca_fda$harmonics)
cat(sprintf("Contra fda::pca.fd (4 primeras): autovalores x N/(N-1) a %.1e relativo, funciones propias a %.1e (salvo signo); fda integra W numéricamente\n",
            max(abs(pca_fda$values[1:4] * sum(ent) / (sum(ent) - 1) - eg$values[1:4]) / eg$values[1:4]),
            max(abs(abs(psi_fda) - abs(Phi_e %*% psi_c[, 1:4])))))
S_todos <- sweep(re$coef, 2, mu_c) %*% W_e %*% psi_c

titulo("5.2b P-19: registros de 1930-1979 y estabilidad del FPCA")
h30 <- rbindlist(lapply(readxl::excel_sheets("datos_curso/P_1930_1979.xlsx"), function(s) as.data.table(readxl::read_excel("datos_curso/P_1930_1979.xlsx", sheet = s))))
h30 <- h30[as.character(STATION_ID) %in% unique(id_est$codigo) & !(MONTH == 2 & DAY == 29) & !is.na(VALUE) & VALUE >= 0 & VALUE <= 1825,
           .(codigo = as.character(STATION_ID), anio = as.integer(YEAR), dia = yday(as.IDate(sprintf("%d-%02d-%02d", YEAR, MONTH, DAY))) - (YEAR %% 4 == 0 & MONTH > 2), P = VALUE)]
h30 <- unique(h30, by = c("codigo", "anio", "dia"))
anios30 <- h30[, .(n_obs = .N, hueco_max = { r <- rle(!(1:365 %in% dia)); if (any(r$values)) max(r$lengths[r$values]) else 0L }), by = .(codigo, anio)]
anios30[, valido := n_obs >= COB_MIN_R5 * 365 & hueco_max <= HUECO_MAX_R5]
h30 <- h30[anios30[valido == TRUE, .(codigo, anio)], on = .(codigo, anio)]
h30[, fila := .GRP, by = .(codigo, anio)]
Y30 <- matrix(NA_real_, max(h30$fila), 365); Y30[cbind(h30$fila, h30$dia)] <- h30$P
C30 <- a_raiz(construir(Y30, T_DIA, K_e)$coef, Phi_e)$coef
eg_30 <- fpca_ent(rbind(re$coef[ent, ], C30), rep(TRUE, sum(ent) + nrow(C30)))
explicada <- function(C, mu, psi, M = 6L) {
  D <- sweep(C, 2, mu); R <- D - (D %*% W_e %*% psi[, 1:M]) %*% t(psi[, 1:M])
  100 * (1 - sum(rowSums((R %*% W_e) * R)) / sum(rowSums((D %*% W_e) * D)))
}
M_ref <- 6L
p19 <- data.table(
  curvas = c("Entrenamiento 1981-1995", "Validación 1996-2000", "Prueba 2001-2003", "Registros 1930-1979 válidos"),
  n = c(sum(ent), sum(val), sum(pru), nrow(C30)),
  pct_explicado_FPCA_entrenamiento = round(c(explicada(re$coef[ent, ], mu_c, psi_c, M_ref), explicada(re$coef[val, ], mu_c, psi_c, M_ref),
                                             explicada(re$coef[pru, ], mu_c, psi_c, M_ref), explicada(C30, mu_c, psi_c, M_ref)), 1))
print(tabla_f(p19, "F5_P19_registros_1930_1979"))
print(tabla_f(anios30[, .(estaciones = uniqueN(codigo), estacion_anios = .N, validos = sum(valido), estaciones_con_validos = uniqueN(codigo[valido]),
                          validos_1970_1979 = sum(valido & anio >= 1970)), by = NULL], "F5_P19_conteo"))
cos_p19 <- svd(t(psi_c[, 1:M_ref]) %*% W_e %*% eg_30$psi[, 1:M_ref])$d
d_mu <- eg_30$mu - mu_c
cat(sprintf("FPCA con 1930-1979 agregados: %.1f %% con M = %d (sin ellos %.1f %%) | cosenos principales entre subespacios: %s | cambio de la media: %.1f %% de su norma\n",
            100 * sum(eg_30$values[1:M_ref]) / sum(eg_30$values), M_ref, 100 * sum(eg$values[1:M_ref]) / sum(eg$values),
            paste(round(cos_p19, 3), collapse = ", "), 100 * sqrt(drop(t(d_mu) %*% W_e %*% d_mu) / drop(t(mu_c) %*% W_e %*% mu_c))))
rm(h30, Y30)

# QUÉ DICEN LOS DATOS (P-19)
#   De las 176 estaciones, 145 tienen registros en 1930-1979 y 140 tienen años que
#   cumplen R5 (>= 90 % de días y ningún hueco > 15): 1 250 curvas, 963 de ellas
#   de los años setenta. Sumarían un 55 % a las 2 257 curvas de entrenamiento.
#   El FPCA de entrenamiento ya las representa: con M = 6 explica el 89,8 % de su
#   variabilidad, cerca del 91,5 % de validación y del 95,2 % de prueba.
#   Si se agregan al ajuste, las primeras cuatro funciones propias no cambian
#   (coseno 1,000), la quinta casi no (0,994) y solo gira algo la sexta (0,912);
#   la media se mueve un 1,7 %. La varianza explicada con M = 6 baja de 94,2 % a
#   93,2 %: esas curvas traen variabilidad de otra época, no estructura nueva.
#   Decisión (P-19): no se usan. La fe de erratas fija 1981-2003, el decodificador
#   sería prácticamente el mismo y habría que depurar 50 años más con R0 a R7.

titulo("5.3 Ajuste y elección de M y lambda con validación")
alt_m <- mean(pares$altitud[ent]); alt_s <- sd(pares$altitud[ent])
Zx <- rs$coef %*% W_s
# Los datos entran como argumentos con valor por defecto (los reales) para que la prueba
# de no fuga pueda reajustar con datos alterados sin tocar los objetos del análisis.
diseno <- function(i, con_sat = TRUE, zx = Zx, alt = pares$altitud, oni = pares$ONI_anual, am = alt_m, ds = alt_s)
  cbind(1, if (con_sat) zx[i, , drop = FALSE], (alt[i] - am) / ds, oni[i])
penal <- function(con_sat = TRUE) { P <- matrix(0, 3 + con_sat * K_s, 3 + con_sat * K_s); if (con_sat) P[1 + 1:K_s, 1 + 1:K_s] <- R_s; P }
ajustar_B <- function(lambda, M, con_sat = TRUE, S = S_todos, ...) {
  Z <- diseno(which(ent), con_sat, ...); solve(crossprod(Z) + lambda * penal(con_sat), crossprod(Z, S[ent, 1:M, drop = FALSE]))
}
predecir_raiz <- function(B, i, con_sat = TRUE, psi = psi_c, ...) {
  M <- ncol(B); sweep((diseno(i, con_sat, ...) %*% B) %*% t(Phi_e %*% psi[, 1:M, drop = FALSE]), 2, drop(Phi_e %*% mu_c), "+")
}
a_mm <- function(R) pmax(R, 0)^2
rmse_l2 <- function(P, i, Y = Y_mm) sqrt(rowMeans((P - Y[i, ])^2))

lambdas <- 10^seq(-2, 10, by = .5)
tabla_seleccion <- function(S = S_todos, Y = Y_mm, ...) rbindlist(lapply(seq_len(K_e), function(M) rbindlist(lapply(lambdas, function(l) {
  B <- ajustar_B(l, M, S = S, ...)
  data.table(M = M, lambda = l, rmse_entrenamiento = mean(rmse_l2(a_mm(predecir_raiz(B, which(ent), psi = psi_todas, ...)), which(ent), Y)),
             rmse_validacion = mean(rmse_l2(a_mm(predecir_raiz(B, which(val), psi = psi_todas, ...)), which(val), Y)))
}))))
# Regla: el menor M cuyo mejor RMSE de validación queda a menos de 0,5 % del mínimo.
elegir_M_lambda <- function(sel) {
  mejor <- sel[, .SD[which.min(rmse_validacion)], by = M]
  M_sel <- mejor[rmse_validacion <= min(rmse_validacion) * 1.005, min(M)]
  list(mejor = mejor, M = M_sel, lambda = mejor[M == M_sel, lambda])
}
sel_lambda <- tabla_seleccion()
tabla_f(sel_lambda, "F5_seleccion_M_lambda")
eleccion <- elegir_M_lambda(sel_lambda); mejor_por_M <- eleccion$mejor
print(mejor_por_M)
LAMBDA <- eleccion$lambda; M <- eleccion$M
psi_c <- psi_c[, 1:M, drop = FALSE]; psi_t <- Phi_e %*% psi_c; S <- S_todos[, 1:M, drop = FALSE]
cat(sprintf("Elegidos: M = %d (%.1f %% de la varianza; el menor M a menos de 0,5 %% del mejor RMSE de validación) y lambda = %s\n",
            M, var_fpca$acumulado[M], format(LAMBDA, scientific = TRUE)))
grafico_f(ggplot(sel_lambda, aes(lambda, rmse_validacion, colour = factor(M))) + geom_line() +
            geom_point(data = mejor_por_M, size = 2) + geom_vline(xintercept = LAMBDA, linetype = 2, colour = "red") + scale_x_log10() +
            labs(x = expression(lambda), y = "RMSE L2 medio en validación (mm/día)", colour = "M",
                 title = "Elección de M y lambda con los años de validación (1996-2000)"),
          "F5_seleccion_M_lambda", ancho = 9, alto = 5)

B <- ajustar_B(LAMBDA, M); B_sin_sat <- ajustar_B(0, M, con_sat = FALSE)
gl_df <- function(Z, P, l) sum(diag(solve(crossprod(Z) + l * P, crossprod(Z))))
cat(sprintf("Grados de libertad efectivos por componente: %.1f de %d parámetros\n", gl_df(diseno(which(ent)), penal(), LAMBDA), ncol(diseno(1))))

titulo("5.4 Superficie beta(s,t) y verificación de la equivalencia")
beta_m <- B[1 + 1:K_s, , drop = FALSE]
BETA <- Phi_s %*% beta_m %*% t(psi_t)
alfa_t <- drop(psi_t %*% B[1, ]); gamma_t <- drop(psi_t %*% B[K_s + 2, ]); delta_t <- drop(psi_t %*% B[K_s + 3, ])
tau_s <- unique(nodos(K_s)); a <- head(tau_s, -1); b <- tail(tau_s, -1)
sq <- as.vector(outer(GL_x, seq_along(a), function(x, k) (b[k] - a[k]) / 2 * x + (a[k] + b[k]) / 2))
wq <- as.vector(outer(GL_w, seq_along(a), function(w, k) (b[k] - a[k]) / 2 * w))
BETA_q <- base_bspline(sq, K_s) %*% beta_m %*% t(psi_t)
i_pru <- which(pru)
x_q <- rs$coef[i_pru, ] %*% t(base_bspline(sq, K_s))
pred_superficie <- sweep(x_q %*% (wq * BETA_q), 2, drop(Phi_e %*% mu_c) + alfa_t, "+") +
  outer((pares$altitud[i_pru] - alt_m) / alt_s, gamma_t) + outer(pares$ONI_anual[i_pru], delta_t)
pred_puntajes <- predecir_raiz(B, i_pru)
dif_equiv <- max(abs(pred_superficie - pred_puntajes))
cat(sprintf("Predicción con la integral de beta(s,t) (cuadratura) vs predicción por puntajes FPCA: diferencia máxima %.2e (criterio <= 1e-4)\n", dif_equiv))
stopifnot(dif_equiv <= 1e-4)

grafico_base("F5_fpca", function() {
  par(mfrow = c(1, 2))
  plot(T_DIA, drop(Phi_e %*% mu_c), type = "l", lwd = 2, xlab = "Día del año", ylab = "raíz(mm/día)", main = "Media de entrenamiento (escala raíz)")
  matplot(T_DIA, psi_t, type = "l", lty = 1, lwd = 2, col = 1:M, xlab = "Día del año", ylab = expression(psi[m](t)),
          main = paste("Funciones propias:", paste(sprintf("%.0f %%", var_fpca$pct[1:M]), collapse = ", ")))
  abline(h = 0, lty = 3)
}, ancho = 12, alto = 4.5)
grafico_base("F5_beta_superficie", function() {
  par(mfrow = c(1, 2))
  sem <- seq(1, 365, by = 7)
  persp(sem, sem, BETA[sem, sem], theta = 35, phi = 30, col = "lightblue", border = "grey50", xlab = "s (satélite)", ylab = "t (estación)",
        zlab = "beta(s,t)", main = expression(hat(beta)(s, t)), ticktype = "simple")
  lim <- max(abs(BETA)) * c(-1, 1)
  image(T_DIA, T_DIA, BETA, zlim = lim, col = hcl.colors(64, "Blue-Red 3"), xlab = "s: día del satélite", ylab = "t: día de la estación", main = "Curvas de nivel")
  contour(T_DIA, T_DIA, BETA, add = TRUE, nlevels = 10, drawlabels = FALSE, col = "grey30"); abline(0, 1, lty = 2)
}, ancho = 13, alto = 5.5)
# Bandas del 95 % de los efectos funcionales. La clase (semana 15) y la tesis del profesor
# (figuras 7.1, 9.1 y 9.2) presentan cada coeficiente con su banda y leen dónde excluye el 0.
# Aquí la banda sale de un bootstrap de ESTACIONES: se remuestrean estaciones de entrenamiento
# (con todos sus años) y se reajusta B con el mismo FPCA, M y lambda. Así la banda respeta que
# los años de una estación no son independientes (la banda sándwich de clase los trata como tales).
B_BOOT_EF <- 500L
filas_est <- split(which(ent), pares$codigo[ent])
efectos_de <- function(Bm) cbind(int_beta = colSums(Phi_s %*% Bm[1 + 1:K_s, , drop = FALSE] %*% t(psi_t)),
                                 gamma = drop(psi_t %*% Bm[K_s + 2, ]), delta = drop(psi_t %*% Bm[K_s + 3, ]))
set.seed(20261119)
boot_ef <- replicate(B_BOOT_EF, {
  f <- unlist(filas_est[sample(length(filas_est), replace = TRUE)], use.names = FALSE)
  Z <- diseno(f); efectos_de(solve(crossprod(Z) + LAMBDA * penal(), crossprod(Z, S[f, , drop = FALSE])))
})
ef_obs <- efectos_de(B); ef_lo <- apply(boot_ef, c(1, 2), quantile, .025); ef_hi <- apply(boot_ef, c(1, 2), quantile, .975)
lejos_borde <- 16:350
bandas_ef <- data.table(efecto = c("Integral de beta sobre s", "Altitud gamma(t)", "ONI delta(t)"),
                        dias_banda_excluye_0 = colSums(ef_lo > 0 | ef_hi < 0),
                        dias_banda_excluye_1 = c(sum(ef_lo[, 1] > 1 | ef_hi[, 1] < 1), NA, NA),
                        dias_banda_excluye_0_lejos_del_borde = colSums((ef_lo > 0 | ef_hi < 0)[lejos_borde, ]),
                        minimo = round(apply(ef_obs, 2, min), 3), maximo = round(apply(ef_obs, 2, max), 3),
                        ancho_medio_banda = round(colMeans(ef_hi - ef_lo), 3))
print(tabla_f(bandas_ef, "F5_efectos_bandas"))
tabla_f(data.table(dia = T_DIA, ef_obs, inf = ef_lo, sup = ef_hi), "F5_efectos_bandas_por_dia")
grafico_base("F5_efectos", function() {
  par(mfrow = c(1, 3))
  titulos <- list(expression(integral(beta(s, t) * ds)), "Efecto de la altitud gamma(t)", "Efecto del ONI delta(t)")
  ejes <- c("", "raíz(mm/día) por 1 DE", "raíz(mm/día) por unidad de ONI"); cols <- c("black", "brown", 4)
  for (j in 1:3) {
    plot(T_DIA, ef_obs[, j], type = "n", ylim = range(c(ef_lo[, j], ef_hi[, j], 0)), xlab = "Día del año", ylab = ejes[j], main = titulos[[j]])
    polygon(c(T_DIA, rev(T_DIA)), c(ef_lo[, j], rev(ef_hi[, j])), col = adjustcolor(cols[j], .2), border = NA)
    lines(T_DIA, ef_obs[, j], lwd = 2, col = cols[j]); abline(h = c(0, if (j == 1) 1), lty = 3)
  }
  mtext(sprintf("Bandas del 95 %% por bootstrap de estaciones (%d remuestreos)", B_BOOT_EF), side = 1, line = -1.2, outer = TRUE, cex = .75)
}, ancho = 14, alto = 4.8)
cat(sprintf("Integral de beta sobre s: máxima el %s, mínima el %s | efecto ONI más negativo el %s (%.3f) | efecto altitud medio %.3f\n",
            fecha_dia(which.max(colSums(BETA))), fecha_dia(which.min(colSums(BETA))), fecha_dia(which.min(delta_t)), min(delta_t), mean(gamma_t)))

titulo("5.5 Evaluación en prueba (2001-2003) con las métricas de la rúbrica")
# Las cuatro métricas de la rúbrica 4.4, por curva:
#   - error L² relativo  = ||ŷ - y|| / ||y||   (y además el error integrado L² en mm/día)
#   - KGE
#   - sesgo del acumulado anual en %
#   - error en la fecha de los dos picos (días)
# IC del 95 % por bootstrap de ESTACIONES (2 000 remuestreos): las curvas de una
# misma estación no son independientes. Prueba de Wilcoxon pareada por estación.
kge <- function(P, Yo) {
  r <- sapply(seq_len(nrow(P)), function(k) cor(P[k, ], Yo[k, ]))
  1 - sqrt((r - 1)^2 + (apply(P, 1, sd) / apply(Yo, 1, sd) - 1)^2 + (rowMeans(P) / rowMeans(Yo) - 1)^2)
}
pico <- function(M, ventana) apply(M[, ventana, drop = FALSE], 1, which.max) + min(ventana) - 1
metricas_de <- function(P, i, modelo) {
  Yo <- Y_mm[i, , drop = FALSE]
  data.table(modelo = modelo, codigo = pares$codigo[i], anio = pares$anio[i], fase = pares$fase[i], particion = pares$particion[i],
             l2_rel = sqrt(rowSums((P - Yo)^2) / rowSums(Yo^2)), rmse = sqrt(rowMeans((P - Yo)^2)), kge = kge(P, Yo),
             sesgo_pct = 100 * (rowMeans(P) - rowMeans(Yo)) / rowMeans(Yo), sesgo_anual = (rowMeans(P) - rowMeans(Yo)) * 365,
             err_pico_1 = abs(pico(P, 60:200) - pico(Yo, 60:200)), err_pico_2 = abs(pico(P, 230:350) - pico(Yo, 230:350)))
}
N_BOOT <- 2000L
boot_ic <- function(v, g) { m <- tapply(v, g, mean); c(mean(m), quantile(replicate(N_BOOT, mean(sample(m, replace = TRUE))), c(.025, .975))) }
boot_media <- function(v, g, dec = 2) { b <- boot_ic(v, g); sprintf(paste0("%.", dec, "f [%.", dec, "f, %.", dec, "f]"), b[1], b[2], b[3]) }
resumir <- function(mc) mc[, .(curvas = .N, estaciones = uniqueN(codigo), L2_relativo = boot_media(l2_rel, codigo), RMSE_L2_mm_dia = boot_media(rmse, codigo),
                               KGE = boot_media(kge, codigo), sesgo_anual_pct = boot_media(sesgo_pct, codigo, 1),
                               error_pico_1_dias = boot_media(err_pico_1, codigo, 1), error_pico_2_dias = boot_media(err_pico_2, codigo, 1))]

modelos <- list(`Satélite sin corregir` = X_mm[i_pru, ], `Regresión funcional` = a_mm(pred_puntajes),
                `Solo altitud y ONI` = a_mm(predecir_raiz(B_sin_sat, i_pru, FALSE)))
recorte_pred <- 100 * mean(pred_puntajes < 0)
metricas_curva <- rbindlist(lapply(names(modelos), function(nm) metricas_de(modelos[[nm]], i_pru, nm)))
metricas <- metricas_curva[, resumir(.SD), by = modelo]
print(tabla_f(metricas, "F5_metricas_prueba"))
tabla_f(metricas_curva, "F5_metricas_prueba_por_curva")
cat(sprintf("Reconstrucción en mm no negativa: %.2f %% de los valores predichos en la escala raíz se recortan en 0\n", recorte_pred))

# ¿La regresión supera al satélite sin corregir? Diferencia pareada por estación.
ancho <- dcast(metricas_curva[modelo != "Solo altitud y ONI"], codigo + anio ~ modelo, value.var = c("l2_rel", "rmse", "kge"))
por_estacion <- metricas_curva[, lapply(.SD, mean), by = .(modelo, codigo), .SDcols = c("l2_rel", "rmse", "kge", "sesgo_pct", "err_pico_1", "err_pico_2")]
pe <- dcast(por_estacion, codigo ~ modelo, value.var = c("l2_rel", "rmse", "kge", "sesgo_pct", "err_pico_1", "err_pico_2"))
comparar <- function(v, mejor_si_baja = TRUE, absoluto = FALSE) {
  a <- pe[[paste0(v, "_Regresión funcional")]]; s <- pe[[paste0(v, "_Satélite sin corregir")]]
  if (absoluto) { a <- abs(a); s <- abs(s) }
  d <- a - s; q <- quantile(replicate(N_BOOT, mean(sample(d, replace = TRUE))), c(.025, .975))
  data.table(metrica = v, satelite = mean(s), regresion = mean(a), diferencia = mean(d), ic_inf = q[1], ic_sup = q[2],
             p_wilcoxon = wilcox.test(a, s, paired = TRUE)$p.value, ic_excluye_0 = q[1] > 0 | q[2] < 0,
             mejora = if (mejor_si_baja) mean(d) < 0 else mean(d) > 0)
}
pareada <- rbind(comparar("l2_rel"), comparar("rmse"), comparar("kge", FALSE), comparar("sesgo_pct", absoluto = TRUE), comparar("err_pico_1"), comparar("err_pico_2"))
pareada[metrica == "sesgo_pct", metrica := "|sesgo_pct|"]
pareada[, c("satelite", "regresion", "diferencia", "ic_inf", "ic_sup") := lapply(.SD, round, 3), .SDcols = c("satelite", "regresion", "diferencia", "ic_inf", "ic_sup")]
pareada[, p_wilcoxon := signif(p_wilcoxon, 2)]
cat("Regresión frente al satélite sin corregir, pareado por estación (", nrow(pe), "estaciones; IC 95 % bootstrap de la diferencia y Wilcoxon):\n")
print(tabla_f(pareada, "F5_regresion_vs_satelite_pareado"))
cat(sprintf("La regresión baja el error L² relativo frente al satélite sin corregir en %.0f %% de las curvas de prueba\n",
            100 * mean(ancho$`l2_rel_Regresión funcional` < ancho$`l2_rel_Satélite sin corregir`)))

# Por fase ENSO. La prueba (2001-2003) no tiene ningún año Niña; por eso se agrega
# la validación (1996-2000), que sí tiene Niña y Niño. Esa parte es ORIENTATIVA:
# la validación se usó para elegir M y lambda (dos números), no para ajustar B.
i_val <- which(val)
mc_val <- rbind(metricas_de(X_mm[i_val, ], i_val, "Satélite sin corregir"), metricas_de(a_mm(predecir_raiz(B, i_val)), i_val, "Regresión funcional"))
por_fase <- rbind(metricas_curva[modelo != "Solo altitud y ONI"], mc_val)[, resumir(.SD), by = .(particion, fase, modelo)]
por_fase[, anios := sapply(seq_len(.N), function(k) paste(sort(unique(pares[particion == por_fase$particion[k] & fase == por_fase$fase[k] & usar, anio])), collapse = " "))]
setcolorder(por_fase, c("particion", "fase", "anios")); setorder(por_fase, -particion, fase, modelo)
print(tabla_f(por_fase, "F5_metricas_por_fase_ENSO"))
por_anio_pru <- metricas_curva[, .(curvas = .N, L2_relativo = round(mean(l2_rel), 3), KGE = round(mean(kge), 2), sesgo_anual_pct = round(mean(sesgo_pct), 1)), by = .(anio, fase, modelo)]
print(tabla_f(por_anio_pru[order(anio, modelo)], "F5_metricas_prueba_por_anio"))

# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ┃ RÚBRICA 4.2 / regla de clase — Contraste con los métodos de regresión vistos en clase
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
titulo("5.5b Contraste con los métodos de regresión vistos en clase")
# La regla del curso: solo métodos vistos en clase. Nuestra regresión (puntajes FPCA de la
# estación sobre ∫φ_j(s)X(s)ds con penalización de rugosidad) es el modelo integral
# funcional-funcional, pero en clase se ajustó de otras dos formas:
#   (a) Taller de repaso, ejercicio 10, paso 3 (Wang et al., 2016): FPCA también de X y
#       regresión de los puntajes de Y sobre los l primeros puntajes de X, B = (ΞᵀΞ)⁻¹ΞᵀH.
#   (b) Taller, ejercicio 9: fda.usc::fregre.basis.fr(), superficie β(s,t) en base tensorial
#       con penalizaciones lambda.s y lambda.t.
# Se ajustan las tres con las mismas curvas √P, el mismo FPCA de la estación (M) y la misma
# regla: parámetros (l, λ_s, λ_t) elegidos con el RMSE de VALIDACIÓN. fregre.basis.fr no
# admite covariables escalares, así que la comparación justa es sin altitud ni ONI; la del
# taller se hace con y sin ellas. Se comparan en prueba con las mismas métricas y bootstrap.
iE_c <- which(ent); iV_c <- which(val)
H_c <- S_todos[, 1:M, drop = FALSE]
curva_de_puntajes <- function(Sc) sweep(Sc %*% t(Phi_e %*% psi_c[, 1:M, drop = FALSE]), 2, drop(Phi_e %*% mu_c), "+")
rmse_val_c <- function(R) mean(rmse_l2(a_mm(R[iV_c, , drop = FALSE]), iV_c))
extras_c <- cbind((pares$altitud - alt_m) / alt_s, pares$ONI_anual)

# (a) FPCA del satélite (solo entrenamiento) y regresión puntajes sobre puntajes
ev_s <- eigen(W_s, symmetric = TRUE)
W12_s <- ev_s$vectors %*% diag(sqrt(ev_s$values)) %*% t(ev_s$vectors); Wm12_s <- ev_s$vectors %*% diag(1 / sqrt(ev_s$values)) %*% t(ev_s$vectors)
mu_x <- colMeans(rs$coef[ent, ])
eg_x <- eigen(W12_s %*% (crossprod(sweep(rs$coef[ent, ], 2, mu_x)) / (sum(ent) - 1)) %*% W12_s, symmetric = TRUE)
psi_x <- Wm12_s %*% eg_x$vectors
Xi <- sweep(rs$coef, 2, mu_x) %*% W_s %*% psi_x                      # puntajes de X en todas las curvas
ajuste_taller <- function(l, con_extras) {
  Z <- cbind(1, Xi[, seq_len(l), drop = FALSE], if (con_extras) extras_c)
  B <- solve(crossprod(Z[iE_c, ]), crossprod(Z[iE_c, ], H_c[iE_c, ]))
  list(B = B, Z = Z, R = curva_de_puntajes(Z %*% B))
}
sel_taller <- rbindlist(lapply(c(TRUE, FALSE), function(ce) rbindlist(lapply(seq_len(K_s), function(l)
  data.table(metodo = "taller PC", con_altitud_ONI = ce, l = l, rmse_validacion = rmse_val_c(ajuste_taller(l, ce)$R))))))
l_con <- sel_taller[con_altitud_ONI == TRUE][which.min(rmse_validacion), l]
l_sin <- sel_taller[con_altitud_ONI == FALSE][which.min(rmse_validacion), l]
# El taller dice "lm() solo para verificar": la B de mínimos cuadrados debe ser la de lm()
aj_t <- ajuste_taller(l_con, TRUE)
dif_lm_taller <- max(abs(aj_t$B - coef(lm(H_c[iE_c, ] ~ aj_t$Z[iE_c, -1]))))

# Nuestra regresión sin altitud ni ONI (λ elegido igual, en validación)
ajuste_nuestra_sin <- function(l) {
  Z <- cbind(1, Zx); P <- matrix(0, ncol(Z), ncol(Z)); P[-1, -1] <- R_s
  curva_de_puntajes(Z %*% solve(crossprod(Z[iE_c, ]) + l * P, crossprod(Z[iE_c, ], H_c[iE_c, ])))
}
sel_nuestra <- data.table(metodo = "nuestra", con_altitud_ONI = FALSE, lambda = lambdas, rmse_validacion = sapply(lambdas, function(l) rmse_val_c(ajuste_nuestra_sin(l))))
LAMBDA_SIN <- sel_nuestra[which.min(rmse_validacion), lambda]

# (b) fda.usc::fregre.basis.fr con las mismas bases B-spline (nodos equiespaciados, como las nuestras)
X_fd <- fdata(rs$coef %*% t(Phi_s), argvals = T_DIA); Y_fd <- fdata(re$coef %*% t(Phi_e), argvals = T_DIA)
base_s_fda <- create.bspline.basis(RANGO, K_s); base_t_fda <- create.bspline.basis(RANGO, K_e)
ajuste_fr <- function(ls, lt) fregre.basis.fr(X_fd[iE_c], Y_fd[iE_c], basis.s = base_s_fda, basis.t = base_t_fda, lambda.s = ls, lambda.t = lt)
rej_fr <- CJ(lambda_s = 10^seq(-2, 14, 2), lambda_t = 10^c(2, 6, 10, 14))
sel_fr <- rej_fr[, .(metodo = "fregre.basis.fr", con_altitud_ONI = FALSE, rmse_validacion = rmse_val_c(predict(ajuste_fr(lambda_s, lambda_t), X_fd)$data)), by = .(lambda_s, lambda_t)]
mejor_fr <- sel_fr[which.min(rmse_validacion)]
fit_fr <- ajuste_fr(mejor_fr$lambda_s, mejor_fr$lambda_t)
stopifnot(mejor_fr$lambda_s < max(rej_fr$lambda_s))   # el óptimo no queda en el borde de la rejilla
tabla_f(rbindlist(list(sel_taller, sel_nuestra, sel_fr), fill = TRUE), "F5_contraste_seleccion_validacion")

modelos_clase <- list(
  `Satélite sin corregir` = X_mm[i_pru, ],
  `Nuestra regresión (con altitud y ONI)` = a_mm(pred_puntajes),
  `Taller PC (con altitud y ONI)` = a_mm(aj_t$R[i_pru, ]),
  `Nuestra regresión (sin altitud ni ONI)` = a_mm(ajuste_nuestra_sin(LAMBDA_SIN)[i_pru, ]),
  `Taller PC (sin altitud ni ONI)` = a_mm(ajuste_taller(l_sin, FALSE)$R[i_pru, ]),
  `fda.usc fregre.basis.fr (sin altitud ni ONI)` = a_mm(predict(fit_fr, X_fd[i_pru])$data))
mc_clase <- rbindlist(lapply(names(modelos_clase), function(nm) metricas_de(modelos_clase[[nm]], i_pru, nm)))
contraste <- mc_clase[, resumir(.SD), by = modelo]
contraste[, parametros := c("", sprintf("M = %d, λ = %.1e", M, LAMBDA), sprintf("M = %d, l = %d", M, l_con),
                            sprintf("M = %d, λ = %.1e", M, LAMBDA_SIN), sprintf("M = %d, l = %d", M, l_sin),
                            sprintf("λs = %.0e, λt = %.0e", mejor_fr$lambda_s, mejor_fr$lambda_t))]
print(tabla_f(contraste, "F5_contraste_metodos_clase"))

# Diferencias pareadas por estación (IC 95 % bootstrap y Wilcoxon), mismas reglas que 5.5
pe_c <- mc_clase[, .(l2_rel = mean(l2_rel), kge = mean(kge)), by = .(modelo, codigo)]
dif_par <- function(a, b, v) {
  x <- pe_c[modelo == a][order(codigo)][[v]]; y <- pe_c[modelo == b][order(codigo)][[v]]; d <- x - y
  q <- quantile(replicate(N_BOOT, mean(sample(d, replace = TRUE))), c(.025, .975))
  data.table(comparacion = paste(a, "−", b), metrica = v, diferencia = round(mean(d), 4), ic_inf = round(q[1], 4), ic_sup = round(q[2], 4),
             p_wilcoxon = signif(wilcox.test(x, y, paired = TRUE)$p.value, 2), ic_excluye_0 = q[1] > 0 | q[2] < 0)
}
pares_clase <- list(c("Nuestra regresión (con altitud y ONI)", "Taller PC (con altitud y ONI)"),
                    c("Nuestra regresión (sin altitud ni ONI)", "Taller PC (sin altitud ni ONI)"),
                    c("Nuestra regresión (sin altitud ni ONI)", "fda.usc fregre.basis.fr (sin altitud ni ONI)"),
                    c("Taller PC (sin altitud ni ONI)", "fda.usc fregre.basis.fr (sin altitud ni ONI)"),
                    c("Nuestra regresión (con altitud y ONI)", "Nuestra regresión (sin altitud ni ONI)"))
contraste_par <- rbindlist(lapply(pares_clase, function(p) rbind(dif_par(p[1], p[2], "l2_rel"), dif_par(p[1], p[2], "kge"))))
print(tabla_f(contraste_par, "F5_contraste_pareado"))

# QUÉ DICEN LOS DATOS (5.5b)
#   Los tres métodos de clase llegan a lo mismo en error: en prueba el L² relativo queda
#   entre 0,43 y 0,44 en todos (satélite sin corregir 0,51), y las diferencias de L² entre
#   nuestra regresión y la del taller no excluyen el 0. Donde sí hay diferencia es en el
#   KGE: la nuestra lo sube unos 0,03 frente a la del taller y frente a fregre.basis.fr
#   (IC que excluyen el 0), porque la penalización de rugosidad conserva mejor la forma y la
#   amplitud de la curva que truncar X en l = 4 componentes. Conclusión: nuestra regresión
#   no es un método aparte; es el modelo integral de clase, y frente a sus dos versiones de
#   clase da el mismo error y mejor KGE. Altitud y ONI no cambian la prueba (IC con el 0;
#   2001-2003 no tiene años Niña), pero en validación bajan el RMSE de 2,42 a 2,29 mm/día,
#   por eso se dejan (y la red los usa).

# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ┃ RÚBRICA 4.1 / 4.4 — Sensibilidad: sin 1985, 1988 y 2000
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
titulo("5.6 Sensibilidad: ¿qué pasa si se quitan 1985, 1988 y 2000?")
# No se quitan (D-47). Esta sección solo mide cuánto cambia la prueba si el modelo
# NO viera esos años: mismo FPCA y mismo M; lambda se vuelve a elegir en validación sin 2000.
ent_todos <- ent; val_todos <- val
ent <- ent_todos & !pares$marcado; val <- val_todos & !pares$marcado
rmse_val_s <- sapply(lambdas, function(l) mean(rmse_l2(a_mm(predecir_raiz(ajustar_B(l, M), which(val))), which(val))))
LAMBDA_S <- lambdas[which.min(rmse_val_s)]; B_sens <- ajustar_B(LAMBDA_S, M)
ent <- ent_todos; val <- val_todos
mc_sens <- metricas_de(a_mm(predecir_raiz(B_sens, i_pru)), i_pru, "Regresión sin 1985, 1988 y 2000")
sens <- rbind(metricas_curva[modelo == "Regresión funcional"][, modelo := "Regresión con todos los años (la que se usa)"], mc_sens)
sensibilidad <- sens[, .(pares_entrenamiento = 0L, lambda = 0, L2_relativo = mean(l2_rel), RMSE_L2 = mean(rmse), KGE = mean(kge), sesgo_anual_pct = mean(sesgo_pct),
                         error_pico_1 = mean(err_pico_1), error_pico_2 = mean(err_pico_2)), by = modelo]
sensibilidad[, `:=`(pares_entrenamiento = c(sum(ent_todos), sum(ent_todos & !pares$marcado)), lambda = c(LAMBDA, LAMBDA_S))]
sens_anio <- dcast(sens[, .(L2 = round(mean(l2_rel), 3), sesgo = round(mean(sesgo_pct), 1)), by = .(modelo, anio)], anio ~ modelo, value.var = c("L2", "sesgo"))
ps <- merge(sens[modelo %like% "todos", .(a = mean(l2_rel)), by = codigo], mc_sens[, .(s = mean(l2_rel)), by = codigo], by = "codigo")
sensibilidad[, c("L2_relativo", "RMSE_L2", "KGE", "sesgo_anual_pct", "error_pico_1", "error_pico_2") := lapply(.SD, round, 3), .SDcols = c("L2_relativo", "RMSE_L2", "KGE", "sesgo_anual_pct", "error_pico_1", "error_pico_2")]
print(tabla_f(sensibilidad, "F5_sensibilidad_anios_marcados"))
print(sens_anio)
cat(sprintf("Con todos los años el error L² relativo en prueba es menor en %d de %d estaciones (Wilcoxon pareado p = %.3g)\n",
            sum(ps$a < ps$s), nrow(ps), wilcox.test(ps$a, ps$s, paired = TRUE)$p.value))

R2_t <- function(P, Yo) 1 - colSums((Yo - P)^2) / colSums(sweep(Yo, 2, colMeans(Yo))^2)
grafico_base("F5_prueba", function() {
  par(mfrow = c(2, 3))
  ej <- i_pru[order(metricas_curva[modelo == "Regresión funcional", rmse])[c(round(length(i_pru) * c(.1, .5, .9)))]]
  for (k in ej) {
    j <- match(k, i_pru)
    matplot(T_DIA, cbind(Y_mm[k, ], modelos[[1]][j, ], modelos[[2]][j, ]), type = "l", lty = c(1, 2, 1), lwd = 2, col = c(1, "darkorange", 4),
            xlab = "Día del año", ylab = "mm/día", main = paste(pares$nombre[k], pares$anio[k]))
    legend("topleft", c("Estación", "Satélite", "Regresión"), col = c(1, "darkorange", 4), lty = c(1, 2, 1), lwd = 2, bty = "n", cex = .8)
  }
  plot(T_DIA, R2_t(modelos[[2]], Y_mm[i_pru, ]), type = "l", lwd = 2, col = 4, ylim = c(-1, 1), xlab = "Día del año", ylab = expression(R^2*(t)),
       main = "Varianza explicada por día (prueba)")
  lines(T_DIA, R2_t(modelos[[1]], Y_mm[i_pru, ]), lwd = 2, col = "darkorange"); abline(h = 0, lty = 3)
  legend("bottomleft", c("Regresión", "Satélite sin corregir"), col = c(4, "darkorange"), lwd = 2, bty = "n")
  boxplot(l2_rel ~ modelo, metricas_curva, ylab = "Error L² relativo", xlab = "", main = "Error por curva (prueba)", cex.axis = .8)
  boxplot(sesgo_pct ~ modelo, metricas_curva, ylab = "% del acumulado anual", xlab = "", main = "Sesgo anual (prueba)", cex.axis = .8); abline(h = 0, lty = 3)
}, ancho = 15, alto = 8)

# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ┃ RÚBRICA 4.2 — Exportación de insumos para la red en Python
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
titulo("5.6b Inferencia y diagnóstico como en la tesis del profesor")
# (a) Prueba F por permutación del aporte del satélite, a la manera de Fperm.fd (clase, semana 15):
#     F(t) = Σ_i (ŷ_i(t) - ȳ̂(t))² / Σ_i (y_i(t) - ŷ_i(t))², en la escala raíz y en entrenamiento.
#     H0: la curva del satélite no aporta nada más allá de la altitud y el ONI. Se permutan las
#     curvas del satélite ENTRE estaciones del mismo año (se conservan el año y su fase ENSO) y
#     se reajusta B con el mismo FPCA, M y lambda. Se reporta el criterio del máximo y el integrado.
Y_raiz <- re$coef %*% t(Phi_e); fil_ent <- which(ent)
F_t_de <- function(P, filas) colSums(sweep(P, 2, colMeans(P))^2) / colSums((Y_raiz[filas, ] - P)^2)
F_obs <- F_t_de(predecir_raiz(B, fil_ent), fil_ent)
anio_ent <- pares$anio[fil_ent]
perm_dentro_anio <- function() { nuevo <- fil_ent; for (a in unique(anio_ent)) { k <- which(anio_ent == a); nuevo[k] <- fil_ent[k][sample.int(length(k))] }; nuevo }
B_FPERM <- 300L
set.seed(20261120)
F_perm <- replicate(B_FPERM, { zxp <- Zx; zxp[fil_ent, ] <- Zx[perm_dentro_anio(), ]
  F_t_de(predecir_raiz(ajustar_B(LAMBDA, M, zx = zxp), fil_ent, zx = zxp), fil_ent) })
c_alfa_F <- apply(F_perm, 1, quantile, .95); c_max_F <- quantile(apply(F_perm, 2, max), .95)
fperm <- data.table(estadistico = c("F máximo", "F integrado (media sobre el año)"), observado = round(c(max(F_obs), mean(F_obs)), 3),
                    critico_95 = round(c(c_max_F, quantile(colMeans(F_perm), .95)), 3),
                    p_permutacion = signif(c((1 + sum(apply(F_perm, 2, max) >= max(F_obs))) / (B_FPERM + 1), (1 + sum(colMeans(F_perm) >= mean(F_obs))) / (B_FPERM + 1)), 2),
                    dias_sobre_el_critico_puntual = c(sum(F_obs > c_alfa_F), NA), permutaciones = B_FPERM)
print(tabla_f(fperm, "F5_permutacion_satelite"))

# (b) Diagnóstico de residuos (tesis, figuras 9.3 y 9.4; clase, correlación de los residuos):
#     residuos estudentizados por día con la DE de entrenamiento, residuos frente a predichos y
#     correlación de los residuos entre días.
res_ent <- Y_raiz[fil_ent, ] - predecir_raiz(B, fil_ent); sig_t <- apply(res_ent, 2, sd)
res_pru <- Y_raiz[i_pru, ] - pred_puntajes; stud_pru <- sweep(res_pru, 2, sig_t, "/")
cor_res <- cor(res_ent)
diag_res <- data.table(medida = c("Residuos estudentizados de prueba fuera de ±2 (%)", "Residuos estudentizados de prueba fuera de ±3 (%)",
                                  "Residuo medio de prueba (raíz de mm/día)", "Correlación de los residuos a 1 día", "Correlación de los residuos a 30 días",
                                  "Correlación de los residuos a 90 días", "Correlación residuo-predicho en prueba"),
                       valor = round(c(100 * mean(abs(stud_pru) > 2), 100 * mean(abs(stud_pru) > 3), mean(res_pru),
                                       mean(cor_res[cbind(1:364, 2:365)]), mean(cor_res[cbind(1:335, 31:365)]), mean(cor_res[cbind(1:275, 91:365)]),
                                       cor(as.vector(res_pru), as.vector(pred_puntajes))), 3))
print(tabla_f(diag_res, "F5_diagnostico_residuos"))

# (c) Traza-variograma de los residuos (tesis, capítulo 8, Giraldo 2009): γ(h) = ½ E||e_i - e_j||²
#     para pares de estaciones del mismo año a distancia h. Si crece con h, los residuos de
#     estaciones cercanas se parecen: hay correlación espacial que el modelo no recoge. Se ajusta
#     el modelo exponencial de la tesis, γ(h) = c0 + σ²(1 - exp(-h/φ)), por mínimos cuadrados ponderados.
coord_par <- est[match(pares$codigo[fil_ent], codigo), .(lat, lon)]
tv <- rbindlist(lapply(sort(unique(anio_ent)), function(a) {
  k <- which(anio_ent == a); E <- res_ent[k, , drop = FALSE]
  G <- E %*% (t(E) * w_trap); nn <- diag(G); D2 <- outer(nn, nn, "+") - 2 * G
  H <- outer(seq_along(k), seq_along(k), function(i, j) haversine(coord_par$lat[k[i]], coord_par$lon[k[i]], coord_par$lat[k[j]], coord_par$lon[k[j]]))
  u <- which(upper.tri(D2), arr.ind = TRUE); data.table(anio = a, h = H[u], g = D2[u] / 2)
}))
tv[, bin := cut(h, c(seq(0, 100, by = 10), 150, 250), include.lowest = TRUE)]
tv_bin <- tv[, .(h = mean(h), gamma = mean(g), pares = .N), by = bin][order(h)]
obj_tv <- function(p) with(tv_bin, sum(pares * (gamma - (p[1] + p[2] * (1 - exp(-h / p[3]))))^2))
fit_tv <- optim(c(min(tv_bin$gamma), diff(range(tv_bin$gamma)), 30), obj_tv, method = "L-BFGS-B", lower = c(0, 0, 1), upper = c(Inf, Inf, 500))$par
tabla_tv <- tv_bin[, .(bin = as.character(bin), h_km = round(h, 1), traza_variograma = round(gamma, 2), pares)]
print(tabla_f(tabla_tv, "F5_traza_variograma_residuos"))
cat(sprintf("Traza-variograma exponencial de los residuos: pepita c0 = %.1f, meseta parcial σ² = %.1f, φ = %.1f km (rango práctico %.0f km); la pepita es %.0f %% de la meseta\n",
            fit_tv[1], fit_tv[2], fit_tv[3], 3 * fit_tv[3], 100 * fit_tv[1] / (fit_tv[1] + fit_tv[2])))
tabla_f(data.table(parametro = c("pepita_c0", "meseta_parcial", "phi_km", "rango_practico_km", "pct_pepita"),
                   valor = round(c(fit_tv, 3 * fit_tv[3], 100 * fit_tv[1] / sum(fit_tv[1:2])), 2)), "F5_traza_variograma_ajuste")

grafico_base("F5_diagnostico", function() {
  par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))
  smoothScatter(as.vector(pred_puntajes), as.vector(res_pru), nrpoints = 0, xlab = "Predicho (raíz de mm/día)", ylab = "Residuo",
                main = "Residuos frente a predichos (prueba)"); abline(h = 0, lty = 2, col = "red")
  q <- apply(stud_pru, 2, quantile, c(.05, .25, .5, .75, .95))
  matplot(T_DIA, t(q), type = "l", lty = c(3, 2, 1, 2, 3), lwd = c(1, 1, 2, 1, 1), col = 1, ylim = c(-4, 4), xlab = "Día del año",
          ylab = "Residuo estudentizado", main = "Residuos estudentizados por día (prueba): cuantiles 5-95 %"); abline(h = c(-2, 2), lty = 3, col = "red")
  contour(T_DIA, T_DIA, cor_res, nlevels = 10, xlab = "Día s", ylab = "Día t", main = "Correlación de los residuos (entrenamiento)")
  with(tv_bin, plot(h, gamma, pch = 19, cex = sqrt(pares) / max(sqrt(pares)) * 2 + .5, ylim = c(0, max(gamma) * 1.1),
                    xlab = "Distancia entre estaciones (km)", ylab = "Traza-variograma", main = "Traza-variograma de los residuos"))
  hh <- seq(0, 250, length.out = 200); lines(hh, fit_tv[1] + fit_tv[2] * (1 - exp(-hh / fit_tv[3])), col = "red", lwd = 2)
}, ancho = 12, alto = 9)

grafico_base("F5_permutacion_satelite", function() {
  matplot(T_DIA, F_perm, type = "l", lty = 1, col = "grey85", ylim = c(0, max(F_obs, F_perm) * 1.05), xlab = "Día del año", ylab = "F(t)",
          main = sprintf("Aporte del satélite: F(t) observado frente a %d permutaciones dentro de cada año", B_FPERM))
  lines(T_DIA, F_obs, lwd = 2); lines(T_DIA, c_alfa_F, lty = 2, col = "steelblue", lwd = 2); abline(h = c_max_F, lty = 2, col = "firebrick", lwd = 2)
  legend("topright", c("Observado", "Permutaciones", "Crítico puntual 95 %", "Crítico del máximo 95 %"), col = c(1, "grey70", "steelblue", "firebrick"),
         lty = c(1, 1, 2, 2), lwd = 2, bty = "n", cex = .85)
}, ancho = 10, alto = 5)

# (d) Curvas observadas y ajustadas de 20 estaciones de un año de prueba (tesis, figura 10.3).
#     Se toma 2002 (año Niño) y 20 estaciones repartidas por cuantiles del error.
k02 <- i_pru[pares$anio[i_pru] == 2002]
err02 <- metricas_curva[modelo == "Regresión funcional" & anio == 2002, l2_rel]
sel20 <- k02[order(err02)][round(seq(1, length(k02), length.out = 20))]
grafico_base("F5_curvas_20_estaciones", function() {
  par(mfrow = c(4, 5), mar = c(2.5, 2.5, 2.5, .5), mgp = c(1.5, .5, 0))
  for (k in sel20) {
    j <- match(k, i_pru)
    matplot(T_DIA, cbind(Y_mm[k, ], modelos[[1]][j, ], modelos[[2]][j, ]), type = "l", lty = c(1, 2, 1), lwd = c(2, 1, 2),
            col = c(1, "darkorange", 4), xlab = "", ylab = "mm/día", main = sprintf("%s (L² %.2f)", substr(pares$nombre[k], 1, 18), err02[match(k, k02)]), cex.main = .85)
  }
}, ancho = 15, alto = 11)

# (e) Mapa del error por estación en prueba (tesis, figuras 10.4 a 10.6).
mapa_err <- merge(metricas_curva[modelo != "Solo altitud y ONI", .(l2 = mean(l2_rel)), by = .(codigo, modelo)][, dcast(.SD, codigo ~ modelo, value.var = "l2")],
                  est[, .(codigo, nombre, lat, lon, vertiente)], by = "codigo")
setnames(mapa_err, c("Regresión funcional", "Satélite sin corregir"), c("l2_reg", "l2_sat"))
mapa_err[, mejora := l2_reg - l2_sat]
tabla_f(mapa_err[order(mejora)], "F5_error_por_estacion_prueba")
cat(sprintf("Prueba, por estación: la regresión baja el error L² relativo en %d de %d estaciones; error medio %.2f en el Cauca y %.2f en el Pacífico\n",
            sum(mapa_err$mejora < 0), nrow(mapa_err), mapa_err[vertiente != "Pacífico", mean(l2_reg)], mapa_err[vertiente == "Pacífico", mean(l2_reg)]))
valle <- sf::st_read("datos_web/limites/deptos_valle_y_vecinos_DANE_2020.gpkg", quiet = TRUE)
grafico_f(ggplot() + geom_sf(data = valle, fill = "grey97", colour = "grey60") +
            geom_point(data = mapa_err, aes(lon, lat, colour = mejora, size = l2_reg)) +
            scale_colour_gradient2(low = "steelblue", mid = "grey80", high = "firebrick", midpoint = 0) + scale_size(range = c(1, 4.5)) +
            coord_sf(xlim = range(mapa_err$lon) + c(-.15, .15), ylim = range(mapa_err$lat) + c(-.15, .15)) +
            labs(x = NULL, y = NULL, colour = "Regresión - satélite\n(error L², < 0 = mejora)", size = "Error L² relativo\nde la regresión",
                 title = "Error por estación en prueba (2001-2003)"),
          "F5_mapa_error", ancho = 8, alto = 7.5)

titulo("5.7 Exportación para la red funcional (Python)")
dir_py <- file.path(SALIDA_F, "para_python"); dir.create(dir_py, showWarnings = FALSE)
fwrite(pares[, .(codigo, anio, particion, usar, celda, altitud, ONI_anual, fase, vertiente)], file.path(dir_py, "pares.csv"))
fwrite(as.data.table(re$coef), file.path(dir_py, "coef_estacion_raiz.csv")); fwrite(as.data.table(rs$coef), file.path(dir_py, "coef_satelite_raiz.csv"))
fwrite(as.data.table(o$coef_est[pares$fila_est, ]), file.path(dir_py, "coef_estacion_mm.csv")); fwrite(as.data.table(o$coef_sat[pares$fila_sat, ]), file.path(dir_py, "coef_satelite_mm.csv"))
for (nm in c("W_e", "W_s", "R_s")) fwrite(as.data.table(get(nm)), file.path(dir_py, paste0(nm, ".csv")), col.names = FALSE)
fwrite(data.table(mu = mu_c, psi_c), file.path(dir_py, "fpca_estaciones.csv"))
fwrite(as.data.table(B), file.path(dir_py, "coeficientes_regresion.csv"), col.names = FALSE)
writeLines(c(sprintf('{"rango": [1, 365], "orden": %d, "K_estaciones": %d, "K_satelite": %d, "M": %d, "lambda": %g,', ORDEN, K_e, K_s, M, LAMBDA),
             sprintf('"nodos_estaciones": [%s], "nodos_satelite": [%s],', paste(nodos(K_e), collapse = ", "), paste(nodos(K_s), collapse = ", ")),
             sprintf('"altitud_media": %.6f, "altitud_de": %.6f, "escala": "raiz de la curva en mm/dia recortada en 0"}', alt_m, alt_s)),
           file.path(dir_py, "especificacion.json"))
cat("Archivos para Python en", dir_py, ":", paste(list.files(dir_py), collapse = ", "), "\n")

# QUÉ DICEN LOS DATOS
#   Escala: al recortar en 0 solo cambia el 2 % de los días de las estaciones y
#   casi ningún día del satélite; la raíz se reproyecta en la misma base con un
#   error de 0,09 (estaciones) y 0,03 (satélite) en raíz de mm/día.
#   Orden de la raíz y el suavizado: suavizar en mm y luego sacar la raíz conserva
#   el acumulado anual (cociente mediano 1,000 frente al dato diario). Sacar la
#   raíz al dato diario y luego suavizar lo reduce a un tercio (0,33), porque la
#   lluvia diaria es intermitente y el cuadrado de la media de raíces no es la
#   media. Como el sesgo del acumulado anual es una métrica de la rúbrica, se
#   usa el primer camino.
#   Datos del modelo: 3 299 pares con todos los años (2 257 de entrenamiento, 672
#   de validación, 370 de prueba). Univalle entra con 996 m (D-46).
#   FPCA (solo entrenamiento): la primera componente explica 62 % (el nivel de
#   lluvia de cada sitio), la segunda 13 % y la tercera 7 %. M y lambda se eligen
#   juntos con validación: el error deja de bajar desde M = 6 (94,2 % de la
#   varianza; el enunciado pide M entre 4 y 6), con lambda = 3,2e8 (10 grados de
#   libertad efectivos de 14).
#   Equivalencia: la predicción con la integral de beta(s,t) por cuadratura y la
#   predicción por puntajes FPCA difieren en 4e-15 (criterio 1e-4). Es la
#   identidad beta(s,t) = sum_m beta_m(s) psi_m(t) que debe reproducir la cabeza
#   lineal de la red.
#   Interpretación: lejos de los bordes, la integral de beta sobre s va de 0,62
#   (diciembre) a 1,03 (1-jun): el modelo traslada el nivel del satélite casi
#   completo a mitad de año y lo amortigua al final. El ONI resta lluvia todo el
#   año y sobre todo a mitad de año (-0,24 raíz(mm/día) por unidad el 2-jul): El Niño profundiza la sequía de
#   junio-julio. La altitud pesa poco una vez que se conoce el píxel.
#   Prueba (2001-2003, 370 curvas de 137 estaciones, IC 95 % por bootstrap de
#   estaciones con 2 000 remuestreos). Frente al satélite sin corregir, la
#   regresión baja el error L² relativo de 0,51 a 0,44, sube el KGE de 0,38 a
#   0,47, baja el sesgo del acumulado anual de +26 % a +7 % y acerca el pico de
#   la primera temporada de 23 a 18 días. Pareado por estación, el IC de la
#   diferencia excluye el cero y el Wilcoxon da p < 0,001 en cinco métricas; en
#   el pico de la segunda temporada la ganancia es de medio día y no es
#   significativa (IC con el cero, p = 0,08). Mejora en 65 % de las curvas. Sin
#   el satélite (solo altitud y ONI) el error queda como el del satélite crudo
#   (0,51) y lejos del 0,44 de la regresión: el satélite sí aporta.
#   Por fase ENSO: la prueba solo tiene años neutros (2001, 2003) y un Niño
#   (2002); no tiene Niña. En 2002 el modelo corrige de más (sesgo -14 % frente a
#   +28 % del satélite). En validación (orientativa) el Niño de 1997 es donde más
#   gana (L² 0,61 a 0,46). En los años Niña 1999-2000 casi empata con el
#   satélite (0,41 a 0,40) y el sesgo empeora (de -0,5 % a -11 %): el modelo baja
#   el nivel del satélite como en un año normal, y el término del ONI no alcanza
#   a devolver lo que el satélite ya perdió en La Niña. Es el punto débil de la
#   referencia lineal y lo que una cabeza no lineal podría mejorar.
#   Sensibilidad: si el modelo no viera 1985, 1988 y 2000, la prueba empeora (L²
#   0,439 frente a 0,434; sesgo 8,9 % frente a 7,4 %) en 79 de 137 estaciones
#   (Wilcoxon p = 0,003). Quitar esos años no ayuda: confirma con
#   datos la indicación del profesor (D-47).
#   La reconstrucción en mm es no negativa y solo 0,02 % de los valores se
#   recorta en 0.
#   Esta regresión es la referencia lineal. La red funcional (codificador con W,
#   decodificador FPCA con estas psi_m y cabezas lineal, KAN y MLP) sigue en el
#   archivo de Python con los insumos de para_python/.

regresion <- list(M = M, lambda = LAMBDA, B = B, mu_c = mu_c, psi_c = psi_c, var_fpca = var_fpca, BETA = BETA, metricas = metricas,
             metricas_curva = metricas_curva, dif_equivalencia = dif_equiv, alt_m = alt_m, alt_s = alt_s)

# =============================================================================
# 6. PRUEBAS Y CIERRE
# =============================================================================
# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ┃ RÚBRICA todas — 13 pruebas automáticas de la parte B
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
titulo("6. Pruebas de la parte funcional")
# Prueba de no fuga. (a) Se reemplaza por ruido TODO lo que no es entrenamiento (curvas
# de estación, píxel, altitud y ONI de validación y prueba) y se reajusta con las mismas
# funciones: media, autovalores, escalas de altitud y B tienen que salir idénticos.
# (b) Se altera solo la prueba y se vuelve a elegir M y lambda: tienen que salir iguales.
# (c) K se eligió con ent_est / ent_sat: se revisa que esas filas sean solo de 1981-1995.
set.seed(7)
ruido_filas <- function(A, filas) { A[filas, ] <- matrix(rnorm(sum(filas) * ncol(A), sd = 3 * sd(A)), sum(filas)); A }
alterar <- function(filas) {
  coef <- ruido_filas(re$coef, filas)
  list(coef = coef, Y = abs(ruido_filas(Y_mm, filas)), zx = ruido_filas(Zx, filas),
       alt = replace(pares$altitud, filas, runif(sum(filas), 0, 4500)), oni = replace(pares$ONI_anual, filas, rnorm(sum(filas), sd = 2)))
}
d <- alterar(!ent)
fp_r <- fpca_ent(d$coef, ent); am_r <- mean(d$alt[ent]); ds_r <- sd(d$alt[ent])
S_r <- sweep(d$coef, 2, fp_r$mu) %*% W_e %*% fp_r$psi
B_r <- ajustar_B(LAMBDA, M, S = S_r, zx = d$zx, alt = d$alt, oni = d$oni, am = am_r, ds = ds_r)
no_fuga_ent <- c(media = max(abs(fp_r$mu - mu_c)), autovalores = max(abs(fp_r$values - var_fpca$autovalor)),
                 escalas = max(abs(c(am_r, ds_r) - c(alt_m, alt_s))), B = max(abs(B_r - B)))
d <- alterar(pru)
S_p <- sweep(d$coef, 2, mu_c) %*% W_e %*% psi_todas
el_p <- elegir_M_lambda(tabla_seleccion(S = S_p, Y = d$Y, zx = d$zx, alt = d$alt, oni = d$oni))
no_fuga_sel <- el_p$M == M && el_p$lambda == LAMBDA
no_fuga_K <- all(id_est$anio[ent_est] <= 1995) && all(id_sat$anio[ent_sat] <= 1995)
cat(sprintf("No fuga: con validación y prueba reemplazadas por ruido, diferencia máxima en media %.1e, autovalores %.1e, escalas %.1e y B %.1e
",
            no_fuga_ent[1], no_fuga_ent[2], no_fuga_ent[3], no_fuga_ent[4]))
cat(sprintf("No fuga: con la prueba alterada se vuelve a elegir M = %d y lambda = %s (original M = %d, lambda = %s)
",
            el_p$M, format(el_p$lambda, scientific = TRUE), M, format(LAMBDA, scientific = TRUE)))
rm(d)
pruebas <- list(
  "El emparejamiento es el mismo del archivo de datos" = identical(pares[order(fila_est), fila_sat], dp$pares[order(fila_est), fila_sat]) && identical(pares[order(fila_est), codigo], dp$pares[order(fila_est), codigo]),
  "Ningún año queda fuera del modelo" = all(pares$usar) && setequal(unique(pares$anio[ent]), 1981:1995) && setequal(unique(pares$anio[val]), 1996:2000) && setequal(unique(pares$anio[pru]), 2001:2003),
  "La base propia coincide con splines::splineDesign" = max(abs(base_bspline(T_DIA, 11) - ref_splines)) < 1e-12,
  "S simétrica, idempotente y tr(S) = K en las dos fuentes" = all(controles[, .(simetria, idempotencia, traza_menos_K)] < 1e-8),
  "W integra el dominio: 1'W1 = 364" = abs(sum(W_e) - 364) < 1e-8 && abs(sum(W_s) - 364) < 1e-8,
  "La media paso a paso coincide con fda::mean.fd" = all(verif_media$vs_fda_mean.fd < 1e-10),
  "Las atípicas coinciden con fda::fbplot" = all(verif_fb),
  "Funciones propias ortonormales con W" = max(abs(t(psi_c) %*% W_e %*% psi_c - diag(M))) < 1e-10,
  "Sin fuga: K, FPCA, escalas y B no cambian si se alteran validación y prueba; M y lambda no cambian si se altera la prueba" =
    no_fuga_K && all(no_fuga_ent < 1e-10) && no_fuga_sel,
  "Equivalencia beta(s,t) vs puntajes FPCA (<= 1e-4)" = dif_equiv <= 1e-4,
  "La reconstrucción en mm no es negativa" = min(modelos[["Regresión funcional"]]) >= 0,
  "Bootstrap con al menos 1000 remuestreos" = N_BOOT >= 1000,
  "Regresión concurrente: b1(t) por año = lm(), ecuaciones normales = mínimos cuadrados directos, sin fuga, rugosidad nula en una recta, b1 lineal con lambda enorme, reconstrucción no negativa y bootstrap suficiente" = all(concurrente$pruebas),
  "ANOVA funcional: SSE_n = (n-k) tr(γ̂) y, con k = 2, F(t) = t² (clase, ejercicio 1)" = all(chequeo_anova < 1e-8),
  "Bandas por bootstrap de estaciones y permutación del satélite con réplicas suficientes" = B_BOOT_EF >= 500 && B_FPERM >= 200 && B_PERM >= 1000,
  "Contraste con los métodos de clase: la B del taller es la de lm() y λ de fregre.basis.fr no queda en el borde" =
    dif_lm_taller < 1e-8 && mejor_fr$lambda_s < max(rej_fr$lambda_s)
)
res_pruebas <- data.table(prueba = names(pruebas), pasa = unlist(pruebas))
print(tabla_f(res_pruebas, "F6_pruebas"))
cat(sprintf("\nPruebas: %d de %d pasan\n", sum(res_pruebas$pasa), nrow(res_pruebas)))
stopifnot(res_pruebas$pasa)
saveRDS(list(construccion = o, eda = eda, comparacion = comparacion, regresion = regresion, pruebas = res_pruebas),
        file.path(SALIDA_F, "resultados_funcional.rds"), compress = "xz")
writeLines(capture.output(sessionInfo()), file.path(SALIDA_F, "sessionInfo.txt"))

# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ┃ RÚBRICA 4.3 / 4.4 — Ablación lineal / KAN / MLP: métricas, IC, Wilcoxon y fase ENSO
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
titulo("7. Evaluación de la red funcional (predicciones de lluvia_valle_red.py)")
# Python entrena la red (PyTorch, lo único permitido en ese apartado) y exporta los coeficientes
# √P predichos para las 3 299 curvas por cada cabeza y semilla. Aquí se evalúan con las MISMAS
# funciones de la sección 5.5 (métricas, bootstrap por estación, wilcox.test), así la regresión y
# las tres cabezas se miden exactamente igual y todo el análisis usa herramientas de clase.
# Orden de ejecución: este archivo → python lluvia_valle_red.py → este archivo otra vez.
ARCH_PRED <- file.path("salidas_python", "predicciones_coef_raiz.csv")
if (!file.exists(ARCH_PRED)) {
  cat("No está", ARCH_PRED, "→ corra python lluvia_valle_red.py y vuelva a correr este archivo para evaluar la red.\n")
} else {
  set.seed(20261114)
  pr <- fread(ARCH_PRED, colClasses = list(character = "codigo"))
  cols_c <- paste0("c", seq_len(K_e))
  # Control 1: las filas son los mismos pares, en el mismo orden
  ok_filas <- pr[, all(codigo == pares$codigo[fila] & anio == pares$anio[fila] & particion == pares$particion[fila])]
  # Control 2: la regresión de R reconstruida en Python da las mismas curvas que aquí
  rp <- pr[modelo == "regresión R (en Python)"][order(fila)]
  dif_py_R <- max(abs(as.matrix(rp[, ..cols_c]) %*% t(Phi_e) - predecir_raiz(B, rp$fila)))
  cat(sprintf("Regresión de R hecha en Python vs en R: diferencia máxima %.1e (√(mm/día))\n", dif_py_R))

  cabezas <- pr[modelo %in% c("lineal", "kan", "mlp")]
  i_vp <- which(val | pru)
  mc_red <- rbindlist(c(
    list(metricas_de(X_mm[i_vp, ], i_vp, "satélite sin corregir")[, semilla := 0L],
         metricas_de(a_mm(predecir_raiz(B, i_vp)), i_vp, "regresión (R)")[, semilla := 0L]),
    lapply(split(cabezas[fila %in% i_vp], by = c("modelo", "semilla")), function(d) {
      d <- d[order(fila)]
      metricas_de(a_mm(as.matrix(d[, ..cols_c]) %*% t(Phi_e)), d$fila, d$modelo[1])[, semilla := d$semilla[1]]
    })))
  tabla_f(mc_red, "F7_metricas_red_por_curva")
  vars_m <- c("l2_rel", "rmse", "kge", "sesgo_pct", "err_pico_1", "err_pico_2")

  # (a) Media ± DE entre semillas: en cada semilla, media de las medias por estación (como en 5.5)
  por_semilla <- mc_red[, lapply(.SD, mean), by = .(particion, modelo, semilla, codigo), .SDcols = vars_m][
    , lapply(.SD, mean), by = .(particion, modelo, semilla), .SDcols = vars_m]
  ablacion <- por_semilla[, c(.(semillas = .N), unlist(lapply(.SD, function(v) list(media = round(mean(v), 3), DE = round(sd(v), 3))), recursive = FALSE)),
                          by = .(particion, modelo), .SDcols = vars_m]
  print(tabla_f(ablacion, "F7_ablacion_media_DE_semillas"))

  # (b) IC 95 % bootstrap por estación, con las métricas de cada curva promediadas en las 3 semillas
  mc_prom <- mc_red[, lapply(.SD, mean), by = .(particion, modelo, codigo, anio, fase), .SDcols = vars_m]
  ic_red <- mc_prom[, resumir(.SD), by = .(particion, modelo)]
  print(tabla_f(ic_red, "F7_ablacion_IC_bootstrap"))

  # (c) Comparaciones pareadas por estación en PRUEBA: una cabeza supera a otra solo si el IC 95 %
  #     de la diferencia excluye el 0 (rúbrica 4.3); Wilcoxon pareado como apoyo.
  pe_red <- mc_prom[particion == "prueba", lapply(.SD, mean), by = .(modelo, codigo), .SDcols = vars_m]
  comparar_red <- function(a, b, v) {
    x <- pe_red[modelo == a][order(codigo)][[v]]; y <- pe_red[modelo == b][order(codigo)][[v]]
    if (v == "sesgo_pct") { x <- abs(x); y <- abs(y) }
    d <- x - y; q <- quantile(replicate(N_BOOT, mean(sample(d, replace = TRUE))), c(.025, .975))
    excl <- q[1] > 0 | q[2] < 0; mejor_a <- if (v == "kge") mean(d) > 0 else mean(d) < 0
    data.table(comparacion = paste(a, "−", b), metrica = if (v == "sesgo_pct") "|sesgo_pct|" else v,
               media_a = round(mean(x), 3), media_b = round(mean(y), 3), diferencia = round(mean(d), 4), ic_inf = round(q[1], 4), ic_sup = round(q[2], 4),
               p_wilcoxon = signif(wilcox.test(x, y, paired = TRUE, exact = FALSE)$p.value, 2), ic_excluye_0 = excl,
               gana = if (!excl) "ninguna (IC con 0)" else if (mejor_a) a else b)
  }
  pares_red <- list(c("kan", "lineal"), c("mlp", "lineal"), c("kan", "mlp"), c("lineal", "regresión (R)"),
                    c("lineal", "satélite sin corregir"), c("kan", "satélite sin corregir"), c("mlp", "satélite sin corregir"))
  comp_red <- rbindlist(lapply(pares_red, function(p) rbindlist(lapply(vars_m, function(v) comparar_red(p[1], p[2], v)))))
  print(tabla_f(comp_red, "F7_comparaciones_pareadas_prueba"))

  # (d) Por fase ENSO (rúbrica 4.4). La prueba no tiene años Niña: la validación se agrega como
  #     ORIENTATIVA (con ella se eligieron λ_w y la parada temprana, no los pesos).
  fase_red <- mc_prom[, resumir(.SD), by = .(particion, fase, modelo)]
  fase_red[, anios := sapply(seq_len(.N), function(k) paste(sort(unique(pares[particion == fase_red$particion[k] & fase == fase_red$fase[k], anio])), collapse = " "))]
  setcolorder(fase_red, c("particion", "fase", "anios")); setorder(fase_red, -particion, fase, modelo)
  print(tabla_f(fase_red, "F7_metricas_por_fase_ENSO"))

  pruebas_red <- list(
    "Las predicciones de Python corresponden a los mismos pares y en el mismo orden" = ok_filas,
    "La regresión de R hecha en Python coincide con la de R (< 1e-8)" = dif_py_R < 1e-8,
    "Tres semillas por cabeza y las 3 299 curvas en cada una" = cabezas[, .N, by = .(modelo, semilla)][, uniqueN(semilla) == 3 && all(N == nrow(pares)), by = modelo][, all(V1)],
    "Métricas finitas en todas las curvas de validación y prueba" = mc_red[, all(is.finite(l2_rel) & is.finite(kge))],
    "Bootstrap con al menos 1000 remuestreos" = N_BOOT >= 1000)
  res_red <- data.table(prueba = names(pruebas_red), pasa = unlist(pruebas_red))
  print(tabla_f(res_red, "F7_pruebas_red"))
  cat(sprintf("\nPruebas de la evaluación de la red: %d de %d pasan\n", sum(res_red$pasa), nrow(res_red)))
  stopifnot(res_red$pasa)
}

# QUÉ DICEN LOS DATOS (sección 7)
#   Control cruzado: la regresión de R rehecha en Python da las mismas curvas (1e-10, la
#   precisión del CSV), así que Python y R miden lo mismo.
#   Prueba 2001-2003 (370 curvas, 137 estaciones; media ± DE de 3 semillas):
#     L² relativo  satélite 0,513 | regresión R 0,437 | lineal 0,440 ± 0,002 | KAN 0,434 ± 0,002 | MLP 0,420 ± 0,010
#     KGE          satélite 0,382 | regresión R 0,467 | lineal 0,454 ± 0,005 | KAN 0,471 ± 0,009 | MLP 0,485 ± 0,020
#   Las tres cabezas mejoran al satélite sin corregir en L², KGE, |sesgo| y pico 1 (IC que
#   excluyen el 0). Entre cabezas (F7_comparaciones_pareadas_prueba.csv):
#     - MLP supera a la lineal en L² (−0,020 [−0,040, −0,002]) y KGE (+0,030 [0,003, 0,058]),
#       por poco; la lineal ubica mejor el pico 1 (0,7 días).
#     - KAN no se distingue de la lineal en L² ni KGE; baja el |sesgo| 4 puntos y ubica peor el
#       pico 2 (2 días).
#     - KAN frente a MLP: iguales en L², KGE, sesgo y pico 1; el MLP mejor en RMSE y pico 2.
#     - La red lineal no se distingue de la regresión de R en L², pero queda algo peor en KGE
#       (−0,013) y sesgo: la penalización de la red (rugosidad de w_k) no es la de R (de β).
#   Por fase ENSO (F7_metricas_por_fase_ENSO.csv) aparece lo más interesante:
#     - Años Niño (1997 en validación, 2002 en prueba): KAN y MLP son las mejores (L² 0,38
#       frente a 0,40 de la lineal en 2002; 0,45 frente a 0,47 en 1997).
#     - Años Niña (1999-2000, validación): la lineal es la mejor (0,39) y las no lineales
#       empeoran (KAN 0,45, MLP 0,42) y sobrestiman la lluvia (+14 % y +7 %).
#   La no linealidad del sesgo existe, pero depende de la fase: lo que KAN y MLP aprenden ayuda
#   en años secos (Niño) y estorba en años húmedos (Niña). Con solo 3 años de prueba y ninguno
#   Niña, la ventaja global de las no lineales es pequeña y no se puede generalizar.

# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ┃ RÚBRICA 4.5 — Interpretación climatológica de la red
# ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
titulo("8. Interpretación de la red (rúbrica 4.5)")
# Se grafica y se lee lo que aprendió la red. Python solo exporta números (pesos del codificador,
# de la cabeza lineal y activaciones de la primera capa de la KAN en el rango de los datos);
# todo el análisis se hace aquí. Se usa la semilla 2026 para las figuras y las tres semillas
# para ver qué tan estable es cada lectura.
ARCH_INT <- file.path("salidas_python", c("pesos_codificador.csv", "cabeza_lineal_pesos.csv", "kan_activaciones_capa1.csv", "kan_importancia_aristas.csv"))
if (!all(file.exists(ARCH_INT)) || !exists("mc_prom")) {
  cat("Faltan las salidas de lluvia_valle_red.py → la sección 8 se omite.\n")
} else {
  SEM <- 2026L
  meses <- data.table(dia = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335), mes = c("E", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"))
  cod <- fread(ARCH_INT[1]); lin <- fread(ARCH_INT[2]); act <- fread(ARCH_INT[3]); imp <- fread(ARCH_INT[4])

  # ---- 8.1 Las w_k(s): qué ventana del año mira cada neurona ----------------------------------
  # w_k(s) = Σ_j a_kj φ_j(s) con la base del satélite. La escala de w_k no importa (z se
  # estandariza), así que se normaliza a max|w_k| = 1. La "ventana" son los días con |w_k| ≥ 0,5.
  segmentos <- function(dias) {
    if (!length(dias)) return("")
    corte <- c(0, which(diff(dias) > 1), length(dias))
    paste(sapply(seq_len(length(corte) - 1), function(k) paste(fecha_dia(dias[corte[k] + 1]), "a", fecha_dia(dias[corte[k + 1]]))), collapse = "; ")
  }
  W_k <- cod[, {
    w <- drop(Phi_s %*% unlist(.SD)); w_n <- w / max(abs(w))
    .(dia = T_DIA, w = w, w_norm = w_n)
  }, by = .(cabeza, semilla, neurona), .SDcols = paste0("a", seq_len(K_s))]
  ventanas <- W_k[, .(dia_max_abs = which.max(abs(w_norm)), signo_en_max = sign(w_norm[which.max(abs(w_norm))]),
                     ventana = segmentos(dia[abs(w_norm) >= .5]), dias_en_ventana = sum(abs(w_norm) >= .5)), by = .(cabeza, semilla, neurona)]
  ventanas[, fecha_max := fecha_dia(dia_max_abs)]
  tabla_f(ventanas, "F8_ventanas_w_k")
  grafico_f(ggplot(W_k[semilla == SEM], aes(dia, w_norm, colour = cabeza)) + geom_hline(yintercept = 0, colour = "grey60") + geom_hline(yintercept = c(-.5, .5), linetype = 3, colour = "grey60") +
              geom_line(linewidth = .7) + facet_wrap(~ paste("neurona", neurona), ncol = 4) +
              scale_x_continuous(breaks = meses$dia, labels = meses$mes) +
              labs(x = "día del año del satélite (s)", y = "w_k(s) normalizada", colour = "cabeza",
                   title = "Pesos del codificador w_k(s): qué parte del año del satélite usa cada neurona (semilla 2026)",
                   subtitle = "Ventana = días con |w_k| ≥ 0,5 (líneas punteadas)"), "F8_pesos_w_k", ancho = 12, alto = 6)

  # ---- 8.2 Cabeza lineal: su β(s,t) frente a la regresión de R y el peso de cada entrada -----
  # La red lineal es la regresión funcional (sección 3 de Python): β(s,t) = Σ_k w_k(s) Σ_m (B_mk/σ_k) ψ_m(t).
  Bl <- dcast(lin[semilla == SEM & entrada != "sesgo"], componente ~ entrada, value.var = "peso")
  ent_nom <- c(paste0("z", 1:8), "altitud", "ONI")
  Bm <- as.matrix(Bl[order(componente), ..ent_nom])                     # M × 10
  cl <- cod[cabeza == "lineal" & semilla == SEM][order(neurona)]
  A_lin <- as.matrix(cl[, paste0("a", seq_len(K_s)), with = FALSE])
  BETA_red <- Phi_s %*% t(A_lin) %*% diag(1 / cl$de_z) %*% t(Bm[, 1:8]) %*% t(psi_t)
  cor_beta <- cor(as.vector(BETA_red), as.vector(BETA))
  # Estabilidad: la β(s,t) de cada semilla frente a la de R y entre semillas. Las neuronas sueltas no
  # son identificables (con cabeza lineal solo importa el espacio que generan las w_k: se pueden
  # rotar sin cambiar la predicción), pero β(s,t), que es lo que predice, sí debe ser estable.
  beta_de <- function(s) {
    B_s <- as.matrix(dcast(lin[semilla == s & entrada != "sesgo"], componente ~ entrada, value.var = "peso")[order(componente), ..ent_nom])
    c_s <- cod[cabeza == "lineal" & semilla == s][order(neurona)]
    Phi_s %*% t(as.matrix(c_s[, paste0("a", seq_len(K_s)), with = FALSE])) %*% diag(1 / c_s$de_z) %*% t(B_s[, 1:8]) %*% t(psi_t)
  }
  betas <- lapply(sort(unique(lin$semilla)), beta_de)
  estab_beta <- data.table(semilla = sort(unique(lin$semilla)), cor_con_regresion_R = sapply(betas, function(b) cor(as.vector(b), as.vector(BETA))),
                           cor_con_semilla_2026 = sapply(betas, function(b) cor(as.vector(b), as.vector(betas[[1]]))))
  print(tabla_f(estab_beta, "F8_estabilidad_beta_lineal"))
  cat(sprintf("β(s,t) de la red lineal (semilla %d) frente a la de la regresión de R: correlación %.3f\n", SEM, cor_beta))
  grafico_base("F8_beta_red_vs_regresion", function() {
    par(mfrow = c(1, 2)); lim <- max(abs(c(BETA, BETA_red))) * c(-1, 1)
    image(T_DIA, T_DIA, BETA, zlim = lim, col = hcl.colors(64, "Blue-Red 3"), xlab = "s: día del satélite", ylab = "t: día de la estación", main = "Regresión funcional (R)")
    image(T_DIA, T_DIA, BETA_red, zlim = lim, col = hcl.colors(64, "Blue-Red 3"), xlab = "s: día del satélite", ylab = "t: día de la estación",
          main = sprintf("Red con cabeza lineal (correlación %.2f)", cor_beta))
  }, ancho = 12, alto = 5)
  # Importancia de cada entrada en la cabeza lineal: norma L² del cambio en la curva √P cuando la
  # entrada sube 1 DE, ||Σ_m B_mi ψ_m(t)||; en la KAN: Σ_j DE de φ_ij(x) sobre entrenamiento.
  imp_lin <- rbindlist(lapply(unique(lin$semilla), function(s) {
    Bs <- as.matrix(dcast(lin[semilla == s & entrada != "sesgo"], componente ~ entrada, value.var = "peso")[order(componente), ..ent_nom])
    data.table(semilla = s, entrada = ent_nom, importancia = sqrt(colSums((psi_t %*% Bs)^2) / 365))
  }))
  imp_kan <- imp[, .(importancia = sum(de_phi_entrenamiento)), by = .(semilla, entrada)]
  importancia <- merge(imp_lin[, .(imp_lineal_media = mean(importancia), imp_lineal_DE = sd(importancia)), by = entrada],
                       imp_kan[, .(imp_kan_media = mean(importancia), imp_kan_DE = sd(importancia)), by = entrada], by = "entrada")
  importancia[, `:=`(rango_lineal = frank(-imp_lineal_media), rango_kan = frank(-imp_kan_media))]
  importancia <- merge(importancia, ventanas[semilla == SEM & cabeza == "lineal", .(entrada = paste0("z", neurona), ventana_lineal = ventana)], by = "entrada", all.x = TRUE)
  importancia <- merge(importancia, ventanas[semilla == SEM & cabeza == "kan", .(entrada = paste0("z", neurona), ventana_kan = ventana)], by = "entrada", all.x = TRUE)
  setorder(importancia, rango_kan)
  print(tabla_f(importancia, "F8_importancia_entradas"))

  # ---- 8.3 Activaciones de la primera capa de la KAN en el rango de los datos ------------------
  # No linealidad de cada arista = 1 − R² de una recta ajustada a φ_ij(x) en el rango p1–p99.
  aristas <- act[, .(no_linealidad = 1 - summary(lm(phi ~ x))$r.squared, rango_phi = diff(range(phi))), by = .(semilla, entrada, oculta)]
  aristas <- merge(aristas, imp, by = c("semilla", "entrada", "oculta"))
  tabla_f(aristas, "F8_kan_aristas")
  a_sem <- aristas[semilla == SEM]
  grafico_f(ggplot(merge(act[semilla == SEM], a_sem[, .(entrada, oculta, de_phi_entrenamiento)], by = c("entrada", "oculta")),
                   aes(x, phi, colour = de_phi_entrenamiento)) + geom_line(linewidth = .8) +
              facet_grid(factor(paste("oculta", oculta)) ~ factor(entrada, levels = ent_nom), scales = "free") +
              scale_colour_viridis_c(name = "importancia\n(DE de φ)") +
              labs(x = "entrada estandarizada (rango p1–p99 de entrenamiento)", y = "φ(x)",
                   title = "Activaciones de la primera capa de la KAN (semilla 2026): cada panel es una arista entrada → neurona oculta"),
            "F8_kan_activaciones", ancho = 15, alto = 7)

  # ---- 8.4 Ajuste simbólico de las activaciones más importantes --------------------------------
  # Candidatas de menor a mayor complejidad; se elige la más simple con R² ≥ 0,99.
  candidatas <- list(
    lineal = list(f = phi ~ x, ini = NULL),
    cuadratica = list(f = phi ~ x + I(x^2), ini = NULL),
    exponencial = list(f = phi ~ a + b * exp(c * x), ini = list(c(a = 0, b = 1, c = .5), c(a = 0, b = -1, c = -.5), c(a = 0, b = 1, c = -.5))),
    tanh = list(f = phi ~ a + b * tanh(c * x + d), ini = list(c(a = 0, b = 1, c = 1, d = 0), c(a = 0, b = -1, c = 1, d = 0), c(a = 0, b = 1, c = 3, d = 1))),
    cubica = list(f = phi ~ x + I(x^2) + I(x^3), ini = NULL),
    seno = list(f = phi ~ a + b * sin(c * x + d), ini = list(c(a = 0, b = 1, c = 1, d = 0), c(a = 0, b = 1, c = 2, d = 1))))
  ajusta <- function(d, nom) {
    cand <- candidatas[[nom]]
    m <- if (is.null(cand$ini)) lm(cand$f, d) else {
      fits <- lapply(cand$ini, function(s0) tryCatch(suppressWarnings(nls(cand$f, d, start = s0, control = nls.control(maxiter = 500, warnOnly = TRUE))), error = function(e) NULL))
      fits <- Filter(Negate(is.null), fits); if (!length(fits)) return(NULL); fits[[which.min(sapply(fits, deviance))]] }
    r2 <- 1 - sum(residuals(m)^2) / sum((d$phi - mean(d$phi))^2)
    data.table(forma = nom, R2 = r2, coeficientes = paste(sprintf("%s = %.3f", names(coef(m)), coef(m)), collapse = ", "))
  }
  top <- a_sem[order(-de_phi_entrenamiento)][1:4]
  simbolico <- rbindlist(lapply(seq_len(nrow(top)), function(k) {
    d <- act[semilla == SEM & entrada == top$entrada[k] & oculta == top$oculta[k]]
    r <- rbindlist(lapply(names(candidatas), function(nm) ajusta(d, nm)))
    r[, `:=`(entrada = top$entrada[k], oculta = top$oculta[k], importancia = top$de_phi_entrenamiento[k], elegida = FALSE)]
    r[which(R2 >= .99)[1], elegida := TRUE]
    if (!any(r$elegida)) r[which.max(R2), elegida := TRUE]
    r
  }))
  setcolorder(simbolico, c("entrada", "oculta", "importancia", "forma", "R2", "elegida"))
  print(tabla_f(simbolico[, R2 := round(R2, 4)], "F8_ajuste_simbolico"))

  # ---- 8.5 Peores casos en prueba ----------------------------------------------------------------
  # Curvas de prueba con mayor error L² relativo medio entre las tres cabezas (promedio de semillas).
  peor <- mc_prom[particion == "prueba" & modelo %in% c("lineal", "kan", "mlp"), .(l2_cabezas = mean(l2_rel)), by = .(codigo, anio, fase)][order(-l2_cabezas)][1:6]
  peor <- merge(peor, est[, .(codigo, nombre, vertiente, altitud)], by = "codigo")[order(-l2_cabezas)]
  peor[, fila := sapply(seq_len(.N), function(k) which(pares$codigo == codigo[k] & pares$anio == anio[k]))]
  peor[, `:=`(P_anual_obs = rowSums(Y_mm[fila, , drop = FALSE]), P_anual_sat = rowSums(X_mm[fila, , drop = FALSE]))]
  curva_cabeza <- function(cab, f) colMeans(do.call(rbind, lapply(split(pr[modelo == cab & fila == f], by = "semilla"), function(d) a_mm(as.matrix(d[, ..cols_c]) %*% t(Phi_e)))))
  peor[, P_anual_kan := sapply(fila, function(f) sum(curva_cabeza("kan", f)))]
  print(tabla_f(peor[, .(codigo, nombre, anio, fase, vertiente, altitud, l2_cabezas = round(l2_cabezas, 2),
                         P_anual_obs = round(P_anual_obs), P_anual_sat = round(P_anual_sat), P_anual_kan = round(P_anual_kan))], "F8_peores_casos"))
  dib <- rbindlist(lapply(seq_len(nrow(peor)), function(k) {
    f <- peor$fila[k]; et <- sprintf("%s %d (%s, %s, %.0f m)", peor$nombre[k], peor$anio[k], peor$fase[k], peor$vertiente[k], peor$altitud[k])
    rbind(data.table(caso = et, dia = T_DIA, serie = "estación", mm = Y_mm[f, ]), data.table(caso = et, dia = T_DIA, serie = "satélite", mm = X_mm[f, ]),
          rbindlist(lapply(c("lineal", "kan", "mlp"), function(cb) data.table(caso = et, dia = T_DIA, serie = cb, mm = curva_cabeza(cb, f)))))
  }))
  grafico_f(ggplot(dib, aes(dia, mm, colour = serie, linewidth = serie == "estación")) + geom_line() + facet_wrap(~ caso, ncol = 2, scales = "free_y") +
              scale_linewidth_manual(values = c(.5, 1.1), guide = "none") + scale_x_continuous(breaks = meses$dia, labels = meses$mes) +
              labs(x = NULL, y = "mm/día", colour = NULL, title = "Los seis peores casos de prueba (mayor error L² relativo de las tres cabezas)"),
            "F8_peores_casos", ancho = 12, alto = 9)
}

# QUÉ DICEN LOS DATOS (sección 8, rúbrica 4.5) — semilla 2026, con la estabilidad entre semillas
#   w_k(s) (F8_pesos_w_k.png, F8_ventanas_w_k.csv):
#     - Cabeza lineal: las w_k son casi rectas (λ_w = 1e8 las suaviza mucho) y cada una pesa el año
#       entero con una tendencia. Una neurona suelta no se puede interpretar: con cabeza lineal solo
#       importa el espacio que generan las 8 w_k (se pueden rotar sin cambiar la predicción), y por
#       eso la importancia de cada neurona cambia mucho entre semillas. Lo que sí es estable es β(s,t):
#       correlación 0,98-0,99 entre semillas y 0,81-0,86 con la regresión de R (F8_estabilidad_beta_lineal.csv).
#     - KAN y MLP (λ_w 1e6-1e7): las w_k contrastan temporadas del régimen bimodal: la 1.ª temporada
#       de lluvias (mar-may) frente a la 2.ª (oct-nov) (neuronas 3 y 6), la temporada seca de
#       mitad de año (jun-ago) frente a la de fin de año (neurona 8), comienzo de año frente a
#       julio-agosto (neurona 5). Es decir, miran la FORMA del año del satélite, no solo su total.
#   Cabeza lineal frente a regresión (F8_beta_red_vs_regresion.png): la β(s,t) de la red tiene la
#     misma estructura que la de R (correlación 0,81) pero más suave en s.
#   Activaciones de la KAN (F8_kan_activaciones.png, F8_kan_aristas.csv):
#     - Las cuatro aristas más importantes (z6, z1, z4 → oculta 3; z6 → oculta 1) son de SATURACIÓN:
#       tanh con R² 0,992-0,997 (F8_ajuste_simbolico.csv; z6→3 queda mejor como cuadrática). Cuando el
#       satélite marca mucha lluvia en esa ventana, la corrección deja de crecer: es la forma de un
#       sesgo que aumenta con la intensidad (el satélite sobrestima más en los meses más lluviosos).
#     - La entrada más no lineal es la ALTITUD (1 − R² de una recta: 0,72 en la mediana de sus
#       aristas; el resto entre 0,14 y 0,36). Su efecto no es monótono: cambia de sentido hacia la
#       media (≈ 1 400 m) y crece fuerte por encima de +2 DE (≈ 2 500 m). Es coherente con la
#       orografía andina, donde la lluvia no aumenta sin fin con la altura (hay un máximo a altitudes
#       medias en las vertientes). La cabeza lineal casi no usa la altitud (importancia 0,04 de 0,56
#       máx.): ese efecto solo aparece cuando se permite la no linealidad.
#     - El ONI entra casi monótono y suave (crece y satura): más ONI (Niño) → la corrección cambia
#       en un sentido estable; coherente con que la ventaja de KAN/MLP sea mayor en años Niño (sección 7).
#   Peores casos (F8_peores_casos.png/.csv): las seis curvas con más error son estaciones que
#   registran muy poca lluvia (500-830 mm al año) donde el satélite marca 1 000-2 400 mm: el fondo
#   seco del valle del Cauca (Guacarí, La Despensa, Aují, 990-1 555 m) y dos estaciones de la cuenca
#   del Dagua (Loboguerrero, Puerto Nuevo; 700-990 m), una zona más seca que la vertiente Pacífica
#   que la rodea. El píxel de 5 km promedia laderas más lluviosas y ninguna cabeza puede bajar tanto
#   desde la curva del satélite: el límite es de representatividad del píxel, no del modelo. Varias
#   de esas curvas caen a 0 al final de diciembre (borde de la base con días faltantes).

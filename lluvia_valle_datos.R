# =============================================================================
# LLUVIA DEL VALLE DEL CAUCA 1981–2003 · ENTENDIMIENTO, DEPURACIÓN Y DATASET FINAL
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
#   8. Depuración (reglas D-36) y dataset final
#   9. Resumen
#
# Ejecutar desde la raíz del repositorio:  Rscript lluvia_valle_datos.R
# (o abrirlo en RStudio con la raíz del repositorio como directorio de trabajo
#  y pulsar Source). Tarda unos minutos, sobre todo la sección 7.
#
# Entradas
#   P_1980_1990.csv, P_1991_2003.csv              lluvia diaria (STATION_ID, YEAR, MONTH, DAY, VALUE)
#   datos_web/estaciones_valle_finales.csv        las 232 estaciones del Valle (D-34/D-35)
#   datos_web/estaciones_valle_del_cauca.xlsx     inventario de los 344 códigos (D-31/D-32)
#   datos_web/limites/deptos_valle_y_vecinos_DANE_2020.gpkg   límites DANE (MGN 2020)
#   eda_fase0/oni_anual_fase.csv                  ONI anual (contexto ENSO)
# Salidas (salidas_R_datos/)
#   lluvia_valle_limpia.rds   dataset final (ver sección 8)
#   tablas/*.csv              tablas de cada sección
#   graficos/*.png            gráficos de cada sección
#
# Librerías: data.table, ggplot2, readxl, sf (EDA). La base B-spline del
# experimento de huecos se programa sin librerías (recursión de Cox–de Boor),
# como pide el curso para la construcción funcional.
# Decisiones y su justificación: DECISIONES.md (D-01 a D-36).
# =============================================================================


# ---- 0. CONFIGURACIÓN -------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(readxl)
  library(sf)
})
sf_use_s2(FALSE)
theme_set(theme_minimal(base_size = 11))
set.seed(20261112)

if (!file.exists("P_1991_2003.csv")) stop("Ejecute el script desde la raíz del repositorio")

ANIOS       <- 1981:2003                    # fe de erratas
PARTICION   <- function(anio) fcase(anio <= 1995, "entrenamiento", anio <= 2000, "validación", default = "prueba")
RADIO_VEC   <- 25                           # km: vecinas para cualquier contraste (H9)
COB_MIN     <- 0.90                         # R5: fracción mínima de días en el año
HUECO_MAX   <- 15                           # R5: hueco máximo permitido (días)
ANIOS_MIN   <- 7                            # R6: años válidos mínimos en entrenamiento

SALIDA <- "salidas_R_datos"
dir.create(file.path(SALIDA, "tablas"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(SALIDA, "graficos"), showWarnings = FALSE, recursive = TRUE)
tabla   <- function(x, nombre) { fwrite(x, file.path(SALIDA, "tablas", paste0(nombre, ".csv"))); invisible(x) }
grafico <- function(p, nombre, ancho = 9, alto = 5) {
  ggsave(file.path(SALIDA, "graficos", paste0(nombre, ".png")), p, width = ancho, height = alto, dpi = 110, bg = "white")
  if (interactive()) print(p)
  invisible(p)
}
titulo <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

dist_km <- function(lat1, lon1, lat2, lon2) {   # distancia de gran círculo
  r <- pi / 180
  a <- sin((lat2 - lat1) * r / 2)^2 + cos(lat1 * r) * cos(lat2 * r) * sin((lon2 - lon1) * r / 2)^2
  12742 * asin(sqrt(a))
}
hueco_mas_largo <- function(obs) { r <- rle(!obs); if (any(r$values)) max(r$lengths[r$values]) else 0L }


# ---- 1. POBLACIÓN: LAS 232 ESTACIONES DEL VALLE -----------------------------
titulo("1. Población")
# Los 344 códigos de los tres archivos de lluvia se ubicaron con fuentes
# oficiales (D-31/D-32). Quedaron 232 del Valle con coordenada: 226 con ubicación
# oficial y 6 identificadas por cruce con el CNE del IDEAM (D-34). La altitud
# solo sale de fuentes oficiales (D-35).

crudo <- rbind(fread("P_1980_1990.csv"), fread("P_1991_2003.csv"))
crudo[, codigo := as.character(STATION_ID)][, STATION_ID := NULL]
crudo <- crudo[YEAR %in% ANIOS]
cobertura_cruda <- crudo[, .(dias_con_dato = .N, anio_ini = min(YEAR), anio_fin = max(YEAR),
                             anios_con_dato = uniqueN(YEAR)), by = codigo]

inventario <- as.data.table(read_excel("datos_web/estaciones_valle_del_cauca.xlsx", sheet = 2))
inventario <- inventario[, .(codigo = format(inventario[[1]], scientific = FALSE, trim = TRUE),
                             depto_inventario = Departamento, nivel_inventario = `Nivel de evidencia`)]

finales <- fread("datos_web/estaciones_valle_finales.csv", encoding = "UTF-8", colClasses = list(character = "codigo"))
finales <- finales[, .(codigo, nombre, municipio, lat, lon, altitud = alt, nivel_c, fuente, fuente_alt, srtm, revisar)]

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
cat("Valle con datos en 1981–2003:", nrow(valle), "de 232 | sin altitud:", valle[is.na(altitud), .N],
    "| coordenada a revisar:", valle[revisar != "", .N], "| a ≤2 km de otro departamento:", valle[km_borde <= 2, .N], "\n")
print(valle[, .N, by = .(nivel_c, vertiente)])

grafico(ggplot(est[grupo == "Valle del Cauca"], aes(lon, lat, colour = dias_con_dato > 0)) +
          geom_sf(data = deptos[deptos$shapeName == "Valle del Cauca", ], inherit.aes = FALSE, fill = NA, colour = "grey50") +
          geom_point(size = 1.8, alpha = .8) +
          scale_colour_manual(values = c(`TRUE` = "steelblue", `FALSE` = "orange"),
                              labels = c(`TRUE` = "con datos 1981–2003", `FALSE` = "sin datos en el periodo"), name = NULL) +
          labs(title = "Las 232 estaciones del Valle y el polígono DANE", x = NULL, y = NULL),
        "1_mapa_estaciones", ancho = 7, alto = 7)

# Matriz de distancias entre las estaciones del Valle con datos
DK <- outer(seq_len(nrow(valle)), seq_len(nrow(valle)),
            function(i, j) dist_km(valle$lat[i], valle$lon[i], valle$lat[j], valle$lon[j]))
dimnames(DK) <- list(valle$codigo, valle$codigo)
vecinas <- function(cod) setdiff(names(which(DK[cod, ] <= RADIO_VEC)), cod)

# QUÉ DICEN LOS DATOS (H1)
#   323 códigos con datos en 1981–2003. De las 232 del Valle, 220 tienen datos en
#   el periodo; las 12 restantes solo tienen datos fuera. Las 232 caen dentro del
#   polígono DANE del Valle (13 a ≤2 km de otro departamento).


# ---- 2. ESTRUCTURA DE LOS ARCHIVOS CRUDOS ------------------------------------
titulo("2. Estructura del archivo")
d <- crudo[codigo %in% valle$codigo, .(codigo, anio = YEAR, mes = MONTH, dia_mes = DAY, P = as.numeric(VALUE))]
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

# Los decimales vienen de las estaciones climatológicas (tipo 02/03 en los
# dígitos 7-8 del código), que leen a 0,1 mm. No son error.
decimales <- d[P > 0, .(pct_decimales = round(100 * mean(P != round(P)), 1), dias_lluvia = .N),
               by = .(tipo = substr(codigo, 7, 8))]
print(tabla(decimales, "2_decimales_por_tipo"))

grafico(ggplot(d[P > 0 & P == round(P)], aes(factor(P %% 10))) + geom_bar() +
          labs(x = "Último dígito del valor (mm, lluvia > 0)", y = "Días",
               title = "Preferencia de dígitos (si la lectura fuera exacta, los 10 dígitos serían parecidos)"),
        "2_preferencia_digitos")

# QUÉ DICEN LOS DATOS (H2)
#   Un faltante es una FILA AUSENTE: no hay NA, negativos ni duplicados, y no
#   aparece el valor imposible (859 375 mm) de otra estación. Hay 832 registros
#   del 29-feb y 5 144 valores con decimales (estaciones climatológicas).


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

orden <- cob[, .(m = mean(pct)), by = codigo][order(m), codigo]
grafico(ggplot(cob, aes(anio, factor(codigo, levels = orden), fill = pct)) + geom_tile() +
          scale_fill_viridis_c(name = "% días\ncon dato") +
          geom_vline(xintercept = c(1995.5, 2000.5), colour = "red", linetype = 2) +
          labs(x = NULL, y = "Estación (ordenada por cobertura)", title = "Cobertura por estación y año",
               subtitle = "Líneas rojas: entrenamiento 1981–1995 | validación 1996–2000 | prueba 2001–2003") +
          theme(axis.text.y = element_blank(), panel.grid = element_blank()),
        "3_cobertura_heatmap", alto = 8)
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


# ---- 4. VALORES: DISTRIBUCIÓN, EXTREMOS Y VECINAS ----------------------------
titulo("4. Valores")
print(d[, .(dias = .N, pct_ceros = round(100 * mean(P == 0), 1), p50_lluvia = median(P[P > 0]),
            p90 = quantile(P[P > 0], .9), p99 = quantile(P[P > 0], .99), p999 = quantile(P[P > 0], .999), max = max(P))])

p99 <- d[P > 0, .(p99 = quantile(P, .99), max = max(P), pct_ceros = 100 * mean(P == 0)), by = codigo]
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

grafico(ggplot(valle[, .(codigo, vertiente)][d[P > 0], on = "codigo"], aes(P, colour = vertiente)) +
          stat_ecdf() + scale_x_log10() +
          labs(x = "Lluvia diaria (mm, días con lluvia, log)", y = "F(x)", colour = "Vertiente",
               title = "Distribución de la lluvia diaria por vertiente"),
        "4_distribucion_lluvia")

# QUÉ DICEN LOS DATOS (H4)
#   60,6 % de los días son ceros; p99 ≈ 76 mm; máximo 600 mm. Solo 5 días superan
#   3·p99 y 150 mm, y en todos las vecinas registran como máximo 78 mm (Pardo 600,
#   Heraclio Uribe 503, La Cascada 421, Costa Rica 283, Venecia 283) → regla R2.


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
print(rz[, .(rachas = .N, de_30_o_mas = sum(largo >= 30), sin_vecinas = sum(is.na(p_lluvia_vecinas)),
             de_30_y_vecinas_llueven_40 = sum(largo >= 30 & p_lluvia_vecinas >= .4, na.rm = TRUE))])
rachas_revisar <- rz[largo >= 30 & p_lluvia_vecinas >= .4]
print(rachas_revisar)
grafico(ggplot(rz[!is.na(p_lluvia_vecinas)], aes(largo, p_lluvia_vecinas, colour = vertiente)) +
          geom_point(alpha = .5) + scale_x_log10() + geom_hline(yintercept = .4, linetype = 2) +
          labs(x = "Largo de la racha de ceros (días, log)", y = "Fracción de días-vecina con lluvia (≤25 km)",
               title = "Rachas ≥15 días de ceros: ¿las vecinas también estaban secas?"),
        "5b_rachas_ceros")

# Revisión a mano de las rachas marcadas: en los mismos meses de otros años,
# ¿cuánto llueve la estación? ¿Y la vecina más cercana durante la racha?
revision <- rbindlist(lapply(seq_len(nrow(rachas_revisar)), function(i) {
  cod <- rachas_revisar$codigo[i]; a <- rachas_revisar$ini[i]; b <- rachas_revisar$fin[i]
  meses <- unique(month(seq(a, b, by = 1)))
  cercana <- names(sort(DK[cod, vecinas(cod)]))[1]
  data.table(codigo = cod, nombre = rachas_revisar$nombre[i], ini = a, fin = b, largo = rachas_revisar$largo[i],
             pct_dias_lluvia_normal = round(100 * d[codigo == cod & month(fecha) %in% meses & !(fecha %between% c(a, b)), mean(P > 0)]),
             vecina_cercana = valle[codigo == cercana, nombre], km = round(DK[cod, cercana], 1),
             dias_lluvia_vecina = d[codigo == cercana & fecha %between% c(a, b), sum(P > 0)],
             mm_vecina = d[codigo == cercana & fecha %between% c(a, b), sum(P)])
}))
print(tabla(revision, "5b_revision_rachas"))

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
#   (se anotó el total mensual) → R3. De 918 rachas de ceros ≥30 días, solo 6
#   tienen vecinas lloviendo ≥40 % (5 de Loboguerrero, 1 de La Despensa) → R4.
#   No hay años copiados; lo más repetido son 11 días seguidos de 1 mm (llovizna).


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
grafico(ggplot(clim[!is.na(altitud)], aes(altitud, P_media, colour = vertiente)) + geom_point() +
          scale_y_log10() + geom_smooth(se = FALSE, method = "loess", formula = y ~ x) +
          labs(x = "Altitud (m)", y = "mm/año (log)", title = "Lluvia anual vs. altitud"),
        "6_lluvia_altitud")

men <- d[buenos, on = .(codigo, anio), nomatch = 0][, .(P = mean(P)), by = .(codigo, mes)]
men <- valle[, .(codigo, vertiente)][men, on = "codigo"]
grafico(ggplot(men, aes(factor(mes), P)) + geom_boxplot(outlier.size = .5) + facet_wrap(~vertiente, scales = "free_y") +
          labs(x = "Mes", y = "mm/día", title = "Ciclo anual: bimodal en la vertiente del Cauca"),
        "6_ciclo_anual")

oni <- fread("eda_fase0/oni_anual_fase.csv")
reg <- anual[valle[, .(codigo, vertiente)], on = "codigo", nomatch = 0]
reg[, anom := P_anual / clim$P_media[match(codigo, clim$codigo)] - 1]
reg <- reg[, .(anom = median(anom)), by = .(vertiente, anio)]
reg <- oni[reg, on = "anio"]
print(reg[, .(r_ONI_anomalia = round(cor(ONI_anual, anom), 2)), by = vertiente])
grafico(ggplot(reg, aes(anio, 100 * anom, fill = fase)) + geom_col() + facet_wrap(~vertiente, ncol = 1) +
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

# R4. Rachas de ceros: no se borran por regla. Revisión a mano (tabla 5b):
#     - Loboguerrero, oct–nov 1995 → faltante: temporada de lluvias, la estación
#       suele llover el 41 % de esos días y Dagua (12 km, mismo cañón) tuvo 12
#       días con lluvia (83 mm).
#     - Las otras 5 (La Despensa ene–feb 1996; Loboguerrero 1982, 1982–83, 1989,
#       1991) caen en temporada seca o en El Niño y la vecina del mismo cañón
#       también estuvo casi seca → se mantienen.
m_cero <- d[codigo == "5301500102" & fecha %between% as.IDate(c("1995-10-27", "1995-11-28")),
            .(codigo, fecha, P, regla = "R4 racha de ceros revisada",
              detalle = "33 días en 0 en temporada de lluvias; Dagua (12 km) 12 días con lluvia")]
stopifnot(nrow(m_cero) == 33, all(m_cero$P == 0))

marcas <- valle[, .(codigo, nombre)][rbind(m_ext, m_mes, m_cero), on = "codigo"]
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

reglas <- c(
  R1 = "Dominio 1–365: se quita el 29 de febrero",
  R2 = "Extremo > 3·p99 y ≥150 mm sin vecina (≤25 km, ±1 día) que supere P/3 → faltante",
  R3 = "Mes con un solo día de lluvia > p99 y (mes con <50 % de días o vecinas con lluvia ≥40 % de los días) → mes faltante",
  R4 = "Rachas de ceros: no se borran por regla; revisión a mano: Loboguerrero oct–nov 1995 → faltante",
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


# ---- 9. RESUMEN ---------------------------------------------------------------
titulo("9. Resumen")
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
writeLines(capture.output(sessionInfo()), file.path(SALIDA, "sessionInfo.txt"))

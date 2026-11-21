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

grafico(ggplot(d[P > 0 & P == round(P)], aes(factor(P %% 10))) + geom_bar() +
          labs(x = "Último dígito del valor (mm, lluvia > 0)", y = "Días",
               title = "Preferencia de dígitos (si la lectura fuera exacta, los 10 dígitos serían parecidos)"),
        "2_preferencia_digitos")

# QUÉ DICEN LOS DATOS (H2)
#   Un faltante es una FILA AUSENTE: no hay NA, negativos ni duplicados, y no
#   aparece el valor imposible (859 375 mm) de otra estación. R0 lo comprueba con
#   un umbral citado (récord OMM, 1 825 mm en 24 h): en el Valle no quita nada. Hay 832 registros
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
print(rz[, .(rachas = .N, de_30_o_mas = sum(largo >= 30), sin_vecinas = sum(is.na(p_lluvia_vecinas)))])
grafico(ggplot(rz[!is.na(p_lluvia_vecinas)], aes(largo, p_lluvia_vecinas, colour = vertiente)) +
          geom_point(alpha = .5) + scale_x_log10() +
          labs(x = "Largo de la racha de ceros (días, log)", y = "Fracción de días-vecina con lluvia (≤25 km)",
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

oni <- fread("datos_web/oni_anual_fase.csv")
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
               title = paste0("Serie péntada a péntada de ", tope_est$nombre[1], " (", tope_est$vertiente[1], "), 1981–2003"),
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
grafico(ggplot(por_est_sat, aes(mm_anual_est, mm_anual_sat, colour = vertiente)) + geom_point() + geom_abline() +
          scale_x_log10() + scale_y_log10() +
          labs(x = "Estación (mm/año, log)", y = "Satélite en su píxel (mm/año, log)",
               title = "El satélite aplana: sobreestima donde llueve poco y subestima donde llueve mucho"),
        "12_sesgo_por_estacion")
grafico(ggplot(por_est_sat, aes(r, fill = vertiente)) + geom_histogram(bins = 30) +
          labs(x = "r entre satélite y estación (péntadas)", y = "Estaciones", title = "Correlación por estación"),
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

# =============================================================================
# 05_limpieza.R
# Aplica las reglas de D-36 y deja el conjunto limpio de estaciones del Valle.
# Periodo: 1981–2003 (fe de erratas). Partición: entrenamiento 1981–1995,
# validación 1996–2000, prueba 2001–2003.
#
# Reglas (D-36):
#   R1. Dominio 1–365: se quita el 29 de febrero (D-14).
#   R2. Extremo sin respaldo → faltante: P > 3·p99 de la estación, P ≥ 150 mm y
#       ninguna vecina a ≤25 km supera P/3 ese día, el anterior o el siguiente.
#   R3. Total mensual anotado en un día → mes entero faltante: el mes tiene un
#       solo día con lluvia y ese valor supera el p99, y además (a) el mes tiene
#       menos de la mitad de sus días observados, o (b) las vecinas a ≤25 km
#       llovieron en ≥40 % de sus días de ese mes (los ceros de la estación serían
#       falsos).
#   R4. Rachas de ceros: no se borran por regla. Las 6 rachas de ≥30 días con
#       vecinas lloviendo (≥40 %) se revisaron a mano; solo una pasa a faltante
#       (Loboguerrero, oct–nov 1995: temporada de lluvias, la estación suele
#       llover el 41 % de esos días y Dagua, a 12 km en el mismo cañón, tuvo 12
#       días con lluvia). Las otras 5 caen en temporada seca o en El Niño y la
#       vecina del mismo cañón también estuvo casi seca: se mantienen.
#   R5. Año válido: ≥90 % de los 365 días con dato y ningún hueco > 15 días.
#   R6. Estación incluida: ≥7 años válidos en entrenamiento.
#   R7. No se imputa (D-13).
#
# Salida: analisis_R/salidas/lluvia_valle_limpia.rds, una lista con
#   $diario     codigo, fecha, anio, dia (1–365), particion, P (NA = faltante),
#               solo estaciones incluidas y años válidos
#   $diario_todo  igual, pero con las 220 estaciones y todos sus años (para el EDA)
#   $anios      codigo, anio, particion, n_obs, hueco_max, valido
#   $estaciones metadatos de las estaciones incluidas
#   $marcas     días que se pasaron a faltante y por qué
#   $reglas     texto de las reglas
# Ejecutar desde la raíz:  Rscript analisis_R/05_limpieza.R
# =============================================================================
suppressPackageStartupMessages(library(data.table))

est <- fread("analisis_R/salidas/estaciones_metadatos.csv", colClasses = list(character = "codigo"))
valle <- est[departamento == "Valle del Cauca" & dias_con_dato > 0]

crudo <- rbind(fread("P_1980_1990.csv"), fread("P_1991_2003.csv"))
crudo[, codigo := as.character(STATION_ID)]
d <- crudo[YEAR %between% c(1981, 2003) & codigo %in% valle$codigo,
           .(codigo, anio = YEAR, mes = MONTH, dia_mes = DAY, P = as.numeric(VALUE))]
rm(crudo)

# --- R1: dominio 1–365 ---------------------------------------------------------
d <- d[!(mes == 2 & dia_mes == 29)]
d[, fecha := as.IDate(sprintf("%d-%02d-%02d", anio, mes, dia_mes))]
d[, dia := yday(fecha) - (anio %% 4 == 0 & mes > 2)]
stopifnot(d[, all(dia %between% c(1, 365))], !anyDuplicated(d[, .(codigo, fecha)]))

# Vecinas a ≤25 km
dist_km <- function(lat1, lon1, lat2, lon2) {
  r <- pi / 180; a <- sin((lat2 - lat1) * r / 2)^2 + cos(lat1 * r) * cos(lat2 * r) * sin((lon2 - lon1) * r / 2)^2
  12742 * asin(sqrt(a))
}
DK <- outer(seq_len(nrow(valle)), seq_len(nrow(valle)),
            function(i, j) dist_km(valle$lat[i], valle$lon[i], valle$lat[j], valle$lon[j]))
dimnames(DK) <- list(valle$codigo, valle$codigo)
vecinas <- function(cod) setdiff(names(which(DK[cod, ] <= 25)), cod)

p99 <- d[P > 0, .(p99 = quantile(P, .99)), by = codigo]
d <- p99[d, on = "codigo"]
setkey(d, codigo, fecha)

# --- R2: extremos sin respaldo -------------------------------------------------
cand <- d[P > 3 * p99 & P >= 150]
# (se sacan los vectores antes: dentro de d[...], `codigo` sería la columna de d)
cand$vecina_max <- mapply(function(cod, f) {
  v <- d[codigo %in% vecinas(cod) & fecha %between% c(f - 1L, f + 1L), P]
  if (length(v)) max(v) else NA_real_
}, cand$codigo, cand$fecha)
m_ext <- cand[!is.na(vecina_max) & vecina_max < P / 3,
              .(codigo, fecha, P, regla = "R2 extremo sin respaldo",
                detalle = sprintf("p99 %.0f; vecina máx %.0f mm", p99, vecina_max))]

# --- R3: totales mensuales en un día --------------------------------------------
mesd <- d[, .(dias_obs = .N, dias_lluvia = sum(P > 0), max = max(P), p99 = p99[1]), by = .(codigo, anio, mes)]
mesd[, dias_mes := fifelse(mes == 2, 28L, c(31L, 28L, 31L, 30L, 31L, 30L, 31L, 31L, 30L, 31L, 30L, 31L)[mes])]
sosp <- mesd[dias_lluvia == 1 & max > p99]
sosp$p_lluvia_vecinas <- mapply(function(cod, a, m) {
  v <- d[codigo %in% vecinas(cod) & anio == a & mes == m, P]
  if (length(v)) mean(v > 0) else NA_real_
}, sosp$codigo, sosp$anio, sosp$mes)
sosp[, mes_vacio := dias_obs < dias_mes / 2]
sosp[, marcar := mes_vacio | (!is.na(p_lluvia_vecinas) & p_lluvia_vecinas >= 0.4)]
m_mes <- d[sosp[marcar == TRUE, .(codigo, anio, mes, dias_obs, p_lluvia_vecinas, mes_vacio)],
           on = .(codigo, anio, mes), nomatch = 0][
  , .(codigo, fecha, P, regla = "R3 total mensual en un día",
      detalle = sprintf("%d días observados; vecinas llovieron %.0f %% de los días%s",
                        dias_obs, 100 * p_lluvia_vecinas, fifelse(mes_vacio, "; mes casi vacío", "")))]

# --- R4: rachas de ceros revisadas a mano ---------------------------------------
m_cero <- d[codigo == "5301500102" & fecha %between% as.IDate(c("1995-10-27", "1995-11-28")),
            .(codigo, fecha, P, regla = "R4 racha de ceros revisada",
              detalle = "33 días en 0 en temporada de lluvias; Dagua (12 km) 12 días con lluvia")]
stopifnot(nrow(m_cero) == 33, all(m_cero$P == 0))

marcas <- rbind(m_ext, m_mes, m_cero)
d[marcas, on = .(codigo, fecha), P := NA_real_]
d <- d[!is.na(P)]   # un faltante es una fila ausente

# --- Rejilla completa estación × año × día ------------------------------------
anios_op <- d[, .(anio = unique(anio)), by = codigo]
g <- anios_op[, .(dia = 1:365), by = .(codigo, anio)]
g <- d[, .(codigo, anio, dia, fecha, P)][g, on = .(codigo, anio, dia)]
g[is.na(fecha), fecha := as.IDate(sprintf("%d-01-01", anio)) + dia - 1L +
                         (anio %% 4 == 0 & dia > 59)]
g[, particion := fcase(anio <= 1995, "entrenamiento", anio <= 2000, "validación", default = "prueba")]
setkey(g, codigo, anio, dia)

# --- R5: años válidos ----------------------------------------------------------
hueco_max <- function(obs) { r <- rle(!obs); if (any(r$values)) max(r$lengths[r$values]) else 0L }
anios <- g[, .(n_obs = sum(!is.na(P)), hueco_max = hueco_max(!is.na(P))), by = .(codigo, anio, particion)]
anios[, valido := n_obs >= 0.9 * 365 & hueco_max <= 15]

# --- R6: estaciones incluidas --------------------------------------------------
porest <- anios[, .(validos_ent = sum(valido & particion == "entrenamiento"),
                    validos_val = sum(valido & particion == "validación"),
                    validos_pru = sum(valido & particion == "prueba")), by = codigo]
porest[, incluida := validos_ent >= 7]
estaciones <- merge(valle, porest, by = "codigo")[incluida == TRUE]
estaciones[, vertiente := fifelse(substr(codigo, 1, 2) == "26", "Cauca (interandina)", "Pacífico")]

diario <- g[anios[valido == TRUE, .(codigo, anio)], on = .(codigo, anio), nomatch = 0][
  codigo %in% estaciones$codigo]

reglas <- c(
  R1 = "Dominio 1–365: se quita el 29 de febrero",
  R2 = "Extremo > 3·p99 y ≥150 mm sin vecina (≤25 km, ±1 día) que supere P/3 → faltante",
  R3 = "Mes con un solo día de lluvia > p99 y (mes con < 50 % de días o vecinas con lluvia ≥40 % de los días) → mes faltante",
  R4 = "Rachas de ceros: no se borran por regla; revisión a mano: Loboguerrero oct–nov 1995 → faltante",
  R5 = "Año válido: ≥90 % de 365 días y hueco máximo ≤15 días",
  R6 = "Estación: ≥7 años válidos en entrenamiento (1981–1995)",
  R7 = "Sin imputación")

saveRDS(list(diario = diario, diario_todo = g, anios = anios, estaciones = estaciones,
             marcas = marcas, meses_sospechosos = sosp, reglas = reglas),
        "analisis_R/salidas/lluvia_valle_limpia.rds", compress = "xz")
fwrite(marcas, "analisis_R/salidas/marcas_limpieza.csv")
fwrite(anios, "analisis_R/salidas/anios_validez.csv")
fwrite(estaciones, "analisis_R/salidas/estaciones_incluidas.csv")

# --- Resumen -------------------------------------------------------------------
cat("Marcas por regla:\n"); print(marcas[, .(dias = .N), by = regla])
cat("R2 extremos a faltante:", nrow(m_ext), "\n"); print(m_ext)
cat("\nR3 meses sospechosos:", nrow(sosp), "| pasados a faltante:", sosp[marcar == TRUE, .N],
    "(", nrow(m_mes), "días )\n")
print(merge(sosp, valle[, .(codigo, nombre)], by = "codigo")[
  , .(codigo, nombre, anio, mes, dias_obs, max, p99 = round(p99), p_vec = round(p_lluvia_vecinas, 2), marcar)])
cat("\nAños válidos:", anios[valido == TRUE, .N], "de", nrow(anios), "\n")
print(anios[, .(validos = sum(valido), total = .N), by = particion])
cat("\nEstaciones con datos:", nrow(valle), "| incluidas (≥7 válidos en entrenamiento):", nrow(estaciones), "\n")
print(estaciones[, .N, by = .(nivel_c)])
cat("Con coordenada a revisar:", estaciones[revisar != "", .N], "\n")
cat("Con ≥1 año válido en validación:", estaciones[validos_val >= 1, .N],
    "| en prueba:", estaciones[validos_pru >= 1, .N], "| en ambos:", estaciones[validos_val >= 1 & validos_pru >= 1, .N], "\n")
cat("Estación-años en $diario:", uniqueN(diario[, .(codigo, anio)]), "| filas:", nrow(diario),
    "| días faltantes dentro de años válidos:", diario[is.na(P), .N], "\n")
print(diario[, .(estacion_anios = uniqueN(paste(codigo, anio))), by = particion])

from pathlib import Path
import json
import textwrap


RUTA = Path(__file__).with_name("01_preprocesamiento_iterativo_desde_crudo.ipynb")


def md(text):
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": textwrap.dedent(text).strip() + "\n",
    }


def code(text):
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": textwrap.dedent(text).strip() + "\n",
    }


with RUTA.open("r", encoding="utf-8") as f:
    nb = json.load(f)

cells = nb["cells"]

# Registro de decisiones actualizado.
cells[21] = code(r'''
registro_decisiones = pd.DataFrame([
    [1, 'Periodo y partición temporal', 'APROBADA', '1981-1995 entrenamiento; 1996-2000 validación; 2001-2003 prueba.'],
    [2, 'Fechas y claves incompletas', 'APROBADA', 'Las fechas imposibles no entran en las curvas; los faltantes se conservan como faltantes.'],
    [3, 'Duplicados diarios', 'APROBADA', 'Los duplicados exactos se consolidan; los conflictivos se invalidan sin promediar.'],
    [4, 'Error confirmado 859375 mm', 'APROBADA', 'Se conserva el original y se convierte solo la copia limpia en NaN.'],
    [5, 'Totales mensuales y digitación', 'APROBADA - ADAPTADA DE JACOBO', 'Reglas relativas; umbrales estimados solo con entrenamiento y congelados.'],
    [6, 'Años copiados entre estaciones', 'APROBADA - ADAPTADA DE JACOBO', 'Coincidencia diaria exacta en 365 días; se excluyen ambos pares sin asignar culpable.'],
    [7, 'Ceros falsos', 'APROBADA - ADAPTADA DE JACOBO', 'Rachas >=30 días y apoyo regional >=40%; se convierten en faltantes, nunca en lluvia satelital.'],
    [8, 'Cobertura y continuidad', 'APROBADA', '90% anual, hueco interior <=30 días y <=15 días faltantes en cada borde.'],
    [9, 'Tratamiento de faltantes', 'APROBADA PARA PREPROCESAMIENTO', 'No se imputa en esta fase; la construcción funcional ajustará con días observados.'],
    [10, 'Mínimo de años válidos', 'APROBADA', 'Cinco años válidos de entrenamiento por estación; luego se aplica la lista congelada.'],
    [11, 'Estaciones reservadas por el docente', 'PENDIENTE', 'Incorporar la lista oficial sin usarla para definir reglas ni ajustar el modelo.'],
], columns=['orden', 'decision', 'estado', 'objetivo'])
display(registro_decisiones)
''')

totales_md = md(r'''
# Decisión 3: adaptar los controles de totales mensuales y digitación de Jacobo

Se reproducen las ideas de la rama de Jacobo, pero se corrige la fuga de información:

- las distribuciones y percentiles se estiman únicamente con **entrenamiento (1981–1995)**;
- los umbrales elegidos se congelan antes de mirar validación y prueba;
- cada cambio afecta solo `precipitacion_limpia_mm`; el valor original y una bandera quedan disponibles;
- las tablas de sensibilidad muestran el costo de cada alternativa en entrenamiento.

Se aplican cuatro controles en orden: mes con un único registro, posible total mensual rodeado de ceros, mes con repetición extrema y error aislado de digitación. No se usa un límite absoluto general que elimine toda lluvia alta.
''')

totales_code = code(r'''
CLAVE_MES = ['station_id', 'year', 'month', 'particion']
MASCARA_TRAIN = datos_trabajo.particion.eq('ENTRENAMIENTO')

def claves_en_tabla(datos, columnas, tabla_claves):
    if tabla_claves.empty:
        return pd.Series(False, index=datos.index)
    izquierda = pd.MultiIndex.from_frame(datos[columnas])
    derecha = pd.MultiIndex.from_frame(tabla_claves[columnas].drop_duplicates())
    return pd.Series(izquierda.isin(derecha), index=datos.index)

def resumen_estado(nombre):
    valido = datos_trabajo.precipitacion_limpia_mm.notna()
    filas = []
    for particion in ORDEN_PARTICIONES + ['TOTAL']:
        mascara = valido if particion == 'TOTAL' else valido & datos_trabajo.particion.eq(particion)
        pares = datos_trabajo.loc[mascara, ['station_id', 'year']].drop_duplicates()
        filas.append({
            'etapa': nombre,
            'particion': particion,
            'registros_validos': int(mascara.sum()),
            'estaciones': int(datos_trabajo.loc[mascara, 'station_id'].nunique()),
            'pares_estacion_anio_con_dato': int(len(pares)),
        })
    return filas

traza_calidad = []
traza_calidad.extend(resumen_estado('00 después del error confirmado 859375'))

# 3A. Mes con un solo registro observado: no describe una curva mensual y se invalida.
resumen_mes = (
    datos_trabajo.groupby(CLAVE_MES, observed=True, dropna=False)
    .agg(registros_mes=('precipitacion_limpia_mm', 'count'))
    .reset_index()
)
meses_unico_registro = resumen_mes.loc[resumen_mes.registros_mes.eq(1), CLAVE_MES].copy()
datos_trabajo['mes_unico_registro'] = claves_en_tabla(datos_trabajo, CLAVE_MES, meses_unico_registro)
datos_trabajo.loc[datos_trabajo.mes_unico_registro, 'precipitacion_limpia_mm'] = np.nan
traza_calidad.extend(resumen_estado('01 mes con un único registro -> NaN'))

# 3B. Posible total mensual cargado en un día y ceros el resto del mes.
# La referencia climática es la mediana del total mensual para la misma estación y mes calendario,
# estimada exclusivamente en entrenamiento después del paso anterior.
def diagnostico_mensual_actual():
    d = datos_trabajo.assign(
        positivo=datos_trabajo.precipitacion_limpia_mm.gt(0),
    )
    return (
        d.groupby(CLAVE_MES, observed=True, dropna=False)
        .agg(
            registros_validos=('precipitacion_limpia_mm', 'count'),
            dias_positivos=('positivo', 'sum'),
            total_mes=('precipitacion_limpia_mm', 'sum'),
            maximo_dia=('precipitacion_limpia_mm', 'max'),
        )
        .reset_index()
    )

resumen_mes = diagnostico_mensual_actual()
medianas_train = (
    resumen_mes.loc[resumen_mes.particion.eq('ENTRENAMIENTO') & resumen_mes.registros_validos.ge(20)]
    .groupby(['station_id', 'month'], observed=True)['total_mes']
    .median()
    .rename('mediana_total_train')
    .reset_index()
)
resumen_mes = resumen_mes.merge(medianas_train, on=['station_id', 'month'], how='left')

sensibilidad_total_mensual = []
for min_registros in [20, 25, 28]:
    for fraccion in [0.25, 0.50, 0.75]:
        umbral = np.maximum(50.0, fraccion * resumen_mes.mediana_total_train)
        seleccion = (
            resumen_mes.particion.eq('ENTRENAMIENTO')
            & resumen_mes.registros_validos.ge(min_registros)
            & resumen_mes.dias_positivos.eq(1)
            & resumen_mes.maximo_dia.ge(umbral)
        )
        sensibilidad_total_mensual.append({
            'min_registros_mes': min_registros,
            'fraccion_mediana_historica': fraccion,
            'meses_marcados_train': int(seleccion.sum()),
            'estaciones_train': int(resumen_mes.loc[seleccion, 'station_id'].nunique()),
            'pares_train': int(resumen_mes.loc[seleccion, ['station_id', 'year']].drop_duplicates().shape[0]),
            'seleccionado': min_registros == 25 and np.isclose(fraccion, .50),
        })
sensibilidad_total_mensual = pd.DataFrame(sensibilidad_total_mensual)

MIN_REGISTROS_TOTAL = 25
FRACCION_MEDIANA_TOTAL = .50
umbral_total = np.maximum(50.0, FRACCION_MEDIANA_TOTAL * resumen_mes.mediana_total_train)
mascara_total = (
    resumen_mes.registros_validos.ge(MIN_REGISTROS_TOTAL)
    & resumen_mes.dias_positivos.eq(1)
    & resumen_mes.maximo_dia.ge(umbral_total)
)
meses_total_mensual = resumen_mes.loc[mascara_total, CLAVE_MES + [
    'registros_validos', 'maximo_dia', 'total_mes', 'mediana_total_train'
]].copy()
datos_trabajo['posible_total_mensual'] = claves_en_tabla(datos_trabajo, CLAVE_MES, meses_total_mensual)
datos_trabajo.loc[datos_trabajo.posible_total_mensual, 'precipitacion_limpia_mm'] = np.nan
traza_calidad.extend(resumen_estado('02 posible total mensual -> mes a NaN'))

# 3C. Mes con repetición anómala de valores extremos relativos a la estación.
# El p99 se aprende con positivos de entrenamiento ya depurados en 3A y 3B.
p99_train = (
    datos_trabajo.loc[
        MASCARA_TRAIN & datos_trabajo.precipitacion_limpia_mm.gt(0),
        ['station_id', 'precipitacion_limpia_mm']
    ]
    .groupby('station_id')['precipitacion_limpia_mm'].quantile(.99)
    .rename('p99_positivo_train')
)
datos_trabajo['p99_positivo_train'] = datos_trabajo.station_id.map(p99_train)

sensibilidad_mes_extremo = []
for multiplicador in [2.5, 3.0, 3.5, 4.0]:
    supera = datos_trabajo.precipitacion_limpia_mm.gt(
        multiplicador * datos_trabajo.p99_positivo_train
    )
    conteo = supera.groupby([
        datos_trabajo.station_id, datos_trabajo.year, datos_trabajo.month
    ]).transform('sum')
    for minimo_dias in [3, 5, 7]:
        filas = MASCARA_TRAIN & conteo.ge(minimo_dias)
        meses = datos_trabajo.loc[filas, ['station_id', 'year', 'month']].drop_duplicates()
        sensibilidad_mes_extremo.append({
            'multiplicador_p99_train': multiplicador,
            'min_dias_extremos_mes': minimo_dias,
            'meses_marcados_train': int(len(meses)),
            'estaciones_train': int(meses.station_id.nunique()),
            'seleccionado': np.isclose(multiplicador, 3.0) and minimo_dias == 5,
        })
sensibilidad_mes_extremo = pd.DataFrame(sensibilidad_mes_extremo)

supera_3p99 = datos_trabajo.precipitacion_limpia_mm.gt(3.0 * datos_trabajo.p99_positivo_train)
conteo_3p99_mes = supera_3p99.groupby([
    datos_trabajo.station_id, datos_trabajo.year, datos_trabajo.month
]).transform('sum')
meses_extremos = datos_trabajo.loc[
    conteo_3p99_mes.ge(5), CLAVE_MES
].drop_duplicates()
datos_trabajo['mes_anomalo_p99'] = claves_en_tabla(datos_trabajo, CLAVE_MES, meses_extremos)
datos_trabajo.loc[datos_trabajo.mes_anomalo_p99, 'precipitacion_limpia_mm'] = np.nan
traza_calidad.extend(resumen_estado('03 mes con >=5 días sobre 3*p99 train -> NaN'))

# 3D. Extremo aislado: regla relativa de Jacobo, con p99 congelado de entrenamiento.
sensibilidad_digitacion_relativa = []
for multiplicador in [3.0, 4.0, 5.0]:
    for piso_mm in [100, 150, 200]:
        seleccion = (
            MASCARA_TRAIN
            & datos_trabajo.precipitacion_limpia_mm.gt(piso_mm)
            & datos_trabajo.precipitacion_limpia_mm.gt(multiplicador * datos_trabajo.p99_positivo_train)
        )
        sensibilidad_digitacion_relativa.append({
            'multiplicador_p99_train': multiplicador,
            'piso_mm': piso_mm,
            'dias_marcados_train': int(seleccion.sum()),
            'estaciones_train': int(datos_trabajo.loc[seleccion, 'station_id'].nunique()),
            'pares_train': int(datos_trabajo.loc[seleccion, ['station_id', 'year']].drop_duplicates().shape[0]),
            'seleccionado': np.isclose(multiplicador, 4.0) and piso_mm == 150,
        })
sensibilidad_digitacion_relativa = pd.DataFrame(sensibilidad_digitacion_relativa)

datos_trabajo['error_digitacion_relativo'] = (
    datos_trabajo.precipitacion_limpia_mm.gt(150)
    & datos_trabajo.precipitacion_limpia_mm.gt(4.0 * datos_trabajo.p99_positivo_train)
)
datos_trabajo.loc[datos_trabajo.error_digitacion_relativo, 'precipitacion_limpia_mm'] = np.nan
traza_calidad.extend(resumen_estado('04 extremo aislado >150 y >4*p99 train -> NaN'))

def resaltar_seleccion(fila):
    color = 'background-color: #d9ead3; font-weight: bold' if bool(fila['seleccionado']) else ''
    return [color] * len(fila)

display(Markdown('**Sensibilidad del posible total mensual (solo entrenamiento):**'))
display(sensibilidad_total_mensual.style.apply(resaltar_seleccion, axis=1))
display(Markdown('**Sensibilidad del mes con repetición extrema (solo entrenamiento):**'))
display(sensibilidad_mes_extremo.style.apply(resaltar_seleccion, axis=1))
display(Markdown('**Sensibilidad del extremo aislado (solo entrenamiento):**'))
display(sensibilidad_digitacion_relativa.style.apply(resaltar_seleccion, axis=1))

impacto_totales = pd.concat([
    meses_unico_registro.assign(regla='Mes con un único registro'),
    meses_total_mensual[CLAVE_MES].assign(regla='Posible total mensual'),
    meses_extremos.assign(regla='Mes con repetición extrema'),
], ignore_index=True)
resumen_impacto_totales = (
    impacto_totales.groupby(['regla', 'particion'], observed=True)
    .agg(meses=('month', 'size'), estaciones=('station_id', 'nunique'))
    .reset_index()
)
resumen_dias_digitacion = (
    datos_trabajo.loc[datos_trabajo.error_digitacion_relativo]
    .groupby('particion', observed=True)
    .agg(dias=('row_id', 'size'), estaciones=('station_id', 'nunique'))
    .reset_index()
)
display(Markdown('**Aplicación congelada a cada partición:**'))
display(resumen_impacto_totales)
display(Markdown('**Días aislados invalidados por la regla relativa:**'))
display(resumen_dias_digitacion)
display(Markdown(
    '**Decisión:** se adopta la estructura de Jacobo, pero las referencias p99 y medianas provienen solo de entrenamiento. '
    'Las lluvias altas que no cumplen estas condiciones se conservan.'
))
''')

copias_md = md(r'''
# Decisión 4: años copiados y posibles ceros falsos

**Años copiados.** Se comparan los 365 valores diarios de cada año, excluyendo el 29 de febrero para usar el mismo dominio anual. Solo una coincidencia exacta y completa entre estaciones distintas activa la regla. Como no puede saberse cuál estación es la original, ambos pares estación–año pasan a faltante.

**Ceros falsos.** Se adapta la regla de Jacobo: una racha de al menos 30 ceros consecutivos se invalida cuando, durante esa racha, al menos 40% de las otras estaciones observadas de la misma macrozona registran lluvia en promedio. El 40% se presenta junto con una tabla de sensibilidad calculada solo en entrenamiento. Es una regla agresiva: no reemplaza el cero por lluvia regional ni satelital; lo convierte en faltante y luego vuelve a evaluar cobertura.
''')

copias_code = code(r'''
import hashlib

# 4A. Años diarios copiados exactamente entre estaciones.
base_copias = datos_trabajo.loc[
    datos_trabajo.fecha.notna()
    & datos_trabajo.station_id.notna()
    & ~((datos_trabajo.fecha.dt.month == 2) & (datos_trabajo.fecha.dt.day == 29)),
    ['station_id', 'year', 'particion', 'fecha', 'precipitacion_original_mm']
].drop_duplicates(['station_id', 'year', 'fecha'], keep=False)

firmas = []
vectores = {}
for (estacion, anio, particion), grupo in base_copias.groupby(
    ['station_id', 'year', 'particion'], observed=True
):
    grupo = grupo.sort_values('fecha')
    if len(grupo) != 365 or grupo.precipitacion_original_mm.isna().any():
        continue
    valores = grupo.precipitacion_original_mm.to_numpy(dtype='<f8')
    firma = hashlib.sha256(valores.tobytes()).hexdigest()
    firmas.append({
        'station_id': estacion, 'year': int(anio), 'particion': particion,
        'firma': firma,
    })
    vectores[(estacion, int(anio))] = valores

firmas = pd.DataFrame(firmas)
grupos_copia = []
if not firmas.empty:
    for (anio, particion, firma), grupo in firmas.groupby(['year', 'particion', 'firma'], observed=True):
        estaciones = sorted(grupo.station_id.unique().tolist())
        if len(estaciones) < 2:
            continue
        referencia = vectores[(estaciones[0], int(anio))]
        iguales = [e for e in estaciones if np.array_equal(vectores[(e, int(anio))], referencia)]
        if len(iguales) >= 2:
            grupos_copia.append({
                'year': int(anio), 'particion': particion,
                'estaciones': ', '.join(iguales), 'n_estaciones': len(iguales),
                'firma_corta': firma[:12],
            })
grupos_copia = pd.DataFrame(grupos_copia)

pares_copiados = []
for fila in grupos_copia.itertuples(index=False):
    for estacion in fila.estaciones.split(', '):
        pares_copiados.append({
            'station_id': estacion, 'year': fila.year, 'particion': fila.particion
        })
pares_copiados = pd.DataFrame(pares_copiados, columns=['station_id', 'year', 'particion']).drop_duplicates()
datos_trabajo['anio_copiado'] = claves_en_tabla(
    datos_trabajo, ['station_id', 'year', 'particion'], pares_copiados
)
datos_trabajo.loc[datos_trabajo.anio_copiado, 'precipitacion_limpia_mm'] = np.nan
traza_calidad.extend(resumen_estado('05 año copiado exacto -> par a NaN'))

display(Markdown('**Grupos de años copiados encontrados por coincidencia exacta:**'))
display(grupos_copia if not grupos_copia.empty else pd.DataFrame({'resultado': ['No se detectaron grupos']}))
display(Markdown('**Pares afectados por partición:**'))
display(
    pares_copiados.groupby('particion').size().rename('pares').reset_index()
    if not pares_copiados.empty else pd.DataFrame({'resultado': ['Ningún par afectado']})
)

# 4B. Ceros falsos: misma aproximación regional de Jacobo, sin usar satélite.
datos_trabajo['macrozona'] = np.where(
    datos_trabajo.station_id.str.startswith('26', na=False), 'CAUCA', 'PACIFICO'
)
validos = datos_trabajo.precipitacion_limpia_mm.notna() & datos_trabajo.fecha.notna()
estado_diario_region = (
    datos_trabajo.loc[validos]
    .assign(lluvia=lambda x: x.precipitacion_limpia_mm.gt(0).astype(int))
    .groupby(['fecha', 'macrozona'], observed=True)
    .agg(estaciones_observadas=('station_id', 'nunique'), estaciones_con_lluvia=('lluvia', 'sum'))
    .reset_index()
)
ceros = datos_trabajo.loc[
    validos & datos_trabajo.precipitacion_limpia_mm.eq(0),
    ['row_id', 'station_id', 'year', 'particion', 'fecha', 'macrozona']
].merge(estado_diario_region, on=['fecha', 'macrozona'], how='left')
ceros['proporcion_otras_con_lluvia'] = (
    ceros.estaciones_con_lluvia / (ceros.estaciones_observadas - 1).replace(0, np.nan)
)

rachas = []
for (estacion, anio, particion, macrozona), grupo in ceros.groupby(
    ['station_id', 'year', 'particion', 'macrozona'], observed=True
):
    grupo = grupo.sort_values('fecha').copy()
    grupo['grupo_racha'] = grupo.fecha.diff().dt.days.ne(1).cumsum()
    for _, racha in grupo.groupby('grupo_racha'):
        if len(racha) < 30:
            continue
        rachas.append({
            'station_id': estacion, 'year': int(anio), 'particion': particion,
            'macrozona': macrozona, 'inicio': racha.fecha.min(), 'fin': racha.fecha.max(),
            'longitud_dias': int(len(racha)),
            'proporcion_regional_lluvia': float(racha.proporcion_otras_con_lluvia.mean()),
        })
rachas_cero = pd.DataFrame(rachas)

sensibilidad_ceros = []
for umbral in [0.20, 0.30, 0.40, 0.50, 0.60]:
    seleccion = (
        rachas_cero.particion.eq('ENTRENAMIENTO')
        & rachas_cero.proporcion_regional_lluvia.ge(umbral)
    ) if not rachas_cero.empty else pd.Series(dtype=bool)
    elegidas = rachas_cero.loc[seleccion] if not rachas_cero.empty else rachas_cero
    sensibilidad_ceros.append({
        'umbral_proporcion_regional': umbral,
        'rachas_marcadas_train': int(len(elegidas)),
        'dias_cero_marcados_train': int(elegidas.longitud_dias.sum()) if len(elegidas) else 0,
        'estaciones_train': int(elegidas.station_id.nunique()) if len(elegidas) else 0,
        'pares_train': int(elegidas[['station_id', 'year']].drop_duplicates().shape[0]) if len(elegidas) else 0,
        'seleccionado': np.isclose(umbral, .40),
    })
sensibilidad_ceros = pd.DataFrame(sensibilidad_ceros)
sensibilidad_ceros['dias_recuperados_vs_umbral_anterior'] = (
    sensibilidad_ceros.dias_cero_marcados_train.shift(1)
    - sensibilidad_ceros.dias_cero_marcados_train
).fillna(0).astype(int)

UMBRAL_CERO_REGIONAL = .40
rachas_cero_marcadas = rachas_cero.loc[
    rachas_cero.proporcion_regional_lluvia.ge(UMBRAL_CERO_REGIONAL)
].copy() if not rachas_cero.empty else rachas_cero.copy()

filas_cero_falso = set()
for racha in rachas_cero_marcadas.itertuples(index=False):
    ids = ceros.loc[
        ceros.station_id.eq(racha.station_id)
        & ceros.year.eq(racha.year)
        & ceros.fecha.between(racha.inicio, racha.fin),
        'row_id'
    ]
    filas_cero_falso.update(ids.astype(int).tolist())
datos_trabajo['cero_falso_regional'] = datos_trabajo.row_id.isin(filas_cero_falso)
datos_trabajo.loc[datos_trabajo.cero_falso_regional, 'precipitacion_limpia_mm'] = np.nan
traza_calidad.extend(resumen_estado('06 racha cero >=30 y apoyo regional >=40% -> NaN'))

display(Markdown('**Sensibilidad del umbral regional, calculada solo en entrenamiento:**'))
display(
    sensibilidad_ceros.style.apply(resaltar_seleccion, axis=1)
    .format({'umbral_proporcion_regional': '{:.0%}'})
)
impacto_ceros_particion = (
    rachas_cero_marcadas.groupby('particion', observed=True)
    .agg(
        rachas=('longitud_dias', 'size'),
        dias_cero_invalidados=('longitud_dias', 'sum'),
        estaciones=('station_id', 'nunique'),
    ).reset_index()
    if not rachas_cero_marcadas.empty else pd.DataFrame({'resultado': ['Ninguna racha marcada']})
)
display(Markdown('**Aplicación del umbral congelado por partición:**'))
display(impacto_ceros_particion)

trazabilidad_calidad = pd.DataFrame(traza_calidad)
trazabilidad_calidad['registros_excluidos_vs_etapa_anterior'] = (
    trazabilidad_calidad.groupby('particion').registros_validos.shift(1)
    - trazabilidad_calidad.registros_validos
).fillna(0).astype(int)
display(Markdown('**Embudo completo antes de evaluar cobertura:**'))
display(
    trazabilidad_calidad.style
    .background_gradient(subset=['registros_excluidos_vs_etapa_anterior'], cmap='Reds')
)
display(Markdown(
    '**Lectura:** validación y prueba aparecen únicamente para medir el impacto de reglas ya congeladas. '
    'No participaron en medianas, percentiles ni en la elección del 40%.'
))
''')

# Insertar las nuevas decisiones después de la decisión 2.
cells[26:26] = [totales_md, totales_code, copias_md, copias_code]

# Renumerar las decisiones posteriores.
for cell in cells:
    if cell.get("cell_type") != "markdown":
        continue
    source = "".join(cell.get("source", []))
    source = source.replace("# Decisión 3: cobertura anual", "# Decisión 5: cobertura anual")
    source = source.replace("# Decisión 4: mínimo de cinco años", "# Decisión 6: mínimo de cinco años")
    cell["source"] = source

# Añadir un cierre explícito con el resultado final por partición.
cells.extend([
    md(r'''
# Resultado final del preprocesamiento terrestre

El resultado conserva los valores originales, documenta cada invalidación y aplica en este orden: errores confirmados, reglas adaptadas de Jacobo, cobertura/continuidad y mínimo de años de entrenamiento. No se imputa precipitación en esta fase.
'''),
    code(r'''
resumen_final = (
    pares_validos_minimo_5.groupby('particion', observed=True)
    .agg(
        estaciones=('station_id', 'nunique'),
        pares_estacion_anio=('year', 'size'),
    )
    .reindex(ORDEN_PARTICIONES)
    .reset_index()
)
resumen_final.loc[len(resumen_final)] = [
    'TOTAL', pares_validos_minimo_5.station_id.nunique(), len(pares_validos_minimo_5)
]
display(Markdown('**Pares elegibles después de todas las decisiones:**'))
display(resumen_final)
display(Markdown(
    '**Salida metodológica:** `datos_trabajo` mantiene la precipitación original, la precipitación limpia y '
    'una bandera por motivo; `pares_validos_minimo_5` identifica las curvas que pueden avanzar a FDA. '
    'La prueba no se usó para seleccionar ninguna regla.'
))
'''),
])

nb["cells"] = cells
with RUTA.open("w", encoding="utf-8") as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)

print(f"Actualizado: {RUTA}")

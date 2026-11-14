from pathlib import Path
import json


ROOT = Path(__file__).resolve().parents[1]
notebooks = list(ROOT.glob("*.ipynb")) + list((ROOT / "KAN").glob("*.ipynb"))

reemplazos = {
    'ROOT / "P_1930_1979.xlsx"': 'ROOT / "Datos" / "P_1930_1979.xlsx"',
    'ROOT / "P_1980_1990.csv"': 'ROOT / "Datos" / "P_1980_1990.csv"',
    'ROOT / "P_1991_2003.csv"': 'ROOT / "Datos" / "P_1991_2003.csv"',
    "ROOT / 'P_1930_1979.xlsx'": "ROOT / 'Datos' / 'P_1930_1979.xlsx'",
    "ROOT / 'P_1980_1990.csv'": "ROOT / 'Datos' / 'P_1980_1990.csv'",
    "ROOT / 'P_1991_2003.csv'": "ROOT / 'Datos' / 'P_1991_2003.csv'",
    'ROOT / "ImagenesSatelitePentadal"': 'ROOT / "Datos" / "ImagenesSatelitePentadal"',
    'ROOT / "ImagenesSatelitePentadal.zip"': 'ROOT / "Datos" / "ImagenesSatelitePentadal.zip"',
    "ROOT / 'ImagenesSatelitePentadal'": "ROOT / 'Datos' / 'ImagenesSatelitePentadal'",
    "ROOT / 'ImagenesSatelitePentadal.zip'": "ROOT / 'Datos' / 'ImagenesSatelitePentadal.zip'",
    "z.extractall(ROOT)": "z.extractall(ROOT / 'Datos')",
    "CATALOGO = ROOT / 'fda_fase1' / 'catalogo_localizacion_docente_preparado.csv'": "RESULTADOS_KAN = ROOT / 'KAN' / 'resultados'\nCATALOGO = RESULTADOS_KAN / 'catalogo_localizacion_docente_preparado.csv'",
    "TRAZABILIDAD = ROOT / 'fda_fase1' / 'pares_estacion_pixel_docente.csv'": "TRAZABILIDAD = RESULTADOS_KAN / 'pares_estacion_pixel_docente.csv'",
    "IDENTIFICADORES = ROOT / 'fda_fase2' / 'identificadores_pares_kan_docente.csv'": "IDENTIFICADORES = RESULTADOS_KAN / 'identificadores_pares_kan_docente.csv'",
    "RESULTADOS = ROOT / 'fda_fase2' / 'piloto_kan_docente'": "RESULTADOS = RESULTADOS_KAN / 'piloto_kan_docente'",
    "fda_fase2/piloto_kan_docente/": "KAN/resultados/piloto_kan_docente/",
}

for ruta in notebooks:
    nb = json.loads(ruta.read_text(encoding="utf-8"))
    cambios = 0
    for celda in nb.get("cells", []):
        fuente = "".join(celda.get("source", []))
        nueva = fuente
        for anterior, posterior in reemplazos.items():
            nueva = nueva.replace(anterior, posterior)

        # Notebooks exploratorios antiguos verifican una lista de rutas relativas.
        if "required = [" in nueva:
            for nombre in [
                "P_1930_1979.xlsx", "P_1980_1990.csv", "P_1991_2003.csv",
                "ImagenesSatelitePentadal.zip",
            ]:
                nueva = nueva.replace(f"'{nombre}'", f"'Datos/{nombre}'")

        # El notebook KAN puede abrirse desde la raíz o desde su propia carpeta.
        if ruta.parent.name == "KAN":
            if "if not (ROOT / 'Datos').exists():" not in nueva:
                nueva = nueva.replace(
                    "ROOT = Path.cwd()\n",
                    "ROOT = Path.cwd()\nif not (ROOT / 'Datos').exists():\n    ROOT = ROOT.parent\n",
                )
            nueva = nueva.replace(
                "%run preparar_emparejamiento_docente.py",
                "script_emparejamiento = ROOT / 'KAN' / 'preparar_emparejamiento_docente.py'\n%run $script_emparejamiento",
            )

        if nueva != fuente:
            celda["source"] = nueva
            cambios += 1
    if cambios:
        ruta.write_text(json.dumps(nb, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"{ruta.relative_to(ROOT)}: {cambios} celdas actualizadas")

# Los tres EDA heredados trabajaban exclusivamente con el ZIP. Ahora leen la
# carpeta extraída para no mantener dos copias de las 2.328 imágenes.
for nombre in [
    "EDA_completo_faltantes_precipitacion.ipynb",
    "EDA_precipitacion_valle_cauca copy.ipynb",
    "EDA_precipitacion_valle_cauca.ipynb",
]:
    ruta = ROOT / nombre
    nb = json.loads(ruta.read_text(encoding="utf-8"))
    for celda in nb["cells"]:
        if celda.get("cell_type") != "code":
            continue
        fuente = "".join(celda.get("source", []))
        if fuente.strip() == "pip install --upgrade pip":
            fuente = "%pip install --upgrade pip"
        fuente = fuente.replace(
            "'Datos/ImagenesSatelitePentadal.zip'",
            "'Datos/ImagenesSatelitePentadal'",
        )
        if "sat_zip = ROOT / 'Datos' / 'ImagenesSatelitePentadal.zip'" in fuente:
            usa_bandera = "PLOT_AVAILABLE" in fuente
            bloque_plot = (
                "if PLOT_AVAILABLE:\n"
                "    fig, ax = plt.subplots(figsize=(11, 3.5))\n"
                "    ax.bar(sat_year.year, sat_year.tiff_count, color=np.where(sat_year.complete_72_expected, 'steelblue', 'darkorange'))\n"
                "    ax.axhline(72, color='black', linestyle='--', linewidth=1, label='72 periodos/año según patrón')\n"
                "    ax.set(title='Cantidad de TIFF por año', xlabel='Año', ylabel='Número de TIFF'); ax.legend(); ax.grid(axis='y', alpha=.2); plt.show()\n"
                if usa_bandera else
                "fig, ax = plt.subplots(figsize=(11, 3.5))\n"
                "ax.bar(sat_year.year, sat_year.tiff_count, color=np.where(sat_year.complete_72_expected, 'steelblue', 'darkorange'))\n"
                "ax.axhline(72, color='black', linestyle='--', linewidth=1, label='72 periodos/año según patrón')\n"
                "ax.set(title='Cantidad de TIFF por año', xlabel='Año', ylabel='Número de TIFF'); ax.legend(); ax.grid(axis='y', alpha=.2); plt.show()\n"
            )
            fuente = f"""sat_dir = ROOT / 'Datos' / 'ImagenesSatelitePentadal'
tiff_paths = sorted(p for p in sat_dir.glob('ppt*.tif') if not p.name.startswith('._'))
tiffs = [p.name for p in tiff_paths]
year_counts = Counter(p.stem[3:7] for p in tiff_paths)
sat_year = pd.DataFrame([{{'year': int(y), 'tiff_count': n, 'complete_72_expected': n == 72}}
                         for y, n in sorted(year_counts.items())])
print('TIFF de datos:', len(tiffs))
display(sat_year)

{bloque_plot}
# Intento opcional de leer metadatos de un TIFF.
try:
    import rasterio
    with rasterio.open(tiff_paths[0]) as ds:
        print({{'archivo': tiff_paths[0].name, 'ancho': ds.width, 'alto': ds.height,
               'bandas': ds.count, 'CRS': str(ds.crs), 'nodata': ds.nodata,
               'tipo': ds.dtypes, 'transform': ds.transform, 'bounds': ds.bounds}})
        band = ds.read(1, masked=True)
        print('mínimo/máximo válidos:', float(band.min()), float(band.max()))
except ImportError:
    print('rasterio no está instalado: no se leyeron metadatos ni valores de píxel.')
except Exception as e:
    print('No se pudo leer el TIFF con rasterio:', repr(e))"""
        celda["source"] = fuente
    ruta.write_text(json.dumps(nb, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"{nombre}: lectura satelital cambiada de ZIP a carpeta")

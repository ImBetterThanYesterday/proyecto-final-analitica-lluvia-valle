from pathlib import Path
import json


ROOT = Path(__file__).resolve().parents[1]

DEPENDENCIAS = {
    "01_EDA_calidad_datos.ipynb": {
        "pandas": "pandas", "numpy": "numpy", "matplotlib": "matplotlib",
        "openpyxl": "openpyxl", "rasterio": "rasterio",
    },
    "01b_valores_unicos.ipynb": {
        "pandas": "pandas", "numpy": "numpy", "matplotlib": "matplotlib",
    },
    "02_datos_externos.ipynb": {
        "pandas": "pandas", "numpy": "numpy", "matplotlib": "matplotlib",
        "rasterio": "rasterio",
    },
    "03_construccion_funcional.ipynb": {
        "pandas": "pandas", "numpy": "numpy", "matplotlib": "matplotlib",
        "rasterio": "rasterio", "scipy": "scipy",
    },
    "04_eda_funcional.ipynb": {
        "pandas": "pandas", "numpy": "numpy", "matplotlib": "matplotlib",
        "scipy": "scipy", "skfda": "scikit-fda",
    },
    "KAN/03_red_funcional_kan_arquitectura.ipynb": {
        "pandas": "pandas", "numpy": "numpy", "matplotlib": "matplotlib",
        "torch": "torch",
    },
    "KAN/04_emparejamiento_docente_y_piloto_kan.ipynb": {
        "pandas": "pandas", "numpy": "numpy", "matplotlib": "matplotlib",
        "openpyxl": "openpyxl", "torch": "torch",
    },
}


def celda_dependencias(paquetes):
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {"tags": ["verificacion-dependencias"]},
        "outputs": [],
        "source": f'''# Verificación reproducible: no reinstala paquetes presentes ni instala sin autorización.
from importlib.util import find_spec
import subprocess
import sys

INSTALAR_FALTANTES = False  # Cambiar a True únicamente si se autoriza instalar en este entorno.
PAQUETES_REQUERIDOS = {paquetes!r}
faltantes = {{modulo: paquete for modulo, paquete in PAQUETES_REQUERIDOS.items()
             if find_spec(modulo) is None}}

if faltantes and INSTALAR_FALTANTES:
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', *faltantes.values()])
    faltantes = {{modulo: paquete for modulo, paquete in PAQUETES_REQUERIDOS.items()
                 if find_spec(modulo) is None}}

if faltantes:
    comando = f"{{sys.executable}} -m pip install {{' '.join(faltantes.values())}}"
    raise ModuleNotFoundError(
        f"Faltan dependencias: {{faltantes}}. No se instalaron automáticamente. "
        f"Autorice la instalación o ejecute: {{comando}}"
    )

versiones = {{modulo: __import__(modulo).__version__
             for modulo in PAQUETES_REQUERIDOS
             if hasattr(__import__(modulo), '__version__')}}
print('Dependencias verificadas:', versiones)
''',
    }


for relativa, paquetes in DEPENDENCIAS.items():
    ruta = ROOT / relativa
    nb = json.loads(ruta.read_text(encoding="utf-8"))
    nb["cells"] = [
        c for c in nb["cells"]
        if "verificacion-dependencias" not in c.get("metadata", {}).get("tags", [])
    ]
    posicion = 1 if nb["cells"] and nb["cells"][0].get("cell_type") == "markdown" else 0
    nb["cells"].insert(posicion, celda_dependencias(paquetes))
    for celda in nb["cells"]:
        if celda.get("cell_type") == "code":
            fuente = "".join(celda.get("source", []))
            fuente = "\n".join(
                linea for linea in fuente.splitlines()
                if not linea.strip().startswith("# %pip install")
            )
            celda["source"] = fuente + ("\n" if fuente else "")
    ruta.write_text(json.dumps(nb, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"Actualizado: {relativa}")

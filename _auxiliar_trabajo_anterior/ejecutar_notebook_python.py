"""Ejecuta un notebook lineal con Python y guarda sus salidas en el archivo.

Este proyecto no depende del servidor Jupyter: captura texto, tablas HTML y figuras.
"""

from pathlib import Path
import base64
import contextlib
import io
import json
import sys
import traceback

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


if len(sys.argv) != 2:
    raise SystemExit("Uso: python ejecutar_notebook_python.py RUTA.ipynb")

ruta = Path(sys.argv[1]).resolve()
nb = json.loads(ruta.read_text(encoding="utf-8"))
salidas_actuales = None


class Markdown:
    def __init__(self, data):
        self.data = str(data)


def display(*objetos):
    for objeto in objetos:
        if isinstance(objeto, Markdown):
            datos = {"text/markdown": objeto.data, "text/plain": objeto.data}
        elif hasattr(objeto, "_repr_html_"):
            try:
                html = objeto._repr_html_()
            except Exception:
                html = None
            datos = {"text/plain": repr(objeto)}
            if html is not None:
                datos["text/html"] = html
        else:
            datos = {"text/plain": repr(objeto)}
        salidas_actuales.append({
            "output_type": "display_data", "data": datos, "metadata": {}
        })


def mostrar_figuras(*args, **kwargs):
    for numero in plt.get_fignums():
        figura = plt.figure(numero)
        buffer = io.BytesIO()
        figura.savefig(buffer, format="png", dpi=110, bbox_inches="tight")
        salidas_actuales.append({
            "output_type": "display_data",
            "data": {"image/png": base64.b64encode(buffer.getvalue()).decode("ascii")},
            "metadata": {},
        })
    plt.close("all")


plt.show = mostrar_figuras
entorno = {
    "__name__": "__main__",
    "display": display,
    "Markdown": Markdown,
}
ejecucion = 0

for indice, celda in enumerate(nb["cells"]):
    if celda.get("cell_type") != "code":
        continue
    ejecucion += 1
    celda["execution_count"] = ejecucion
    celda["outputs"] = []
    salidas_actuales = celda["outputs"]
    fuente = "".join(celda.get("source", []))
    lineas = []
    for linea in fuente.splitlines():
        if linea.lstrip().startswith("%"):
            continue
        if linea.strip().startswith("from IPython.display import"):
            continue
        lineas.append(linea)
    fuente_ejecutable = "\n".join(lineas)
    primera = next((x.strip() for x in lineas if x.strip()), "")
    print(f"Celda {indice}: {primera[:80]}", flush=True)
    salida_texto = io.StringIO()
    try:
        with contextlib.redirect_stdout(salida_texto), contextlib.redirect_stderr(salida_texto):
            exec(compile(fuente_ejecutable, f"{ruta.name}:celda_{indice}", "exec"), entorno)
    except Exception as exc:
        texto = salida_texto.getvalue()
        if texto:
            celda["outputs"].insert(0, {"output_type": "stream", "name": "stdout", "text": texto})
        celda["outputs"].append({
            "output_type": "error",
            "ename": type(exc).__name__,
            "evalue": str(exc),
            "traceback": traceback.format_exc().splitlines(),
        })
        ruta.write_text(json.dumps(nb, ensure_ascii=False, indent=1), encoding="utf-8")
        raise
    texto = salida_texto.getvalue()
    if texto:
        celda["outputs"].insert(0, {"output_type": "stream", "name": "stdout", "text": texto})
    ruta.write_text(json.dumps(nb, ensure_ascii=False, indent=1), encoding="utf-8")

print(f"Notebook ejecutado y guardado: {ruta}")

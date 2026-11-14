"""Ejecutor mínimo de notebooks usando el kernel local disponible.

Se usa porque este equipo tiene ipykernel/jupyter_client, pero no nbconvert.
Conserva stream, tablas HTML, texto, errores e imágenes en las salidas del .ipynb.
"""

from pathlib import Path
import json
import sys
from queue import Empty

from jupyter_client import KernelManager


if len(sys.argv) != 2:
    raise SystemExit("Uso: python ejecutar_notebook_local.py RUTA.ipynb")

ruta = Path(sys.argv[1]).resolve()
nb = json.loads(ruta.read_text(encoding="utf-8"))

km = KernelManager(kernel_name=nb.get("metadata", {}).get("kernelspec", {}).get("name", "python3"))
km.start_kernel(cwd=str(ruta.parent))
kc = km.client()
kc.start_channels()
kc.wait_for_ready(timeout=60)

contador = 0
try:
    for indice, celda in enumerate(nb["cells"]):
        if celda.get("cell_type") != "code":
            continue
        fuente = "".join(celda.get("source", []))
        contador += 1
        celda["execution_count"] = contador
        celda["outputs"] = []
        primera = next((x.strip() for x in fuente.splitlines() if x.strip()), "")
        print(f"Celda {indice}: {primera[:80]}", flush=True)
        msg_id = kc.execute(fuente, stop_on_error=True)
        fallo = None
        while True:
            try:
                msg = kc.get_iopub_msg(timeout=1800)
            except Empty as exc:
                raise TimeoutError(f"La celda {indice} superó el tiempo límite") from exc
            if msg.get("parent_header", {}).get("msg_id") != msg_id:
                continue
            tipo = msg["header"]["msg_type"]
            contenido = msg["content"]
            if tipo == "status" and contenido.get("execution_state") == "idle":
                break
            if tipo == "stream":
                celda["outputs"].append({
                    "output_type": "stream",
                    "name": contenido["name"],
                    "text": contenido["text"],
                })
            elif tipo in {"display_data", "execute_result"}:
                salida = {
                    "output_type": tipo,
                    "data": contenido.get("data", {}),
                    "metadata": contenido.get("metadata", {}),
                }
                if tipo == "execute_result":
                    salida["execution_count"] = contador
                celda["outputs"].append(salida)
            elif tipo == "error":
                fallo = contenido
                celda["outputs"].append({
                    "output_type": "error",
                    "ename": contenido.get("ename", "Error"),
                    "evalue": contenido.get("evalue", ""),
                    "traceback": contenido.get("traceback", []),
                })
        ruta.write_text(json.dumps(nb, ensure_ascii=False, indent=1), encoding="utf-8")
        if fallo:
            raise RuntimeError(
                f"Celda {indice}: {fallo.get('ename')}: {fallo.get('evalue')}"
            )
finally:
    kc.stop_channels()
    km.shutdown_kernel(now=True)

print(f"Notebook ejecutado y guardado: {ruta}")

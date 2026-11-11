"""Inventario conservador del código R y cruce con dibujos reales del DOCX.

No ejecuta cálculos R. Un patrón gráfico nuevo no reconocido falla cerrado.
El log previo se contrasta como evidencia, nunca como fuente única del inventario.
"""
from pathlib import Path
from zipfile import ZipFile
from collections import Counter
import csv
import hashlib
import json
import re
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "resultados_R"
NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main",
      "a": "http://schemas.openxmlformats.org/drawingml/2006/main",
      "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inventario():
    codigo = (ROOT / "proyecto_funcional.R").read_text()
    construccion = ["estaciones", "satelite"]
    eda = ["estaciones", "satelite_tierra"]
    # Estas listas son el dominio explícito de los bucles del recorrido vigente.
    for literal in ('c("ECM", "VC", "VCG")', 'c("varianza", "de", "cv")',
                    'c("covarianza", "correlacion")', 'c("ONI", "grupos")',
                    'c("estaciones", "satelite")', '"satelite_tierra", oni)'):
        assert literal in codigo, f"Revisar dominios del inventario R: {literal}"
    patrones = {
        '"criterio_", crit, "_", nombre, ".pdf"':
            ([f"criterio_{c}" for c in ("ECM", "VC", "VCG")], construccion),
        '"seleccion_K_", nombre, ".pdf"': (["seleccion_K"], construccion),
        '"curva_banda_residuos_", nombre, ".pdf"': (["curva_banda_residuos"], construccion),
        '"residual_", nombre, ".pdf"': (["residual"], construccion),
        '"media_", nombre, ".pdf"': (["media"], eda),
        '"EDA_", nombre, ".pdf"': (["EDA"], eda),
        '"boxplot_", nombre, ".pdf"': (["boxplot"], eda),
        'variable, "_", nombre, ".pdf"': (["varianza", "de", "cv", "covarianza", "correlacion"], eda),
        'tipo, "_", nombre, ".pdf"': (["ONI", "grupos"], eda),
        '"contexto_", nombre, ".pdf"': (["contexto"], eda),
        '"observados_", fuente, ".pdf"': (["observados"], construccion),
        '"bases_", fuente, ".pdf"': (["bases"], construccion),
        '"derivada_", nombre, ".pdf"': (["derivada"], eda),
    }
    llamadas = list(re.finditer(r'exportar_grafico\(paste0\(([^\n]+?)\), function\(\)', codigo))
    assert len(llamadas) == len(re.findall(r'exportar_grafico\(', codigo)), "Exportación R nueva sin inventariar"
    filas = {}
    for llamada in llamadas:
        patron = llamada.group(1)
        assert patron in patrones, f"Patrón gráfico R nuevo: {patron}"
        prefijos, fuentes = patrones[patron]
        for prefijo in prefijos:
            for fuente in fuentes:
                nombre = f"{prefijo}_{fuente}"
                linea = str(codigo[:llamada.start()].count("\n") + 1)
                if nombre in filas:
                    filas[nombre]["linea_R"] += ";" + linea
                    continue
                componentes = ([f"criterio_{c}_{fuente}" for c in ("ECM", "VC", "VCG")]
                               if prefijo == "seleccion_K" else
                               [f"grupos_{fuente}", f"ONI_{fuente}"] if prefijo == "contexto" else [])
                filas[nombre] = {"figura_R": nombre, "linea_R": linea, "fuente": fuente,
                                 "componentes": componentes}
    log = (OUT / "prueba_rscript_directo.log").read_text()
    evidencia = re.findall(r"Figura guardada: resultados_R/(\S+)\.pdf", log)
    assert "Recorrido terminado" in log and Counter(evidencia) == Counter(filas), "Log y código R no concuerdan"
    return filas


def dibujos(docx):
    with ZipFile(docx) as z:
        assert z.testzip() is None
        doc = ET.fromstring(z.read("word/document.xml"))
        relaciones = {r.get("Id"): r.get("Target") for r in ET.fromstring(z.read("word/_rels/document.xml.rels"))}
        hashes = []
        for b in doc.findall(".//a:blip", NS):
            target = relaciones[b.get("{" + NS["r"] + "}embed")]
            hashes.append(hashlib.sha256(z.read("word/" + target)).hexdigest())
        textos = ["".join(t.text or "" for t in p.findall(".//w:t", NS)) for p in doc.findall(".//w:p", NS)]
        pies = [t for t in textos if re.match(r"Figuras? \d", t)]
        media = [n for n in z.namelist() if n.startswith("word/media/")]
        return {"drawings": len(doc.findall(".//w:drawing", NS)), "blips": len(hashes),
                "media": len(media), "media_unicos": len({hashlib.sha256(z.read(n)).hexdigest() for n in media}),
                "hashes_orden": hashes, "pies": pies, "sha256_docx": sha(docx)}


def correspondencia(figuras):
    invent = inventario()
    assert len(figuras) == len(set(figuras)), "Figura repetida en el generador"
    relevantes = {n for n, r in invent.items() if not r["componentes"]}
    assert set(figuras) == relevantes, f"Omisiones={relevantes-set(figuras)}; extras={set(figuras)-relevantes}"
    filas = []
    for orden, (nombre, r) in enumerate(invent.items(), 1):
        png, pdf = OUT / (nombre + ".png"), OUT / (nombre + ".pdf")
        assert png.is_file() and pdf.is_file(), nombre
        numeros = [figuras.index(c) + 1 for c in r["componentes"]] if r["componentes"] else [figuras.index(nombre) + 1]
        filas.append({"orden_inventario": orden, "figura_R": nombre, "linea_R": r["linea_R"],
                      "fuente": r["fuente"], "png": str(png.relative_to(ROOT)), "pdf": str(pdf.relative_to(ROOT)),
                      "sha256_png": sha(png), "sha256_pdf": sha(pdf),
                      "figura_Word": ";".join(map(str, numeros)),
                      "estado": "compatibilidad_no_reinsertar" if r["componentes"] else "incluida",
                      "paneles": len(r["componentes"]) or 1,
                      "equivalentes": ";".join(r["componentes"]),
                      "evidencia": "codigo_R_y_log_previo"})
    with (OUT / "correspondencia_figuras_R_Word.csv").open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(filas[0]))
        writer.writeheader(); writer.writerows(filas)
    return filas


def baseline():
    destino = OUT / "auditoria_imagenes_antes.json"
    if destino.exists():
        raise FileExistsError("La evidencia previa ya existe; no se sobrescribe")
    invent = inventario()
    doc = dibujos(ROOT / "INFORME_RECONSTRUCCION_EDA_ACTUAL.docx")
    mapeo = {sha(OUT / (n + ".png")): n for n in invent}
    nombres = [mapeo[h] for h in doc["hashes_orden"]]
    doc.update({"figuras_en_orden": nombres,
                "omitidas_individuales": sorted(n for n, r in invent.items() if not r["componentes"] and n not in nombres),
                "compuestas_no_incluidas": {n: r["componentes"] for n, r in invent.items() if r["componentes"]},
                "duplicadas": {n: c for n, c in Counter(nombres).items() if c > 1},
                "sha256_R": sha(ROOT / "proyecto_funcional.R"),
                "sha256_fuentes": {p.name: sha(p) for p in [ROOT / "generar_informe_actual.py", ROOT / "validar_informe_actual.py"]},
                "sha256_salidas_R": {p.name: sha(p) for p in OUT.iterdir() if p.suffix in (".csv", ".rds", ".pdf", ".png") and not p.name.startswith(("pagina-", "INFORME_"))},
                "reporte_previo": "28 no prueba cobertura; el validador fijaba 28",
                "conteo_manual_usuaria": "Word 27 / ejecución R 30; no reproducible en los artefactos actuales"})
    destino.write_text(json.dumps(doc, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({k: v for k, v in doc.items() if k in ("drawings", "media", "media_unicos", "omitidas_individuales", "duplicadas")}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    baseline()

"""Valida el entregable real, sus cifras, imágenes y PDF ya renderizado."""
from pathlib import Path
from zipfile import ZipFile
import csv
import hashlib
import json
import struct
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "resultados_R"
NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main",
      "wp": "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing",
      "a": "http://schemas.openxmlformats.org/drawingml/2006/main",
      "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}


def csv_real(nombre):
    with (OUT / nombre).open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def celdas(fila):
    return ["".join(c.itertext()) for c in fila.findall("w:tc", NS)]


def decimal(texto):
    return float(texto.replace(".", "").replace(",", "."))


def main():
    manifiesto = json.loads((OUT / "manifiesto_informe_actual.json").read_text())
    for nombre, esperado in manifiesto["archivos"].items():
        assert hashlib.sha256((ROOT / nombre).read_bytes()).hexdigest() == esperado, nombre
    with ZipFile(ROOT / manifiesto["documento"]) as z:
        assert z.testzip() is None
        doc = ET.fromstring(z.read("word/document.xml"))
        texto = " ".join(t.text or "" for t in doc.findall(".//w:t", NS))
        assert "desde la primera línea" in texto and "Source ejecuta todo" in texto
        assert not any(obsoleto in texto for obsoleto in ("expression", "eval(clase", "principal()", "sin su coma"))
        tablas = doc.findall(".//w:tbl", NS)
        assert len(tablas) == 16
        relaciones = ET.fromstring(z.read("word/_rels/document.xml.rels"))
        destino = {r.attrib["Id"]: r.attrib["Target"] for r in relaciones}
        for inline in doc.findall(".//wp:inline", NS):
            embed = inline.find(".//a:blip", NS).get("{" + NS["r"] + "}embed")
            datos = z.read("word/" + destino[embed])
            w, h = struct.unpack(">II", datos[16:24])
            ext = inline.find("wp:extent", NS)
            assert abs(int(ext.get("cx")) / int(ext.get("cy")) / (w / h) - 1) < .001
        medias = [z.read(n) for n in z.namelist() if n.startswith("word/media/")]
        usados = [ROOT / n for n in manifiesto["archivos"] if n.endswith(".png")]
        assert len(medias) == len(usados) == manifiesto["figuras"]
        assert {hashlib.sha256(b).hexdigest() for b in medias} == {
            hashlib.sha256(p.read_bytes()).hexdigest() for p in usados}
        # Toda figura exigida por la rúbrica debe estar incluida (no sólo exportada).
        figmap = {r["figura_R"]: r for r in csv_real("manifiesto_figuras.csv")}
        requeridas = [f"criterio_{c}_{f}" for c in ("ECM", "VC", "VCG")
                      for f in ("estaciones", "satelite")]
        requeridas += [f"{v}_{f}" for v in ("varianza", "de", "cv")
                       for f in ("estaciones", "satelite_tierra")]
        requeridas += [f"residual_{f}" for f in ("estaciones", "satelite")]
        requeridas += ["diagrama_flujo_general", "diagrama_decision_K",
                       "diagrama_decision_construccion", "diagrama_decision_NA",
                       "diagrama_decision_centralidad", "diagrama_decision_temporadas",
                       "diagrama_decision_fuentes"]
        for nombre in requeridas:
            assert figmap[nombre]["estado"] == "incluida", f"Figura exigida fuera del Word: {nombre}"
        # Cifras leídas del Word, contrastadas con los CSV de la ejecución.
        cobertura = csv_real("cobertura.csv")
        for fila, referencia in zip(tablas[0].findall("w:tr", NS)[1:], cobertura):
            cells = celdas(fila)
            assert cells[0] == referencia["fuente"]
            assert decimal(cells[1]) == int(referencia["curvas"])
            anual = csv_real("cobertura_anual_" + cells[0] + ".csv")
            assert decimal(cells[-1]) == sum(int(r["valores_faltantes"]) for r in anual)
        minimos = tablas[3].findall("w:tr", NS)[1:]
        assert len(minimos) == 6
        for fila in minimos:
            fuente, criterio, k, valor = celdas(fila)
            mejor = min(csv_real(f"seleccion_K_{fuente}.csv"), key=lambda r: float(r[criterio]))
            assert int(k) == int(mejor["K"])
            assert abs(decimal(valor) - float(mejor[criterio])) < .00051
    for entrada in csv_real("integridad_entradas.csv"):
        assert hashlib.md5((ROOT / entrada["archivo"]).read_bytes()).hexdigest() == entrada["md5"]
    # El PDF es opcional: en este Mac no hay LibreOffice/pdflatex para renderizarlo.
    # Se valida sólo si existe una versión más reciente que el DOCX.
    pdf = OUT / "INFORME_RECONSTRUCCION_EDA_ACTUAL.pdf"
    paginas_archivo = OUT / "inspeccion_paginas.json"
    pdf_estado = "no_regenerado_sin_LibreOffice"
    paginas_n = None
    if pdf.exists() and paginas_archivo.exists() and \
            pdf.stat().st_mtime >= (ROOT / manifiesto["documento"]).stat().st_mtime:
        paginas = json.loads(paginas_archivo.read_text())
        paginas_n = len(paginas)
        for i, p in enumerate(paginas, 1):
            assert p["caracteres"] > 500, f"Página vacía/casi vacía: {i}"
            assert f"Construcción y EDA · {i}" in p["texto"]
            assert (OUT / f"pagina-{i:02d}.png").stat().st_size > 10000
        pdf_estado = "actual_y_renderizado"
    trazabilidad = list(csv.DictReader((ROOT / "trazabilidad_metodos_clase.csv").open()))
    assert len(trazabilidad) == 38 and len({r["id"] for r in trazabilidad}) == 38
    resultado = {"estado": "OK", "paginas": paginas_n, "imagenes": manifiesto["figuras"], "tablas": 16,
                 "trazabilidad": 38, "manifiesto_sha256": "OK", "proporciones_imagenes": "OK",
                 "figuras_exigidas_incluidas": "OK",
                 "cobertura_y_minimos_word_vs_csv": "OK", "entradas_preservadas": "OK",
                  "pdf": pdf_estado, "instrucciones_script_secuencial": "OK"}
    (OUT / "validacion_informe_actual.json").write_text(json.dumps(resultado, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(resultado, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

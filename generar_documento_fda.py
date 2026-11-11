"""Genera el Word de la guía y cuatro láminas de árboles de decisiones.

Requiere Graphviz (dot) y pandoc. Ejecutar desde cualquier carpeta:
    python generar_documento_fda.py
La fuente editable del contenido es GUIA_RECONSTRUCCION_EDA.md.
"""

from pathlib import Path
import re
import subprocess
import tempfile
import json
import struct
from zipfile import ZipFile, ZIP_DEFLATED
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parent
W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
ET.register_namespace("w", W)


def diagramas():
    # Cada decisión tiene dos ramas reales. Las alternativas no ensayadas no se
    # presentan como perdedoras de una comparación que todavía no se ha hecho.
    arboles = [
        {
            "archivo":"diagrama_reconstruccion_funcional", "titulo":"RECONSTRUCCIÓN A · Datos y familia de curvas",
            "nodos":{
                "d0":("decision","¿Están verificadas las unidades\ny el calendario del producto?"),
                "p0":("pendiente","R01 · Procedencia pendiente\nSe encontró documentación CHIRPS,\npero no prueba identidad del TIFF.\nConversión actual: condicionada."),
                "d1":("decision","¿El año tiene observaciones\npara cubrir todo su dominio?"),
                "p1":("pendiente","R06 · No tratarlo como año completo\n2013: 24 de 72 archivos.\nSolicitar los 48 restantes; no\ninventar mayo–diciembre.\n1981–2012 sí: 72/año."),
                "d2":("decision","¿Se impone y justifica\nun cierre periódico?"),
                "a2":("alternativa","Fourier · Alternativa de clase\nPermite representación periódica.\nNo adoptada ni comparada aún;\nno se afirma que sea peor.\nC: ejercicio, línea 471."),
                "d3":("decision","¿Se requieren derivadas\ncontinuas de la curva?"),
                "a3":("alternativa","B-spline de orden 2\nLineal por tramos; pendiente\nconstante por tramos. No elegida\npara el análisis actual de derivadas.\nE: 354–364, comparación de órdenes."),
                "ok":("elegida","R02 · B-spline de orden 4\nAdoptada: curva con derivadas\ncontinuas y antecedente de lluvia.\nPhi se programa manualmente.\nE: 101–113; PDF EDA p. 43.")},
            "aristas":[("d0","d1","Sí"),("d0","p0","No"),("d1","d2","Sí"),("d1","p1","No"),("d2","a2","Sí"),("d2","d3","No"),("d3","ok","Sí"),("d3","a3","No")]
        },
        {
            "archivo":"arbol_reconstruccion_ajuste", "titulo":"RECONSTRUCCIÓN B · Ajuste, flexibilidad y huecos",
            "nodos":{
                "d0":("decision","¿Se adopta una base amplia\ny penalización de rugosidad?"),
                "alt":("alternativa","Penalización con lambda por VCG\nAlternativa de clase pendiente\nde contraste; no se ha probado\nque supere al ajuste actual.\nGrados de libertad: tr(S_lambda).\nNo exigir S²=S. C: 421–473."),
                "mco":("elegida","R03 · MCO con K variable\nBase actual: ajuste manual y\npropiedades de proyección claras.\nSe conserva, sin atribuirle\nsuperioridad sobre penalización.\nC: 83–105."),
                "d1":("decision","¿Elegir K solo porque\nbaja el ECM de ajuste?"),
                "noecm":("rechazada","No se elige ese criterio solo\nPremia el ajuste de los datos\nusados y puede sobreajustar.\nC: 304–305."),
                "vcg":("elegida","R04 · Mostrar ECM, VC y VCG\nElegir por VCG: K=6 / K=11.\nResultados descriptivos existentes;\nrevisar sensibilidad de derivadas.\nC: 215–305."),
                "d2":("decision","¿La curva tiene\nobservaciones faltantes?"),
                "gap":("elegida","R05–R08 · Solo filas observadas\nPhi_o y y_o; comprobar rango y S_o.\nNo poner ceros ni inventar registros.\nBanda sigma²I condicionada;\nguardar coeficientes e IDs.\nT: 129–143; C: 350–387."),
                "full":("elegida","R05–R08 · Usar toda la rejilla\nComprobar S²=S y tr(S)=K.\nResiduos y banda sigma²I puntual;\nno extender a bandas simultáneas.\nGuardar coeficientes e IDs.\nC: 95–105 y 350–387.")},
            "aristas":[("d0","alt","Sí"),("d0","mco","No"),("mco","d1","Continuar"),("d1","noecm","Sí"),("d1","vcg","No"),("vcg","d2","Continuar"),("d2","gap","Sí"),("d2","full","No")]
        },
        {
            "archivo":"diagrama_eda_funcional", "titulo":"EDA A · Resumen, dispersión y centralidad",
            "nodos":{
                "d0":("decision","¿Basta presentar\nsolo una media funcional?"),
                "solo":("rechazada","Media sola: insuficiente\nLa rúbrica también pide dispersión,\nFM, MBD y boxplot. No se omiten.\nWord: criterio 2."),
                "res":("elegida","E01 · Media + resúmenes robustos\nMedia por coeficientes y librería.\nMediana y recortada describen\ncentralidad sin depender solo\nde valores grandes. E: 155–173,\n445–489; PDF pp. 20, 60."),
                "d1":("decision","¿Queremos la DE exacta\nen la rejilla evaluada?"),
                "std":("alternativa","std.fd: versión suavizada\nAlternativa mostrada en clase;\nno confundirla con la raíz directa\nde varianza ni declarar igualdad.\nE: 187–192; PDF p. 24."),
                "de":("elegida","E02 · Raíz directa de varianza\nAdoptada para ambas fuentes.\nCV complementa: requiere media\ndistinta de cero; ±2DE no es IC.\nE: 180–196."),
                "d2":("decision","¿La pregunta de centralidad\nse basa en bandas de curvas?"),
                "mbd":("elegida","E05–E06 · MBD y boxplot\nAdoptados: región central 50 %,\nexpansión 1,5; empates inclusivos.\nAtípica no equivale a error.\nE: 493–518; PDF pp. 60–62."),
                "fm":("elegida","E05 · FM como complemento\nCentralidad punto a punto;\nmediana y recortada adoptadas.\nModal: otra alternativa de clase\nno aplicada; no sustituye MBD\nexigida. E: 449–489.")},
            "aristas":[("d0","solo","Sí"),("d0","res","No"),("res","d1","Continuar"),("d1","de","Sí"),("d1","std","No"),("de","d2","Continuar"),("d2","mbd","Sí"),("d2","fm","No")]
        },
        {
            "archivo":"arbol_eda_interpretacion", "titulo":"EDA B · Asociación, derivadas e interpretación",
            "nodos":{
                "d0":("decision","¿Hay sitios y años\nemparejados entre fuentes?"),
                "cross":("alternativa","Comparación funcional pareada\nDiferencias y correlación cruzada.\nMétodos de clase; fase posterior,\nno implementados sin metadatos.\nE: 313–323 y 406–438."),
                "cov":("elegida","E03 · Covarianza por fuente\nAdoptada: Phi Sigma_coef Phi'.\nDiagonal = varianza; correlación\nno demuestra causalidad ni sesgo.\nE: 253–289."),
                "d1":("decision","¿La curva derivada representa\nprecipitación acumulada?"),
                "ac":("alternativa","Derivada de acumulada\nSe interpreta como intensidad\nen la unidad temporal usada.\nEs el ejemplo de clase, pero no\nla escala actual en mm/día.\nE: 351–364."),
                "rate":("elegida","E04 · Derivada de tasa mm/día\nAdoptada: Dphi por coeficientes.\nResultado en mm/día²:\ncambio de tasa, no acumulado.\nE: 351–364; adaptación de escala."),
                "d2":("decision","¿El marcador queda en el borde\no no puede identificarse?"),
                "border":("pendiente","E04 · No dar fecha confiable\nFin_2=365 en 4.185/4.426 curvas.\nDeclarar limitación y estudiar\nregla y sensibilidad con clase.\nPendiente; no resuelto al dibujar."),
                "inside":("pendiente","E04/E07 · Marcador candidato\nEstar dentro del año no lo valida.\nJustificar regla, revisar estabilidad\ny alcance climatológico.\nNo inferir sesgo ENSO sin pares.\nResponsabilidad del equipo.")},
            "aristas":[("d0","cross","Sí"),("d0","cov","No"),("cov","d1","Continuar"),("d1","ac","Sí"),("d1","rate","No"),("rate","d2","Continuar"),("d2","border","Sí"),("d2","inside","No")]
        },
    ]
    colores={"decision":"#e8effa","elegida":"#e3f3eb","alternativa":"#eef0f4","pendiente":"#fff0d3","rechazada":"#fbe5e6"}
    for arbol in arboles:
        lineas=['digraph decisiones {',
            'graph [rankdir=TB, bgcolor="white", pad="0.18", nodesep="0.22", ranksep="0.30", fontname="Arial", fontsize=11, labelloc=t];',
            'node [fontname="Arial", fontsize=11, style="filled,rounded", color="#546577", margin="0.13,0.10"];',
            'edge [fontname="Arial", fontsize=10, arrowsize=0.65, color="#526678"];',
            'label='+json.dumps(arbol['titulo']+'\nVerde: adoptado · Gris: alternativa · Naranja: pendiente · Rojo: insuficiente\nC: Construcción FDA · E: Análisis Exploratorio · T: t-Funcional',ensure_ascii=False)+';']
        interiores={a for a,_,_ in arbol['aristas']}
        for ident,(tipo,texto) in arbol['nodos'].items():
            shape='diamond' if tipo=='decision' else 'box'
            grupo=', group="tronco"' if ident in interiores else ''
            lineas.append(f'{ident} [shape={shape}, fillcolor="{colores[tipo]}", label={json.dumps(texto,ensure_ascii=False)}{grupo}];')
        for a,b,etiqueta in arbol['aristas']:
            peso=100 if b in interiores else 1
            lineas.append(f'{a} -> {b} [label={json.dumps(etiqueta,ensure_ascii=False)}, weight={peso}];')
        for ident,(tipo,_) in arbol['nodos'].items():
            if tipo=='decision':
                assert sorted(e[2] for e in arbol['aristas'] if e[0]==ident)==['No','Sí'], ident
        lineas.append('}')
        dot=ROOT/(arbol['archivo']+'.dot');dot.write_text('\n'.join(lineas)+'\n')
        subprocess.run(['dot','-Tpng','-Gdpi=190',str(dot),'-o',str(dot.with_suffix('.png'))],check=True)
        subprocess.run(['dot','-Tsvg',str(dot),'-o',str(dot.with_suffix('.svg'))],check=True)
    (ROOT/'arboles_decision_fda.json').write_text(json.dumps(arboles,indent=2,ensure_ascii=False)+'\n')


def formatear_word(path):
    """Ajustes de legibilidad sobre el DOCX generado por pandoc."""
    with ZipFile(path) as z:
        partes = {n:z.read(n) for n in z.namelist()}
    styles = ET.fromstring(partes["word/styles.xml"])
    for style in styles.findall(f"{{{W}}}style"):
        nombre = style.get(f"{{{W}}}styleId", "")
        if nombre in ("Normal","BodyText","FirstParagraph","Compact"):
            rp = style.find(f"{{{W}}}rPr")
            if rp is None: rp=ET.SubElement(style,f"{{{W}}}rPr")
            fonts = rp.find(f"{{{W}}}rFonts")
            if fonts is None: fonts=ET.SubElement(rp,f"{{{W}}}rFonts")
            for key in ("ascii","hAnsi","cs"): fonts.set(f"{{{W}}}{key}","Calibri")
            sz=rp.find(f"{{{W}}}sz")
            if sz is None: sz=ET.SubElement(rp,f"{{{W}}}sz")
            sz.set(f"{{{W}}}val","22")
    partes["word/styles.xml"]=ET.tostring(styles,encoding="utf-8",xml_declaration=True)
    tree=ET.fromstring(partes["word/document.xml"])
    # Saltar antes del encabezado, no con párrafos vacíos que pueden generar
    # una página en blanco al terminar una figura de página completa.
    nuevos=("Árbol R-A", "Árbol R-B", "Árbol E-A", "Árbol E-B", "3. Reconstrucción:", "4. EDA funcional:", "8. Preguntas para el profesor")
    for par in tree.findall(f".//{{{W}}}p"):
        titulo=''.join(n.text or '' for n in par.iter(f"{{{W}}}t"))
        if titulo.startswith(nuevos):
            pr=par.find(f"{{{W}}}pPr")
            if pr is None:pr=ET.Element(f"{{{W}}}pPr");par.insert(0,pr)
            ET.SubElement(pr,f"{{{W}}}pageBreakBefore")
    for fila in tree.findall(f".//{{{W}}}tr"):
        pr=fila.find(f"{{{W}}}trPr")
        if pr is None: pr=ET.SubElement(fila,f"{{{W}}}trPr")
        ET.SubElement(pr,f"{{{W}}}cantSplit")
    for sec in tree.findall(f".//{{{W}}}sectPr"):
        tam=sec.find(f"{{{W}}}pgSz")
        if tam is not None:
            tam.set(f"{{{W}}}w","11906"); tam.set(f"{{{W}}}h","16838")
        marg=sec.find(f"{{{W}}}pgMar")
        if marg is not None:
            for key in ("top","bottom","left","right"): marg.set(f"{{{W}}}{key}","1080")
    partes["word/document.xml"]=ET.tostring(tree,encoding="utf-8",xml_declaration=True)
    with ZipFile(path,"w",compression=ZIP_DEFLATED) as z:
        for n,b in partes.items(): z.writestr(n,b)


def main():
    diagramas()
    texto=(ROOT/"GUIA_RECONSTRUCCION_EDA.md").read_text()
    # Conservar proporciones y encajar cada árbol en una página Word legible.
    for nombre in ('diagrama_reconstruccion_funcional','arbol_reconstruccion_ajuste','diagrama_eda_funcional','arbol_eda_interpretacion'):
        ancho,alto=struct.unpack('>II',(ROOT/(nombre+'.png')).read_bytes()[16:24])
        ancho_pulgadas=min(6.3,8.0*ancho/alto)
        texto=re.sub(r'(\]\('+re.escape(nombre)+r'\.png\))\{[^}]+\}',r'\g<1>{width='+f'{ancho_pulgadas:.3f}in'+'}',texto)
    texto=re.sub(r"```mermaid\n.*?\n```", "", texto, flags=re.S)
    texto=texto.replace("La secuencia completa también se conserva en formato Mermaid para la consulta técnica:", "")
    destino=ROOT/"GUIA_RECONSTRUCCION_EDA.docx"
    with tempfile.TemporaryDirectory(prefix="guia-fda-") as temp:
        fuente=Path(temp)/"guia.md"; fuente.write_text(texto)
        subprocess.run(["pandoc",str(fuente),"--from=markdown","--resource-path="+str(ROOT),"-M","lang=es-CO","-o",str(destino)],check=True)
    formatear_word(destino)
    with ZipFile(destino) as z:
        assert z.testzip() is None
        for n in z.namelist():
            if n.endswith(".xml"): ET.fromstring(z.read(n))
        imagenes=[n for n in z.namelist() if n.startswith("word/media/")]
        assert len(imagenes)==4, imagenes
    print(f"Word generado y XML/ZIP verificados: {destino}; diagramas incrustados: {len(imagenes)}")


if __name__=="__main__":
    main()

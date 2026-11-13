"""Comprueba las funciones reales de 03/04 desde los archivos preparados.

No ejecuta descargas, depuración ni modelos. Los resultados anteriores se leen;
solo genera verificacion_fda_clase.json y el diagrama PNG/SVG en la raíz.
Dependencias: numpy, pandas, matplotlib. Ejecución: python verificar_fda_clase.py
"""

import ast
import hashlib
import itertools
import json
from pathlib import Path
import platform
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ROOT = Path(__file__).resolve().parent


def funciones_notebook(nombre, nombres, contexto):
    """Carga solo definiciones seleccionadas: no ejecuta celdas ni sus salidas."""
    ruta = ROOT / nombre
    if not ruta.exists():
        ruta = ROOT / "archivo_historico_notebooks" / nombre
    nb = json.loads(ruta.read_text())
    encontrados = set()
    for i, celda in enumerate(nb["cells"], 1):
        if celda["cell_type"] != "code":
            continue
        arbol = ast.parse("".join(celda["source"]))
        for nodo in arbol.body:
            if isinstance(nodo, ast.FunctionDef) and nodo.name in nombres:
                exec(compile(ast.Module(body=[nodo], type_ignores=[]), f"{nombre}:celda{i}", "exec"), contexto)
                encontrados.add(nodo.name)
    if encontrados != set(nombres):
        raise ValueError(f"Definiciones no encontradas: {set(nombres) - encontrados}")
    return contexto


def diagrama():
    izquierda = [
        ("R01 · Dominio y unidades", "365 días; tasa en mm/día. Centros de péntada.\nPor qué: compartir dominio sin inventar mediciones.", "EDA R 47–61 / 127–133 · calendario: adaptación"),
        ("R02 · Base B-spline y Phi", "Orden 4 para curvas y derivadas continuas.\nPhi: valores de las bases en cada tiempo.", "Construcción R 73–80; EDA R 101–113"),
        ("R03 · Mínimos cuadrados", "Resolver (Phi'Phi)c = Phi'y.\nPor qué: estimar coeficientes paso a paso.", "Construcción R 83–105"),
        ("R04 · ECM, VC y VCG", "Elegir K por VCG, mostrando los tres criterios.\nPor qué: ECM de ajuste solo premia flexibilidad.", "Construcción R 215–305"),
        ("R05 · Verificar S", "Simetría, S²=S y tr(S)=K, sin penalización.\nPor qué: comprobar la proyección calculada.", "Construcción R 95–102"),
        ("R06 · Ajustar con observados", "Una máscara y un ajuste por curva incompleta.\nPor qué: ausencia no equivale a lluvia cero.", "t-Funcional R 129–143: observados, allí penalizados"),
        ("R07 · Residuos y bandas", "Curva ±2 EE bajo sigma² I.\nPor qué: incertidumbre puntual, no simultánea.", "Construcción R 350–387"),
        ("R08 · Guardar la construcción", "Coeficientes, IDs, años, dominio y unidades.\nPor qué: reconstruir y reutilizar en el EDA.", "Construcción R 83–105 · formatos: proyecto"),
    ]
    derecha = [
        ("E01 · Media funcional", "Promedio de coeficientes y contraste con librería.\nPor qué: describir y comprobar la curva media.", "EDA R 155–173 · PDF p. 20"),
        ("E02 · Varianza, DE y CV", "Dispersión absoluta y relativa a lo largo del año.\nMedia ±2 DE no es banda de confianza.", "EDA R 180–196 · PDF p. 24"),
        ("E03 · Covarianza y correlación", "Asociación entre épocas, para cada fuente.\nComprobar diagonal de covarianza = varianza.", "EDA R 253–289 · PDF p. 30"),
        ("E04 · Derivadas", "Cambios de tasa, picos y limitaciones de borde.\nFechas de temporada: regla adaptada, no receta.", "EDA R 351–364 · PDF p. 43"),
        ("E05 · Profundidad y centralidad", "FM y MBD; mediana y media recortada.\nPor qué: ordenar curvas del centro a la periferia.", "EDA R 445–489 / 493–518 · PDF pp. 53 / 60"),
        ("E06 · Boxplot funcional", "Región central 50 %, expansión 1,5 y atípicas.\nIdentificar estación–año; no borrar por defecto.", "EDA R 493–518 · PDF pp. 60–62"),
        ("E07 · Interpretación y límites", "Describir bimodalidad, dispersión y atípicas.\nSin pares geográficos no se demuestra sesgo.", "Word criterio 2 · ejemplo pareado: EDA R 60–61"),
        ("Registro y verificación", "Fuente de clase + decisión + código + resultado.\nListo solo para la comprobación efectuada.", "Guía, matriz CSV y verificación JSON del proyecto"),
    ]
    fig = plt.figure(figsize=(16, 18), facecolor="#f5f7fb")
    ax = fig.add_axes([0, 0, 1, 1]); ax.set(xlim=(0, 16), ylim=(0, 18)); ax.axis("off")
    ax.text(.65, 17.35, "DE REGISTROS A CURVAS, Y DE CURVAS A INTERPRETACIÓN", fontsize=17, weight="bold", color="#16334c")
    ax.text(.65, 16.90, "Reconstrucción y EDA funcional · métodos de clase y adaptaciones explícitas", fontsize=13, color="#415a70")
    for x, titulo, pasos, color in [(.65, "01 / RECONSTRUCCIÓN", izquierda, "#1c6572"), (8.25, "02 / EDA FUNCIONAL", derecha, "#6353a0")]:
        ax.text(x, 16.35, titulo, fontsize=13, weight="bold", color=color)
        for i, (nombre, motivo, fuente) in enumerate(pasos):
            y = 14.64 - i * 1.72
            ax.add_patch(FancyBboxPatch((x, y), 7.05, 1.46, boxstyle="round,pad=0.08", facecolor="white", edgecolor=color, linewidth=1.1))
            ax.text(x+.18, y+1.10, nombre, fontsize=12, weight="bold", color=color)
            ax.text(x+.18, y+.62, motivo, fontsize=10.2, va="center", linespacing=1.5, color="#253d50")
            ax.text(x+.18, y+.12, "Visto en: " + fuente, fontsize=8.5, color="#506678")
            if i < 7:
                ax.annotate("", xy=(x+3.52, y-.24), xytext=(x+3.52, y-.10), arrowprops={"arrowstyle": "->", "color": color})
    ax.text(.65, .62, "Fuentes en ../analitica/: Construcción R = Script_ConstrucciÃ³n_FDA.R; EDA R = Script_AnalisisExploratorio.R.", fontsize=9, color="#415a70")
    ax.text(.65, .29, "Las referencias numéricas son líneas de scripts y páginas físicas de PDF. Detalles, estados y límites: GUIA_RECONSTRUCCION_EDA.md.", fontsize=9, color="#415a70")
    for extension in ("png", "svg"):
        fig.savefig(ROOT / f"diagrama_reconstruccion_eda.{extension}", dpi=150, facecolor=fig.get_facecolor())
    plt.close(fig)


def verificar():
    spec = json.loads((ROOT / "fda_fase1/base_especificacion.json").read_text())
    base_contexto = {"np": np, "pd": pd, "ORDEN": spec["orden"], "RANGO": tuple(spec["rango"])}
    construccion = funciones_notebook("03_construccion_funcional.ipynb", ["nodos", "base_bspline", "ajuste_mco", "criterios", "banda"], dict(base_contexto))
    eda = funciones_notebook("04_eda_funcional.ipynb", ["nodos", "base_en_nodos", "base_bspline", "var_paso", "cov_paso", "deriv_base", "rangos", "profundidades", "boxplot_funcional"], dict(base_contexto))
    controles = []

    def control(nombre, error, tolerancia, detalle=""):
        error = float(error)
        controles.append({"control": nombre, "estado": "OK" if np.isfinite(error) and error <= tolerancia else "REVISAR", "error": error, "tolerancia": tolerancia, "detalle": detalle})

    archivos = [p for carpeta in ("fda_fase1", "fda_fase2") for p in (ROOT / carpeta).iterdir() if p.is_file()]
    antes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in archivos}
    curvas = pd.read_csv(ROOT / "fda_fase1/curvas_estaciones_mm.csv.gz", dtype={"STATION_ID": str})
    ce = pd.read_csv(ROOT / "fda_fase1/coef_estaciones.csv", dtype={"STATION_ID": str})
    e = np.load(ROOT / "fda_fase1/coeficientes_estaciones.npz", allow_pickle=False)
    s = np.load(ROOT / "fda_fase1/coeficientes_satelite.npz", allow_pickle=False)
    Y = curvas.iloc[:, 2:].to_numpy(float); C = ce.iloc[:, 2:].to_numpy(float)
    t = np.arange(1, 366, dtype=float)
    K = spec["K_estaciones"]; Ks = spec["K_satelite"]
    Phi = construccion["base_bspline"](t, tuple(spec["rango"]), K)
    control("R02: partición de la unidad", np.max(np.abs(Phi.sum(1)-1)), 1e-12)
    control("R02: misma base en 03 y 04", np.max(np.abs(Phi-eda["base_bspline"](t, K))), 1e-12)
    control("R08: claves estaciones CSV/NPZ", not (curvas.iloc[:, :2].equals(ce.iloc[:, :2]) and np.array_equal(ce.STATION_ID, e["station_id"].astype(str)) and np.array_equal(ce.YEAR, e["anio"])), 0)
    control("R08: coeficientes CSV/NPZ", np.max(np.abs(C-e["coef"])), 1e-10)
    error_c = 0.; max_traza = 0.; max_condicion = 0.
    for i, y in enumerate(Y):
        # Resolver por curva con las observaciones reales, sin alterar las matrices.
        m = np.isfinite(y); P = Phi[m]
        c = np.linalg.solve(P.T @ P, P.T @ y[m])
        error_c = max(error_c, float(np.max(np.abs(c-C[i]))))
        if not m.all():
            max_condicion = max(max_condicion, float(np.linalg.cond(P.T @ P)))
    control("R03/R06: reconstrucción de todas las curvas terrestres", error_c, 1e-9)
    for etiqueta, y in [("completa", Y[np.isfinite(Y).all(1)][0]), ("incompleta", Y[~np.isfinite(Y).all(1)][0])]:
        c, mapa, H, m = construccion["ajuste_mco"](y, Phi)
        control(f"R05: S simétrica {etiqueta}", np.max(np.abs(H-H.T)), 1e-10)
        control(f"R05: S idempotente {etiqueta}", np.max(np.abs(H@H-H)), 1e-10)
        control(f"R05: traza {etiqueta}", abs(np.trace(H)-K), 1e-10)
        ee = construccion["banda"](Phi, mapa, np.eye(m.sum()))
        control(f"R07: EE no negativo y finito {etiqueta}", not (np.isfinite(ee).all() and (ee >= 0).all()), 0)
    completa = Y[np.isfinite(Y).all(1)][0]
    c, mapa, H, _ = construccion["ajuste_mco"](completa, Phi)
    vc = np.mean(((completa-H@completa)/(1-np.diag(H)))**2)
    pred = [Phi[j]@np.linalg.lstsq(np.delete(Phi,j,0), np.delete(completa,j), rcond=None)[0] for j in range(365)]
    control("R04: VC atajo vs dejar uno fuera", abs(vc-np.mean((completa-pred)**2)), 1e-8)
    for fuente, k in [("estaciones", K), ("satelite", Ks)]:
        sel = pd.read_csv(ROOT / f"fda_fase1/seleccion_K_{fuente}.csv")
        control(f"R04: K {fuente} coincide con mínimo VCG guardado", abs(int(sel.loc[sel.VCG.idxmin(), "K"])-k), 0, "Comprueba tabla exportada; no repite toda la búsqueda de K.")
    sat = pd.read_csv(ROOT / "fda_fase1/curvas_satelite_mm.csv.gz")
    ts = np.array([float(c[1:]) for c in sat.columns[3:]])
    Ps = construccion["base_bspline"](ts, tuple(spec["rango"]), Ks)
    Cs = sat.iloc[:,3:].to_numpy(float) @ np.linalg.solve(Ps.T@Ps, Ps.T).T
    control("R08: claves satélite CSV/NPZ", not np.array_equal(sat.iloc[:,:3].to_numpy(), np.column_stack([s["anio"],s["fila"],s["col"]])), 0)
    control("R03: coeficientes satélite desde CSV redondeado", np.max(np.abs(Cs-s["coef"])), 1e-4, "CSV a 5 decimales y coeficientes NPZ float32; tolerancia por serialización.")
    Hs = Ps@np.linalg.solve(Ps.T@Ps,Ps.T)
    control("R05: S satélite idempotente", np.max(np.abs(Hs@Hs-Hs)), 1e-10)
    control("R05: traza satélite", abs(np.trace(Hs)-Ks), 1e-10)
    medias = pd.read_csv(ROOT / "fda_fase2/media_de_mediana.csv")
    tierra = pd.read_csv(ROOT / "fda_fase2/coef_satelite_tierra.csv.gz")
    indice = pd.MultiIndex.from_arrays([s["anio"],s["fila"],s["col"]])
    posiciones = indice.get_indexer(pd.MultiIndex.from_frame(tierra.iloc[:,:3]))
    if (posiciones < 0).any():
        raise ValueError("Curvas terrestres satelitales sin clave en el NPZ")
    for fuente, coef, nbase in [("est", C, K), ("sat", s["coef"][posiciones].astype(float), Ks)]:
        B = eda["base_bspline"](t, nbase); F = coef@B.T
        mu = B@coef.mean(0); v = eda["var_paso"](F, mu); cov = eda["cov_paso"](coef,B)
        control(f"E01: media por coeficientes vs curvas {fuente}", np.max(np.abs(mu-F.mean(0))), 1e-9)
        control(f"E01: media exportada {fuente}", np.max(np.abs(mu-medias[f"media_{fuente}"])), 1e-9)
        control(f"E02: DE exportada {fuente}", np.max(np.abs(np.sqrt(v)-medias[f"de_{fuente}"])), 1e-9)
        control(f"E03: diagonal covarianza vs varianza {fuente}", np.max(np.abs(np.diag(cov)-v)), 1e-9)
        D = eda["deriv_base"](t, nbase)
        control(f"E04: derivada de media vs media de derivadas {fuente}", np.max(np.abs(D@coef.mean(0)-(coef@D.T).mean(0))), 1e-9)
        if fuente == "est":
            Fe = F
    fm, mbd = eda["profundidades"](Fe)
    profundidad = pd.read_csv(ROOT / "fda_fase2/profundidad_y_fechas_estaciones.csv", dtype={"STATION_ID": str})
    control("E05: claves de profundidades", not profundidad.iloc[:,:2].equals(ce.iloc[:,:2]), 0)
    control("E05: FM de estaciones exportada", np.max(np.abs(fm.mean(1)-profundidad.FM)), 1e-10)
    control("E05: MBD de estaciones exportada", np.max(np.abs(mbd.mean(1)-profundidad.MBD)), 1e-10)
    bx = eda["boxplot_funcional"](Fe)
    at = pd.read_csv(ROOT / "fda_fase2/atipicas_todas_las_estaciones.csv", dtype={"STATION_ID": str})
    esperadas = set(map(tuple, at[["STATION_ID","YEAR"]].to_numpy()))
    obtenidas = set(map(tuple, ce.loc[bx["atip"], ["STATION_ID","YEAR"]].to_numpy()))
    control("E06: atípicas identificadas por estación y año", len(esperadas.symmetric_difference(obtenidas)), 0)
    # Contraste de la definición de clase, incluido un caso con empates.
    def mbd_por_pares(A):
        pares = list(itertools.combinations(range(len(A)),2))
        return np.array([np.mean([((np.minimum(A[j],A[k])<=a)&(a<=np.maximum(A[j],A[k]))).mean() for j,k in pares]) for a in A])
    pendientes = []
    for nombre, A in [("curvas reales", Fe[:20]), ("empates", np.array([[0.,0.,1.],[0.,0.,1.],[1.,0.,2.],[2.,0.,2.]]))]:
        diferencia = float(np.max(np.abs(eda["profundidades"](A)[1].mean(1)-mbd_por_pares(A))))
        control(f"E05: MBD vs definición por pares ({nombre})", diferencia, 1e-10)
    trazas = pd.read_csv(ROOT / "trazabilidad_metodos_clase.csv")
    control("Trazabilidad: 15 IDs únicos", abs(len(trazas.id.unique())-15), 0)
    for fuente in trazas.fuente_clase.unique():
        control(f"Trazabilidad: existe {fuente}", not (ROOT.parent/"analitica"/fuente).is_file(), 0)
    despues = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in archivos}
    control("Integridad: datos y resultados previos intactos", antes != despues, 0)
    pendientes.extend([
        {"id":"E04", "estado":"PENDIENTE", "motivo":"Validar regla adaptada de temporadas; fin_2 de borde no confiable", "fin_2_dia365":int((profundidad.fin_2==365).sum()), "fin_2_ausente":int(profundidad.fin_2.isna().sum())},
        {"id":"E07", "estado":"PENDIENTE", "motivo":"Representatividad de 1981–2013 y sesgo pareado no demostrados"},
        {"id":"R07", "estado":"PENDIENTE", "motivo":"Banda estacional heredada sin receta localizada en clase; banda sigma²I sí respaldada"},
    ])
    errores = sum(c["estado"]!="OK" for c in controles)
    resultado = {
        "fecha_utc":datetime.now(timezone.utc).isoformat(),
        "alcance":"Recálculo desde matrices preparadas; no ejecución completa desde crudos. No se repite el contraste con scikit-fda: sus salidas guardadas son evidencia previa.",
        "versiones":{"python":platform.python_version(),"numpy":np.__version__,"pandas":pd.__version__,"matplotlib":matplotlib.__version__},
        "estado":"COMPROBACIONES_OK_CON_PENDIENTES" if errores==0 else "REVISAR_COMPROBACIONES",
        "conteos":{"estaciones":int(ce.STATION_ID.nunique()),"curvas_estaciones":len(Y),"curvas_completas":int(np.isfinite(Y).all(1).sum()),"curvas_incompletas":int((~np.isfinite(Y).all(1)).sum()),"curvas_satelite":len(sat),"curvas_satelite_tierra":len(tierra),"K_estaciones":K,"K_satelite":Ks,"condicion_maxima_con_huecos":max_condicion,"puntos_negativos_estaciones":int((Fe<0).sum()),"minimo_estaciones_mm_dia":float(Fe.min())},
        "controles":controles,"pendientes":pendientes,"sha256_resultados_leidos":antes,
    }
    (ROOT/"verificacion_fda_clase.json").write_text(json.dumps(resultado,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    diagrama()
    print(json.dumps({"estado":resultado["estado"],"controles_OK":len(controles)-errores,"controles_a_revisar":errores,"pendientes":pendientes,"conteos":resultado["conteos"]},indent=2,ensure_ascii=False))
    return int(errores>0)


if __name__ == "__main__":
    raise SystemExit(verificar())

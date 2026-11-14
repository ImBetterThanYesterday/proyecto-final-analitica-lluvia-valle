"""Audita que las salidas EDA/FDA de Jacobo estén listas para la KAN."""

from pathlib import Path
from importlib.util import find_spec
import json

import pandas as pd


KAN_ROOT = Path(__file__).resolve().parent
ROOT = KAN_ROOT.parent
E0 = ROOT / "eda_fase0"
F1 = ROOT / "fda_fase1"
RESULTADOS = KAN_ROOT / "resultados"
RESULTADOS.mkdir(exist_ok=True)

embudo = pd.read_csv(E0 / "embudo_depuracion_preliminar.csv")
pares = pd.read_csv(
    E0 / "pares_estacion_anio_validos_preliminar.csv",
    dtype={"STATION_ID": "string"},
)
cobertura = pd.read_csv(
    E0 / "cobertura_estacion_anio_1981_2003_tras_depuracion.csv",
    dtype={"STATION_ID": "string"},
)
curvas_tierra = pd.read_csv(
    F1 / "curvas_estaciones_mm.csv.gz",
    usecols=["STATION_ID", "YEAR"],
    dtype={"STATION_ID": "string"},
)
curvas_satelite = pd.read_csv(
    F1 / "curvas_satelite_mm.csv.gz",
    usecols=["anio", "fila", "col"],
)
pares_kan = pd.read_csv(
    RESULTADOS / "identificadores_pares_kan_docente.csv",
    dtype={"station_id": "string"},
)

pares["YEAR"] = pares["YEAR"].astype(int)
pares["particion"] = pd.cut(
    pares["YEAR"], bins=[1980, 1995, 2000, 2003],
    labels=["entrenamiento", "validacion", "prueba"],
)

assert not pares.duplicated(["STATION_ID", "YEAR"]).any()
assert set(pares["YEAR"]) == set(range(1981, 2004))
assert not pares["particion"].isna().any()

control_cobertura = pares.merge(
    cobertura[["STATION_ID", "YEAR", "cob", "hueco", "meses20", "ok"]],
    on=["STATION_ID", "YEAR"], how="left", validate="one_to_one",
)
assert control_cobertura["ok"].fillna(False).all()
assert control_cobertura["cob"].ge(.90).all()
assert control_cobertura["hueco"].le(30).all()
assert control_cobertura["meses20"].eq(12).all()

anios_train = (
    pares.loc[pares["particion"].eq("entrenamiento")]
    .groupby("STATION_ID")["YEAR"].nunique()
)
assert anios_train.min() >= 7

claves_finales = set(map(tuple, pares[["STATION_ID", "YEAR"]].to_numpy()))
claves_curvas = set(map(tuple, curvas_tierra[["STATION_ID", "YEAR"]].to_numpy()))
assert claves_finales == claves_curvas

assert not curvas_satelite.duplicated(["anio", "fila", "col"]).any()
sat_por_anio = curvas_satelite.groupby("anio").size()
assert set(sat_por_anio.index) == set(range(1981, 2004))
assert sat_por_anio.nunique() == 1

particiones = pares.groupby("particion", observed=True).size().astype(int).to_dict()
kan_particiones = pares_kan.groupby("particion").size().astype(int).to_dict()

reporte = {
    "eda": {
        "estaciones_finales": int(pares.STATION_ID.nunique()),
        "pares_finales": int(len(pares)),
        "pares_por_particion": {str(k): int(v) for k, v in particiones.items()},
        "minimo_anios_entrenamiento_observado": int(anios_train.min()),
        "reglas_continuidad_verificadas": {
            "cobertura_minima": float(control_cobertura.cob.min()),
            "hueco_maximo": int(control_cobertura.hueco.max()),
            "meses_con_20_dias_minimo": int(control_cobertura.meses20.min()),
        },
        "embudo": embudo.to_dict(orient="records"),
    },
    "fda": {
        "curvas_terrestres": int(len(curvas_tierra)),
        "curvas_satelitales": int(len(curvas_satelite)),
        "pixeles_satelitales": int(curvas_satelite[["fila", "col"]].drop_duplicates().shape[0]),
        "anios_satelitales": int(curvas_satelite.anio.nunique()),
        "pares_eda_iguales_a_curvas_terrestres": True,
    },
    "emparejamiento_kan": {
        "estaciones": int(pares_kan.station_id.nunique()),
        "pares": int(len(pares_kan)),
        "pares_por_particion": {str(k): int(v) for k, v in kan_particiones.items()},
        "pixeles": int(pares_kan[["pixel_row", "pixel_col"]].drop_duplicates().shape[0]),
    },
    "dependencias_faltantes_entorno_actual": [
        modulo for modulo in ["rasterio", "skfda"] if find_spec(modulo) is None
    ],
}

(RESULTADOS / "auditoria_entradas_jacobo.json").write_text(
    json.dumps(reporte, ensure_ascii=False, indent=2), encoding="utf-8"
)
print(json.dumps(reporte, ensure_ascii=True, indent=2))

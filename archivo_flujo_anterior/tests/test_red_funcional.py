"""Pruebas de la fase 4 (05_red_funcional.py).

Componentes: base KAN, rugosidad, FPCA y paridad de parámetros (sintéticas, rápidas).
Resultados: equivalencia lineal, ausencia de fuga, no negatividad y formato de entrega
(requieren haber ejecutado `python 05_red_funcional.py`).
Uso:  python -m unittest tests/test_red_funcional.py -v
"""

import importlib.util
import json
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import torch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("rf", ROOT / "05_red_funcional.py")
rf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rf)
F4 = ROOT / "fda_fase4"


class Componentes(unittest.TestCase):
    def test_bases_kan_particion_de_la_unidad_en_el_rango(self):
        capa = rf.KANLayer(3, 2, G=5, rango=(-3.0, 3.0))
        x = torch.linspace(-2.999, 2.999, 400)[:, None].repeat(1, 3)
        B = capa.bases(x)
        self.assertLess(float((B.sum(-1) - 1).abs().max()), 1e-12)
        self.assertGreaterEqual(float(B.min()), 0.0)

    def test_activacion_kan_tiene_la_forma_del_enunciado(self):
        torch.manual_seed(0)
        capa = rf.KANLayer(2, 1, G=4)
        x = torch.randn(50, 2)
        spl = torch.einsum("nig,iog->nio", capa.bases(x), capa.coef)
        esperado = (capa.w_b * torch.nn.functional.silu(x).unsqueeze(-1) + capa.w_s * spl).sum(1)
        self.assertLess(float((capa(x) - esperado).abs().max()), 1e-14)

    def test_rugosidad_anula_funciones_lineales_y_es_semidefinida(self):
        R = rf.rugosidad(11)
        Phi = rf.base_bspline(rf.DIAS, (1, 365), 11)
        a = np.linalg.lstsq(Phi, 2.0 + 0.01 * rf.DIAS, rcond=None)[0]     # recta exacta en la base
        self.assertLess(abs(a @ R @ a), 1e-12)
        ev = np.linalg.eigvalsh(R)
        self.assertGreater(ev.min(), -1e-12)
        self.assertEqual(int((ev < 1e-10).sum()), 2)

    def test_gram_coincide_con_la_fase_1(self):
        W = np.load(ROOT / "fda_fase1" / "matrices_gram.npz")
        self.assertLess(np.abs(rf.gram(8) - W["W_est"]).max(), 1e-10)
        self.assertLess(np.abs(rf.gram(11) - W["W_sat"]).max(), 1e-10)

    def test_fpca_ortonormal_en_L2_y_ordenada(self):
        rng = np.random.default_rng(1)
        W = rf.gram(8)
        C = rng.normal(size=(300, 8)) @ rng.normal(size=(8, 8))
        mu, D, lam = rf.fpca(C, W)
        self.assertTrue(np.allclose(D.T @ W @ D, np.eye(8), atol=1e-10))
        self.assertTrue(np.all(np.diff(lam) <= 1e-12))
        # la varianza total en L² se conserva: traza(W Σ) = Σ λ
        S = np.cov(C - mu, rowvar=False, bias=True)
        self.assertAlmostEqual(np.trace(W @ S), lam.sum(), places=8)


class Resultados(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not (F4 / "resumen_fase4.json").exists():
            raise unittest.SkipTest("Ejecute primero python 05_red_funcional.py")
        cls.res = json.loads((F4 / "resumen_fase4.json").read_text(encoding="utf-8"))

    def test_kan_sintetica_rmse(self):
        self.assertLessEqual(self.res["KAN_sintetica"]["RMSE_prueba"], 1e-2)

    def test_equivalencia_lineal_con_regresion_funcional(self):
        eq = pd.read_csv(F4 / "tablas" / "equivalencia_lineal_regresion_funcional.csv")
        self.assertEqual(len(eq), 3)
        self.assertTrue((eq.max_abs_diferencia_sqrt_mm <= 1e-4).all())

    def test_parametros_kan_mlp_dentro_de_10_por_ciento(self):
        p = pd.read_csv(F4 / "tablas" / "parametros_cabezas.csv").set_index("cabeza")
        self.assertLessEqual(abs(p.loc["kan", "param_total"] - p.loc["mlp", "param_total"]) / p.loc["kan", "param_total"], .10)

    def test_tres_semillas_por_cabeza(self):
        e = pd.read_csv(F4 / "tablas" / "entrenamiento_por_semilla.csv")
        self.assertEqual(e.groupby("cabeza").semilla.nunique().to_dict(), {"kan": 3, "lineal": 3, "mlp": 3})

    def test_sin_fuga_entre_particiones(self):
        m = pd.read_csv(F4 / "tablas" / "metricas_por_curva.csv", dtype={"STATION_ID": str})
        self.assertTrue((m[m.particion == "validacion"].YEAR.between(1996, 2000)).all())
        self.assertTrue((m[m.particion == "prueba"].YEAR.between(2001, 2003)).all())
        self.assertEqual(self.res["estandarizacion_altitud"]["fuente"], "solo entrenamiento")
        self.assertIn("entrenamiento", self.res["FPCA"]["criterio"])

    def test_reconstruccion_no_negativa_y_formato_de_entrega(self):
        archivos = list((F4 / "predicciones").glob("predicciones_prueba_*.csv"))
        self.assertEqual(len(archivos), 1)
        p = pd.read_csv(archivos[0])
        self.assertEqual(list(p.columns), ["station_id", "year", "day", "value_mm"])
        self.assertGreaterEqual(p.value_mm.min(), 0)
        self.assertTrue(p.day.between(1, 365).all() and p.year.between(2001, 2003).all())
        self.assertEqual(len(p) % 365, 0)

    def test_bootstrap_con_al_menos_1000_remuestreos_y_wilcoxon(self):
        self.assertGreaterEqual(rf.N_BOOT, 1000)
        c = pd.read_csv(F4 / "tablas" / "comparaciones_pareadas_bootstrap_wilcoxon.csv")
        self.assertTrue({"IC95_inf", "IC95_sup", "wilcoxon_p"} <= set(c.columns))
        self.assertTrue({"Niña", "Neutro", "Niño", "todas"} >= set(c.fase_ENSO))

    def test_ajuste_simbolico_reporta_r2(self):
        s = pd.read_csv(F4 / "tablas" / "ajuste_simbolico_resumen.csv")
        self.assertTrue(s.R2.between(-1, 1).all() and len(s) >= 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)

"""Regresión de MBD: definición inclusiva de bandas del material de clase."""

import itertools
import unittest

import numpy as np

from verificar_fda_clase import funciones_notebook


class ProfundidadBandasTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.funciones = funciones_notebook(
            "04_eda_funcional.ipynb", ["rangos", "profundidades"], {"np": np}
        )

    def test_curvas_identicas_pertenecen_a_todas_las_bandas(self):
        _, mbd = self.funciones["profundidades"](np.zeros((4, 5)))
        np.testing.assert_allclose(mbd, 1., atol=1e-14)

    def test_empates_contra_definicion_inclusiva_por_pares(self):
        Y = np.array([[0., 0., 1.], [0., 0., 1.], [1., 0., 2.], [2., 0., 2.]])
        pares = list(itertools.combinations(range(len(Y)), 2))
        esperada = np.array([
            np.mean([(np.minimum(Y[j], Y[k]) <= y) & (y <= np.maximum(Y[j], Y[k]))
                     for j, k in pares], axis=0) for y in Y
        ])
        np.testing.assert_allclose(self.funciones["profundidades"](Y)[1], esperada, atol=1e-14)

    def test_sin_empates_contra_definicion_por_pares(self):
        Y = np.random.default_rng(2026).normal(size=(12, 9))
        pares = list(itertools.combinations(range(len(Y)), 2))
        esperada = np.array([
            np.mean([(np.minimum(Y[j], Y[k]) <= y) & (y <= np.maximum(Y[j], Y[k]))
                     for j, k in pares], axis=0) for y in Y
        ])
        np.testing.assert_allclose(self.funciones["profundidades"](Y)[1], esperada, atol=1e-14)


if __name__ == "__main__":
    unittest.main()

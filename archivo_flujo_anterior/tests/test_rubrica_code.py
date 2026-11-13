"""Source-level regression checks for rubric-critical behavior.

These checks do not replace numerical verification; they guard that the production
constructor stays package-independent and that EDA conventions remain explicit.
"""

from pathlib import Path
import json
import unittest


ROOT = Path(__file__).resolve().parents[1]


def notebook_sources(name: str) -> str:
    notebook = json.loads((ROOT / name).read_text(encoding="utf-8"))
    return "\n".join(
        "".join(cell.get("source", []))
        for cell in notebook.get("cells", [])
    )


class RubricSourceTests(unittest.TestCase):
    def test_construction_entrypoint_has_no_package_based_calculation_or_workspace_wipe(self):
        source = (ROOT / "03_construccion_funcional.R").read_text(encoding="utf-8")
        self.assertNotIn("splines::", source)
        self.assertNotIn("lsfit(", source)
        self.assertNotIn("rm(list = ls())", source)
        self.assertNotIn("rm(list = ls())", (ROOT / "04_eda_funcional.R").read_text(encoding="utf-8"))

    def test_training_only_k_selection_is_applied_to_both_sources(self):
        source = (ROOT / "03_construccion_funcional.R").read_text(encoding="utf-8")
        self.assertIn("ent_est <- completas_est & anio_de(est$id) <= ENTREN_FIN", source)
        self.assertIn("ent_sat <- completas_sat & anio_de(sat$id) <= ENTREN_FIN", source)

    def test_cross_language_checks_fail_for_both_sources_and_gram(self):
        source = (ROOT / "tests/validar_salidas_cruzadas.R").read_text(encoding="utf-8")
        self.assertIn("coef_satelite_python.csv.gz", source)
        self.assertIn("stopifnot", source)
        self.assertIn("error_sat", source)
        self.assertIn("error_w_sat", source)

    def test_eda_uses_teacher_scaled_fraiman_muniz_depth(self):
        r_source = (ROOT / "04_eda_funcional.R").read_text(encoding="utf-8")
        py_source = notebook_sources("04_eda_funcional.ipynb")
        self.assertIn("2 * (1 - abs(.5 - r / N)) - 1", r_source)
        self.assertIn("fm_pdf_t = 1 - np.abs(.5 - r / N)", py_source)
        self.assertIn("fm_t = 2 * fm_pdf_t - 1", py_source)

    def test_library_mean_verification_is_required_for_both_sources(self):
        r_source = (ROOT / "04_eda_funcional.R").read_text(encoding="utf-8")
        py_source = notebook_sources("04_eda_funcional.ipynb")
        self.assertIn('if (requireNamespace("fda", quietly = TRUE))', r_source)
        self.assertIn("err_media_e", r_source)
        self.assertIn("err_media_s", r_source)
        self.assertIn("PENDIENTE: fda no está instalado", r_source)
        self.assertIn("assert err_fm < 1e-10", py_source)
        self.assertIn("assert err_mbd < 1e-10", py_source)

    def test_fm_trim_includes_quantile_ties_and_boxplot_ties_use_class_order(self):
        r_source = (ROOT / "04_eda_funcional.R").read_text(encoding="utf-8")
        py_source = notebook_sources("04_eda_funcional.ipynb")
        self.assertIn("pr$FM >= quantile(pr$FM, .25", r_source)
        self.assertIn("order(-prof, seq_along(prof))[seq_len(ceiling(nrow(Y) / 2))]", r_source)
        self.assertIn("y <= inf | y >= sup", r_source)
        self.assertIn("np.lexsort((np.arange(len(prof)), -prof))", py_source)
        self.assertIn("incluidas = prof >= corte", py_source)
        self.assertIn("Y <= inf", py_source)
        self.assertIn("Y >= sup", py_source)

    def test_unpaired_enso_results_are_not_called_satellite_bias(self):
        py_source = notebook_sources("04_eda_funcional.ipynb").lower()
        self.assertNotIn("sesgo satelital depende de la fase", py_source)
        self.assertNotIn("this justifies including the oni", py_source)
        self.assertIn("fechas candidatas", py_source)


if __name__ == "__main__":
    unittest.main(verbosity=2)

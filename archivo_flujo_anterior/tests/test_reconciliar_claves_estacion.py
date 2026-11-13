import unittest
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from reconciliar_claves_estacion import cargar_claves  # noqa: E402


class ReconcileStationKeysTests(unittest.TestCase):
    def test_local_main_and_jacobo_keys_match_the_audit_counts(self):
        main, jacobo = cargar_claves(ROOT, "main")
        self.assertEqual(len(main), 4426)
        self.assertEqual(len(jacobo), 4485)
        self.assertEqual(len(main & jacobo), 4360)
        self.assertEqual(len(main - jacobo), 66)
        self.assertEqual(len(jacobo - main), 125)

    def test_diff_rows_keep_both_station_and_year_as_the_key(self):
        main, jacobo = cargar_claves(ROOT, "main")
        only_main = sorted(main - jacobo)
        only_jacobo = sorted(jacobo - main)
        self.assertTrue(all(len(k) == 2 and k[0] and k[1] for k in only_main + only_jacobo))
        self.assertEqual(len(set(only_main + only_jacobo)), 191)

    def test_every_key_difference_is_explained_by_a_rule(self):
        from reconciliar_claves_estacion import reconcile
        records, summary = reconcile(ROOT, "main")
        solo_j = [r for r in records if r["pertenece"] == "solo_jacobo"]
        solo_m = [r for r in records if r["pertenece"] == "solo_main"]
        self.assertEqual(len(solo_j), 125)
        self.assertTrue(all(r["diagnostico"].startswith("P6_cambio_de_regla") for r in solo_j))
        self.assertTrue(all(r["anios_validos_entren"] >= 7 and r["anios_validos_1981_2003"] < 10 for r in solo_j))
        p6 = [r for r in solo_m if r["diagnostico"].startswith("P6")]
        p5 = [r for r in solo_m if r["diagnostico"].startswith("P5")]
        self.assertEqual((len(p5), len(p6)), (10, 56))
        self.assertTrue(all(r["anios_validos_entren"] < 7 for r in p6))
        self.assertFalse([r for r in records if "revisar" in r["diagnostico"]])


if __name__ == "__main__":
    unittest.main(verbosity=2)

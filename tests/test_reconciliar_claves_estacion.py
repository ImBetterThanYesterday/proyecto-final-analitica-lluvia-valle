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


if __name__ == "__main__":
    unittest.main(verbosity=2)

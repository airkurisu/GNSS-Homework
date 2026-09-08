"""Regression checks for the optional offline reference-analysis tools."""
import math
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"tools"))
from compare_precise import enu, gpstime, load_sinex, load_sp3, sp3_velocity


class ReferenceToolsTests(unittest.TestCase):
    def test_gps_origin_and_week(self):
        self.assertEqual(gpstime(["1980", "1", "6", "0", "0", "0"]), 0)
        self.assertEqual(gpstime(["2024", "6", "7", "0", "0", "0"]), 2317*604800+432000)

    def test_sinex_estimate_not_apriori(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"station.snx"
            path.write_text("+SOLUTION/APRIORI\n 1 STAX JFNG A 1 24:159:43200 m 2 999 0\n-SOLUTION/APRIORI\n"
                            "+SOLUTION/ESTIMATE\n 1 STAX JFNG A 1 24:159:43200 m 2 1.0 0.001\n"
                            " 2 STAY JFNG A 1 24:159:43200 m 2 2.0 0.002\n"
                            " 3 STAZ JFNG A 1 24:159:43200 m 2 3.0 0.003\n-SOLUTION/ESTIMATE\n")
            result = load_sinex(path, "JFNG")
            self.assertEqual(result["xyz_m"], [1, 2, 3])
            with self.assertRaises(ValueError):
                load_sinex(path, "ABCD")

    def test_sp3_units_and_time_system(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"orbit.sp3"
            text = "%c M cc GPS\n"
            for i in range(9):
                text += f"* 2024 6 7 0 {i*5} 0\nPG01{20000:14.6f}{10000:14.6f}{15000:14.6f}{1:14.6f}\n"
            path.write_text(text)
            positions, epochs, _ = load_sp3(path)
            self.assertEqual(positions[("G01", epochs[0])], [20000000, 10000000, 15000000])
            path.write_text(text.replace("GPS", "UTC"))
            with self.assertRaises(ValueError):
                load_sp3(path)

    def test_velocity_polynomial_and_no_extrapolation(self):
        records = {("G01", t): [2*t+t**3, -3*t, 7.0] for t in range(-4, 5)}
        v = sp3_velocity(records, "G01", 0, 1)
        for actual, expected in zip(v, [2, -3, 0]):
            self.assertAlmostEqual(actual, expected, places=12)
        self.assertIsNone(sp3_velocity(records, "G01", -4, 1))

    def test_enu_equator(self):
        self.assertEqual(enu([1, 2, 3], [6378137, 0, 0]), [2, 3, 1])


if __name__ == "__main__":
    unittest.main()

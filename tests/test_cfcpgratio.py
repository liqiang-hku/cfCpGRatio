from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SOURCE_DIR = PROJECT_ROOT / "src"
if not (SOURCE_DIR / "cfcpgratio.py").is_file():
    raise FileNotFoundError("Expected the main script at src/cfcpgratio.py")
sys.path.insert(0, str(SOURCE_DIR))

from cfcpgratio import calculate_cfcpgratio, load_marker_indices


class TestCfCpGRatio(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tempdir.cleanup)
        self.path = Path(self.tempdir.name)
        atlas = pd.DataFrame(
            [
                ["chr1", 1, 4, 1, 4, "CellA", "marker_a", "U", 0, 0, 0],
                ["chr1", 5, 8, 4, 6, "CellB", "marker_b", "U", 0, 0, 0],
                ["chr1", 9, 12, 6, 7, "CellC", "marker_c", "U", 0, 0, 0],
            ],
            columns=[
                "chr", "start", "end", "startCpG", "endCpG", "target",
                "name", "direction", "CellA", "CellB", "CellC",
            ],
        )
        self.atlas_path = self.path / "atlas.tsv"
        atlas.to_csv(self.atlas_path, sep="\t", index=False)
        self.lbeta_path = self.path / "sample.lbeta"
        np.asarray(
            [[2, 10], [4, 10], [7, 10], [5, 5], [1, 5], [0, 0]],
            dtype="<u2",
        ).tofile(self.lbeta_path)

    def test_counts_and_scaling(self):
        names, indices, counts = load_marker_indices(self.atlas_path)
        self.assertEqual(names, ["CellA", "CellB", "CellC"])
        result = calculate_cfcpgratio(self.lbeta_path, names, indices, counts)
        rows = result.set_index("CellType")
        a, b, c = rows.loc["CellA"], rows.loc["CellB"], rows.loc["CellC"]
        self.assertEqual(int(a.Total_CpG_observations), 30)
        self.assertEqual(int(a.Unmethylated_CpG_observations), 17)
        self.assertEqual(int(a.Methylated_CpG_observations), 13)
        self.assertAlmostEqual(a.cfCpGRatio_raw, 17 / 30)
        self.assertEqual(int(b.Total_CpG_observations), 10)
        self.assertEqual(int(b.Unmethylated_CpG_observations), 4)
        self.assertAlmostEqual(b.cfCpGRatio_raw, 4 / 10)
        self.assertAlmostEqual(a.cfCpGRatio_scaled, 17 / 29)
        self.assertAlmostEqual(b.cfCpGRatio_scaled, 12 / 29)
        self.assertTrue(np.isnan(c.cfCpGRatio_raw))
        self.assertEqual(c.cfCpGRatio_scaled, 0.0)
        self.assertAlmostEqual(float(rows.cfCpGRatio_scaled.sum()), 1.0)

    def test_reject_invalid_methylation_counts(self):
        np.asarray(
            [[12, 10], [4, 10], [7, 10], [5, 5], [1, 5], [0, 0]],
            dtype="<u2",
        ).tofile(self.lbeta_path)
        names, indices, counts = load_marker_indices(self.atlas_path)
        with self.assertRaisesRegex(ValueError, "Methylated counts exceed"):
            calculate_cfcpgratio(self.lbeta_path, names, indices, counts)

    def test_command_line(self):
        output = self.path / "result.tsv"
        completed = subprocess.run(
            [
                sys.executable,
                str(SOURCE_DIR / "cfcpgratio.py"),
                "--lbeta", str(self.lbeta_path),
                "--atlas", str(self.atlas_path),
                "--output", str(output),
            ],
            capture_output=True,
            text=True,
        )
        self.assertEqual(completed.returncode, 0, msg=completed.stderr)
        self.assertTrue(output.is_file())
        result = pd.read_csv(output, sep="\t")
        self.assertEqual(len(result), 3)
        self.assertAlmostEqual(float(result.cfCpGRatio_scaled.sum()), 1.0)


if __name__ == "__main__":
    unittest.main()

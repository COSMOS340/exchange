import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("pc", ROOT / "Monocle2_PC.py")
pc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pc)


def metadata(n=20):
    return pd.DataFrame({
        "cell_id": [f"test_{i}" for i in range(n)],
        "tissue_group": ["Primary_CRC"] * 12 + ["Normal_adjacent"] * (n - 12),
        "trajectory_eligible": True, "cnv_summary_eligible": True, "identity_review": False,
        "known_CD45_enriched": False, "CNV_relative_band": "Low", "common_gene_UMI": 4,
        "coarse_state": ["C0"] * 6 + ["C1"] * (n - 6), "specimen_id": [f"S{i}" for i in range(n)],
        "sample_id": [f"L{i}" for i in range(n)], "subject_key": [f"P{i}" for i in range(n)],
        "dataset": "synthetic_test_only"})


class HandoffTests(unittest.TestCase):
    def test_exact_scope(self):
        self.assertEqual(int(pc.crc_mask(metadata()).sum()), 12)

    def test_duplicates_rejected(self):
        m = metadata(); m.loc[1, "cell_id"] = m.loc[0, "cell_id"]
        with self.assertRaises(ValueError):
            pc.crc_mask(m)

    def test_wrong_CNV_label_rejected(self):
        m = metadata(); m.loc[0, "CNV_relative_band"] = "low"
        with self.assertRaises(ValueError):
            pc.crc_mask(m)

    def test_ineligible_not_silently_deleted(self):
        m = metadata(); m.loc[0, "trajectory_eligible"] = False
        with self.assertRaises(ValueError):
            pc.crc_mask(m)

    def test_root_does_not_depend_on_CNV(self):
        m = metadata().iloc[:12].copy(); h = np.arange(12 * 30).reshape(12, 30)
        a = pc.choose_roots(m, h)[:2]
        m["CNV_relative_band"] = "High"
        self.assertEqual(a, pc.choose_roots(m, h)[:2])

    def test_source_engine_transformation(self):
        source = (ROOT / "reference/Overall_Monocle2_generated.R").read_text()
        engine = pc.build_engine(source)
        self.assertNotIn("patch_blocked_dispersion", engine)
        self.assertNotIn("patch_sparse_size_factors", engine)
        self.assertIn("memory_patch();patch_root_fallback()", engine)
        self.assertIn("residualModelFormulaStr='~dataset'", engine)
        with self.assertRaises(ValueError):
            pc.build_engine("unknown source")

    def test_manifest_modified_file_and_traversal(self):
        with tempfile.TemporaryDirectory() as d:
            folder = Path(d); (folder / "x").write_text("abc")
            records = pc.manifest(folder); pc.verify(folder, records)
            (folder / "x").write_text("abd")
            with self.assertRaises(ValueError):
                pc.verify(folder, records)
            records[0]["path"] = "../outside"
            with self.assertRaises(ValueError):
                pc.verify(folder, records)

    def test_complete_preparation_synthetic_only(self):
        with tempfile.TemporaryDirectory() as d:
            folder = Path(d); source = folder / "data/Overall"; source.mkdir(parents=True)
            m = metadata(); m.to_csv(source / "metadata.csv.gz", index=False)
            pd.DataFrame({"gene": ["A", "B", "C", "D"], "ordering": True}).to_csv(source / "genes.csv", index=False)
            x = sp.csc_matrix(np.ones((4, 20), dtype=np.int32))
            with h5py.File(source / "counts.h5", "w") as h:
                for name, value in [("data", x.data), ("indices", x.indices), ("indptr", x.indptr), ("shape", x.shape)]:
                    h.create_dataset(name, data=value)
            for name in ["Harmony", "PCA"]:
                np.save(source / f"{name}.npy", np.arange(600).reshape(20, 30))
            (folder / "reference").mkdir()
            (folder / "reference/Overall_Monocle2_generated.R").write_text(
                (ROOT / "reference/Overall_Monocle2_generated.R").read_text())
            (folder / "config.json").write_text((ROOT / "config.json").read_text())
            pc.save_json(folder / "INPUT_MANIFEST.json", pc.manifest(folder / "data"))
            with patch.object(pc, "ROOT", folder), patch.object(pc, "EXPECTED_SOURCE", 20), patch.object(pc, "EXPECTED_CRC", 12):
                pc.prepare()
                self.assertEqual(json.loads((folder / "work/preparation_summary.json").read_text())["CRC_cells"], 12)
                with h5py.File(folder / "work/CRC/counts.h5") as h:
                    self.assertEqual(list(h["shape"][:]), [4, 12])
                with self.assertRaises(RuntimeError):
                    pc.prepare()


if __name__ == "__main__":
    unittest.main()

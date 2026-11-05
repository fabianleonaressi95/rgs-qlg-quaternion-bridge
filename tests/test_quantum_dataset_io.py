import tempfile
import unittest
from pathlib import Path

from experiments.quantum_models.dataset_io import resolve_dataset_path

class DatasetPathTests(unittest.TestCase):
    def test_prefers_explicit_environment_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            explicit = root / "input.csv"
            explicit.write_text("Kingdom\n", encoding="utf-8")
            result = resolve_dataset_path(root / "experiment" / "base.py", environ={"CODON_TRAIN_CSV": str(explicit)}, cwd=root)
            self.assertEqual(result, explicit.resolve())

    def test_finds_dataset_next_to_module(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            module_dir = root / "experiment"
            module_dir.mkdir()
            data = module_dir / "codon_train.csv"
            data.write_text("Kingdom\n", encoding="utf-8")
            result = resolve_dataset_path(module_dir / "base.py", environ={}, cwd=root)
            self.assertEqual(result, data.resolve())

    def test_missing_dataset_has_actionable_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(FileNotFoundError, "CODON_TRAIN_CSV"):
                resolve_dataset_path(Path(tmp) / "experiment" / "base.py", environ={}, cwd=tmp)

    def test_invalid_explicit_path_does_not_fall_back(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "codon_train.csv").write_text("Kingdom\n", encoding="utf-8")
            with self.assertRaisesRegex(FileNotFoundError, "points to a missing file"):
                resolve_dataset_path(root / "experiment" / "base.py", environ={"CODON_TRAIN_CSV": str(root / "missing.csv")}, cwd=root)

if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import hashlib
import io
import os
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest.mock import patch

from experiments.shared_io import atomic_write_bytes, extract_zip_safely


class SharedIOTests(unittest.TestCase):
    def test_atomic_write_replaces_complete_file(self):
        with tempfile.TemporaryDirectory() as folder:
            target = Path(folder) / "payload.bin"
            target.write_bytes(b"old")
            atomic_write_bytes(b"new-complete", target)
            self.assertEqual(target.read_bytes(), b"new-complete")
            self.assertEqual(sorted(p.name for p in Path(folder).iterdir()), ["payload.bin"])

    def test_zip_extracts_selected_file(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            archive = root / "input.zip"
            with zipfile.ZipFile(archive, "w") as z:
                z.writestr("nested/data.csv", "a,b\n1,2\n")
                z.writestr("ignored.bin", b"x")
            result = extract_zip_safely(archive, root / "out", expected_suffixes=(".csv",))
            self.assertEqual(len(result), 1)
            self.assertEqual(result[0].read_text(), "a,b\n1,2\n")

    def test_zip_traversal_is_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            archive = root / "bad.zip"
            with zipfile.ZipFile(archive, "w") as z:
                z.writestr("../escape.csv", "bad")
            with self.assertRaises(ValueError):
                extract_zip_safely(archive, root / "out", expected_suffixes=(".csv",))
            self.assertFalse((root / "escape.csv").exists())

    def test_zip_size_limit_is_enforced(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            archive = root / "big.zip"
            with zipfile.ZipFile(archive, "w") as z:
                z.writestr("large.csv", "12345")
            with self.assertRaises(ValueError):
                extract_zip_safely(archive, root / "out", max_total_size=4)


if __name__ == "__main__":
    unittest.main()

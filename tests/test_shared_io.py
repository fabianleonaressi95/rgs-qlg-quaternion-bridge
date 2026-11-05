from __future__ import annotations

import gzip
import hashlib
import io
import os
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest.mock import patch

from experiments.shared_io import atomic_write_bytes, download_to_path, extract_zip_safely


class SharedIOTests(unittest.TestCase):
    def test_download_checks_size_and_sha_before_publish(self):
        payload = b"validated payload"
        digest = hashlib.sha256(payload).hexdigest()
        with tempfile.TemporaryDirectory() as folder:
            target = Path(folder) / "download.bin"
            target.write_bytes(b"previous-good")
            with patch("experiments.shared_io.urllib.request.urlopen", return_value=io.BytesIO(payload)):
                result = download_to_path(
                    "https://example.invalid/data", target,
                    expected_size=len(payload), expected_sha256=digest,
                )
            self.assertEqual(result, digest)
            self.assertEqual(target.read_bytes(), payload)

    def test_download_mismatch_preserves_existing_file(self):
        payload = b"short response"
        with tempfile.TemporaryDirectory() as folder:
            target = Path(folder) / "download.bin"
            target.write_bytes(b"previous-good")
            with patch("experiments.shared_io.urllib.request.urlopen", return_value=io.BytesIO(payload)):
                with self.assertRaises(ValueError):
                    download_to_path("https://example.invalid/data", target, expected_size=999)
            self.assertEqual(target.read_bytes(), b"previous-good")
            self.assertEqual(sorted(p.name for p in Path(folder).iterdir()), ["download.bin"])

    def test_gzip_download_is_validated_before_publish(self):
        payload = gzip.compress(b"fasta-data")
        with tempfile.TemporaryDirectory() as folder:
            target = Path(folder) / "genome.gz"
            with patch("experiments.shared_io.urllib.request.urlopen", return_value=io.BytesIO(payload)):
                download_to_path("https://example.invalid/genome.gz", target)
            self.assertEqual(gzip.decompress(target.read_bytes()), b"fasta-data")

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

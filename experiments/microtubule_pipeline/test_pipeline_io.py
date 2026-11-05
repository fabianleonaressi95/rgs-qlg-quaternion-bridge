import gzip
import hashlib
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from io import BytesIO

from pipeline_io import download_to_path


class _Response(BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


class DownloadToPathTests(unittest.TestCase):
    def test_streams_and_publishes_atomically(self):
        payload = b"scientific-data"
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "dataset.bin"
            with patch("pipeline_io.urllib.request.urlopen", return_value=_Response(payload)) as opener:
                digest = download_to_path("https://example.invalid/data", target)
            self.assertEqual(target.read_bytes(), payload)
            self.assertEqual(digest, hashlib.sha256(payload).hexdigest())
            self.assertEqual(opener.call_args.kwargs["timeout"], 60.0)

    def test_checksum_mismatch_does_not_replace_existing_file(self):
        payload = b"new"
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "dataset.bin"
            target.write_bytes(b"old")
            with patch("pipeline_io.urllib.request.urlopen", return_value=_Response(payload)):
                with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                    download_to_path("https://example.invalid/data", target, expected_sha256="0" * 64)
            self.assertEqual(target.read_bytes(), b"old")

    def test_gzip_stream_is_validated(self):
        payload = gzip.compress(b"valid gzip payload")
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "dataset.gz"
            with patch("pipeline_io.urllib.request.urlopen", return_value=_Response(payload)):
                download_to_path("https://example.invalid/data.gz", target)
            self.assertEqual(gzip.decompress(target.read_bytes()), b"valid gzip payload")

    def test_invalid_timeout_rejected(self):
        with self.assertRaises(ValueError):
            download_to_path("https://example.invalid/data", "unused", timeout=0)


if __name__ == "__main__":
    unittest.main()

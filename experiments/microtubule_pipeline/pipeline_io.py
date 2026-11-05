"""Shared, side-effect-safe I/O helpers for the microtubule experiments.

Downloads are streamed to a temporary sibling file and atomically published only
after transfer and format checks complete. Callers may pin a SHA-256 digest when
one is available from an independently trusted source.
"""
from __future__ import annotations

import gzip
import hashlib
import os
import shutil
import tempfile
import urllib.request
from pathlib import Path
from typing import Union

PathLike = Union[str, os.PathLike[str]]


def download_to_path(
    url: str,
    destination: PathLike,
    *,
    timeout: float = 60.0,
    expected_sha256: str | None = None,
    chunk_size: int = 1024 * 1024,
) -> str:
    """Download *url* atomically and return its SHA-256 hex digest.

    If expected_sha256 is supplied, a mismatch raises ValueError and the
    destination is left untouched. Gzip files are fully read before publish,
    which validates their stream CRC/trailer. A missing checksum is not
    represented as cryptographic verification.
    """
    if timeout <= 0:
        raise ValueError("timeout must be positive")
    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive")
    if expected_sha256 is not None:
        expected_sha256 = expected_sha256.strip().lower()
        if len(expected_sha256) != 64 or any(c not in "0123456789abcdef" for c in expected_sha256):
            raise ValueError("expected_sha256 must be a 64-character hexadecimal digest")

    target = Path(destination).expanduser()
    target.parent.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256()
    fd, temporary_name = tempfile.mkstemp(
        prefix=f".{target.name}.", suffix=".part", dir=str(target.parent)
    )
    os.close(fd)
    temporary = Path(temporary_name)

    try:
        request = urllib.request.Request(
            url, headers={"User-Agent": "rgs-qlg-quaternion-bridge/1.0"}
        )
        with urllib.request.urlopen(request, timeout=timeout) as response:
            with temporary.open("wb") as output:
                while True:
                    block = response.read(chunk_size)
                    if not block:
                        break
                    output.write(block)
                    digest.update(block)
                output.flush()
                os.fsync(output.fileno())

        actual_sha256 = digest.hexdigest()
        if expected_sha256 is not None and actual_sha256 != expected_sha256:
            raise ValueError(
                f"SHA-256 mismatch for {url}: expected {expected_sha256}, "
                f"received {actual_sha256}"
            )

        if target.suffix.lower() == ".gz":
            with gzip.open(temporary, "rb") as compressed:
                while compressed.read(chunk_size):
                    pass

        os.replace(temporary, target)
        return actual_sha256
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass

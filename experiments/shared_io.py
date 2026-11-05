"""Shared safe I/O primitives for repository experiments.

Downloads are streamed to sibling temporary files and published atomically.
Checksums/sizes are enforced before publication; ZIP extraction rejects path
traversal and symlink entries.
"""
from __future__ import annotations

import gzip
import hashlib
import os
import shutil
import stat
import tempfile
import urllib.request
import zipfile
from pathlib import Path
from typing import Iterable


def atomic_write_bytes(data: bytes, destination: str | os.PathLike[str]) -> None:
    target = Path(destination).expanduser()
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f".{target.name}.", suffix=".part", dir=target.parent)
    tmp = Path(name)
    try:
        with os.fdopen(fd, "wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(tmp, target)
    finally:
        tmp.unlink(missing_ok=True)


def download_to_path(
    url: str,
    destination: str | os.PathLike[str],
    *,
    timeout: float = 60.0,
    expected_sha256: str | None = None,
    expected_size: int | None = None,
    chunk_size: int = 1024 * 1024,
) -> str:
    """Stream a URL, validate optional size/SHA-256, then atomically publish."""
    if timeout <= 0 or chunk_size <= 0:
        raise ValueError("timeout and chunk_size must be positive")
    if expected_size is not None and expected_size < 0:
        raise ValueError("expected_size must be non-negative")
    if expected_sha256 is not None:
        expected_sha256 = expected_sha256.strip().lower()
        if len(expected_sha256) != 64 or any(c not in "0123456789abcdef" for c in expected_sha256):
            raise ValueError("expected_sha256 must be a 64-character hexadecimal digest")

    target = Path(destination).expanduser()
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f".{target.name}.", suffix=".part", dir=target.parent)
    os.close(fd)
    tmp = Path(name)
    digest = hashlib.sha256()
    size = 0
    try:
        request = urllib.request.Request(url, headers={"User-Agent": "rgs-qlg-quaternion-bridge/1.0"})
        with urllib.request.urlopen(request, timeout=timeout) as response, tmp.open("wb") as output:
            while True:
                block = response.read(chunk_size)
                if not block:
                    break
                output.write(block)
                digest.update(block)
                size += len(block)
            output.flush()
            os.fsync(output.fileno())

        actual = digest.hexdigest()
        if expected_size is not None and size != expected_size:
            raise ValueError(f"size mismatch for {url}: expected {expected_size}, received {size}")
        if expected_sha256 is not None and actual != expected_sha256:
            raise ValueError(f"SHA-256 mismatch for {url}: expected {expected_sha256}, received {actual}")
        if target.suffix.lower() == ".gz":
            with gzip.open(tmp, "rb") as stream:
                while stream.read(chunk_size):
                    pass
        os.replace(tmp, target)
        return actual
    finally:
        tmp.unlink(missing_ok=True)


def decompress_gzip_to_path(source: str | os.PathLike[str], destination: str | os.PathLike[str]) -> None:
    """Validate and decompress gzip data without exposing a partial output."""
    target = Path(destination).expanduser()
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f".{target.name}.", suffix=".part", dir=target.parent)
    tmp = Path(name)
    try:
        with os.fdopen(fd, "wb") as output, gzip.open(source, "rb") as compressed:
            shutil.copyfileobj(compressed, output, length=1024 * 1024)
            output.flush()
            os.fsync(output.fileno())
        os.replace(tmp, target)
    finally:
        tmp.unlink(missing_ok=True)


def extract_zip_safely(
    archive: str | os.PathLike[str],
    destination: str | os.PathLike[str],
    *,
    expected_suffixes: Iterable[str] = (".vcf", ".csv", ".txt", ".xlsx"),
    max_total_size: int = 5 * 1024**3,
) -> list[Path]:
    """Extract selected regular files; reject traversal, symlinks and oversized archives."""
    root = Path(destination).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    suffixes = tuple(s.lower() for s in expected_suffixes)
    extracted: list[Path] = []
    total = 0
    with zipfile.ZipFile(archive, "r") as bundle:
        for info in bundle.infolist():
            if info.is_dir() or not info.filename.lower().endswith(suffixes):
                continue
            mode = info.external_attr >> 16
            if stat.S_ISLNK(mode):
                raise ValueError(f"ZIP symlink entry is not allowed: {info.filename}")
            member = Path(info.filename)
            if member.is_absolute() or any(part in ("..", "") for part in member.parts):
                raise ValueError(f"unsafe ZIP member path: {info.filename}")
            target = (root / member).resolve()
            if target != root and root not in target.parents:
                raise ValueError(f"ZIP member escapes destination: {info.filename}")
            total += info.file_size
            if total > max_total_size:
                raise ValueError("ZIP selected content exceeds max_total_size")
            target.parent.mkdir(parents=True, exist_ok=True)
            fd, name = tempfile.mkstemp(prefix=f".{target.name}.", suffix=".part", dir=target.parent)
            tmp = Path(name)
            try:
                with os.fdopen(fd, "wb") as output, bundle.open(info, "r") as source:
                    shutil.copyfileobj(source, output, length=1024 * 1024)
                    output.flush()
                    os.fsync(output.fileno())
                if tmp.stat().st_size != info.file_size:
                    raise ValueError(f"truncated ZIP member: {info.filename}")
                os.replace(tmp, target)
            finally:
                tmp.unlink(missing_ok=True)
            extracted.append(target)
    return extracted

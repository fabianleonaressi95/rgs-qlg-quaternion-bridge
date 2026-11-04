#!/usr/bin/env python3
"""Build a reproducible metadata/provenance manifest for repository .npy files.

The manifest describes file/container metadata only. It does not certify the
scientific meaning, units, coordinate convention, or correctness of an array.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import math
import re
import struct
from pathlib import Path
from typing import Any

import numpy as np


def _read_npy_header(path: Path) -> tuple[tuple[int, ...], Any, bool, str, int]:
    with path.open("rb") as stream:
        magic = stream.read(6)
        if magic != b"\x93NUMPY":
            raise ValueError("invalid .npy magic")
        version_raw = stream.read(2)
        if len(version_raw) != 2:
            raise ValueError("truncated .npy version")
        version = (version_raw[0], version_raw[1])
        if version == (1, 0):
            length_raw = stream.read(2)
            if len(length_raw) != 2:
                raise ValueError("truncated v1 header length")
            header_len = struct.unpack("<H", length_raw)[0]
        elif version in ((2, 0), (3, 0)):
            length_raw = stream.read(4)
            if len(length_raw) != 4:
                raise ValueError("truncated v2/v3 header length")
            header_len = struct.unpack("<I", length_raw)[0]
        else:
            raise ValueError(f"unsupported .npy version {version}")
        raw_header = stream.read(header_len)
        if len(raw_header) != header_len:
            raise ValueError("truncated .npy header")
        encoding = "utf-8" if version == (3, 0) else "latin1"
        header = ast.literal_eval(raw_header.decode(encoding).strip())
        if not isinstance(header, dict):
            raise ValueError("header is not a dictionary")
        shape = tuple(int(n) for n in header["shape"])
        if any(n < 0 for n in shape):
            raise ValueError("negative dimension in shape")
        dtype = np.dtype(header["descr"])
        fortran_order = bool(header["fortran_order"])
        data_offset = 6 + 2 + len(length_raw) + header_len
        return shape, dtype, fortran_order, f"{version[0]}.{version[1]}", data_offset


def _git_blob_sha1(data: bytes) -> str:
    payload = b"blob " + str(len(data)).encode("ascii") + b"\0" + data
    return hashlib.sha1(payload).hexdigest()


def _python_mentions(root: Path, basename: str) -> list[str]:
    # These are leads for manual provenance review, not proof of generation.
    needle = basename.encode("utf-8")
    hits: list[str] = []
    for source in sorted(root.rglob("*.py")):
        if ".git" in source.parts:
            continue
        try:
            if needle in source.read_bytes():
                hits.append(source.relative_to(root).as_posix())
        except OSError:
            continue
    return hits


def build_manifest(root: Path) -> dict[str, Any]:
    data_root = root / "data"
    files = sorted(data_root.rglob("*.npy")) if data_root.exists() else []
    records: list[dict[str, Any]] = []
    for path in files:
        rel = path.relative_to(root).as_posix()
        raw = path.read_bytes()
        sha256 = hashlib.sha256(raw).hexdigest()
        try:
            shape, dtype, fortran_order, version, offset = _read_npy_header(path)
            elements = math.prod(shape) if shape else 1
            payload_bytes = len(raw) - offset
            expected_payload_bytes = elements * dtype.itemsize
            if dtype.hasobject:
                payload_check = "object_dtype_size_not_derivable"
            else:
                payload_check = (
                    "match" if payload_bytes == expected_payload_bytes else "mismatch"
                )
            header_status = "parsed"
            dtype_text = dtype.str if dtype.fields is None else repr(dtype.descr)
            shape_value: list[int] = list(shape)
        except Exception as exc:  # keep the inventory complete on malformed files
            shape_value = []
            dtype_text = None
            fortran_order = None
            version = None
            offset = None
            payload_bytes = None
            expected_payload_bytes = None
            payload_check = "not_checked"
            header_status = f"error: {type(exc).__name__}: {exc}"

        records.append({
            "path": rel,
            "size_bytes": len(raw),
            "sha256": sha256,
            "git_blob_sha1": _git_blob_sha1(raw),
            "npy_version": version,
            "shape": shape_value,
            "dtype": dtype_text,
            "fortran_order": fortran_order,
            "header_status": header_status,
            "data_offset_bytes": offset,
            "payload_bytes": payload_bytes,
            "expected_payload_bytes": expected_payload_bytes,
            "payload_size_check": payload_check,
            "units": None,
            "coordinate_convention": None,
            "semantic_role": "unassigned_pending_domain_review",
            "producer_candidates": _python_mentions(root, path.name),
            "provenance_status": "candidate_references_only_not_confirmed",
        })

    by_hash: dict[str, list[str]] = {}
    for record in records:
        by_hash.setdefault(record["sha256"], []).append(record["path"])
    for record in records:
        group = by_hash[record["sha256"]]
        record["duplicate_paths"] = group if len(group) > 1 else []

    return {
        "schema_version": 1,
        "generated_by": "tools/build_npy_manifest.py",
        "scientific_validation": "not_performed",
        "provenance_policy": (
            "Filename mentions are candidate references only; units, semantic role, "
            "coordinate convention and generating script require human confirmation."
        ),
        "array_count": len(records),
        "duplicate_content_groups": sum(1 for paths in by_hash.values() if len(paths) > 1),
        "arrays": records,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    root = args.repo_root.resolve()
    output = args.output or (root / "data" / "manifest" / "npy_manifest.json")
    manifest = build_manifest(root)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Wrote {manifest['array_count']} array records to {output}")
    print(f"Exact-content duplicate groups: {manifest['duplicate_content_groups']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

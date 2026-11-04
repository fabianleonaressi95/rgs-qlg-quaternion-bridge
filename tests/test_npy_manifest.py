from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np


def test_manifest_records_shape_dtype_hash_and_duplicates(tmp_path: Path) -> None:
    root = tmp_path
    (root / "data" / "sample").mkdir(parents=True)
    array = np.arange(12, dtype=np.float64).reshape(3, 4)
    first = root / "data" / "sample" / "A.npy"
    second = root / "data" / "sample" / "A-copy.npy"
    np.save(first, array)
    second.write_bytes(first.read_bytes())

    tool = Path(__file__).resolve().parents[1] / "tools" / "build_npy_manifest.py"
    output = root / "manifest.json"
    subprocess.run(
        [sys.executable, str(tool), "--repo-root", str(root), "--output", str(output)],
        check=True,
    )
    manifest = json.loads(output.read_text(encoding="utf-8"))
    assert manifest["array_count"] == 2
    assert manifest["duplicate_content_groups"] == 1
    by_path = {item["path"]: item for item in manifest["arrays"]}
    a = by_path["data/sample/A.npy"]
    assert a["shape"] == [3, 4]
    assert a["dtype"] == "<f8"
    assert len(a["sha256"]) == 64
    assert a["payload_size_check"] == "match"
    assert a["duplicate_paths"] == [
        "data/sample/A-copy.npy",
        "data/sample/A.npy",
    ]


def test_manifest_marks_semantics_as_unassigned(tmp_path: Path) -> None:
    (tmp_path / "data").mkdir()
    np.save(tmp_path / "data" / "unknown.npy", np.zeros((2,), dtype=np.int16))
    tool = Path(__file__).resolve().parents[1] / "tools" / "build_npy_manifest.py"
    output = tmp_path / "manifest.json"
    subprocess.run(
        [sys.executable, str(tool), "--repo-root", str(tmp_path), "--output", str(output)],
        check=True,
    )
    record = json.loads(output.read_text(encoding="utf-8"))["arrays"][0]
    assert record["units"] is None
    assert record["coordinate_convention"] is None
    assert record["provenance_status"] == "candidate_references_only_not_confirmed"

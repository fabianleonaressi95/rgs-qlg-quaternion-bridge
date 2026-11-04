import importlib.util
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CORE = ROOT / "experiments" / "rgs_psu_core"


def _load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_rgs_psu_core_sources_parse():
    import ast

    for path in CORE.glob("*.py"):
        if path.name in {
            "RGS_PSU_closure.py",
            "RGS_PSU_Equidistribution_Proof.py",
        }:
            ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


def test_600_cell_coordinate_count_and_unit_norms():
    module = _load_module(
        CORE / "RGS_PSU_Equidistribution_Proof.py",
        "rgs_psu_equidistribution_audit",
    )
    vertices = module.build_600_cell_vertices()
    assert vertices.shape == (120, 4)
    assert np.allclose(np.linalg.norm(vertices, axis=1), 1.0, atol=1e-12)
    assert np.unique(np.round(vertices, 12), axis=0).shape[0] == 120


def test_600_cell_finite_moments_match_s3_haar():
    module = _load_module(
        CORE / "RGS_PSU_Equidistribution_Proof.py",
        "rgs_psu_equidistribution_audit",
    )
    report = module.moment_audit()
    assert report["checks"]["center_of_mass_zero"]
    assert report["checks"]["second_moment_matches_S3_Haar"]
    assert report["checks"]["fourth_moment_matches_S3_Haar"]
    assert report["continuous_haar_equidistribution_claim"] == "NOT_ESTABLISHED"


def test_600_cell_quotient_claim_is_not_misidentified():
    module = _load_module(
        CORE / "RGS_PSU_Equidistribution_Proof.py",
        "rgs_psu_equidistribution_audit",
    )
    report = module.moment_audit()
    assert "one orbit point" in report["quotient_warning"]

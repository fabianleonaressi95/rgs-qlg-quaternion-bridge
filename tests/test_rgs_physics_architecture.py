import pytest

from src.rgs_physics_architecture import (
    ARCHITECTURE_SCHEMA,
    audit_core_conventions,
    audit_physics_architecture,
    physics_architecture_manifest,
)


def test_core_metric_clifford_and_q6_conventions_are_compatible():
    report = audit_core_conventions()
    assert report["passed"] is True
    assert report["metric_signature"] == [-1, 1, 1, 1]
    assert report["metric_match_residual"] < 1e-13
    assert report["gamma_match_max_residual"] < 1e-13
    assert report["gamma5_match_residual"] < 1e-13
    assert report["clifford_max_residual"] < 1e-12
    assert report["q6"]["vertices"] == 64
    assert report["q6"]["edges"] == 192
    assert report["q6"]["faces"] == 240
    assert report["q6"]["max_abs_B1B2"] == 0


def test_architecture_registry_keeps_distinct_physical_sectors_separate():
    manifest = physics_architecture_manifest()
    by_name = {item["sector"]: item for item in manifest["sectors"]}
    assert by_name["Dirac_Einstein_Cartan"]["status"] == "ACTION_DERIVED_LOCAL_OPERATORS"
    assert by_name["TEGR"]["status"] == "SEPARATE_GEOMETRIC_FORMULATION"
    assert by_name["SU3_QCD"]["status"] == "REPRESENTATION_AND_FLAVOUR_COMPONENTS_ONLY"
    assert by_name["MHD_fluids_superfluids_elasticity"]["status"] == "EFFECTIVE_LIMITS_NOT_YET_UNIFIED"
    assert manifest["schema"] == ARCHITECTURE_SCHEMA


def test_architecture_audit_does_not_claim_physical_unification():
    report = audit_physics_architecture()
    assert report["all_core_conventions_compatible"] is True
    assert report["unification_status"] == "CORE_CONVENTIONS_COMPATIBLE_BUT_PHYSICAL_COUPLINGS_OPEN"
    assert "spinor_covariant_derivative" in report["open_coupling_gates"]
    assert "Noether_and_Bianchi" in report["open_coupling_gates"]


def test_audit_rejects_invalid_tolerance():
    with pytest.raises(ValueError):
        audit_core_conventions(tolerance=0.0)

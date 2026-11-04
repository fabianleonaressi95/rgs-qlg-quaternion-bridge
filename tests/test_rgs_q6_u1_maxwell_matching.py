import math
import pytest

from src.rgs_q6_u1_maxwell_matching import q6_u1_maxwell_matching_report


@pytest.mark.parametrize(
    "a_s,a_t,g2,d",
    [
        (1.0, 0.05, 1.0, 6),
        (0.2, 0.01, 0.003, 6),
        (0.7, 0.03, 0.4, 3),
        (1.3, 0.2, 2.1, 4),
    ],
)
def test_wilson_projected_hamiltonian_matches_d5_maxwell_coefficients(a_s, a_t, g2, d):
    report = q6_u1_maxwell_matching_report(
        spatial_spacing=a_s,
        temporal_spacing=a_t,
        g_squared=g2,
        spatial_dim=d,
        hodge_eigenvalue=12.0,
    )
    assert report["all_checks_passed"], report["checks"]
    assert report["electric_hessian_transfer"] == pytest.approx(
        report["electric_hessian_D5"], rel=1e-12
    )
    assert report["magnetic_hessian_transfer"] == pytest.approx(
        report["magnetic_hessian_D5"], rel=1e-12
    )
    assert report["frequency_transfer"] == pytest.approx(
        math.sqrt(12.0) / a_s, rel=1e-12
    )


def test_matching_explicitly_does_not_claim_integer_lattice_embedding():
    report = q6_u1_maxwell_matching_report()
    assert any("integer sublattice" in item for item in report["not_proved"])
    assert "fixed-flux" in report["scope"]


@pytest.mark.parametrize("kwargs", [
    {"spatial_spacing": 0.0},
    {"temporal_spacing": -1.0},
    {"g_squared": 0.0},
    {"spatial_dim": 1},
])
def test_matching_rejects_invalid_parameters(kwargs):
    with pytest.raises(ValueError):
        q6_u1_maxwell_matching_report(**kwargs)

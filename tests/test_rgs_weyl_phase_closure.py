import numpy as np
import pytest

from src.rgs_weyl_phase_closure import (
    audit_original_q6_c6_edge_weights, audit_weyl_helicity,
    build_icosahedral_c6, helicity_projectors,
)

@pytest.mark.parametrize("momentum", [
    np.array([1.0, 0.0, 0.0]),
    np.array([0.3, -0.7, 1.2]),
    np.array([-2.0, 0.4, 0.1]),
])
def test_massless_weyl_helicity_projectors(momentum):
    report = audit_weyl_helicity(momentum)
    assert report["passed"]
    assert max(v for k, v in report.items() if k != "passed") < 1e-12

def test_helicity_is_undefined_at_zero_momentum():
    with pytest.raises(ValueError, match="zero spatial momentum"):
        helicity_projectors(np.zeros(3))

def test_icosahedral_c6_quadratic_identity_and_zero_diagonal():
    data = build_icosahedral_c6()
    c6 = data["C6"]
    assert c6.shape == (6, 6)
    assert data["C6_square_residual"] < 1e-12
    assert data["C6_diagonal_max"] == 0.0
    assert np.allclose(c6, c6.T)

def test_original_q6_weight_rule_is_uniform_not_c6_anisotropic():
    report = audit_original_q6_c6_edge_weights()
    assert report["q6_vertices"] == 64
    assert report["q6_edges"] == 192
    assert report["all_weights_equal_one"]
    assert report["weight_std"] == 0.0

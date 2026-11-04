import numpy as np
import pytest

from src.rgs_q6_integer_a5_d5_measure import (
    audit_q6_integer_a5_d5_measure,
    d5_bergman_volume_density,
    d5_type_iv_norm_polynomial,
)


def test_exact_rational_q6_mode_space_has_five_integer_flux_generators():
    report = audit_q6_integer_a5_d5_measure()
    assert report["rational_eigenspace_dimension"] == 5
    assert report["integer_mode_basis_shape"] == [192, 5]
    assert report["integer_mode_basis_rank"] == 5
    assert report["integer_gauss_failures"] == 0
    assert report["integer_hodge_failures"] == 0


def test_full_a5_action_preserves_integer_q6_mode_intersection():
    report = audit_q6_integer_a5_d5_measure()
    assert report["A5_group_order"] == 60
    assert report["integer_action_failures"] == 0
    assert report["A5_orthogonality_max"] < 1e-8


def test_a5_preserves_d5_type_iv_norm_and_bergman_density():
    report = audit_q6_integer_a5_d5_measure()
    assert report["A5_type_IV_norm_invariance_max"] < 1e-10
    z = np.array([.01+.02j, -.02+.01j, .015j, .005-.004j, .003+.006j])
    assert d5_bergman_volume_density(z) > 0
    assert d5_bergman_volume_density(np.zeros(5, dtype=complex)) == pytest.approx(1e5)
    assert d5_type_iv_norm_polynomial(np.zeros(5, dtype=complex)) == pytest.approx(1.0)


def test_d5_bergman_density_rejects_points_outside_domain():
    with pytest.raises(ValueError):
        d5_bergman_volume_density(np.array([1.0+0j, 0j, 0j, 0j, 0j]))

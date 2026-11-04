import numpy as np
import pytest

from src.rgs_ec_contact_tetrad_variation import (
    contact_lagrangian_density, contact_tetrad_variation,
    finite_difference_tetrad_gradient,
)


def test_tetrad_variation_matches_central_finite_difference():
    tetrad = np.array([
        [1.10, 0.08, 0.00, 0.00],
        [0.02, 0.93, 0.04, 0.00],
        [0.00, 0.03, 1.07, 0.02],
        [0.01, 0.00, 0.02, 0.98],
    ])
    current = np.array([0.4, -0.2, 0.3, 0.1])
    result = contact_tetrad_variation(tetrad, current, 1.7)
    numeric = finite_difference_tetrad_gradient(tetrad, current, 1.7)
    assert np.allclose(result["tetrad_action_gradient"], numeric,
                       rtol=2e-7, atol=2e-9)


def test_explicit_hilbert_tensor_is_metric_times_contact_lagrangian():
    tetrad = np.diag([1.2, 0.9, 1.1, 1.05])
    current = np.array([0.0, 0.3, -0.2, 0.1])
    result = contact_tetrad_variation(tetrad, current, 2.0)
    expected = result["contact_lagrangian"] * result["metric_covariant"]
    assert np.allclose(result["hilbert_stress_covariant_explicit"], expected)
    assert np.allclose(result["hilbert_stress_covariant_explicit"],
                       result["hilbert_stress_covariant_explicit"].T)
    assert result["full_on_shell_dirac_stress_tensor"] is None
    assert result["fluid_equation_of_state"] is None


def test_lorentz_scalar_is_frame_invariant_under_spatial_rotation():
    current = np.array([0.5, 0.2, -0.4, 0.1])
    angle = 0.37
    rotation = np.array([
        [1, 0, 0, 0],
        [0, np.cos(angle), -np.sin(angle), 0],
        [0, np.sin(angle), np.cos(angle), 0],
        [0, 0, 0, 1],
    ])
    rotated = rotation @ current
    assert contact_lagrangian_density(current, 1.3) == pytest.approx(
        contact_lagrangian_density(rotated, 1.3), abs=1e-14)


def test_rejects_singular_or_wrong_orientation_tetrads():
    current = np.array([0.1, 0.2, 0.3, 0.4])
    with pytest.raises(ValueError):
        contact_tetrad_variation(np.zeros((4, 4)), current, 1.0)
    with pytest.raises(ValueError):
        contact_tetrad_variation(-np.eye(4), current, 1.0)

import numpy as np
import pytest

from src.rgs_ec_contact_tetrad_variation import contact_tetrad_stress_energy


def test_fixed_state_tetrad_variation_gives_metric_times_contact_lagrangian():
    # M^{ab}=<J^a J^b>, isotropic rest-frame sample with timelike component.
    m = np.diag([2.0, 0.5, 0.5, 0.5])
    out = contact_tetrad_stress_energy(m, 16.0 / 3.0)
    assert out["mean_contact_lagrangian"] == pytest.approx(-0.5)
    assert np.allclose(out["stress_energy_tangent"], np.diag([0.5, -0.5, -0.5, -0.5]))


def test_spacelike_second_moment_has_explicit_vacuum_form_stress():
    m = np.diag([0.0, 1.0, 1.0, 1.0])
    out = contact_tetrad_stress_energy(m, 16.0 / 3.0)
    assert out["mean_contact_lagrangian"] == pytest.approx(3.0)
    assert np.allclose(out["stress_energy_tangent"], np.diag([-3.0, 3.0, 3.0, 3.0]))
    assert out["interaction_energy_density_explicit"] == pytest.approx(-3.0)
    assert out["interaction_pressure_explicit"] == pytest.approx(3.0)
    assert out["equation_of_state_explicit"] == -1.0
    assert out["state_response_included"] is False


def test_contact_tetrad_variation_rejects_invalid_inputs():
    with pytest.raises(ValueError):
        contact_tetrad_stress_energy(np.zeros((3, 3)), 1.0)
    with pytest.raises(ValueError):
        contact_tetrad_stress_energy(np.zeros((4, 4)), 0.0)
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


def test_spin_ensemble_isotropy_contract_exposes_time_space_covariance():
    from src.rgs_ec_spin_ensemble import spatial_isotropy_residual

    result = spatial_isotropy_residual(np.diag([0.2, 0.5, 0.5, 0.5]))
    assert result["spatial_isotropy_residual"] == pytest.approx(0.0)
    assert result["time_space_covariance_norm"] == pytest.approx(0.0)

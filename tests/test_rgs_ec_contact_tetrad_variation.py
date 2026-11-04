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

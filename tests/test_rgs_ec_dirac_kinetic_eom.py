import numpy as np
import pytest

from src.rgs_ec_algebraic_torsion import (
    dirac_axial_current, dirac_gamma_matrices_minus_plus,
)
from src.rgs_ec_dirac_kinetic_eom import (
    dirac_kinetic_euler_operator, effective_ec_dirac_eom_residual,
)
from src.rgs_ec_nonlinear_dirac_eom import effective_dirac_contact_sources


def test_kinetic_euler_operator_is_gamma_contraction_with_action_sign():
    derivative = np.arange(16, dtype=float).reshape(4, 4).astype(complex)
    gamma = dirac_gamma_matrices_minus_plus()
    expected_dirac = sum(gamma[a] @ derivative[:, a] for a in range(4))
    result = dirac_kinetic_euler_operator(derivative)
    assert np.allclose(result["dirac_operator"], expected_dirac)
    assert np.allclose(result["kinetic_euler_operator"], -expected_dirac)
    assert result["action_variation"].startswith("delta_S_D/delta_barpsi=")
    assert result["mass_term_included"] is False


def test_kinetic_operator_rejects_invalid_derivative_shape():
    with pytest.raises(ValueError):
        dirac_kinetic_euler_operator(np.zeros((4, 3)))


def test_full_local_ec_dirac_residual_adds_contact_source():
    psi = np.array([0.7 + 0.2j, -0.3 + 0.5j, 0.4 - 0.1j, 0.2 + 0.6j])
    current = dirac_axial_current(psi)
    derivative = np.zeros((4, 4), dtype=complex)
    derivative[:, 0] = np.array([0.1j, 0.2, -0.3j, 0.4])
    kinetic = dirac_kinetic_euler_operator(derivative)["kinetic_euler_operator"]
    source = effective_dirac_contact_sources(psi, current, 0.8)
    result = effective_ec_dirac_eom_residual(derivative, psi, current, 0.8)
    assert np.allclose(result["residual"], kinetic + source["dL_dbarpsi"])
    assert result["spin_connection_supplied_externally"] is True
    assert result["mass_term_included"] is False
    assert result["conservation_status"].startswith("NOT_ESTABLISHED")


def test_zero_derivative_and_zero_spinor_give_zero_local_residual():
    psi = np.zeros(4, dtype=complex)
    result = effective_ec_dirac_eom_residual(
        np.zeros((4, 4), dtype=complex), psi, np.zeros(4), 1.0
    )
    assert np.allclose(result["residual"], 0.0)
    assert result["on_shell"] is True


def test_kinetic_euler_operator_is_linear_in_covariant_derivative():
    rng = np.random.default_rng(20261104)
    left = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    right = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    a, b = 0.3 - 0.2j, -0.7 + 0.1j
    combined = dirac_kinetic_euler_operator(a * left + b * right)
    separate_left = dirac_kinetic_euler_operator(left)
    separate_right = dirac_kinetic_euler_operator(right)
    expected = (
        a * separate_left["kinetic_euler_operator"]
        + b * separate_right["kinetic_euler_operator"]
    )
    assert np.allclose(combined["kinetic_euler_operator"], expected)

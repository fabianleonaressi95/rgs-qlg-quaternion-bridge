import numpy as np
import pytest

from src.rgs_ec_algebraic_torsion import dirac_axial_current
from src.rgs_ec_full_dirac_stress import (
    axial_contact_hilbert_tensor,
    dirac_kinetic_hilbert_from_bilinears,
    dirac_kinetic_hilbert_from_spinor,
    full_dirac_hilbert_tensor,
    pure_spinor_second_moment,
)


def test_bilinear_kinetic_hilbert_tensor_is_symmetric_and_real():
    right = 1j * np.array([
        [1.0, 0.2, 0.0, 0.0],
        [0.2, 2.0, 0.1, 0.0],
        [0.0, 0.1, 3.0, 0.4],
        [0.0, 0.0, 0.4, 4.0],
    ])
    result = dirac_kinetic_hilbert_from_bilinears(right, np.zeros((4, 4)))
    assert np.allclose(result["tensor"], -0.5 * right.imag)
    assert result["symmetry_residual"] < 1e-14
    assert result["reality_residual"] == pytest.approx(0.0)


def test_bilinear_kinetic_rejects_nonreal_physical_tensor():
    with pytest.raises(ArithmeticError):
        dirac_kinetic_hilbert_from_bilinears(np.ones((4, 4)), np.zeros((4, 4)))


def test_spinor_builder_zero_derivative_has_zero_kinetic_tensor():
    psi = np.array([1.0 + 0.2j, 0.3j, -0.4, 0.7 + 0.1j])
    result = dirac_kinetic_hilbert_from_spinor(psi, np.zeros((4, 4)))
    assert np.allclose(result["tensor"], 0.0)
    assert result["reality_residual"] < 1e-14


def test_contact_tensor_matches_fixed_state_variation():
    moment = np.diag([2.0, 0.5, 0.5, 0.5])
    result = axial_contact_hilbert_tensor(moment, 16.0 / 3.0)
    assert result["contact_lagrangian"] == pytest.approx(-0.5)
    assert np.allclose(result["tensor"], np.diag([0.5, -0.5, -0.5, -0.5]))


def test_full_tensor_adds_kinetic_and_contact_without_claiming_closure():
    kinetic = np.diag([1.0, 2.0, 3.0, 4.0])
    moment = np.diag([2.0, 0.5, 0.5, 0.5])
    result = full_dirac_hilbert_tensor(kinetic, moment, 16.0 / 3.0)
    assert np.allclose(
        result["total_hilbert_tensor"],
        result["kinetic_hilbert_tensor"] + result["contact_hilbert_tensor"],
    )
    assert result["dirac_eom_enforced"] is False
    assert result["conservation_status"].startswith("REQUIRES_")


def test_pure_spinor_second_moment_is_outer_product_of_axial_current():
    psi = np.array([1.0 + 0.2j, 0.3j, -0.4, 0.7 + 0.1j])
    moment = pure_spinor_second_moment(psi)
    current = dirac_axial_current(psi)
    assert np.allclose(moment, np.outer(current, current))

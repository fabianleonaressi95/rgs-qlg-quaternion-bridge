import numpy as np
import pytest

from src.rgs_spin_connection import (
    dirac_covariant_derivative,
    levi_civita_spin_connection,
    orthonormal_dirac_derivative,
    spinor_connection_from_lorentz_connection,
)


def test_constant_minkowski_coframe_has_zero_spin_connection():
    result = levi_civita_spin_connection(np.eye(4), np.zeros((4, 4, 4)))
    assert result["passed"] is True
    assert np.allclose(result["metric"], np.diag([-1.0, 1.0, 1.0, 1.0]))
    assert np.linalg.norm(result["omega_mixed"]) == pytest.approx(0.0)
    assert result["torsion_free_residual"] < 1e-13


def test_spherical_minkowski_coframe_has_torsion_free_connection():
    radius = 2.3
    theta = 0.8
    coframe = np.diag([1.0, 1.0, radius, radius * np.sin(theta)])
    derivatives = np.zeros((4, 4, 4))
    # de[a, coordinate_index, derivative_coordinate]
    derivatives[2, 2, 1] = 1.0
    derivatives[3, 3, 1] = np.sin(theta)
    derivatives[3, 3, 2] = radius * np.cos(theta)

    result = levi_civita_spin_connection(coframe, derivatives)
    assert result["passed"] is True
    assert result["antisymmetry_residual"] < 1e-12
    assert result["tetrad_postulate_residual"] < 1e-12
    assert result["torsion_free_residual"] < 1e-12
    assert np.linalg.norm(result["omega_lower"]) > 0.1


def test_spinor_derivative_uses_spin_connection_and_frame_conversion():
    omega = np.zeros((4, 4, 4))
    omega[0, 1, 2] = 0.2
    omega[1, 0, 2] = -0.2
    spin_connection = spinor_connection_from_lorentz_connection(omega)
    psi = np.array([1.0, 0.2j, -0.3, 0.5 + 0.1j])
    partial = np.arange(16).reshape(4, 4).astype(complex) / 10.0
    coordinate = dirac_covariant_derivative(partial, psi, spin_connection)
    expected = np.column_stack([
        partial[:, mu] + spin_connection[mu] @ psi for mu in range(4)
    ])
    assert np.allclose(coordinate, expected)
    inverse = np.array([[1.0, 0.0, 0.0, 0.0],
                        [0.0, 0.5, 0.0, 0.0],
                        [0.0, 0.0, 2.0, 0.0],
                        [0.0, 0.0, 0.0, 1.5]])
    frame = orthonormal_dirac_derivative(coordinate, inverse)
    assert np.allclose(frame, coordinate @ inverse)


def test_spin_connection_rejects_singular_coframe():
    with pytest.raises(ValueError):
        levi_civita_spin_connection(np.zeros((4, 4)), np.zeros((4, 4, 4)))


def test_spinor_derivative_rejects_wrong_shapes():
    with pytest.raises(ValueError):
        dirac_covariant_derivative(np.zeros((4, 3)), np.ones(4), np.zeros((4, 4, 4)))


def test_frame_derivative_feeds_existing_ec_dirac_euler_operator():
    from src.rgs_ec_dirac_kinetic_eom import dirac_kinetic_euler_operator

    omega = np.zeros((4, 4, 4))
    omega[0, 1, 2] = 0.15
    omega[1, 0, 2] = -0.15
    spin_connection = spinor_connection_from_lorentz_connection(omega)
    psi = np.array([0.4 + 0.1j, -0.2j, 0.7, 0.3 - 0.2j])
    partial = np.arange(16, dtype=float).reshape(4, 4).astype(complex) / 20.0
    coordinate = dirac_covariant_derivative(partial, psi, spin_connection)
    inverse = np.diag([1.0, 0.8, 1.2, 0.9])
    frame = orthonormal_dirac_derivative(coordinate, inverse)

    result = dirac_kinetic_euler_operator(frame)
    expected = -sum(
        __import__("src.rgs_ec_algebraic_torsion", fromlist=["dirac_gamma_matrices_minus_plus"])
        .dirac_gamma_matrices_minus_plus()[a] @ frame[:, a]
        for a in range(4)
    )
    assert frame.shape == (4, 4)
    assert np.allclose(result["kinetic_euler_operator"], expected)
    assert result["mass_term_included"] is False

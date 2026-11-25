import numpy as np
import pytest

from src.pellis_su2_fibre import (
    EPSILON,
    conjugate_intertwining_residual,
    doubled_edge_covariance_residual,
    edge_covariance_residual,
    fibre_exchange_matrix,
    su2_matrix,
)


def test_epsilon_is_unitary_skew_and_squares_to_minus_identity():
    identity = np.eye(2, dtype=complex)
    assert np.array_equal(EPSILON.T, -EPSILON)
    assert np.allclose(EPSILON.conj().T @ EPSILON, identity, atol=1e-14)
    assert np.allclose(EPSILON @ EPSILON, -identity, atol=1e-14)


@pytest.mark.parametrize("theta,phase", [
    (0.0, 0.0), (0.2, 0.7), (1.1, 2.4), (np.pi / 2, -1.2),
])
def test_su2_equivariance(theta, phase):
    a = np.cos(theta) * np.exp(1j * phase)
    b = np.sin(theta) * np.exp(-0.3j * phase)
    u = su2_matrix(a, b)
    assert conjugate_intertwining_residual(u) < 1e-12
    assert np.allclose(EPSILON @ u @ np.linalg.inv(EPSILON),
                       np.conjugate(u), atol=1e-12)


def test_fibre_exchange_is_self_adjoint_unitary_involution():
    j = fibre_exchange_matrix()
    identity = np.eye(4, dtype=complex)
    assert np.allclose(j.conj().T, j, atol=1e-14)
    assert np.allclose(j @ j, identity, atol=1e-14)


def test_edge_covariance_on_random_fibre_data():
    rng = np.random.default_rng(20261125)
    for _ in range(25):
        z = rng.normal(size=2) + 1j * rng.normal(size=2)
        z /= np.linalg.norm(z)
        u = su2_matrix(z[0], z[1])
        source = rng.normal(size=2) + 1j * rng.normal(size=2)
        target = rng.normal(size=2) + 1j * rng.normal(size=2)
        assert edge_covariance_residual(u, source, target) < 1e-12


def test_doubled_edge_exchange_covariance():
    rng = np.random.default_rng(120)
    for _ in range(25):
        z = rng.normal(size=2) + 1j * rng.normal(size=2)
        z /= np.linalg.norm(z)
        u = su2_matrix(z[0], z[1])
        vectors = [rng.normal(size=2) + 1j * rng.normal(size=2) for _ in range(4)]
        assert doubled_edge_covariance_residual(u, *vectors) < 1e-12


def test_su2_parameter_validation():
    with pytest.raises(ValueError):
        su2_matrix(1.0, 1.0)

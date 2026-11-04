import numpy as np
import pytest

from src.rgs_q6_euclidean_time_wilson import build_q6_euclidean_time_data
from src.rgs_q6_u1_character_electric_flux import (
    audit_q6_u1_character_electric_flux,
    build_q6_integer_cycle_basis,
    q6_u1_bessel_coefficients,
    q6_u1_character_kernel,
    q6_u1_electric_weights,
    q6_u1_integer_gauss_law,
)


def test_character_expansion_reconstructs_wilson_boltzmann_factor():
    beta = 2.3
    theta = np.linspace(-np.pi, np.pi, 257)
    reconstructed = q6_u1_character_kernel(beta, theta, 18)
    assert np.allclose(reconstructed, np.exp(beta * np.cos(theta)), rtol=1e-11, atol=1e-12)


def test_bessel_coefficients_are_positive_and_symmetric_in_integer_order():
    coeff = q6_u1_bessel_coefficients(1.7, 12)
    assert coeff.shape == (13,)
    assert np.all(coeff > 0)
    assert q6_u1_electric_weights(1.7, 12) == pytest.approx(coeff)


def test_fundamental_integer_cycle_basis_has_exact_gauss_law():
    data = build_q6_euclidean_time_data(2)
    basis = build_q6_integer_cycle_basis(data)
    assert basis.shape == (129, 192)
    assert np.issubdtype(basis.dtype, np.integer)
    divergence = np.asarray(data["D0"], dtype=np.int64).T @ basis.T
    assert np.array_equal(divergence, np.zeros((64, 129), dtype=np.int64))
    assert np.linalg.matrix_rank(basis.astype(float)) == 129


def test_integer_gauss_law_rejects_noninteger_flux():
    with pytest.raises(ValueError, match="integer-valued"):
        q6_u1_integer_gauss_law(np.full(192, 0.5))


def test_integer_gauss_law_accepts_cycle_flux():
    basis = build_q6_integer_cycle_basis()
    assert np.array_equal(q6_u1_integer_gauss_law(3 * basis[7]), np.zeros(64, dtype=int))


def test_character_audit_passes():
    report = audit_q6_u1_character_electric_flux()
    assert report["all_checks_passed"], report["checks"]
    assert report["integer_cycle_rank"] == 129


@pytest.mark.parametrize("beta", [0.0, -1.0, np.nan])
def test_invalid_beta_rejected(beta):
    with pytest.raises(ValueError):
        q6_u1_bessel_coefficients(beta, 4)


@pytest.mark.parametrize("nmax", [-1, 1.5])
def test_invalid_fourier_cutoff_rejected(nmax):
    with pytest.raises(ValueError):
        q6_u1_bessel_coefficients(1.0, nmax)

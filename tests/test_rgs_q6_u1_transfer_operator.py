import numpy as np
import pytest

from src.rgs_q6_u1_character_electric_flux import build_q6_integer_cycle_basis
from src.rgs_q6_u1_transfer_operator import (
    audit_q6_u1_transfer_operator,
    q6_u1_anisotropic_wilson_couplings,
    q6_u1_electric_transfer_eigenvalues,
    q6_u1_electric_hamiltonian_eigenvalues,
    q6_u1_spatial_wilson_potential,
    q6_u1_symmetric_transfer_kernel,
)


def test_q6_integer_flux_transfer_is_diagonal_positive_and_gauss_restricted():
    basis = build_q6_integer_cycle_basis()
    states = np.vstack([np.zeros((1, 192), dtype=np.int64), basis[:8], -basis[:8]])
    values = q6_u1_electric_transfer_eigenvalues(states, 12.0)
    assert values.shape == (17,)
    assert values[0] == pytest.approx(1.0)
    assert np.all(values > 0.0)
    assert np.all(values <= 1.0 + 1e-12)


def test_electric_hamiltonian_approaches_quadratic_flux_energy():
    basis = build_q6_integer_cycle_basis()
    state = basis[0:1]
    beta, at = 100.0, 0.1
    exact = q6_u1_electric_hamiltonian_eigenvalues(state, beta, at)[0]
    leading = np.sum(state[0] ** 2) / (2.0 * beta * at)
    assert abs(exact - leading) / leading < 0.025


def test_q6_spatial_wilson_potential_is_gauge_invariant():
    from src.rgs_q6_euclidean_time_wilson import build_q6_euclidean_time_data
    data = build_q6_euclidean_time_data(2)
    rng = np.random.default_rng(20261104)
    links = rng.normal(size=192)
    gauge = rng.normal(size=64)
    moved = links + data["D0"].T @ gauge
    assert q6_u1_spatial_wilson_potential(links, 0.8, data) == pytest.approx(
        q6_u1_spatial_wilson_potential(moved, 0.8, data), abs=1e-10
    )


def test_symmetric_transfer_kernel_is_hermitian():
    basis = build_q6_integer_cycle_basis()
    states = np.vstack([np.zeros((1, 192), dtype=np.int64), basis[:12], -basis[:12]])
    rng = np.random.default_rng(7)
    a, b = rng.normal(size=192) * 0.1, rng.normal(size=192) * 0.1
    ab = q6_u1_symmetric_transfer_kernel(b, a, states, 9.0, 0.5)
    ba = q6_u1_symmetric_transfer_kernel(a, b, states, 9.0, 0.5)
    assert ab == pytest.approx(np.conjugate(ba), abs=1e-11)


def test_anisotropic_couplings_have_correct_d3_and_d6_scaling():
    d3 = q6_u1_anisotropic_wilson_couplings(0.2, 0.05, 0.003, spatial_dim=3)
    assert d3["beta_time"] == pytest.approx(0.2 / (0.003 * 0.05))
    assert d3["beta_space"] == pytest.approx(0.05 / (0.003 * 0.2))
    d6 = q6_u1_anisotropic_wilson_couplings(0.2, 0.05, 0.003, spatial_dim=6)
    assert d6["beta_time"] == pytest.approx(0.2**4 / (0.003 * 0.05))
    assert d6["beta_space"] == pytest.approx(0.05 * 0.2**2 / 0.003)


def test_transfer_operator_audit_passes():
    report = audit_q6_u1_transfer_operator()
    assert report["all_checks_passed"], report["checks"]


def test_transfer_inputs_reject_nonpositive_couplings():
    with pytest.raises(ValueError):
        q6_u1_electric_transfer_eigenvalues(np.zeros((1, 192), dtype=int), 0.0)
    with pytest.raises(ValueError):
        q6_u1_anisotropic_wilson_couplings(0.2, 0.05, -1.0, spatial_dim=6)

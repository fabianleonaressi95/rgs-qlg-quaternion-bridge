import numpy as np
import pytest

from src.q6_bisheet_entanglement import (
    analytic_q6_entropy,
    analytic_schmidt_spectrum,
    audit_q6_bisheet_entanglement,
    normalized_bisheet_state,
    von_neumann_entropy,
)


def test_q6_bisheet_maximally_entangled_at_zero_tau():
    C, rho = normalized_bisheet_state(0.0)
    assert C.shape == (64, 64)
    assert np.isclose(von_neumann_entropy(rho), np.log(64), atol=1e-12)
    assert np.allclose(rho, np.eye(64) / 64, atol=1e-12)


def test_q6_bisheet_entropy_matches_exact_spectrum():
    report = audit_q6_bisheet_entanglement()
    assert report["pass"]
    assert report["max_entropy_residual"] < 1e-10
    assert report["max_schmidt_spectrum_residual"] < 1e-10


def test_entanglement_decreases_with_heat_kernel_scale():
    values = [analytic_q6_entropy(t) for t in (0.0, 0.1, 0.25, 0.5, 1.0, 2.0)]
    assert all(a > b for a, b in zip(values, values[1:]))


def test_schmidt_probabilities_are_normalized():
    for tau in (0.0, 0.2, 1.0):
        p = analytic_schmidt_spectrum(tau)
        assert p.shape == (64,)
        assert np.all(p >= 0)
        assert np.isclose(p.sum(), 1.0)


@pytest.mark.parametrize("tau", [-1.0, float("nan"), float("inf")])
def test_invalid_tau_rejected(tau):
    with pytest.raises(ValueError):
        normalized_bisheet_state(tau)

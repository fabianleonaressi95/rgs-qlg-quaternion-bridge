import numpy as np
import pytest

from src.q6_complex import build_q6_complex
from src.q6_gibbs_dynamics import (
    integrated_stress_correlation,
    sample_q6_gibbs_stress,
    stress_autocorrelation,
)


def test_q6_gibbs_sampler_returns_finite_trace_and_valid_acceptance():
    trace = sample_q6_gibbs_stress(
        sweeps=3, burn_in=1, thin=1, beta=0.7,
        proposal_width=0.08, reference_volume=3.0, seed=77,
    )
    assert trace.stress.shape == (3,)
    assert trace.action.shape == (3,)
    assert np.all(np.isfinite(trace.stress))
    assert np.all(np.isfinite(trace.action))
    assert 0.0 <= trace.acceptance_rate <= 1.0
    assert trace.proposed == (3 + 1) * 192


def test_q6_gibbs_sampler_is_reproducible_for_fixed_seed():
    kwargs = dict(
        sweeps=4, burn_in=1, beta=0.5,
        proposal_width=0.06, reference_volume=2.0, seed=901,
    )
    first = sample_q6_gibbs_stress(**kwargs)
    second = sample_q6_gibbs_stress(**kwargs)
    assert np.array_equal(first.action, second.action)
    assert np.array_equal(first.stress, second.stress)
    assert first.accepted == second.accepted


def test_sampler_rejects_wrong_link_graph():
    from src.holonomies import SU2Lattice
    with pytest.raises(ValueError):
        sample_q6_gibbs_stress(
            sweeps=1, initial_links=SU2Lattice([0, 1]),
        )


def test_autocorrelation_is_normalized_and_finite():
    values = np.array([1.0, -1.0, 0.5, -0.5, 1.5, -1.5])
    corr = stress_autocorrelation(values)
    assert corr[0] == 1.0
    assert np.all(np.isfinite(corr))


def test_integrated_correlation_is_algorithmic_not_physical():
    values = np.array([0.2, -0.4, 0.7, -0.1, 0.3, -0.5, 0.8])
    result = integrated_stress_correlation(values, max_lag=3, algorithmic_step=0.25)
    assert result["correlation"].shape == (4,)
    assert result["lags_algorithmic"][-1] == 0.75
    assert "not physical" in result["status"]


def test_constant_trace_has_no_defined_normalized_autocorrelation():
    with pytest.raises(ValueError):
        stress_autocorrelation(np.ones(8))


def test_face_count_and_link_count_match_q6_contract():
    q6 = build_q6_complex()
    assert len(q6.edges) == 192
    assert len(q6.faces) == 240

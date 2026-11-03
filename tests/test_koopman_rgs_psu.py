import numpy as np
import pytest

from src.koopman_rgs_psu import (
    barycentric_koopman_features, build_cantor_q6_markov_chain,
    constant_frequency_koopman_generator, edmd_koopman_audit,
    koopman_markov_audit, koopman_semigroup_residual,
    build_q6_rgs_geometry_lookup,
)

def test_cantor_parry_q6_chain_is_stochastic_and_stationary():
    P, mu, meta = build_cantor_q6_markov_chain()
    assert P.shape == (192, 192)
    assert np.min(P) >= 0
    assert np.max(np.abs(P.sum(axis=1) - 1.0)) < 1e-12
    assert np.max(np.abs(P.T @ mu - mu)) < 1e-12
    assert abs(mu.sum() - 1.0) < 1e-14
    assert abs(meta['psi']**3 - meta['psi']**2 - 1.0) < 1e-12
    assert np.max(np.abs(meta['parry_transition'].sum(axis=1) - 1.0)) < 1e-12

def test_markov_koopman_and_weighted_adjoint_audit():
    P, mu, _ = build_cantor_q6_markov_chain()
    r = koopman_markov_audit(P, mu)
    assert r['dimension'] == 192
    assert r['row_stochastic_residual'] < 1e-12
    assert r['stationarity_residual'] < 1e-12
    assert r['constant_observable_residual'] < 1e-12
    assert r['weighted_adjoint_identity_residual'] < 1e-12
    assert r['spectral_radius'] <= 1.0 + 1e-10

def test_markov_audit_rejects_nonstochastic_matrix():
    with pytest.raises(ValueError):
        koopman_markov_audit([[0.2, 0.2], [0.0, 1.0]])

def test_edmd_recovers_exact_affine_koopman_on_linear_system():
    x = np.linspace(-1.0, 1.0, 31)
    X = np.column_stack((np.ones_like(x), x))
    Y = np.column_stack((np.ones_like(x), 2.0*x + 1.0))
    r = edmd_koopman_audit(X, Y)
    K = np.asarray(r['koopman_matrix'])
    expected = np.array([[1.0, 1.0], [0.0, 2.0]])
    assert np.allclose(K, expected, atol=1e-12)
    assert r['weighted_relative_residual'] < 1e-12
    assert r['rank'] == 2

def test_edmd_reports_nonclosure_for_insufficient_dictionary():
    x = np.linspace(-1.0, 1.0, 41)
    r = edmd_koopman_audit(x[:, None], (x**2)[:, None])
    assert r['weighted_relative_residual'] > 1e-3

def test_barycentric_feature_map_has_exact_centered_moire_observables():
    X = np.array([[[0., 0.], [1., 0.], [0., 1.]],
                  [[0., 0.], [2., 0.], [0., 2.]]])
    W = np.array([[1., 2., 1.], [1., 1., 2.]])
    dk = np.array([0.2, -0.1])
    F, names = barycentric_koopman_features(X, W, dk, delta_omega=0.3,
                                             times=[0.0, 1.0], phase0=0.2)
    assert F.shape == (2, len(names))
    assert names[0] == '1' and np.all(F[:, 0] == 1.0)
    assert np.all(np.isfinite(F))
    cr, ci = F[:, names.index('moire_char_real')], F[:, names.index('moire_char_imag')]
    cp, sp = F[:, names.index('cos_phase')], F[:, names.index('sin_phase')]
    exact = 1.0 - (cp*cr - sp*ci)
    assert np.allclose(F[:, names.index('moire_average_exact')], exact, atol=1e-13)

def test_constant_frequency_koopman_generator_is_skew_and_semigroup_consistent():
    L = constant_frequency_koopman_generator(1.37)
    assert np.allclose(L + L.T, 0.0, atol=1e-14)
    assert koopman_semigroup_residual(L, 0.07, steps=9) < 1e-10

def test_invalid_edmd_inputs_are_rejected():
    with pytest.raises(ValueError):
        edmd_koopman_audit([[1.0], [2.0]], [[1.0]])
    with pytest.raises(ValueError):
        edmd_koopman_audit([[1.0], [2.0]], [[1.0], [2.0]], sample_weights=[1.0, -1.0])

def test_q6_rgs_lookup_reindexes_archived_big_endian_rows_exactly():
    canonical = np.array([[(g >> a) & 1 for a in range(6)] for g in range(64)])
    # Simulate the archive's lexicographic row order (MSB-first).
    archived_states = canonical[np.argsort(canonical @ (2 ** np.arange(6))[::-1])]
    # Coordinates are deliberately nonlinear in the six bits.
    X = np.column_stack((
        np.sum(archived_states, axis=1),
        archived_states[:, 0] * archived_states[:, 3] + 0.25 * archived_states[:, 5],
        np.sum(archived_states * (2 ** np.arange(6)), axis=1) ** 2 / 4096.0,
    ))
    result = build_q6_rgs_geometry_lookup(X, archived_states)
    lookup = np.asarray(result["positions_by_q6_integer"])
    ids = archived_states @ (2 ** np.arange(6))
    assert np.allclose(lookup[ids], X, atol=0.0)
    assert result["geometry_map"].startswith("exact 64-state lookup")
    assert result["affine_relative_residual"] > 1e-3


def test_q6_rgs_lookup_rejects_duplicate_or_missing_states():
    states = np.array([[(g >> a) & 1 for a in range(6)] for g in range(64)])
    states[-1] = states[0]
    with pytest.raises(ValueError):
        build_q6_rgs_geometry_lookup(np.zeros((64, 3)), states)

import numpy as np
import pytest

from src.rgs_weyl_c6_q6 import (
    build_q6_c6_quadratic_operator,
    positive_c6_constitutive_matrix,
    q6_c6_energy,
    q6_directional_gradient_matrix,
)


def test_directional_gradient_has_six_hypercube_edges_per_vertex():
    d = q6_directional_gradient_matrix()
    assert d.shape == (384, 64)
    assert np.all(np.sum(d != 0, axis=1) == 2)
    assert np.all(np.sum(d, axis=1) == 0)


def test_positive_c6_constitutive_matrix_uses_offdiagonal_channels():
    h, report = positive_c6_constitutive_matrix(beta=0.2)
    assert report["strictly_positive"]
    assert report["C6_offdiagonal_norm"] > 0
    assert np.linalg.eigvalsh(h).min() > 0


def test_q6_operator_is_symmetric_positive_and_has_constant_kernel():
    k, report = build_q6_c6_quadratic_operator(beta=0.2)
    assert k.shape == (64, 64)
    assert report["symmetry_residual"] < 1e-12
    assert report["constant_mode_residual"] < 1e-11
    assert np.linalg.eigvalsh(k).min() > -1e-10
    assert np.linalg.matrix_rank(k, tol=1e-9) == 63


def test_energy_is_nonnegative_and_zero_on_constants():
    rng = np.random.default_rng(20261104)
    f = rng.normal(size=64)
    assert q6_c6_energy(f, beta=0.15) >= 0
    assert q6_c6_energy(np.ones(64), beta=0.15) == pytest.approx(0.0, abs=1e-10)


def test_axis_relabelling_is_covariant_when_mapping_is_relabelled():
    rng = np.random.default_rng(42)
    f = rng.normal(size=64)
    identity = tuple(range(6))
    permutation = (2, 5, 1, 4, 0, 3)
    e0 = q6_c6_energy(f, beta=0.13, axis_map=identity)
    e1 = q6_c6_energy(f, beta=0.13, axis_map=permutation)
    # Reordering channels alone is a different assignment; covariance is tested
    # directly at the constitutive-matrix level under simultaneous relabelling.
    h0, _ = positive_c6_constitutive_matrix(beta=0.13, axis_map=identity)
    h1, _ = positive_c6_constitutive_matrix(beta=0.13, axis_map=permutation)
    p = np.eye(6)[list(permutation)]
    assert np.allclose(h1, p @ h0 @ p.T, atol=1e-12)
    assert np.isfinite(e0) and np.isfinite(e1)


def test_invalid_axis_map_and_nonpositive_parameters_rejected():
    with pytest.raises(ValueError):
        positive_c6_constitutive_matrix(axis_map=(0, 1, 2, 3, 4, 4))
    with pytest.raises(ValueError):
        positive_c6_constitutive_matrix(epsilon=0.0)

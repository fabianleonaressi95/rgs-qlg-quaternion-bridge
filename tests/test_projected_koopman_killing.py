import numpy as np
import pytest

from src.koopman_rgs_psu import (
    discrete_euclidean_killing_audit,
    project_q6_laplacian_to_rgs_tangent,
)


def test_projected_q6_graph_diffusion_is_symmetric_negative_semidefinite():
    rng = np.random.default_rng(2718)
    Q, _ = np.linalg.qr(rng.normal(size=(192, 188)))
    L6 = np.zeros((64, 64))
    for g in range(64):
        for bit in range(6):
            h = g ^ (1 << bit)
            L6[g, g] += 1.0
            L6[g, h] -= 1.0
    result = project_q6_laplacian_to_rgs_tangent(Q, L6)
    G = result["generator"]
    assert G.shape == (188, 188)
    assert np.linalg.norm(G - G.T) < 1e-11
    assert np.linalg.eigvalsh(G)[-1] < 1e-10
    assert result["tangent_orthogonality_residual"] < 1e-11
    assert result["projected_min_eigenvalue"] >= -1e-10
    assert result["projected_positive_modes"] == 188


def test_projected_q6_graph_diffusion_rejects_nonorthonormal_frame():
    with pytest.raises(ValueError):
        project_q6_laplacian_to_rgs_tangent(np.ones((192, 188)))


def test_discrete_euclidean_rotations_are_edgewise_killing():
    X = np.array([[(g >> a) & 1 for a in range(3)] for g in range(64)], dtype=float)
    edges = np.array([(g, g ^ (1 << a)) for g in range(64) for a in range(3) if not (g >> a) & 1])
    result = discrete_euclidean_killing_audit(X, edges)
    assert result["passed"]
    assert result["max_killing_residual"] < 1e-12
    assert result["edges"] == 96


def test_discrete_killing_rejects_invalid_edge_endpoint():
    X = np.zeros((64, 3))
    with pytest.raises(ValueError):
        discrete_euclidean_killing_audit(X, [[0, 64]])

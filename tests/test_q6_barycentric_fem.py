import numpy as np
import pytest

from src.q6_barycentric_fem import (
    audit_q6_barycentric_fem,
    build_q6_barycentric_fem,
    q6_dirichlet_energy,
    q6_edge_gradients,
    q6_hypercube_coordinates,
)


def test_q6_hypercube_embedding_and_fem_matrix_invariants():
    fem = build_q6_barycentric_fem()
    report = audit_q6_barycentric_fem(fem)
    assert report["vertices"] == 64
    assert report["edges"] == 192
    assert report["embedding_dimension"] == 6
    assert report["edge_length_min"] == pytest.approx(1.0)
    assert report["edge_length_max"] == pytest.approx(1.0)
    assert np.allclose(fem.barycentre, np.full(6, 0.5))
    assert np.allclose(fem.covariance, 0.25 * np.eye(6))
    assert report["mass_positive_definite"]
    assert report["stiffness_semidefinite"]
    assert report["stiffness_kernel_dimension"] == 1
    assert report["stiffness_constant_mode_residual"] < 1e-12
    assert report["stiffness_graph_factorization_residual"] < 1e-12
    assert report["boundary_squared_exact_zero"]


def test_q6_coordinate_fields_have_expected_edge_gradients_and_energy():
    fem = build_q6_barycentric_fem()
    coordinates = q6_hypercube_coordinates()
    for axis in range(6):
        gradient = q6_edge_gradients(fem, coordinates[:, axis])
        assert gradient.shape == (192, 6)
        assert np.allclose(gradient[:, axis], 1.0)
        assert np.allclose(np.delete(gradient, axis, axis=1), 0.0)
        assert q6_dirichlet_energy(fem, coordinates[:, axis]) == pytest.approx(16.0)


def test_translation_invariance_and_anisotropic_embedding():
    x = q6_hypercube_coordinates()
    base = build_q6_barycentric_fem(x)
    shifted = build_q6_barycentric_fem(x + np.array([3., -2., 1., 0., 4., -1.]))
    assert np.allclose(base.mass, shifted.mass)
    assert np.allclose(base.stiffness, shifted.stiffness)
    scaled = build_q6_barycentric_fem(x @ np.diag([1., 2., 3., 4., 5., 6.]))
    assert scaled.edge_lengths.min() == pytest.approx(1.0)
    assert scaled.edge_lengths.max() == pytest.approx(6.0)
    assert np.linalg.eigvalsh(scaled.mass).min() > 0.0


def test_degenerate_q6_embedding_is_rejected():
    x = np.zeros((64, 6))
    with pytest.raises(ValueError, match="degenerate"):
        build_q6_barycentric_fem(x)

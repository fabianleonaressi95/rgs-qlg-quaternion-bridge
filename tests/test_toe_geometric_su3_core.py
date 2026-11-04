import numpy as np
import pytest

from src.barycentric_geometry import (
    audit_simplex, barycentric_coordinates, p1_element_matrices,
    simplex_geometry, weighted_barycentre,
)
from src.su3_flavour import su3_three_fundamental_audit


def test_tetrahedron_barycentric_affine_and_fem_invariants():
    vertices = np.array([[0., 0., 0.], [2., 0., 0.], [0., 3., 0.], [0., 0., 4.]])
    report = audit_simplex(vertices)
    assert report["volume"] == pytest.approx(4.0)
    assert report["partition_of_unity_residual"] < 1e-12
    assert report["gradient_sum_residual"] < 1e-12
    assert report["mass_positive_definite"]
    assert report["stiffness_positive_semidefinite"]
    assert report["stiffness_constant_mode_residual"] < 1e-12


def test_barycentric_reconstructs_points_and_weighted_centre():
    vertices = np.array([[0., 0.], [2., 0.], [0., 2.]])
    point = np.array([0.5, 0.75])
    lam = barycentric_coordinates(vertices, point)
    assert np.allclose(lam @ vertices, point)
    assert np.sum(lam) == pytest.approx(1.0)
    centre, centred, covariance = weighted_barycentre(vertices, [1, 2, 1])
    assert np.allclose(centred.T @ np.array([.25, .5, .25]), 0.0)
    assert np.allclose(centre, np.array([1.0, .5]))
    assert np.linalg.eigvalsh(covariance).min() >= -1e-12


def test_degenerate_simplex_rejected():
    with pytest.raises(ValueError, match="degenerate"):
        simplex_geometry(np.array([[0., 0.], [1., 0.], [2., 0.]]))


def test_su3_three_fundamental_representation_content():
    report = su3_three_fundamental_audit()
    assert report["dimension_decomposition_exact"]
    assert report["total_casimir_eigenvalue_multiplicities"] == {
        "0.0": 1, "3.0": 16, "6.0": 10
    }
    assert report["gell_mann_trace_orthogonality_residual"] < 1e-12
    assert report["mixed_projector_idempotence_residual"] < 1e-12
    assert report["total_casimir_commutator_residual"] < 1e-10

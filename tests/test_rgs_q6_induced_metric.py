import numpy as np

from src.rgs_q6_induced_metric import (
    audit_rgs_q6_induced_metric,
    build_rgs_q6_induced_metric,
)


def test_rgs_q6_metric_uses_hash_verified_seed_and_canonical_cells():
    data = build_rgs_q6_induced_metric()
    assert data.coordinates.shape == (64, 3)
    assert data.bit_coordinates.shape == (64, 6)
    assert data.edges.shape == (192, 2)
    assert data.faces.shape == (240, 4)
    assert np.all(data.edge_lengths > 0.0)


def test_local_six_dimensional_pullback_is_explicitly_degenerate():
    data = build_rgs_q6_induced_metric()
    for jac, gram in zip(data.local_jacobians, data.local_pullback_grams):
        assert jac.shape == (3, 6)
        assert gram.shape == (6, 6)
        assert np.linalg.matrix_rank(jac, tol=1e-10) <= 3
        assert np.linalg.matrix_rank(gram, tol=1e-9) <= 3
        assert np.allclose(gram, jac.T @ jac, atol=1e-12)


def test_face_gram_and_quadrilateral_area_are_consistent():
    data = build_rgs_q6_induced_metric()
    assert data.face_grams.shape == (240, 2, 2)
    assert np.all(np.linalg.eigvalsh(data.face_grams) >= -1e-12)
    assert np.all(data.face_parallelogram_areas >= 0.0)
    assert np.all(data.face_triangle_areas >= 0.0)
    assert np.all(data.face_closure_defects >= 0.0)


def test_metric_weighted_graph_p1_mass_and_stiffness():
    data = build_rgs_q6_induced_metric()
    mass_eigs = np.linalg.eigvalsh(data.edge_mass)
    stiffness_eigs = np.linalg.eigvalsh(data.edge_stiffness)
    assert mass_eigs[0] > 0.0
    assert stiffness_eigs[0] >= -1e-10
    assert np.linalg.norm(data.edge_stiffness @ np.ones(64)) < 1e-10
    assert np.allclose(data.edge_mass, data.edge_mass.T, atol=1e-12)
    assert np.allclose(data.edge_stiffness, data.edge_stiffness.T, atol=1e-12)


def test_audit_keeps_coordinate_oh_separate_from_q6_equivariance():
    report = audit_rgs_q6_induced_metric()
    assert report["counts"] == {"vertices": 64, "edges": 192, "square_faces": 240}
    assert report["checks"]["coordinate_Oh_action_verified"]
    assert report["pullback_6d"]["rank_max"] <= 3
    assert report["evidence"]["continuum_hodge_mosco"] == "OPEN"
    assert report["pass"]

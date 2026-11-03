"""Tests for Q6 spectral and homological formalization."""
import numpy as np
from src.q6_complex import build_q6_complex
from src.q6_spectral_audit import (
    antipodal_vertex_map, audit_q6_formalization, expected_q6_spectrum,
    graph_laplacian, homology_audit, spectral_audit,
)


def test_q6_laplacian_exact_hypercube_spectrum():
    result = spectral_audit()
    assert result["pass"]
    assert result["observed_multiplicities"] == expected_q6_spectrum()
    assert result["max_abs_spectral_error"] < 1e-10


def test_graph_cycle_rank_differs_from_cellular_first_betti():
    result = homology_audit()
    assert result["pass"]
    assert result["graph_cycle_rank_beta1"] == 129
    assert result["cellular_betti_b1"] == 0
    assert result["cellular_betti_b2"] == 111


def test_antipodal_complement_is_graph_automorphism():
    q6 = build_q6_complex()
    p = antipodal_vertex_map()
    assert np.array_equal(p[p], np.arange(64))
    edges = {tuple(e) for e in q6.edges}
    assert all(tuple(sorted((int(p[u]), int(p[v])))) in edges for u, v in q6.edges)


def test_combined_audit_passes_but_keeps_geometric_identification_open():
    report = audit_q6_formalization()
    assert report["full_pass"]
    assert report["vertex_count"] == 64
    assert report["edge_count"] == 192
    assert report["face_count"] == 240
    assert report["boundary_squared_exact_zero"]
    assert report["geometric_embedding_Oh_identification"].startswith("OPEN")
    assert report["gauge_dynamics_on_Q6"].startswith("OPEN")


def test_graph_laplacian_is_symmetric_positive_semidefinite():
    L = graph_laplacian()
    assert np.array_equal(L, L.T)
    assert np.linalg.eigvalsh(L.astype(float)).min() > -1e-10

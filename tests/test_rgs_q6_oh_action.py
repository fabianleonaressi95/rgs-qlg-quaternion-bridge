"""Tests for the induced O_h action on the recovered RGS seed."""
import numpy as np

from src.rgs_ellis_refinement import load_rgs_seed_data
from src.rgs_q6_oh_action import (
    audit_rgs_q6_oh,
    induced_vertex_permutation,
    signed_permutation_matrices,
)


def test_signed_permutation_group_has_order_48():
    group = signed_permutation_matrices()
    assert len(group) == 48
    assert len({tuple(g.ravel()) for g in group}) == 48
    assert all(np.array_equal(g.T @ g, np.eye(3, dtype=int)) for g in group)


def test_all_oh_matrices_act_as_unique_coordinate_permutations():
    seed = load_rgs_seed_data()
    for matrix in signed_permutation_matrices():
        permutation, audit = induced_vertex_permutation(seed.coordinates, matrix)
        assert np.array_equal(np.sort(permutation), np.arange(64))
        assert audit["max_coordinate_residual"] < 1e-10
        assert audit["minimum_match_gap"] > 1e-10


def test_induced_oh_action_is_faithful_and_respects_group_product():
    report = audit_rgs_q6_oh()
    assert report["group"]["order"] == 48
    assert report["group"]["faithful_on_64_vertices"]
    assert report["group"]["action_homomorphism_exact"]
    assert report["antipodal"]["involution_exact"]
    assert report["antipodal"]["fixed_point_count"] == 0
    assert report["pass"]


def test_q6_equivariance_is_reported_separately_from_coordinate_oh():
    report = audit_rgs_q6_oh()
    q6 = report["q6"]
    assert q6["vertices"] == 64
    assert q6["edges"] == 192
    assert q6["square_faces"] == 240
    assert 0 <= q6["group_elements_preserving_edges"] <= 48
    assert 0 <= q6["group_elements_preserving_faces"] <= 48
    # Coordinate O_h correctness must not be conflated with Q6-label equivariance.
    assert "induced_Q6_graph_action" in report["status"]

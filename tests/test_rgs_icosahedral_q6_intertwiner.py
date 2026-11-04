import numpy as np
import pytest
from src.rgs_icosahedral_q6_intertwiner import (
    audit_icosahedral_q6_intertwiner,
    icosahedral_signed_axis_representation,
    q6_vertex_action,
)


def test_proper_icosahedral_group_has_60_signed_axis_actions():
    data = icosahedral_signed_axis_representation()
    assert data["rotations"].shape == (60, 3, 3)
    assert data["signed_permutations"].shape == (60, 6, 6)
    for s in data["signed_permutations"]:
        assert np.array_equal(np.sum(np.abs(s), axis=0), np.ones(6))
        assert np.array_equal(np.sum(np.abs(s), axis=1), np.ones(6))


def test_signed_axis_action_preserves_c6_and_is_group_representation():
    report = audit_icosahedral_q6_intertwiner()
    assert report["C6_equivariance_residual"] < 1e-8
    assert report["group_homomorphism_residual"] < 1e-8
    assert report["all_checks_passed"]


def test_signed_icosahedral_actions_are_q6_graph_automorphisms():
    for s in icosahedral_signed_axis_representation()["signed_permutations"]:
        image = [q6_vertex_action(v, s) for v in range(64)]
        assert len(set(image)) == 64
        for v in range(64):
            for axis in range(6):
                assert (image[v] ^ image[v ^ (1 << (5-axis))]).bit_count() == 1


def test_q6_vertex_action_rejects_invalid_inputs():
    with pytest.raises(ValueError):
        q6_vertex_action(0, np.zeros((6, 6), dtype=int))
    with pytest.raises(ValueError):
        q6_vertex_action(64, np.eye(6, dtype=int))

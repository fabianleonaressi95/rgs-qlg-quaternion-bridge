from itertools import permutations, product

import pytest

from src.gec6_a5_root_cayley_torus import (
    act_permutation, cayley_neighbors, finite_a5_audit, from_root_coordinates,
    root_step, root_steps, to_root_coordinates,
)


def test_root_coordinate_round_trip():
    for q in product(range(-2, 3), repeat=5):
        assert to_root_coordinates(from_root_coordinates(q)) == q


def test_all_30_roots_are_distinct_and_generate_a5():
    roots = root_steps()
    assert len(roots) == 30
    assert len(set(roots)) == 30
    assert {root_step(i, i + 1) for i in range(5)} == {
        tuple(1 if k == i else 0 for k in range(5)) for i in range(5)
    }


def test_s6_preserves_root_system():
    roots = set(root_steps())
    for sigma in permutations(range(6)):
        assert {act_permutation(r, sigma) for r in roots} == roots


def test_s6_action_is_faithful_on_root_set():
    roots = root_steps()
    identity = tuple(range(6))
    for sigma in permutations(range(6)):
        images = tuple(act_permutation(r, sigma) for r in roots)
        if sigma != identity:
            assert images != roots


def test_s6_action_remains_faithful_modulo_three():
    roots = root_steps()
    identity = tuple(range(6))
    for sigma in permutations(range(6)):
        if sigma == identity:
            continue
        assert any(
            tuple(x % 3 for x in act_permutation(r, sigma))
            != tuple(x % 3 for x in r)
            for r in roots
        )


def test_finite_quotient_is_simple_30_regular():
    m = 3
    for v in product(range(m), repeat=5):
        neighbors = cayley_neighbors(v, m)
        assert len(neighbors) == 30
        assert len(set(neighbors)) == 30
        assert v not in neighbors


def test_finite_quotient_audit_counts():
    result = finite_a5_audit(3)
    assert result["vertex_count"] == 243
    assert result["undirected_edge_count"] == 3645
    assert result["root_action_checks"] == 21600
    assert result["nonidentity_permutations_detected"] == 719
    assert result["vertex_degree_checks"] == 243


def test_invalid_modulus_rejected():
    with pytest.raises(ValueError):
        finite_a5_audit(2)

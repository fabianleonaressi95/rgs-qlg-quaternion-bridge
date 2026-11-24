from itertools import combinations, permutations

from src.q6_prism_exterior_square_map import (
    PAIRS,
    audit,
    edge_cocycle,
    permute_pair_vector,
    square_curvature,
    wedge_basis,
)


def test_all_240_square_curvatures_are_coordinate_wedges():
    count = 0
    for i, j in PAIRS:
        mask = (1 << i) | (1 << j)
        for fixed in range(64):
            if fixed & mask:
                continue
            assert square_curvature(i, j, fixed) == wedge_basis(i, j)
            count += 1
    assert count == 240


def test_all_s6_coordinate_permutations_act_naturally_on_wedge_basis():
    for perm in permutations(range(6)):
        for i, j in PAIRS:
            assert permute_pair_vector(wedge_basis(i, j), perm) == wedge_basis(perm[i], perm[j])


def test_edge_cochain_has_expected_support():
    # Direction 4, with bits 0 and 2 set, gives e0^e4 + e2^e4.
    alpha = edge_cocycle((1 << 0) | (1 << 2), 4)
    expected = [0] * 15
    expected[PAIRS.index((0, 4))] = 1
    expected[PAIRS.index((2, 4))] = 1
    assert tuple(expected) == alpha


def test_audit():
    result = audit()
    assert result["source_cycle_rank"] == 129
    assert result["target_H1_rank"] == 15
    assert result["all_Q6_square_faces_checked"] == 240
    assert result["S6_equivariance_checks"] == 720 * 15
    assert result["A6_equivariant"]
    assert result["surjective_on_added_H1_module"]

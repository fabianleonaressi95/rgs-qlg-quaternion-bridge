from itertools import combinations, permutations

from src.q6_gec6_a5_quotient import (
    audit, edge_cocycle, permute_vector, square_curl, unit_difference,
)


def test_all_240_q6_square_curls_are_a5_roots():
    count = 0
    for i, j in combinations(range(6), 2):
        for fixed in range(64):
            if not ((fixed >> i) & 1) and not ((fixed >> j) & 1):
                assert square_curl(i, j, fixed) == unit_difference(i, j)
                count += 1
    assert count == 240


def test_root_vectors_lie_in_a5():
    for i, j in combinations(range(6), 2):
        vector = unit_difference(i, j)
        assert sum(vector) == 0
        assert sum(abs(x) for x in vector) == 2


def test_s6_equivariance_of_root_assignment():
    for sigma in permutations(range(6)):
        for i, j in combinations(range(6), 2):
            assert permute_vector(unit_difference(i, j), sigma) == unit_difference(
                sigma[i], sigma[j]
            )


def test_edge_cocycle_values_have_zero_coordinate_sum():
    for vertex in range(64):
        for direction in range(6):
            assert sum(edge_cocycle(vertex, direction)) == 0


def test_audit_counts():
    result = audit()
    assert result["square_curl_checks"] == 240
    assert result["S6_root_equivariance_checks"] == 10800
    assert result["kernel_rank"] == 124

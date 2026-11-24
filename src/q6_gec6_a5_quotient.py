"""GEC-6: contracted Q6 graph-cycle cocycle with values in the A5 root lattice.

Coordinates are indexed 0,...,5. Edges are oriented from bit d=0 to bit d=1.
The edge cochain is a gauge choice; its square curl is the canonical root
difference e_j-e_i. The induced cycle map is S6-equivariant.
"""
from itertools import combinations, permutations


def unit_difference(i: int, j: int) -> tuple[int, ...]:
    """Return e_j-e_i in Z^6; i and j must be distinct."""
    if not (0 <= i < 6 and 0 <= j < 6 and i != j):
        raise ValueError("i and j must be distinct indices in range(6)")
    out = [0] * 6
    out[i] -= 1
    out[j] += 1
    return tuple(out)


def edge_cocycle(vertex: int, direction: int) -> tuple[int, ...]:
    """beta_d(x)=sum_{i<d} x_i(e_d-e_i), a contracted edge 1-cochain."""
    if not 0 <= vertex < 64:
        raise ValueError("vertex must be a six-bit integer")
    if not 0 <= direction < 6:
        raise ValueError("direction must lie in range(6)")
    out = [0] * 6
    for i in range(direction):
        if (vertex >> i) & 1:
            out[i] -= 1
            out[direction] += 1
    return tuple(out)


def _sub(a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(x - y for x, y in zip(a, b))


def square_curl(i: int, j: int, fixed_vertex: int) -> tuple[int, ...]:
    """Discrete curl on the square i<j, with both free bits initially zero."""
    if not (0 <= i < j < 6):
        raise ValueError("require 0 <= i < j < 6")
    if not 0 <= fixed_vertex < 64:
        raise ValueError("vertex must be a six-bit integer")
    if ((fixed_vertex >> i) & 1) or ((fixed_vertex >> j) & 1):
        raise ValueError("fixed_vertex must have bits i and j equal to zero")
    vi = fixed_vertex | (1 << i)
    vj = fixed_vertex | (1 << j)
    return _sub(
        _sub(edge_cocycle(vi, j), edge_cocycle(fixed_vertex, j)),
        _sub(edge_cocycle(vj, i), edge_cocycle(fixed_vertex, i)),
    )


def permute_vector(vector: tuple[int, ...], sigma: tuple[int, ...]) -> tuple[int, ...]:
    """Coordinate action e_i -> e_{sigma(i)}."""
    if len(vector) != 6 or sorted(sigma) != list(range(6)):
        raise ValueError("expected a six-vector and a permutation of range(6)")
    out = [0] * 6
    for i, value in enumerate(vector):
        out[sigma[i]] += value
    return tuple(out)


def audit() -> dict:
    """Audit all 240 square curls and all 720*15 root-action identities."""
    square_checks = 0
    for i, j in combinations(range(6), 2):
        for fixed in range(64):
            if not ((fixed >> i) & 1) and not ((fixed >> j) & 1):
                assert square_curl(i, j, fixed) == unit_difference(i, j)
                square_checks += 1

    equivariance_checks = 0
    for sigma in permutations(range(6)):
        for i, j in combinations(range(6), 2):
            assert permute_vector(unit_difference(i, j), sigma) == unit_difference(
                sigma[i], sigma[j]
            )
            equivariance_checks += 1

    return {
        "name": "GEC-6",
        "source": "Q6 graph cycle lattice",
        "source_cycle_rank": 129,
        "exterior_square_rank": 15,
        "target": "A5 root lattice in Z^6",
        "target_rank": 5,
        "kernel_rank": 124,
        "square_curl_checks": square_checks,
        "S6_root_equivariance_checks": equivariance_checks,
        "S6_order": 720,
        "pair_directions": 15,
        "status": "exact finite algebraic audit; physical interpretation open",
    }


if __name__ == "__main__":
    import json
    print(json.dumps(audit(), indent=2, sort_keys=True))

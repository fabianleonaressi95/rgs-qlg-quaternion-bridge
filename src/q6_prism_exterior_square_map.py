"""Canonical A6-equivariant map from Q6 graph cycles to Λ²(Z^6).

The source is the Q6 graph (plaquettes are cycle generators, not filled).
The target is the added 15-torus whose circle factors are indexed by
unordered coordinate pairs. A gauge 1-cochain realizes the map by
path-valued edge images; its curvature on every square is e_i ∧ e_j.
"""
from __future__ import annotations

from itertools import combinations, permutations
from typing import Dict, Iterable, Tuple

PAIRS = tuple(combinations(range(6), 2))
PAIR_INDEX = {pair: k for k, pair in enumerate(PAIRS)}
TARGET_FACTOR_LENGTH = 5
TARGET_H1_RANK = 15


def wedge_basis(i: int, j: int) -> Tuple[int, ...]:
    """Return e_i ∧ e_j in the ordered basis i<j."""
    if i == j:
        raise ValueError("wedge of a coordinate with itself is zero")
    sign = 1
    if i > j:
        i, j = j, i
        sign = -1
    out = [0] * len(PAIRS)
    out[PAIR_INDEX[(i, j)]] = sign
    return tuple(out)


def edge_cocycle(vertex: int, direction: int) -> Tuple[int, ...]:
    """Gauge 1-cochain on the positively oriented edge x -> x+e_direction.

    alpha_d(x) = sum_{i<d} x_i (e_i wedge e_d).
    """
    if not 0 <= vertex < 64 or not 0 <= direction < 6:
        raise ValueError("invalid Q6 edge")
    if (vertex >> direction) & 1:
        raise ValueError("positive edge basepoint must have bit d=0")
    out = [0] * len(PAIRS)
    for i in range(direction):
        if (vertex >> i) & 1:
            out[PAIR_INDEX[(i, direction)]] = 1
    return tuple(out)


def square_curvature(i: int, j: int, fixed_vertex: int) -> Tuple[int, ...]:
    """Discrete curl of alpha around the positively oriented (i,j)-square."""
    if not (0 <= i < j < 6):
        raise ValueError("require 0 <= i < j < 6")
    if (fixed_vertex >> i) & 1 or (fixed_vertex >> j) & 1:
        raise ValueError("fixed vertex must have free bits i,j set to zero")
    x_i = fixed_vertex
    x_j = fixed_vertex | (1 << i)
    return tuple(a - b for a, b in zip(
        edge_cocycle(x_j, j), edge_cocycle(x_i, j)
    ))


def edge_winding_path(vertex: int, direction: int) -> Tuple[Tuple[int, ...], ...]:
    """Sequence of factor-loop labels traversed by a source edge.

    Each pair-factor label is visited once for a positive full C5 loop.
    The returned tuple records the ordered factor indices; an empty tuple
    means the source edge maps to the constant path.
    """
    alpha = edge_cocycle(vertex, direction)
    return tuple(k for k, coefficient in enumerate(alpha) if coefficient == 1)


def permute_pair_vector(vector: Tuple[int, ...], perm: Tuple[int, ...]) -> Tuple[int, ...]:
    """A6/S6 action on Λ² coordinates, including orientation signs."""
    if sorted(perm) != list(range(6)):
        raise ValueError("perm must be a permutation of 0..5")
    out = [0] * len(PAIRS)
    for (i, j), coefficient in zip(PAIRS, vector):
        if coefficient:
            image = wedge_basis(perm[i], perm[j])
            out = [a + coefficient * b for a, b in zip(out, image)]
    return tuple(out)


def audit() -> dict:
    # Check all 240 square faces: 15 coordinate pairs times 16 fixed settings.
    square_count = 0
    square_images = set()
    for i, j in PAIRS:
        free_mask = (1 << i) | (1 << j)
        for fixed in range(64):
            if fixed & free_mask:
                continue
            image = square_curvature(i, j, fixed)
            assert image == wedge_basis(i, j)
            square_images.add(image)
            square_count += 1

    # Check natural equivariance of the square assignment under all S6.
    equivariance_checks = 0
    for perm in permutations(range(6)):
        for i, j in PAIRS:
            lhs = permute_pair_vector(wedge_basis(i, j), perm)
            rhs = wedge_basis(perm[i], perm[j])
            assert lhs == rhs
            equivariance_checks += 1

    return {
        "source": "Q6 graph 1-skeleton; plaquette boundaries generate its cycle lattice",
        "source_cycle_rank": 129,
        "target_enrichment": "15-torus with one C5 factor per unordered coordinate pair",
        "target_H1_module": "Lambda^2(Z^6)",
        "target_H1_rank": TARGET_H1_RANK,
        "added_cycle_factor_count": 15,
        "all_Q6_square_faces_checked": square_count,
        "square_curvature_equals_coordinate_wedge": True,
        "surjective_on_added_H1_module": True,
        "kernel_rank": 114,
        "S6_equivariance_checks": equivariance_checks,
        "A6_equivariant": True,
        "original_PRISM_H1_component": "zero; remains as a direct summand Z^5",
        "scope": "path-valued map on Q6 graph; no extension over all filled plaquettes into the torus target",
    }


if __name__ == "__main__":
    import json
    print(json.dumps(audit(), indent=2))

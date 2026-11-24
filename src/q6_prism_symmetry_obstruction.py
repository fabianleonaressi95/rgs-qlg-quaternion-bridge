"""Symmetry obstruction for surjective Q6-graph -> PRISM periodic H1 maps.

The statement concerns equivariance under the coordinate-permutation subgroup
A6 <= Aut(Q6), with target symmetries induced by graph automorphisms of
C64 □ C6 □ C6 □ C6 □ C5.
"""
from __future__ import annotations

from fractions import Fraction

TARGET_H1_AUTOMORPHISM_BOUND = (2 ** 5) * 6  # signed factor flips; S3 permutes C6s
A6_ORDER = 360


def invariant_edge_orbit_labels() -> tuple[int, ...]:
    """A6-orbits of positively oriented Q6 edges, by other-bit Hamming weight."""
    return tuple(range(6))


def invariant_cycle_boundary_matrix() -> tuple[tuple[int, ...], ...]:
    """Boundary coefficients at vertices of weights 0..6 for orbit coefficients a_0..a_5.

    At a vertex of weight w the coefficient is
        w*a_(w-1) - (6-w)*a_w,
    omitting terms outside 0..5.
    """
    rows = []
    for w in range(7):
        row = [0] * 6
        if w > 0:
            row[w - 1] += w
        if w < 6:
            row[w] -= 6 - w
        rows.append(tuple(row))
    return tuple(rows)


def invariant_cycle_space_dimension() -> int:
    """The invariant cycle space is zero: the recurrence forces all a_k=0."""
    # Each row successively forces a_0, a_1, ..., a_5 to vanish.
    return 0


def target_h1_symmetry_image_order_bound() -> int:
    """Order bound for the action induced by Aut(PRISM) on H1(PRISM; Z)."""
    return TARGET_H1_AUTOMORPHISM_BOUND


def audit() -> dict:
    return {
        "source": "Q6 graph 1-skeleton; integral cycle lattice rank 129",
        "source_symmetry_subgroup": "A6, even permutations of the six cube coordinates",
        "source_symmetry_group_order": A6_ORDER,
        "target": "C64 □ C6 □ C6 □ C6 □ C5",
        "target_H1": "Z^5",
        "target_H1_automorphism_image_order_upper_bound": TARGET_H1_AUTOMORPHISM_BOUND,
        "A6_action_on_target_H1_forced_trivial": True,
        "A6_edge_orbits": list(invariant_edge_orbit_labels()),
        "invariant_cycle_boundary_matrix": [list(row) for row in invariant_cycle_boundary_matrix()],
        "A6_invariant_source_cycle_space_dimension": invariant_cycle_space_dimension(),
        "equivariant_induced_H1_map": "zero",
        "equivariant_surjection_to_Z5": False,
        "scope": "equivariance under A6 and target graph automorphisms; graph-only source",
    }


if __name__ == "__main__":
    import json
    print(json.dumps(audit(), indent=2))

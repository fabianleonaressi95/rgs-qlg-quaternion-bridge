"""Trigram-block-symmetric path-valued map Q6^(1) -> PRISM^(1).

Bits 0..2 are the lower trigram; bits 3..5 are the upper trigram.
The construction is equivariant under S3(lower) x S3(upper), with S3(lower)
permuting the three C6 target factors. It is graph-only, not a filled-cube map.
"""
from __future__ import annotations

from itertools import permutations
from typing import Iterable

LENGTHS = (64, 6, 6, 6, 5)
BASE = (0, 0, 0, 0, 0)


def q6_edges() -> list[tuple[int, int, int]]:
    """Canonical positively oriented edges (lower endpoint, upper endpoint, bit)."""
    return [
        (v, v | (1 << d), d)
        for v in range(64)
        for d in range(6)
        if not (v >> d) & 1
    ]


def upper_weight(v: int) -> int:
    return ((v >> 3) & 0b111).bit_count()


def edge_winding(v: int, direction: int) -> tuple[int, int, int, int, int]:
    """Winding vector of the edge path, in factor order (C64,C6^3,C5)."""
    if direction >= 3:
        return (0, 0, 0, 0, 0)
    k = upper_weight(v)
    out = [0, 0, 0, 0, 0]
    out[direction + 1] = k
    out[0] = int(k >= 2)
    out[4] = int(k == 3)
    return tuple(out)


def edge_path(v: int, direction: int) -> tuple[tuple[int, ...], ...]:
    """A based edge path realizing edge_winding, with positive factor loops."""
    winding = edge_winding(v, direction)
    point = list(BASE)
    path = [tuple(point)]
    for factor, count in enumerate(winding):
        for _ in range(count):
            for _step in range(LENGTHS[factor]):
                point[factor] = (point[factor] + 1) % LENGTHS[factor]
                path.append(tuple(point))
    return tuple(path)


def square_winding(lower_direction: int, upper_direction: int,
                   other_upper_mask: int) -> tuple[int, int, int, int, int]:
    """Image of an oriented mixed square; mask uses the two other upper bits.

    Orientation: lower edge, upper edge, reversed lower edge, reversed upper edge.
    The upper-direction edges have constant image.
    """
    if lower_direction not in range(3) or upper_direction not in range(3, 6):
        raise ValueError("expected one lower and one upper direction")
    if other_upper_mask not in range(4):
        raise ValueError("other_upper_mask must encode two bits")
    k = other_upper_mask.bit_count()
    out = [0, 0, 0, 0, 0]
    out[lower_direction + 1] = -1
    out[0] = int(k >= 2) - int(k + 1 >= 2)
    out[4] = int(k == 3) - int(k + 1 == 3)
    return tuple(out)


def selected_unimodular_square_images() -> list[tuple[int, ...]]:
    """Five square images forming a Z-basis of target H1."""
    # Three k=0 squares give the C6 generators; k=1 and k=2
    # squares in lower direction 0 then add the C64 and C5 generators.
    return [
        square_winding(0, 3, 0),
        square_winding(1, 3, 0),
        square_winding(2, 3, 0),
        square_winding(0, 3, 1),
        square_winding(0, 3, 3),
    ]


def permute_bits(v: int, permutation: tuple[int, int, int, int, int, int]) -> int:
    out = 0
    for old in range(6):
        if (v >> old) & 1:
            out |= 1 << permutation[old]
    return out


def target_permute(point: tuple[int, ...],
                   lower_permutation: tuple[int, int, int]) -> tuple[int, ...]:
    out = [0] * 5
    out[0], out[4] = point[0], point[4]
    for old in range(3):
        out[lower_permutation[old] + 1] = point[old + 1]
    return tuple(out)


def symmetry_audit() -> dict:
    perms = list(permutations(range(3)))
    checked = 0
    for pl in perms:
        for pu in perms:
            p = pl + tuple(3 + x for x in pu)
            for v, _w, d in q6_edges():
                vp = permute_bits(v, p)
                dp = p[d]
                lhs = tuple(target_permute(x, pl) for x in edge_path(v, d))
                rhs = edge_path(vp, dp)
                if lhs != rhs:
                    raise AssertionError("path assignment is not S3 x S3 equivariant")
                checked += 1
    return {
        "source_group": "S3(lower trigram) x S3(upper trigram)",
        "source_group_order": 36,
        "target_action": "S3 permutes the three equal C6 factors; C64 and C5 fixed",
        "edge_path_equivariance_checks": checked,
        "equivariant": True,
    }


def audit() -> dict:
    edges = q6_edges()
    paths = [edge_path(v, d) for v, _w, d in edges]
    selected = selected_unimodular_square_images()
    # The five rows are the images of five source square cycles.
    import sympy as sp
    matrix = sp.Matrix.hstack(*(sp.Matrix(x) for x in selected))
    return {
        "source": "Q6 graph 1-skeleton with lower/upper trigram bit blocks",
        "source_vertices": 64,
        "source_edges": len(edges),
        "target": "C64 x C6^3 x C5",
        "target_H1": "Z^5",
        "lower_trigram_bits": [0, 1, 2],
        "upper_trigram_bits": [3, 4, 5],
        "edge_paths_are_based_loops": all(p[0] == BASE and p[-1] == BASE for p in paths),
        "nonconstant_source_edge_paths": sum(len(p) > 1 for p in paths),
        "total_target_path_steps": sum(len(p) - 1 for p in paths),
        "selected_square_images": [list(x) for x in selected],
        "selected_minor_determinant": int(matrix.det()),
        "induced_H1_map_is_surjective": abs(int(matrix.det())) == 1,
        "equivariance": symmetry_audit(),
        "scope": "integral graph-homology surjection; no extension over all filled Q6 plaquettes is asserted",
    }


if __name__ == "__main__":
    import json
    print(json.dumps(audit(), indent=2))

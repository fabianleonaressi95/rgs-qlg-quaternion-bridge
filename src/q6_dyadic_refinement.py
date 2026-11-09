"""Exact dyadic subdivision of the canonical Q6 cubical 2-skeleton.

This is a geometric/combinatorial reference refinement, NOT a Kathara/Krist-All
identification. Coordinates are in {0,1/2,1}^6, stored as integers {0,1,2}.
The construction preserves the source Q6 orientation conventions.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, product
import numpy as np

from src.q6_complex import build_q6_complex


@dataclass(frozen=True)
class Q6DyadicRefinement:
    vertices: tuple[tuple[int, ...], ...]
    edges: tuple[tuple[tuple[int, ...], tuple[int, ...]], ...]
    faces: tuple[tuple[tuple[int, ...], ...], ...]
    boundary_1: np.ndarray
    boundary_2: np.ndarray
    S0: np.ndarray
    S1: np.ndarray
    S2: np.ndarray


def _bits(v: int) -> tuple[int, ...]:
    return tuple(int(x) for x in format(v, "06b"))


def _fine_edge_endpoints(axis: int, segment: int, fixed: tuple[int, ...]):
    """Expand five fixed coordinates into a positive-axis fine edge."""
    lo = [None] * 6
    hi = [None] * 6
    k = 0
    for j in range(6):
        if j == axis:
            lo[j], hi[j] = segment, segment + 1
        else:
            lo[j] = hi[j] = fixed[k]
            k += 1
    return tuple(lo), tuple(hi)


def _fine_face(axes: tuple[int, int], fixed: tuple[int, ...],
               sa: int, sb: int):
    """Return a positively oriented (+a,+b) fine square boundary."""
    a, b = axes
    base = [None] * 6
    k = 0
    for j in range(6):
        if j not in axes:
            base[j] = fixed[k]
            k += 1
    p00 = list(base)
    p10 = list(base)
    p11 = list(base)
    p01 = list(base)
    p00[a], p10[a], p11[a], p01[a] = sa, sa + 1, sa + 1, sa
    p00[b], p10[b], p11[b], p01[b] = sb, sb, sb + 1, sb + 1
    return tuple(map(tuple, (p00, p10, p11, p01)))


def build_q6_dyadic_refinement() -> Q6DyadicRefinement:
    """Build fine cells and exact integer chain maps S0,S1,S2.

    Coarse Q6 is the 2-skeleton with 64 vertices, 192 edges and 240 squares.
    Fine cells are all grid cells in its standard midpoint subdivision:
    496 vertices, 1344 edges and 960 squares.
    """
    coarse = build_q6_complex()

    vertices = tuple(
        p for p in product(range(3), repeat=6)
        if sum(x == 1 for x in p) <= 2
    )
    vertex_index = {p: i for i, p in enumerate(vertices)}

    # Every fine edge has one varying coordinate and at most one other
    # coordinate fixed at the midpoint. Its orientation is increasing.
    edge_set = set()
    for axis in range(6):
        other = tuple(j for j in range(6) if j != axis)
        for segment in (0, 1):
            for fixed in product((0, 2), repeat=5):
                edge_set.add(_fine_edge_endpoints(axis, segment, fixed))
            for midpoint_axis in other:
                remaining = tuple(j for j in other if j != midpoint_axis)
                for fixed4 in product((0, 2), repeat=4):
                    assignment = dict(zip(remaining, fixed4))
                    assignment[midpoint_axis] = 1
                    fixed = tuple(assignment[j] for j in other)
                    edge_set.add(_fine_edge_endpoints(axis, segment, fixed))
    edges = tuple(sorted(edge_set))
    edge_index = {e: i for i, e in enumerate(edges)}

    # Each fine square is positively oriented by its ordered coordinate axes.
    faces_list = []
    for a, b in combinations(range(6), 2):
        other = tuple(j for j in range(6) if j not in (a, b))
        for fixed in product((0, 2), repeat=4):
            for sa, sb in product((0, 1), repeat=2):
                faces_list.append(_fine_face((a, b), fixed, sa, sb))
    faces = tuple(faces_list)

    b1 = np.zeros((len(vertices), len(edges)), dtype=np.int64)
    for k, (u, v) in enumerate(edges):
        b1[vertex_index[u], k] = -1
        b1[vertex_index[v], k] = 1

    b2 = np.zeros((len(edges), len(faces)), dtype=np.int64)
    for face_id, face in enumerate(faces):
        for source, target in zip(face, face[1:] + face[:1]):
            key = (source, target) if source < target else (target, source)
            b2[edge_index[key], face_id] += 1 if source < target else -1

    S0 = np.zeros((len(vertices), 64), dtype=np.int64)
    for v in coarse.vertices:
        S0[vertex_index[tuple(2*x for x in _bits(v))], v] = 1

    S1 = np.zeros((len(edges), len(coarse.edges)), dtype=np.int64)
    for parent_id, (u, v) in enumerate(coarse.edges):
        bu, bv = _bits(u), _bits(v)
        axis = next(j for j in range(6) if bu[j] != bv[j])
        midpoint = tuple(bu[j] + bv[j] for j in range(6))
        pu = tuple(2*x for x in bu)
        pv = tuple(2*x for x in bv)
        for child in ((pu, midpoint), (midpoint, pv)):
            key = (child[0], child[1]) if child[0] < child[1] else (child[1], child[0])
            sign = 1 if child[0] < child[1] else -1
            S1[edge_index[key], parent_id] += sign

    fine_face_index = {face: i for i, face in enumerate(faces)}
    S2 = np.zeros((len(faces), len(coarse.faces)), dtype=np.int64)
    for parent_id, parent in enumerate(coarse.faces):
        v0, va, vab, vb = parent
        bits0 = _bits(v0)
        ax = next(j for j in range(6) if _bits(va)[j] != bits0[j])
        bx = next(j for j in range(6) if _bits(vb)[j] != bits0[j])
        fixed = tuple(2*bits0[j] for j in range(6) if j not in (ax, bx))
        for sa, sb in product((0, 1), repeat=2):
            child = _fine_face((ax, bx), fixed, sa, sb)
            S2[fine_face_index[child], parent_id] = 1

    return Q6DyadicRefinement(
        vertices=vertices, edges=edges, faces=faces,
        boundary_1=b1, boundary_2=b2, S0=S0, S1=S1, S2=S2,
    )



def _edge_axis(edge: tuple[tuple[int, ...], tuple[int, ...]]) -> int:
    """Return the unique coordinate direction of an axis-parallel edge."""
    u, v = edge
    return next(j for j in range(6) if u[j] != v[j])


def q6_axis_direction_transfer(refinement: Q6DyadicRefinement | None = None) -> np.ndarray:
    """Return the averaged 6x6 parent-axis to child-axis count matrix.

    Each coarse Q6 direction contains 32 edges. The matrix entry (i,j)
    is the number of child edges in direction i per parent edge in
    direction j. For this axis-preserving midpoint subdivision it is 2 I_6.
    """
    coarse = build_q6_complex()
    fine = refinement or build_q6_dyadic_refinement()
    counts = np.zeros((6, 6), dtype=np.int64)
    for parent_id, (u, v) in enumerate(coarse.edges):
        parent_axis = next(j for j in range(6)
                           if _bits(u)[j] != _bits(v)[j])
        for child_id in np.flatnonzero(fine.S1[:, parent_id]):
            child_axis = _edge_axis(fine.edges[int(child_id)])
            counts[child_axis, parent_axis] += abs(int(fine.S1[child_id, parent_id]))
    return counts // 32

def audit_q6_dyadic_refinement() -> dict:
    """Return exact counts, map shapes, chain residuals and subdivision sums."""
    coarse = build_q6_complex()
    fine = build_q6_dyadic_refinement()
    r1 = fine.boundary_1 @ fine.S1 - fine.S0 @ coarse.boundary_1
    r2 = fine.boundary_2 @ fine.S2 - fine.S1 @ coarse.boundary_2
    return {
        "coarse_counts": [len(coarse.vertices), len(coarse.edges), len(coarse.faces)],
        "fine_counts": [len(fine.vertices), len(fine.edges), len(fine.faces)],
        "map_shapes": [list(fine.S0.shape), list(fine.S1.shape), list(fine.S2.shape)],
        "boundary_squared_fine_max_abs": int(np.max(np.abs(fine.boundary_1 @ fine.boundary_2))),
        "chain_map_S1_residual_max_abs": int(np.max(np.abs(r1))),
        "chain_map_S2_residual_max_abs": int(np.max(np.abs(r2))),
        "all_parent_edges_have_two_children": bool(np.all(np.sum(np.abs(fine.S1), axis=0) == 2)),
        "all_parent_faces_have_four_children": bool(np.all(np.sum(np.abs(fine.S2), axis=0) == 4)),
        "axis_direction_transfer": q6_axis_direction_transfer(fine).tolist(),
        "all_exact": bool(
            np.count_nonzero(r1) == 0 and np.count_nonzero(r2) == 0
            and np.count_nonzero(fine.boundary_1 @ fine.boundary_2) == 0
        ),
    }

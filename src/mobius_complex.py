"""Canonical six-square Möbius-band subcomplex embedded in Q6.

The six face patterns are fixed project data. A '*' denotes a free bit.
For each face, the two free bits are ordered cyclically as 00, 10, 11, 01;
this is essential: lexicographic product order would create crossed diagonals.

Vertex labels retain Q6's six-bit big-endian convention. Edges are oriented
from the smaller integer label to the larger one. Each face boundary is
oriented by its cyclic vertex list.
"""
from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import dataclass

import numpy as np


MOBIUS_FACE_PATTERNS = (
    ("*", "*", "1", "1", "0", "0"),
    ("*", "0", "1", "*", "0", "0"),
    ("0", "0", "*", "*", "0", "0"),
    ("*", "1", "1", "*", "0", "0"),
    ("*", "1", "*", "0", "0", "0"),
    ("0", "*", "*", "0", "0", "0"),
)
_CYCLIC_FREE_BITS = ((0, 0), (1, 0), (1, 1), (0, 1))


@dataclass(frozen=True)
class MobiusComplex:
    vertices: tuple[int, ...]
    edges: tuple[tuple[int, int], ...]
    faces: tuple[tuple[int, int, int, int], ...]
    boundary_1: np.ndarray
    boundary_2: np.ndarray


def _face_vertices(pattern: tuple[str, ...]) -> tuple[int, int, int, int]:
    free_positions = [i for i, bit in enumerate(pattern) if bit == "*"]
    if len(free_positions) != 2:
        raise ValueError("Every Möbius face pattern must contain two wildcards.")
    result = []
    for free_bits in _CYCLIC_FREE_BITS:
        bits = list(pattern)
        for position, bit in zip(free_positions, free_bits):
            bits[position] = str(bit)
        result.append(int("".join(bits), 2))
    return tuple(result)


def build_mobius_complex() -> MobiusComplex:
    """Build the agreed 12-vertex, 18-edge, 6-square Möbius subcomplex."""
    faces = tuple(_face_vertices(p) for p in MOBIUS_FACE_PATTERNS)
    vertices = tuple(sorted({v for face in faces for v in face}))
    edge_set = {
        (min(u, v), max(u, v))
        for face in faces
        for u, v in zip(face, face[1:] + face[:1])
    }
    edges = tuple(sorted(edge_set))
    vertex_index = {v: i for i, v in enumerate(vertices)}
    edge_index = {e: i for i, e in enumerate(edges)}

    b1 = np.zeros((len(vertices), len(edges)), dtype=np.int8)
    for k, (u, v) in enumerate(edges):
        b1[vertex_index[u], k] = -1
        b1[vertex_index[v], k] = 1

    b2 = np.zeros((len(edges), len(faces)), dtype=np.int8)
    for j, face in enumerate(faces):
        for u, v in zip(face, face[1:] + face[:1]):
            edge = (min(u, v), max(u, v))
            b2[edge_index[edge], j] += 1 if (u, v) == edge else -1

    return MobiusComplex(vertices, edges, faces, b1, b2)


def _boundary_components(complex_: MobiusComplex) -> tuple[tuple[int, ...], ...]:
    edge_face_count = np.count_nonzero(complex_.boundary_2, axis=1)
    boundary_edges = [
        complex_.edges[i] for i, count in enumerate(edge_face_count) if count == 1
    ]
    adjacency: dict[int, set[int]] = defaultdict(set)
    for u, v in boundary_edges:
        adjacency[u].add(v)
        adjacency[v].add(u)

    components = []
    unseen = set(adjacency)
    while unseen:
        root = min(unseen)
        queue = [root]
        unseen.remove(root)
        component = []
        while queue:
            u = queue.pop()
            component.append(u)
            for v in adjacency[u] & unseen:
                unseen.remove(v)
                queue.append(v)
        components.append(tuple(sorted(component)))
    return tuple(sorted(components))


def _orientation_propagation(complex_: MobiusComplex) -> tuple[bool, int]:
    """Try to orient all faces consistently across shared interior edges.

    Returns (orientable, number_of_orientation_conflicts).
    """
    edge_faces: dict[tuple[int, int], list[tuple[int, int]]] = defaultdict(list)
    for face_id, face in enumerate(complex_.faces):
        for u, v in zip(face, face[1:] + face[:1]):
            edge = (min(u, v), max(u, v))
            sign = 1 if (u, v) == edge else -1
            edge_faces[edge].append((face_id, sign))

    adjacency: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for incident in edge_faces.values():
        if len(incident) == 2:
            (i, si), (j, sj) = incident
            relative_orientation = -si * sj
            adjacency[i].append((j, relative_orientation))
            adjacency[j].append((i, relative_orientation))

    orientation: dict[int, int] = {}
    conflicts = 0
    for root in range(len(complex_.faces)):
        if root in orientation:
            continue
        orientation[root] = 1
        queue = deque([root])
        while queue:
            i = queue.popleft()
            for j, relative in adjacency[i]:
                required = orientation[i] * relative
                if j not in orientation:
                    orientation[j] = required
                    queue.append(j)
                elif orientation[j] != required:
                    conflicts += 1
    return conflicts == 0, conflicts


def audit_mobius_complex() -> dict[str, int | bool | list[list[int]]]:
    """Return reproducible combinatorial, homological and orientability data."""
    complex_ = build_mobius_complex()
    b1 = complex_.boundary_1.astype(np.int64)
    b2 = complex_.boundary_2.astype(np.int64)
    chain_residual = b1 @ b2
    edge_face_count = np.count_nonzero(b2, axis=1)
    boundary_edges = int(np.count_nonzero(edge_face_count == 1))
    interior_edges = int(np.count_nonzero(edge_face_count == 2))
    orientable, conflicts = _orientation_propagation(complex_)
    rank_b1 = int(np.linalg.matrix_rank(b1.astype(float)))
    rank_b2 = int(np.linalg.matrix_rank(b2.astype(float)))
    beta0 = len(complex_.vertices) - rank_b1
    beta1 = len(complex_.edges) - rank_b1 - rank_b2
    beta2 = len(complex_.faces) - rank_b2
    components = _boundary_components(complex_)
    return {
        "vertices": len(complex_.vertices),
        "edges": len(complex_.edges),
        "faces": len(complex_.faces),
        "boundary_edges": boundary_edges,
        "interior_edges": interior_edges,
        "boundary_components": len(components),
        "boundary_component_vertices": [list(c) for c in components],
        "euler_characteristic": len(complex_.vertices) - len(complex_.edges) + len(complex_.faces),
        "rank_boundary_1": rank_b1,
        "rank_boundary_2": rank_b2,
        "betti_0": beta0,
        "betti_1": beta1,
        "betti_2": beta2,
        "boundary_squared_max_abs": int(np.max(np.abs(chain_residual))),
        "boundary_squared_exact_zero": bool(np.count_nonzero(chain_residual) == 0),
        "orientable": orientable,
        "orientation_conflicts": conflicts,
    }

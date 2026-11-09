"""Source-transcribed Kathara 12-center / 15-primary-line graph.

The 15 source line labels are represented as:
  L1: central vertical polyline 1-2-5-8-11-12
  L2: left vertical polyline 3-6-9
  L3: right vertical polyline 4-7-10
  L4..L15: the twelve cross-links transcribed from the published grid.

Every elementary directed edge is oriented from the lower to the higher
center label. Faces are computed from the planar rotation system induced by
the displayed integer coordinates. This is a combinatorial transcription of
the diagram, not a physical model.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from math import atan2
from typing import Iterable

COORDS: dict[int, tuple[int, int]] = {
    1: (0, 0), 2: (0, 1), 3: (-1, 2), 4: (1, 2),
    5: (0, 3), 6: (-1, 4), 7: (1, 4), 8: (0, 5),
    9: (-1, 6), 10: (1, 6), 11: (0, 7), 12: (0, 8),
}
PRIMARY_LINES: dict[int, tuple[int, ...]] = {
    1: (1, 2, 5, 8, 11, 12),
    2: (3, 6, 9),
    3: (4, 7, 10),
    4: (1, 3),
    5: (2, 3),
    6: (1, 4),
    7: (2, 4),
    8: (9, 12),
    9: (9, 11),
    10: (10, 12),
    11: (10, 11),
    12: (8, 9),
    13: (8, 10),
    14: (5, 6),
    15: (5, 7),
}


@dataclass(frozen=True)
class KatharaGraphAudit:
    vertices: tuple[int, ...]
    directed_edges: tuple[tuple[int, int], ...]
    primary_lines: dict[int, tuple[int, ...]]
    faces: tuple[tuple[int, ...], ...]
    bounded_faces: tuple[tuple[int, ...], ...]
    euler_characteristic: int
    connected: bool
    planar_straight_line_embedding: bool
    acyclic_orientation: bool


def elementary_edges(lines: dict[int, tuple[int, ...]] = PRIMARY_LINES) -> tuple[tuple[int, int], ...]:
    """Expand primary polylines into unique directed elementary segments."""
    undirected: set[tuple[int, int]] = set()
    for path in lines.values():
        if len(path) < 2:
            raise ValueError("each primary line must contain at least two centers")
        if any(v not in COORDS for v in path):
            raise ValueError("unknown Kathara center in primary line")
        for a, b in zip(path, path[1:]):
            if a == b:
                raise ValueError("self-loop in primary line")
            undirected.add(tuple(sorted((a, b))))
    return tuple(sorted(undirected))


def _orient(a: tuple[int, int], b: tuple[int, int], c: tuple[int, int]) -> int:
    return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])


def _proper_crossing(e1: tuple[int, int], e2: tuple[int, int]) -> bool:
    if set(e1) & set(e2):
        return False
    a, b = COORDS[e1[0]], COORDS[e1[1]]
    c, d = COORDS[e2[0]], COORDS[e2[1]]
    o1, o2 = _orient(a, b, c), _orient(a, b, d)
    o3, o4 = _orient(c, d, a), _orient(c, d, b)
    return o1*o2 < 0 and o3*o4 < 0


def _face_walk(edges: tuple[tuple[int, int], ...]) -> tuple[tuple[int, ...], ...]:
    """Enumerate oriented face boundary walks from the planar rotation system."""
    adj: dict[int, list[int]] = defaultdict(list)
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    for u in adj:
        x, y = COORDS[u]
        adj[u].sort(key=lambda v: atan2(COORDS[v][1]-y, COORDS[v][0]-x))

    seen: set[tuple[int, int]] = set()
    faces: list[tuple[int, ...]] = []
    for u in sorted(adj):
        for v in adj[u]:
            if (u, v) in seen:
                continue
            face: list[int] = []
            a, b = u, v
            while (a, b) not in seen:
                seen.add((a, b))
                face.append(a)
                neighbors = adj[b]
                reverse_index = neighbors.index(a)
                c = neighbors[(reverse_index - 1) % len(neighbors)]
                a, b = b, c
            faces.append(tuple(face))
    return tuple(faces)


def audit_kathara_15_line_graph() -> KatharaGraphAudit:
    """Build and audit the 12-center graph and its bounded planar faces."""
    if set(PRIMARY_LINES) != set(range(1, 16)):
        raise AssertionError("the source transcription must contain line IDs 1..15")
    edges = elementary_edges()
    vertices = tuple(sorted(COORDS))
    crossing_pairs = [
        (a, b) for i, a in enumerate(edges) for b in edges[i+1:]
        if _proper_crossing(a, b)
    ]
    # Connectivity by graph traversal.
    adj: dict[int, set[int]] = {v: set() for v in vertices}
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
    reached = set()
    stack = [vertices[0]]
    while stack:
        u = stack.pop()
        if u not in reached:
            reached.add(u)
            stack.extend(adj[u] - reached)
    connected = len(reached) == len(vertices)
    faces = _face_walk(edges)
    # The unbounded face is the unique face with maximal boundary length here.
    outer = max(faces, key=len)
    bounded = tuple(face for face in faces if face != outer)
    # Orientation by increasing center label is a DAG orientation.
    indegree = {v: 0 for v in vertices}
    for u, v in edges:
        indegree[v] += 1
    acyclic = all(u < v for u, v in edges)
    euler = len(vertices) - len(edges) + len(faces)
    return KatharaGraphAudit(
        vertices=vertices,
        directed_edges=edges,
        primary_lines=dict(PRIMARY_LINES),
        faces=faces,
        bounded_faces=bounded,
        euler_characteristic=euler,
        connected=connected,
        planar_straight_line_embedding=not crossing_pairs,
        acyclic_orientation=acyclic,
    )

"""Cellular chain complex induced by filling the ten bounded Kathara faces.

The source diagram supplies a graph; declaring its bounded planar regions to
be 2-cells is an additional mathematical convention.
"""
from __future__ import annotations
from dataclasses import dataclass
from src.kathara_15_line_graph import audit_kathara_15_line_graph, elementary_edges

BOUNDED_FACES = (
    (1, 2, 3), (1, 4, 2), (2, 4, 7, 5), (2, 5, 6, 3),
    (5, 7, 10, 8), (5, 8, 9, 6), (8, 10, 11), (8, 11, 9),
    (9, 11, 12), (10, 12, 11),
)

def _rank_q(matrix):
    a = [[int(x) for x in row] for row in matrix]
    if not a:
        return 0
    m, n, rank = len(a), len(a[0]), 0
    for col in range(n):
        pivot = next((r for r in range(rank, m) if a[r][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        p = a[rank][col]
        for r in range(rank + 1, m):
            q = a[r][col]
            if q:
                for c in range(col + 1, n):
                    numerator = p*a[r][c] - q*a[rank][c]
                    a[r][c] = numerator // _gcd(p, q)
                a[r][col] = 0
        rank += 1
        if rank == m:
            break
    return rank

def _gcd(a, b):
    a, b = abs(a), abs(b)
    while b:
        a, b = b, a % b
    return a or 1

def _transpose(a):
    return [list(row) for row in zip(*a)]

def _matmul(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]

def _add(a, b):
    return [[x+y for x, y in zip(ar, br)] for ar, br in zip(a, b)]

@dataclass(frozen=True)
class KatharaCellularAudit:
    vertices: tuple
    edges: tuple
    faces: tuple
    b1: tuple
    b2: tuple
    b1b2_zero: bool
    rank_b1: int
    rank_b2: int
    betti: tuple
    euler_characteristic: int
    laplacian_0: tuple
    laplacian_1: tuple
    laplacian_2: tuple

def audit_kathara_cellular_complex():
    graph = audit_kathara_15_line_graph()
    expected = {frozenset(f) for f in BOUNDED_FACES}
    actual = {frozenset(f) for f in graph.bounded_faces}
    if expected != actual:
        raise AssertionError("bounded-face inventory mismatch")
    vertices, edges = graph.vertices, elementary_edges()
    edge_index = {e: j for j, e in enumerate(edges)}
    b1 = [[0]*len(edges) for _ in vertices]
    for j, (u, v) in enumerate(edges):
        b1[u-1][j], b1[v-1][j] = -1, 1
    b2 = [[0]*len(BOUNDED_FACES) for _ in edges]
    for j, face in enumerate(BOUNDED_FACES):
        for u, v in zip(face, face[1:] + face[:1]):
            e = tuple(sorted((u, v)))
            if e not in edge_index:
                raise AssertionError(f"face boundary uses absent edge {e}")
            b2[edge_index[e]][j] = 1 if (u, v) == e else -1
    product = _matmul(b1, b2)
    if any(x for row in product for x in row):
        raise AssertionError("B1 @ B2 is not zero")
    b1t, b2t = _transpose(b1), _transpose(b2)
    l0 = _matmul(b1, b1t)
    l1 = _add(_matmul(b1t, b1), _matmul(b2, b2t))
    l2 = _matmul(b2t, b2)
    r1, r2 = _rank_q(b1t), _rank_q(b2t)
    betti = (len(vertices)-r1, len(edges)-r1-r2, len(BOUNDED_FACES)-r2)
    chi = len(vertices)-len(edges)+len(BOUNDED_FACES)
    if chi != betti[0]-betti[1]+betti[2]:
        raise AssertionError("Euler-Poincare identity failed")
    freeze = lambda a: tuple(tuple(row) for row in a)
    return KatharaCellularAudit(
        vertices, edges, BOUNDED_FACES, freeze(b1), freeze(b2), True,
        r1, r2, betti, chi, freeze(l0), freeze(l1), freeze(l2))

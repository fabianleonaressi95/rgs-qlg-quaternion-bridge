"""A coordinate-defined, deliberately coarse chain map K_Kathara -> Q3.

The map is a candidate bridge, not a canonical identification. It sends the
Kathara planar embedding to the square face b2=0 of the trigram cube. All
incidence matrices are integral. Run the accompanying tests before treating
the implementation as computationally verified.
"""
from __future__ import annotations

from dataclasses import dataclass

from src.kathara_15_line_graph import COORDS
from src.kathara_cellular_complex import (
    BOUNDED_FACES, audit_kathara_cellular_complex,
)


def _transpose(a):
    return [list(row) for row in zip(*a)]


def _matmul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def _zero(a):
    return all(x == 0 for row in a for x in row)


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
                    a[r][c] = (p * a[r][c] - q * a[rank][c]) // _gcd(p, q)
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


def _edge_index():
    # Q3 edges oriented from bit k=0 to bit k=1.
    edges = []
    for k in range(3):
        for v in range(8):
            if not (v & (1 << k)):
                edges.append((v, v | (1 << k)))
    return tuple(edges), {e: i for i, e in enumerate(edges)}


def _face_boundaries():
    # Six square faces: choose the varying axes i<j and fix the remaining bit.
    faces = []
    for i in range(3):
        for j in range(i + 1, 3):
            k = 3 - i - j
            for fixed in (0, 1):
                base = fixed << k
                faces.append((base, base | (1 << i),
                              base | (1 << i) | (1 << j),
                              base | (1 << j)))
    return tuple(faces)


def _boundary_matrix_1(vertices, edges):
    index = {v: i for i, v in enumerate(vertices)}
    out = [[0] * len(edges) for _ in vertices]
    for j, (u, v) in enumerate(edges):
        out[index[u]][j] = -1
        out[index[v]][j] = 1
    return out


def _boundary_matrix_2(edges, faces):
    edge_index = {e: i for i, e in enumerate(edges)}
    out = [[0] * len(faces) for _ in edges]
    for j, face in enumerate(faces):
        for u, v in zip(face, face[1:] + face[:1]):
            e = (min(u, v), max(u, v))
            out[edge_index[e]][j] += 1 if (u, v) == e else -1
    return out


def _map_vertex(v):
    x, y = COORDS[v]
    # Two geometric half-plane tests; third trigram bit is fixed to zero.
    bit0 = int(x > 0)       # right half-plane
    bit1 = int(y >= 4)       # upper half-plane relative to AzurA=(0,4)
    bit2 = 0
    return bit0 | (bit1 << 1) | (bit2 << 2)


def _path_chain(start, finish, edges, edge_index):
    """Canonical path: flip differing bits in order 0,1,2."""
    chain = [0] * len(edges)
    current = start
    for k in range(3):
        mask = 1 << k
        if (current & mask) != (finish & mask):
            nxt = current ^ mask
            e = (min(current, nxt), max(current, nxt))
            chain[edge_index[e]] += 1 if current < nxt else -1
            current = nxt
    if current != finish:
        raise AssertionError("bit-flip path did not reach its endpoint")
    return chain


@dataclass(frozen=True)
class KatharaQ3ChainMapAudit:
    source_vertex_map: tuple
    target_vertices: tuple
    target_edges: tuple
    target_faces: tuple
    S0: tuple                 # C0(K) -> C0(Q3), shape 8 x 12
    S1: tuple                 # C1(K) -> C1(Q3), shape 12 x 21
    S2: tuple                 # C2(K) -> C2(Q3), shape 6 x 10
    chain_identity_0: bool    # B1(Q3) S1 = S0 B1(K)
    chain_identity_1: bool    # B2(Q3) S2 = S1 B2(K)
    rank_S0: int
    rank_S1: int
    rank_S2: int
    image_uses_one_square_face: bool


def audit_kathara_to_q3_chain_map():
    source = audit_kathara_cellular_complex()
    qv = tuple(range(8))
    qe, qedge_index = _edge_index()
    qf = _face_boundaries()
    qB1 = _boundary_matrix_1(qv, qe)
    qB2 = _boundary_matrix_2(qe, qf)

    # S0: one-hot vertex map.
    S0 = [[0] * len(source.vertices) for _ in qv]
    mapped = tuple(_map_vertex(v) for v in source.vertices)
    for j, target_v in enumerate(mapped):
        S0[target_v][j] = 1

    # S1: each source edge maps to a canonical Q3 edge path.
    S1 = [[0] * len(source.edges) for _ in qe]
    for j, (u, v) in enumerate(source.edges):
        path = _path_chain(_map_vertex(u), _map_vertex(v), qe, qedge_index)
        for i, value in enumerate(path):
            S1[i][j] = value

    # S2: each mapped face cycle lies in the b2=0 square. Its coefficient
    # is read from the target edge 0 -> 1, which fixes the square orientation.
    S2 = [[0] * len(source.faces) for _ in qf]
    square_index = qf.index((0, 1, 3, 2))
    reference_edge = qedge_index[(0, 1)]
    source_B2 = [list(row) for row in source.b2]
    mapped_face_chains = _matmul(S1, source_B2)
    for j in range(len(source.faces)):
        coefficient = mapped_face_chains[reference_edge][j]
        S2[square_index][j] = coefficient

    lhs0 = _matmul(qB1, S1)
    rhs0 = _matmul(S0, [list(row) for row in source.b1])
    lhs1 = _matmul(qB2, S2)
    rhs1 = mapped_face_chains
    identity0, identity1 = _zero([
        [x - y for x, y in zip(a, b)] for a, b in zip(lhs0, rhs0)
    ]), _zero([
        [x - y for x, y in zip(a, b)] for a, b in zip(lhs1, rhs1)
    ])

    # This construction intentionally factors through a single square face.
    nonzero_face_rows = [
        i for i, row in enumerate(S2) if any(row)
    ]
    return KatharaQ3ChainMapAudit(
        mapped, qv, qe, qf,
        tuple(tuple(row) for row in S0),
        tuple(tuple(row) for row in S1),
        tuple(tuple(row) for row in S2),
        identity0, identity1,
        _rank_q(_transpose(S0)), _rank_q(_transpose(S1)),
        _rank_q(_transpose(S2)),
        nonzero_face_rows == [square_index],
    )

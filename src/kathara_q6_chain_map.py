"""Exact integral Kathara -> Q6 chain map candidate.

The 12 Kathara centres are encoded by their already-declared integer drawing
coordinates: x-column (-2,0,+2) -> 0,1,2 and y-level 0..8 -> its integer
four-bit coordinate. The resulting Q6 label is (column << 4) | level.

This is an explicit chart, not a theorem that the chart is selected by the
RGS metric. It uses no distance threshold or nearest-neighbour assignment.
Edges are sent to root-routed Q6 1-chains P(v)-P(u), where P(g) is the
canonical bit-flip path from 0 to g. Consequently every source cycle maps
to zero as a 1-chain and R2=0 is a valid (degenerate) cellular chain map.
The 10 Kathara 2-cells are therefore collapsed; this is not an embedding.
"""
from __future__ import annotations

from dataclasses import dataclass

from src.kathara_cellular_complex import audit_kathara_cellular_complex
from src.kathara_15_line_graph import COORDS, PRIMARY_LINES
from src.q6_complex import build_q6_complex
from src.rgs_ellis_refinement import load_rgs_seed_data


def _column_code(x: int) -> int:
    table = {-2: 0, 0: 1, 2: 2}
    if x not in table:
        raise ValueError(f"unsupported Kathara x-column {x}")
    return table[x]


def kathara_vertex_to_q6_label(vertex: int) -> int:
    """Exact coordinate code; Q6 integer label is six-bit big-endian."""
    if vertex not in COORDS:
        raise ValueError(f"unknown Kathara vertex {vertex}")
    x, y = COORDS[vertex]
    if not 0 <= y <= 8:
        raise ValueError(f"unsupported Kathara y-level {y}")
    return (_column_code(x) << 4) | int(y)


def _root_path_chain(label: int, edges, edge_index) -> list[int]:
    """Integral chain of the canonical 0 -> label path, flipping axes 0..5."""
    chain = [0] * len(edges)
    current = 0
    for axis in range(6):
        mask = 1 << (5 - axis)
        if label & mask:
            nxt = current ^ mask
            edge = (min(current, nxt), max(current, nxt))
            chain[edge_index[edge]] += 1 if current < nxt else -1
            current = nxt
    if current != label:
        raise AssertionError("root path failed to reach its Q6 label")
    return chain


def _matmul(a, b):
    if not a or not b:
        return []
    return [[sum(a[i][k] * b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def _is_zero(a):
    return all(value == 0 for row in a for value in row)


@dataclass(frozen=True)
class KatharaQ6ChainMapAudit:
    source_vertices: tuple
    source_edges: tuple
    source_faces: tuple
    target_vertices: tuple
    target_edges: tuple
    target_faces: tuple
    vertex_labels: tuple
    R0: tuple
    R1: tuple
    R2: tuple
    chain_identity_0: bool
    chain_identity_1: bool
    injective_on_vertices: bool
    collapsed_source_2_cells: int
    line_edge_membership: tuple
    rgs_coordinates_of_image: tuple
    scope: str


def audit_kathara_to_q6_chain_map() -> KatharaQ6ChainMapAudit:
    """Construct and exactly check R0/R1/R2 over the integers."""
    source = audit_kathara_cellular_complex()
    target = build_q6_complex()
    vertex_index = {v: i for i, v in enumerate(source.vertices)}
    edge_index = {e: i for i, e in enumerate(target.edges)}

    labels = tuple(kathara_vertex_to_q6_label(v) for v in source.vertices)
    if len(set(labels)) != len(labels):
        raise AssertionError("coordinate code is not injective on Kathara centres")

    # R0: one-hot map C0(K) -> C0(Q6).
    R0 = [[0] * len(source.vertices) for _ in target.vertices]
    for j, label in enumerate(labels):
        R0[label][j] = 1

    # Root paths P(label), then edge (u,v) -> P(v)-P(u).
    root_paths = {
        label: _root_path_chain(label, target.edges, edge_index)
        for label in set(labels)
    }
    R1 = [[0] * len(source.edges) for _ in target.edges]
    for j, (u, v) in enumerate(source.edges):
        pu, pv = root_paths[labels[vertex_index[u]]], root_paths[labels[vertex_index[v]]]
        for i in range(len(target.edges)):
            R1[i][j] = pv[i] - pu[i]

    # This chain map deliberately factors through the Q6 1-skeleton.
    # Every closed source 1-cycle telescopes to zero, hence R2=0.
    R2 = [[0] * len(source.faces) for _ in target.faces]

    lhs0 = _matmul([list(row) for row in target.boundary_1], R1)
    rhs0 = _matmul(R0, [list(row) for row in source.b1])
    lhs1 = _matmul([list(row) for row in target.boundary_2], R2)
    rhs1 = _matmul(R1, [list(row) for row in source.b2])
    identity0 = _is_zero([[x-y for x, y in zip(a, b)]
                          for a, b in zip(lhs0, rhs0)])
    identity1 = _is_zero([[x-y for x, y in zip(a, b)]
                          for a, b in zip(lhs1, rhs1)])
    if not identity0 or not identity1:
        raise AssertionError("constructed Kathara -> Q6 chain identities failed")

    # Retain source's 15-line hierarchy as incidence metadata, not weights.
    membership = []
    for line_id, path in sorted(PRIMARY_LINES.items()):
        segments = tuple(tuple(sorted(pair)) for pair in zip(path, path[1:]))
        membership.append((line_id, segments))

    # Pull back the archived RGS geometry at the selected Q6 labels.
    seed = load_rgs_seed_data()  # verifies archived array hashes
    image_coordinates = tuple(
        tuple(float(x) for x in seed.coordinates[label]) for label in labels
    )

    return KatharaQ6ChainMapAudit(
        source.vertices, source.edges, source.faces,
        target.vertices, target.edges, target.faces,
        labels,
        tuple(tuple(row) for row in R0),
        tuple(tuple(row) for row in R1),
        tuple(tuple(row) for row in R2),
        identity0, identity1, len(set(labels)) == len(labels),
        len(source.faces), tuple(membership), image_coordinates,
        "VALID INTEGRAL CHAIN MAP; injective on vertices, not a cellular embedding; "
        "all source 2-cells collapse. Coordinate code is explicit but not RGS-canonical."
    )

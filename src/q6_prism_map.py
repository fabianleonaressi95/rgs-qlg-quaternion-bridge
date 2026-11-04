"""Explicit partial Q6-to-PRISM graph map.

The inverse-Gray-code map embeds a Hamiltonian cycle of Q6 into the C64
factor of PRISM, with all other PRISM coordinates fixed at zero. It is
injective and preserves exactly that 64-edge Hamiltonian cycle, not all
192 Q6 edges. This distinction is part of the result.
"""
from __future__ import annotations


def gray_code(index: int) -> int:
    if not 0 <= int(index) < 64:
        raise ValueError("Gray-code index must lie in 0..63.")
    return int(index) ^ (int(index) >> 1)


def inverse_gray_code(value: int) -> int:
    if not 0 <= int(value) < 64:
        raise ValueError("Gray-code value must lie in 0..63.")
    result = int(value)
    shift = result >> 1
    while shift:
        result ^= shift
        shift >>= 1
    return result


def q6_vertex_to_prism(vertex: int) -> tuple[int, int, int, int, int]:
    """Inject a Q6 label into PRISM coordinates (C64,C6,C6,C6,C5)."""
    if not 0 <= int(vertex) < 64:
        raise ValueError("Q6 vertex must lie in 0..63.")
    return (inverse_gray_code(int(vertex)), 0, 0, 0, 0)


def audit_q6_prism_map() -> dict[str, object]:
    vertices = range(64)
    image = [q6_vertex_to_prism(v) for v in vertices]
    mapped_cycle_edges = 0
    for i in range(64):
        j = (i + 1) % 64
        u, v = gray_code(i), gray_code(j)
        if (u ^ v).bit_count() == 1:
            mapped_cycle_edges += 1
    return {
        "audit": "Q6_TO_PRISM_PARTIAL_GRAPH_EMBEDDING",
        "map": "v -> (inverse_gray(v),0,0,0,0)",
        "injective": len(set(image)) == 64,
        "image_vertices": len(set(image)),
        "q6_vertices": 64,
        "q6_edges": 192,
        "preserved_hamiltonian_cycle_edges": mapped_cycle_edges,
        "q6_edge_preservation_fraction": mapped_cycle_edges / 192.0,
        "preserves_all_q6_adjacencies": mapped_cycle_edges == 192,
        "target_graph": "C64 x C6 x C6 x C6 x C5",
        "scope": "injective vertex map preserving a spanning 64-cycle; not a Q6 graph isomorphism or full edge embedding",
    }

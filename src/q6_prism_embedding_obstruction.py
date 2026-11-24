"""Exact graph-theoretic obstruction and minimal product extension for Q6.

The existing PRISM graph is C64 □ C6 □ C6 □ C6 □ C5 (five cycle
factors). Q6 cannot embed injectively as an edge-preserving subgraph into
this product. The proof is recorded in the accompanying research note.

A six-factor product does contain an induced Q6: use {0,1} in each of
six distinct cycle factors. A deliberately collapsed parity homomorphism
into the existing product is also provided as a baseline, not as a useful
geometric identification.
"""
from __future__ import annotations
from itertools import product
from typing import Sequence

EXISTING_PRISM_CYCLES = (64, 6, 6, 6, 5)
SIX_FACTOR_EXTENSION = (64, 6, 6, 6, 6, 5)


def q6_vertices() -> tuple[int, ...]:
    return tuple(range(64))


def q6_edges() -> tuple[tuple[int, int], ...]:
    return tuple((v, v ^ (1 << d)) for v in range(64)
                 for d in range(6) if not (v >> d) & 1)


def has_cycle_factor_embedding_obstruction(cube_dimension: int,
                                           factor_count: int) -> bool:
    """Necessary local condition: each pair of cube directions needs distinct factors."""
    if cube_dimension < 0 or factor_count < 0:
        raise ValueError("dimensions and factor counts must be nonnegative")
    return factor_count < cube_dimension


def embed_q6_in_six_factor_product(vertex: int,
                                   cycle_lengths: Sequence[int] = SIX_FACTOR_EXTENSION
                                   ) -> tuple[int, ...]:
    """Induced Q6 embedding: bit d selects 0 or 1 in factor d."""
    if not isinstance(vertex, int) or vertex not in range(64):
        raise ValueError("Q6 vertex must be an integer in 0..63")
    lengths = tuple(cycle_lengths)
    if len(lengths) < 6 or any(n < 5 for n in lengths):
        raise ValueError("need at least six cycle factors, each of length >=5")
    return tuple((vertex >> d) & 1 if d < 6 else 0 for d in range(len(lengths)))


def product_adjacent(a: Sequence[int], b: Sequence[int],
                     cycle_lengths: Sequence[int]) -> bool:
    """Adjacency in a Cartesian product of cycles."""
    lengths = tuple(cycle_lengths)
    if len(a) != len(lengths) or len(b) != len(lengths):
        raise ValueError("coordinate count must equal number of cycle factors")
    changed = [i for i, (x, y) in enumerate(zip(a, b)) if x != y]
    if len(changed) != 1:
        return False
    i = changed[0]
    return (a[i] - b[i]) % lengths[i] in (1, lengths[i] - 1)


def parity_collapse_to_existing_prism(vertex: int) -> tuple[int, ...]:
    """A valid but highly non-injective graph homomorphism Q6 -> existing PRISM."""
    if not isinstance(vertex, int) or vertex not in range(64):
        raise ValueError("Q6 vertex must be an integer in 0..63")
    return (vertex.bit_count() % 2, 0, 0, 0, 0)


def audit() -> dict:
    return {
        "existing_prism_cycle_factors": list(EXISTING_PRISM_CYCLES),
        "existing_factor_count": len(EXISTING_PRISM_CYCLES),
        "q6_dimension": 6,
        "injective_edge_preserving_embedding": False,
        "obstruction": "each pair of incident Q6 directions spans a square; in a product of cycles of length >=5 every simple 4-cycle is a rectangle using two distinct factors; therefore the six incident Q6 edges at a vertex require six distinct factors",
        "minimal_factor_count": 6,
        "six_factor_extension": list(SIX_FACTOR_EXTENSION),
        "six_factor_embedding_is_induced": True,
        "parity_homomorphism": {
            "exists": True,
            "image_vertices": 2,
            "fiber_sizes": [32, 32],
            "interpretation": "valid graph homomorphism but collapses all six Q6 directions to one edge; not an embedding",
        },
        "scope": "graph-theoretic result only; no physical or transport equivalence",
    }

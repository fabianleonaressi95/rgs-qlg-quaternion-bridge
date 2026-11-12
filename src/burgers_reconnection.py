"""Local Burgers-conserving reconnection for equal-Burgers strands.

The operation is deliberately narrow: two disjoint oriented edges with the
same Burgers vector exchange their heads.  It preserves the full C_0-valued
boundary exactly.  Unequal-Burgers reactions require a separately specified
junction reaction and are rejected here.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Hashable, Iterable

Vertex = Hashable
Burgers = tuple[int, ...]


@dataclass(frozen=True)
class OrientedBurgersEdge:
    tail: Vertex
    head: Vertex
    burgers: Burgers


def boundary_flux(
    edges: Iterable[OrientedBurgersEdge],
) -> dict[Vertex, Burgers]:
    """Compute ∂₁c as a vertex -> Burgers-vector map."""
    totals: dict[Vertex, list[int]] = defaultdict(list)
    dimension: int | None = None
    for edge in edges:
        if not edge.burgers:
            raise ValueError("Burgers vectors must have at least one component.")
        if dimension is None:
            dimension = len(edge.burgers)
        elif len(edge.burgers) != dimension:
            raise ValueError("All Burgers vectors must have the same dimension.")
        if edge.tail == edge.head:
            raise ValueError("Self-loop edges are not admissible in this local rule.")
        if edge.tail not in totals:
            totals[edge.tail] = [0] * len(edge.burgers)
        if edge.head not in totals:
            totals[edge.head] = [0] * len(edge.burgers)
        for i, value in enumerate(edge.burgers):
            totals[edge.tail][i] -= value
            totals[edge.head][i] += value
    return {vertex: tuple(vector) for vertex, vector in totals.items()}


def reconnect_equal_burgers(
    first: OrientedBurgersEdge,
    second: OrientedBurgersEdge,
) -> tuple[OrientedBurgersEdge, OrientedBurgersEdge]:
    """Exchange heads of two disjoint strands carrying the same Burgers vector.

    If first=(a,b,B) and second=(c,d,B), return (a,d,B),(c,b,B).
    The local boundary ∂₁c is unchanged exactly.  This is a graph-level
    bookkeeping rule, not a physical criterion for when reconnection occurs.
    """
    if len({first.tail, first.head, second.tail, second.head}) != 4:
        raise ValueError("The two strands must have four distinct endpoints.")
    if first.burgers != second.burgers:
        raise ValueError(
            "Head exchange preserves Burgers boundary only for equal Burgers vectors."
        )
    if not first.burgers:
        raise ValueError("Burgers vectors must have at least one component.")

    before = boundary_flux((first, second))
    reconnected = (
        OrientedBurgersEdge(first.tail, second.head, first.burgers),
        OrientedBurgersEdge(second.tail, first.head, second.burgers),
    )
    after = boundary_flux(reconnected)
    if before != after:  # defensive invariant check
        raise ArithmeticError("Internal error: reconnection changed Burgers boundary.")
    return reconnected


def junction_is_burgers_conservative(
    incoming: Iterable[Burgers],
    outgoing: Iterable[Burgers],
) -> bool:
    """Check Σ b_in = Σ b_out for a declared oriented junction."""
    incoming = tuple(tuple(v) for v in incoming)
    outgoing = tuple(tuple(v) for v in outgoing)
    all_vectors = incoming + outgoing
    if not all_vectors:
        return True
    dimension = len(all_vectors[0])
    if any(len(v) != dimension for v in all_vectors):
        raise ValueError("All Burgers vectors at a junction must share a dimension.")
    return tuple(sum(v[i] for v in incoming) for i in range(dimension)) == tuple(
        sum(v[i] for v in outgoing) for i in range(dimension)
    )

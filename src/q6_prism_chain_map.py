"""An explicit surjective cubical-chain map Q6 -> a Q5 corner of PRISM.

PRISM has cycle factors C64 x C6 x C6 x C6 x C5. Its binary corner
{0,1}^5 is an induced cubical Q5. The vertex map is
    (v0,v1,v2,v3,v4,v5) -> (v0 XOR v5,v1,v2,v3,v4).
It identifies the two Q6 directions 0 and 5. The chain map below sends
a cell containing both directions to zero; otherwise it sends it to the
corresponding target cell with the orientation sign forced by the fold.

This is an algebraic chain map into a contractible corner subcomplex,
not a homotopy equivalence onto the full PRISM torus.
"""
from __future__ import annotations
from itertools import product
from typing import Dict, Iterable, Tuple

Cell = Tuple[Tuple[int, ...], Tuple[int, ...]]
Chain = Dict[Cell, int]
PRISM_CYCLE_LENGTHS = (64, 6, 6, 6, 5)


def all_cells(dimension: int) -> Iterable[Cell]:
    """Enumerate cubical cells by states 0=fixed zero, 1=fixed one, 2=free."""
    for states in product((0, 1, 2), repeat=dimension):
        base = tuple(0 if s == 2 else s for s in states)
        directions = tuple(i for i, s in enumerate(states) if s == 2)
        yield base, directions


def add_term(chain: Chain, cell: Cell, coefficient: int) -> None:
    if coefficient:
        chain[cell] = chain.get(cell, 0) + coefficient
        if chain[cell] == 0:
            del chain[cell]


def boundary(cell: Cell) -> Chain:
    """Oriented cubical boundary for directions in increasing order."""
    base, directions = cell
    result: Chain = {}
    for j, d in enumerate(directions):
        remaining = tuple(x for x in directions if x != d)
        upper = list(base)
        upper[d] = 1
        add_term(result, (tuple(upper), remaining), (-1) ** j)
        add_term(result, (base, remaining), -((-1) ** j))
    return result


def vertex_map(vertex: Tuple[int, ...]) -> Tuple[int, ...]:
    if len(vertex) != 6 or any(x not in (0, 1) for x in vertex):
        raise ValueError("vertex must be a six-bit tuple")
    return (vertex[0] ^ vertex[5], vertex[1], vertex[2], vertex[3], vertex[4])


def prism_corner_vertex(vertex: Tuple[int, ...]) -> Tuple[int, ...]:
    """Place a Q5 vertex in PRISM factor coordinates (C64,C6,C6,C6,C5)."""
    return vertex_map(vertex)


def map_cell(cell: Cell) -> Tuple[Cell, int] | None:
    """Map one oriented Q6 cell to a signed Q5 cell, or to zero."""
    base, directions = cell
    if len(base) != 6 or any(x not in (0, 1) for x in base):
        raise ValueError("invalid Q6 cell base")
    if tuple(sorted(set(directions))) != directions or any(d not in range(6) for d in directions):
        raise ValueError("cell directions must be strictly increasing in 0..5")
    if any(base[d] != 0 for d in directions):
        raise ValueError("free coordinates must be zero in the canonical base")

    # Both source directions fold to the same target coordinate: this cell
    # loses one dimension and is sent to zero.
    if 0 in directions and 5 in directions:
        return None

    has_fold_direction = 0 in directions or 5 in directions
    fixed_fold_coordinate = base[0] ^ base[5]
    target_directions = tuple(sorted(
        ([0] if has_fold_direction else []) +
        [d for d in directions if 1 <= d <= 4]
    ))

    target_base = [0] * 5
    if 0 not in target_directions:
        target_base[0] = fixed_fold_coordinate
    for d in range(1, 5):
        if d not in target_directions:
            target_base[d] = base[d]

    # If the varying source coordinate is 0 or 5, its image is negatively
    # oriented exactly when the other folded bit is 1.
    sign = 1
    if 0 in directions:
        sign *= (-1) ** base[5]
    elif 5 in directions:
        sign *= (-1) ** base[0]
        # Source direction 5 moves to target direction 0, crossing all
        # remaining ordered directions in the orientation wedge.
        sign *= (-1) ** (len(directions) - 1)

    return (tuple(target_base), target_directions), sign


def map_chain(chain: Chain) -> Chain:
    result: Chain = {}
    for cell, coefficient in chain.items():
        mapped = map_cell(cell)
        if mapped is not None:
            target, sign = mapped
            add_term(result, target, coefficient * sign)
    return result


def chain_map_identity_holds(cell: Cell) -> bool:
    """Check ∂F(cell) = F∂(cell) exactly over integer coefficients."""
    mapped = map_cell(cell)
    lhs: Chain = {}
    if mapped is not None:
        target, sign = mapped
        for face, coefficient in boundary(target).items():
            add_term(lhs, face, sign * coefficient)
    rhs = map_chain(boundary(cell))
    return lhs == rhs


def audit() -> dict:
    cells = tuple(all_cells(6))
    failures = [cell for cell in cells if not chain_map_identity_holds(cell)]
    image_vertices = {
        vertex_map(tuple(bits))
        for bits in product((0, 1), repeat=6)
    }
    mapped_cells = {
        map_cell(cell)[0]
        for cell in cells
        if map_cell(cell) is not None
    }
    return {
        "source": "cubical Q6",
        "target": "binary-corner cubical Q5 inside C64 x C6 x C6 x C6 x C5",
        "vertex_map": "(v0,v1,v2,v3,v4,v5) -> (v0 XOR v5,v1,v2,v3,v4)",
        "source_cell_count": len(cells),
        "target_cell_count": 3 ** 5,
        "source_vertices": 64,
        "target_vertices_reached": len(image_vertices),
        "vertex_map_fiber_sizes": [2] * len(image_vertices),
        "source_edges": 192,
        "target_edges_reached": sum(1 for c in mapped_cells if len(c[1]) == 1),
        "source_cells_collapsed": sum(1 for c in cells if map_cell(c) is None),
        "chain_map_identity_checked_on_all_cells": len(failures) == 0,
        "chain_map_failures": len(failures),
        "surjective_on_target_cells": len(mapped_cells) == 3 ** 5,
        "induced_homology": {
            "H0": "isomorphism (both cubical cubes are connected; vertex map is surjective)",
            "higher_degrees": "zero on homology because source and target cubes are contractible",
        },
        "scope": "explicit integral cubical chain map; not a homotopy equivalence to full PRISM and not a transport theorem",
    }

"""Path-valued graph map Q6^(1) -> PRISM^(1) inducing Z^129 -> Z^5 surjectively.

This is an explicit existence construction, not a canonical geometric model:
all 64 source vertices map to the PRISM basepoint; tree edges map to
constant paths; five chosen non-tree edges map to one full coordinate
cycle each; all remaining non-tree edges map to constant paths.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, List, Tuple

LENGTHS = (64, 6, 6, 6, 5)
Vertex = int
Edge = Tuple[int, int]


def q6_edges() -> List[Edge]:
    return sorted(
        (v, v ^ (1 << d))
        for v in range(64)
        for d in range(6)
        if ((v >> d) & 1) == 0
    )


def parent(v: Vertex) -> Vertex:
    if v == 0:
        raise ValueError("root has no parent")
    return v ^ (1 << (v.bit_length() - 1))


def spanning_tree_edges() -> set[Edge]:
    return {tuple(sorted((v, parent(v)))) for v in range(1, 64)}


def fundamental_chords() -> List[Edge]:
    tree = spanning_tree_edges()
    return [e for e in q6_edges() if e not in tree]


def factor_loop(factor: int) -> List[Tuple[int, ...]]:
    """Positive unit steps around one PRISM factor, based at the origin."""
    if factor not in range(5):
        raise ValueError("factor must be in 0..4")
    path = []
    for k in range(LENGTHS[factor] + 1):
        point = [0] * 5
        point[factor] = k % LENGTHS[factor]
        path.append(tuple(point))
    return path


@dataclass(frozen=True)
class EdgeImage:
    source_edge: Edge
    factor: int | None
    path: Tuple[Tuple[int, ...], ...]


def edge_images() -> Dict[Edge, EdgeImage]:
    chords = fundamental_chords()
    selected = {chords[i]: i for i in range(5)}
    result = {}
    base = (0, 0, 0, 0, 0)
    for edge in q6_edges():
        if edge in selected:
            factor = selected[edge]
            path = tuple(factor_loop(factor))
        else:
            factor = None
            path = (base,)
        result[edge] = EdgeImage(edge, factor, path)
    return result


def fundamental_cycle_image(chord: Edge) -> Tuple[int, ...]:
    """Image of the fundamental cycle: chord plus its tree return path."""
    chords = fundamental_chords()
    if chord not in chords:
        raise ValueError("edge is not a non-tree chord")
    selected = {chords[i]: i for i in range(5)}
    if chord in selected:
        out = [0] * 5
        out[selected[chord]] = 1
        return tuple(out)
    return (0, 0, 0, 0, 0)


def audit() -> dict:
    edges = q6_edges()
    tree = spanning_tree_edges()
    chords = fundamental_chords()
    images = edge_images()
    matrix_columns = [fundamental_cycle_image(c) for c in chords]
    nonzero = sum(any(x) for x in matrix_columns)
    return {
        "source_graph": "Q6 1-skeleton only",
        "target": "1-skeleton of C64 x C6 x C6 x C6 x C5",
        "source_vertices": 64,
        "source_edges": len(edges),
        "spanning_tree_edges": len(tree),
        "fundamental_chords": len(chords),
        "source_graph_cycle_rank": len(edges) - 64 + 1,
        "target_periodic_H1_rank": 5,
        "all_source_vertices_map_to": [0, 0, 0, 0, 0],
        "selected_chords": [list(chords[i]) for i in range(5)],
        "selected_chord_images": [list(fundamental_cycle_image(chords[i])) for i in range(5)],
        "induced_map_rank": 5,
        "induced_map_is_surjective": True,
        "nonzero_fundamental_cycle_images": nonzero,
        "target_factor_loop_lengths": list(LENGTHS),
        "total_target_steps_on_selected_edges": sum(LENGTHS),
        "all_other_source_edges_map_to_constant_paths": True,
        "extends_over_all_240_Q6_plaquettes": False,
        "interpretation": "constructive path-valued graph map; existence proof only, not a canonical Q6/PRISM geometric or dynamical correspondence",
    }


if __name__ == "__main__":
    import json
    print(json.dumps(audit(), indent=2))

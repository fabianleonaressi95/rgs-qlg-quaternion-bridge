from itertools import product
from src.q6_prism_chain_map import (
    all_cells, audit, boundary, chain_map_identity_holds, map_cell,
    map_chain, vertex_map, PRISM_CYCLE_LENGTHS,
)

def test_vertex_map_is_surjective_two_to_one():
    images = [vertex_map(v) for v in product((0, 1), repeat=6)]
    assert len(set(images)) == 32
    assert sorted(images.count(x) for x in set(images)) == [2] * 32

def test_each_q6_edge_maps_to_a_prism_corner_edge():
    edges = [c for c in all_cells(6) if len(c[1]) == 1]
    assert len(edges) == 192
    mapped = [map_cell(c) for c in edges]
    assert all(m is not None and len(m[0][1]) == 1 for m in mapped)
    assert len({m[0] for m in mapped}) == 80

def test_every_cell_satisfies_chain_map_identity_over_integers():
    cells = tuple(all_cells(6))
    assert len(cells) == 3 ** 6 == 729
    assert all(chain_map_identity_holds(c) for c in cells)

def test_map_is_surjective_on_all_target_cells():
    images = {map_cell(c)[0] for c in all_cells(6) if map_cell(c) is not None}
    assert len(images) == 3 ** 5 == 243

def test_cells_with_both_folded_directions_are_collapsed():
    for base in product((0, 1), repeat=6):
        cell = (tuple(0 if d in (0, 5) else base[d] for d in range(6)), (0, 5))
        assert map_cell(cell) is None

def test_audit_and_prism_factor_order():
    result = audit()
    assert PRISM_CYCLE_LENGTHS == (64, 6, 6, 6, 5)
    assert result["target_vertices_reached"] == 32
    assert result["chain_map_identity_checked_on_all_cells"]
    assert result["chain_map_failures"] == 0
    assert result["surjective_on_target_cells"]

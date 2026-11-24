from src.q6_prism_path_valued_homology import (
    LENGTHS, audit, edge_images, fundamental_chords,
    fundamental_cycle_image, q6_edges, spanning_tree_edges,
)

def test_q6_graph_counts_and_fundamental_basis():
    assert len(q6_edges()) == 192
    assert len(spanning_tree_edges()) == 63
    assert len(fundamental_chords()) == 129

def test_five_selected_fundamental_cycles_map_to_standard_basis():
    chords = fundamental_chords()
    for i in range(5):
        expected = tuple(1 if j == i else 0 for j in range(5))
        assert fundamental_cycle_image(chords[i]) == expected

def test_remaining_fundamental_cycles_map_to_zero():
    chords = fundamental_chords()
    assert all(fundamental_cycle_image(c) == (0, 0, 0, 0, 0)
               for c in chords[5:])

def test_each_source_edge_has_a_based_target_path():
    images = edge_images()
    assert len(images) == 192
    for edge_image in images.values():
        assert edge_image.path[0] == (0, 0, 0, 0, 0)
        assert edge_image.path[-1] == (0, 0, 0, 0, 0)

def test_selected_paths_are_exact_factor_cycles():
    chords = fundamental_chords()
    images = edge_images()
    for i in range(5):
        path = images[chords[i]].path
        assert len(path) - 1 == LENGTHS[i]
        assert len(set(path[:-1])) == LENGTHS[i]

def test_audit_rank_and_scope():
    result = audit()
    assert result["induced_map_rank"] == 5
    assert result["induced_map_is_surjective"]
    assert result["extends_over_all_240_Q6_plaquettes"] is False

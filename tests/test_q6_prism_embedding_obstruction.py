from src.q6_prism_embedding_obstruction import (
    EXISTING_PRISM_CYCLES, SIX_FACTOR_EXTENSION, audit, q6_edges,
    embed_q6_in_six_factor_product, product_adjacent,
    parity_collapse_to_existing_prism, has_cycle_factor_embedding_obstruction,
)

def test_existing_prism_has_five_factors_and_cannot_host_q6_embedding():
    assert len(EXISTING_PRISM_CYCLES) == 5
    assert has_cycle_factor_embedding_obstruction(6, 5)
    assert not has_cycle_factor_embedding_obstruction(6, 6)

def test_six_factor_extension_contains_induced_q6():
    images = {v: embed_q6_in_six_factor_product(v) for v in range(64)}
    assert len(set(images.values())) == 64
    for v, w in q6_edges():
        assert product_adjacent(images[v], images[w], SIX_FACTOR_EXTENSION)
    # No non-edge becomes an edge: distinct cube vertices differ in >=1
    # selected binary coordinates, and a product edge changes exactly one.
    for v in range(64):
        for w in range(v + 1, 64):
            if (v ^ w).bit_count() != 1:
                assert not product_adjacent(images[v], images[w], SIX_FACTOR_EXTENSION)

def test_parity_collapse_is_a_homomorphism_but_not_an_embedding():
    images = {v: parity_collapse_to_existing_prism(v) for v in range(64)}
    assert len(set(images.values())) == 2
    assert sorted(sum(1 for x in images.values() if x == p)
                  for p in set(images.values())) == [32, 32]
    for v, w in q6_edges():
        assert product_adjacent(images[v], images[w], EXISTING_PRISM_CYCLES)

def test_audit_states_exact_scope():
    result = audit()
    assert result["injective_edge_preserving_embedding"] is False
    assert result["minimal_factor_count"] == 6
    assert result["six_factor_embedding_is_induced"] is True
    assert result["parity_homomorphism"]["image_vertices"] == 2

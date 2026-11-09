from src.kathara_q6_chain_map import audit_kathara_to_q6_chain_map


def test_kathara_q6_chain_identities_are_exact():
    a = audit_kathara_to_q6_chain_map()
    assert a.chain_identity_0
    assert a.chain_identity_1


def test_map_dimensions_and_vertex_injectivity():
    a = audit_kathara_to_q6_chain_map()
    assert len(a.R0) == 64 and all(len(row) == 12 for row in a.R0)
    assert len(a.R1) == 192 and all(len(row) == 21 for row in a.R1)
    assert len(a.R2) == 240 and all(len(row) == 10 for row in a.R2)
    assert a.injective_on_vertices
    assert len(set(a.vertex_labels)) == 12


def test_map_does_not_claim_to_embed_kathara_2_cells():
    a = audit_kathara_to_q6_chain_map()
    assert a.collapsed_source_2_cells == 10
    assert "not a cellular embedding" in a.scope
    assert "not RGS-canonical" in a.scope


def test_fifteen_primary_lines_are_retained_as_metadata():
    a = audit_kathara_to_q6_chain_map()
    assert len(a.line_edge_membership) == 15

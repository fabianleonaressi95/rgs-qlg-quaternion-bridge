from src.kathara_q3_chain_map import audit_kathara_to_q3_chain_map


def test_chain_map_commutes_with_both_boundary_operators():
    a = audit_kathara_to_q3_chain_map()
    assert a.chain_identity_0
    assert a.chain_identity_1


def test_map_is_integral_and_has_expected_shapes():
    a = audit_kathara_to_q3_chain_map()
    assert len(a.S0) == 8 and all(len(r) == 12 for r in a.S0)
    assert len(a.S1) == 12 and all(len(r) == 21 for r in a.S1)
    assert len(a.S2) == 6 and all(len(r) == 10 for r in a.S2)
    assert all(x in (-1, 0, 1) for row in a.S0 for x in row)


def test_projection_is_explicitly_coarse_and_uses_one_square():
    a = audit_kathara_to_q3_chain_map()
    assert set(a.source_vertex_map) <= {0, 1, 2, 3}
    assert a.image_uses_one_square_face

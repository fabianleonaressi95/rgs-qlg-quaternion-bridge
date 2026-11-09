from src.kathara_cellular_complex import audit_kathara_cellular_complex

def test_chain_dimensions_and_boundary_squared_zero():
    a = audit_kathara_cellular_complex()
    assert (len(a.vertices), len(a.edges), len(a.faces)) == (12, 21, 10)
    assert a.b1b2_zero

def test_exact_ranks_betti_and_euler_characteristic():
    a = audit_kathara_cellular_complex()
    assert (a.rank_b1, a.rank_b2) == (11, 10)
    assert a.betti == (1, 0, 0)
    assert a.euler_characteristic == 1

def test_hodge_laplacian_dimensions():
    a = audit_kathara_cellular_complex()
    assert len(a.laplacian_0) == 12 and all(len(r) == 12 for r in a.laplacian_0)
    assert len(a.laplacian_1) == 21 and all(len(r) == 21 for r in a.laplacian_1)
    assert len(a.laplacian_2) == 10 and all(len(r) == 10 for r in a.laplacian_2)

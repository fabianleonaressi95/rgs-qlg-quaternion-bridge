import numpy as np
import pytest
from src.geb_formal_layer import (
    audit_geb, bach_spectral_audit, bach_substitution_matrix,
    connected_double_cover_criterion, escher_sheet_involution,
    godel_finite_encoding, strange_loop_residual, torus_audit,
)

def test_bach_supergolden_matrix_identity():
    a=bach_substitution_matrix()
    assert np.array_equal(np.linalg.matrix_power(a,3)-np.linalg.matrix_power(a,2),np.eye(3,dtype=int))
    r=bach_spectral_audit()
    assert r["spectral_radius_residual"]<1e-12
    assert r["matrix_polynomial_residual"]<1e-12

def test_escher_triangulated_torus_chain_complex():
    r=torus_audit(4)
    assert (r["V"],r["E"],r["F"],r["euler_characteristic"])==(16,48,32,0)
    assert r["betti_numbers_from_ranks"]==[1,2,1]
    assert r["boundary_squared_residual"]==0
    assert r["period_matrix"]==[[1.0,0.0],[0.0,1.0]]
    assert max(r["closed_cochain_residuals"])<1e-12
    assert r["hodge_zero_mode_counts"]==[1,2,1]

def test_escher_involution():
    for s in (0,1):
        assert escher_sheet_involution(*escher_sheet_involution(3,2,s))==(3,2,s)

def test_connected_cover_requires_nonzero_class():
    assert not connected_double_cover_criterion((0,0))["connected"]
    assert connected_double_cover_criterion((1,0))["connected"]

def test_godel_hash_is_canonical_but_not_incompleteness():
    a=godel_finite_encoding({"b":2,"a":1})
    b=godel_finite_encoding({"a":1,"b":2})
    assert a["sha256"]==b["sha256"]
    assert "no incompleteness theorem" in a["interpretation"]

def test_strange_loop_residual():
    assert strange_loop_residual(np.array([1.,2.]),lambda x:x)==0.
    assert strange_loop_residual(np.array([1.,2.]),lambda x:.5*x)>0.

def test_audit_status_and_open_proofs():
    r=audit_geb()
    assert r["status"].startswith("finite combinatorial")
    assert len(r["open_obligations"])>=6

def test_invalid_cover_class():
    with pytest.raises(ValueError): connected_double_cover_criterion((0,2))

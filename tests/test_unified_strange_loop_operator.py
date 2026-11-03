import numpy as np
from src.unified_strange_loop_operator import (
    audit_unified_operator, build_strange_loop_operator,
    escher_lifted_torus_laplacian, torus_q6_embedding,
)

def test_operator_is_symmetric_psd_and_constant_mode_kernel():
    H,_=build_strange_loop_operator()
    assert H.shape==(384,384)
    assert np.linalg.norm(H-H.T)<1e-12
    assert np.linalg.norm(H@np.ones(384))<1e-11
    assert np.linalg.eigvalsh(H)[0]>-1e-10

def test_torus_gray_embedding_and_lift():
    ids=torus_q6_embedding()
    assert len(ids)==16 and len(set(ids.tolist()))==16
    L=escher_lifted_torus_laplacian()
    assert np.linalg.norm(L-L.T)<1e-12
    assert np.linalg.norm(L@np.ones(64))<1e-12

def test_unified_audit_and_prism_constant_intertwiner():
    r=audit_unified_operator()
    assert r["dimension"]==384
    assert r["constant_mode_residual"]<1e-10
    assert r["prism_intertwiner"]["rank"]==1
    assert r["prism_intertwiner"]["intertwining_residual"]<1e-10
    assert "only" in r["prism_intertwiner"]["interpretation"]

def test_invalid_weights():
    import pytest
    with pytest.raises(ValueError):
        build_strange_loop_operator(kappa=-1)

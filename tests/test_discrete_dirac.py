import numpy as np
import pytest
from src.discrete_dirac import (assemble_spatial_dirac,curved_clifford_audit,curved_gamma_matrices,discrete_dirac_spectrum,polar_spin_links,static_tetrad_4d)
from src.clifford_qlg import gamma_matrices

def test_static_tetrad_and_curved_clifford_relations():
    rng=np.random.default_rng(12); e=[]
    for _ in range(5):
        a=rng.normal(size=(3,3)); e.append(a.T@a+np.eye(3))
    e=np.asarray(e); tetrad=static_tetrad_4d(e)
    assert tetrad.shape==(5,4,4) and np.allclose(tetrad[:,0,0],1)
    audit=curved_clifford_audit(e); assert audit["passed"] and audit["max_clifford_residual"]<2e-10

def test_curved_gamma_uses_inverse_coframe():
    e=np.broadcast_to(np.diag([2.,3.,4.]),(2,3,3)).copy(); curved=curved_gamma_matrices(e); gs=gamma_matrices()
    for mu,scale in enumerate([2.,3.,4.]): assert np.allclose(curved[0,mu+1],gs[mu+1]/scale)

def test_polar_spin_links_are_unitary():
    e=np.array([np.eye(3),np.diag([1.2,.9,1.1])]); links,res=polar_spin_links(e,[(0,1)])
    assert res<1e-12; s=links[(0,1)]
    assert s.shape==(2,2) and np.allclose(s.conj().T@s,np.eye(2),atol=1e-12) and np.isclose(np.linalg.det(s),1,atol=1e-12)

def test_discrete_spatial_dirac_is_hermitian_and_has_real_spectrum():
    x=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
    e=np.broadcast_to(np.eye(3),(4,3,3)).copy(); edges=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
    d,audit=assemble_spatial_dirac(x,e,edges)
    assert audit["hermiticity_residual"]<1e-12
    spec=discrete_dirac_spectrum(d); assert len(spec["eigenvalues"])==16 and np.all(np.isfinite(spec["eigenvalues"]))

def test_dirac_operator_scales_under_coordinate_dilation():
    x=np.array([[0.,0.,0.],[1.,0.,0.]]); e=np.broadcast_to(np.eye(3),(2,3,3)).copy()
    d1,_=assemble_spatial_dirac(x,e,[(0,1)]); d2,_=assemble_spatial_dirac(2*x,e,[(0,1)])
    assert np.allclose(d2,d1/2)

def test_nonpositive_vierbein_rejected():
    x=np.array([[0.,0.,0.],[1.,0.,0.]]); e=np.array([np.eye(3),np.diag([-1.,1.,1.])])
    with pytest.raises(ValueError): assemble_spatial_dirac(x,e,[(0,1)])

def test_curved_clifford_for_nonsymmetric_coframe():
    e=np.array([[[1.2,.2,.1],[.0,1.1,.3],[.1,.0,.9]]])
    assert np.linalg.det(e[0])>0
    audit=curved_clifford_audit(e)
    assert audit["passed"]
    assert audit["max_clifford_residual"]<2e-10

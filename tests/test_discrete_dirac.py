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

def test_compact_kernel_field_transform_is_gauge_covariant():
    from src.discrete_dirac import compact_kernel_transform
    rng=np.random.default_rng(41)
    def unitary():
        a=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
        q,r=np.linalg.qr(a)
        return q
    f=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
    p=unitary(); k=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
    g0,g1=unitary(),unitary()
    kernels={(0,1):k}; transports={(0,1):p}
    y,support=compact_kernel_transform(f,[(0,1)],kernels,transports)
    fp=np.stack([g0@f[0],g1@f[1]])
    pp={(0,1):g1@p@g0.conj().T}
    kp={(0,1):g1@k@g1.conj().T}
    yp,supportp=compact_kernel_transform(fp,[(0,1)],kp,pp)
    assert support==supportp=={(0,1)}
    assert np.allclose(yp[1],g1@y[1],atol=2e-12)


def _su2_phase(theta):
    return np.array([[np.cos(theta), 1j*np.sin(theta)],
                     [1j*np.sin(theta), np.cos(theta)]], dtype=complex)


def test_wilson_laplacian_is_hermitian_positive_and_scale_covariant():
    from src.discrete_dirac import assemble_wilson_laplacian
    x=np.array([[0.,0.,0.],[1.,0.,0.],[1.,1.,0.]])
    e=np.broadcast_to(np.eye(3),(3,3,3)).copy()
    edges=[(0,1),(1,2),(0,2)]
    w,a=assemble_wilson_laplacian(x,e,edges,wilson_r=1.25)
    assert a["hermiticity_residual"]<1e-12
    assert a["minimum_eigenvalue"]>-1e-12
    w2,_=assemble_wilson_laplacian(3*x,e,edges,wilson_r=1.25)
    assert np.allclose(w2,w/3,atol=1e-12)


def test_wilson_laplacian_is_locally_gauge_covariant():
    from src.discrete_dirac import assemble_wilson_laplacian
    x=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]])
    e=np.broadcast_to(np.eye(3),(3,3,3)).copy()
    edges=[(0,1),(0,2),(1,2)]
    links={(0,1):_su2_phase(.17),(0,2):_su2_phase(-.31),(1,2):_su2_phase(.43)}
    gauges=[_su2_phase(.11),_su2_phase(-.29),_su2_phase(.37)]
    links_g={(i,j):gauges[j]@links[(i,j)]@gauges[i].conj().T for i,j in edges}
    l,_=assemble_wilson_laplacian(x,e,edges,spin_links=links,wilson_r=.8)
    lg,_=assemble_wilson_laplacian(x,e,edges,spin_links=links_g,wilson_r=.8)
    big_g=np.zeros((12,12),dtype=complex)
    for i,g in enumerate(gauges):
        g4=np.block([[g,np.zeros((2,2))],[np.zeros((2,2)),g]])
        big_g[4*i:4*i+4,4*i:4*i+4]=g4
    assert np.allclose(lg,big_g@l@big_g.conj().T,atol=2e-12)


def test_wilson_term_lifts_pi_mode_on_periodic_1d_ring():
    from src.discrete_dirac import assemble_wilson_laplacian
    n=16
    x=np.column_stack([np.arange(n,dtype=float),np.zeros(n),np.zeros(n)])
    e=np.broadcast_to(np.eye(3),(n,3,3)).copy()
    edges=[(i,i+1) for i in range(n-1)]+[(0,n-1)]
    lengths={(i,j):1.0 for i,j in edges}
    w,a=assemble_wilson_laplacian(x,e,edges,wilson_r=1.0,edge_lengths=lengths)
    # A normalized plane wave on the periodic ring is an exact eigenvector.
    def rayleigh(k):
        phase=np.exp(2j*np.pi*k*np.arange(n)/n)
        v=np.zeros(4*n,dtype=complex)
        v[0::4]=phase/np.sqrt(n)
        return float(np.real(np.vdot(v,w@v)))
    assert abs(rayleigh(0))<1e-12
    assert np.isclose(rayleigh(n//2),2.0,atol=1e-12)
    assert a["minimum_eigenvalue"]>-1e-12


def test_wilson_term_breaks_exact_chiral_anticommutation():
    from src.discrete_dirac import assemble_spatial_dirac, chirality_matrix
    x=np.array([[0.,0.,0.],[1.,0.,0.]])
    e=np.broadcast_to(np.eye(3),(2,3,3)).copy()
    edges=[(0,1)]
    d0,_=assemble_spatial_dirac(x,e,edges)
    dw,a=assemble_spatial_dirac(x,e,edges,wilson_r=1.0)
    g5=chirality_matrix()
    assert np.allclose(g5,g5.conj().T,atol=1e-12)
    assert np.allclose(g5@g5,np.eye(4),atol=1e-12)
    gamma5=np.kron(np.eye(2),g5)
    naive=float(np.linalg.norm(d0@gamma5+gamma5@d0,ord="fro"))
    wilson=float(np.linalg.norm(dw@gamma5+gamma5@dw,ord="fro"))
    assert naive<1e-12
    assert wilson>1e-6
    assert a["hermiticity_residual"]<1e-12
    assert a["wilson_audit"]["minimum_eigenvalue"]>-1e-12


def test_assemble_spatial_dirac_rejects_negative_wilson_parameter():
    x=np.array([[0.,0.,0.],[1.,0.,0.]])
    e=np.broadcast_to(np.eye(3),(2,3,3)).copy()
    with pytest.raises(ValueError):
        assemble_spatial_dirac(x,e,[(0,1)],wilson_r=-1.0)


def test_wilson_laplacian_rejects_nonunitary_spin_links():
    from src.discrete_dirac import assemble_wilson_laplacian
    x=np.array([[0.,0.,0.],[1.,0.,0.]])
    e=np.broadcast_to(np.eye(3),(2,3,3)).copy()
    with pytest.raises(ValueError,match="unitary|SU\\(2\\)"):
        assemble_wilson_laplacian(x,e,[(0,1)],spin_links={(0,1):2*np.eye(2)})

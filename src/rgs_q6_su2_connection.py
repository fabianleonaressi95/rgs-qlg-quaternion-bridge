"""SU(2) parallel transport and a gauge-covariant Q6 Dirac-Kaehler operator.

Links live on canonically oriented Q6 edges. For e=(a,b), a<b, U_e maps
the fibre at b to the fibre at a. Under vertex gauge transformations
G_v, U_ab -> G_a U_ab G_b^dagger. The covariant coboundary is
(D_U psi)_e = U_ab psi_b - psi_a.

The resulting block operator [[0,D_U^dagger],[D_U,0]] is Hermitian and
gauge covariant. It is a finite graph Dirac-Kaehler operator, not a
Clifford-spin Dirac operator or a continuum fermion theory.
"""
from __future__ import annotations
import numpy as np
from src.rgs_q6_equivariant_cochains import build_q6_cochain_complex


def _su2_from_quaternion(q):
    q=np.asarray(q,dtype=float)
    if q.shape[-1]!=4: raise ValueError("quaternions must have four real components")
    norm=np.linalg.norm(q,axis=-1,keepdims=True)
    if np.any(norm<=1e-14) or not np.all(np.isfinite(q)): raise ValueError("invalid quaternion")
    a,b,c,d=(q/norm[...,None] if norm.ndim==q.ndim-1 else q/norm).T if q.ndim==2 else (q/norm).T
    out=np.empty((len(q),2,2),dtype=complex)
    out[:,0,0]=a+1j*b; out[:,0,1]=c+1j*d
    out[:,1,0]=-c+1j*d; out[:,1,1]=a-1j*b
    return out


def random_su2_matrices(count, seed=0):
    """Deterministic Haar-distributed SU(2) samples via normalized Gaussians."""
    if not isinstance(count,(int,np.integer)) or count<1: raise ValueError("count must be positive")
    rng=np.random.default_rng(seed)
    q=rng.normal(size=(int(count),4))
    q/=np.linalg.norm(q,axis=1,keepdims=True)
    out=np.empty((int(count),2,2),dtype=complex)
    a,b,c,d=q.T
    out[:,0,0]=a+1j*b; out[:,0,1]=c+1j*d
    out[:,1,0]=-c+1j*d; out[:,1,1]=a-1j*b
    return out


def validate_su2_field(matrices, count=None, atol=1e-10):
    u=np.asarray(matrices,dtype=complex)
    if u.ndim!=3 or u.shape[1:]!=(2,2) or (count is not None and len(u)!=count):
        raise ValueError("SU(2) field must have shape (N,2,2)")
    eye=np.eye(2)
    unit=max((float(np.linalg.norm(x.conj().T@x-eye)) for x in u),default=0.0)
    det=max((float(abs(np.linalg.det(x)-1.0)) for x in u),default=0.0)
    if not np.all(np.isfinite(u)) or unit>atol or det>atol:
        raise ValueError("matrices must lie in SU(2)")
    return {"unitarity_residual":unit,"determinant_residual":det}


def identity_link_field(complex_data=None):
    c=build_q6_cochain_complex() if complex_data is None else complex_data
    return np.repeat(np.eye(2,dtype=complex)[None,:,:],len(c["edges"]),axis=0)


def covariant_coboundary(links, complex_data=None):
    """Build D_U: C0(Q6,C2)->C1(Q6,C2), shape (384,128)."""
    c=build_q6_cochain_complex() if complex_data is None else complex_data
    u=np.asarray(links,dtype=complex)
    validate_su2_field(u,len(c["edges"]))
    d=np.zeros((2*len(c["edges"]),2*64),dtype=complex)
    eye=np.eye(2,dtype=complex)
    for e,(a,b,_direction) in enumerate(c["edges"]):
        rows=slice(2*e,2*e+2)
        d[rows,2*a:2*a+2]=-eye
        d[rows,2*b:2*b+2]=u[e]
    return d


def dirac_kaehler_operator(links, complex_data=None):
    """Hermitian off-diagonal graph Dirac operator on vertex plus edge spinors."""
    d=covariant_coboundary(links,complex_data)
    out=np.zeros((d.shape[1]+d.shape[0],)*2,dtype=complex)
    n=d.shape[1]
    out[:n,n:]=d.conj().T
    out[n:,:n]=d
    return out


def gauge_transform_links(links, gauge, complex_data=None):
    """Apply U_ab -> G_a U_ab G_b^dagger."""
    c=build_q6_cochain_complex() if complex_data is None else complex_data
    u=np.asarray(links,dtype=complex); g=np.asarray(gauge,dtype=complex)
    validate_su2_field(u,len(c["edges"])); validate_su2_field(g,64)
    transformed=np.empty_like(u)
    for e,(a,b,_d) in enumerate(c["edges"]):
        transformed[e]=g[a]@u[e]@g[b].conj().T
    return transformed


def gauge_covariance_matrices(gauge, complex_data=None):
    c=build_q6_cochain_complex() if complex_data is None else complex_data
    g=np.asarray(gauge,dtype=complex); validate_su2_field(g,64)
    p0=np.zeros((128,128),dtype=complex)
    p1=np.zeros((384,384),dtype=complex)
    for v in range(64): p0[2*v:2*v+2,2*v:2*v+2]=g[v]
    for e,(tail,_head,_d) in enumerate(c["edges"]):
        p1[2*e:2*e+2,2*e:2*e+2]=g[tail]
    return p0,p1


def plaquette_holonomies(links, complex_data=None):
    """Return oriented SU(2) holonomy for each elementary square."""
    c=build_q6_cochain_complex() if complex_data is None else complex_data
    u=np.asarray(links,dtype=complex); validate_su2_field(u,len(c["edges"]))
    out=[]
    for _v,_i,_j,cycle in c["faces"]:
        h=np.eye(2,dtype=complex)
        for a,b in zip(cycle,cycle[1:]+cycle[:1]):
            idx=c["edge_index"][(a,b)] if a<b else c["edge_index"][(b,a)]
            step=u[idx] if a<b else u[idx].conj().T
            h=h@step
        out.append(h)
    return np.asarray(out)


def audit_q6_su2_connection(seed=20261104):
    """Audit SU(2), gauge covariance, Hermiticity and Wilson-loop invariance."""
    c=build_q6_cochain_complex()
    links=random_su2_matrices(192,seed=seed)
    gauge=random_su2_matrices(64,seed=seed+1)
    link_report=validate_su2_field(links,192)
    gauge_report=validate_su2_field(gauge,64)
    d=covariant_coboundary(links,c)
    dk=dirac_kaehler_operator(links,c)
    links_g=gauge_transform_links(links,gauge,c)
    dg=covariant_coboundary(links_g,c)
    dkg=dirac_kaehler_operator(links_g,c)
    p0,p1=gauge_covariance_matrices(gauge,c)
    cov_d=float(np.linalg.norm(dg-p1@d@p0.conj().T))
    p=np.zeros((512,512),dtype=complex); p[:128,:128]=p0; p[128:,128:]=p1
    cov_dk=float(np.linalg.norm(dkg-p@dk@p.conj().T))
    herm=float(np.linalg.norm(dk-dk.conj().T))
    eig=np.linalg.eigvalsh(dk)
    min_eig=float(eig[0])
    h=plaquette_holonomies(links,c)
    hg=plaquette_holonomies(links_g,c)
    traces=np.trace(h,axis1=1,axis2=2)
    traces_g=np.trace(hg,axis1=1,axis2=2)
    trace_res=float(np.max(np.abs(traces-traces_g)))
    # Holonomy transforms by conjugation at the first vertex of each oriented face.
    conj_res=0.0
    for k,(_v,_i,_j,cycle) in enumerate(c["faces"]):
        expected=gauge[cycle[0]]@h[k]@gauge[cycle[0]].conj().T
        conj_res=max(conj_res,float(np.linalg.norm(hg[k]-expected)))
    return {
      "vertices":64,"edges":192,"plaquettes":240,"fiber_dimension":2,
      "link_unitarity_residual":link_report["unitarity_residual"],
      "link_determinant_residual":link_report["determinant_residual"],
      "gauge_unitarity_residual":gauge_report["unitarity_residual"],
      "covariant_coboundary_shape":list(d.shape),
      "dirac_kaehler_shape":list(dk.shape),
      "dirac_kaehler_hermiticity_residual":herm,
      "dirac_kaehler_min_eigenvalue":min_eig,
      "gauge_covariance_DU_residual":cov_d,
      "gauge_covariance_DK_residual":cov_dk,
      "wilson_trace_gauge_invariance_residual":trace_res,
      "plaquette_conjugation_residual":conj_res,
      "mean_wilson_trace_real":float(np.mean(traces.real)),
      "all_checks_passed":bool(max(link_report["unitarity_residual"],link_report["determinant_residual"],
        gauge_report["unitarity_residual"],gauge_report["determinant_residual"],herm,cov_d,cov_dk,
        trace_res,conj_res)<1e-10 and min_eig>-1e-10)
    }

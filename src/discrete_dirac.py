"""Weiler-compatible vierbein and discrete spatial Dirac construction.

Finite-graph scope only: no continuum limit or homogeneous-space assumption.
"""
from __future__ import annotations
from collections.abc import Mapping, Sequence
import numpy as np
from src.clifford_qlg import ETA, gamma_matrices
from src.representations_2i import quaternion_to_fundamental_su2

def _validate_geometry(coordinates, vierbein):
    x=np.asarray(coordinates,dtype=float); e=np.asarray(vierbein,dtype=float)
    if x.ndim!=2 or x.shape[1]!=3: raise ValueError("coordinates must have shape (N,3)")
    if e.shape!=(len(x),3,3): raise ValueError("vierbein must have shape (N,3,3)")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(e)): raise ValueError("geometry must be finite")
    for i,ei in enumerate(e):
        if np.linalg.det(ei)<=0 or np.linalg.cond(ei)>1e12: raise ValueError(f"vierbein[{i}] must be invertible and orientation-preserving")
    return x,e

def static_tetrad_4d(vierbein):
    e=np.asarray(vierbein,dtype=float)
    if e.ndim!=3 or e.shape[1:]!=(3,3): raise ValueError("vierbein must have shape (N,3,3)")
    out=np.zeros((len(e),4,4)); out[:,0,0]=1.; out[:,1:,1:]=e
    return out

def curved_gamma_matrices(vierbein,gammas=None):
    e=np.asarray(vierbein,dtype=float); gs=tuple(gamma_matrices() if gammas is None else gammas)
    if e.ndim!=3 or e.shape[1:]!=(3,3) or len(gs)!=4: raise ValueError("expected (N,3,3) vierbein and four gamma matrices")
    out=np.zeros((len(e),4,4,4),dtype=complex); out[:,0]=gs[0]
    for i,ei in enumerate(e):
        inv=np.linalg.inv(ei)
        for mu in range(3): out[i,mu+1]=sum(inv[mu,a]*gs[a+1] for a in range(3))
    return out

def curved_clifford_audit(vierbein,tol=2e-10):
    e=np.asarray(vierbein,dtype=float); curved=curved_gamma_matrices(e); mc=0.; mm=0.
    for i,ei in enumerate(e):
        tetrad=np.eye(4); tetrad[1:,1:]=ei
        gc=tetrad.T@ETA@tetrad; gi=np.linalg.inv(gc)
        mm=max(mm,float(np.linalg.norm(gc@gi-np.eye(4))))
        for mu in range(4):
            for nu in range(4):
                lhs=curved[i,mu]@curved[i,nu]+curved[i,nu]@curved[i,mu]
                mc=max(mc,float(np.linalg.norm(lhs-2*gi[mu,nu]*np.eye(4),ord="fro")))
    return {"vertices":len(e),"max_clifford_residual":mc,"max_metric_inverse_residual":mm,"passed":max(mc,mm)<=tol,"time_leg":"static e^0_0=1"}

def _rotation_to_quaternion(r):
    r=np.asarray(r,dtype=float)
    if r.shape!=(3,3) or np.linalg.norm(r.T@r-np.eye(3))>1e-7: raise ValueError("rotation must be orthogonal")
    if np.linalg.det(r)<0: raise ValueError("reflection requires Pin, not SU(2)")
    tr=float(np.trace(r))
    if tr>0:
        s=2*np.sqrt(max(0.,tr+1)); q=np.array([.25*s,(r[2,1]-r[1,2])/s,(r[0,2]-r[2,0])/s,(r[1,0]-r[0,1])/s])
    else:
        k=int(np.argmax(np.diag(r)))
        if k==0:
            s=2*np.sqrt(max(0.,1+r[0,0]-r[1,1]-r[2,2])); q=np.array([(r[2,1]-r[1,2])/s,.25*s,(r[0,1]+r[1,0])/s,(r[0,2]+r[2,0])/s])
        elif k==1:
            s=2*np.sqrt(max(0.,1+r[1,1]-r[0,0]-r[2,2])); q=np.array([(r[0,2]-r[2,0])/s,(r[0,1]+r[1,0])/s,.25*s,(r[1,2]+r[2,1])/s])
        else:
            s=2*np.sqrt(max(0.,1+r[2,2]-r[0,0]-r[1,1])); q=np.array([(r[1,0]-r[0,1])/s,(r[0,2]+r[2,0])/s,(r[1,2]+r[2,1])/s,.25*s])
    q/=np.linalg.norm(q)
    return -q if q[0]<0 else q

def polar_spin_links(vierbein,edges):
    """Candidate spin transport: lift the polar rotation of E_j E_i^{-1}."""
    e=np.asarray(vierbein,dtype=float); links={}; residual=0.
    for i,j in edges:
        if not (0<=i<j<len(e)): raise ValueError("edges must be canonical i<j")
        m=e[j]@np.linalg.inv(e[i]); u,_,vh=np.linalg.svd(m); r=u@vh
        if np.linalg.det(r)<0: u[:,-1]*=-1; r=u@vh
        residual=max(residual,float(np.linalg.norm(r.T@r-np.eye(3))))
        links[(i,j)]=quaternion_to_fundamental_su2(_rotation_to_quaternion(r))
    return links,residual

def compact_kernel_transform(features,edges,kernels:Mapping,transporters:Mapping):
    """Finite-support field transform y_target += K_target,source P_target,source f_source."""
    f=np.asarray(features,dtype=complex)
    if f.ndim!=2: raise ValueError("features must have shape (N,C)")
    outdim=next(iter(kernels.values())).shape[0] if kernels else f.shape[1]
    y=np.zeros((len(f),outdim),dtype=complex); support=set()
    for a,b in edges:
        for source,target in ((a,b),(b,a)):
            key=(source,target)
            if key not in kernels: continue
            if key not in transporters: raise KeyError(f"missing transporter {key}")
            k=np.asarray(kernels[key],dtype=complex); p=np.asarray(transporters[key],dtype=complex)
            if k.shape[1]!=p.shape[0] or p.shape[1]!=f.shape[1]: raise ValueError("kernel/transporter/feature dimensions disagree")
            y[target]+=k@p@f[source]; support.add(key)
    return y,support


def assemble_wilson_laplacian(
    coordinates, vierbein, edges, spin_links=None, wilson_r=1.0, edge_lengths=None
):
    """Build the positive, gauge-covariant graph Wilson term.

    For each canonical edge i<j with weight w_ij=r/(2 length_ij):
      L[j,i] -= w_ij S_ji; L[i,j] -= w_ij S_ji^dagger
      L[i,i] += w_ij I;     L[j,j] += w_ij I.
    The returned operator is volume-normalized by V^(-1/2) on both sides.
    edge_lengths optionally supplies physical lengths for periodic edges
    whose endpoint coordinate difference is not the periodic geodesic length.
    """
    x,e=_validate_geometry(coordinates,vierbein)
    r=float(wilson_r)
    if not np.isfinite(r) or r<0:
        raise ValueError("wilson_r must be finite and non-negative")
    n=len(x)
    if spin_links is None:
        spin_links,_=polar_spin_links(e,edges)
    raw=np.zeros((4*n,4*n),dtype=complex)
    identity=np.eye(4,dtype=complex)
    min_length=np.inf
    for i,j in edges:
        if not (0<=i<j<n):
            raise ValueError("edges must be canonical i<j")
        if edge_lengths is None:
            length=float(np.linalg.norm(x[j]-x[i]))
        else:
            length=float(edge_lengths.get((i,j),np.linalg.norm(x[j]-x[i])))
        if not np.isfinite(length) or length<=1e-14:
            raise ValueError(f"edge {(i,j)} must have a positive finite length")
        min_length=min(min_length,length)
        sji=np.asarray(spin_links[(i,j)],dtype=complex)
        if sji.shape==(2,2):
            sji=np.block([[sji,np.zeros((2,2))],[np.zeros((2,2)),sji]])
        if sji.shape!=(4,4):
            raise ValueError("spin link must be 2x2 SU(2) or 4x4")
        si=slice(4*i,4*i+4); sj=slice(4*j,4*j+4)
        weight=r/(2.0*length)
        raw[si,si]+=weight*identity
        raw[sj,sj]+=weight*identity
        raw[sj,si]-=weight*sji
        raw[si,sj]-=weight*sji.conj().T
    volumes=np.linalg.det(e)
    v4=np.repeat(volumes,4)
    invsqrt=1.0/np.sqrt(v4)
    operator=invsqrt[:,None]*raw*invsqrt[None,:]
    herm=float(np.linalg.norm(operator-operator.conj().T,ord="fro"))
    eigmin=float(np.linalg.eigvalsh(operator)[0]) if operator.size else 0.0
    return operator,{
        "vertices":n,"edges":len(edges),"wilson_r":r,
        "min_edge_length":None if not np.isfinite(min_length) else float(min_length),
        "volume_min":float(volumes.min()),"volume_max":float(volumes.max()),
        "hermiticity_residual":herm,"minimum_eigenvalue":eigmin,
        "status":"GAUGE_COVARIANT_WILSON_LAPLACIAN"
    }


def chirality_matrix():
    """Return gamma_5=i gamma^0 gamma^1 gamma^2 gamma^3 in this convention."""
    gs=gamma_matrices()
    return 1j*gs[0]@gs[1]@gs[2]@gs[3]


def assemble_spatial_dirac(coordinates,vierbein,edges,spin_links=None,local_spatial_gamma=None,wilson_r=0.0):
    """Hermitian volume-weighted nearest-neighbour spatial Dirac candidate."""
    x,e=_validate_geometry(coordinates,vierbein); n=len(x); gs=gamma_matrices()
    if local_spatial_gamma is None: local=np.broadcast_to(np.asarray(gs[1:]),(n,3,4,4)).copy()
    else:
        local=np.asarray(local_spatial_gamma,dtype=complex)
        if local.shape!=(n,3,4,4): raise ValueError("local_spatial_gamma must have shape (N,3,4,4)")
    if spin_links is None: spin_links,_=polar_spin_links(e,edges)
    h=np.zeros((4*n,4*n),dtype=complex)
    for i,j in edges:
        dx=x[j]-x[i]; length=float(np.linalg.norm(dx))
        if length<=1e-14: raise ValueError("edge endpoints must be distinct")
        vi=e[i]@dx; vj=e[j]@dx; ni=vi/np.linalg.norm(vi); nj=vj/np.linalg.norm(vj)
        ci=sum(ni[a]*local[i,a] for a in range(3)); cj=sum(nj[a]*local[j,a] for a in range(3))
        sji=np.asarray(spin_links[(i,j)],dtype=complex)
        if sji.shape==(2,2): sji=np.block([[sji,np.zeros((2,2))],[np.zeros((2,2)),sji]])
        if sji.shape!=(4,4): raise ValueError("spin link must be 2x2 SU(2) or 4x4")
        block=.5j/length*(cj@sji+sji@ci)
        si=slice(4*i,4*i+4); sj=slice(4*j,4*j+4)
        h[sj,si]+=block; h[si,sj]+=block.conj().T
    volumes=np.linalg.det(e); v4=np.repeat(volumes,4); invsqrt=1/np.sqrt(v4)
    d=invsqrt[:,None]*h*invsqrt[None,:]
    r=float(wilson_r)
    if not np.isfinite(r) or r<0: raise ValueError("wilson_r must be finite and non-negative")
    wilson_audit=None
    if r>0:
        wilson,wilson_audit=assemble_wilson_laplacian(x,e,edges,spin_links=spin_links,wilson_r=r)
        d=d+wilson
    return d,{"vertices":n,"edges":len(edges),"spinor_dimension":4,"volume_min":float(volumes.min()),"volume_max":float(volumes.max()),"hermiticity_residual":float(np.linalg.norm(d-d.conj().T,ord="fro")),"operator_norm":float(np.linalg.norm(d,ord=2)),"spin_links_derived_from":"right-polar factor of E_j E_i^{-1}","wilson_r":r,"wilson_audit":wilson_audit,"status":"CANDIDATE_DISCRETE_SPATIAL_DIRAC_WITH_WILSON_TERM" if r>0 else "CANDIDATE_DISCRETE_SPATIAL_DIRAC"}

def discrete_dirac_spectrum(operator):
    d=np.asarray(operator,dtype=complex); residual=float(np.linalg.norm(d-d.conj().T,ord="fro"))
    if residual>1e-9: raise ValueError(f"operator is not Hermitian: {residual}")
    values=np.linalg.eigvalsh(d)
    return {"eigenvalues":values,"min":float(values[0]),"max":float(values[-1]),"near_zero_count_1e-9":int(np.sum(np.abs(values)<1e-9)),"spectral_pairing_residual":float(np.max(np.abs(values+values[::-1])))}

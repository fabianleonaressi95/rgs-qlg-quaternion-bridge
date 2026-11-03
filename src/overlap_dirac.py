"""Finite-dimensional overlap construction from the repository's Wilson pair.

The repository stores a Hermitian Hamiltonian-like naive operator H0 with
{Gamma5,H0}=0. The Euclidean gamma5-Hermitian Wilson kernel used here is
D_E = i H0 + W, and H_W = Gamma5 (a D_E - rho I).
The overlap operator is D_ov = (rho/a) [I + Gamma5 sign(H_W)].

This is an explicit algebraic conversion, not a claim that the Q6 operator
already has a continuum Euclidean Dirac interpretation.
"""
from __future__ import annotations
import numpy as np


def assemble_overlap_from_wilson(
    naive_hermitian, wilson_laplacian, gamma5, lattice_spacing=1.0,
    rho=1.0, zero_tolerance=1e-11
):
    h0=np.asarray(naive_hermitian,dtype=complex)
    w=np.asarray(wilson_laplacian,dtype=complex)
    g=np.asarray(gamma5,dtype=complex)
    a=float(lattice_spacing); rho=float(rho)
    if h0.ndim!=2 or h0.shape[0]!=h0.shape[1] or w.shape!=h0.shape or g.shape!=h0.shape:
        raise ValueError("naive, Wilson and Gamma5 matrices must be equal-size square matrices")
    if not np.isfinite(a) or a<=0 or not np.isfinite(rho) or rho<=0:
        raise ValueError("lattice_spacing and rho must be positive and finite")
    if np.linalg.norm(h0-h0.conj().T)>1e-9 or np.linalg.norm(w-w.conj().T)>1e-9:
        raise ValueError("naive and Wilson operators must be Hermitian")
    if np.linalg.norm(g-g.conj().T)>1e-9 or np.linalg.norm(g@g-np.eye(len(g)))>1e-9:
        raise ValueError("Gamma5 must be a Hermitian involution")
    naive_chiral=float(np.linalg.norm(h0@g+g@h0,ord="fro"))
    wilson_comm=float(np.linalg.norm(w@g-g@w,ord="fro"))
    if naive_chiral>1e-8 or wilson_comm>1e-8:
        raise ValueError("overlap kernel requires {Gamma5,H0}=0 and [Gamma5,W]=0")
    d_e=1j*h0+w
    kernel=g@(a*d_e-rho*np.eye(len(g)))
    kernel_herm=float(np.linalg.norm(kernel-kernel.conj().T,ord="fro"))
    if kernel_herm>1e-8:
        raise ValueError(f"gamma5 Wilson kernel is not Hermitian: {kernel_herm}")
    values,vectors=np.linalg.eigh((kernel+kernel.conj().T)/2)
    min_gap=float(np.min(np.abs(values)))
    if min_gap<=zero_tolerance:
        raise ValueError("Wilson kernel has a zero/near-zero eigenvalue; sign is undefined")
    sign=(vectors*np.sign(values))@vectors.conj().T
    identity=np.eye(len(g),dtype=complex)
    overlap=(rho/a)*(identity+g@sign)
    gw=float(np.linalg.norm(g@overlap+overlap@g-(a/rho)*overlap@g@overlap,ord="fro"))
    gamma5_herm=float(np.linalg.norm(overlap.conj().T-g@overlap@g,ord="fro"))
    return overlap,{
        "dimension":len(g),
        "lattice_spacing":a,
        "rho":rho,
        "wilson_kernel_hermiticity_residual":kernel_herm,
        "wilson_kernel_minimum_absolute_eigenvalue":min_gap,
        "ginsparg_wilson_residual":gw,
        "gamma5_hermiticity_residual":gamma5_herm,
        "naive_chiral_residual":naive_chiral,
        "wilson_gamma5_commutator_residual":wilson_comm,
        "status":"FINITE_MATRIX_OVERLAP_GW_OPERATOR",
        "continuum_locality_and_index_theorem":"OPEN"
    }


def overlap_index(overlap, gamma5, lattice_spacing=1.0, rho=1.0,
                  zero_tolerance=1e-8, identity_tolerance=1e-7):
    """Compute the finite-matrix GW index trace and resolve exact zero-mode chirality.

    Convention: index = Tr[Gamma5 (I - a D_ov/(2 rho))].
    The returned nearest integer is explicitly labelled a candidate unless
    the trace is real and within identity_tolerance of that integer.
    """
    d=np.asarray(overlap,dtype=complex); g=np.asarray(gamma5,dtype=complex)
    a=float(lattice_spacing); rho=float(rho)
    if d.ndim!=2 or d.shape[0]!=d.shape[1] or g.shape!=d.shape:
        raise ValueError("overlap and Gamma5 must be equal-size square matrices")
    if not np.isfinite(a) or a<=0 or not np.isfinite(rho) or rho<=0:
        raise ValueError("lattice_spacing and rho must be positive and finite")
    eye=np.eye(len(d),dtype=complex)
    if np.linalg.norm(g-g.conj().T)>1e-9 or np.linalg.norm(g@g-eye)>1e-9:
        raise ValueError("Gamma5 must be a Hermitian involution")
    gw=float(np.linalg.norm(g@d+d@g-(a/rho)*d@g@d,ord="fro"))
    gh=float(np.linalg.norm(d.conj().T-g@d@g,ord="fro"))
    if gw>identity_tolerance or gh>identity_tolerance:
        raise ValueError("index trace requires a GW operator with gamma5-hermiticity")
    raw=np.trace(g@(eye-(a/(2.0*rho))*d))
    real=float(np.real(raw)); imag=float(abs(np.imag(raw)))
    nearest=int(np.rint(real))
    integer_residual=float(abs(real-nearest))
    vals,vecs=np.linalg.eig(d)
    zero_ids=np.flatnonzero(np.abs(vals)<=zero_tolerance)
    if len(zero_ids):
        z=vecs[:,zero_ids]
        q,_=np.linalg.qr(z)
        restricted=q.conj().T@g@q
        chiralities=np.linalg.eigvalsh((restricted+restricted.conj().T)/2)
        zero_chirality_sum=float(np.sum(chiralities))
        chirality_eigenvalues=chiralities.tolist()
    else:
        zero_chirality_sum=0.0
        chirality_eigenvalues=[]
    return {
        "index_trace_real":real,
        "index_trace_imag_abs":imag,
        "nearest_integer_candidate":nearest,
        "integer_residual":integer_residual,
        "integer_within_tolerance":bool(imag<=identity_tolerance and integer_residual<=identity_tolerance),
        "zero_mode_count":int(len(zero_ids)),
        "zero_mode_chirality_eigenvalues":chirality_eigenvalues,
        "zero_mode_chirality_sum":zero_chirality_sum,
        "gw_residual":gw,
        "gamma5_hermiticity_residual":gh,
        "convention":"Tr[Gamma5 (I - a D_ov/(2 rho))]",
        "status":"FINITE_MATRIX_INDEX_TRACE; CONTINUUM_INDEX_THEOREM_OPEN"
    }


def overlap_kernel_spectral_flow(naive_hermitian, wilson_laplacian, gamma5,
                                 lattice_spacing, rho_values):
    """Track Hermitian H_W(rho)=Gamma5(a(iD0+W)-rho I) across rho.

    This is a finite-dimensional spectral-flow diagnostic. A sampled scan
    cannot certify all crossings between samples or a continuum index theorem.
    """
    h0=np.asarray(naive_hermitian,dtype=complex)
    w=np.asarray(wilson_laplacian,dtype=complex)
    g=np.asarray(gamma5,dtype=complex)
    a=float(lattice_spacing); rhos=np.asarray(rho_values,dtype=float)
    if h0.ndim!=2 or h0.shape[0]!=h0.shape[1] or w.shape!=h0.shape or g.shape!=h0.shape:
        raise ValueError("naive, Wilson and Gamma5 matrices must have equal square shape")
    if not np.isfinite(a) or a<=0 or rhos.ndim!=1 or len(rhos)<2 or not np.all(np.isfinite(rhos)):
        raise ValueError("positive spacing and at least two finite rho samples are required")
    if np.any(rhos<=0):
        raise ValueError("all rho samples must be positive")
    if np.linalg.norm(h0-h0.conj().T)>1e-9 or np.linalg.norm(w-w.conj().T)>1e-9:
        raise ValueError("naive and Wilson operators must be Hermitian")
    if np.linalg.norm(g-g.conj().T)>1e-9 or np.linalg.norm(g@g-np.eye(len(g)))>1e-9:
        raise ValueError("Gamma5 must be a Hermitian involution")
    if np.linalg.norm(h0@g+g@h0)>1e-8 or np.linalg.norm(w@g-g@w)>1e-8:
        raise ValueError("spectral flow requires {Gamma5,D0}=0 and [Gamma5,W]=0")
    de=1j*h0+w
    rows=[]
    for rho in rhos:
        hw=g@(a*de-rho*np.eye(len(g)))
        herm=float(np.linalg.norm(hw-hw.conj().T,ord="fro"))
        if herm>1e-8:
            raise ValueError("sampled Wilson kernel is not Hermitian")
        ev=np.linalg.eigvalsh((hw+hw.conj().T)/2)
        rows.append({
            "rho":float(rho),
            "minimum_absolute_eigenvalue":float(np.min(np.abs(ev))),
            "negative_eigenvalue_count":int(np.count_nonzero(ev<0)),
            "positive_eigenvalue_count":int(np.count_nonzero(ev>0)),
            "eigenvalues":ev.tolist()
        })
    crossings=[]
    for left,right in zip(rows[:-1],rows[1:]):
        if left["negative_eigenvalue_count"]!=right["negative_eigenvalue_count"]:
            crossings.append({"rho_left":left["rho"],"rho_right":right["rho"],
                              "negative_count_left":left["negative_eigenvalue_count"],
                              "negative_count_right":right["negative_eigenvalue_count"]})
    return {
        "samples":rows,
        "sampled_sign_count_changes":crossings,
        "status":"FINITE_SAMPLED_SPECTRAL_FLOW; BETWEEN_SAMPLE_CROSSINGS_NOT_RESOLVED"
    }


def overlap_locality_profile(overlap, edges, vertex_count, internal_dimension=4):
    """Summarize block norms by graph distance; not a uniform locality bound."""
    d=np.asarray(overlap,dtype=complex)
    n=int(vertex_count); k=int(internal_dimension)
    if n<1 or k<1 or d.shape!=(n*k,n*k):
        raise ValueError("operator dimension must equal vertex_count*internal_dimension")
    adjacency=[set() for _ in range(n)]
    for i,j in edges:
        i=int(i); j=int(j)
        if i==j or i<0 or j<0 or i>=n or j>=n:
            raise ValueError("edges must join distinct valid vertices")
        adjacency[i].add(j); adjacency[j].add(i)
    distances=np.full((n,n),-1,dtype=int)
    for root in range(n):
        distances[root,root]=0; queue=[root]
        for u in queue:
            for v in adjacency[u]:
                if distances[root,v]<0:
                    distances[root,v]=distances[root,u]+1; queue.append(v)
    if np.any(distances<0):
        raise ValueError("graph must be connected for a single-shell locality profile")
    shells={}
    for i in range(n):
        for j in range(n):
            if i==j: continue
            shell=int(distances[i,j])
            bi=d[k*i:k*(i+1),k*j:k*(j+1)]
            shells.setdefault(shell,[]).append(float(np.linalg.norm(bi,ord="fro")))
    profile={str(s):{"count":len(v),"mean":float(np.mean(v)),
                     "max":float(np.max(v))} for s,v in sorted(shells.items())}
    usable=[(s,v) for s,v in sorted(shells.items()) if np.mean(v)>0]
    if len(usable)>=2:
        slope=float(np.polyfit([x[0] for x in usable],
                               np.log([np.mean(x[1]) for x in usable]),1)[0])
        factor=float(np.exp(slope))
    else:
        slope=None; factor=None
    return {
        "distance":"shortest-path graph distance",
        "shells":profile,
        "log_mean_fit_slope":slope,
        "empirical_factor_per_step":factor,
        "status":"FINITE_GRAPH_DIAGNOSTIC; NO_UNIFORM_EXPONENTIAL_LOCALITY_CLAIM"
    }

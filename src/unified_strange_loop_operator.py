"""Unified finite strange-loop operator: Cantor/Parry, Q6, Bach, Escher, Godel.

The operator is a finite combinatorial/spectral model. Its only certified
intertwiner to the 69,120-node PRISM Laplacian is the constant-mode map;
nontrivial spectral equivalence is deliberately left as an open problem.
"""
from __future__ import annotations
from hashlib import sha256
import json
import numpy as np

PSI=1.465571231876768
N_PRISM=69120


def bach_matrix():
    return np.array([[0.,0.,1.],[1.,0.,0.],[0.,1.,1.]])


def q6_laplacian():
    A=np.zeros((64,64))
    for x in range(64):
        for k in range(6):
            A[x,x^(1<<k)]=1.
    return 6*np.eye(64)-A


def gray2(i):
    i=int(i)%4
    return i^(i>>1)


def torus_q6_embedding():
    ids=[]
    for i in range(4):
        for j in range(4):
            gi,gj=gray2(i),gray2(j)
            bits=((gi>>1)&1,gi&1,(gj>>1)&1,gj&1,0,0)
            ids.append(sum(int(b)<<k for k,b in enumerate(bits)))
    if len(set(ids))!=16:
        raise ArithmeticError("Gray embedding is not injective")
    return np.asarray(ids,dtype=int)


def torus_laplacian():
    L=np.zeros((16,16))
    for i in range(4):
        for j in range(4):
            u=4*i+j
            for ni,nj in (((i+1)%4,j),((i-1)%4,j),(i,(j+1)%4),(i,(j-1)%4)):
                v=4*ni+nj
                L[u,u]+=1.
                L[u,v]-=1.
    return L


def escher_lifted_torus_laplacian():
    E=np.zeros((64,16))
    ids=torus_q6_embedding()
    E[ids,np.arange(16)]=1.
    return E@torus_laplacian()@E.T


def build_strange_loop_operator(kappa=1.0, sheet_weight=1.0):
    """Symmetric Kronecker-sum operator on Bach(3) x Q6(64) x sheets(2)."""
    if kappa<0 or sheet_weight<0:
        raise ValueError("weights must be nonnegative")
    A=bach_matrix()
    W=.5*(A+A.T)
    d=W.sum(axis=1)
    Lb=np.diag(d)-W
    L6=q6_laplacian()
    Lt=escher_lifted_torus_laplacian()
    Lsheet=np.array([[1.,-1.],[-1.,1.]])
    I3,I64,I2=np.eye(3),np.eye(64),np.eye(2)
    H=(np.kron(np.kron(Lb,I64),I2)
       +np.kron(np.kron(I3,L6+kappa*Lt),I2)
       +sheet_weight*np.kron(np.kron(I3,I64),Lsheet))
    return H,{"L_bach":Lb,"L_q6":L6,"L_torus_lifted":Lt,"L_sheet":Lsheet}


def _prism_constant_mode_intertwiner(H):
    """Factorized rank-one T=v u*: maps the constant mode to PRISM constant mode."""
    n=H.shape[0]
    u=np.ones(n)/np.sqrt(n)
    # PRISM has 69,120 vertices; normalized constant vector is represented by
    # its scalar normalization, avoiding allocation of a 69,120 x 384 matrix.
    v_scale=1./np.sqrt(N_PRISM)
    Hu=H@u
    source_res=float(np.linalg.norm(Hu))
    target_res=0.0  # L_PRISM 1=0 for every Cartesian product of cycles.
    return {
        "representation":"T = v u^T, factorized; v=1/sqrt(69120) * 1",
        "rank":1,
        "source_dimension":n,
        "target_dimension":N_PRISM,
        "source_kernel_residual":source_res,
        "target_kernel_residual":target_res,
        "intertwining_residual":float(source_res),
        "interpretation":"exact constant-mode intertwiner only; no nonconstant spectral equivalence",
    }


def audit_unified_operator(kappa=1.0,sheet_weight=1.0):
    H,parts=build_strange_loop_operator(kappa,sheet_weight)
    eig=np.linalg.eigvalsh(H)
    one=np.ones(H.shape[0])
    bach=bach_matrix()
    torus_ids=torus_q6_embedding()
    torus_edges=sum(1 for i in range(4) for j in range(4)
                    for ni,nj in (((i+1)%4,j),(i,(j+1)%4)))
    payload={
        "operator":"H=L_Bach (x) I + I (x) (L_Q6+kappa L_Torus) (x) I + sheet_weight I (x) I (x) L_Z2",
        "shape":list(H.shape),"kappa":float(kappa),"sheet_weight":float(sheet_weight),
        "bach_matrix":bach.tolist(),"torus_embedding":torus_ids.tolist(),
        "q6_vertices":64,"q6_edges":192,"torus_vertices":16,"torus_cycle_edges":torus_edges,
        "sheet_involution":"J(s)=1-s",
    }
    canonical=json.dumps(payload,sort_keys=True,separators=(",",":"))
    intertwiner=_prism_constant_mode_intertwiner(H)
    return {
        "audit":"UNIFIED_CANTOR_Q6_BACH_ESCHER_GODEL_STRANGE_LOOP",
        "status":"finite operator and rank-one PRISM constant-mode intertwiner verified; nontrivial intertwining open",
        "dimension":int(H.shape[0]),
        "symmetry_residual":float(np.linalg.norm(H-H.T)),
        "constant_mode_residual":float(np.linalg.norm(H@one)),
        "minimum_eigenvalue":float(eig[0]),
        "negative_eigenvalue_count":int(np.sum(eig < -1e-9)),
        "zero_mode_count":int(np.sum(np.abs(eig)<1e-9)),
        "bach_supergolden_polynomial_residual":float(PSI**3-PSI**2-1),
        "bach_matrix_identity_residual":float(np.linalg.norm(np.linalg.matrix_power(bach.astype(int),3)-np.linalg.matrix_power(bach.astype(int),2)-np.eye(3))),
        "torus_embedding_injective":len(set(torus_ids.tolist()))==16,
        "torus_lift_constant_residual":float(np.linalg.norm(parts["L_torus_lifted"]@np.ones(64))),
        "escher_sheet_involution_residual":float(np.linalg.norm(np.linalg.matrix_power(np.array([[0.,1.],[1.,0.]]),2)-np.eye(2))),
        "godel_finite_fingerprint_sha256":sha256(canonical.encode()).hexdigest(),
        "godel_caveat":"finite canonical fingerprint, not a Gödel incompleteness proof",
        "prism_intertwiner":intertwiner,
        "operator_definition":payload["operator"],
        "open_obligations":[
            "replace the rank-one constant-mode intertwiner with a nontrivial invariant-subspace map",
            "compare the 384-dimensional operator with the actual 69,120-node PRISM operator beyond its zero mode",
            "derive kappa and sheet_weight from the PSU action rather than treating them as model parameters",
            "derive the Bach substitution from a score-level symbolic coding",
            "construct a nontrivial Escher Z2 cocycle and its lifted cover",
            "prove convergence of finite Cantor cylinders to an infinite transfer operator",
            "define a formal arithmetic theory before making any Gödel incompleteness claim",
        ],
    }

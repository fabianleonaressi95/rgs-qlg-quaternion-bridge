"""Coupled Cantor--Q6--supergolden transfer operator.

Symbolic model only: its pressure and dimensions use the declared symbolic
metric. A geometric RGS/PSU embedding and Q6 intertwiner remain open.
"""
from __future__ import annotations
from math import log, sqrt
from typing import Iterable
import numpy as np

PHI=(1.0+sqrt(5.0))/2.0
PSI=1.465571231876768

def supergolden_matrix():
    return np.array([[1.,0.,1.],[1.,0.,0.],[0.,1.,0.]])

def perron_data():
    A=supergolden_matrix()
    r=np.array([1.,1./PSI,1./PSI**2])
    P=A*r[np.newaxis,:]/(PSI*r[:,np.newaxis])
    vals,vecs=np.linalg.eig(P.T)
    k=int(np.argmin(np.abs(vals-1.)))
    pi=np.real(vecs[:,k])
    if pi.sum()<0: pi=-pi
    pi=pi/pi.sum()
    return {"A":A,"psi":PSI,"perron_vector":r,"P":P,"stationary":pi,
            "perron_residual":float(np.linalg.norm(A@r-PSI*r)),
            "row_sum_residual":float(np.linalg.norm(P.sum(axis=1)-1.)),
            "stationary_residual":float(np.linalg.norm(pi@P-pi))}

def golden_branch_probabilities(phi=PHI):
    phi=float(phi)
    if not np.isfinite(phi) or phi<=1.: raise ValueError("phi must be finite and > 1")
    q_phi=float(np.exp(-1./phi**2)/phi)
    r=perron_data()["perron_vector"]
    gm=float(np.prod(r)**(1./3.))
    q=q_phi*r/gm
    p1=q/(1.+q); p0=1./(1.+q)
    return {"phi":phi,"q_phi":q_phi,"geometric_mean_perron":gm,
            "q_state":q,"p0":p0,"p1":p1}

def pressure_matrix(s,phi=PHI):
    s=float(s)
    if not np.isfinite(s): raise ValueError("s must be finite")
    d=perron_data(); z=golden_branch_probabilities(phi)
    A,P=d["A"],d["P"]; M=np.zeros((3,3))
    for i in range(3):
        branches=z["p0"][i]**s+z["p1"][i]**s
        for j in range(3):
            if A[i,j]:
                trans=1. if s==0. else P[i,j]**s
                M[i,j]=trans*branches
    return M

def pressure(s,phi=PHI):
    rho=float(np.max(np.abs(np.linalg.eigvals(pressure_matrix(s,phi)))))
    if not np.isfinite(rho) or rho<=0.: raise ArithmeticError("invalid pressure spectral radius")
    return float(np.log(rho))

def generalized_dimension(s,phi=PHI):
    s=float(s); den=log(3.*PSI)
    if abs(s-1.)>1e-7: return pressure(s,phi)/((1.-s)*den)
    h=1e-5
    deriv=(pressure(1.+h,phi)-pressure(1.-h,phi))/(2.*h)
    return float(-deriv/den)

def cylinder_probability(states: Iterable[int],bits: Iterable[int],phi=PHI):
    states=tuple(int(x) for x in states); bits=tuple(int(x) for x in bits)
    if len(states)!=len(bits)+1: raise ValueError("states length must be bits length + 1")
    if not states or any(x not in (0,1,2) for x in states): raise ValueError("invalid state path")
    if any(x not in (0,1) for x in bits): raise ValueError("bits must be binary")
    d=perron_data(); z=golden_branch_probabilities(phi)
    value=float(d["stationary"][states[0]])
    for k,b in enumerate(bits):
        i,j=states[k],states[k+1]; trans=d["P"][i,j]
        if trans==0.: return 0.
        value*=float(trans*(z["p0"][i] if b==0 else z["p1"][i]))
    return float(value)

def q6_laplacian():
    A=np.zeros((64,64))
    for v in range(64):
        for bit in range(6): A[v,v^(1<<bit)]=1.
    return 6.*np.eye(64)-A

def q6_spectrum_residual():
    actual=np.linalg.eigvalsh(q6_laplacian())
    expected=np.sort(np.concatenate([np.full(int(__import__("math").comb(6,j)),2.*j) for j in range(7)]))
    return float(np.max(np.abs(actual-expected)))

def intertwiner_residual(T,L_source,L_target):
    T=np.asarray(T); L_source=np.asarray(L_source); L_target=np.asarray(L_target)
    if T.ndim!=2 or L_source.shape!=(T.shape[1],T.shape[1]): raise ValueError("bad source dimensions")
    if L_target.shape!=(T.shape[0],T.shape[0]): raise ValueError("bad target dimensions")
    num=np.linalg.norm(T@L_source-L_target@T,ord="fro")
    return float(num/max(1.,np.linalg.norm(T@L_source,ord="fro")))

def audit_cantor_q6_supergolden():
    d=perron_data(); z=golden_branch_probabilities(); ss=[0.,.5,1.,2.,3.]
    radii={str(s):float(np.max(np.abs(np.linalg.eigvals(pressure_matrix(s))))) for s in ss}
    dims={str(s):generalized_dimension(s) for s in ss}
    return {
      "audit":"CANTOR_Q6_SUPERGOLDEN_COUPLED_TRANSFER",
      "status":"symbolic model; physical RGS-PSU realization remains open",
      "supergolden":{"polynomial_residual":float(PSI**3-PSI**2-1.),
        "matrix_identity_residual":float(np.linalg.norm(np.linalg.matrix_power(d["A"],3)-np.linalg.matrix_power(d["A"],2)-np.eye(3))),
        "perron_root":PSI,"perron_vector":d["perron_vector"].tolist(),
        "transition_matrix":d["P"].tolist(),"stationary_distribution":d["stationary"].tolist(),
        "row_sum_residual":d["row_sum_residual"],"stationary_residual":d["stationary_residual"]},
      "golden_measure":{"q_phi":z["q_phi"],"q_state":z["q_state"].tolist(),
        "p0":z["p0"].tolist(),"p1":z["p1"].tolist()},
      "pressure_spectral_radii":radii,"generalized_dimensions":dims,
      "q6":{"vertices":64,"edges":192,"degree":6,"laplacian_spectrum_residual":q6_spectrum_residual()},
      "metric":"synchronous symbolic metric ell_n=(3*psi)^(-n)",
      "open_obligations":["prove separation or control overlaps for a geometric embedding",
        "construct an intertwiner with the Q6 Laplacian","derive psi^(-n) from the RGS-PSU action",
        "establish finite-level operator convergence","derive physical observables rather than identify symbolic parameters with them"]}

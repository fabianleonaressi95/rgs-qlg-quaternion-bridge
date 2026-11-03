"""Unweighted Dirac-Kahler operator on the Q6 two-skeleton.

Cochain differential is assembled from the canonical Q6 incidence matrices:
d0 = boundary_1.T : C^0 -> C^1
d1 = boundary_2.T : C^1 -> C^2
D_K = d + d^dagger on C^0 direct-sum C^1 direct-sum C^2.

This is a finite combinatorial operator. No metric Hodge star, Whitney map,
continuum spinor equivalence, or physical fermion interpretation is assumed.
"""
from __future__ import annotations
import numpy as np
from src.q6_complex import build_q6_complex


def assemble_q6_dirac_kahler():
    complex_=build_q6_complex()
    d0=complex_.boundary_1.astype(float).T
    d1=complex_.boundary_2.astype(float).T
    if np.count_nonzero(d1@d0):
        raise ValueError("Q6 incidence data violate d1*d0=0")
    n0,n1,n2=d0.shape[1],d0.shape[0],d1.shape[0]
    total=n0+n1+n2
    operator=np.zeros((total,total),dtype=float)
    s0=slice(0,n0); s1=slice(n0,n0+n1); s2=slice(n0+n1,total)
    operator[s1,s0]=d0
    operator[s0,s1]=d0.T
    operator[s2,s1]=d1
    operator[s1,s2]=d1.T
    grading=np.diag(np.concatenate((np.ones(n0),-np.ones(n1),np.ones(n2))))
    lap0=d0.T@d0
    lap1=d0@d0.T+d1.T@d1
    lap2=d1@d1.T
    expected=np.zeros_like(operator)
    expected[s0,s0]=lap0
    expected[s1,s1]=lap1
    expected[s2,s2]=lap2
    herm=float(np.linalg.norm(operator-operator.T,ord="fro"))
    square=float(np.linalg.norm(operator@operator-expected,ord="fro"))
    chiral=float(np.linalg.norm(operator@grading+grading@operator,ord="fro"))
    return operator,{
        "cochain_dimensions":[n0,n1,n2],
        "total_dimension":total,
        "d_squared_residual":float(np.linalg.norm(d1@d0,ord="fro")),
        "hermiticity_residual":herm,
        "square_to_hodge_laplacian_residual":square,
        "degree_grading_anticommutator_residual":chiral,
        "status":"UNWEIGHTED_FINITE_Q6_DIRAC_KAHLER"
    }

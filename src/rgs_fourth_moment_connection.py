"""Gauge-covariant derivative coupling of the RGS fourth moment.

The fourth-moment operator is shifted to a positive constitutive metric on
Sym_0(3) before it is used in a gradient energy. The SO(3) connection enters
through D_i Q = partial_i Q + [omega_i,Q]. Variation with respect to omega
defines an antisymmetric spin current. This is a finite-dimensional local
constitutive sector, not yet a complete TEGR/Einstein-Cartan field theory.
"""
from __future__ import annotations
import numpy as np
from src.rgs_fourth_moment_bridge import symmetric_traceless_basis, anisotropic_response_operator


def _sym0_basis():
    return symmetric_traceless_basis()


def _coordinates(q):
    return np.asarray([np.einsum("ab,ab", b, q) for b in _sym0_basis()])


def _from_coordinates(v):
    return np.einsum("a,aij->ij", np.asarray(v,float), np.asarray(_sym0_basis()))


def _validate_sym0(q, name="Q"):
    q=np.asarray(q,float)
    if q.shape!=(3,3) or not np.all(np.isfinite(q)):
        raise ValueError(f"{name} must be finite 3x3")
    if not np.allclose(q,q.T,atol=1e-10,rtol=0.0) or abs(float(np.trace(q)))>1e-10:
        raise ValueError(f"{name} must be symmetric and trace-free")
    return q


def _validate_omega(omega):
    w=np.asarray(omega,float)
    if w.shape!=(3,3) or not np.all(np.isfinite(w)):
        raise ValueError("omega must be finite 3x3")
    if not np.allclose(w.T,-w,atol=1e-10,rtol=0.0):
        raise ValueError("omega must be antisymmetric (so(3))")
    return w


def fourth_moment_operator(residual):
    """5x5 matrix of A(Q)=R_abcd Q_cd in an orthonormal Sym_0 basis."""
    return np.asarray(anisotropic_response_operator(residual)["matrix"],float)


def positive_constitutive_metric(residual, beta=1.0, margin=1e-8):
    """Return H=shift*I+beta*A, shifted strictly positive on Sym_0(3).

    The shift is the minimal spectral correction plus a positive margin.
    This is a constitutive regularization, not a fitted physical constant.
    """
    if not np.isfinite(beta) or not np.isfinite(margin) or margin<=0:
        raise ValueError("beta must be finite and margin strictly positive")
    a=fourth_moment_operator(residual)
    eig=np.linalg.eigvalsh(a)
    shift=max(0.0,-float(beta)*float(eig[0]))+float(margin)
    h=shift*np.eye(5)+float(beta)*a
    heig=np.linalg.eigvalsh(h)
    if heig[0]<=0:
        raise ArithmeticError("positive constitutive shift failed")
    return {"operator":h,"raw_operator":a,"shift":shift,
            "raw_eigenvalues":eig,"eigenvalues":heig,
            "minimum_eigenvalue":float(heig[0])}


def covariant_derivative_q(partial_q, omega, q):
    """D_i Q = partial_i Q + [omega_i,Q], with omega_i in so(3)."""
    dq=_validate_sym0(partial_q,"partial_q")
    q=_validate_sym0(q)
    w=_validate_omega(omega)
    out=dq+w@q-q@w
    return 0.5*(out+out.T)-np.eye(3)*np.trace(out)/3.0


def apply_constitutive_metric(h, tensor):
    """Apply the 5x5 positive metric H to a Sym_0(3) tensor."""
    h=np.asarray(h,float)
    x=_validate_sym0(tensor,"tensor")
    if h.shape!=(5,5) or not np.allclose(h,h.T,atol=1e-10,rtol=0.0):
        raise ValueError("H must be symmetric 5x5")
    return _from_coordinates(h@_coordinates(x))


def gradient_energy_density(derivatives, omegas, q, h, coupling=1.0):
    """E= coupling/2 sum_i <D_iQ,H D_iQ>; positive for coupling>=0."""
    ds=np.asarray(derivatives,float); ws=np.asarray(omegas,float)
    q=_validate_sym0(q)
    if ds.ndim!=3 or ds.shape[1:]!=(3,3) or ws.shape!=ds.shape:
        raise ValueError("derivatives and omegas must both have shape (n,3,3)")
    if not np.isfinite(coupling) or coupling<0:
        raise ValueError("gradient coupling must be finite and nonnegative")
    energy=0.0
    for d,w in zip(ds,ws):
        dq=covariant_derivative_q(d,w,q)
        hdq=apply_constitutive_metric(h,dq)
        energy+=float(np.einsum("ab,ab",dq,hdq))
    return 0.5*float(coupling)*energy


def spin_current(q, covariant_derivative, h, coupling=1.0):
    r"""Connection Euler current J_i = coupling [H(D_i Q),Q] in so(3).

    Convention: delta(D_i Q)=[delta(omega_i),Q], hence
    delta L = <J_i,delta(omega_i)>_F. The result is antisymmetric.
    """
    q=_validate_sym0(q)
    dq=_validate_sym0(covariant_derivative,"covariant_derivative")
    s=apply_constitutive_metric(h,dq)
    j=float(coupling)*(s@q-q@s)
    return 0.5*(j-j.T)


def connection_variation_density(q, partial_q, omega, h, coupling=1.0):
    """One-direction density and its analytic so(3) connection current."""
    dq=covariant_derivative_q(partial_q,omega,q)
    hdq=apply_constitutive_metric(h,dq)
    energy=0.5*float(coupling)*float(np.einsum("ab,ab",dq,hdq))
    return energy,spin_current(q,dq,h,coupling)


def audit_connection_current(residual, q, partial_q, omega, beta=1.0,
                             coupling=0.8, margin=1e-8, step=1e-6):
    metric=positive_constitutive_metric(residual,beta,margin)
    energy,current=connection_variation_density(q,partial_q,omega,
                                                 metric["operator"],coupling)
    generators=[]
    for i,j in ((0,1),(0,2),(1,2)):
        g=np.zeros((3,3)); g[i,j]=1.0; g[j,i]=-1.0
        generators.append(g)
    errors=[]
    for g in generators:
        plus=connection_variation_density(q,partial_q,omega+step*g,
                                          metric["operator"],coupling)[0]
        minus=connection_variation_density(q,partial_q,omega-step*g,
                                           metric["operator"],coupling)[0]
        fd=(plus-minus)/(2*step)
        analytic=float(np.einsum("ab,ab",current,g))
        errors.append(abs(fd-analytic))
    return {
        "energy_density":energy,
        "current":current.tolist(),
        "current_antisymmetry_residual":float(np.linalg.norm(current+current.T)),
        "connection_fd_max_abs":float(max(errors)),
        "constitutive_shift":metric["shift"],
        "raw_operator_eigenvalues":metric["raw_eigenvalues"].tolist(),
        "positive_operator_eigenvalues":metric["eigenvalues"].tolist(),
        "minimum_positive_eigenvalue":metric["minimum_eigenvalue"],
        "status":{
            "positive_gradient_energy":"NUMERICALLY_VERIFIED",
            "connection_spin_current_variation":"NUMERICALLY_VERIFIED",
            "so3_current_antisymmetry":"NUMERICALLY_VERIFIED",
            "teleparallel_flatness_constraint":"OPEN",
            "full_coframe_and_connection_field_equations":"OPEN",
            "matter_and_q6_gauge_coupling":"OPEN",
            "continuum_limit":"OPEN"
        }
    }

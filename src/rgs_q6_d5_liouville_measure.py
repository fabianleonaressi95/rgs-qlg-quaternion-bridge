"""Exact pullback law between Q6 five-mode Liouville measure and D5 Bergman volume."""
from __future__ import annotations
import math
import numpy as np

def d5_generic_norm(z):
    z=np.asarray(z,dtype=complex)
    if z.shape!=(5,) or not np.all(np.isfinite(z)):
        raise ValueError("z must be finite with shape (5,)")
    w=np.dot(z,z)
    return float(1-2*np.vdot(z,z).real+abs(w)**2), w

def d5_bergman_metric(z):
    n,w=d5_generic_norm(z)
    z=np.asarray(z,dtype=complex)
    if n<=0 or abs(w)>=1:
        raise ValueError("z outside D5")
    dn=-2*np.conjugate(z)+2*z*np.conjugate(w)
    ddn=-2*np.eye(5,dtype=complex)+4*np.outer(z,np.conjugate(z))
    g=-5*(ddn/n-np.outer(dn,np.conjugate(dn))/n**2)
    return (g+g.conj().T)/2

def q6_maxwell_to_d5_chart(q,p,beta=1.0,eigenvalue=12.0):
    q,p=np.asarray(q,float),np.asarray(p,float)
    if q.shape!=(5,) or p.shape!=(5,) or not np.all(np.isfinite(q)) or not np.all(np.isfinite(p)):
        raise ValueError("q and p must be finite vectors of shape (5,)")
    if not math.isfinite(beta) or beta<=0 or not math.isfinite(eigenvalue) or eigenvalue<=0:
        raise ValueError("beta and eigenvalue must be finite and positive")
    # H = 1/2 [beta*lambda*q^2 + beta^(-1)*p^2].
    # The flow frequency is Omega=sqrt(lambda), and the canonical
    # reciprocal scale is kappa=beta*Omega (not sqrt(beta*lambda)).
    omega=math.sqrt(eigenvalue)
    kappa=beta*omega
    return (math.sqrt(kappa)*q+1j*p/math.sqrt(kappa))/math.sqrt(20.)

def q6_d5_bergman_liouville_report(q=None,p=None,beta=1.,eigenvalue=12.):
    if q is None: q=np.array([.07,-.03,.02,.01,-.04])
    if p is None: p=np.array([-.02,.05,.01,-.06,.03])
    z=q6_maxwell_to_d5_chart(q,p,beta,eigenvalue)
    n,w=d5_generic_norm(z)
    if n<=0 or abs(w)>=1: raise ValueError("selected point outside D5")
    g=d5_bergman_metric(z)
    detg=float(np.linalg.det(g).real)
    expected=10.**5/n**5
    jac=20.**-5
    # dVol_B=(i g_ij dz_i wedge dbar(z_j))^5/5!
    # hence dVol_B=2^5 det(g) d^10(Re z, Im z).
    pullback=float(2.**5*detg*jac)
    predicted=float(n**-5)
    g0=d5_bergman_metric(np.zeros(5,dtype=complex))
    omega=math.sqrt(eigenvalue)
    hamiltonian=0.5*(beta*eigenvalue*np.dot(q,q)+np.dot(p,p)/beta)
    complex_energy=10.0*omega*float(np.vdot(z,z).real)
    checks={
      "inside_D5":n>0 and abs(w)<1,
      "maxwell_energy_identity":math.isclose(hamiltonian,complex_energy,rel_tol=2e-12,abs_tol=2e-14),
      "positive_Bergman_metric":np.linalg.eigvalsh(g).min()>0,
      "determinant_identity":math.isclose(detg,expected,rel_tol=2e-9,abs_tol=1e-8),
      "chart_jacobian":math.isclose(jac,1/(32*10**5),rel_tol=1e-14),
      "measure_pullback":math.isclose(pullback,predicted,rel_tol=2e-9,abs_tol=1e-8),
      "origin_metric_10I":np.linalg.norm(g0-10*np.eye(5))<1e-12,
      "origin_volume_equals_Liouville":math.isclose(float(2**5*np.linalg.det(g0).real*jac),1.,rel_tol=1e-12)
    }
    return {
      "audit":"Q6_D5_BERGMAN_LIOUVILLE_PULLBACK_LAW",
      "frequency":omega,"canonical_complex_scale":beta*omega,
      "hamiltonian":hamiltonian,"complex_energy":complex_energy,"generic_norm":n,
      "metric_determinant":detg,"expected_metric_determinant":expected,
      "real_chart_jacobian":jac,"bergman_volume_pullback_density":pullback,
      "predicted_density":predicted,"exact_law":"Phi^*(dVol_B)=N(Phi(q,p))^(-5) dmu_L",
      "origin_metric":g0,"checks":checks,"all_checks_passed":bool(all(checks.values())),
      "scope":"Exact volume identity for this explicit chart on its D5 preimage; symplectic equality holds at the origin, not asserted globally."
    }


def q6_maxwell_oscillator_flow(q,p,time,beta=1.0,eigenvalue=12.0):
    """Exact isotropic five-mode Hamiltonian flow for H=(beta*lambda*q²+p²/beta)/2."""
    q,p=np.asarray(q,float),np.asarray(p,float)
    if q.shape!=(5,) or p.shape!=(5,) or not np.all(np.isfinite(q)) or not np.all(np.isfinite(p)):
        raise ValueError("q and p must be finite vectors of shape (5,)")
    if not math.isfinite(beta) or beta<=0 or not math.isfinite(eigenvalue) or eigenvalue<=0 or not math.isfinite(time):
        raise ValueError("beta, eigenvalue and time must be finite; beta and eigenvalue positive")
    omega=math.sqrt(eigenvalue)
    c,s=math.cos(omega*time),math.sin(omega*time)
    return c*q+s*p/(beta*omega), c*p-s*beta*omega*q


def q6_d5_liouville_flow_report(q=None,p=None,time=0.37,beta=1.7,eigenvalue=12.0):
    """Audit common-phase flow and preservation of N^(-5) dmu_L."""
    if q is None: q=np.array([.07,-.03,.02,.01,-.04])
    if p is None: p=np.array([-.02,.05,.01,-.06,.03])
    q,p=np.asarray(q,float),np.asarray(p,float)
    z0=q6_maxwell_to_d5_chart(q,p,beta,eigenvalue)
    qt,pt=q6_maxwell_oscillator_flow(q,p,time,beta,eigenvalue)
    zt=q6_maxwell_to_d5_chart(qt,pt,beta,eigenvalue)
    omega=math.sqrt(eigenvalue)
    expected=np.exp(-1j*omega*time)*z0
    n0,_=d5_generic_norm(z0)
    nt,_=d5_generic_norm(zt)
    checks={
      "common_phase_flow":np.linalg.norm(zt-expected)<2e-12,
      "generic_norm_preserved":math.isclose(n0,nt,rel_tol=2e-12,abs_tol=2e-14),
      "bergman_density_preserved":math.isclose(n0**-5,nt**-5,rel_tol=2e-11,abs_tol=2e-13),
      "inside_D5_initial":n0>0,
      "inside_D5_final":nt>0,
    }
    return {"audit":"Q6_D5_LIOUVILLE_WEIGHTED_FLOW","frequency":omega,
            "time":float(time),"generic_norm_initial":n0,"generic_norm_final":nt,
            "bergman_density_initial":n0**-5,"bergman_density_final":nt**-5,
            "checks":checks,"all_checks_passed":bool(all(checks.values())),
            "scope":"Exact isotropic five-mode flow; no claim for anisotropic or arbitrary Hamiltonians."}

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
    omega=math.sqrt(beta*eigenvalue)
    return (math.sqrt(omega)*q+1j*p/math.sqrt(omega))/math.sqrt(20.)

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
    checks={
      "inside_D5":n>0 and abs(w)<1,
      "positive_Bergman_metric":np.linalg.eigvalsh(g).min()>0,
      "determinant_identity":math.isclose(detg,expected,rel_tol=2e-9,abs_tol=1e-8),
      "chart_jacobian":math.isclose(jac,1/(32*10**5),rel_tol=1e-14),
      "measure_pullback":math.isclose(pullback,predicted,rel_tol=2e-9,abs_tol=1e-8),
      "origin_metric_10I":np.linalg.norm(g0-10*np.eye(5))<1e-12,
      "origin_volume_equals_Liouville":math.isclose(float(2**5*np.linalg.det(g0).real*jac),1.,rel_tol=1e-12)
    }
    return {
      "audit":"Q6_D5_BERGMAN_LIOUVILLE_PULLBACK_LAW",
      "frequency":math.sqrt(beta*eigenvalue),"generic_norm":n,
      "metric_determinant":detg,"expected_metric_determinant":expected,
      "real_chart_jacobian":jac,"bergman_volume_pullback_density":pullback,
      "predicted_density":predicted,"exact_law":"Phi^*(dVol_B)=N(Phi(q,p))^(-5) dmu_L",
      "origin_metric":g0,"checks":checks,"all_checks_passed":bool(all(checks.values())),
      "scope":"Exact volume identity for this explicit chart on its D5 preimage; symplectic equality holds at the origin, not asserted globally."
    }

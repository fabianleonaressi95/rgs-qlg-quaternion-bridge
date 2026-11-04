"""Q6 hypercube curvature law, Green--Kubo response, and supergolden scaling.

Finite identities and conditional response formulas are implemented separately
from any claim of a physical Yang--Mills/PSU continuum theory.
"""
from __future__ import annotations
from dataclasses import dataclass
from math import sqrt
import numpy as np
from src.q6_geometric_hodge import build_unit_cube_complex

@dataclass(frozen=True)
class GreenKuboResult:
    eta: float
    correlation_integral: float
    eigenvalues: np.ndarray
    modal_amplitudes: np.ndarray
    modal_integrals: np.ndarray
    zero_mode_amplitude: float
    convergent: bool

def supergolden_constant(tol: float = 1e-15, max_iter: int = 256) -> float:
    """Positive real root of x^3-x^2-1, computed by bisection."""
    if tol <= 0 or max_iter < 1: raise ValueError("tol must be positive and max_iter >= 1")
    lo, hi = 1.0, 2.0
    for _ in range(max_iter):
        mid = 0.5*(lo+hi)
        if mid**3-mid**2-1.0 > 0: hi=mid
        else: lo=mid
        if hi-lo <= tol: break
    return 0.5*(lo+hi)

def supergolden_companion_matrix() -> np.ndarray:
    """Companion matrix with characteristic polynomial x^3-x^2-1."""
    return np.array([[1.,0.,1.],[1.,0.,0.],[0.,1.,0.]])

def supergolden_geometric_sum(a0: float, n_terms: int | None = None) -> float:
    """Sum a_n=a0/psi^n; None returns the infinite sum."""
    if not np.isfinite(a0): raise ValueError("a0 must be finite")
    psi=supergolden_constant()
    if n_terms is None: return float(a0*psi/(psi-1.))
    if not isinstance(n_terms,int) or n_terms<0: raise ValueError("n_terms must be non-negative or None")
    if n_terms==0: return 0.
    return float(a0*(1.-psi**(-n_terms))/(1.-1./psi))

def supergolden_refinement_scales(scale0: float, levels: int) -> np.ndarray:
    """Candidate nested scales scale_n=scale0*psi**(-n), not a physical law."""
    if not np.isfinite(scale0) or scale0<=0: raise ValueError("scale0 must be finite and positive")
    if not isinstance(levels,int) or levels<0: raise ValueError("levels must be non-negative")
    return scale0*supergolden_constant()**(-np.arange(levels,dtype=float))

def supergolden_spectral_ratio_audit(eigenvalues: np.ndarray) -> dict[str,object]:
    """Compare positive spectral ratios with psi^2 and phi^2; no fit is imposed."""
    values=np.asarray(eigenvalues,dtype=float)
    if values.ndim!=1 or not np.all(np.isfinite(values)): raise ValueError("eigenvalues must be finite 1D")
    positive=np.sort(values[values>0])
    if positive.size<2: raise ValueError("at least two positive eigenvalues are required")
    ratios=positive/positive[0]
    psi,phi=supergolden_constant(),(1.+sqrt(5.))/2.
    return {"positive_eigenvalue_count":int(positive.size),"ratios_to_minimum":ratios.tolist(),
        "supergolden_squared":float(psi**2),"golden_squared":float(phi**2),
        "nearest_supergolden_squared":float(ratios[np.argmin(np.abs(ratios-psi**2))]),
        "nearest_golden_squared":float(ratios[np.argmin(np.abs(ratios-phi**2))]),
        "supergolden_absolute_mismatch":float(np.min(np.abs(ratios-psi**2))),
        "golden_absolute_mismatch":float(np.min(np.abs(ratios-phi**2))),
        "interpretation":"diagnostic only; proximity is not a theorem or evidence of universality"}

def q6_closed_two_cochain_projector() -> np.ndarray:
    """Orthogonal projector onto ker(B3.T) in the 240-dimensional 2-cochain space."""
    cube=build_unit_cube_complex(6)
    b3=cube.boundary[3].astype(float)
    # SVD with an explicit relative cutoff avoids promoting roundoff-sized
    # singular values of this rank-deficient integer incidence matrix.
    left, singular, _ = np.linalg.svd(b3, full_matrices=False)
    cutoff = max(b3.shape) * np.finfo(float).eps * singular[0] * 100.0
    rank = int(np.count_nonzero(singular > cutoff))
    if rank != 111:
        raise ArithmeticError(f"unexpected rank(B3)={rank}; expected 111 for the full six-cube")
    image_basis = left[:, :rank]
    p = np.eye(b3.shape[0]) - image_basis @ image_basis.T
    return 0.5*(p+p.T)

def q6_hypercube_cosine_law(fluxes: np.ndarray,w: float=2.) -> dict[str,object]:
    """Evaluate Δμ4=8 w^4 Σ(cos Ωp-1) and F_Q6=-Δμ4/(8w^4)."""
    omega=np.asarray(fluxes,dtype=float)
    if omega.shape!=(240,) or not np.all(np.isfinite(omega)): raise ValueError("fluxes must be finite length-240")
    if not np.isfinite(w) or w<=0: raise ValueError("w must be finite and positive")
    delta=float(8.*w**4*np.sum(np.cos(omega)-1.))
    action=float(np.sum(1.-np.cos(omega)))
    quadratic,quartic=float(.5*np.dot(omega,omega)),float(-np.sum(omega**4)/24.)
    return {"plaquette_count":240,"w":float(w),"delta_mu4":delta,"F_Q6":action,
        "F_Q6_from_delta_mu4":float(-delta/(8.*w**4)),"quadratic_term":quadratic,
        "quartic_term":quartic,"taylor_through_quartic":float(quadratic+quartic),
        "taylor_remainder":float(action-quadratic-quartic),
        "identity_residual":float(action+delta/(8.*w**4)),
        "status":"exact algebraic identity under the supplied plaquette-flux convention"}

def q6_projected_quadratic_flux_action(fluxes: np.ndarray) -> dict[str,object]:
    """Return 1/2 Ω.T P_cl Ω, the closed-sector quadratic diagnostic."""
    omega=np.asarray(fluxes,dtype=float)
    if omega.shape!=(240,) or not np.all(np.isfinite(omega)): raise ValueError("fluxes must be finite length-240")
    p=q6_closed_two_cochain_projector()
    return {"projector":p,"projector_idempotence_residual":float(np.linalg.norm(p@p-p)),
        "projector_symmetry_residual":float(np.linalg.norm(p-p.T)),"projected_flux":p@omega,
        "projected_quadratic_action":float(.5*omega@p@omega),"rank_projector":int(np.linalg.matrix_rank(p)),
        "status":"closed-2-cochain projection; distinct from the unprojected cosine action"}

def green_kubo_ou(relaxation: np.ndarray,covariance: np.ndarray,stress: np.ndarray,
                  temperature: float=1.,volume: float=1.,atol: float=1e-10) -> GreenKuboResult:
    """Integrate C(t)=b.T exp(-A t) Sigma b for symmetric stable OU relaxation."""
    a,sigma,b=np.asarray(relaxation,float),np.asarray(covariance,float),np.asarray(stress,float)
    if a.ndim!=2 or a.shape[0]!=a.shape[1]: raise ValueError("relaxation must be square")
    n=a.shape[0]
    if sigma.shape!=(n,n) or b.shape!=(n,): raise ValueError("input shapes do not match")
    if not (np.all(np.isfinite(a)) and np.all(np.isfinite(sigma)) and np.all(np.isfinite(b))):
        raise ValueError("inputs must be finite")
    if temperature<=0 or volume<=0 or atol<=0: raise ValueError("temperature, volume, atol must be positive")
    if not np.allclose(a,a.T,atol=atol) or not np.allclose(sigma,sigma.T,atol=atol):
        raise ValueError("relaxation and covariance must be symmetric")
    evals,vectors=np.linalg.eigh(.5*(a+a.T))
    if evals[0]<-atol: raise ValueError("relaxation must be positive semidefinite")
    if np.linalg.eigvalsh(.5*(sigma+sigma.T))[0]<-atol: raise ValueError("covariance must be positive semidefinite")
    amplitudes=(vectors.T@b)*(vectors.T@sigma@b)
    zero=evals<=atol
    zero_amplitude=float(np.sum(amplitudes[zero]))
    convergent=abs(zero_amplitude)<=atol
    modal=np.zeros_like(evals)
    modal[~zero]=amplitudes[~zero]/evals[~zero]
    integral=float(np.sum(modal)) if convergent else float("inf")
    return GreenKuboResult(float(integral/(temperature*volume)),integral,evals,amplitudes,
                           modal,zero_amplitude,bool(convergent))

def green_kubo_single_exponential(amplitude:float,rate:float,temperature:float=1.,volume:float=1.) -> float:
    """For C(t)=amplitude exp(-rate t), eta=amplitude/(rate kB T V)."""
    if rate<=0 or temperature<=0 or volume<=0: raise ValueError("rate, temperature, volume must be positive")
    return float(amplitude/(rate*temperature*volume))

def shear_stress_from_energy(energy_at_strain,strain_step:float=1e-5) -> float:
    """Central difference for Sigma=-dH/dgamma at gamma=0."""
    if strain_step<=0 or not np.isfinite(strain_step): raise ValueError("strain_step must be positive")
    return float(-(energy_at_strain(strain_step)-energy_at_strain(-strain_step))/(2.*strain_step))

def green_kubo_unitary_equivalence_residuals(unitary,h_source,stress_source,state_source,
                                              h_target,stress_target,state_target) -> dict[str,float]:
    """Audit U Hs U†=Ht, U Πs U†=Πt, U ρs U†=ρt."""
    u=np.asarray(unitary,complex)
    hs,ht,ps,pt,rs,rt=(np.asarray(x) for x in
        (h_source,h_target,stress_source,stress_target,state_source,state_target))
    if u.ndim!=2 or u.shape[0]!=u.shape[1]: raise ValueError("unitary must be square")
    n=u.shape[0]
    if any(x.shape!=(n,n) for x in (hs,ht,ps,pt,rs,rt)): raise ValueError("all matrices must match unitary")
    return {"unitarity":float(np.linalg.norm(u.conj().T@u-np.eye(n))),
        "hamiltonian_intertwining":float(np.linalg.norm(u@hs@u.conj().T-ht)),
        "stress_intertwining":float(np.linalg.norm(u@ps@u.conj().T-pt)),
        "state_intertwining":float(np.linalg.norm(u@rs@u.conj().T-rt))}

"""Exact Hessian identity for finite conditional sums of the compact Q6 Wilson action.

The routine evaluates the derivative of -log(sum_j a_j exp(-V(q,eta_j)))
for an explicitly supplied finite fibre ensemble. It is exact for that
finite sum, not a claim that a finite list exhausts the rank-124 lattice.
"""
from __future__ import annotations

import math
import numpy as np
from scipy.special import logsumexp


def compact_wilson_fiber_hessian(
    plaquette_mode_matrix,
    fibre_plaquette_angles,
    log_base_weights,
    beta_space: float,
) -> dict:
    """Compute the exact finite-fibre Hessian at q=0.

    Parameters
    ----------
    plaquette_mode_matrix:
        M = D1 B, shape (n_plaquettes, n_retained_modes).
    fibre_plaquette_angles:
        Rows eta_j = D1 U xi_j (or any explicitly specified fibre
        plaquette configurations), shape (n_fibre_points, n_plaquettes).
    log_base_weights:
        Logarithms of positive measure weights a_j, one per fibre point.
        These may include multiplicities and/or electric transfer weights.
    beta_space:
        Positive spatial Wilson coupling.

    For theta_j(q)=eta_j+M q and
    V_j(q)=beta sum_f (1-cos(theta_jf)), the exact identity is
      Hess(-log Z)|0 = E[V''|0] - Cov(V'|0,V'|0),
    with expectation under weights proportional to a_j exp(-V_j(0)).
    """
    beta = float(beta_space)
    m = np.asarray(plaquette_mode_matrix, dtype=float)
    eta = np.asarray(fibre_plaquette_angles, dtype=float)
    loga = np.asarray(log_base_weights, dtype=float)
    if not math.isfinite(beta) or beta <= 0:
        raise ValueError("beta_space must be finite and positive")
    if m.ndim != 2 or eta.ndim != 2 or eta.shape[1] != m.shape[0]:
        raise ValueError("incompatible plaquette mode matrix and fibre angles")
    if loga.shape != (eta.shape[0],) or eta.shape[0] < 1:
        raise ValueError("one finite log weight is required per fibre point")
    if not (np.all(np.isfinite(m)) and np.all(np.isfinite(eta))
            and np.all(np.isfinite(loga))):
        raise ValueError("all inputs must be finite")

    # Evaluate the compact Wilson potential and its first two q-derivatives.
    theta = eta
    sin_theta = np.sin(theta)
    cos_theta = np.cos(theta)
    v0 = beta * np.sum(1.0 - cos_theta, axis=1)
    grad = beta * np.einsum("jf,fi->ji", sin_theta, m)
    hess_each = beta * np.einsum("jf,fi,fk->jik", cos_theta, m, m)

    log_terms = loga - v0
    log_z = float(logsumexp(log_terms))
    probs = np.exp(log_terms - log_z)
    mean_grad = np.einsum("j,ji->i", probs, grad)
    mean_hess = np.einsum("j,jik->ik", probs, hess_each)
    centered = grad - mean_grad
    covariance = np.einsum("j,ji,jk->ik", probs, centered, centered)
    effective_hessian = mean_hess - covariance
    effective_hessian = 0.5 * (effective_hessian + effective_hessian.T)

    return {
        "audit": "Q6_COMPACT_WILSON_FINITE_FIBRE_HESSIAN",
        "beta_space": beta,
        "fibre_points": int(eta.shape[0]),
        "plaquettes": int(m.shape[0]),
        "retained_modes": int(m.shape[1]),
        "log_partition_at_zero": log_z,
        "probability_sum_residual": float(abs(probs.sum() - 1.0)),
        "mean_gradient": mean_grad.tolist(),
        "mean_hessian": mean_hess.tolist(),
        "gradient_covariance": covariance.tolist(),
        "effective_hessian": effective_hessian.tolist(),
        "eigenvalues": np.linalg.eigvalsh(effective_hessian).tolist(),
        "identity": "H_eff = E[V''|q=0] - Cov(V'|q=0,V'|q=0)",
        "scope": (
            "Exact derivative of the explicitly supplied finite conditional "
            "fibre sum. It is not the complete rank-124 integer Wilson marginal."
        ),
        "complete_fibre_sum_evaluated": False,
    }


def finite_fibre_log_partition(q, plaquette_mode_matrix,
                               fibre_plaquette_angles, log_base_weights,
                               beta_space: float) -> float:
    """Stable log partition for finite-difference validation of the identity."""
    q = np.asarray(q, dtype=float)
    m = np.asarray(plaquette_mode_matrix, dtype=float)
    eta = np.asarray(fibre_plaquette_angles, dtype=float)
    loga = np.asarray(log_base_weights, dtype=float)
    theta = eta + (m @ q)[None, :]
    v = float(beta_space) * np.sum(1.0 - np.cos(theta), axis=1)
    return float(logsumexp(loga - v))

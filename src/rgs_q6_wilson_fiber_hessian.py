"""Quadratic Wilson fibre marginalization onto the five Q6-D5 modes.

This computes the exact Gaussian (quadratic Wilson) Schur complement.
It does not claim to integrate the nonlinear compact Wilson theory.
"""
from __future__ import annotations

import math
import numpy as np
from scipy.linalg import null_space

from src.rgs_q6_equivariant_cochains import build_q6_cochain_complex
from src.rgs_q6_d5_phase_space import build_q6_d5_phase_space_sector


def q6_wilson_fiber_hessian_report(beta_space: float = 1.0,
                                  eigenvalue: float = 12.0) -> dict:
    """Integrate the 124 discarded Coulomb modes in the quadratic Wilson action.

    In Coulomb gauge C=ker(D0.T), split C=ran(B) orthogonally into the
    selected five-mode eigenspace B and its 124-dimensional complement U.
    For S2(x)=beta_space/2*x.T*L1*x, the fibre Gaussian integral is exact.
    The returned Schur complement is H_BB-H_BU H_UU^{-1} H_UB.
    """
    beta = float(beta_space)
    lam = float(eigenvalue)
    if not math.isfinite(beta) or beta <= 0 or not math.isfinite(lam) or lam <= 0:
        raise ValueError("beta_space and eigenvalue must be finite and positive")

    complex_data = build_q6_cochain_complex()
    d0 = np.asarray(complex_data["D0"], dtype=float)
    d1 = np.asarray(complex_data["D1"], dtype=float)
    l1 = d1.T @ d1

    sector = build_q6_d5_phase_space_sector(beta=1.0, eigenvalue=lam)
    b = np.asarray(sector["selected_mode_basis"], dtype=float)

    # Orthonormal basis of the 129-dimensional Coulomb tangent.
    _, singular, vh = np.linalg.svd(d0.T, full_matrices=True)
    rank = int(np.sum(singular > 1e-10))
    c = vh[rank:].T
    # Null space inside Coulomb gauge of the five retained coordinates.
    fibre_coordinates = null_space(b.T @ c, rcond=1e-11)
    u = c @ fibre_coordinates

    h = beta * l1
    hbb = b.T @ h @ b
    hbu = b.T @ h @ u
    huu = u.T @ h @ u
    hub = hbu.T
    schur = hbb - hbu @ np.linalg.solve(huu, hub)

    # Exact Gaussian integration over eta contributes only a q-independent
    # determinant. Its q-Hessian is therefore zero.
    sign, logdet_huu = np.linalg.slogdet(huu)
    if sign <= 0:
        raise ArithmeticError("discarded-mode Hessian is not positive definite")
    gaussian_log_normalizer = (
        0.5 * u.shape[1] * math.log(2.0 * math.pi) - 0.5 * logdet_huu
    )

    target = beta * lam * np.eye(5)
    checks = {
        "coulomb_dimension_129": c.shape == (192, 129),
        "fibre_dimension_124": u.shape == (192, 124),
        "retained_basis_orthonormal": np.linalg.norm(b.T @ b - np.eye(5)) < 1e-9,
        "fibre_basis_orthonormal": np.linalg.norm(u.T @ u - np.eye(124)) < 1e-9,
        "retained_fibre_orthogonal": np.linalg.norm(b.T @ u) < 1e-9,
        "both_in_coulomb_gauge": max(np.linalg.norm(d0.T @ b),
                                     np.linalg.norm(d0.T @ u)) < 1e-8,
        "fibre_hessian_positive": sign > 0,
        "cross_hessian_vanishes": np.linalg.norm(hbu) < 1e-8,
        "schur_equals_retained_hessian": np.linalg.norm(schur - hbb) < 1e-8,
        "retained_hessian_is_beta_lambda_identity":
            np.linalg.norm(schur - target) < 1e-8,
    }
    return {
        "audit": "Q6_WILSON_QUADRATIC_FIBRE_MARGINALIZATION",
        "beta_space": beta,
        "hodge_eigenvalue": lam,
        "coulomb_dimension": int(c.shape[1]),
        "retained_dimension": int(b.shape[1]),
        "integrated_fibre_dimension": int(u.shape[1]),
        "retained_hessian": hbb.tolist(),
        "cross_hessian_frobenius": float(np.linalg.norm(hbu)),
        "fibre_hessian_min_eigenvalue":
            float(np.linalg.eigvalsh(huu).min()),
        "fibre_hessian_logdet": float(logdet_huu),
        "gaussian_log_normalizer_q_independent": float(gaussian_log_normalizer),
        "schur_complement": schur.tolist(),
        "expected_schur_complement": target.tolist(),
        "schur_residual": float(np.linalg.norm(schur - target)),
        "effective_quadratic_action":
            "S_eff^(2)(q) = (beta_space*lambda/2)*q.T@q + constant",
        "checks": checks,
        "all_checks_passed": bool(all(checks.values())),
        "proved_within_scope": [
            "the 129-dimensional Coulomb tangent splits orthogonally into 5 retained and 124 discarded real modes",
            "the quadratic Wilson Hessian has zero retained-discarded cross block",
            "Gaussian integration of the 124 discarded modes changes only a q-independent determinant",
            "the retained Maxwell quadratic coefficient is exactly beta_space*lambda",
        ],
        "not_proved": [
            "integration of the full nonlinear compact Wilson cosine action",
            "a canonical compact integer-fibre measure on the 124-dimensional kernel",
            "a Wilson-derived Hua cutoff or Bergman normalization",
            "a microscopic prediction of the U(1) charge or fine-structure constant",
            "a 3+1-dimensional continuum limit",
        ],
        "scope": (
            "Exact finite-dimensional Gaussian marginalization of the quadratic "
            "spatial Wilson action on the Q6 Coulomb tangent. This is not the "
            "full finite-beta compact Wilson path integral."
        ),
    }

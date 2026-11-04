"""Tensorially correct interface from the RGS fourth moment to Maxwell response.

The fully symmetric rank-four residual cannot couple directly to two
antisymmetric electromagnetic field strengths.  It first acts on an
independent STF(3) order parameter Q, producing a symmetric-traceless
constitutive anisotropy A(Q)=R:Q.  The resulting Maxwell response is a
candidate material constitutive law, not a derived alpha prediction.
"""
from __future__ import annotations

import numpy as np


def _validate_stf(q: np.ndarray, tol: float = 1e-10) -> np.ndarray:
    q = np.asarray(q, dtype=float)
    if q.shape != (3, 3) or not np.all(np.isfinite(q)):
        raise ValueError("Q must be a finite 3x3 matrix")
    if not np.allclose(q, q.T, atol=tol, rtol=0.0):
        raise ValueError("Q must be symmetric")
    if abs(float(np.trace(q))) > tol:
        raise ValueError("Q must be trace-free")
    return q


def fourth_moment_anisotropy(residual: np.ndarray) -> np.ndarray:
    """Validate R_ijkl and return its induced STF operator matrix."""
    r = np.asarray(residual, dtype=float)
    if r.shape != (3, 3, 3, 3) or not np.all(np.isfinite(r)):
        raise ValueError("R must be a finite rank-four 3x3x3x3 tensor")
    permutations = ((1, 0, 2, 3), (0, 2, 1, 3), (0, 1, 3, 2),
                    (3, 2, 1, 0))
    if any(not np.allclose(r, r.transpose(p), atol=1e-9, rtol=0.0)
           for p in permutations):
        raise ValueError("R must be fully symmetric")
    double_trace = np.einsum("iikl->kl", r)
    if not np.allclose(double_trace, 0.0, atol=1e-8, rtol=0.0):
        raise ValueError("R must have vanishing pair contraction (STF rank four)")
    return r


def maxwell_constitutive_response(
    residual: np.ndarray,
    order_parameter: np.ndarray,
    electric_baseline: float = 1.0,
    magnetic_baseline: float = 1.0,
    electric_coupling: float = 0.0,
    magnetic_coupling: float = 0.0,
) -> dict:
    """Build positive electric/magnetic response matrices from R:Q.

    L_EM = 1/2 E_i (Z_E)_ij E_j - 1/2 B_i (Z_B)_ij B_j.
    Positivity is checked, not silently repaired by clipping eigenvalues.
    """
    r = fourth_moment_anisotropy(residual)
    q = _validate_stf(order_parameter)
    ze = float(electric_baseline)
    zb = float(magnetic_baseline)
    ge = float(electric_coupling)
    gb = float(magnetic_coupling)
    if not np.all(np.isfinite([ze, zb, ge, gb])) or ze <= 0 or zb <= 0:
        raise ValueError("baselines must be positive and all coefficients finite")

    a = np.einsum("ijkl,kl->ij", r, q)
    a = 0.5 * (a + a.T)
    if abs(float(np.trace(a))) > 1e-8:
        raise ArithmeticError("R:Q failed to land in Sym_0(3)")
    z_e = ze * np.eye(3) + ge * a
    z_b = zb * np.eye(3) + gb * a
    eig_e = np.linalg.eigvalsh(z_e)
    eig_b = np.linalg.eigvalsh(z_b)
    if eig_e[0] <= 0 or eig_b[0] <= 0:
        raise ValueError("constitutive matrices are not positive definite")

    return {
        "anisotropy_operator": a,
        "electric_matrix": z_e,
        "magnetic_matrix": z_b,
        "electric_eigenvalues": eig_e,
        "magnetic_eigenvalues": eig_b,
        "electric_coupling_bound": (
            float(ze / np.linalg.norm(a, 2)) if np.linalg.norm(a, 2) else np.inf
        ),
        "magnetic_coupling_bound": (
            float(zb / np.linalg.norm(a, 2)) if np.linalg.norm(a, 2) else np.inf
        ),
        "status": "POSITIVE_CONSTITUTIVE_CANDIDATE",
    }


def maxwell_density(electric, magnetic, electric_matrix, magnetic_matrix) -> float:
    """Evaluate the real 3+1 Maxwell density in the chosen constitutive frame."""
    e = np.asarray(electric, dtype=float)
    b = np.asarray(magnetic, dtype=float)
    ze = np.asarray(electric_matrix, dtype=float)
    zb = np.asarray(magnetic_matrix, dtype=float)
    if e.shape != (3,) or b.shape != (3,) or ze.shape != (3, 3) or zb.shape != (3, 3):
        raise ValueError("E,B must be 3-vectors and response matrices 3x3")
    return float(0.5 * e @ ze @ e - 0.5 * b @ zb @ b)


def direct_antisymmetric_contraction(residual: np.ndarray, field_strength: np.ndarray) -> float:
    """Return R_ijkl F_ij F_kl; it vanishes for symmetric R and antisymmetric F."""
    r = fourth_moment_anisotropy(residual)
    f = np.asarray(field_strength, dtype=float)
    if f.shape != (3, 3) or not np.allclose(f, -f.T, atol=1e-10, rtol=0.0):
        raise ValueError("F must be a 3x3 antisymmetric tensor")
    return float(np.einsum("ijkl,ij,kl", r, f, f))


def audit_fourth_moment_maxwell_bridge(tolerance: float = 1e-9) -> dict:
    """Audit the real 64-site RGS seed and the tensorial Maxwell interface."""
    from src.rgs_ellis_refinement import load_rgs_seed_data
    from src.rgs_fourth_moment_bridge import fourth_moment_decomposition
    from src.rgs_q6_oh_action import signed_permutation_matrices

    seed = load_rgs_seed_data()
    parts = fourth_moment_decomposition(seed.coordinates)
    r = fourth_moment_anisotropy(parts["M4_anisotropic"])
    q = np.diag([1.0, -1.0, 0.0]) / np.sqrt(2.0)
    a = np.einsum("ijkl,kl->ij", r, q)
    norm_a = float(np.linalg.norm(a, 2))
    coupling = 0.25 / max(norm_a, 1.0)
    response = maxwell_constitutive_response(
        r, q, electric_coupling=coupling, magnetic_coupling=-coupling
    )

    covariance_residuals = []
    for g in signed_permutation_matrices():
        rg = np.einsum("ia,jb,kc,ld,abcd->ijkl", g, g, g, g, r)
        qg = g @ q @ g.T
        ag = np.einsum("ijkl,kl->ij", rg, qg)
        covariance_residuals.append(float(np.linalg.norm(ag - g @ a @ g.T)))

    f = np.array([[0.0, 0.7, -0.2], [-0.7, 0.0, 0.4], [0.2, -0.4, 0.0]])
    no_go = direct_antisymmetric_contraction(r, f)
    report = {
        "audit": "RGS_FOURTH_MOMENT_MAXWELL_CONSTITUTIVE_BRIDGE",
        "seed_shape": list(seed.coordinates.shape),
        "M2_equals_80I_residual": float(
            np.linalg.norm(parts["M2"] - 80.0 * np.eye(3))
        ),
        "fourth_moment_anisotropy_norm": float(np.linalg.norm(r)),
        "stf_order_parameter": q.tolist(),
        "induced_stf_response": a.tolist(),
        "electric_eigenvalues": response["electric_eigenvalues"].tolist(),
        "magnetic_eigenvalues": response["magnetic_eigenvalues"].tolist(),
        "electric_coupling_bound": response["electric_coupling_bound"],
        "magnetic_coupling_bound": response["magnetic_coupling_bound"],
        "Oh_group_order": len(covariance_residuals),
        "Oh_covariance_max_residual": max(covariance_residuals),
        "direct_R_F_F_contraction": no_go,
        "checks": {
            "positive_electric_response": bool(response["electric_eigenvalues"][0] > 0),
            "positive_magnetic_response": bool(response["magnetic_eigenvalues"][0] > 0),
            "Oh_covariance": max(covariance_residuals) < tolerance,
            "direct_antisymmetric_contraction_no_go": abs(no_go) < tolerance,
            "induced_response_is_tracefree": abs(float(np.trace(a))) < tolerance,
        },
        "status": {
            "64_site_fourth_moment": "NUMERICALLY_VERIFIED_BY_EXISTING_AUDIT",
            "R_to_STF_response": "EXACT_TENSOR_CONSTRUCTION",
            "Maxwell_constitutive_law": "CANDIDATE",
            "Q6_A4_restricted_intertwiner_matrix": "OPEN",
            "Wilson_derived_order_parameter_measure": "OPEN",
            "scalar_Maxwell_normalization_and_alpha": "NOT_DERIVED",
        },
    }
    report["pass"] = bool(all(report["checks"].values()))
    return report

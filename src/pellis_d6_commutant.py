"""Full D6 commutant for the three-mode factor of the Fibonacci six-orbit.

D6 covariance alone permits two independent mode sectors: the invariant
axis (dimension one) and its orthogonal plane (dimension two). The special
M_parallel == M_perp case is the stronger threefold-degenerate ansatz.
"""
from __future__ import annotations
import numpy as np
from src.fibonacci_d6_channel_intertwiner import d6_channel_intertwiner


def mode_isotypic_projectors() -> tuple[np.ndarray, np.ndarray]:
    """Project onto span(1,1,1) and its orthogonal two-dimensional plane."""
    u = np.ones(3, dtype=float) / np.sqrt(3.0)
    p_axis = np.outer(u, u)
    return p_axis.astype(complex), (np.eye(3) - p_axis).astype(complex)


def d6_covariant_radial_transfer(parallel: np.ndarray, transverse: np.ndarray) -> np.ndarray:
    """Build M_parallel⊗I_sheet⊗P_axis + M_transverse⊗I_sheet⊗P_plane."""
    mp, mt = np.asarray(parallel, dtype=complex), np.asarray(transverse, dtype=complex)
    if mp.shape != (2, 2) or mt.shape != (2, 2):
        raise ValueError("parallel and transverse matrices must both be 2x2")
    if not np.all(np.isfinite(mp)) or not np.all(np.isfinite(mt)):
        raise ValueError("transfer matrices must be finite")
    p0, p2 = mode_isotypic_projectors()
    eye_sheet = np.eye(2, dtype=complex)
    return np.kron(mp, np.kron(eye_sheet, p0)) + np.kron(mt, np.kron(eye_sheet, p2))


def audit_d6_commutant(parallel: np.ndarray, transverse: np.ndarray, tolerance: float = 1e-12) -> dict[str, float | bool]:
    """Check D6 covariance and whether the stronger triplet degeneracy holds."""
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("tolerance must be positive and finite")
    data = d6_channel_intertwiner()
    transfer = d6_covariant_radial_transfer(parallel, transverse)
    r = np.kron(np.eye(2, dtype=complex), data["r_channel"])
    s = np.kron(np.eye(2, dtype=complex), data["s_channel"])
    residuals = {
        "rotation_commutator": float(np.linalg.norm(transfer @ r - r @ transfer)),
        "reflection_commutator": float(np.linalg.norm(transfer @ s - s @ transfer)),
    }
    mp, mt = np.asarray(parallel, dtype=complex), np.asarray(transverse, dtype=complex)
    return {
        **residuals,
        "triplet_degenerate": bool(np.linalg.norm(mp - mt) <= tolerance),
        "max_residual": max(residuals.values()),
        "passed": max(residuals.values()) <= tolerance,
    }

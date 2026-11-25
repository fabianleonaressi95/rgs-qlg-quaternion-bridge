"""Finite operator-level audit for the conditional D6 x two-sheet radial lift.

This module proves algebraic covariance of the declared tensor lift
M(E) ⊗ I_6 and parity boundary operator d(E)P_+ + b(E)P_-.
It does not assert that a physical Pellis/GSR model has these channels.
"""
from __future__ import annotations

import numpy as np

from src.fibonacci_d6_channel_intertwiner import d6_channel_intertwiner


def sheet_parity_projectors() -> tuple[np.ndarray, np.ndarray]:
    """Orthogonal projectors for the central D6 half-turn r^3."""
    sigma_x = np.array([[0.0, 1.0], [1.0, 0.0]])
    exchange = np.kron(sigma_x, np.eye(3))
    eye = np.eye(6)
    return (eye + exchange) / 2, (eye - exchange) / 2


def lifted_radial_transfer(matrix: np.ndarray) -> np.ndarray:
    """Lift a 2x2 Cauchy transfer identically over six channels."""
    m = np.asarray(matrix, dtype=complex)
    if m.shape != (2, 2) or not np.all(np.isfinite(m)):
        raise ValueError("matrix must be a finite 2x2 array")
    return np.kron(m, np.eye(6, dtype=complex))


def parity_boundary_operator(b: complex, d: complex) -> np.ndarray:
    """Return B = d P_+ + b P_- on sheet ⊗ mode channels."""
    if not np.isfinite(b) or not np.isfinite(d):
        raise ValueError("b and d must be finite")
    p_plus, p_minus = sheet_parity_projectors()
    return d * p_plus + b * p_minus


def audit_d6_radial_lift(matrix: np.ndarray, b: complex, d: complex) -> dict[str, float | int | bool]:
    """Check D6 covariance and the threefold parity-kernel prediction.

    The returned kernel dimensions are algebraic dimensions of the finite
    channel boundary matrix, not multiplicities of a physical spectrum.
    """
    data = d6_channel_intertwiner()
    r, s = data["r_channel"], data["s_channel"]
    transfer = lifted_radial_transfer(matrix)
    boundary = parity_boundary_operator(b, d)
    eye2 = np.eye(2, dtype=complex)
    r_lift, s_lift = np.kron(eye2, r), np.kron(eye2, s)
    residuals = {
        "transfer_rotation_covariance": float(np.linalg.norm(transfer @ r_lift - r_lift @ transfer)),
        "transfer_reflection_covariance": float(np.linalg.norm(transfer @ s_lift - s_lift @ transfer)),
        "boundary_rotation_covariance": float(np.linalg.norm(boundary @ r - r @ boundary)),
        "boundary_reflection_covariance": float(np.linalg.norm(boundary @ s - s @ boundary)),
        "boundary_self_adjoint_when_real": float(
            np.linalg.norm(boundary - boundary.conj().T)
        ) if np.isreal(b) and np.isreal(d) else 0.0,
    }
    return {
        **residuals,
        "symmetric_kernel_dimension": int(3 if d == 0 and b != 0 else (6 if d == 0 and b == 0 else 0)),
        "antisymmetric_kernel_dimension": int(3 if b == 0 and d != 0 else (6 if d == 0 and b == 0 else 0)),
        "max_residual": max(residuals.values()),
        "passed": max(residuals.values()) <= 1e-12,
    }

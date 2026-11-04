"""Discrete frame alignment and Spin(4)/Clifford transport on Q6.

This module defines a *candidate edgewise frame-alignment connection*.
It is not the Levi-Civita connection of a varying metric and is not a
coframe-derived TEGR spin connection unless the edgewise metric-compatibility
conditions documented below are satisfied.

Convention: coframes are 4x4 matrices with internal index first,
g_v = e_v.T @ e_v (Euclidean signature), and
R_{v<-u} is defined by the right polar factor of e_v @ inv(e_u).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import expm, logm

from src.q6_complex import Q6Complex, build_q6_complex
from src.unified_discrete_action import UnifiedQ6Fields


_PAULI = (
    np.array([[0, 1], [1, 0]], dtype=complex),
    np.array([[0, -1j], [1j, 0]], dtype=complex),
    np.array([[1, 0], [0, -1]], dtype=complex),
)
_I2 = np.eye(2, dtype=complex)


def euclidean_gamma_matrices() -> tuple[np.ndarray, ...]:
    """Return Hermitian 4x4 Euclidean gamma matrices with {gamma_a,gamma_b}=2δab."""
    s1, s2, s3 = _PAULI
    return (
        np.kron(s1, _I2),
        np.kron(s2, _I2),
        np.kron(s3, s1),
        np.kron(s3, s2),
    )


@dataclass(frozen=True)
class EdgeSpinTransport:
    source: int
    target: int
    vector_rotation: np.ndarray
    spin_transport: np.ndarray
    metric_source: np.ndarray
    metric_target: np.ndarray
    metric_mismatch: float
    coframe_alignment_residual: float
    logarithm_residual: float
    clifford_residual: float


def _validate_coframe(e: np.ndarray) -> np.ndarray:
    e = np.asarray(e, dtype=float)
    if e.shape != (4, 4) or not np.all(np.isfinite(e)):
        raise ValueError("each coframe must be a finite real 4x4 matrix")
    det = float(np.linalg.det(e))
    if abs(det) <= 1e-12:
        raise ValueError("coframe must be non-degenerate")
    if det < 0.0:
        raise ValueError(
            "Spin(4) transport requires a fixed orientation sector: det(e)>0"
        )
    return e


def _gamma_spin_lift(rotation: np.ndarray) -> tuple[np.ndarray, float, float]:
    """Lift R in SO(4) to Spin(4), using S gamma_c S† = R_bc gamma_b."""
    r = np.asarray(rotation, dtype=float)
    if r.shape != (4, 4):
        raise ValueError("rotation must be 4x4")
    if not np.allclose(r.T @ r, np.eye(4), atol=2e-9, rtol=2e-9):
        raise ValueError("vector transport is not orthogonal")
    if not np.isclose(np.linalg.det(r), 1.0, atol=2e-9):
        raise ValueError("vector transport must lie in SO(4)")

    omega_complex = logm(r.astype(complex))
    imag_residual = float(np.linalg.norm(omega_complex.imag))
    omega = omega_complex.real
    skew_residual = float(np.linalg.norm(omega + omega.T))
    if imag_residual > 2e-8 or skew_residual > 2e-7:
        raise ValueError("SO(4) logarithm is not a resolved real skew generator")
    omega = 0.5 * (omega - omega.T)

    gamma = euclidean_gamma_matrices()
    generator = np.zeros((4, 4), dtype=complex)
    for a in range(4):
        for b in range(4):
            generator += 0.25 * omega[a, b] * (gamma[a] @ gamma[b])
    spin = expm(generator)

    # The spin representation must be unitary and intertwine the vector rotation.
    if not np.allclose(spin.conj().T @ spin, np.eye(4), atol=2e-8, rtol=2e-8):
        raise ArithmeticError("Spin(4) lift is not unitary")
    residual = 0.0
    for c in range(4):
        target = sum((r[b, c] * gamma[b] for b in range(4)),
                     np.zeros((4, 4), dtype=complex))
        residual = max(
            residual,
            float(np.linalg.norm(spin @ gamma[c] @ spin.conj().T - target)),
        )
    if residual > 3e-7:
        raise ArithmeticError(
            f"Spin/Clifford intertwining failed (residual={residual:.3e})"
        )
    return spin, residual, float(np.linalg.norm(expm(omega) - r))


def edge_frame_spin_transport(
    coframe_source: np.ndarray,
    coframe_target: np.ndarray,
    source: int = 0,
    target: int = 1,
) -> EdgeSpinTransport:
    """Construct polar frame alignment and its Spin(4) lift.

    Let A=e_target @ inv(e_source)=R P be the right polar decomposition.
    R is the closest orientation-preserving orthogonal factor in Frobenius
    norm. This is a frame-alignment prescription, not a metric connection
    for arbitrary independently varying metrics.

    If g_target == g_source, then A is exactly orthogonal and R=A, so this
    prescription is an exact local-frame rotation on that edge.
    """
    e0 = _validate_coframe(coframe_source)
    e1 = _validate_coframe(coframe_target)
    relative = e1 @ np.linalg.inv(e0)
    left, singular_values, right_t = np.linalg.svd(relative)
    rotation = left @ right_t
    if np.linalg.det(rotation) < 0.0:
        raise ValueError("coframe pair induces an orientation-reversing alignment")
    if np.min(singular_values) <= 1e-12:
        raise ValueError("relative coframe map is singular")

    spin, clifford_residual, log_residual = _gamma_spin_lift(rotation)
    g0 = e0.T @ e0
    g1 = e1.T @ e1
    mismatch = float(np.linalg.norm(g1 - g0))
    alignment = float(np.linalg.norm(e1 - rotation @ e0))
    return EdgeSpinTransport(
        source=int(source),
        target=int(target),
        vector_rotation=rotation,
        spin_transport=spin,
        metric_source=g0,
        metric_target=g1,
        metric_mismatch=mismatch,
        coframe_alignment_residual=alignment,
        logarithm_residual=log_residual,
        clifford_residual=clifford_residual,
    )


def build_q6_spin_transports(
    fields: UnifiedQ6Fields,
    complex_: Q6Complex | None = None,
) -> dict[tuple[int, int], EdgeSpinTransport]:
    """Build candidate edgewise Spin(4) transports for all canonical Q6 edges."""
    q6 = build_q6_complex() if complex_ is None else complex_
    f = fields.validated(len(q6.vertices))
    return {
        (u, v): edge_frame_spin_transport(f.coframe[u], f.coframe[v], u, v)
        for u, v in q6.edges
    }


def clifford_algebra_residual() -> float:
    """Maximum Frobenius residual in the Euclidean Clifford anticommutators."""
    gamma = euclidean_gamma_matrices()
    identity = np.eye(4, dtype=complex)
    residual = 0.0
    for a in range(4):
        for b in range(4):
            expected = 2.0 * identity if a == b else np.zeros((4, 4), complex)
            residual = max(
                residual,
                float(np.linalg.norm(gamma[a] @ gamma[b] + gamma[b] @ gamma[a] - expected)),
            )
    return residual


def reverse_spin_transport(transport: EdgeSpinTransport) -> np.ndarray:
    """Return the inverse edge transport; for unitary Spin(4), this is S†."""
    return transport.spin_transport.conj().T

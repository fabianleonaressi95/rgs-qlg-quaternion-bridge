"""Exact D6 channel factorization of a reversible Fibonacci six-orbit.

The orbit basis is e_k = e_{T^k p}, k=0,...,5.  The target basis is
|sheet,mode>, ordered sheet-major with mode j=0,1,2.  This is a finite
representation-theoretic intertwiner; it is not, by itself, a physical
Pellis transfer operator.
"""
from __future__ import annotations

import numpy as np


def d6_channel_intertwiner() -> dict[str, np.ndarray]:
    """Return the exact permutation intertwiner C^6 -> C^2_sheet tensor C^3_mode.

    The index map is
      k=2j       -> (sheet=0, mode=j),
      k=2j+3 mod6 -> (sheet=1, mode=j).
    Thus the half-turn r^3 is sigma_x tensor I_3 and r^2 is I_2 tensor C_3.
    Matrices act on column vectors.
    """
    r = np.zeros((6, 6), dtype=complex)
    s = np.zeros((6, 6), dtype=complex)
    phi = np.zeros((6, 6), dtype=complex)
    for k in range(6):
        r[(k + 1) % 6, k] = 1.0
        s[(-k) % 6, k] = 1.0
        sheet = k % 2
        mode = k // 2 if sheet == 0 else ((k - 3) // 2) % 3
        phi[3 * sheet + mode, k] = 1.0

    c3 = np.zeros((3, 3), dtype=complex)
    f3 = np.zeros((3, 3), dtype=complex)
    for j in range(3):
        c3[(j + 1) % 3, j] = 1.0
        f3[(-j) % 3, j] = 1.0

    zero = np.zeros((3, 3), dtype=complex)
    r_channel = np.block([[zero, c3 @ c3], [c3.conj().T, zero]])
    s_channel = np.block([[f3, zero], [zero, f3]])
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
    half_turn_channel = np.kron(sigma_x, np.eye(3, dtype=complex))
    mode_rotation_channel = np.kron(np.eye(2, dtype=complex), c3)

    return {
        "Phi": phi,
        "r_orbit": r,
        "s_orbit": s,
        "r_channel": r_channel,
        "s_channel": s_channel,
        "C3": c3,
        "F3": f3,
        "S_channel": half_turn_channel,
        "r2_channel": mode_rotation_channel,
    }


def audit_d6_channel_intertwiner(tolerance: float = 1e-12) -> dict[str, float | bool]:
    """Check unitarity, D6 relations, and all claimed channel factorizations."""
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("tolerance must be positive and finite")
    data = d6_channel_intertwiner()
    phi = data["Phi"]
    r = data["r_orbit"]
    s = data["s_orbit"]
    rc = data["r_channel"]
    sc = data["s_channel"]
    eye6 = np.eye(6, dtype=complex)
    eye3 = np.eye(3, dtype=complex)
    r2 = rc @ rc
    r3 = r2 @ rc
    residuals = {
        "unitarity": float(np.linalg.norm(phi.conj().T @ phi - eye6, ord="fro")),
        "rotation_intertwining": float(np.linalg.norm(phi @ r - rc @ phi, ord="fro")),
        "reflection_intertwining": float(np.linalg.norm(phi @ s - sc @ phi, ord="fro")),
        "r6_identity": float(np.linalg.norm(np.linalg.matrix_power(rc, 6) - eye6, ord="fro")),
        "s2_identity": float(np.linalg.norm(sc @ sc - eye6, ord="fro")),
        "srs_relation": float(np.linalg.norm(sc @ rc @ sc - rc.conj().T, ord="fro")),
        "half_turn_factorization": float(np.linalg.norm(r3 - data["S_channel"], ord="fro")),
        "three_mode_factorization": float(np.linalg.norm(r2 - data["r2_channel"], ord="fro")),
        "reflection_factorization": float(np.linalg.norm(sc - np.kron(np.eye(2), data["F3"]), ord="fro")),
        "c3_order": float(np.linalg.norm(np.linalg.matrix_power(data["C3"], 3) - eye3, ord="fro")),
    }
    return {
        **residuals,
        "max_residual": max(residuals.values()),
        "passed": max(residuals.values()) <= tolerance,
    }


def triangle_to_dipole_intertwiner() -> dict[str, np.ndarray]:
    """Identify the three-point permutation module with the l=1 O(3) module.

    This is a candidate geometric realization of each three-mode sheet as
    the dipole (l=1) angular sector of a spherical resonator. It does not
    assert that Pellis's GSR operator selects this sector.
    """
    W = np.array([
        [1 / np.sqrt(3), 1 / np.sqrt(3), 1 / np.sqrt(3)],
        [2 / np.sqrt(6), -1 / np.sqrt(6), -1 / np.sqrt(6)],
        [0.0, 1 / np.sqrt(2), -1 / np.sqrt(2)],
    ], dtype=float)
    theta = 2 * np.pi / 3
    rotation = np.array([
        [1.0, 0.0, 0.0],
        [0.0, np.cos(theta), -np.sin(theta)],
        [0.0, np.sin(theta), np.cos(theta)],
    ])
    reflection = np.diag([1.0, 1.0, -1.0])
    c3 = np.array([
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ])
    f3 = np.array([
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 1.0, 0.0],
    ])
    return {
        "W3": W,
        "C3": c3,
        "F3": f3,
        "dipole_rotation": rotation,
        "dipole_reflection": reflection,
    }


def audit_triangle_to_dipole_intertwiner(tolerance: float = 1e-12) -> dict[str, float | bool]:
    """Audit the real orthogonal equivalence with the spherical l=1 module."""
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("tolerance must be positive and finite")
    data = triangle_to_dipole_intertwiner()
    w = data["W3"]
    residuals = {
        "orthogonality": float(np.linalg.norm(w @ w.T - np.eye(3), ord="fro")),
        "rotation_intertwining": float(np.linalg.norm(
            w @ data["C3"] @ w.T - data["dipole_rotation"], ord="fro")),
        "reflection_intertwining": float(np.linalg.norm(
            w @ data["F3"] @ w.T - data["dipole_reflection"], ord="fro")),
    }
    return {
        **residuals,
        "max_residual": max(residuals.values()),
        "passed": max(residuals.values()) <= tolerance,
    }

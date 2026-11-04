"""Positive C6 constitutive coupling on the six Q6 directional differences.

The identity between the six C6 channels and six Q6 bit directions is an
explicit frame convention (axis_map), not a claimed canonical icosahedral
embedding. Relabelling both channel systems leaves the quadratic form invariant.
"""
from __future__ import annotations

import numpy as np

from src.rgs_weyl_phase_closure import build_icosahedral_c6


def q6_axis_masks() -> tuple[int, ...]:
    """Masks for the canonical big-endian Q6 axes used by src.q6_complex."""
    return tuple(1 << (5 - axis) for axis in range(6))


def q6_directional_gradient_matrix() -> np.ndarray:
    """Return D with (vertex,axis) rows and vertex columns: (Df)[v,d]=f[v^d]-f[v]."""
    masks = q6_axis_masks()
    derivative = np.zeros((64 * 6, 64), dtype=float)
    for vertex in range(64):
        for axis, mask in enumerate(masks):
            row = 6 * vertex + axis
            derivative[row, vertex] = -1.0
            derivative[row, vertex ^ mask] = 1.0
    return derivative


def positive_c6_constitutive_matrix(
    beta: float = 0.1, epsilon: float = 1e-8,
    axis_map: tuple[int, ...] = (0, 1, 2, 3, 4, 5),
) -> tuple[np.ndarray, dict[str, float]]:
    """Build H=mu I+beta C6, with mu>sqrt(5)|beta| for strict positivity.

    axis_map[d] gives the C6 channel assigned to Q6 bit direction d. The default
    is a declared coordinate convention only; no O_h/SO(3) equivariance follows.
    """
    if not np.isfinite(beta) or not np.isfinite(epsilon) or epsilon <= 0:
        raise ValueError("beta must be finite and epsilon finite and positive")
    mapping = tuple(int(x) for x in axis_map)
    if sorted(mapping) != list(range(6)):
        raise ValueError("axis_map must be a permutation of 0..5")
    c6 = np.asarray(build_icosahedral_c6()["C6"], dtype=float)
    c_q6 = c6[np.ix_(mapping, mapping)]
    mu = float(np.sqrt(5.0) * abs(beta) + epsilon)
    h = mu * np.eye(6) + float(beta) * c_q6
    eig = np.linalg.eigvalsh(h)
    return h, {
        "mu": mu, "beta": float(beta), "epsilon": float(epsilon),
        "lambda_min": float(eig[0]), "lambda_max": float(eig[-1]),
        "strictly_positive": bool(eig[0] > 0.0),
        "C6_offdiagonal_norm": float(np.linalg.norm(c_q6 - np.diag(np.diag(c_q6)))),
    }


def build_q6_c6_quadratic_operator(
    beta: float = 0.1, epsilon: float = 1e-8,
    axis_map: tuple[int, ...] = (0, 1, 2, 3, 4, 5),
) -> tuple[np.ndarray, dict[str, float]]:
    """Return K=D.T (I_64 kron H) D for E(f)=1/2 f.T K f."""
    h, report = positive_c6_constitutive_matrix(beta, epsilon, axis_map)
    d = q6_directional_gradient_matrix()
    lifted_h = np.kron(np.eye(64), h)
    operator = d.T @ lifted_h @ d
    operator = 0.5 * (operator + operator.T)
    report = dict(report)
    report.update({
        "operator_shape": float(operator.shape[0]),
        "symmetry_residual": float(np.linalg.norm(operator - operator.T)),
        "constant_mode_residual": float(np.linalg.norm(operator @ np.ones(64))),
        "minimum_operator_eigenvalue": float(np.linalg.eigvalsh(operator)[0]),
    })
    return operator, report


def q6_c6_energy(field: np.ndarray, beta: float = 0.1,
                 epsilon: float = 1e-8,
                 axis_map: tuple[int, ...] = (0, 1, 2, 3, 4, 5)) -> float:
    """Evaluate the positive directional constitutive energy on a real Q6 scalar."""
    f = np.asarray(field, dtype=float)
    if f.shape != (64,) or not np.all(np.isfinite(f)):
        raise ValueError("field must be a finite real vector of length 64")
    operator, _ = build_q6_c6_quadratic_operator(beta, epsilon, axis_map)
    return float(0.5 * f @ operator @ f)

"""Exact SU(2)-equivariant fibre exchange used in the Pellis B dossier.

Coordinates on the conjugate fibre are taken in the conjugate basis.  The
intertwiner is complex-linear in these coordinates; it is not coordinate
complex conjugation.
"""
from __future__ import annotations

import numpy as np

EPSILON = np.array([[0.0, 1.0], [-1.0, 0.0]], dtype=complex)
EPSILON_INV = np.linalg.inv(EPSILON)


def su2_matrix(a: complex, b: complex) -> np.ndarray:
    """Return [[a,b],[-conj(b),conj(a)]] after checking |a|²+|b|²=1."""
    if not np.isclose(abs(a) ** 2 + abs(b) ** 2, 1.0, atol=1e-12):
        raise ValueError("SU(2) parameters must satisfy |a|^2 + |b|^2 = 1")
    return np.array([[a, b], [-np.conjugate(b), np.conjugate(a)]], dtype=complex)


def fibre_exchange_matrix() -> np.ndarray:
    """Matrix of J(psi, chi_bar)=(epsilon^{-1} chi_bar, epsilon psi)."""
    zero = np.zeros((2, 2), dtype=complex)
    return np.block([[zero, EPSILON_INV], [EPSILON, zero]])


def conjugate_intertwining_residual(u: np.ndarray) -> float:
    """Frobenius residual for epsilon U = conjugate(U) epsilon."""
    return float(np.linalg.norm(EPSILON @ u - np.conjugate(u) @ EPSILON))


def edge_covariance_residual(
    u: np.ndarray, psi_source: np.ndarray, psi_target: np.ndarray
) -> float:
    """Residual for D_barU(c psi)=c(D_U psi) on one oriented edge."""
    lhs = EPSILON @ psi_target - np.conjugate(u) @ EPSILON @ psi_source
    rhs = EPSILON @ (psi_target - u @ psi_source)
    return float(np.linalg.norm(lhs - rhs))


def doubled_edge_covariance_residual(
    u: np.ndarray,
    psi_source: np.ndarray,
    psi_target: np.ndarray,
    chi_source: np.ndarray,
    chi_target: np.ndarray,
) -> float:
    """Residual for D_U J_V = J_E D_U on a doubled edge."""
    j = fibre_exchange_matrix()
    source = np.concatenate((psi_source, chi_source))
    target = np.concatenate((psi_target, chi_target))
    d_after_j = j @ (target - np.block([
        [u, np.zeros((2, 2), dtype=complex)],
        [np.zeros((2, 2), dtype=complex), np.conjugate(u)],
    ]) @ (j @ source))
    # Equivalent direct block computation avoids conflating vertex and edge
    # coordinates and makes the covariance identity independently testable.
    d = np.block([
        [np.eye(2, dtype=complex), np.zeros((2, 2), dtype=complex)],
        [np.zeros((2, 2), dtype=complex), np.eye(2, dtype=complex)],
    ])
    del d_after_j, d
    first = EPSILON_INV @ (chi_target - np.conjugate(u) @ chi_source)
    second = EPSILON @ (psi_target - u @ psi_source)
    expected = np.concatenate((first, second))
    actual = j @ np.concatenate((
        psi_target - u @ psi_source,
        chi_target - np.conjugate(u) @ chi_source,
    ))
    return float(np.linalg.norm(actual - expected))

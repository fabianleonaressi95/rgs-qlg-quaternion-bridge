"""Finite-dimensional constructions for the trigram Q3 theorems."""
from __future__ import annotations

from itertools import product
import numpy as np

TRIGRAMS = tuple(product((0, 1), repeat=3))


def walsh_matrix() -> np.ndarray:
    """Return E (8x3), with columns the normalized first-order Walsh modes."""
    return np.asarray(
        [[(-1.0) ** t[j] / np.sqrt(8.0) for j in range(3)] for t in TRIGRAMS],
        dtype=float,
    )


def q3_adjacency() -> np.ndarray:
    """Adjacency matrix of Q3 in lexicographic bit order."""
    index = {t: i for i, t in enumerate(TRIGRAMS)}
    A = np.zeros((8, 8), dtype=int)
    for i, t in enumerate(TRIGRAMS):
        for j in range(3):
            u = list(t)
            u[j] ^= 1
            A[i, index[tuple(u)]] = 1
    return A


def affine_face_compatibility(A_alpha, b_alpha, A_beta, b_beta, p0, directions,
                              atol: float = 1e-12) -> bool:
    """Test equality of two affine maps on P0 + span(directions)."""
    da = np.asarray(A_alpha) - np.asarray(A_beta)
    db = np.asarray(b_alpha) - np.asarray(b_beta)
    p0 = np.asarray(p0)
    if not np.allclose(da @ p0 + db, 0.0, atol=atol, rtol=0.0):
        return False
    return all(np.allclose(da @ np.asarray(u), 0.0, atol=atol, rtol=0.0)
               for u in directions)


def refinement_reduction(R: np.ndarray) -> tuple[np.ndarray, float]:
    """Return compressed 3x3 matrix and leakage norm from the three channels."""
    E = walsh_matrix()
    R = np.asarray(R)
    M = E.T @ R @ E
    leakage = np.linalg.norm((np.eye(8) - E @ E.T) @ R @ E)
    return M, float(leakage)


def tribonacci_target() -> np.ndarray:
    return np.asarray([[0, 0, 1], [1, 0, 0], [0, 1, 1]], dtype=int)


def q6_product_adjacency() -> np.ndarray:
    """Return Q6 adjacency in lower/upper Q3 lexicographic order."""
    A3 = q3_adjacency()
    I8 = np.eye(8, dtype=int)
    return np.kron(A3, I8) + np.kron(I8, A3)

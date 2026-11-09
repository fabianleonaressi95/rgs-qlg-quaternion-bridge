"""Exact obstruction checks for a trigram Q3 refinement transfer."""
from __future__ import annotations
from itertools import product, permutations
import numpy as np

TRIGRAMS = tuple(product((0, 1), repeat=3))


def walsh_first_order() -> np.ndarray:
    return np.asarray(
        [[(-1.0) ** t[j] / np.sqrt(8.0) for j in range(3)] for t in TRIGRAMS],
        dtype=float,
    )


def q3_adjacency() -> np.ndarray:
    idx = {t: i for i, t in enumerate(TRIGRAMS)}
    A = np.zeros((8, 8), dtype=float)
    for i, t in enumerate(TRIGRAMS):
        for j in range(3):
            s = list(t)
            s[j] ^= 1
            A[i, idx[tuple(s)]] = 1.0
    return A


def translation_permutation(g: tuple[int, int, int]) -> np.ndarray:
    idx = {t: i for i, t in enumerate(TRIGRAMS)}
    T = np.zeros((8, 8), dtype=float)
    for i, t in enumerate(TRIGRAMS):
        T[idx[tuple(t[j] ^ g[j] for j in range(3))], i] = 1.0
    return T


def coordinate_permutation_matrix(p: tuple[int, int, int]) -> np.ndarray:
    P = np.zeros((3, 3), dtype=float)
    for old, new in enumerate(p):
        P[new, old] = 1.0
    return P


def tribonacci_target() -> np.ndarray:
    return np.asarray([[0., 0., 1.], [1., 0., 0.], [0., 1., 1.]])


def q6_adjacency() -> np.ndarray:
    A3 = q3_adjacency()
    I8 = np.eye(8)
    return np.kron(A3, I8) + np.kron(I8, A3)


def q6_lower_channel_compression() -> np.ndarray:
    E = walsh_first_order()
    u0 = np.ones(8) / np.sqrt(8.0)
    F = np.kron(E, u0[:, None])
    return F.T @ q6_adjacency() @ F


def symmetric_channel_form(M: np.ndarray, atol: float = 1e-12) -> bool:
    M = np.asarray(M)
    alpha = float(np.trace(M) / 3.0)
    off = M[~np.eye(3, dtype=bool)]
    beta = float(np.mean(off))
    target = alpha*np.eye(3) + beta*(np.ones((3, 3))-np.eye(3))
    return bool(np.allclose(M, target, atol=atol, rtol=0.0))


def algebraic_realization(M: np.ndarray, lam: float = 0.0) -> np.ndarray:
    """Diagnostic realization only; deliberately not a geometric refinement."""
    E = walsh_first_order()
    return E @ np.asarray(M) @ E.T + lam * (np.eye(8) - E @ E.T)


def audit() -> dict:
    E, A3, Mpsi = walsh_first_order(), q3_adjacency(), tribonacci_target()
    A6_compressed = q6_lower_channel_compression()
    Rdiag = algebraic_realization(Mpsi, lam=0.25)
    recovered = E.T @ Rdiag @ E
    leakage = np.linalg.norm((np.eye(8)-E@E.T) @ Rdiag @ E)
    Rconv = 2*np.eye(8) + A3
    translation_residual = max(
        np.linalg.norm(Rconv @ translation_permutation(g)
                       - translation_permutation(g) @ Rconv)
        for g in TRIGRAMS
    )
    Mconv = E.T @ Rconv @ E
    return {
        "q3_walsh_orthogonality_residual": float(np.linalg.norm(E.T@E-np.eye(3))),
        "q3_adjacency_channel_residual": float(np.linalg.norm(A3@E-E)),
        "translation_commutator_max": float(translation_residual),
        "translation_invariant_compression": Mconv.tolist(),
        "translation_compression_offdiagonal_max": float(np.max(np.abs(Mconv-np.diag(np.diag(Mconv))))),
        "q6_lower_channel_compression": A6_compressed.tolist(),
        "q6_compression_equals_4I_residual": float(np.linalg.norm(A6_compressed-4*np.eye(3))),
        "tribonacci_is_translation_diagonal": bool(np.allclose(Mpsi, np.diag(np.diag(Mpsi)))),
        "tribonacci_has_full_S3_commutant_form": symmetric_channel_form(Mpsi),
        "diagnostic_realization_compression_residual": float(np.linalg.norm(recovered-Mpsi)),
        "diagnostic_realization_leakage": float(leakage),
        "warning": "The diagnostic algebraic realization is circular, not a geometric derivation.",
    }

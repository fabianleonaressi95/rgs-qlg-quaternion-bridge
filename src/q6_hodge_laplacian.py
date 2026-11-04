"""Weighted combinatorial Hodge Laplacians on the canonical Q6 square complex.

Positive diagonal cochain inner products define the adjoints of the coboundary
maps. Returned matrices are symmetric representatives in orthonormalized
cochain coordinates. Unit weights recover the standard combinatorial
Hodge Laplacians. No continuum metric or dual-cell construction is assumed.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from src.q6_complex import Q6Complex, build_q6_complex


@dataclass(frozen=True)
class Q6HodgeLaplacians:
    weights_0: np.ndarray
    weights_1: np.ndarray
    weights_2: np.ndarray
    laplacian_0: np.ndarray
    laplacian_1: np.ndarray
    laplacian_2: np.ndarray


def _positive_weights(values: np.ndarray | None, size: int, name: str) -> np.ndarray:
    if values is None:
        return np.ones(size, dtype=float)
    weights = np.asarray(values, dtype=float)
    if weights.shape != (size,):
        raise ValueError(f"{name} must have shape ({size},)")
    if not np.all(np.isfinite(weights)) or np.any(weights <= 0.0):
        raise ValueError(f"{name} must contain finite, strictly positive values")
    return weights.copy()


def build_q6_hodge_laplacians(
    weights_0: np.ndarray | None = None,
    weights_1: np.ndarray | None = None,
    weights_2: np.ndarray | None = None,
    complex_: Q6Complex | None = None,
) -> Q6HodgeLaplacians:
    """Build weighted Hodge Laplacians in degrees 0, 1 and 2.

    With d0=B1.T and d1=B2.T, Hk=diag(weights_k) defines the cochain
    inner products and d_k^*=H_k^{-1} d_k.T H_{k+1}. The symmetric
    representatives are H_k^(1/2) Delta_k H_k^(-1/2).
    """
    q6 = build_q6_complex() if complex_ is None else complex_
    n0, n1, n2 = len(q6.vertices), len(q6.edges), len(q6.faces)
    h0 = _positive_weights(weights_0, n0, "weights_0")
    h1 = _positive_weights(weights_1, n1, "weights_1")
    h2 = _positive_weights(weights_2, n2, "weights_2")

    b1 = np.asarray(q6.boundary_1, dtype=np.int64)
    b2 = np.asarray(q6.boundary_2, dtype=np.int64)
    if b1.shape != (n0, n1) or b2.shape != (n1, n2):
        raise ValueError("incidence matrix shapes are inconsistent with the cell counts")
    if np.count_nonzero(b1 @ b2):
        raise ValueError("invalid chain complex: B1 @ B2 must vanish")

    d0, d1 = b1.astype(float).T, b2.astype(float).T
    ih0, ih1, ih2 = 1/np.sqrt(h0), 1/np.sqrt(h1), 1/np.sqrt(h2)
    sh1, sh2 = np.sqrt(h1), np.sqrt(h2)

    # Delta_0 = H0^(-1/2) d0.T H1 d0 H0^(-1/2).
    gram0 = d0.T @ (h1[:, None] * d0)
    lap0 = ih0[:, None] * gram0 * ih0[None, :]

    # Delta_1 = H1^(1/2) d0 H0^(-1) d0.T H1^(1/2)
    #         + H1^(-1/2) d1.T H2 d1 H1^(-1/2).
    lower = d0 @ (ih0**2[:, None] * d0.T)
    upper = d1.T @ (h2[:, None] * d1)
    lap1 = sh1[:, None] * lower * sh1[None, :] + ih1[:, None] * upper * ih1[None, :]

    # Delta_2 = H2^(1/2) d1 H1^(-1) d1.T H2^(1/2).
    gram2 = d1 @ (ih1**2[:, None] * d1.T)
    lap2 = sh2[:, None] * gram2 * sh2[None, :]

    return Q6HodgeLaplacians(
        weights_0=h0, weights_1=h1, weights_2=h2,
        laplacian_0=(lap0 + lap0.T)/2,
        laplacian_1=(lap1 + lap1.T)/2,
        laplacian_2=(lap2 + lap2.T)/2,
    )


def audit_q6_hodge_laplacians(
    hodge: Q6HodgeLaplacians | None = None,
    complex_: Q6Complex | None = None,
    atol: float = 1e-9,
) -> dict[str, object]:
    """Audit symmetry, positivity and harmonic dimensions against homology."""
    q6 = build_q6_complex() if complex_ is None else complex_
    model = build_q6_hodge_laplacians(complex_=q6) if hodge is None else hodge
    matrices = (model.laplacian_0, model.laplacian_1, model.laplacian_2)
    eigenvalues = tuple(np.linalg.eigvalsh(matrix) for matrix in matrices)
    sizes = (len(q6.vertices), len(q6.edges), len(q6.faces))
    ranks = (
        int(np.linalg.matrix_rank(q6.boundary_1.astype(float))),
        int(np.linalg.matrix_rank(q6.boundary_2.astype(float))),
    )
    expected_betti = (
        sizes[0] - ranks[0],
        sizes[1] - ranks[0] - ranks[1],
        sizes[2] - ranks[1],
    )
    harmonic_dimensions = tuple(
        int(np.count_nonzero(np.abs(values) <= atol)) for values in eigenvalues
    )
    symmetry_residuals = tuple(
        float(np.linalg.norm(matrix - matrix.T)) for matrix in matrices
    )
    minimum_eigenvalues = tuple(float(values[0]) for values in eigenvalues)
    chain_residual = q6.boundary_1.astype(np.int64) @ q6.boundary_2.astype(np.int64)
    pass_ = (
        np.all(np.asarray(symmetry_residuals) <= atol)
        and np.all(np.asarray(minimum_eigenvalues) >= -atol)
        and harmonic_dimensions == expected_betti
        and np.count_nonzero(chain_residual) == 0
    )
    return {
        "cell_counts": {"C0": sizes[0], "C1": sizes[1], "C2": sizes[2]},
        "boundary_ranks": {"rank_B1": ranks[0], "rank_B2": ranks[1]},
        "expected_betti_numbers": expected_betti,
        "observed_harmonic_dimensions": harmonic_dimensions,
        "minimum_eigenvalues": minimum_eigenvalues,
        "symmetry_residuals": symmetry_residuals,
        "boundary_squared_exact_zero": bool(np.count_nonzero(chain_residual) == 0),
        "pass": bool(pass_),
        "scope": (
            "finite weighted combinatorial Hodge theory on the Q6 square complex; "
            "no metric-derived dual cells, continuum limit, or physical field dynamics"
        ),
    }

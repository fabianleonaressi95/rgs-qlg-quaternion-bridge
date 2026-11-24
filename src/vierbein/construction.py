"""Vierbein construction from RGS covariance matrices.

Historical reference values (to be verified):
    ℓ² = 7.035514834154695
    c₀ = 43.96495130838031

Construction:
    D_ij = X_j - X_i
    W_ij = exp(-‖D_ij‖² / ℓ²)
    C_i = Σ_{j≠i} W_ij D_ij D_ij^T
    E_i = √(C_i / c₀)  (principal symmetric positive square root)
"""

import numpy as np
from numpy.typing import NDArray
from typing import Tuple, Dict
import warnings


FloatArray = NDArray[np.float64]


def compute_distances(
    X: FloatArray,
) -> FloatArray:
    """Compute distance matrix D_ij = X_j - X_i.

    Args:
        X: (N, 3) array of vertex positions

    Returns:
        D: (N, N, 3) array where D[i,j] = X[j] - X[i]
    """
    N = X.shape[0]
    D = np.zeros((N, N, 3), dtype=np.float64)

    for i in range(N):
        for j in range(N):
            D[i, j] = X[j] - X[i]

    return D


def compute_weights(
    D: FloatArray,
    ell_squared: float,
) -> FloatArray:
    """Compute Gaussian weights W_ij = exp(-‖D_ij‖² / ℓ²).

    Args:
        D: (N, N, 3) distance array
        ell_squared: ℓ² parameter

    Returns:
        W: (N, N) weight matrix
    """
    N = D.shape[0]
    W = np.zeros((N, N), dtype=np.float64)

    for i in range(N):
        for j in range(N):
            if i != j:
                dist_sq = np.sum(D[i, j] ** 2)
                W[i, j] = np.exp(-dist_sq / ell_squared)

    return W


def compute_covariance_matrices(
    D: FloatArray,
    W: FloatArray,
) -> FloatArray:
    """Compute covariance matrices C_i = Σ_{j≠i} W_ij D_ij D_ij^T.

    Args:
        D: (N, N, 3) distance array
        W: (N, N) weight matrix

    Returns:
        C: (N, 3, 3) array of covariance matrices
    """
    N = D.shape[0]
    C = np.zeros((N, 3, 3), dtype=np.float64)

    for i in range(N):
        for j in range(N):
            if i != j:
                w_ij = W[i, j]
                d_ij = D[i, j].reshape(3, 1)
                C[i] += w_ij * (d_ij @ d_ij.T)

    return C


def compute_vierbein(
    C: FloatArray,
    c0: float,
    tol: float = 1e-12,
) -> Tuple[FloatArray, Dict[str, any]]:
    """Compute vierbein E_i = √(C_i / c₀).

    Uses principal symmetric positive square root via SVD.

    Args:
        C: (N, 3, 3) array of covariance matrices
        c0: Normalization constant
        tol: Tolerance for positive definiteness check

    Returns:
        E: (N, 3, 3) array of vierbein matrices
        diagnostics: Dictionary with diagnostic information
    """
    N = C.shape[0]
    E = np.zeros((N, 3, 3), dtype=np.float64)

    min_eigenvalue = np.inf
    max_condition_number = 0.0
    positive_definite_count = 0

    for i in range(N):
        # Compute symmetric positive square root via SVD
        # C / c0 = U Σ V^T, with U = V for symmetric matrix
        # √(C/c0) = U √Σ U^T
        Ci_scaled = C[i] / c0

        # SVD
        U, s, Vt = np.linalg.svd(Ci_scaled, hermitian=True)

        # Check positive definiteness
        if np.all(s > tol):
            positive_definite_count += 1
            min_eigenvalue = min(min_eigenvalue, np.min(s))

            # Compute square root
            sqrt_s = np.sqrt(s)
            E[i] = U @ np.diag(sqrt_s) @ Vt

            # Condition number
            cond = np.max(s) / np.min(s)
            max_condition_number = max(max_condition_number, cond)
        else:
            warnings.warn(f"C[{i}] / c0 has non-positive eigenvalues: {s}")
            # Still compute using absolute values (may not be physical)
            sqrt_s = np.sqrt(np.maximum(s, 0))
            E[i] = U @ np.diag(sqrt_s) @ Vt

    diagnostics = {
        "positive_definite_count": positive_definite_count,
        "total_matrices": N,
        "min_eigenvalue": float(min_eigenvalue),
        "max_condition_number": float(max_condition_number),
        "c0_value": float(c0),
    }

    return E, diagnostics


def build_vierbein(
    X: FloatArray,
    ell_squared: float = 7.035514834154695,
    c0: float = 43.96495130838031,
    tol: float = 1e-12,
) -> Tuple[FloatArray, Dict[str, any]]:
    """Complete vierbein construction pipeline.

    Args:
        X: (N, 3) array of vertex positions
        ell_squared: ℓ² parameter (historical value)
        c0: Normalization constant (historical value)
        tol: Tolerance for positive definiteness

    Returns:
        E: (N, 3, 3) array of vierbein matrices
        metadata: Dictionary with all intermediate results and diagnostics
    """
    # Step 1: Compute distances
    D = compute_distances(X)

    # Step 2: Compute weights
    W = compute_weights(D, ell_squared)

    # Step 3: Compute covariance matrices
    C = compute_covariance_matrices(D, W)

    # Step 4: Compute vierbein
    E, vierbein_diagnostics = compute_vierbein(C, c0, tol)

    metadata = {
        "ell_squared": float(ell_squared),
        "c0": float(c0),
        "num_vertices": int(X.shape[0]),
        "distance_matrix_shape": tuple(D.shape),
        "weight_matrix_shape": tuple(W.shape),
        "covariance_shape": tuple(C.shape),
        "vierbein_shape": tuple(E.shape),
        **vierbein_diagnostics,
    }

    return E, metadata


def reconstruct_metric(
    E: FloatArray,
    metric_signature: Tuple[int, int, int, int] = (1, -1, -1, -1),
) -> FloatArray:
    """Reconstruct metric from vierbein: g_μν = η_AB e^A_μ e^B_ν.

    For spatial vierbein (3D), we use spatial metric δ_ij.

    Args:
        E: (N, 3, 3) array of vierbein matrices
        metric_signature: Signature for 4D extension (not used for 3D)

    Returns:
        g: (N, 3, 3) array of metric tensors
    """
    N = E.shape[0]
    g = np.zeros((N, 3, 3), dtype=np.float64)

    # For 3D spatial case, η_AB = δ_AB
    for i in range(N):
        g[i] = E[i].T @ E[i]

    return g


def verify_metric_identity(
    C: FloatArray,
    c0: float,
    g: FloatArray,
    tol: float = 1e-14,
) -> Dict[str, float]:
    """Verify the identity g = C / c0.

    Args:
        C: (N, 3, 3) covariance matrices
        c0: Normalization constant
        g: (N, 3, 3) reconstructed metrics
        tol: Tolerance for identity check

    Returns:
        diagnostics: Dictionary with max residual and pass/fail status
    """
    N = C.shape[0]
    max_residual = 0.0

    for i in range(N):
        expected = C[i] / c0
        residual = np.linalg.norm(g[i] - expected, ord="fro")
        max_residual = max(max_residual, residual)

    return {
        "max_residual": float(max_residual),
        "tolerance": float(tol),
        "identity_holds": max_residual <= tol,
    }

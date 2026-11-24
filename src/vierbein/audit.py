"""Audit functions for vierbein construction and Killing symmetries."""

import numpy as np
from numpy.typing import NDArray
from typing import Dict, Tuple, List
from .construction import build_vierbein, reconstruct_metric, verify_metric_identity
from .killing import KillingOperator, find_killing_vectors, compute_lie_bracket


FloatArray = NDArray[np.float64]


def audit_vierbein(
    X: FloatArray,
    ell_squared: float = 7.035514834154695,
    c0: float = 43.96495130838031,
    tol: float = 1e-12,
) -> Dict[str, any]:
    """Complete audit of vierbein construction.

    Args:
        X: (N, 3) vertex positions
        ell_squared: ℓ² parameter
        c0: Normalization constant
        tol: Tolerance for checks

    Returns:
        audit_results: Dictionary with all audit results
    """
    # Build vierbein
    E, metadata = build_vierbein(X, ell_squared, c0, tol)

    # Reconstruct metric
    g = reconstruct_metric(E)

    # Verify metric identity
    # Need C matrices - recompute for verification
    from .construction import (
        compute_distances,
        compute_weights,
        compute_covariance_matrices,
    )

    D = compute_distances(X)
    W = compute_weights(D, ell_squared)
    C = compute_covariance_matrices(D, W)

    identity_check = verify_metric_identity(C, c0, g, tol)

    # Additional vierbein checks
    N = X.shape[0]
    norm_checks = []
    det_checks = []
    condition_numbers = []

    for i in range(N):
        # Check norm (should be unit for proper vierbein)
        E_norm = np.linalg.norm(E[i], ord="fro")
        norm_checks.append(E_norm)

        # Check determinant
        det = np.linalg.det(E[i])
        det_checks.append(det)

        # Condition number
        if np.linalg.det(E[i]) != 0:
            cond = np.linalg.cond(E[i])
            condition_numbers.append(cond)

    audit_results = {
        "vierbein_construction": {
            "positive_definite_count": metadata["positive_definite_count"],
            "total_matrices": metadata["total_matrices"],
            "min_eigenvalue": metadata["min_eigenvalue"],
            "max_condition_number": metadata["max_condition_number"],
        },
        "metric_identity": identity_check,
        "vierbein_norms": {
            "mean": float(np.mean(norm_checks)),
            "std": float(np.std(norm_checks)),
            "min": float(np.min(norm_checks)),
            "max": float(np.max(norm_checks)),
        },
        "vierbein_determinants": {
            "mean": float(np.mean(det_checks)),
            "std": float(np.std(det_checks)),
            "min": float(np.min(det_checks)),
            "max": float(np.max(det_checks)),
        },
        "condition_numbers": {
            "mean": float(np.mean(condition_numbers)) if condition_numbers else None,
            "max": float(np.max(condition_numbers)) if condition_numbers else None,
        },
        "num_vertices": int(N),
        "ell_squared": float(ell_squared),
        "c0": float(c0),
    }

    return audit_results


def audit_metric(
    g: FloatArray,
    X: FloatArray,
    tol: float = 1e-10,
) -> Dict[str, any]:
    """Audit reconstructed metric properties.

    Args:
        g: (N, 3, 3) metric tensors
        X: (N, 3) vertex positions
        tol: Tolerance for checks

    Returns:
        audit_results: Dictionary with metric audit results
    """
    N = g.shape[0]

    # Check symmetry
    symmetry_errors = []
    for i in range(N):
        symmetry_error = np.linalg.norm(g[i] - g[i].T, ord="fro")
        symmetry_errors.append(symmetry_error)

    # Check positive definiteness
    min_eigenvalues = []
    positive_definite_count = 0

    for i in range(N):
        eigvals = np.linalg.eigvalsh(g[i])
        min_eigenvalues.append(np.min(eigvals))
        if np.all(eigvals > tol):
            positive_definite_count += 1

    # Check signature (should be + + + for spatial metric)
    signature_checks = []
    for i in range(N):
        eigvals = np.linalg.eigvalsh(g[i])
        positive = np.sum(eigvals > tol)
        negative = np.sum(eigvals < -tol)
        signature_checks.append((positive, negative))

    audit_results = {
        "symmetry": {
            "max_error": float(np.max(symmetry_errors)),
            "mean_error": float(np.mean(symmetry_errors)),
            "all_symmetric": float(np.max(symmetry_errors)) < tol,
        },
        "positive_definiteness": {
            "count": positive_definite_count,
            "total": N,
            "min_eigenvalue": float(np.min(min_eigenvalues)),
            "all_positive": positive_definite_count == N,
        },
        "signature": {
            "positive_counts": [s[0] for s in signature_checks],
            "negative_counts": [s[1] for s in signature_checks],
            "expected_positive": 3,
            "expected_negative": 0,
        },
    }

    return audit_results


def audit_killing_symmetries(
    X: FloatArray,
    g: FloatArray,
    adjacency: FloatArray = None,
    tol: float = 1e-8,
) -> Dict[str, any]:
    """Audit Killing vector symmetries of the metric.

    Args:
        X: (N, 3) vertex positions
        g: (N, 3, 3) metric tensors
        adjacency: (N, N) adjacency matrix
        tol: Tolerance for Killing equation

    Returns:
        audit_results: Dictionary with Killing symmetry audit results
    """
    # Find Killing vectors
    killing_fields, metadata = find_killing_vectors(X, g, adjacency, tol=tol)

    N = X.shape[0]
    num_killing = len(killing_fields)

    # If we found Killing vectors, compute their algebra
    lie_brackets = []
    if num_killing > 0:
        for i in range(num_killing):
            for j in range(i + 1, num_killing):
                bracket = compute_lie_bracket(
                    killing_fields[i], killing_fields[j], X
                )
                lie_brackets.append({
                    "indices": (i, j),
                    "bracket_norm": float(np.linalg.norm(bracket)),
                })

    # Check if Killing vectors form a closed algebra
    algebra_closed = True
    if num_killing > 1:
        for bracket_info in lie_brackets:
            if bracket_info["bracket_norm"] > tol:
                algebra_closed = False
                break

    audit_results = {
        "num_killing_vectors": num_killing,
        "killing_residuals": metadata["killing_residuals"],
        "tolerance": float(tol),
        "algebra_closed": algebra_closed,
        "lie_brackets": lie_brackets if num_killing > 1 else [],
        "expected_dimension": 6,  # For Euclidean 3D space
        "dimension_matches_expected": num_killing == 6,
    }

    return audit_results


def print_audit_report(
    vierbein_audit: Dict[str, any],
    metric_audit: Dict[str, any],
    killing_audit: Dict[str, any],
):
    """Print formatted audit report."""
    print("=" * 70)
    print("RGS-PSU VIERBEIN AND KILLING SYMMETRY AUDIT REPORT")
    print("=" * 70)
    print()

    print("VIERBEIN CONSTRUCTION")
    print("-" * 70)
    vc = vierbein_audit["vierbein_construction"]
    print(f"Positive definite matrices: {vc['positive_definite_count']}/{vc['total_matrices']}")
    print(f"Min eigenvalue: {vc['min_eigenvalue']:.2e}")
    print(f"Max condition number: {vc['max_condition_number']:.2e}")
    print()

    print("METRIC IDENTITY (g = C/c0)")
    print("-" * 70)
    mi = vierbein_audit["metric_identity"]
    print(f"Max residual: {mi['max_residual']:.2e}")
    print(f"Identity holds: {mi['identity_holds']}")
    print()

    print("METRIC PROPERTIES")
    print("-" * 70)
    mp = metric_audit["symmetry"]
    print(f"Max symmetry error: {mp['max_error']:.2e}")
    print(f"All symmetric: {mp['all_symmetric']}")
    pd = metric_audit["positive_definiteness"]
    print(f"Positive definite: {pd['count']}/{pd['total']}")
    print(f"Min eigenvalue: {pd['min_eigenvalue']:.2e}")
    print()

    print("KILLING SYMMETRIES")
    print("-" * 70)
    print(f"Number of Killing vectors found: {killing_audit['num_killing_vectors']}")
    print(f"Expected (Euclidean 3D): {killing_audit['expected_dimension']}")
    print(f"Algebra closed: {killing_audit['algebra_closed']}")
    if killing_audit["killing_residuals"]:
        print(f"Killing residuals: {[f'{r:.2e}' for r in killing_audit['killing_residuals']]}")
    print()

    print("=" * 70)

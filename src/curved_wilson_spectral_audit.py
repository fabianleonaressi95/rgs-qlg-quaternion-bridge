"""Finite curved-graph spectral diagnostics for the Wilson kernel.

This module deliberately reports finite-matrix diagnostics only. It does not
infer a Brillouin zone, a continuum species count, or refinement convergence.
"""
from __future__ import annotations
import numpy as np
from src.discrete_dirac import (
    assemble_euclidean_wilson_dirac,
    assemble_wilson_laplacian,
    chirality_matrix,
)


def curved_graph_wilson_spectral_audit(
    coordinates, vierbein, edges, spin_links=None, wilson_r=1.0,
    physical_cutoff=1e-8, tolerance=1e-9,
):
    """Return finite-graph Wilson-kernel and Euclidean-Dirac diagnostics."""
    cutoff, tol = float(physical_cutoff), float(tolerance)
    if not np.isfinite(cutoff) or cutoff < 0:
        raise ValueError("physical_cutoff must be finite and non-negative")
    if not np.isfinite(tol) or tol <= 0:
        raise ValueError("tolerance must be finite and positive")
    w, wa = assemble_wilson_laplacian(
        coordinates, vierbein, edges, spin_links=spin_links, wilson_r=wilson_r
    )
    de, da = assemble_euclidean_wilson_dirac(
        coordinates, vierbein, edges, spin_links=spin_links, wilson_r=wilson_r
    )
    # W is Hermitian positive semidefinite by its edgewise construction.
    wvals = np.linalg.eigvalsh(w)
    singular = np.linalg.svd(de, compute_uv=False)
    gamma5 = np.kron(np.eye(len(np.asarray(coordinates))), chirality_matrix())
    gauge_covariance_residual = None
    link_unitarity_residual = 0.0
    if spin_links is not None:
        for link in spin_links.values():
            q = np.asarray(link, dtype=complex)
            if q.shape == (2, 2):
                link_unitarity_residual = max(
                    link_unitarity_residual,
                    float(np.linalg.norm(q.conj().T @ q - np.eye(2), ord="fro")),
                    float(abs(np.linalg.det(q) - 1.0)),
                )
            elif q.shape == (4, 4):
                link_unitarity_residual = max(
                    link_unitarity_residual,
                    float(np.linalg.norm(q.conj().T @ q - np.eye(4), ord="fro")),
                )
            else:
                raise ValueError("spin links must be 2x2 SU(2) or 4x4 unitary")
    w_herm = float(np.linalg.norm(w - w.conj().T, ord="fro"))
    de_gamma5 = float(np.linalg.norm(
        de.conj().T - gamma5 @ de @ gamma5, ord="fro"
    ))
    report = {
        "scope": "FINITE_CURVED_GRAPH_ONLY",
        "vertices": int(len(np.asarray(coordinates))),
        "spinor_dimension": 4,
        "matrix_dimension": int(w.shape[0]),
        "wilson_r": float(wilson_r),
        "physical_cutoff": cutoff,
        "wilson_kernel": {
            "minimum_eigenvalue": float(wvals[0]),
            "maximum_eigenvalue": float(wvals[-1]),
            "near_zero_multiplicity": int(np.sum(wvals <= tol)),
            "multiplicity_below_cutoff": int(np.sum(wvals <= cutoff)),
            "hermiticity_residual": w_herm,
            "assembler_minimum_eigenvalue": wa["minimum_eigenvalue"],
        },
        "euclidean_dirac": {
            "minimum_singular_value": float(singular[-1]),
            "maximum_singular_value": float(singular[0]),
            "singular_values_below_cutoff": int(np.sum(singular <= cutoff)),
            "gamma5_hermiticity_residual": de_gamma5,
            "assembler_gamma5_hermiticity_residual": da["gamma5_hermiticity_residual"],
        },
        "link_unitarity_residual": link_unitarity_residual,
        "checks": {
            "wilson_kernel_hermitian": w_herm <= tol,
            "wilson_kernel_positive_semidefinite": float(wvals[0]) >= -tol,
            "euclidean_gamma5_hermiticity": de_gamma5 <= tol,
            "provided_links_unitary": link_unitarity_residual <= tol,
        },
        "limitations": [
            "finite graph spectrum is not a Brillouin-zone species count",
            "a single graph does not establish uniform estimates over refinement",
            "no curved continuum or strong-resolvent convergence is asserted",
        ],
    }
    report["all_finite_checks_passed"] = all(report["checks"].values())
    return report

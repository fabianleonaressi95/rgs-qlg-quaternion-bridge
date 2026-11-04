"""Exact tensor-product refinement of the Q6 coexact five-mode sector.

This module concerns the cubical product complex on [0, 1]^6 with m
subintervals per coordinate and the standard tensor-product cochain inner
product (absolute/Neumann boundary convention). It is not a claim about
Wilson-measure convergence or a Whitney mass-matrix discretization.
"""
from __future__ import annotations

import math
import numpy as np


def cubical_cell_counts(m: int) -> dict[str, int]:
    """Numbers of k-cells for k=0,1,2 in the m^6 cubical subdivision."""
    _validate_m(m)
    return {
        "vertices": (m + 1) ** 6,
        "edges": 6 * m * (m + 1) ** 5,
        "square_faces": 15 * m**2 * (m + 1) ** 4,
    }


def interval_hodge_eigenvalue(j: int, m: int) -> float:
    """Nonzero interval Hodge eigenvalue; j=0 denotes the scalar zero mode."""
    _validate_m(m)
    if not isinstance(j, int) or not 0 <= j <= m:
        raise ValueError("j must be an integer in [0, m]")
    if j == 0:
        return 0.0
    return 2.0 - 2.0 * math.cos(j * math.pi / (m + 1))


def helmert_coexact_basis() -> np.ndarray:
    """Orthonormal basis of {x in R^6 : sum(x)=0}; shape (6,5)."""
    out = np.zeros((6, 5), dtype=float)
    for col in range(5):
        k = col + 1
        scale = math.sqrt(k * (k + 1))
        out[:k, col] = 1.0 / scale
        out[k, col] = -k / scale
    return out


def q6_fundamental_coexact_sector(m: int = 1) -> dict:
    """Exact fundamental joint spectral block and its coexact part.

    The joint index (1,1,1,1,1,1) has six 1-form component vectors and
    eigenvalue 6*mu_1. The gradient of the corresponding scalar product
    eigenfunction is the diagonal vector (1,...,1)/sqrt(6). Its orthogonal
    complement is the rank-five coexact sector.
    """
    _validate_m(m)
    mu = interval_hodge_eigenvalue(1, m)
    basis = helmert_coexact_basis()
    exact_unit = np.ones(6, dtype=float) / math.sqrt(6.0)
    checks = {
        "basis_shape_6_by_5": basis.shape == (6, 5),
        "basis_orthonormal": float(np.linalg.norm(basis.T @ basis - np.eye(5))) < 1e-13,
        "orthogonal_to_exact_gradient": float(np.linalg.norm(basis.T @ exact_unit)) < 1e-13,
        "rank_five": int(np.linalg.matrix_rank(basis)) == 5,
    }
    return {
        "subdivision": m,
        "mesh_spacing": 1.0 / m,
        "cell_counts": cubical_cell_counts(m),
        "joint_index": (1, 1, 1, 1, 1, 1),
        "interval_eigenvalue": mu,
        "full_one_form_eigenvalue": 6.0 * mu,
        "exact_rank": 1,
        "coexact_rank": 5,
        "rescaled_coexact_eigenvalue": (m * m) * 6.0 * mu,
        "continuum_limit": 6.0 * math.pi**2,
        "coexact_coordinate_basis": basis,
        "checks": checks,
        "all_checks_passed": all(checks.values()),
        "scope": (
            "Exact tensor-product combinatorial Hodge result with the "
            "standard product inner product and absolute/Neumann boundary "
            "convention. It identifies a rank-five coexact joint block; "
            "it does not establish isolation from all other blocks having "
            "the same total eigenvalue."
        ),
    }


def _validate_m(m: int) -> None:
    if not isinstance(m, int) or isinstance(m, bool) or m < 1:
        raise ValueError("m must be a positive integer")

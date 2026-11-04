"""Finite cochain complex and combinatorial Hodge Laplacians.

Boundary convention: B1:C1->C0 and B2:C2->C1. Cochain differentials
are d0=B1.T and d1=B2.T, so d1@d0=0 iff B1@B2=0.
The returned Laplacians use unit cell weights; weighted Hodge stars are
an explicit future extension, not silently inferred from geometry.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class DiscreteDifferentialOperators:
    d0: np.ndarray
    d1: np.ndarray
    laplacian_0: np.ndarray
    laplacian_1: np.ndarray
    laplacian_2: np.ndarray


def build_discrete_differential_operators(
    boundary_1: np.ndarray, boundary_2: np.ndarray
) -> DiscreteDifferentialOperators:
    b1 = np.asarray(boundary_1, dtype=float)
    b2 = np.asarray(boundary_2, dtype=float)
    if b1.ndim != 2 or b2.ndim != 2 or b1.shape[1] != b2.shape[0]:
        raise ValueError("Boundary matrices must have compatible C2->C1->C0 shapes.")
    if not np.all(np.isfinite(b1)) or not np.all(np.isfinite(b2)):
        raise ValueError("Boundary matrices must be finite.")
    chain = b1 @ b2
    if np.max(np.abs(chain), initial=0.0) > 0.0:
        raise ValueError("Boundary maps do not form a chain complex: B1 @ B2 != 0.")
    d0, d1 = b1.T, b2.T
    l0 = d0.T @ d0
    l1 = d0 @ d0.T + d1.T @ d1
    l2 = d1 @ d1.T
    return DiscreteDifferentialOperators(d0, d1, l0, l1, l2)


def audit_discrete_differential_operators(
    boundary_1: np.ndarray, boundary_2: np.ndarray, atol: float = 1e-12
) -> dict[str, object]:
    ops = build_discrete_differential_operators(boundary_1, boundary_2)
    nilpotency = float(np.linalg.norm(ops.d1 @ ops.d0, ord="fro"))
    symmetry = max(
        float(np.linalg.norm(x - x.T, ord="fro"))
        for x in (ops.laplacian_0, ops.laplacian_1, ops.laplacian_2)
    )
    min_eigenvalue = min(
        float(np.linalg.eigvalsh(x).min())
        for x in (ops.laplacian_0, ops.laplacian_1, ops.laplacian_2)
    )
    return {
        "audit": "FINITE_DISCRETE_DE_RHAM_HODGE",
        "cochain_dimensions": [
            int(ops.laplacian_0.shape[0]),
            int(ops.laplacian_1.shape[0]),
            int(ops.laplacian_2.shape[0]),
        ],
        "d1_d0_nilpotency_fro": nilpotency,
        "laplacian_symmetry_max_fro": symmetry,
        "laplacian_min_eigenvalue": min_eigenvalue,
        "passed": bool(nilpotency <= atol and symmetry <= atol and min_eigenvalue >= -atol),
        "scope": "finite unit-weight combinatorial complex; not a continuum or Mosco theorem",
    }

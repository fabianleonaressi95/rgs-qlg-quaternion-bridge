"""Uniform stability diagnostics for the flat periodic P1 refinement family.

The analytic coercivity and compactness statements are documented in
docs/FLAT_TORUS_G5_STABILITY.md. Numerical eigenvalues are diagnostics,
not substitutes for the Rellich/Poincare arguments.
"""
from __future__ import annotations

from typing import Any

import numpy as np

from src.flat_torus_refinement import build_flat_torus_mesh
from src.flat_torus_embedding import assemble_consistent_mass0


def audit_g5_level(
    cells_per_axis: int, dimension: int = 2, side_length: float = 1.0
) -> dict[str, Any]:
    """Measure shape regularity and the first generalized P1 eigenvalue."""
    mesh = build_flat_torus_mesh(cells_per_axis, dimension, side_length)
    mass = assemble_consistent_mass0(mesh)
    stiffness = (mesh.stiffness0 + mesh.stiffness0.T) / 2.0

    # M = L L^T. Transform K u = lambda M u to the symmetric
    # standard problem L^{-1} K L^{-T} y = lambda y.
    chol = np.linalg.cholesky(mass)
    inverse_chol = np.linalg.inv(chol)
    standard = inverse_chol @ stiffness @ inverse_chol.T
    eigenvalues = np.linalg.eigvalsh((standard + standard.T) / 2.0)
    positive = eigenvalues[eigenvalues > 1e-9]
    if positive.size == 0:
        raise RuntimeError("no positive generalized stiffness eigenvalue found")

    h = side_length / cells_per_axis
    continuum_poincare = (2.0 * np.pi / side_length) ** 2
    lambda1 = float(positive[0])
    return {
        "cells_per_axis": int(cells_per_axis),
        "dimension": int(dimension),
        "spacing": float(h),
        "max_simplex_condition": float(mesh.max_simplex_condition),
        "max_edge_over_h": float(mesh.max_edge_length / h),
        "constant_mode_stiffness_residual": float(
            np.linalg.norm(mesh.stiffness0 @ np.ones(mesh.vertices.shape[0]))
        ),
        "first_positive_generalized_eigenvalue": lambda1,
        "continuum_poincare_lower_bound": float(continuum_poincare),
        "poincare_bound_holds": bool(
            lambda1 >= continuum_poincare * (1.0 - 1e-8)
        ),
        "mass_positive_definite": bool(np.linalg.eigvalsh(mass)[0] > 0.0),
        "h1_norm_identity": "||J_h u||_L2^2 + q_h(u)",
    }


def flat_torus_g5_audit(
    levels: tuple[int, ...] = (4, 8, 16),
    dimension: int = 2,
    side_length: float = 1.0,
) -> dict[str, Any]:
    """Audit G5 on at least three regular flat-torus refinement levels."""
    ns = tuple(int(n) for n in levels)
    if len(ns) < 3 or any(b <= a for a, b in zip(ns, ns[1:])):
        raise ValueError("levels must contain at least three strictly increasing integers")
    rows = [audit_g5_level(n, dimension, side_length) for n in ns]
    conditions = np.asarray([row["max_simplex_condition"] for row in rows])
    diameters = np.asarray([row["max_edge_over_h"] for row in rows])
    shape_spread = float(np.max(conditions) - np.min(conditions))
    diameter_spread = float(np.max(diameters) - np.min(diameters))
    return {
        "audit": "FLAT_TORUS_G5_UNIFORM_STABILITY",
        "geometry": "flat periodic torus with exact Euclidean P1 metric",
        "levels": list(ns),
        "dimension": int(dimension),
        "side_length": float(side_length),
        "uniform_shape_regular_diagnostic": bool(shape_spread < 1e-12),
        "shape_condition_spread": shape_spread,
        "uniform_scaled_diameter_diagnostic": bool(diameter_spread < 1e-12),
        "scaled_diameter_spread": diameter_spread,
        "poincare_bound_all_levels": bool(all(row["poincare_bound_holds"] for row in rows)),
        "mass_positive_all_levels": bool(all(row["mass_positive_definite"] for row in rows)),
        "refinements": rows,
        "analytic_compactness_basis": (
            "Rellich compactness: J_h u_h are bounded in H1(T_L^d) "
            "when their exact L2 norms and P1 energies are bounded."
        ),
        "scope_warning": (
            "G5 is established only for this flat periodic P1 family; "
            "no uniform estimate for a curved RGS-derived family is claimed."
        ),
    }

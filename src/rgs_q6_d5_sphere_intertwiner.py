"""A5 intertwiner from the Q6 five-mode cotangent sector to R^3.

The five-dimensional coexact Hodge eigenspace is an irreducible A5 module V5.
Its exterior square contains a 3-dimensional irreducible summand. This module
constructs the unique (up to sign and scale) intertwiner from Lambda^2(V5)
to the actual geometric 3D rotation representation. Applying it to q wedge p
defines an A5-equivariant, SO(2)-phase-invariant relative-direction vector.

This is a candidate nonlinear bridge on the D5 tangent/cotangent sector. It
does not yet identify the finite Q6 vertex 5-design with the pushforward of
the continuous phase-space measure.
"""
from __future__ import annotations

import numpy as np

from src.rgs_q6_d5_phase_space import build_q6_d5_phase_space_sector
from src.rgs_icosahedral_q6_intertwiner import icosahedral_signed_axis_representation


def _wedge2_matrix(rho: np.ndarray) -> np.ndarray:
    """Matrix of Lambda^2(rho) in the lexicographic basis e_i wedge e_j."""
    r = np.asarray(rho, dtype=float)
    if r.shape != (5, 5):
        raise ValueError("rho must have shape (5,5)")
    pairs = [(i, j) for i in range(5) for j in range(i + 1, 5)]
    out = np.empty((10, 10), dtype=float)
    for a, (i, j) in enumerate(pairs):
        for b, (k, ell) in enumerate(pairs):
            out[a, b] = r[i, k] * r[j, ell] - r[i, ell] * r[j, k]
    return out


def _intertwiner_constraints(rotation_matrices, wedge_matrices):
    columns = []
    for a in range(3):
        for b in range(10):
            basis = np.zeros((3, 10), dtype=float)
            basis[a, b] = 1.0
            columns.append(basis)
    blocks = []
    for rotation, wedge in zip(rotation_matrices, wedge_matrices):
        block = np.column_stack([
            (basis @ wedge - rotation @ basis).ravel()
            for basis in columns
        ])
        blocks.append(block)
    return np.vstack(blocks)


def build_q6_d5_relative_direction_intertwiner() -> dict:
    """Compute the A5-equivariant Lambda^2(V5)->R3 map by nullspace SVD."""
    phase = build_q6_d5_phase_space_sector()
    mode_basis = np.asarray(phase["selected_mode_basis"], dtype=float)
    group = icosahedral_signed_axis_representation()
    from src.rgs_q6_equivariant_cochains import (
        build_q6_cochain_complex,
        induced_cochain_actions,
    )
    complex_data = build_q6_cochain_complex()
    rotations = np.asarray(group["rotations"], dtype=float)
    signed = np.asarray(group["signed_permutations"])
    rho_matrices = []
    wedge_matrices = []
    for s in signed:
        _p0, p1, _p2 = induced_cochain_actions(s, complex_data)
        rho = mode_basis.T @ p1 @ mode_basis
        rho_matrices.append(rho)
        wedge_matrices.append(_wedge2_matrix(rho))
    constraints = _intertwiner_constraints(rotations, wedge_matrices)
    _u, singular_values, vh = np.linalg.svd(constraints, full_matrices=False)
    nullity = int(np.count_nonzero(singular_values < 1e-8))
    if nullity != 1:
        raise ArithmeticError(
            f"expected a unique geometric-3 intertwiner, nullity={nullity}; "
            f"smallest singular values={singular_values[-4:]}"
        )
    intertwiner = vh[-1].reshape(3, 10)
    # Schur's lemma gives T T^T = c I; normalize c=1.
    gram = intertwiner @ intertwiner.T
    scale = float(np.trace(gram) / 3.0)
    if scale <= 1e-14:
        raise ArithmeticError("intertwiner has zero norm")
    intertwiner /= np.sqrt(scale)
    # Fix the otherwise arbitrary global sign deterministically.
    first = np.flatnonzero(np.abs(intertwiner.ravel()) > 1e-10)
    if len(first) and intertwiner.ravel()[first[0]] < 0:
        intertwiner = -intertwiner

    residual = max(
        float(np.linalg.norm(intertwiner @ wedge - rotation @ intertwiner))
        for rotation, wedge in zip(rotations, wedge_matrices)
    )
    orthogonality = float(np.linalg.norm(
        intertwiner @ intertwiner.T - np.eye(3)
    ))
    rank = int(np.linalg.matrix_rank(intertwiner, tol=1e-9))
    checks = {
        "unique_intertwiner_nullity_one": nullity == 1,
        "intertwiner_rank_three": rank == 3,
        "intertwiner_rows_orthonormal": orthogonality < 1e-8,
        "A5_intertwining_residual": residual < 1e-8,
        "source_five_mode_sector_verified": bool(phase["all_checks_passed"]),
    }
    return {
        "audit": "Q6_D5_WEDGE2_TO_GEOMETRIC_R3_INTERTWINER",
        "source_irrep_dimension": 5,
        "exterior_square_dimension": 10,
        "target_geometric_irrep_dimension": 3,
        "A5_group_order": 60,
        "intertwiner": intertwiner,
        "geometric_rotations": rotations,
        "A5_mode_representation": np.asarray(rho_matrices),
        "intertwiner_nullity": nullity,
        "constraint_singular_values": singular_values,
        "intertwining_residual": residual,
        "row_orthogonality_residual": orthogonality,
        "intertwiner_rank": rank,
        "all_checks_passed": bool(all(checks.values())),
        "checks": checks,
        "scope": (
            "A5-equivariant projection Lambda^2(V5)->R3; a nonlinear map "
            "from the 10D cotangent sector via q wedge p, not a linear "
            "identification V5=R3 and not yet a measure-preserving D5 map"
        ),
    }


def relative_direction_linear_map_in_p(q, intertwiner):
    """Return M(q) such that T(q wedge p)=M(q) p."""
    q = np.asarray(q, dtype=float)
    T = np.asarray(intertwiner, dtype=float)
    if q.shape != (5,) or T.shape != (3, 10):
        raise ValueError("q must have shape (5,) and intertwiner shape (3,10)")
    if not np.all(np.isfinite(q)) or not np.all(np.isfinite(T)):
        raise ValueError("q and intertwiner must be finite")
    matrix = np.empty((3, 5), dtype=float)
    for j in range(5):
        basis = np.eye(5)[j]
        wedge = np.asarray([
            q[i]*basis[k] - q[k]*basis[i]
            for i in range(5) for k in range(i+1, 5)
        ])
        matrix[:, j] = T @ wedge
    return matrix


def relative_direction_from_cotangent(q, p, intertwiner=None):
    """Return T(q wedge p) and its unit direction when nonzero.

    Under the canonical SO(2) oscillator phase rotation of (q,p), q wedge p
    is invariant. Under A5, the returned vector transforms in the geometric
    three-dimensional rotation representation.
    """
    q = np.asarray(q, dtype=float)
    p = np.asarray(p, dtype=float)
    if q.shape != (5,) or p.shape != (5,):
        raise ValueError("q and p must have shape (5,)")
    if not np.all(np.isfinite(q)) or not np.all(np.isfinite(p)):
        raise ValueError("q and p must be finite")
    T = (
        np.asarray(build_q6_d5_relative_direction_intertwiner()["intertwiner"])
        if intertwiner is None else np.asarray(intertwiner, dtype=float)
    )
    if T.shape != (3, 10):
        raise ValueError("intertwiner must have shape (3,10)")
    wedge = np.asarray(
        [q[i]*p[j] - q[j]*p[i] for i in range(5) for j in range(i+1, 5)],
        dtype=float,
    )
    vector = T @ wedge
    norm = float(np.linalg.norm(vector))
    direction = None if norm <= 1e-12 else vector / norm
    return {"wedge_coordinates": wedge, "vector": vector,
            "norm": norm, "direction": direction}

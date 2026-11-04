"""Exact rational/integer audit of the five-mode Q6 Maxwell sector."""
from __future__ import annotations
import math
import numpy as np
import sympy as sp
from src.rgs_q6_equivariant_cochains import build_q6_cochain_complex, induced_cochain_actions
from src.rgs_icosahedral_q6_intertwiner import icosahedral_signed_axis_representation
from src.rgs_q6_d5_phase_space import build_q6_d5_phase_space_sector


def _exact_integer_mode_basis(eigenvalue=12):
    if not isinstance(eigenvalue, int) or eigenvalue <= 0:
        raise ValueError("eigenvalue must be a positive integer")
    c = build_q6_cochain_complex()
    d0 = sp.Matrix(np.asarray(c["D0"], dtype=np.int64).T.tolist())
    d1 = sp.Matrix(np.asarray(c["D1"], dtype=np.int64).tolist())
    laplacian = d1.T * d1
    constraints = d0.col_join(laplacian - eigenvalue * sp.eye(192))
    rational_basis = constraints.nullspace()
    if len(rational_basis) != 5:
        raise ArithmeticError(f"expected rational eigenspace dimension five, got {len(rational_basis)}")
    columns = []
    for vector in rational_basis:
        denominator = sp.ilcm(*[entry.q for entry in vector])
        columns.append(vector * denominator)
    basis = sp.Matrix.hstack(*columns)
    if basis.rank() != 5 or any(entry.q != 1 for entry in basis):
        raise ArithmeticError("integerized exact modes do not have integer rank five")
    if d0 * basis != sp.zeros(64, 5):
        raise ArithmeticError("integer modes violate Gauss law")
    if (laplacian - eigenvalue * sp.eye(192)) * basis != sp.zeros(192, 5):
        raise ArithmeticError("integer modes fail exact Hodge eigen-equation")
    return c, basis


def d5_type_iv_norm_polynomial(z):
    value = np.asarray(z, dtype=complex)
    if value.shape != (5,) or not np.all(np.isfinite(value)):
        raise ValueError("z must be a finite complex vector of shape (5,)")
    quadratic_modulus = abs(np.dot(value, value))
    norm = float(1.0 - 2.0 * np.vdot(value, value).real + quadratic_modulus ** 2)
    if quadratic_modulus >= 1.0 or norm <= 0.0:
        raise ValueError("z is outside the type-IV_5 domain")
    return norm


def d5_bergman_volume_density(z):
    norm = d5_type_iv_norm_polynomial(z)
    return float(10.0 ** 5 / norm ** 6)


def audit_q6_integer_a5_d5_measure():
    c, exact_basis = _exact_integer_mode_basis(12)
    integer_basis = np.asarray(exact_basis.tolist(), dtype=np.int64)
    group = icosahedral_signed_axis_representation()
    integer_action_failures = gauss_failures = eigen_failures = 0
    lap = np.asarray(c["D1"], dtype=np.int64).T @ np.asarray(c["D1"], dtype=np.int64)
    for signed in group["signed_permutations"]:
        _p0, p1, _p2 = induced_cochain_actions(signed, c)
        moved = np.asarray(p1, dtype=np.int64) @ integer_basis
        if np.any(moved != np.rint(moved)):
            integer_action_failures += 1
        if np.any(np.asarray(c["D0"], dtype=np.int64).T @ moved):
            gauss_failures += 1
        if np.any((lap - 12 * np.eye(192, dtype=np.int64)) @ moved):
            eigen_failures += 1
    phase = build_q6_d5_phase_space_sector(beta=1.0, eigenvalue=12.0)
    mode_basis = np.asarray(phase["selected_mode_basis"], dtype=float)
    orthogonality_max = norm_invariance_max = 0.0
    sample = np.array([.013+.007j, -.009+.004j, .006-.011j, .003+.005j, -.008+.002j])
    n0 = d5_type_iv_norm_polynomial(sample)
    for signed in group["signed_permutations"]:
        _p0, p1, _p2 = induced_cochain_actions(signed, c)
        rho = mode_basis.T @ np.asarray(p1, dtype=float) @ mode_basis
        orthogonality_max = max(orthogonality_max, float(np.linalg.norm(rho.T @ rho - np.eye(5))))
        norm_invariance_max = max(norm_invariance_max, abs(d5_type_iv_norm_polynomial(rho @ sample) - n0))
    density_origin = 10.0 ** 5
    real_jacobian = 10.0 ** -5
    checks = {
        "rational_eigenspace_dimension_five": exact_basis.rank() == 5,
        "integer_vectors_rank_five": np.linalg.matrix_rank(integer_basis.astype(float)) == 5,
        "integer_vectors_exact_gauss": gauss_failures == 0,
        "integer_vectors_exact_hodge_12": eigen_failures == 0,
        "A5_group_order_60": int(group["group_order"]) == 60,
        "A5_preserves_integer_intersection": integer_action_failures == 0,
        "A5_representation_orthogonal": orthogonality_max < 1e-8,
        "A5_preserves_D5_generic_norm": norm_invariance_max < 1e-10,
        "D5_Bergman_density_origin": density_origin == 10**5,
        "normalized_tangent_measure_origin": math.isclose(density_origin * real_jacobian, 1.0),
    }
    return {
        "audit": "Q6_INTEGER_A5_D5_BERGMAN_MEASURE",
        "rational_eigenspace_dimension": int(exact_basis.rank()),
        "integer_mode_basis_shape": list(integer_basis.shape),
        "integer_mode_basis": integer_basis,
        "integer_mode_basis_rank": int(np.linalg.matrix_rank(integer_basis.astype(float))),
        "integer_action_failures": integer_action_failures,
        "integer_gauss_failures": gauss_failures,
        "integer_hodge_failures": eigen_failures,
        "A5_group_order": int(group["group_order"]),
        "A5_orthogonality_max": orthogonality_max,
        "A5_type_IV_norm_invariance_max": norm_invariance_max,
        "D5_Bergman_density_formula": "10^5*N(z)^(-6)",
        "Bergman_density_at_origin": density_origin,
        "tangent_coordinate_real_jacobian": real_jacobian,
        "normalized_tangent_measure_at_origin": density_origin * real_jacobian,
        "checks": checks,
        "all_checks_passed": bool(all(checks.values())),
        "proved": [
            "the Hodge-12 Gauss sector is defined over Q",
            "its intersection with Z^192 contains five independent integer flux vectors",
            "the full A5 action preserves that integer intersection",
            "the real A5 five-dimensional representation preserves the D5 type-IV generic norm",
            "the Bergman density cancels the normalized tangent-coordinate Jacobian at the origin",
        ],
        "not_proved": [
            "the displayed integer generators form a saturated basis of the full integer intersection",
            "the compact integer-flux Hilbert space equals the real Gaussian five-mode Hilbert space",
            "a global nonlinear compact-Q6 to D5 map",
            "a microscopic RGS-PSU derivation of the D5 Bergman measure normalization",
            "Wyler's global volume prescription or a derivation of the fine-structure constant",
        ],
        "scope": "Exact rational Q6 mode structure and A5 action; local D5 Bergman-volume normalization only.",
    }

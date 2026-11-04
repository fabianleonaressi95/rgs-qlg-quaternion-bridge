"""A5-invariant angular measure induced by the Q6/D5 cotangent sphere.

The map (q,p) -> T(q wedge p) is A5-equivariant. The normalized map from
the unit phase sphere S^9 therefore pushes its invariant surface measure to
an A5-invariant probability measure on S^2, provided its zero set has
surface measure zero. The A5 harmonic character calculation shows that
every such measure has the uniform S^2 moments through degree five.
"""
from __future__ import annotations

import math
import numpy as np

from src.rgs_q6_d5_sphere_intertwiner import (
    build_q6_d5_relative_direction_intertwiner,
    relative_direction_linear_map_in_p,
)


_A5_CLASS_SIZES = (1, 15, 20, 12, 12)
_A5_ROTATION_ANGLES = (
    0.0,
    math.pi,
    2.0 * math.pi / 3.0,
    2.0 * math.pi / 5.0,
    4.0 * math.pi / 5.0,
)


def so3_spherical_harmonic_character(degree: int, angle: float) -> float:
    """Character of the degree-l real SO(3) spherical-harmonic module."""
    if not isinstance(degree, (int, np.integer)) or degree < 0:
        raise ValueError("degree must be a nonnegative integer")
    if not np.isfinite(angle):
        raise ValueError("angle must be finite")
    if abs(math.sin(angle / 2.0)) < 1e-14:
        return float(2 * int(degree) + 1)
    return float(
        math.sin((int(degree) + 0.5) * angle) / math.sin(angle / 2.0)
    )


def a5_invariant_harmonic_multiplicity(degree: int) -> int:
    """Dimension of the A5-fixed subspace in degree-l spherical harmonics.

    Uses the class equation (1/60) sum_C |C| chi_l(C). The five classes
    have rotation angles 0, pi, 2pi/3, 2pi/5 and 4pi/5.
    """
    value = sum(
        size * so3_spherical_harmonic_character(degree, angle)
        for size, angle in zip(_A5_CLASS_SIZES, _A5_ROTATION_ANGLES)
    ) / 60.0
    rounded = int(round(value))
    if abs(value - rounded) > 2e-10:
        raise ArithmeticError(
            f"A5 character average is not integral at degree {degree}: {value}"
        )
    return rounded


def a5_invariant_harmonic_spectrum(max_degree: int = 12) -> list[int]:
    """Return fixed-space dimensions for degrees 0,...,max_degree."""
    if not isinstance(max_degree, (int, np.integer)) or max_degree < 0:
        raise ValueError("max_degree must be a nonnegative integer")
    return [
        a5_invariant_harmonic_multiplicity(degree)
        for degree in range(int(max_degree) + 1)
    ]


def build_q6_d5_a5_angular_measure_audit(max_degree: int = 12) -> dict:
    """Audit the exact symmetry consequence for the normalized S9 pushforward.

    Let X=(Q,P) be distributed by normalized surface measure on S^9 in
    R^5 x R^5. The A5 action is diag(rho_g,rho_g), hence preserves this
    measure. The map F(Q,P)=T(Q wedge P) is equivariant. Away from F=0,
    n=F/||F|| is equivariant too, so its pushforward is A5-invariant.

    The invariant-harmonic multiplicities are zero in degrees 1,...,5.
    Consequently every A5-invariant probability measure on S^2, including
    this pushforward whenever defined almost everywhere, agrees with uniform
    S^2 on all polynomial moments of total degree at most five. Degree six
    is the first possible anisotropic harmonic and has a one-dimensional
    A5-fixed subspace. Symmetry alone does not determine its coefficient.
    """
    if max_degree < 6:
        raise ValueError("max_degree must be at least six for the first anisotropy audit")

    bridge = build_q6_d5_relative_direction_intertwiner()
    rho = np.asarray(bridge["A5_mode_representation"], dtype=float)
    rotations = np.asarray(bridge["geometric_rotations"], dtype=float)
    T = np.asarray(bridge["intertwiner"], dtype=float)

    if rho.shape != (60, 5, 5) or rotations.shape != (60, 3, 3):
        raise ArithmeticError("expected 60 A5 matrices in dimensions five and three")

    orthogonal_residual = max(
        float(np.linalg.norm(matrix.T @ matrix - np.eye(5)))
        for matrix in rho
    )
    geometric_orthogonal_residual = max(
        float(np.linalg.norm(matrix.T @ matrix - np.eye(3)))
        for matrix in rotations
    )
    determinant_residual = max(
        abs(float(np.linalg.det(matrix)) - 1.0)
        for matrix in rotations
    )
    intertwining_residual = float(bridge["intertwining_residual"])
    witness_q = np.arange(1.0, 6.0)
    witness_matrix = relative_direction_linear_map_in_p(witness_q, T)
    witness_rank = int(np.linalg.matrix_rank(witness_matrix, tol=1e-9))

    multiplicities = a5_invariant_harmonic_spectrum(max_degree)
    low_degree_zero = multiplicities[1:6] == [0, 0, 0, 0, 0]
    first_nontrivial_degree_six = multiplicities[6] == 1
    checks = {
        "A5_group_order_60": len(rho) == 60,
        "five_mode_action_orthogonal": orthogonal_residual < 2e-8,
        "geometric_action_SO3": (
            geometric_orthogonal_residual < 2e-8
            and determinant_residual < 2e-8
        ),
        "wedge_intertwiner_equivariant": intertwining_residual < 2e-8,
        "nonzero_polynomial_map_witness": witness_rank > 0,
        "no_invariant_harmonics_degrees_1_to_5": low_degree_zero,
        "first_invariant_harmonic_degree_6_is_unique": first_nontrivial_degree_six,
        "source_bridge_checks_passed": bool(bridge["all_checks_passed"]),
    }
    return {
        "audit": "Q6_D5_A5_ANGULAR_PUSHFORWARD_MEASURE_V1",
        "phase_space": "unit sphere S^9 in R^5_q x R^5_p",
        "phase_measure": "normalized O(10)-invariant surface probability measure",
        "map": "n(q,p)=T(q wedge p)/||T(q wedge p)|| where nonzero",
        "group_order": 60,
        "A5_invariant_spherical_harmonic_multiplicities": multiplicities,
        "first_possible_angular_anisotropy_degree": 6,
        "degree_six_invariant_dimension": multiplicities[6],
        "five_mode_orthogonality_residual": orthogonal_residual,
        "geometric_SO3_orthogonality_residual": geometric_orthogonal_residual,
        "geometric_determinant_residual": determinant_residual,
        "wedge_intertwining_residual": intertwining_residual,
        "zero_locus_witness_q": witness_q,
        "zero_locus_witness_matrix_rank": witness_rank,
        "zero_locus_has_surface_measure_zero": bool(witness_rank > 0),
        "checks": checks,
        "all_checks_passed": bool(all(checks.values())),
        "theorem": (
            "The map T(q wedge p) is a nonzero polynomial map (the recorded "
            "rank witness proves nontriviality), so its common zero locus has "
            "S9 surface measure zero. Its normalized pushforward is therefore "
            "an A5-invariant probability measure. Since the A5-fixed subspace "
            "of degree-l spherical harmonics vanishes for 1<=l<=5, its angular "
            "moments through degree five equal those of uniform S2."
        ),
        "proved_by_finite_group_character": [
            "A5-invariant harmonic multiplicities are zero for degrees 1 through 5",
            "the degree-six invariant harmonic subspace is one-dimensional",
            "an A5-invariant angular probability measure has no anisotropic harmonic below degree six",
        ],
        "conditional_or_open": [
            "a complete global classification of the zero locus T(q wedge p)=0",
            "the degree-six coefficient of the continuous pushforward measure",
            "equality of the continuous pushforward measure with the 32-point Q6 measure",
            "a measure-preserving global embedding into the D5 Lie ball",
        ],
        "scope_warning": (
            "The result is a symmetry theorem for the normalized S9 phase-shell "
            "measure, not a derivation of the physical Liouville measure on an "
            "unbounded phase space and not a continuum QED or alpha prediction."
        ),
    }

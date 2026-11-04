"""Character-theoretic proof of the A5 type of the Q6 top coexact mode.

The key calculation avoids choosing an eigenbasis inside the degenerate
five-dimensional eigenspace. Spectral projectors are polynomials in the
integer Hodge matrices, so their equivariant traces are basis-independent.
"""
from __future__ import annotations

import math
import numpy as np

from src.rgs_q6_equivariant_cochains import (
    build_q6_cochain_complex,
    induced_cochain_actions,
)
from src.rgs_icosahedral_q6_intertwiner import (
    icosahedral_signed_axis_representation,
)

SQRT5 = math.sqrt(5.0)
PHI = (1.0 + SQRT5) / 2.0
CLASS_SIZES = np.array([1.0, 15.0, 20.0, 12.0, 12.0])
CLASS_TRACES = np.array([3.0, -1.0, 0.0, PHI, (1.0 - SQRT5) / 2.0])
IRREPS = {
    "1": np.array([1.0, 1.0, 1.0, 1.0, 1.0]),
    "3": np.array([3.0, -1.0, 0.0, PHI, (1.0 - SQRT5) / 2.0]),
    "3_prime": np.array([3.0, -1.0, 0.0, (1.0 - SQRT5) / 2.0, PHI]),
    "4": np.array([4.0, 0.0, 1.0, -1.0, -1.0]),
    "5": np.array([5.0, 1.0, -1.0, 0.0, 0.0]),
}


def _spectral_projector(matrix: np.ndarray, eigenvalue: float,
                        spectrum: tuple[float, ...]) -> np.ndarray:
    """Lagrange polynomial projector for a symmetric matrix with known spectrum."""
    a = np.asarray(matrix, dtype=float)
    out = np.eye(a.shape[0])
    for other in spectrum:
        if not math.isclose(other, eigenvalue, rel_tol=0.0, abs_tol=1e-12):
            out = out @ ((a - other * np.eye(a.shape[0])) /
                         (eigenvalue - other))
    return 0.5 * (out + out.T)


def _class_index(rotation_trace: float) -> int:
    t = float(rotation_trace)
    index = int(np.argmin(np.abs(CLASS_TRACES - t)))
    if abs(float(CLASS_TRACES[index]) - t) > 1e-7:
        raise ArithmeticError(f"unrecognized A5 rotation trace {t}")
    return index


def q6_top_coexact_a5_character_report(tolerance: float = 2e-7) -> dict:
    """Recover the exact A5 character of the lambda=12 coexact eigenspace.

    On C^1(Q6), Delta_1 has spectrum {2,4,6,8,10,12}; its lambda=12
    eigenspace has dimension six. Exactly one of those six modes is exact,
    namely d0 applied to the unique lambda=12 eigenvector of Delta_0.
    Subtracting its character gives the coexact five-dimensional character.

    For each group element g:
      chi_coex(g) = tr(P_12^(1) rho_1(g)) - tr(P_12^(0) rho_0(g)).
    Both projectors are exact Lagrange polynomials in the Hodge matrices.
    """
    c = build_q6_cochain_complex()
    d0 = np.asarray(c["D0"], dtype=float)
    d1 = np.asarray(c["D1"], dtype=float)
    delta0 = d0.T @ d0
    delta1 = d0 @ d0.T + d1.T @ d1

    p1 = _spectral_projector(delta1, 12.0, (2., 4., 6., 8., 10., 12.))
    p0 = _spectral_projector(delta0, 12.0, (0., 2., 4., 6., 8., 10., 12.))
    rank_p1 = int(round(float(np.trace(p1))))
    rank_p0 = int(round(float(np.trace(p0))))

    group = icosahedral_signed_axis_representation()
    class_values = [[] for _ in range(5)]
    full_values = [[] for _ in range(5)]
    exact_values = [[] for _ in range(5)]
    max_equivariance = 0.0
    max_projector_idempotence = max(
        float(np.linalg.norm(p1 @ p1 - p1)),
        float(np.linalg.norm(p0 @ p0 - p0)),
    )
    for rotation, signed in zip(group["rotations"], group["signed_permutations"]):
        p0g, p1g, _p2g = induced_cochain_actions(signed, c)
        ci = _class_index(np.trace(rotation))
        full_trace = float(np.trace(p1 @ p1g))
        exact_trace = float(np.trace(p0 @ p0g))
        coexact_trace = full_trace - exact_trace
        full_values[ci].append(full_trace)
        exact_values[ci].append(exact_trace)
        class_values[ci].append(coexact_trace)
        max_equivariance = max(
            max_equivariance,
            float(np.linalg.norm(delta0 @ p0g - p0g @ delta0)),
            float(np.linalg.norm(delta1 @ p1g - p1g @ delta1)),
            float(np.linalg.norm(p0 @ p0g - p0g @ p0)),
            float(np.linalg.norm(p1 @ p1g - p1g @ p1)),
        )

    character = np.array([np.mean(v) for v in class_values])
    full_character = np.array([np.mean(v) for v in full_values])
    exact_character = np.array([np.mean(v) for v in exact_values])
    expected = IRREPS["5"]
    char_residual = float(np.max(np.abs(character - expected)))
    inner_product = float(np.dot(CLASS_SIZES, character * character) / 60.0)
    multiplicities = {
        name: float(np.dot(CLASS_SIZES, character * chi) / 60.0)
        for name, chi in IRREPS.items()
    }
    checks = {
        "full_lambda12_dimension_six": rank_p1 == 6,
        "exact_lambda12_dimension_one": rank_p0 == 1,
        "coexact_lambda12_dimension_five": rank_p1 - rank_p0 == 5,
        "character_is_5_1_minus1_0_0": char_residual < tolerance,
        "character_inner_product_is_one": abs(inner_product - 1.0) < tolerance,
        "irreducible_five_occurs_once": abs(multiplicities["5"] - 1.0) < tolerance,
        "other_irreps_absent": all(
            abs(multiplicities[name]) < tolerance
            for name in ("1", "3", "3_prime", "4")
        ),
        "projectors_are_idempotent": max_projector_idempotence < 1e-7,
        "group_action_commutes_with_hodge_projectors": max_equivariance < 1e-7,
    }
    return {
        "audit": "Q6_TOP_COEXACT_A5_CHARACTER_THEOREM",
        "group_order": int(group["group_order"]),
        "class_sizes": CLASS_SIZES.astype(int).tolist(),
        "class_rotation_traces": CLASS_TRACES.tolist(),
        "lambda12_full_one_cochain_dimension": rank_p1,
        "lambda12_exact_one_cochain_dimension": rank_p0,
        "lambda12_coexact_dimension": rank_p1 - rank_p0,
        "full_lambda12_character": full_character,
        "exact_lambda12_character": exact_character,
        "coexact_lambda12_character": character,
        "expected_irreducible_five_character": expected,
        "character_inner_product": inner_product,
        "irreducible_multiplicities": multiplicities,
        "residuals": {
            "character": char_residual,
            "projector_idempotence": max_projector_idempotence,
            "group_hodge_commutator": max_equivariance,
        },
        "checks": checks,
        "all_checks_passed": bool(all(checks.values())),
        "proof_method": (
            "Basis-independent character extraction by Lagrange spectral "
            "projectors: tr(P12 on C1 rho1(g)) minus tr(P12 on C0 rho0(g)); "
            "the exact lambda-12 1-form is d0 of the unique lambda-12 "
            "vertex eigenmode."
        ),
        "logical_conclusion": (
            "The top coexact Q6 eigenspace is, as a real A5-module, one "
            "copy of the irreducible five-dimensional representation V5."
        ),
        "scope": (
            "Finite-complex representation theorem. It identifies the "
            "A5-module type, not a preferred basis, a global D5 embedding, "
            "a measure pushforward, or a physical alpha normalization."
        ),
    }

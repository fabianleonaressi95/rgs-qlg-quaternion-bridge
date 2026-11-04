"""Exact integer-fiber quotient of Q6 fluxes onto the five D5 modes.

The top coexact eigenspace is computed over Q from the integer Hodge matrix.
Its rational mode map sends the rank-129 integral cycle lattice to a
rank-five lattice. The kernel intersected with Z^129 has rank 124.
The pointwise Wilson weight does not descend to this quotient; the positive
fiber-summed weight does, by summing over the complete integer kernel.
"""
from __future__ import annotations

import math
import numpy as np
import sympy as sp

from src.rgs_q6_euclidean_time_wilson import build_q6_euclidean_time_data
from src.rgs_q6_u1_character_electric_flux import build_q6_integer_cycle_basis
from src.rgs_q6_u1_transfer_operator import (
    q6_u1_electric_transfer_eigenvalues,
    validate_q6_integer_flux_states,
)


def _integerize_rational_vector(vector: sp.Matrix) -> sp.Matrix:
    den = sp.ilcm(*[int(x.q) for x in vector])
    values = [int(x * den) for x in vector]
    common = int(sp.igcd(*values)) if any(values) else 1
    return sp.Matrix([x // common for x in values])


def build_q6_d5_rational_integer_quotient(data=None):
    """Return exact rational D5 mode map and independent integer kernel vectors.

    The 5D eigenspace is ker(D1.T*D1 - 12 I) over Q. Its basis need not be
    orthonormal: any rational basis gives the same kernel and quotient.
    """
    d = build_q6_euclidean_time_data(2) if data is None else data
    d0 = sp.Matrix(np.asarray(d["D0"], dtype=np.int64).tolist())
    d1 = sp.Matrix(np.asarray(d["D1"], dtype=np.int64).tolist())
    cycles_np = build_q6_integer_cycle_basis(d).astype(np.int64)
    cycles = sp.Matrix(cycles_np.tolist())
    hodge = d1.T * d1
    top_basis = (hodge - 12 * sp.eye(192)).nullspace()
    if len(top_basis) != 5:
        raise ArithmeticError(f"expected rational top eigenspace dimension 5, got {len(top_basis)}")
    v = sp.Matrix.hstack(*top_basis)
    if d0.T * v != sp.zeros(64, 5):
        raise ArithmeticError("top eigenspace is not in the Coulomb cycle space")
    if hodge * v != 12 * v:
        raise ArithmeticError("exact top-mode eigen-equation failed")

    mode_map = cycles * v  # 129 x 5, rational and full column rank
    if mode_map.rank() != 5:
        raise ArithmeticError("cycle-to-D5 rational mode map does not have rank five")

    rational_kernel = mode_map.T.nullspace()
    if len(rational_kernel) != 124:
        raise ArithmeticError("expected a 124-dimensional rational kernel")
    integer_witnesses = [_integerize_rational_vector(x) for x in rational_kernel]
    kernel_rows = sp.Matrix.vstack(*[x.T for x in integer_witnesses])
    if kernel_rows.rank() != 124 or kernel_rows * mode_map != sp.zeros(124, 5):
        raise ArithmeticError("integer kernel witnesses failed exact verification")

    # Find a small-support exact witness, using the first prefix with rank five.
    witness = None
    support_size = None
    for width in range(6, 32):
        prefix = mode_map[:width, :]
        if prefix.rank() == 5:
            local_kernel = prefix.T.nullspace()
            if local_kernel:
                local = _integerize_rational_vector(local_kernel[0])
                witness = sp.zeros(129, 1)
                for i, value in enumerate(local):
                    witness[i] = value
                support_size = int(sum(1 for x in witness if x != 0))
                break
    if witness is None:
        witness = integer_witnesses[0]
        support_size = int(sum(1 for x in witness if x != 0))
    if (witness.T * mode_map) != sp.zeros(1, 5):
        raise ArithmeticError("integer kernel witness does not project to zero")

    return {
        "cycles": cycles_np,
        "top_mode_basis_rational": v,
        "mode_map_rational": mode_map,
        "integer_kernel_witnesses": kernel_rows,
        "kernel_witness": np.asarray(witness, dtype=np.int64).reshape(129),
        "kernel_witness_support": support_size,
        "data": d,
    }


def q6_d5_fiber_summed_electric_weight(
    representative,
    kernel_flux_coefficients,
    beta_time: float,
    data=None,
) -> float:
    """Finite positive partial sum over supplied integer kernel translations.

    For a complete enumeration of K_Z this becomes the exact effective
    quotient weight. This routine intentionally accepts only an explicit
    finite subset; it never labels a truncation as the complete fiber sum.
    """
    beta = float(beta_time)
    if not math.isfinite(beta) or beta <= 0:
        raise ValueError("beta_time must be finite and positive")
    d = build_q6_euclidean_time_data(2) if data is None else data
    c = np.asarray(representative)
    k = np.asarray(kernel_flux_coefficients)
    if c.shape != (129,) or not np.all(c == np.rint(c)):
        raise ValueError("representative must be an integer vector of shape (129,)")
    if k.ndim != 2 or k.shape[1] != 129 or not np.all(k == np.rint(k)):
        raise ValueError("kernel_flux_coefficients must be an integer matrix with 129 columns")
    states = (c[None, :] + k).astype(np.int64) @ build_q6_integer_cycle_basis(d)
    states = validate_q6_integer_flux_states(states, d)
    weights = q6_u1_electric_transfer_eigenvalues(states, beta, d)
    return float(np.sum(weights))


def audit_q6_d5_integer_fiber_quotient(beta_time=80.0):
    """Exact quotient-rank audit and a witness that Wilson weight needs fibers."""
    beta = float(beta_time)
    if not math.isfinite(beta) or beta <= 0:
        raise ValueError("beta_time must be finite and positive")
    q = build_q6_d5_rational_integer_quotient()
    cycles = q["cycles"]
    mode_map = q["mode_map_rational"]
    kernel_rows = q["integer_kernel_witnesses"]
    witness = q["kernel_witness"]
    d = q["data"]

    flux = witness @ cycles
    states = np.vstack((np.zeros(192, dtype=np.int64), np.asarray(flux, dtype=np.int64)))
    states = validate_q6_integer_flux_states(states, d)
    weights = q6_u1_electric_transfer_eigenvalues(states, beta, d)
    witness_weight = float(weights[1])
    vacuum_weight = float(weights[0])

    # The quotient is the image lattice in Q^5; Z^129/K_Z is torsion-free
    # because it embeds into the additive group Q^5.
    checks = {
        "cycle_lattice_rank_129": cycles.shape == (129, 192),
        "exact_top_eigenspace_dimension_5": q["top_mode_basis_rational"].shape == (192, 5),
        "exact_mode_map_rank_5": mode_map.rank() == 5,
        "exact_integer_kernel_rank_124": kernel_rows.rank() == 124,
        "all_kernel_witnesses_project_to_zero": kernel_rows * mode_map == sp.zeros(124, 5),
        "kernel_witness_is_nonzero_integer": bool(np.any(witness != 0)),
        "kernel_witness_satisfies_gauss_law": bool(
            np.max(np.abs(np.asarray(d["D0"], dtype=np.int64).T @ flux)) == 0
        ),
        "vacuum_weight_is_one": math.isclose(vacuum_weight, 1.0, rel_tol=1e-13),
        "wilson_weight_not_constant_on_kernel_fiber": 0.0 < witness_weight < 1.0,
    }
    return {
        "audit": "Q6_D5_EXACT_INTEGER_FIBER_QUOTIENT",
        "quotient": {
            "source_integer_cycle_rank": 129,
            "target_rational_mode_rank": 5,
            "integer_kernel_rank": 124,
            "quotient_rank": 5,
            "quotient_group": "Z^129 / (Z^129 intersect ker_Q(Pi_5))",
            "isomorphism": "the image of the rational mode map as a rank-five subgroup of Q^5",
            "mode_map_basis": "exact rational basis of ker(D1.T D1 - 12 I); basis-independent kernel",
            "integer_kernel_witness_count": int(kernel_rows.rows),
            "witness_support_size": int(q["kernel_witness_support"]),
        },
        "fiber_weight_test": {
            "beta_time": beta,
            "vacuum_transfer_weight": vacuum_weight,
            "nonzero_kernel_translation_transfer_weight": witness_weight,
            "strict_noninvariance": bool(witness_weight < vacuum_weight),
            "meaning": "the pointwise Wilson weight does not descend to the five-mode quotient",
        },
        "effective_quotient_weight": {
            "definition": "W_eff([c]) = sum_{k in K_Z} product_e I_|((c+k)C)_e|(beta_t)/I_0(beta_t)",
            "positivity": "strictly positive",
            "convergence": (
                "For fixed beta, I_m(beta)/I_0(beta) <= (beta/2)^m/m!. "
                "The injective integer cycle map has a positive norm lower bound, "
                "so the factorial edge-flux decay makes the sum over the rank-124 "
                "integer kernel absolutely convergent."
            ),
            "numerically_evaluated": False,
            "reason": "the full 124-dimensional fiber sum requires a separate controlled lattice enumeration",
        },
        "checks": checks,
        "all_checks_passed": bool(all(checks.values())),
        "conclusion": (
            "The exact rational D5 mode map defines a rank-five quotient of the "
            "integer Q6 cycle lattice with rank-124 integer kernel. The Wilson "
            "weight is not constant on kernel fibers, so pointwise quotienting is "
            "invalid. The mathematically correct positive quotient weight is the "
            "sum over the complete integer fiber; its convergence follows from "
            "factorial Bessel decay."
        ),
        "scope": {
            "established": [
                "exact rational top eigenspace and mode-map rank",
                "rank-124 integer kernel and rank-five torsion-free quotient",
                "explicit integer kernel translation changing the Wilson weight",
                "definition and absolute convergence of the fiber-summed electric weight",
            ],
            "open": [
                "controlled numerical evaluation of the full rank-124 fiber sum",
                "spatial Wilson potential and canonical q-sector on the quotient",
                "full ten-real-dimensional D5 phase-space measure",
                "identification of a Wilson-induced Hua generic-norm cutoff",
                "global Q6-to-D5 map and any physical 3+1D reduction",
            ],
        },
    }

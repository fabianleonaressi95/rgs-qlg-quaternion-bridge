"""Exact integral splitting certificate for the Q6 annihilator/gauge tori.

All lattice bases and pairings are computed over ZZ/QQ. No floating-point
nullspace is used. The routine certifies the abstract five-torus quotient
by constructing dual unimodular character/cocharacter bases.
"""
from __future__ import annotations

import numpy as np
import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp

from src.rgs_q6_d5_integer_fiber_quotient import (
    build_q6_d5_rational_integer_quotient,
)
from src.rgs_q6_d5_saturated_kernel_bounds import (
    build_q6_d5_saturated_integer_kernel,
)


def _snf(matrix: sp.Matrix):
    """Return exact Smith form and unimodular transforms U,V with U*A*V=S."""
    dm = DomainMatrix.from_Matrix(sp.Matrix(matrix)).convert_to(sp.ZZ)
    smith, left, right = smith_normal_decomp(dm)
    s, u, v = smith.to_Matrix(), left.to_Matrix(), right.to_Matrix()
    a = sp.Matrix(matrix)
    if u * a * v != s:
        raise ArithmeticError("exact Smith identity U*A*V=S failed")
    return s, u, v


def _rank_from_smith(s: sp.Matrix) -> int:
    return sum(1 for i in range(min(s.rows, s.cols)) if s[i, i] != 0)


def build_q6_integral_gauge_quotient_splitting(data=None) -> dict:
    """Construct the exact H/G cocharacter complement and dual character basis.

    H = ker(L : T^192 -> T^124), G = im(D : T^64 -> T^192).
    The output includes a Z-basis of ker_Z(L), a gauge-adapted unimodular
    basis in that lattice, and the dual five character vectors.
    """
    q = build_q6_d5_rational_integer_quotient(data)
    sat = build_q6_d5_saturated_integer_kernel(data)
    d = q["data"] if data is None else data

    cycles = sp.Matrix(np.asarray(q["cycles"], dtype=np.int64).tolist())
    d0 = sp.Matrix(np.asarray(d["D0"], dtype=np.int64).tolist())
    if cycles.shape != (129, 192) or d0.shape != (192, 64):
        raise ArithmeticError("unexpected Q6 cycle/incidence matrix dimensions")
    if cycles * d0 != sp.zeros(129, 64):
        raise ArithmeticError("cycle-gradient identity C*D0=0 failed")

    kernel_basis = sp.Matrix(sat["kernel_basis"])
    lattice = kernel_basis * cycles
    if lattice.shape != (124, 192) or lattice.rank() != 124:
        raise ArithmeticError("annihilator lattice must have rank 124")
    s_l, u_l, v_l = _snf(lattice)
    rank_l = _rank_from_smith(s_l)
    if rank_l != 124:
        raise ArithmeticError("annihilator lattice rank is not 124")
    invariants_l = [abs(int(s_l[i, i])) for i in range(rank_l)]
    if invariants_l != [1] * 124:
        raise ArithmeticError("annihilator row lattice is not primitive")

    # The last 68 columns of V_L are a saturated Z-basis of ker_Z(L).
    q_h = v_l[:, rank_l:]
    if q_h.shape != (192, 68) or lattice * q_h != sp.zeros(124, 68):
        raise ArithmeticError("Smith kernel basis for H failed")

    # Express the 64 incidence columns in the integral H-cocharacter basis.
    # Select 68 independent rows of Q_H and solve the resulting square system.
    pivot_rows = list(q_h.T.rref()[1])
    if len(pivot_rows) != 68:
        raise ArithmeticError("H cocharacter basis does not have rank 68")
    q_pivot = q_h[pivot_rows, :]
    d_pivot = d0[pivot_rows, :]
    x_g = q_pivot.inv() * d_pivot
    if any(entry.q != 1 for entry in x_g):
        raise ArithmeticError("gauge generators are not integral in H coordinates")
    x_g = sp.Matrix(x_g)
    if q_h * x_g != d0:
        raise ArithmeticError("gauge incidence map does not reconstruct exactly")
    if x_g.rank() != 63:
        raise ArithmeticError("gauge cocharacter rank is not 63")

    s_g, u_g, v_g = _snf(x_g)
    rank_g = _rank_from_smith(s_g)
    invariants_g = [abs(int(s_g[i, i])) for i in range(rank_g)]
    if rank_g != 63 or invariants_g != [1] * 63:
        raise ArithmeticError("gauge cocharacter lattice is not primitive in H")

    # U_G^{-1} is a unimodular basis in H coordinates, adapted to im(X_G).
    u_g_inv = u_g.inv()
    if u_g * x_g * v_g != s_g:
        raise ArithmeticError("gauge Smith identity failed on reconstruction")
    adapted_coordinates = u_g_inv
    adapted_edge_basis = q_h * adapted_coordinates
    gauge_edge_basis = adapted_edge_basis[:, :63]
    quotient_cocharacters = adapted_edge_basis[:, 63:]
    if lattice * adapted_edge_basis != sp.zeros(124, 68):
        raise ArithmeticError("adapted cocharacter basis leaves ker(L)")
    if q_h * x_g != d0:
        raise ArithmeticError("gauge inclusion verification failed")

    # A = M^T maps Z^129 -> Z^5 after clearing a common denominator.
    # Its Smith right transform splits the cycle-coordinate lattice into
    # five quotient representatives plus the saturated rank-124 kernel.
    a = sp.Matrix(sat["integer_map"])
    s_a, u_a, v_a = _snf(a)
    rank_a = _rank_from_smith(s_a)
    if rank_a != 5:
        raise ArithmeticError("five-mode integer map rank is not five")
    cycle_quotient_representatives = v_a[:, :5]
    cycle_kernel_basis = v_a[:, 5:]
    if a * cycle_kernel_basis != sp.zeros(5, 124):
        raise ArithmeticError("Smith kernel of five-mode map failed")
    if cycle_quotient_representatives.shape != (129, 5):
        raise ArithmeticError("five quotient representatives have wrong shape")

    # Cycle representatives are characters of H/G. Their pairing with the
    # five complementary cocharacters must be a unimodular 5x5 matrix.
    quotient_characters = cycles.T * cycle_quotient_representatives
    pairing = quotient_characters.T * quotient_cocharacters
    if pairing.shape != (5, 5) or abs(int(pairing.det())) != 1:
        raise ArithmeticError("quotient character/cocharacter pairing is not unimodular")

    checks = {
        "cycle_gradient_identity": cycles * d0 == sp.zeros(129, 64),
        "annihilator_rank_124": rank_l == 124,
        "annihilator_saturated": invariants_l == [1] * 124,
        "H_cocharacter_rank_68": q_h.rank() == 68,
        "gauge_inclusion_exact": q_h * x_g == d0,
        "gauge_rank_63": rank_g == 63,
        "gauge_saturated_in_H": invariants_g == [1] * 63,
        "adapted_basis_rank_68": adapted_coordinates.det() in (1, -1),
        "adapted_basis_annihilated_by_L": lattice * adapted_edge_basis == sp.zeros(124, 68),
        "quotient_character_pairing_unimodular": abs(int(pairing.det())) == 1,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"integral splitting certificate failed: {checks}")

    return {
        "cycles": cycles,
        "incidence": d0,
        "annihilator_lattice": lattice,
        "H_cocharacter_basis": q_h,
        "gauge_coordinates_in_H": x_g,
        "gauge_adapted_coordinates": adapted_coordinates,
        "gauge_adapted_edge_basis": adapted_edge_basis,
        "gauge_edge_basis": gauge_edge_basis,
        "quotient_cocharacters": quotient_cocharacters,
        "cycle_quotient_representatives": cycle_quotient_representatives,
        "quotient_characters": quotient_characters,
        "quotient_pairing": pairing,
        "smith_invariants_H_constraints": invariants_l,
        "smith_invariants_gauge_inclusion": invariants_g,
        "checks": checks,
    }


def audit_q6_integral_gauge_quotient_splitting(data=None) -> dict:
    """Return a JSON-compatible exact certificate for H/G ~= T^5."""
    out = build_q6_integral_gauge_quotient_splitting(data)
    return {
        "audit": "Q6_INTEGRAL_GAUGE_QUOTIENT_SPLITTING",
        "groups": {
            "annihilator": "connected T^68",
            "vertex_gauge": "connected T^63",
            "quotient": "T^5",
            "quotient_character_lattice": "Z^129 / K ~= image of the exact five-mode integer map",
        },
        "matrix_shapes": {
            "H_cocharacter_basis": list(out["H_cocharacter_basis"].shape),
            "gauge_coordinates_in_H": list(out["gauge_coordinates_in_H"].shape),
            "gauge_adapted_edge_basis": list(out["gauge_adapted_edge_basis"].shape),
            "quotient_cocharacters": list(out["quotient_cocharacters"].shape),
            "quotient_characters": list(out["quotient_characters"].shape),
            "quotient_pairing": list(out["quotient_pairing"].shape),
        },
        "smith_invariants_H_constraints": out["smith_invariants_H_constraints"],
        "smith_invariants_gauge_inclusion": out["smith_invariants_gauge_inclusion"],
        "quotient_pairing": [[int(x) for x in out["quotient_pairing"].row(i)]
                             for i in range(5)],
        "checks": out["checks"],
        "all_checks_passed": all(out["checks"].values()),
        "evidence_status": {
            "proved_by_exact_integer_arithmetic": [
                "H has a saturated rank-68 cocharacter lattice",
                "the 63-dimensional gauge cocharacter lattice is primitive in H",
                "the five quotient cocharacters extend the gauge lattice to a unimodular basis",
                "the five-mode quotient characters pair unimodularly with the five quotient cocharacters",
            ],
            "not_claimed": [
                "Wilson fibre marginal equals the D5 Hua density",
                "a continuum or physical 3+1D limit",
            ],
        },
    }

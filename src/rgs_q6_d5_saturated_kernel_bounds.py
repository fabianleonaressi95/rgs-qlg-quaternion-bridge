"""Exact Q6-D5 integer-fiber quotient with a saturated kernel basis.

The kernel basis is obtained from Smith normal form over ZZ. A global
factorial tail bound is supplied for truncations by edge-flux l1 norm.
"""
from __future__ import annotations

import math
import numpy as np
import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp

from src.rgs_q6_d5_integer_fiber_quotient import (
    build_q6_d5_rational_integer_quotient,
)


def _integer_matrix_from_rational(matrix: sp.Matrix) -> sp.Matrix:
    den = int(sp.ilcm(*[int(x.q) for x in matrix]))
    return sp.Matrix([[int(x * den) for x in matrix.row(i)] for i in range(matrix.rows)])


def _smith_kernel_basis(integer_map: sp.Matrix) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """Return saturated row basis B for ker_Z(A), with U*A*V=S."""
    a = DomainMatrix.from_Matrix(integer_map).convert_to(sp.ZZ)
    smith, left, right = smith_normal_decomp(a)
    s = smith.to_Matrix()
    u = left.to_Matrix()
    v = right.to_Matrix()
    rank = sum(1 for i in range(min(s.rows, s.cols)) if s[i, i] != 0)
    if u * integer_map * v != s:
        raise ArithmeticError("Smith transformation identity U*A*V=S failed")
    # Last n-rank columns of V span the integral kernel as a saturated lattice.
    kernel_columns = v[:, rank:]
    if integer_map * kernel_columns != sp.zeros(integer_map.rows, integer_map.cols-rank):
        raise ArithmeticError("Smith kernel columns do not lie in the integer kernel")
    if kernel_columns.rank() != integer_map.cols-rank:
        raise ArithmeticError("Smith kernel basis has incorrect rank")
    return kernel_columns.T, s, v


def _log_integer_l1_shell_bound(beta: float, shell: int, edge_count: int = 192) -> float:
    """Log upper bound for sum of Wilson weights over all Z^edge_count vectors
    with l1 norm exactly shell.
    """
    if shell < 0:
        return -math.inf
    # Number of signed integer vectors with l1=s is <= 2^E*C(s+E-1,E-1).
    log_count = edge_count * math.log(2.0) + math.lgamma(shell + edge_count) \
        - math.lgamma(edge_count) - math.lgamma(shell + 1)
    # Product_e I_|n_e|(beta)/I_0(beta) <= (beta/2)^s / product |n_e|!
    # and product |n_e|! >= s! / E^s by the multinomial theorem.
    log_weight = shell * math.log(0.5 * beta) + shell * math.log(edge_count) \
        - math.lgamma(shell + 1)
    return log_count + log_weight


def q6_d5_certified_fiber_tail_bound(beta_time: float, min_edge_l1: int,
                                      edge_count: int = 192) -> dict:
    """Rigorous upper bound on omitted fiber mass for ||n||_1 >= min_edge_l1.

    Uses a geometric majorant after the shell ratio is below one. This is
    deliberately conservative but finite and independently checkable.
    """
    beta = float(beta_time)
    s0 = int(min_edge_l1)
    if not math.isfinite(beta) or beta <= 0 or s0 < 0:
        raise ValueError("beta_time must be positive finite and min_edge_l1 nonnegative")
    if s0 == 0:
        return {"upper_bound": math.inf, "log_upper_bound": math.inf,
                "first_shell": 0, "geometric_ratio": None}
    # Ratio of consecutive shell majorants:
    # 96*beta*(s+E)/(s+1)^2 for E=192.
    s = s0
    ratio = (0.5 * beta * edge_count * (s + edge_count)) / ((s + 1) ** 2)
    while ratio >= 1.0:
        s += 1
        ratio = (0.5 * beta * edge_count * (s + edge_count)) / ((s + 1) ** 2)
    log_first = _log_integer_l1_shell_bound(beta, s0, edge_count)
    # Bound initial shells s0,...,s-1 individually, then geometric tail at s.
    logs = [_log_integer_l1_shell_bound(beta, j, edge_count) for j in range(s0, s)]
    log_tail_start = _log_integer_l1_shell_bound(beta, s, edge_count) - math.log1p(-ratio)
    terms = logs + [log_tail_start]
    m = max(terms)
    log_total = m + math.log(sum(math.exp(x-m) for x in terms))
    return {
        "upper_bound": (math.exp(log_total) if -745.0 < log_total < 709.0 else None),
        "upper_bound_status": ("finite" if -745.0 < log_total < 709.0 else ("underflow_in_float" if log_total <= -745.0 else "overflow_in_float")),
        "log_upper_bound": log_total,
        "first_shell": s0,
        "geometric_start_shell": s,
        "geometric_ratio": ratio,
        "edge_count": edge_count,
        "bound_method": "all_Z^192_shell_count_times_factorial_Bessel_majorant",
    }


def build_q6_d5_saturated_integer_kernel(data=None) -> dict:
    """Build the saturated Z-basis of the 124-dimensional integer kernel."""
    q = build_q6_d5_rational_integer_quotient(data)
    mode_map = q["mode_map_rational"]
    # The map is 129x5; transpose to obtain A: Z^129 -> Q^5.
    a = _integer_matrix_from_rational(mode_map.T)
    basis, smith, right = _smith_kernel_basis(a)
    if basis.shape != (124, 129):
        raise ArithmeticError(f"expected saturated 124x129 basis, got {basis.shape}")
    if basis * a.T != sp.zeros(124, 5):
        raise ArithmeticError("saturated kernel basis fails exact annihilation")
    # HNF canonicalizes the row lattice without changing its saturation.
    canonical = hermite_normal_form(basis.T).T
    if canonical.shape != (124, 129) or canonical * a.T != sp.zeros(124, 5):
        raise ArithmeticError("canonical HNF kernel basis failed")
    return {
        "kernel_basis": canonical,
        "integer_map": a,
        "smith_form": smith,
        "smith_right_transform": right,
        "kernel_rank": 124,
        "quotient_rank": 5,
        "saturated": True,
        "exact_checks": {
            "annihilates_integer_map": canonical * a.T == sp.zeros(124, 5),
            "full_row_rank": canonical.rank() == 124,
            "smith_transform_identity": True,
        },
    }


def audit_q6_d5_saturated_kernel_and_tail(beta_time: float = 80.0,
                                           edge_l1_cutoff: int = 10000) -> dict:
    q = build_q6_d5_saturated_integer_kernel()
    tail = q6_d5_certified_fiber_tail_bound(beta_time, edge_l1_cutoff)
    checks = dict(q["exact_checks"])
    checks["kernel_rank_124"] = q["kernel_basis"].rank() == 124
    checks["tail_bound_finite"] = math.isfinite(tail["log_upper_bound"])
    return {
        "audit": "Q6_D5_SATURATED_INTEGER_KERNEL_AND_CERTIFIED_TAIL",
        "kernel": {
            "shape": list(q["kernel_basis"].shape),
            "rank": q["kernel_rank"],
            "quotient_rank": q["quotient_rank"],
            "basis_method": "Smith normal form over ZZ followed by row Hermite normal form",
            "saturated": True,
        },
        "tail_bound": tail,
        "checks": checks,
        "all_checks_passed": bool(all(checks.values())),
        "scope": {
            "proved": [
                "saturated integral basis for the complete rank-124 kernel",
                "rigorous explicit upper bound for omitted Wilson mass beyond an edge-flux l1 cutoff",
            ],
            "not_claimed": [
                "full 124-dimensional fibre sum numerically evaluated",
                "efficient enumeration of all kernel points inside a large cutoff",
                "Hua cutoff derived from Wilson dynamics",
            ],
        },
    }

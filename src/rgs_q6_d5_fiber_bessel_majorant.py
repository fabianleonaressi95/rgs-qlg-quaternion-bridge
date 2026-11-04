"""Analytic majorants for the Q6-D5 Wilson integer-fiber marginal.

This module does NOT enumerate the rank-124 fiber. It supplies an explicit,
finite logarithmic upper bound for its complete positive Wilson sum. The
bound is intentionally conservative and is reported separately from finite
partial sums; it must not be mistaken for a numerical evaluation of W_eff.
"""
from __future__ import annotations

import math
import numpy as np
import sympy as sp
from scipy.special import ive, logsumexp

from src.rgs_q6_d5_integer_fiber_quotient import (
    build_q6_d5_rational_integer_quotient,
    q6_d5_fiber_summed_electric_weight,
)


def _pivot_minor_inverse_norm(cycles: np.ndarray):
    """Return exact ||A^-1||_1 for a nonsingular 129-column minor A.

    Coefficients are row vectors: k = n_P A^-1. The induced matrix
    1-norm (maximum absolute column sum) therefore controls ||k||_infinity
    by ||n_P||_infinity.
    """
    c = sp.Matrix(np.asarray(cycles, dtype=np.int64).tolist())
    if c.rows != 129 or c.rank() != 129:
        raise ArithmeticError("cycle coefficient matrix must have exact row rank 129")
    _, pivots = c.rref()
    if len(pivots) != 129:
        raise ArithmeticError("could not select 129 independent edge columns")
    a = c[:, list(pivots)]
    inv = a.inv()
    column_sums = [sum(abs(inv[i, j]) for i in range(129)) for j in range(129)]
    norm_one = max(column_sums)
    return norm_one, tuple(int(i) for i in pivots)


def q6_d5_fiber_bessel_majorant(beta_time: float, representative=None, quotient=None):
    r"""Bound the complete fiber sum by a convergent one-dimensional series.

    For n_e integer and beta>0,
       I_|n_e|(beta)/I_0(beta) <= (beta/2)^|n_e|/|n_e|!.
    If N=||n||_1, multiplication and the multinomial theorem give
       prod_e ratio_e <= (96 beta)^N/N!   (192 edges).

    Let C be the 129x192 integral cycle matrix and choose an invertible
    129x129 column minor A. For n=(c+k)C and k in Z^129,
       ||k||_infty <= ||A^-1||_1 N + ||c||_infty.
    Hence the number of fiber points with flux norm N is no greater than
       (2(M N + ||c||_infty)+1)^129.
    Summing this box bound gives an explicit convergent majorant.
    """
    beta = float(beta_time)
    if not math.isfinite(beta) or beta <= 0:
        raise ValueError("beta_time must be finite and positive")
    q = build_q6_d5_rational_integer_quotient() if quotient is None else quotient
    c0 = np.zeros(129, dtype=np.int64) if representative is None else np.asarray(representative)
    if c0.shape != (129,) or not np.all(np.isfinite(c0)) or not np.all(c0 == np.rint(c0)):
        raise ValueError("representative must be a finite integer vector of shape (129,)")
    c0 = c0.astype(np.int64)
    cycles = np.asarray(q["cycles"], dtype=np.int64)
    M, pivot_columns = _pivot_minor_inverse_norm(cycles)
    c_norm = int(np.max(np.abs(c0)))
    rank = 129
    A = 96.0 * beta

    def log_shell_bound(N: int) -> float:
        # Use the real-valued upper bound 2R+1 >= 2 floor(R)+1.
        # Keeping this smooth bound is essential for the ratio-test tail.
        radius = float(M * N + c_norm)
        log_count = rank * math.log(2.0 * radius + 1.0)
        if N == 0:
            return log_count
        return log_count + N * math.log(A) - math.lgamma(N + 1)

    # For x_N=2(MN+c)+1, B_(N+1)/B_N is bounded by
    # (x_(N+1)/x_N)^129 * A/(N+1), a decreasing function of N.
    n0 = 0
    while True:
        x0 = 2.0 * (float(M * n0) + c_norm) + 1.0
        x1 = 2.0 * (float(M * (n0 + 1)) + c_norm) + 1.0
        ratio = (x1 / x0) ** rank * A / (n0 + 1)
        if ratio < 1.0:
            break
        n0 += 1
        if n0 > 2_000_000:
            raise ArithmeticError("failed to find a geometric-tail starting index")

    finite_logs = np.asarray([log_shell_bound(n) for n in range(n0)], dtype=float)
    log_finite = float(logsumexp(finite_logs)) if finite_logs.size else -math.inf
    log_tail = log_shell_bound(n0) - math.log1p(-ratio)
    log_total = float(np.logaddexp(log_finite, log_tail))
    return {
        "beta_time": beta,
        "cycle_rank": 129,
        "kernel_rank": 124,
        "edge_count": 192,
        "pivot_columns": list(pivot_columns),
        "exact_inverse_one_norm": str(M),
        "representative_infinity_norm": c_norm,
        "bessel_product_bound": "(96*beta_time)^N/N!, N=||flux||_1",
        "coefficient_box_count_bound": "(2*(M*N+||c||_infinity)+1)^129",
        "geometric_tail_start_N": n0,
        "geometric_tail_ratio_bound": float(ratio),
        "log_complete_fiber_sum_upper_bound": log_total,
        "complete_fiber_sum_upper_bound": (
            float(math.exp(log_total)) if log_total < math.log(np.finfo(float).max) else None
        ),
        "upper_bound_overflowed_float": bool(log_total >= math.log(np.finfo(float).max)),
        "status": "ANALYTIC_CONVERGENT_MAJORANT_NOT_FULL_FIBER_EVALUATION",
    }


def audit_q6_d5_fiber_majorant(beta_time=2.0):
    """Check the Bessel coefficient inequality and produce the fiber majorant."""
    beta = float(beta_time)
    if not math.isfinite(beta) or beta <= 0:
        raise ValueError("beta_time must be finite and positive")
    orders = np.arange(0, 81, dtype=int)
    ratios = ive(orders, beta) / ive(0, beta)
    rhs = np.asarray([(beta / 2.0) ** int(m) / math.factorial(int(m)) for m in orders])
    residual = float(np.max(ratios - rhs))
    inequality_ok = bool(np.all(ratios <= rhs * (1.0 + 5e-13)))
    q = build_q6_d5_rational_integer_quotient()
    majorant = q6_d5_fiber_bessel_majorant(beta, quotient=q)
    checks = {
        "Bessel_factorial_bound_orders_0_to_80": inequality_ok,
        "majorant_is_finite_in_log_space": math.isfinite(
            majorant["log_complete_fiber_sum_upper_bound"]
        ),
        "geometric_tail_ratio_below_one": majorant["geometric_tail_ratio_bound"] < 1.0,
        "exact_cycle_minor_selected": len(majorant["pivot_columns"]) == 129,
    }
    return {
        "audit": "Q6_D5_FIBER_BESSEL_MAJORANT",
        "beta_time": beta,
        "bessel_inequality_max_residual": residual,
        "majorant": majorant,
        "checks": checks,
        "all_checks_passed": bool(all(checks.values())),
        "scope": {
            "proved_by_analytic_bound": [
                "factorial Bessel decay",
                "explicit coefficient-box counting bound from an exact cycle minor",
                "convergence of the complete integer-fiber sum",
                "a finite logarithmic upper bound for the full marginal",
            ],
            "not_computed": [
                "complete rank-124 fiber sum",
                "sharp or practically narrow numerical enclosure",
                "spatial Wilson marginal",
                "D5 Hua generic-norm cutoff",
            ],
        },
    }

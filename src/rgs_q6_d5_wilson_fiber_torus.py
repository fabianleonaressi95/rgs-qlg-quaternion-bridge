"""Exact Q6-D5 Wilson fiber marginal on its Smith-parameterized annihilator.

The Fourier identity and Smith data are exact. Scrambled Sobol output is
only a diagnostic estimate; it is never reported as a certified interval.
"""
from __future__ import annotations

import math
import numpy as np
import sympy as sp
from scipy.special import i0e
from scipy.stats import qmc
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp

from src.rgs_q6_d5_integer_fiber_quotient import build_q6_d5_rational_integer_quotient
from src.rgs_q6_d5_saturated_kernel_bounds import build_q6_d5_saturated_integer_kernel


def smith_torus_parameterization(lattice_basis):
    """Return exact SNF coordinates for the annihilator of a full-row lattice."""
    lattice = sp.Matrix(lattice_basis)
    if lattice.rows >= lattice.cols or lattice.rank() != lattice.rows:
        raise ValueError("lattice basis must have full row rank and fewer rows than columns")
    dm = DomainMatrix.from_Matrix(lattice).convert_to(sp.ZZ)
    diagonal, left, right = smith_normal_decomp(dm)
    smith, u, v = diagonal.to_Matrix(), left.to_Matrix(), right.to_Matrix()
    if u * lattice * v != smith:
        raise ArithmeticError("invalid Smith decomposition")
    rank = lattice.rows
    factors = [abs(int(smith[i, i])) for i in range(rank)]
    if any(value == 0 for value in factors):
        raise ArithmeticError("zero Smith invariant")
    return {
        "smith_form": smith,
        "left_transform": u,
        "right_transform": v,
        "invariants": factors,
        "rank": rank,
        "ambient_dimension": lattice.cols,
        "torus_dimension": lattice.cols - rank,
    }


def q6_d5_wilson_annihilator_torus_integral(
    representative=None, beta_time: float = 2.0, power: int = 8,
    scrambles: int = 4, seed: int = 20261104, quotient=None,
) -> dict:
    """Evaluate a reproducible QMC diagnostic for the exact fiber integral."""
    beta = float(beta_time)
    if not math.isfinite(beta) or beta <= 0:
        raise ValueError("beta_time must be finite and positive")
    if not isinstance(power, int) or not 4 <= power <= 20:
        raise ValueError("power must be an integer in [4, 20]")
    if not isinstance(scrambles, int) or not 2 <= scrambles <= 64:
        raise ValueError("scrambles must be an integer in [2, 64]")

    quotient = build_q6_d5_rational_integer_quotient() if quotient is None else quotient
    saturated = build_q6_d5_saturated_integer_kernel()
    cycles = sp.Matrix(np.asarray(quotient["cycles"], dtype=np.int64).tolist())
    lattice = saturated["kernel_basis"] * cycles
    snf = smith_torus_parameterization(lattice)

    c = np.zeros(129, dtype=np.int64) if representative is None else np.asarray(representative)
    if c.shape != (129,) or not np.all(np.isfinite(c)) or not np.all(c == np.rint(c)):
        raise ValueError("representative must be an integer vector of shape (129,)")
    c = c.astype(np.int64)
    n0 = np.asarray((sp.Matrix(1, 129, list(map(int, c))) * cycles).tolist(),
                    dtype=np.float64).reshape(192)

    rank, ambient = snf["rank"], snf["ambient_dimension"]
    free_dim = snf["torus_dimension"]
    factors = snf["invariants"]
    discrete_indices = [i for i, factor in enumerate(factors) if factor > 1]
    qmc_dim = len(discrete_indices) + free_dim
    v_exact = snf["right_transform"]
    max_v = max(abs(int(value)) for value in v_exact)
    if max_v > 2**52 or max(factors) > 2**52:
        raise ArithmeticError("Smith coordinates exceed exact float sampling range")
    v = np.asarray(v_exact.tolist(), dtype=np.float64)
    d = np.asarray(factors, dtype=np.int64)

    log_i0 = beta + math.log(float(i0e(beta)))
    log_scale = ambient * (beta - log_i0)
    estimates = []
    for scramble in range(scrambles):
        points = qmc.Sobol(d=qmc_dim, scramble=True, seed=seed + scramble).random_base2(power)
        z = np.zeros((len(points), ambient), dtype=np.float64)
        for column, smith_index in enumerate(discrete_indices):
            count = int(d[smith_index])
            z[:, smith_index] = 2.0 * math.pi * np.floor(
                points[:, column] * count
            ) / count
        z[:, rank:] = 2.0 * math.pi * points[:, len(discrete_indices):]
        x = np.remainder(z @ v.T, 2.0 * math.pi)
        amplitude = np.exp(np.minimum(beta * np.cos(x).sum(axis=1) - log_scale, 0.0))
        phase = x @ n0
        estimates.append(np.mean(amplitude * (np.cos(phase) - 1j * np.sin(phase))))

    estimates = np.asarray(estimates, dtype=np.complex128)
    mean = complex(estimates.mean())
    real_se = float(estimates.real.std(ddof=1) / math.sqrt(scrambles))
    imag_se = float(estimates.imag.std(ddof=1) / math.sqrt(scrambles))
    log_abs = log_scale + math.log(abs(mean)) if abs(mean) else None
    real_estimate = None
    if mean.real > 0 and log_scale + math.log(mean.real) < math.log(np.finfo(float).max):
        real_estimate = math.exp(log_scale + math.log(mean.real))

    return {
        "audit": "Q6_D5_WILSON_ANNIHILATOR_TORUS_FIBER_MARGINAL",
        "exact_identity": {
            "formula": "W_eff([c]) = integral_H exp(beta*sum cos(x_e))/I0(beta)^192 * exp(-i*(cC).x) dH",
            "kernel_rank": rank,
            "ambient_dimension": ambient,
            "continuous_torus_dimension": free_dim,
            "smith_invariants": [str(value) for value in factors],
            "finite_component_count": str(math.prod(factors)),
            "smith_nontrivial_indices": discrete_indices,
            "right_transform_max_abs_entry": max_v,
        },
        "numerical_qmc": {
            "beta_time": beta,
            "power": power,
            "samples_per_scramble": 2**power,
            "scrambles": scrambles,
            "seed": seed,
            "qmc_dimension": qmc_dim,
            "discrete_coordinates": len(discrete_indices),
            "continuous_coordinates": free_dim,
            "scaled_complex_mean": [mean.real, mean.imag],
            "per_scramble_scaled_means": [[z.real, z.imag] for z in estimates],
            "scaled_real_standard_error": real_se,
            "scaled_imag_standard_error": imag_se,
            "log_absolute_estimate": log_abs,
            "real_estimate": real_estimate,
            "status": "QMC_ESTIMATE_NOT_CERTIFIED_INTERVAL",
        },
        "scope": {
            "proved": [
                "the complete integer-fiber sum is the Fourier coefficient on the exact annihilator",
                "the annihilator dimension is ambient_dimension minus kernel_rank",
                "Smith invariants determine its finite component group",
            ],
            "open": [
                "a rigorous numerical error enclosure",
                "a derived pushforward to the D5 Hua measure",
                "a physical 3+1D reduction or fine-structure constant prediction",
            ],
        },
    }

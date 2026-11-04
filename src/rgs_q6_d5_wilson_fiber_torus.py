"""Exact annihilator-torus representation of the Q6-D5 Wilson fiber sum.

The rank-124 fiber lattice is embedded in Z^192 by the cycle map. Its
annihilator is a (possibly disconnected) subgroup of T^192 of dimension 68.
The Fourier identity below is exact. Sobol integration is only diagnostic.
"""
from __future__ import annotations
import math
import numpy as np
import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp
from scipy.special import i0e
from scipy.stats import qmc

from src.rgs_q6_d5_integer_fiber_quotient import build_q6_d5_rational_integer_quotient
from src.rgs_q6_d5_saturated_kernel_bounds import build_q6_d5_saturated_integer_kernel


def smith_torus_parameterization(lattice_basis):
    """For row lattice L, return SNF data parameterizing {x: L x=0 mod 2pi}."""
    lmat = sp.Matrix(lattice_basis)
    if lmat.rows >= lmat.cols or lmat.rank() != lmat.rows:
        raise ValueError("lattice basis must have full row rank and fewer rows than columns")
    dm = DomainMatrix.from_Matrix(lmat).convert_to(sp.ZZ)
    smith, left, right = smith_normal_decomp(dm)
    s, u, v = smith.to_Matrix(), left.to_Matrix(), right.to_Matrix()
    if u * lmat * v != s:
        raise ArithmeticError("Smith identity U*L*V=S failed")
    rank = lmat.rows
    invariants = [abs(int(s[i, i])) for i in range(rank)]
    if any(d == 0 for d in invariants):
        raise ArithmeticError("unexpected zero Smith invariant")
    # smith_normal_decomp returns unimodular transforms by construction;
    # avoid an expensive determinant of the full 192x192 transform.
    return {
        "smith_form": s, "left_transform": u, "right_transform": v,
        "invariants": invariants, "rank": rank,
        "ambient_dimension": lmat.cols,
        "torus_dimension": lmat.cols-rank,
    }


def q6_d5_wilson_annihilator_torus_integral(
    representative=None, beta_time: float = 2.0, power: int = 10,
    scrambles: int = 8, seed: int = 20261104, quotient=None,
) -> dict:
    r"""Estimate the exact annihilator-torus integral for the full fiber sum.

    If L = K_Z C subset Z^192 and n0=c C, then
      sum_{k in K_Z} w((c+k)C)
        = integral_{H_L} exp(beta sum_e cos(x_e))/I0(beta)^192
                    exp(-i n0.x) dmu_Haar(x),
    where H_L={x in T^192: ell.x=0 mod 2pi for all ell in L}.
    H_L has dimension 68; its finite components are encoded by SNF invariants.
    """
    beta = float(beta_time)
    if not math.isfinite(beta) or beta <= 0:
        raise ValueError("beta_time must be finite and positive")
    if not isinstance(power, int) or power < 4 or power > 20:
        raise ValueError("power must be an integer between 4 and 20")
    if not isinstance(scrambles, int) or scrambles < 2 or scrambles > 64:
        raise ValueError("scrambles must be an integer between 2 and 64")

    q = build_q6_d5_rational_integer_quotient() if quotient is None else quotient
    sat = build_q6_d5_saturated_integer_kernel()
    cycles = sp.Matrix(np.asarray(q["cycles"], dtype=np.int64).tolist())
    kernel_basis = sat["kernel_basis"]
    lattice = kernel_basis * cycles
    snf = smith_torus_parameterization(lattice)

    c = np.zeros(129, dtype=np.int64) if representative is None else np.asarray(representative)
    if c.shape != (129,) or not np.all(np.isfinite(c)) or not np.all(c == np.rint(c)):
        raise ValueError("representative must be an integer vector of shape (129,)")
    c = c.astype(np.int64)
    n0 = np.asarray((sp.Matrix(1, 129, list(map(int, c))) * cycles).tolist(),
                    dtype=np.float64).reshape(192)

    v_exact = snf["right_transform"]
    max_v = max(abs(int(x)) for x in v_exact)
    if max_v > 2**52:
        raise ArithmeticError("Smith transform exceeds exact float integer range; QMC parameterization unsafe")
    v = np.asarray(v_exact.tolist(), dtype=np.float64)
    if max(snf["invariants"]) > 2**52:
        raise ArithmeticError("Smith invariant exceeds exact float sampling range")
    d = np.asarray(snf["invariants"], dtype=np.int64)
    rank, ambient = snf["rank"], snf["ambient_dimension"]
    free_dim = snf["torus_dimension"]
    log_i0 = beta + math.log(float(i0e(beta)))
    log_scale = ambient * (beta-log_i0)
    means = []
    for j in range(scrambles):
        sampler = qmc.Sobol(d=ambient, scramble=True, seed=seed+j)
        u = sampler.random_base2(power)
        discrete = np.floor(u[:, :rank] * d[None, :]).astype(np.int64)
        z = np.empty((u.shape[0], ambient), dtype=np.float64)
        z[:, :rank] = 2.0*math.pi*discrete/d[None, :]
        z[:, rank:] = 2.0*math.pi*u[:, rank:]
        x = np.remainder(z @ v.T, 2.0*math.pi)
        log_amp = beta*np.sum(np.cos(x), axis=1)-log_scale
        amp = np.exp(np.minimum(log_amp, 0.0))
        phase = x @ n0
        vals = amp*(np.cos(phase)-1j*np.sin(phase))
        means.append(np.mean(vals))
    means = np.asarray(means, dtype=np.complex128)
    mean = complex(np.mean(means))
    re_se = float(np.std(means.real, ddof=1)/math.sqrt(scrambles))
    im_se = float(np.std(means.imag, ddof=1)/math.sqrt(scrambles))
    log_abs = log_scale+math.log(abs(mean)) if abs(mean)>0 else -math.inf
    real_est = None
    if mean.real>0 and log_scale+math.log(mean.real)<math.log(np.finfo(float).max):
        real_est = math.exp(log_scale+math.log(mean.real))
    return {
        "audit": "Q6_D5_WILSON_ANNIHILATOR_TORUS_FIBER_MARGINAL",
        "exact_identity": {
            "formula": "W_eff([c]) = integral_{H_L} exp(beta*sum_e cos(x_e))/I0(beta)^192 * exp(-i*(c*C).x) dmu_Haar(x)",
            "fiber_lattice": "L=K_Z*C subset Z^192",
            "annihilator": "H_L={x in T^192: ell dot x=0 mod 2pi for every ell in L}",
            "kernel_rank": rank,
            "ambient_dimension": ambient,
            "continuous_torus_dimension": free_dim,
            "smith_invariants": [str(x) for x in snf["invariants"]],
            "finite_component_count": str(math.prod(snf["invariants"])),
            "right_transform_max_abs_entry": max_v,
            "normalization": "uniform over all finite components and normalized Haar on T^68",
        },
        "numerical_qmc": {
            "beta_time": beta, "power": power, "samples_per_scramble": 2**power,
            "scrambles": scrambles, "seed": seed, "log_scale_removed": log_scale,
            "scaled_complex_mean": [mean.real, mean.imag],
            "scaled_real_standard_error": re_se,
            "scaled_imag_standard_error": im_se,
            "log_absolute_estimate": log_abs, "real_estimate": real_est,
            "relative_scramble_error": re_se/abs(mean.real) if mean.real else None,
            "status": "QMC_ESTIMATE_NOT_CERTIFIED_INTERVAL",
        },
        "scope": {
            "proved_exactly": [
                "complete rank-124 fiber sum equals a Fourier coefficient integrated over the exact annihilator subgroup",
                "annihilator subgroup dimension is 68",
                "Smith invariants encode all finite connected components",
            ],
            "not_claimed": [
                "rigorous numerical error interval",
                "identification with D5 Hua measure or cutoff",
                "physical 3+1D reduction or fine-structure constant prediction",
            ],
        },
    }

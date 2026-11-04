"""Exact five-torus Fourier representation of the Q6-D5 Wilson fiber sum.

For K = ker_Z(A) with A an integer 5x129 matrix, Fourier orthogonality
turns the complete sum over c+K into a five-dimensional torus integral.
The QMC routine is diagnostic numerical integration, not a rigorous
error enclosure.
"""
from __future__ import annotations
import math
import numpy as np
import sympy as sp
from scipy.special import i0e
from scipy.stats import qmc
from src.rgs_q6_d5_integer_fiber_quotient import build_q6_d5_rational_integer_quotient


def _integerize_matrix(matrix: sp.Matrix) -> sp.Matrix:
    den = int(sp.ilcm(*[int(x.q) for x in matrix]))
    return sp.Matrix([[int(matrix[i, j] * den) for j in range(matrix.cols)]
                      for i in range(matrix.rows)])


def _torus_frequency_data(integer_map, cycles):
    a = sp.Matrix(integer_map)
    c = sp.Matrix(np.asarray(cycles, dtype=np.int64).tolist())
    if a.cols != c.rows or a.rank() != 5 or c.rank() != c.rows:
        raise ValueError("expected rank-five integer map and full-row-rank cycle matrix")
    if any(not x.is_Integer for x in a):
        raise ValueError("integer_map must have integer entries")
    return a, c, a * c


def q6_d5_wilson_fiber_torus_integral(
    representative=None, beta_time: float = 2.0, power: int = 12,
    scrambles: int = 8, seed: int = 20261104, quotient=None,
) -> dict:
    r"""Estimate W_eff([c]) using its exact T^5 Fourier-integral identity.

    W_eff([c]) = integral_T5 exp(beta sum_e cos((A C)_e dot t))
                 / I_0(beta)^E * exp(-i (c A^T).t) dmu_Haar(t).
    The identity is exact; only the QMC integration is approximate.
    """
    beta = float(beta_time)
    if not math.isfinite(beta) or beta <= 0:
        raise ValueError("beta_time must be finite and positive")
    if not isinstance(power, int) or power < 4 or power > 22:
        raise ValueError("power must be an integer between 4 and 22")
    if not isinstance(scrambles, int) or scrambles < 2 or scrambles > 64:
        raise ValueError("scrambles must be an integer between 2 and 64")
    q = build_q6_d5_rational_integer_quotient() if quotient is None else quotient
    a = _integerize_matrix(q["mode_map_rational"].T)
    a, cmat, frequency = _torus_frequency_data(a, q["cycles"])
    n = cmat.rows
    c = np.zeros(n, dtype=np.int64) if representative is None else np.asarray(representative)
    if c.shape != (n,) or not np.all(np.isfinite(c)) or not np.all(c == np.rint(c)):
        raise ValueError(f"representative must be an integer vector of shape ({n},)")
    c = c.astype(np.int64)
    freq = np.asarray(frequency.tolist(), dtype=np.float64)
    character = np.asarray((sp.Matrix(1, n, list(map(int, c))) * a.T).tolist(),
                           dtype=np.float64).reshape(5)
    edges = cmat.cols
    log_i0 = beta + math.log(float(i0e(beta)))
    log_scale = edges * (beta - log_i0)
    scaled_means = []
    sample_count = 2 ** power
    for j in range(scrambles):
        sampler = qmc.Sobol(d=5, scramble=True, seed=seed + j)
        points = sampler.random_base2(power)
        angles = (2.0 * math.pi) * points
        edge_angles = angles @ freq
        log_scaled = beta * np.sum(np.cos(edge_angles), axis=1) - log_scale
        amplitude = np.exp(np.minimum(log_scaled, 0.0))
        phase = angles @ character
        values = amplitude * (np.cos(phase) - 1j * np.sin(phase))
        scaled_means.append(np.mean(values))
    means = np.asarray(scaled_means, dtype=np.complex128)
    mean = complex(np.mean(means))
    re_se = float(np.std(means.real, ddof=1) / math.sqrt(scrambles))
    im_se = float(np.std(means.imag, ddof=1) / math.sqrt(scrambles))
    scaled_abs = abs(mean)
    log_abs_estimate = log_scale + math.log(scaled_abs) if scaled_abs > 0 else -math.inf
    real_estimate = None
    if mean.real > 0 and log_scale + math.log(mean.real) < math.log(np.finfo(float).max):
        real_estimate = math.exp(log_scale + math.log(mean.real))
    return {
        "audit": "Q6_D5_WILSON_FIVE_TORUS_FOURIER_MARGINAL",
        "exact_identity": {
            "formula": "W_eff([c]) = integral_T5 exp(beta*sum_e cos((A*C)_e dot t))/I0(beta)^E * exp(-i*(c*A.T) dot t) dmu_Haar(t)",
            "integer_map_shape": list(a.shape),
            "cycle_matrix_shape": list(cmat.shape),
            "frequency_matrix_shape": list(frequency.shape),
            "kernel_characterization": "K_Z = {k in Z^129 : k*A.T = 0}",
            "normalization": "normalized Haar measure (2*pi)^-5 dt",
        },
        "numerical_qmc": {
            "beta_time": beta, "representative_linf": int(np.max(np.abs(c))),
            "power": power, "samples_per_scramble": sample_count,
            "scrambles": scrambles, "seed": seed, "log_scale_removed": log_scale,
            "scaled_complex_mean": [mean.real, mean.imag],
            "scaled_real_standard_error": re_se,
            "scaled_imag_standard_error": im_se,
            "log_absolute_estimate": log_abs_estimate,
            "real_estimate": real_estimate,
            "relative_scramble_error": re_se / abs(mean.real) if mean.real != 0 else None,
            "status": "QMC_ESTIMATE_NOT_CERTIFIED_INTERVAL",
        },
        "scope": {
            "proved_exactly": [
                "complete rank-124 integer-fiber sum equals a normalized five-torus Fourier coefficient",
                "formula is independent of representative modulo K_Z",
                "no integer-fiber truncation appears in the integral identity",
            ],
            "not_claimed": [
                "rigorous numerical error interval for QMC estimate",
                "identification with D5 Hua measure or cutoff",
                "physical 3+1D reduction or fine-structure constant prediction",
            ],
        },
    }

"""Independent scrambled-Sobol cross-check of the Q6-D5 degree-six coefficient."""
from __future__ import annotations

import math
import numpy as np
from scipy.special import ndtri
from scipy.stats import qmc

from src.rgs_icosahedral_q6_intertwiner import icosahedral_signed_axis_representation
from src.rgs_q6_d5_a5_degree_six import (
    _invariant_harmonic_values,
    _uniform_sphere_harmonic_rms,
    build_q6_d5_degree_six_coefficient_audit,
)
from src.rgs_q6_d5_sphere_intertwiner import build_q6_d5_relative_direction_intertwiner
from src.rgs_q6_relative_sphere import build_q6_a5_relative_sphere_map


def _sobol_direction_batch(power: int, seed: int, intertwiner: np.ndarray) -> np.ndarray:
    """Map one independently scrambled Sobol net in (0,1)^10 to S^2."""
    sampler = qmc.Sobol(d=10, scramble=True, seed=int(seed))
    uniform = sampler.random_base2(power)
    eps = np.finfo(float).eps
    gaussian = ndtri(np.clip(uniform, eps, 1.0 - eps))
    q, p = gaussian[:, :5], gaussian[:, 5:]
    pairs = [(i, j) for i in range(5) for j in range(i + 1, 5)]
    wedge = np.column_stack([
        q[:, i] * p[:, j] - q[:, j] * p[:, i] for i, j in pairs
    ])
    vectors = wedge @ np.asarray(intertwiner, dtype=float).T
    norms = np.linalg.norm(vectors, axis=1)
    if np.any(~np.isfinite(norms)) or np.any(norms <= np.finfo(float).tiny):
        raise ArithmeticError("Sobol sample reached the zero locus or became non-finite")
    return vectors / norms[:, None]


def build_q6_d5_degree_six_sobol_crosscheck(
    power: int = 12,
    scrambles: int = 8,
    seed: int = 20261105,
) -> dict:
    """Compare independent scrambled-Sobol and Gaussian estimates.

    Each Sobol replicate is a separately scrambled 2**power point net.
    The across-scramble standard error is an empirical randomized-QMC error
    estimate, not a deterministic quadrature error bound.
    """
    if not isinstance(power, int) or power < 8:
        raise ValueError("power must be an integer >= 8")
    if not isinstance(scrambles, int) or scrambles < 4:
        raise ValueError("scrambles must be an integer >= 4")
    if not isinstance(seed, int):
        raise ValueError("seed must be an integer")

    group = icosahedral_signed_axis_representation()
    rotations = np.asarray(group["rotations"], dtype=float)
    bridge = build_q6_d5_relative_direction_intertwiner()
    if not bridge["all_checks_passed"]:
        raise ArithmeticError("Q6-D5 intertwiner audit failed")
    intertwiner = np.asarray(bridge["intertwiner"], dtype=float)
    directions = np.asarray(
        build_q6_a5_relative_sphere_map()["unique_directions"], dtype=float
    )
    if directions.shape != (32, 3):
        raise ArithmeticError("finite Q6 design must contain 32 directions")

    rms = _uniform_sphere_harmonic_rms(rotations)
    if not np.isfinite(rms) or rms <= 1e-14:
        raise ArithmeticError("degree-six invariant harmonic has zero norm")
    finite_raw = _invariant_harmonic_values(directions, rotations)
    finite = float(np.mean(finite_raw) / rms)

    seed_rng = np.random.default_rng(seed)
    replicate_means = []
    for _ in range(scrambles):
        scramble_seed = int(seed_rng.integers(0, 2**31 - 1))
        sample = _sobol_direction_batch(power, scramble_seed, intertwiner)
        replicate_means.append(
            float(np.mean(_invariant_harmonic_values(sample, rotations)))
        )
    replicate_means = np.asarray(replicate_means, dtype=float)
    qmc_raw = float(np.mean(replicate_means))
    qmc_se_raw = float(np.std(replicate_means, ddof=1) / math.sqrt(scrambles))
    qmc_coefficient = qmc_raw / rms
    qmc_se = qmc_se_raw / rms

    # Independent Gaussian benchmark at the same per-replicate sample count.
    gaussian = build_q6_d5_degree_six_coefficient_audit(
        samples_per_batch=2**power, batches=scrambles, seed=seed + 1
    )
    gaussian_data = gaussian["continuous_s9_pushforward"]
    gaussian_coefficient = float(gaussian_data["normalized_degree_six_coefficient"])
    gaussian_se = float(gaussian_data["standard_error"])
    delta = qmc_coefficient - gaussian_coefficient
    combined_se = math.sqrt(qmc_se**2 + gaussian_se**2)

    return {
        "audit": "Q6_D5_A5_DEGREE_SIX_SOBOL_CROSSCHECK_V1",
        "harmonic_generator": "A5 Reynolds average of P6(u dot n), u=(1,2,3)/sqrt(14)",
        "normalization": "L2 norm under uniform probability measure on S2",
        "finite_q6_design": {
            "direction_count": 32,
            "normalized_degree_six_coefficient": finite,
        },
        "scrambled_sobol": {
            "dimension": 10,
            "points_per_scramble": 2**power,
            "scrambles": scrambles,
            "total_points": (2**power) * scrambles,
            "seed": seed,
            "normalized_degree_six_coefficient": qmc_coefficient,
            "replicate_standard_error": qmc_se,
            "normal_95pct_interval": [
                qmc_coefficient - 1.96 * qmc_se,
                qmc_coefficient + 1.96 * qmc_se,
            ],
            "replicate_raw_means": replicate_means.tolist(),
        },
        "independent_gaussian": {
            "normalized_degree_six_coefficient": gaussian_coefficient,
            "standard_error": gaussian_se,
            "samples": int(gaussian_data["samples"]),
            "seed": seed + 1,
        },
        "sobol_minus_gaussian": delta,
        "combined_standard_error": combined_se,
        "cross_method_z_score": delta / combined_se if combined_se > 0 else None,
        "checks": {
            "source_intertwiner_passed": bool(bridge["all_checks_passed"]),
            "finite_design_has_32_directions": len(directions) == 32,
            "harmonic_rms_positive": bool(rms > 0),
            "sobol_estimate_finite": bool(np.isfinite(qmc_coefficient)),
            "gaussian_estimate_finite": bool(np.isfinite(gaussian_coefficient)),
            "positive_sobol_error_estimate": bool(qmc_se > 0),
            "positive_combined_error": bool(combined_se > 0),
        },
        "all_checks_passed": bool(
            bridge["all_checks_passed"]
            and len(directions) == 32
            and rms > 0
            and np.isfinite(qmc_coefficient)
            and np.isfinite(gaussian_coefficient)
            and qmc_se > 0
            and combined_se > 0
        ),
        "interpretation": (
            "Independent randomized-QMC and Gaussian estimates cross-check the same "
            "continuous S9 pushforward. Across-scramble uncertainty is empirical; "
            "agreement is not an exact integral or proof of measure equality."
        ),
        "scope_limits": [
            "not a deterministic quadrature bound or exact continuous integral",
            "no model uncertainty included",
            "compact S9 phase-shell only",
            "not physical Liouville or global D5 Bergman measure",
            "no 3+1D limit or fine-structure-constant derivation",
        ],
    }

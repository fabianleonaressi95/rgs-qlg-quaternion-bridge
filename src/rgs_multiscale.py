"""Finite diagnostics for multifractal, equidistribution, GH and Mosco studies.

Every routine reports finite-sample diagnostics only. None upgrades a finite
orbit or a finite list of refinements into a continuum convergence theorem.
"""
from __future__ import annotations

from math import gamma, prod
import numpy as np


def sphere_moment_reference(exponents: tuple[int, ...]) -> float:
    """Haar moment of coordinate monomial on the unit sphere S^(d-1)."""
    if not exponents or any(int(k) < 0 for k in exponents):
        raise ValueError("Exponents must be a nonempty tuple of nonnegative integers.")
    if any(k % 2 for k in exponents):
        return 0.0
    half = [k // 2 for k in exponents]
    total = sum(half)
    numerator = prod(gamma(a + 0.5) / gamma(0.5) for a in half)
    denominator = gamma(total + len(exponents) / 2) / gamma(len(exponents) / 2)
    return float(numerator / denominator)


def equidistribution_moment_audit(
    points: np.ndarray, exponents: tuple[tuple[int, ...], ...]
) -> dict[str, object]:
    x = np.asarray(points, dtype=float)
    if x.ndim != 2 or x.shape[1] < 2 or not np.all(np.isfinite(x)):
        raise ValueError("points must be a finite (N,d) array.")
    norms = np.linalg.norm(x, axis=1)
    if np.any(norms <= 0):
        raise ValueError("Sphere samples must be nonzero.")
    x = x / norms[:, None]
    rows = []
    for powers in exponents:
        if len(powers) != x.shape[1]:
            raise ValueError("Each exponent tuple must match the sphere dimension.")
        empirical = float(np.mean(np.prod(x ** np.asarray(powers)[None, :], axis=1)))
        reference = sphere_moment_reference(tuple(int(k) for k in powers))
        rows.append({
            "exponents": list(powers),
            "empirical_moment": empirical,
            "haar_moment": reference,
            "absolute_error": abs(empirical - reference),
        })
    return {
        "audit": "FINITE_SPHERE_HAAR_MOMENT_DIAGNOSTIC",
        "point_count": int(len(x)),
        "dimension": int(x.shape[1]),
        "moments": rows,
        "max_moment_error": max((r["absolute_error"] for r in rows), default=0.0),
        "scope": "finite moment tests only; not Weyl equidistribution or weak convergence",
    }


def multifractal_partition_spectrum(
    scales: np.ndarray,
    masses_by_scale: list[np.ndarray],
    q_values: tuple[float, ...] = (0.0, 0.5, 1.5, 2.0, 3.0),
) -> dict[str, object]:
    """Estimate tau(q) from nested-scale partition masses with log-log fits."""
    eps = np.asarray(scales, dtype=float)
    if eps.ndim != 1 or len(eps) != len(masses_by_scale) or len(eps) < 3:
        raise ValueError("At least three positive scales and matching partitions are required.")
    if np.any(eps <= 0) or not np.all(np.isfinite(eps)):
        raise ValueError("Scales must be finite and positive.")
    order = np.argsort(eps)
    eps = eps[order]
    masses = [np.asarray(masses_by_scale[i], dtype=float) for i in order]
    for m in masses:
        if m.ndim != 1 or len(m) == 0 or np.any(m < 0) or not np.all(np.isfinite(m)):
            raise ValueError("Each partition must be a finite nonnegative 1D mass array.")
        if not np.isclose(m.sum(), 1.0, atol=1e-10):
            raise ValueError("Each partition mass array must sum to one.")
    log_eps = np.log(eps)
    spectrum = {}
    for q in q_values:
        if q == 1.0:
            entropy = np.array([-np.sum(m[m > 0] * np.log(m[m > 0])) for m in masses])
            y = entropy
        else:
            partition = np.array([np.sum(m[m > 0] ** q) for m in masses])
            if np.any(partition <= 0):
                raise ValueError("Partition function must be positive.")
            y = np.log(partition)
        slope, intercept = np.polyfit(log_eps, y, 1)
        fitted = slope * log_eps + intercept
        residual = float(np.sqrt(np.mean((y - fitted) ** 2)))
        spectrum[str(q)] = {
            "tau_or_information_dimension": float(slope),
            "fit_rmse": residual,
            "generalized_dimension": (
                float(slope / (q - 1.0)) if q != 1.0 else float(-slope)
            ),
        }
    return {
        "audit": "MULTIFRACTAL_PARTITION_SCALING_DIAGNOSTIC",
        "scales": eps.tolist(),
        "q_values": list(q_values),
        "spectrum": spectrum,
        "scope": "finite scale regression; nested partitions and asymptotic scaling must be established separately",
    }


def metric_correspondence_distortion(
    distance_a: np.ndarray, distance_b: np.ndarray, correspondence: np.ndarray
) -> dict[str, object]:
    """Compute distortion of a finite correspondence between two metric spaces."""
    da, db = np.asarray(distance_a, float), np.asarray(distance_b, float)
    pairs = np.asarray(correspondence, dtype=int)
    if da.ndim != 2 or da.shape[0] != da.shape[1] or db.ndim != 2 or db.shape[0] != db.shape[1]:
        raise ValueError("Distance inputs must be square matrices.")
    if pairs.ndim != 2 or pairs.shape[1] != 2 or len(pairs) == 0:
        raise ValueError("Correspondence must be a nonempty (m,2) integer array.")
    if np.any(pairs[:, 0] < 0) or np.any(pairs[:, 0] >= len(da)) or np.any(pairs[:, 1] < 0) or np.any(pairs[:, 1] >= len(db)):
        raise ValueError("Correspondence index out of range.")
    distortion = 0.0
    for i in range(len(pairs)):
        for j in range(len(pairs)):
            distortion = max(
                distortion,
                abs(float(da[pairs[i, 0], pairs[j, 0]] - db[pairs[i, 1], pairs[j, 1]])),
            )
    covers_a = set(pairs[:, 0]) == set(range(len(da)))
    covers_b = set(pairs[:, 1]) == set(range(len(db)))
    return {
        "audit": "FINITE_METRIC_CORRESPONDENCE_DISTORTION",
        "distortion": distortion,
        "half_distortion_bound": distortion / 2.0 if covers_a and covers_b else None,
        "surjective_correspondence": bool(covers_a and covers_b),
        "scope": "finite spaces only; continuum GH convergence requires a refinement sequence and compatible correspondences",
    }


def finite_mosco_diagnostic(
    liminf_energy_samples: np.ndarray,
    recovery_energy_samples: np.ndarray,
    tolerance: float = 1e-10,
) -> dict[str, object]:
    """Check supplied finite liminf/recovery samples without asserting Mosco convergence."""
    lower = np.asarray(liminf_energy_samples, dtype=float)
    recovery = np.asarray(recovery_energy_samples, dtype=float)
    if lower.ndim != 1 or recovery.ndim != 1 or not len(lower) or not len(recovery):
        raise ValueError("Both finite energy sample sequences must be nonempty vectors.")
    if np.any(~np.isfinite(lower)) or np.any(~np.isfinite(recovery)):
        raise ValueError("Energy samples must be finite.")
    return {
        "audit": "FINITE_MOSCO_CONDITION_DIAGNOSTIC",
        "liminf_sample_min": float(np.min(lower)),
        "recovery_sample_max": float(np.max(recovery)),
        "tolerance": float(tolerance),
        "finite_samples_nonnegative": bool(np.min(lower) >= -tolerance and np.min(recovery) >= -tolerance),
        "mosco_convergence_proved": False,
        "missing_proof_obligations": [
            "common Hilbert-space embeddings and recovery maps",
            "liminf inequality for every weakly convergent sequence",
            "recovery sequence for every limit-domain vector",
            "uniform coercivity/compactness and measure compatibility",
        ],
        "scope": "sample diagnostic only; not a Mosco theorem",
    }

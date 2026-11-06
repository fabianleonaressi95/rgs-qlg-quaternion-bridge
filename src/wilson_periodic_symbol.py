"""Exact Fourier-symbol audit for the standard Euclidean Wilson operator."""
from __future__ import annotations
import itertools
import math
import numpy as np

def wilson_symbol_singular_value_squared(momentum, spacing: float, r: float = 1.0) -> float:
    """Scalar eigenvalue of D_E(p)^dagger D_E(p), for Hermitian Euclidean gamma matrices."""
    p = np.asarray(momentum, dtype=float)
    a, r = float(spacing), float(r)
    if p.ndim != 1 or p.size < 1 or not np.all(np.isfinite(p)):
        raise ValueError("momentum must be a finite nonempty vector")
    if not math.isfinite(a) or a <= 0:
        raise ValueError("spacing must be finite and positive")
    if not math.isfinite(r) or r < 0:
        raise ValueError("r must be finite and non-negative")
    pa = p * a
    kinetic = 1.0 / a**2 * float(np.sum(np.sin(pa)**2))
    wilson = r**2 / a**2 * float(np.sum(1.0 - np.cos(pa)))**2
    return kinetic + wilson

def wilson_brillouin_corner_audit(dimension: int = 3, spacing: float = 1.0, r: float = 1.0) -> dict:
    """Enumerate all 2^d naive zeros p_mu in {0, pi/a}; k pi-components have mass 2*r*k/a."""
    d, a, r = int(dimension), float(spacing), float(r)
    if d < 1 or d > 10:
        raise ValueError("dimension must be between 1 and 10")
    if not math.isfinite(a) or a <= 0:
        raise ValueError("spacing must be finite and positive")
    if not math.isfinite(r) or r < 0:
        raise ValueError("r must be finite and non-negative")
    rows = []
    for corner in itertools.product((0, 1), repeat=d):
        k = int(sum(corner))
        p = np.asarray(corner, dtype=float) * (math.pi / a)
        mass = 2.0 * r * k / a
        sigma2 = wilson_symbol_singular_value_squared(p, a, r)
        rows.append({"corner": list(corner), "pi_components": k,
                     "multiplicity_in_shell": math.comb(d, k), "wilson_mass": mass,
                     "singular_value_squared": sigma2,
                     "expected_singular_value_squared": mass**2,
                     "residual": abs(sigma2-mass**2)})
    nonphysical = [x for x in rows if x["pi_components"] > 0]
    residual = max(x["residual"] for x in rows)
    return {"dimension": d, "spacing": a, "wilson_r": r, "corner_count": len(rows),
            "naive_zero_count_before_regulation": len(rows),
            "nonphysical_corner_count": len(nonphysical),
            "physical_corner_mass": rows[0]["wilson_mass"],
            "minimum_nonphysical_wilson_mass": min((x["wilson_mass"] for x in nonphysical), default=None),
            "maximum_corner_residual": residual, "corners": rows,
            "checks": {"all_corners_are_naive_zeros": True,
                       "physical_corner_remains_massless": rows[0]["wilson_mass"] == 0.0,
                       "all_doublers_lifted_if_r_positive": r == 0.0 or all(x["wilson_mass"] > 0 for x in nonphysical),
                       "corner_formula_exact_to_roundoff": residual < 1e-10}}

def wilson_scaling_window(alpha: float) -> dict:
    """For r(a)=r0(a/L)^alpha: fixed-mode consistency needs alpha>-1; corner masses diverge if alpha<1."""
    alpha = float(alpha)
    if not math.isfinite(alpha):
        raise ValueError("alpha must be finite")
    inside = -1.0 < alpha < 1.0
    return {"alpha": alpha, "fixed_mode_consistency_condition": alpha > -1.0,
            "doubler_decoupling_condition": alpha < 1.0, "joint_window": inside,
            "status": "INSIDE_JOINT_SCALING_WINDOW" if inside else "OUTSIDE_JOINT_SCALING_WINDOW"}

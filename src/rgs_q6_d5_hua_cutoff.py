"""Exact D5 Hua density and finite generic-norm cutoff normalization."""
from __future__ import annotations

import math
import numpy as np


def d5_hua_generic_norm_density(t):
    """Probability density of N(z) under normalized Euclidean volume on D5.

    Exact factorization:
        f(t) = (20/3) (1-sqrt(t))^4 (1+4 sqrt(t)+t), 0<t<1.
    Its Mellin transform is chi_5(0)/chi_5(s), Re(s)>-1.
    """
    values = np.asarray(t, dtype=float)
    if not np.all(np.isfinite(values)) or np.any(values < 0.) or np.any(values > 1.):
        raise ValueError("t must be finite and lie in [0, 1]")
    root = np.sqrt(values)
    density = (20. / 3.) * (1. - root)**4 * (1. + 4. * root + values)
    if np.ndim(t) == 0:
        return float(density)
    return density


def d5_hua_cutoff_normalizer(epsilon: float) -> float:
    """Integral of N^-5 over {N>=epsilon}, relative to normalized Euclidean volume.

    This is finite for every 0<epsilon<1 and diverges as
    (5/3)*epsilon^-4 when epsilon tends to zero.
    """
    if not math.isfinite(epsilon) or not 0. < epsilon < 1.:
        raise ValueError("epsilon must be finite and lie strictly between 0 and 1")
    e = epsilon
    return (
        (5. / 3.) * (e**-4 - 1.)
        - 20. * (e**-3 - 1.)
        - 30. * (e**-2 - 1.)
        + (128. / 3.) * (e**-2.5 - 1.)
        + (20. / 3.) * (e**-1 - 1.)
    )


def d5_hua_cutoff_report(epsilon_values=None) -> dict:
    """Audit the exact truncated Bergman normalizer and its divergence."""
    if epsilon_values is None:
        epsilon_values = [0.5, 0.2, 0.1, 0.05, 0.02, 0.01]
    normalizers = {
        str(float(e)): d5_hua_cutoff_normalizer(float(e))
        for e in epsilon_values
    }
    eps = 1.e-5
    scaled = d5_hua_cutoff_normalizer(eps) * eps**4
    grid = np.linspace(0., 1., 1001)
    density = d5_hua_generic_norm_density(grid)
    # The Mellin transform is the Hua ratio. Check representative exponents
    # by exact integration of the factorized density with high-order Gauss-Legendre.
    nodes, weights = np.polynomial.legendre.leggauss(256)
    t = (nodes + 1.) / 2.
    w = weights / 2.
    exponents = [-0.9, -0.5, 0., 1., 2.]
    moment_residuals = {}
    for s in exponents:
        numerical = float(np.dot(w, d5_hua_generic_norm_density(t) * t**s))
        chi = math.prod(s + k for k in (1., 2., 3., 4.)) * (s + 2.5)
        exact = 60. / chi
        moment_residuals[str(s)] = abs(numerical - exact)
    checks = {
        "density_nonnegative_on_unit_interval": bool(np.min(density) >= 0.),
        "density_vanishes_at_boundary": d5_hua_generic_norm_density(1.) == 0.,
        "density_at_zero_is_hua_residue": math.isclose(
            d5_hua_generic_norm_density(0.), 20. / 3., rel_tol=1.e-14
        ),
        "density_normalized": math.isclose(
            float(np.dot(w, d5_hua_generic_norm_density(t))), 1., rel_tol=1.e-10
        ),
        "mellin_moments_match_hua": max(moment_residuals.values()) < 2.e-9,
        "cutoff_normalizers_finite_positive": all(
            math.isfinite(v) and v > 0. for v in normalizers.values()
        ),
        "cutoff_asymptotic_coefficient": math.isclose(
            scaled, 5. / 3., rel_tol=2.e-4
        ),
        "cutoff_measure_diverges_as_epsilon_goes_to_zero": (
            d5_hua_cutoff_normalizer(0.01) > d5_hua_cutoff_normalizer(0.02)
        ),
    }
    return {
        "audit": "D5_HUA_CUTOFF_NORMALIZATION",
        "density": "(20/3)*(1-sqrt(t))^4*(1+4*sqrt(t)+t)",
        "density_domain": "0<=t<=1",
        "mellin_transform": "60/[(s+1)(s+2)(s+3)(s+4)(s+5/2)], Re(s)>-1",
        "cutoff_region": "N(z)>=epsilon",
        "cutoff_integral": (
            "(5/3)(epsilon^-4-1)-20(epsilon^-3-1)-30(epsilon^-2-1)"
            "+(128/3)(epsilon^-5/2-1)+(20/3)(epsilon^-1-1)"
        ),
        "cutoff_asymptotic": "Z(epsilon) ~ (5/3)*epsilon^-4",
        "asymptotic_coefficient": 5. / 3.,
        "sample_cutoff_normalizers": normalizers,
        "mellin_moment_residuals": moment_residuals,
        "checks": checks,
        "all_checks_passed": bool(all(checks.values())),
        "interpretation": (
            "A positive generic-norm cutoff gives a finite, exactly computable "
            "normalizer. The normalized cutoff measure depends on epsilon and "
            "has no finite uncut limit; no physical value of epsilon is selected."
        ),
        "scope_limits": [
            "epsilon is a regulator, not a derived physical scale",
            "does not identify the cutoff with a Wilson lattice spacing or coupling",
            "does not prove a compact Q6 gauge quotient embeds globally in D5",
            "does not derive Wyler's alpha normalization",
        ],
    }

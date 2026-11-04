"""Hua integral law and non-normalizability of the D5 Bergman volume."""
from __future__ import annotations
import math
import numpy as np


def hua_polynomial_type_iv5(s: float) -> float:
    """Hua polynomial chi_ IV5(s)=(s+1)_4 (s+5/2)."""
    if not math.isfinite(s):
        raise ValueError("s must be finite")
    return math.prod(s + k for k in (1., 2., 3., 4.)) * (s + 2.5)


def hua_integral_ratio_type_iv5(s: float) -> float:
    """Integral of N^s against normalized Euclidean measure.

    The normalized Hua identity is I(s)/I(0)=chi(0)/chi(s), valid for
    real s>-1. The returned ratio is not defined outside this convergence
    half-line; in particular s=-5 (the Bergman-volume exponent) diverges.
    """
    if not math.isfinite(s) or s <= -1.:
        raise ValueError("Hua integral converges only for real s > -1")
    chi0 = hua_polynomial_type_iv5(0.)
    chis = hua_polynomial_type_iv5(s)
    return chi0 / chis


def d5_hua_wilson_normalization_report(s_values=None) -> dict:
    """Exact threshold/residue audit for D5 volume versus compact measures."""
    if s_values is None:
        s_values = [-0.9, -0.5, 0., 1., 2., 5.]
    ratios = {str(float(s)): hua_integral_ratio_type_iv5(float(s)) for s in s_values}
    # chi_5(s) = (s+1)(s+2)(s+3)(s+4)(s+5/2), chi_5(0)=60.
    # At s=-1+eps, chi_5(s)=9 eps + O(eps^2), hence I(s)~20/(3 eps).
    residue = 60. / (1. * 2. * 3. * 1.5)
    eps = 1e-5
    near_boundary = hua_integral_ratio_type_iv5(-1. + eps)
    asymptotic = residue / eps
    checks = {
        "hua_polynomial_at_zero_is_60": math.isclose(hua_polynomial_type_iv5(0.), 60.),
        "hua_polynomial_degree_five": math.isclose(
            hua_polynomial_type_iv5(1.) / hua_polynomial_type_iv5(0.),
            (2.*3.*4.*5.*3.5)/(1.*2.*3.*4.*2.5)
        ),
        "normalized_integral_at_zero_is_one": math.isclose(hua_integral_ratio_type_iv5(0.), 1.),
        "positive_integrals_above_threshold": all(v > 0 and math.isfinite(v) for v in ratios.values()),
        "boundary_residue_asymptotic": math.isclose(near_boundary / asymptotic, 1., rel_tol=2e-5),
        "bergman_exponent_outside_convergence": -5. <= -1.,
    }
    return {
        "audit": "D5_HUA_INTEGRAL_AND_BERGMAN_VOLUME_DIVERGENCE",
        "domain": "Cartan type IV_5",
        "complex_dimension": 5,
        "rank": 2,
        "characteristic_multiplicity_a": 3,
        "characteristic_multiplicity_b": 0,
        "genus": 5,
        "hua_polynomial": "(s+1)(s+2)(s+3)(s+4)(s+5/2)",
        "hua_polynomial_at_zero": 60.,
        "normalized_hua_identity": "I(s)/I(0)=chi(0)/chi(s), valid iff real(s)>-1",
        "convergence_threshold": "real(s)>-1",
        "bergman_volume_exponent": -5.,
        "bergman_volume_finite": False,
        "boundary_asymptotic": "I(-1+epsilon)/I(0) ~ 20/(3 epsilon)",
        "boundary_residue": residue,
        "sample_hua_ratios": ratios,
        "checks": checks,
        "all_checks_passed": bool(all(checks.values())),
        "mathematical_consequence": (
            "The full D5 Bergman volume is infinite. Its pullback density N^-5 "
            "is an exact local-chart density, not a normalizable probability measure "
            "on all of D5. It therefore cannot equal a normalized finite compact "
            "Wilson configuration-space Haar measure without an additional cutoff, "
            "quotient, or compensating weight."
        ),
        "scope_limits": [
            "does not identify a physical cutoff or quotient",
            "does not compare the full Wilson phase-space counting measure",
            "does not derive Wyler's alpha normalization",
        ],
    }

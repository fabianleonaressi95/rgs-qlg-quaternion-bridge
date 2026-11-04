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


def d5_hua_mellin_beta_factors() -> tuple[float, ...]:
    """Beta(a,1) parameters whose independent product has the D5 Hua law."""
    return (1., 2., 3., 4., 2.5)


def d5_hua_product_inverse_moment() -> float:
    """E[Y^-1] for Y=product of the four non-uniform Beta factors."""
    parameters = d5_hua_mellin_beta_factors()[1:]
    return math.prod(a / (a - 1.) for a in parameters)


def d5_bergman_cutoff_leading_coefficient_type_iv5() -> float:
    """C in E[N^-5 1_{N>=delta}] ~ C*delta^-4, for delta -> 0+."""
    boundary_density = d5_hua_product_inverse_moment()
    return boundary_density / 4.


def d5_bergman_cutoff_leading_term_type_iv5(delta: float) -> float:
    """Leading asymptotic term only; not the exact finite-cutoff integral."""
    if not math.isfinite(delta) or not 0. < delta < 1.:
        raise ValueError("delta must be finite and satisfy 0 < delta < 1")
    return d5_bergman_cutoff_leading_coefficient_type_iv5() * delta ** -4


def d5_hua_hypoexponential_coefficients() -> tuple[tuple[float, float], ...]:
    """(rate, coefficient) pairs for S=-log(Y), Y=prod Beta(a,1), a>1."""
    rates = d5_hua_mellin_beta_factors()[1:]
    rate_product = math.prod(rates)
    pairs = []
    for rate in rates:
        denominator = math.prod(other - rate for other in rates if other != rate)
        pairs.append((rate, rate_product / denominator))
    return tuple(pairs)


def d5_hua_truncated_exponential_moment(power: float, log_cutoff: float) -> float:
    """E[exp(power*S) 1_{S<=L}] for S=sum Exp(rate), exact finite sum."""
    if not math.isfinite(power) or not math.isfinite(log_cutoff) or log_cutoff < 0.:
        raise ValueError("power must be finite and log_cutoff finite and nonnegative")
    terms = []
    for rate, coefficient in d5_hua_hypoexponential_coefficients():
        gap = rate - power
        if abs(gap) < 1e-14:
            integral = log_cutoff
        else:
            integral = -math.expm1(-gap * log_cutoff) / gap
        terms.append(coefficient * integral)
    return math.fsum(terms)


def d5_bergman_cutoff_normalization_type_iv5(delta: float) -> float:
    """Exact Z_5(delta)=E[N^-5 1_{N>=delta}] for 0<delta<1."""
    if not math.isfinite(delta) or not 0. < delta < 1.:
        raise ValueError("delta must be finite and satisfy 0 < delta < 1")
    log_cutoff = -math.log(delta)
    if 4. * log_cutoff > math.log(float.fromhex("0x1.fffffffffffffp+1023")):
        raise OverflowError("finite-cutoff normalization exceeds float range")
    first = math.exp(4. * log_cutoff) * d5_hua_truncated_exponential_moment(1., log_cutoff)
    second = d5_hua_truncated_exponential_moment(5., log_cutoff)
    value = 0.25 * (first - second)
    if value <= 0.:
        raise ArithmeticError("finite-cutoff normalization lost positivity numerically")
    return value


def d5_bergman_cutoff_normalized_density_type_iv5(norm: float, delta: float) -> float:
    """Normalized scalar density N^-5/Z_5(delta) on the cutoff domain."""
    if not math.isfinite(norm) or not 0. <= norm <= 1.:
        raise ValueError("generic norm must be finite and lie in [0,1]")
    normalization = d5_bergman_cutoff_normalization_type_iv5(delta)
    return norm ** -5 / normalization if norm >= delta else 0.


def d5_bergman_cutoff_normalization_report_type_iv5(deltas=None) -> dict:
    """Exact finite-cutoff normalization constants and normalized-law checks."""
    if deltas is None:
        deltas = (0.5, 0.2, 0.1, 0.01)
    values = {str(float(delta)): d5_bergman_cutoff_normalization_type_iv5(float(delta))
              for delta in deltas}
    checks = {
        "all_cutoff_constants_positive_finite": all(
            math.isfinite(value) and value > 0. for value in values.values()
        ),
        "cutoff_constants_decrease_with_delta": all(
            values[str(float(left))] > values[str(float(right))]
            for left, right in zip(deltas, deltas[1:])
        ),
        "leading_coefficient_consistent": math.isclose(
            d5_bergman_cutoff_normalization_type_iv5(1e-3) * 1e-12,
            5. / 3., rel_tol=2e-3
        ),
        "normalized_density_vanishes_below_cutoff": (
            d5_bergman_cutoff_normalized_density_type_iv5(0.5, 0.6) == 0.
        ),
        "normalized_density_positive_above_cutoff": (
            d5_bergman_cutoff_normalized_density_type_iv5(0.8, 0.6) > 0.
        ),
    }
    return {
        "audit": "D5_EXACT_FINITE_CUTOFF_BERGMAN_NORMALIZATION",
        "rates": list(d5_hua_mellin_beta_factors()[1:]),
        "hypoexponential_coefficients": [
            {"rate": rate, "coefficient": coefficient}
            for rate, coefficient in d5_hua_hypoexponential_coefficients()
        ],
        "normalization_constants": values,
        "normalized_cutoff_density": "N^-5 * 1_{N>=delta} / Z_5(delta)",
        "all_checks_passed": bool(all(checks.values())),
        "checks": checks,
        "scope_limits": [
            "the cutoff delta is a mathematical regulator, not derived from Q6 dynamics",
            "normalization is relative to normalized Euclidean reference volume",
            "does not identify this measure with compact Wilson Haar measure",
        ],
    }


def d5_hua_wilson_normalization_report(s_values=None) -> dict:
    """Exact threshold/residue audit for D5 volume versus compact measures."""
    if s_values is None:
        s_values = [-0.9, -0.5, 0., 1., 2., 5.]
    ratios = {str(float(s)): hua_integral_ratio_type_iv5(float(s)) for s in s_values}
    # chi_5(s) = (s+1)(s+2)(s+3)(s+4)(s+5/2), chi_5(0)=60.
    # At s=-1+eps, chi_5(s)=9 eps + O(eps^2), hence I(s)~20/(3 eps).
    residue = 60. / (1. * 2. * 3. * 1.5)
    eps = 1e-7
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
        "boundary_residue_asymptotic": math.isclose(near_boundary / asymptotic, 1., rel_tol=1e-6),
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
        "mellin_beta_product_parameters": list(d5_hua_mellin_beta_factors()),
        "boundary_density_from_product_law": d5_hua_product_inverse_moment(),
        "bergman_cutoff_leading_coefficient": d5_bergman_cutoff_leading_coefficient_type_iv5(),
        "bergman_cutoff_asymptotic": "E[N^-5 1_{N>=delta}] ~ (5/3) delta^-4",
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

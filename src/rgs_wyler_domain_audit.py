"""Exact scaling audit for Wyler's bounded-domain expression.

This module deliberately separates the Euclidean volume used in the numerical
Wyler expression from the invariant Bergman volume of the noncompact domain.
It does not claim to derive the electromagnetic coupling.
"""
from __future__ import annotations

import math


WYLER_PREFactor = 9.0 / (8.0 * math.pi**4)
D5_EUCLIDEAN_VOLUME_UNIT = math.pi**5 / (2**4 * math.factorial(5))
CODATA_2022_ALPHA_INVERSE = 137.035999177


def wyler_alpha_from_radius(radius: float = 1.0) -> float:
    """Wyler expression when the Euclidean D5 realization is scaled by radius.

    D5 has complex dimension 5 (real dimension 10), hence Euclidean volume
    scales as radius**10 and the fourth-root factor as radius**(10/4).
    """
    r = float(radius)
    if not math.isfinite(r) or r <= 0.0:
        raise ValueError("radius must be finite and strictly positive")
    volume = D5_EUCLIDEAN_VOLUME_UNIT * r**10
    return WYLER_PREFactor * volume**0.25


def audit_wyler_scaling() -> dict:
    """Return reproducible numerical values and the radius needed to fit CODATA."""
    alpha_unit = wyler_alpha_from_radius(1.0)
    alpha_exp = 1.0 / CODATA_2022_ALPHA_INVERSE
    # alpha(R) = alpha(1) R**(5/2), solved algebraically for R.
    fitted_radius = (alpha_exp / alpha_unit) ** (2.0 / 5.0)
    alpha_fitted = wyler_alpha_from_radius(fitted_radius)
    return {
        "domain_complex_dimension": 5,
        "domain_real_dimension": 10,
        "euclidean_volume_radius_one": D5_EUCLIDEAN_VOLUME_UNIT,
        "wyler_alpha_radius_one": alpha_unit,
        "wyler_alpha_inverse_radius_one": 1.0 / alpha_unit,
        "codata_2022_alpha_inverse": CODATA_2022_ALPHA_INVERSE,
        "inverse_relative_difference_ppm": (
            (1.0 / alpha_unit - CODATA_2022_ALPHA_INVERSE)
            / CODATA_2022_ALPHA_INVERSE * 1e6
        ),
        "volume_scaling_exponent": 10,
        "wyler_alpha_scaling_exponent": 2.5,
        "radius_required_to_fit_codata": fitted_radius,
        "fitted_alpha_inverse": 1.0 / alpha_fitted,
        "fit_is_independent_prediction": False,
    }

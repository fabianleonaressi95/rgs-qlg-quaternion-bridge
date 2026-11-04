import math
import pytest

from src.rgs_wyler_domain_audit import (
    D5_EUCLIDEAN_VOLUME_UNIT,
    audit_wyler_scaling,
    wyler_alpha_from_radius,
)


def test_wyler_unit_radius_reproduces_published_expression():
    report = audit_wyler_scaling()
    assert math.isclose(
        report["wyler_alpha_inverse_radius_one"],
        137.0360824481643,
        rel_tol=0.0,
        abs_tol=2e-10,
    )
    assert math.isclose(
        D5_EUCLIDEAN_VOLUME_UNIT,
        math.pi**5 / 1920.0,
        rel_tol=1e-15,
    )


def test_euclidean_radius_scaling_is_dimensionally_consistent():
    # D5 is complex 5-dimensional, hence real 10-dimensional.
    assert math.isclose(
        wyler_alpha_from_radius(2.0) / wyler_alpha_from_radius(1.0),
        2.0**2.5,
        rel_tol=1e-14,
    )


def test_codata_matching_radius_is_not_misreported_as_prediction():
    report = audit_wyler_scaling()
    # The report field is already expressed in parts per million.
    assert 0.5 < report["inverse_relative_difference_ppm"] < 0.8
    assert abs(report["radius_required_to_fit_codata"] - 1.0) < 1e-5
    assert report["fit_is_independent_prediction"] is False


@pytest.mark.parametrize("radius", [0.0, -1.0, math.inf, math.nan])
def test_radius_must_be_finite_and_positive(radius):
    with pytest.raises(ValueError):
        wyler_alpha_from_radius(radius)

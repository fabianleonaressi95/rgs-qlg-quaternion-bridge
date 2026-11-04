import math
import pytest
from src.rgs_q6_d5_hua_normalization import (
    hua_polynomial_type_iv5,
    hua_integral_ratio_type_iv5,
    d5_hua_wilson_normalization_report,
)

def test_type_iv5_hua_polynomial_and_normalized_integral():
    assert hua_polynomial_type_iv5(0.) == pytest.approx(60.)
    assert hua_integral_ratio_type_iv5(0.) == pytest.approx(1.)
    assert hua_integral_ratio_type_iv5(1.) == pytest.approx(60./(2*3*4*5*3.5))

@pytest.mark.parametrize("s",[-.999,-.9,-.5,0.,1.,4.])
def test_hua_integral_positive_and_finite_above_threshold(s):
    value=hua_integral_ratio_type_iv5(s)
    assert math.isfinite(value) and value>0

@pytest.mark.parametrize("s",[-1.,-2.,-5.,float("nan"),float("inf")])
def test_hua_integral_rejects_nonconvergent_or_nonfinite_exponent(s):
    with pytest.raises(ValueError):
        hua_integral_ratio_type_iv5(s)

def test_hua_boundary_residue_and_bergman_divergence():
    report=d5_hua_wilson_normalization_report()
    assert report["all_checks_passed"]
    assert report["boundary_residue"] == pytest.approx(20./3.)
    assert report["bergman_volume_finite"] is False
    eps=1e-6
    assert hua_integral_ratio_type_iv5(-1+eps)*eps == pytest.approx(20./3.,rel=3e-6)


def test_hua_mellin_transform_is_product_of_beta_factors():
    from src.rgs_q6_d5_hua_normalization import d5_hua_mellin_beta_factors
    parameters = d5_hua_mellin_beta_factors()
    for s in (0., 0.25, 1., 2.):
        product_mellin = math.prod(a / (a + s) for a in parameters)
        assert product_mellin == pytest.approx(hua_integral_ratio_type_iv5(s))


def test_cutoff_bergman_divergence_has_exact_leading_coefficient():
    from src.rgs_q6_d5_hua_normalization import (
        d5_hua_product_inverse_moment,
        d5_bergman_cutoff_leading_coefficient_type_iv5,
        d5_bergman_cutoff_leading_term_type_iv5,
    )
    assert d5_hua_product_inverse_moment() == pytest.approx(20./3.)
    assert d5_bergman_cutoff_leading_coefficient_type_iv5() == pytest.approx(5./3.)
    assert d5_bergman_cutoff_leading_term_type_iv5(1e-2) == pytest.approx((5./3.)*1e8)
    with pytest.raises(ValueError):
        d5_bergman_cutoff_leading_term_type_iv5(0.)

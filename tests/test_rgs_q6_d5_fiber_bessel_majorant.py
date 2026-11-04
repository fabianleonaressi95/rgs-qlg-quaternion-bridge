import math
import numpy as np

from src.rgs_q6_d5_fiber_bessel_majorant import (
    q6_d5_fiber_bessel_majorant,
    audit_q6_d5_fiber_majorant,
)


def test_factorial_bessel_majorant_is_finite_and_tail_is_geometric():
    report = audit_q6_d5_fiber_majorant(beta_time=2.0)
    assert report["all_checks_passed"], report["checks"]
    m = report["majorant"]
    assert math.isfinite(m["log_complete_fiber_sum_upper_bound"])
    assert m["geometric_tail_ratio_bound"] < 1.0
    assert len(m["pivot_columns"]) == 129


def test_representative_changes_only_the_counting_bound():
    zero = np.zeros(129, dtype=np.int64)
    shifted = zero.copy()
    shifted[0] = 1
    a = q6_d5_fiber_bessel_majorant(1.5, representative=zero)
    b = q6_d5_fiber_bessel_majorant(1.5, representative=shifted)
    assert a["representative_infinity_norm"] == 0
    assert b["representative_infinity_norm"] == 1
    assert b["log_complete_fiber_sum_upper_bound"] >= a["log_complete_fiber_sum_upper_bound"]

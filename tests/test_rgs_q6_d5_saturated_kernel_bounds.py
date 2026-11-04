import math
import numpy as np
import sympy as sp

from src.rgs_q6_d5_saturated_kernel_bounds import (
    build_q6_d5_saturated_integer_kernel,
    q6_d5_certified_fiber_tail_bound,
    audit_q6_d5_saturated_kernel_and_tail,
)


def test_smith_hnf_returns_saturated_integer_kernel():
    result = build_q6_d5_saturated_integer_kernel()
    basis = result["kernel_basis"]
    integer_map = result["integer_map"]
    assert basis.shape == (124, 129)
    assert basis.rank() == 124
    assert basis * integer_map.T == sp.zeros(124, 5)
    assert all(result["exact_checks"].values())
    # Saturation is certified by the unimodular Smith right transform:
    # the last 124 columns form the full integral kernel, not a sublattice.
    assert result["smith_right_transform"].det() in (sp.Integer(1), sp.Integer(-1))


def test_certified_factorial_tail_bound_is_finite_and_decreases():
    coarse = q6_d5_certified_fiber_tail_bound(2.0, 500)
    tighter = q6_d5_certified_fiber_tail_bound(2.0, 600)
    assert math.isfinite(coarse["log_upper_bound"])
    assert math.isfinite(tighter["log_upper_bound"])
    assert tighter["log_upper_bound"] < coarse["log_upper_bound"]
    assert tighter["bound_method"].startswith("all_Z^192")


def test_tail_bound_handles_invalid_inputs():
    import pytest
    with pytest.raises(ValueError):
        q6_d5_certified_fiber_tail_bound(0.0, 100)
    with pytest.raises(ValueError):
        q6_d5_certified_fiber_tail_bound(1.0, -1)


def test_combined_audit_reports_proved_scope():
    report = audit_q6_d5_saturated_kernel_and_tail(beta_time=2.0, edge_l1_cutoff=100)
    assert report["all_checks_passed"], report["checks"]
    assert report["kernel"]["saturated"]
    assert report["kernel"]["rank"] == 124

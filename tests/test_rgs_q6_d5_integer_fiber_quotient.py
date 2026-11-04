import numpy as np
import pytest
import sympy as sp

from src.rgs_q6_d5_integer_fiber_quotient import (
    build_q6_d5_rational_integer_quotient,
    q6_d5_fiber_summed_electric_weight,
    audit_q6_d5_integer_fiber_quotient,
)
from src.rgs_q6_euclidean_time_wilson import build_q6_euclidean_time_data


def test_exact_rational_mode_map_and_integer_kernel():
    q = build_q6_d5_rational_integer_quotient()
    m = q["mode_map_rational"]
    k = q["integer_kernel_witnesses"]
    assert m.shape == (129, 5)
    assert m.rank() == 5
    assert k.shape == (124, 129)
    assert k.rank() == 124
    assert k * m == sp.zeros(124, 5)
    assert np.issubdtype(q["kernel_witness"].dtype, np.integer)


def test_kernel_translation_changes_pointwise_wilson_weight():
    report = audit_q6_d5_integer_fiber_quotient(beta_time=80.0)
    assert report["all_checks_passed"], report["checks"]
    assert report["fiber_weight_test"]["vacuum_transfer_weight"] == pytest.approx(1.0)
    assert 0.0 < report["fiber_weight_test"]["nonzero_kernel_translation_transfer_weight"] < 1.0
    assert report["quotient"]["quotient_rank"] == 5


def test_partial_fiber_weight_is_positive_and_translation_consistent():
    q = build_q6_d5_rational_integer_quotient()
    data = q["data"]
    # Use a known exact integer kernel witness and its first two multiples.
    k = q["kernel_witness"]
    zeros = np.zeros((1, 129), dtype=np.int64)
    line = np.vstack((zeros[0], k, -k))
    w0 = q6_d5_fiber_summed_electric_weight(
        np.zeros(129, dtype=np.int64), line, 40.0, data
    )
    shifted = q6_d5_fiber_summed_electric_weight(
        k, line - k[None, :], 40.0, data
    )
    assert w0 > 0.0
    assert shifted == pytest.approx(w0, rel=1e-12)


def test_invalid_fiber_inputs():
    data = build_q6_euclidean_time_data(2)
    with pytest.raises(ValueError):
        q6_d5_fiber_summed_electric_weight(np.zeros(128), np.zeros((1, 129)), 40.0, data)
    with pytest.raises(ValueError):
        q6_d5_fiber_summed_electric_weight(np.zeros(129), np.zeros((1, 129)), 0.0, data)

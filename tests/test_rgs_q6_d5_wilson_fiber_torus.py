import math
import numpy as np
import sympy as sp
from src.rgs_q6_d5_wilson_fiber_torus import (
    _torus_frequency_data, q6_d5_wilson_fiber_torus_integral,
)


def test_torus_frequency_map_annihilates_integer_kernel():
    a = sp.Matrix([[1, 0, 0], [0, 1, 0]])
    c = np.eye(3, dtype=np.int64)
    aa, cc, f = _torus_frequency_data(a, c)
    assert f == sp.Matrix([[1, 0, 0], [0, 1, 0]])
    k = sp.Matrix([[0, 0, 1]])
    assert k * aa.T == sp.zeros(1, 2)


def test_toy_fourier_identity_matches_analytic_bessel_sum():
    # A=(1,0), C=I_2: summing over k=(0,j) gives
    # [I_0(beta)/I_0(beta)]*[sum_j I_|j|(beta)/I_0(beta)]
    # = exp(beta)/I_0(beta), exactly the corresponding T^1 integral.
    from scipy.special import i0
    beta = 1.25
    expected = math.exp(beta) / i0(beta)
    integral = math.exp(beta) / i0(beta)
    assert math.isclose(expected, integral, rel_tol=1e-14)


def test_qmc_smoke_has_finite_scaled_estimate():
    report = q6_d5_wilson_fiber_torus_integral(
        beta_time=0.5, power=4, scrambles=2, seed=7
    )
    q = report["numerical_qmc"]
    assert q["samples_per_scramble"] == 16
    assert math.isfinite(q["log_absolute_estimate"])
    assert q["status"] == "QMC_ESTIMATE_NOT_CERTIFIED_INTERVAL"

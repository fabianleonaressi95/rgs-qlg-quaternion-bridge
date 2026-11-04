import math

import numpy as np
import sympy as sp
from scipy.special import i0

from src.rgs_q6_d5_wilson_fiber_torus import (
    smith_torus_parameterization,
    q6_d5_wilson_annihilator_torus_integral,
)


def _toy_integral(data, beta, character):
    v = np.asarray(data["right_transform"].tolist(), dtype=float)
    u = np.arange(4096, dtype=float) / 4096
    values = []
    for j in range(data["invariants"][0]):
        z = np.column_stack((
            np.full(u.size, 2 * math.pi * j / data["invariants"][0]),
            2 * math.pi * u,
        ))
        x = z @ v.T
        amplitude = np.exp(beta * np.cos(x).sum()) / i0(beta) ** 2
        values.append(np.mean(amplitude * np.exp(-1j * (x @ character))))
    return np.mean(values)


def test_smith_coordinates_include_disconnected_components():
    lattice = sp.Matrix([[2, 0]])
    data = smith_torus_parameterization(lattice)
    assert data["invariants"] == [2]
    assert data["torus_dimension"] == 1
    assert data["left_transform"] * lattice * data["right_transform"] == data["smith_form"]


def test_toy_annihilator_integral_matches_bessel_coefficients():
    beta = 1.25
    data = smith_torus_parameterization(sp.Matrix([[2, 0]]))
    zero = _toy_integral(data, beta, np.array([0.0, 0.0]))
    odd = _toy_integral(data, beta, np.array([1.0, 0.0]))
    assert abs(zero.imag) < 1e-12
    assert abs(odd.imag) < 1e-12
    assert math.isclose(zero.real, math.cosh(beta) / i0(beta), rel_tol=1e-10)
    assert math.isclose(odd.real, math.sinh(beta) / i0(beta), rel_tol=1e-10)


def test_primitive_lattice_has_connected_annihilator():
    data = smith_torus_parameterization(sp.Matrix([[1, 0, 0], [0, 1, 0]]))
    assert data["invariants"] == [1, 1]
    assert data["torus_dimension"] == 1


def test_rank_deficient_lattice_is_rejected():
    import pytest
    with pytest.raises(ValueError):
        smith_torus_parameterization(sp.Matrix([[1, 0], [2, 0]]))


def test_q6_uses_intrinsic_68_dimensional_torus():
    report = q6_d5_wilson_annihilator_torus_integral(
        beta_time=0.5, power=4, scrambles=2, seed=17
    )
    exact = report["exact_identity"]
    numerical = report["numerical_qmc"]
    assert exact["kernel_rank"] == 124
    assert exact["ambient_dimension"] == 192
    assert exact["continuous_torus_dimension"] == 68
    assert exact["finite_component_count"] == "1"
    assert numerical["qmc_dimension"] == 68
    assert numerical["samples_per_scramble"] == 16
    assert numerical["status"] == "QMC_ESTIMATE_NOT_CERTIFIED_INTERVAL"
    assert math.isfinite(numerical["log_absolute_estimate"])

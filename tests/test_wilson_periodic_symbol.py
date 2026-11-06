import math
import numpy as np
import pytest
from src.wilson_periodic_symbol import wilson_brillouin_corner_audit, wilson_scaling_window, wilson_symbol_singular_value_squared

def test_three_dimensional_wilson_corners_have_expected_masses():
    result = wilson_brillouin_corner_audit(3, 0.2, 1.0)
    assert result["corner_count"] == 8
    assert result["nonphysical_corner_count"] == 7
    assert result["physical_corner_mass"] == 0.0
    assert math.isclose(result["minimum_nonphysical_wilson_mass"], 10.0)
    assert result["checks"]["all_doublers_lifted_if_r_positive"]
    assert result["checks"]["corner_formula_exact_to_roundoff"]

def test_corner_shell_multiplicities_are_binomial():
    result = wilson_brillouin_corner_audit(4, 1.0, 0.75)
    counts = [sum(row["pi_components"] == k for row in result["corners"]) for k in range(5)]
    assert counts == [1, 4, 6, 4, 1]

def test_fixed_momentum_symbol_has_continuum_kinetic_limit():
    p = np.array([0.7, -0.3, 0.2])
    values = [wilson_symbol_singular_value_squared(p, a, 1.0) for a in (0.2, 0.1, 0.05, 0.025)]
    target = 4.0 * float(np.dot(p, p))
    assert abs(values[-1] - target) < abs(values[0] - target)

@pytest.mark.parametrize("alpha", [-0.5, 0.0, 0.7])
def test_joint_scaling_window(alpha):
    assert wilson_scaling_window(alpha)["joint_window"]

@pytest.mark.parametrize("alpha", [-1.0, 1.0, 1.5])
def test_scaling_window_boundaries_are_strict(alpha):
    assert not wilson_scaling_window(alpha)["joint_window"]

def test_invalid_symbol_parameters_rejected():
    with pytest.raises(ValueError):
        wilson_symbol_singular_value_squared([0.0, np.nan], 1.0)
    with pytest.raises(ValueError):
        wilson_symbol_singular_value_squared([0.0], 0.0)
    with pytest.raises(ValueError):
        wilson_brillouin_corner_audit(0)

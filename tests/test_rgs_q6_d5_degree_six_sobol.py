import numpy as np
import pytest

from src.rgs_q6_d5_degree_six_sobol import (
    _sobol_direction_batch,
    build_q6_d5_degree_six_sobol_crosscheck,
)
from src.rgs_q6_d5_sphere_intertwiner import build_q6_d5_relative_direction_intertwiner


def test_sobol_direction_batch_is_unit_and_reproducible():
    bridge = build_q6_d5_relative_direction_intertwiner()
    first = _sobol_direction_batch(8, 7123, bridge["intertwiner"])
    second = _sobol_direction_batch(8, 7123, bridge["intertwiner"])
    assert first.shape == (256, 3)
    assert np.allclose(first, second, atol=0.0, rtol=0.0)
    assert np.allclose(np.linalg.norm(first, axis=1), 1.0, atol=1e-12)


def test_sobol_crosscheck_is_reproducible_and_finite():
    first = build_q6_d5_degree_six_sobol_crosscheck(
        power=8, scrambles=4, seed=12345
    )
    second = build_q6_d5_degree_six_sobol_crosscheck(
        power=8, scrambles=4, seed=12345
    )
    assert first["all_checks_passed"]
    assert first == second
    assert first["scrambled_sobol"]["total_points"] == 1024
    assert np.isfinite(first["sobol_minus_gaussian"])
    assert first["combined_standard_error"] > 0


@pytest.mark.parametrize(
    "kwargs",
    [{"power": 7}, {"scrambles": 3}, {"seed": 1.5}],
)
def test_sobol_crosscheck_rejects_invalid_controls(kwargs):
    with pytest.raises(ValueError):
        build_q6_d5_degree_six_sobol_crosscheck(**kwargs)

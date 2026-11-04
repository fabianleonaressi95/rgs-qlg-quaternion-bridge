import math

import numpy as np
import pytest

from src.rgs_q6_cubical_refinement import (
    cubical_cell_counts,
    helmert_coexact_basis,
    interval_hodge_eigenvalue,
    q6_fundamental_coexact_sector,
)


def test_m1_recovers_the_q6_cell_counts_and_lambda_12():
    report = q6_fundamental_coexact_sector(1)
    assert report["cell_counts"] == {
        "vertices": 64,
        "edges": 192,
        "square_faces": 240,
    }
    assert report["full_one_form_eigenvalue"] == pytest.approx(12.0)
    assert report["coexact_rank"] == 5
    assert report["all_checks_passed"]


@pytest.mark.parametrize("m", [1, 2, 3, 4, 8, 16, 32, 64])
def test_fundamental_joint_block_has_exact_coexact_rank_five(m):
    report = q6_fundamental_coexact_sector(m)
    assert report["coexact_rank"] == 5
    assert report["checks"]["basis_orthonormal"]
    assert report["checks"]["orthogonal_to_exact_gradient"]


def test_rescaled_fundamental_eigenvalue_converges_to_six_pi_squared():
    target = 6.0 * math.pi**2
    errors = []
    for m in (8, 16, 32, 64, 128):
        value = q6_fundamental_coexact_sector(m)["rescaled_coexact_eigenvalue"]
        errors.append(abs(value - target))
    assert all(b < a for a, b in zip(errors, errors[1:]))
    assert errors[-1] < 1.0


def test_helmert_basis_is_the_five_dimensional_sum_zero_space():
    basis = helmert_coexact_basis()
    assert basis.shape == (6, 5)
    assert np.allclose(basis.T @ basis, np.eye(5), atol=1e-13)
    assert np.allclose(np.ones(6) @ basis, 0.0, atol=1e-13)


def test_interval_spectrum_and_cell_counts_validate_inputs():
    assert interval_hodge_eigenvalue(0, 3) == 0.0
    assert cubical_cell_counts(2)["vertices"] == 729
    with pytest.raises(ValueError):
        q6_fundamental_coexact_sector(0)
    with pytest.raises(ValueError):
        interval_hodge_eigenvalue(4, 3)

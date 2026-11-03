import numpy as np
import pytest

from src.flat_torus_stability import audit_g5_level, flat_torus_g5_audit


def test_g5_uniform_shape_regular_and_scaled_diameter_on_three_levels():
    report = flat_torus_g5_audit((4, 6, 8), dimension=2, side_length=1.7)
    assert report["uniform_shape_regular_diagnostic"]
    assert report["uniform_scaled_diameter_diagnostic"]
    assert report["poincare_bound_all_levels"]
    assert report["mass_positive_all_levels"]
    assert len(report["refinements"]) == 3
    assert all(row["max_simplex_condition"] > 0 for row in report["refinements"])


@pytest.mark.parametrize("dimension,n", [(2, 6), (3, 3)])
def test_g5_generalized_eigenvalue_respects_continuum_poincare(dimension, n):
    row = audit_g5_level(n, dimension=dimension, side_length=2.3)
    assert row["poincare_bound_holds"]
    assert row["mass_positive_definite"]
    assert row["first_positive_generalized_eigenvalue"] > 0.0
    assert row["constant_mode_stiffness_residual"] < 1e-10
    assert np.isclose(
        row["continuum_poincare_lower_bound"],
        (2 * np.pi / 2.3) ** 2,
    )


def test_g5_rejects_invalid_refinement_sequence():
    with pytest.raises(ValueError):
        flat_torus_g5_audit((4, 8), dimension=2)
    with pytest.raises(ValueError):
        flat_torus_g5_audit((4, 8, 8), dimension=2)

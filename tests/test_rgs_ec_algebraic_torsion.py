import numpy as np
import pytest

from src.rgs_ec_algebraic_torsion import (
    DERIVED_CARTAN_COEFFICIENT, DIRAC_EC_CONVENTION_ID,
    axial_torsion_from_current, axial_torsion_from_dirac_spinor,
    dirac_axial_current, dirac_gamma_matrices_minus_plus, gamma5_matrix,
    cartan_algebraic_residual, torsion_quadratic_invariant,
)


def test_minus_plus_gamma_matrices_obey_clifford_algebra():
    gamma = dirac_gamma_matrices_minus_plus()
    eta = np.diag([-1.0, 1.0, 1.0, 1.0])
    residual = 0.0
    for a in range(4):
        for b in range(4):
            residual = max(residual, np.linalg.norm(
                gamma[a] @ gamma[b] + gamma[b] @ gamma[a]
                - 2.0 * eta[a, b] * np.eye(4)
            ))
    assert residual < 1e-14
    assert np.linalg.norm(gamma5_matrix() @ gamma5_matrix() - np.eye(4)) < 1e-14


def test_dirac_axial_current_is_real_and_spinor_derived():
    psi = np.array([1+0.2j, -0.3+0.7j, 0.4-0.1j, 0.8+0.5j])
    j5 = dirac_axial_current(psi)
    assert j5.shape == (4,)
    assert np.all(np.isfinite(j5))
    assert np.linalg.norm(dirac_axial_current(np.zeros(4))) == 0.0


def test_action_derived_cartan_coefficient_and_spinor_to_torsion_pipeline():
    psi = np.array([0.2+0.1j, 0.5-0.3j, -0.4+0.2j, 0.7+0.6j])
    result = axial_torsion_from_dirac_spinor(psi, 8.0, units_id="test-units")
    assert result["cartan_coefficient"] == DERIVED_CARTAN_COEFFICIENT == -0.5
    assert result["convention_id"] == DIRAC_EC_CONVENTION_ID
    assert result["totally_antisymmetric"]
    check = cartan_algebraic_residual(
        result["torsion"], result["spinor_axial_current"], 8.0,
        units_id="test-units",
    )
    assert check["equation_satisfied"]
    assert check["residual_norm"] < 1e-13


def test_torsion_invariant_finite_and_zero_for_zero_current():
    zero = axial_torsion_from_current(np.zeros(4), 1.0, units_id="u")["torsion"]
    assert torsion_quadratic_invariant(zero) == 0.0
    nonzero = axial_torsion_from_current([0.2, -0.1, 0.3, 0.4], 1.0)
    assert np.isfinite(torsion_quadratic_invariant(nonzero["torsion"]))


@pytest.mark.parametrize("kwargs", [{"convention_id": ""}, {"units_id": ""}])
def test_missing_convention_or_units_are_rejected(kwargs):
    args = {"axial_current": np.ones(4), "gravitational_coupling": 1.0}
    args.update(kwargs)
    with pytest.raises(ValueError):
        axial_torsion_from_current(**args)

import numpy as np
import pytest

from src.rgs_ec_algebraic_torsion import (
    DEFAULT_CONVENTION_ID, axial_torsion_from_current,
    audit_algebraic_axial_torsion, cartan_algebraic_residual,
    torsion_quadratic_invariant,
)


def test_axial_torsion_is_totally_antisymmetric_and_solves_parameterized_equation():
    j5 = np.array([0.4, -0.2, 0.1, 0.3])
    result = axial_torsion_from_current(
        j5, 8.0, cartan_coefficient=-0.25,
        convention_id=DEFAULT_CONVENTION_ID, units_id="test-units",
    )
    assert result["totally_antisymmetric"]
    assert result["antisymmetry_residual"] < 1e-13
    check = cartan_algebraic_residual(
        result["torsion"], j5, 8.0, cartan_coefficient=-0.25,
        convention_id=DEFAULT_CONVENTION_ID, units_id="test-units",
    )
    assert check["equation_satisfied"]
    assert check["residual_norm"] < 1e-13


def test_response_is_linear_in_axial_current_and_kappa():
    j5 = np.array([0.2, 0.3, -0.1, 0.5])
    a = axial_torsion_from_current(j5, 2.0, units_id="u")["torsion"]
    b = axial_torsion_from_current(3.0*j5, 2.0, units_id="u")["torsion"]
    c = axial_torsion_from_current(j5, 6.0, units_id="u")["torsion"]
    assert np.allclose(b, 3.0*a)
    assert np.allclose(c, 3.0*a)


def test_torsion_invariant_is_finite_and_zero_for_zero_current():
    zero = axial_torsion_from_current(
        np.zeros(4), 1.0, units_id="test-units"
    )["torsion"]
    assert torsion_quadratic_invariant(zero) == 0.0
    report = audit_algebraic_axial_torsion()
    assert np.isfinite(report["torsion_quadratic_invariant"])
    assert report["algebraic_equation_satisfied"]
    assert report["propagating_torsion_dofs"] == 0


@pytest.mark.parametrize("kwargs", [{"convention_id": ""}, {"units_id": ""}])
def test_missing_convention_or_units_are_rejected(kwargs):
    args = {"axial_current": np.ones(4), "gravitational_coupling": 1.0}
    args.update(kwargs)
    with pytest.raises(ValueError):
        axial_torsion_from_current(**args)

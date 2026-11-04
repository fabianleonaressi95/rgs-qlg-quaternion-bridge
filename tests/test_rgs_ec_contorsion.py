import numpy as np
import pytest

from src.rgs_ec_algebraic_torsion import (
    DERIVED_CARTAN_COEFFICIENT, axial_torsion_from_current,
)
from src.rgs_ec_contorsion import (
    audit_dirac_torsion_identities, contorsion_from_dirac_torsion,
    reconstruct_torsion_from_contorsion, torsion_quadratic_invariant,
)


def test_contorsion_reconstructs_totally_antisymmetric_torsion():
    j5 = np.array([0.2, -0.1, 0.3, 0.4])
    kappa = 2.3
    torsion = axial_torsion_from_current(j5, kappa)["torsion"]
    result = contorsion_from_dirac_torsion(torsion)
    reconstructed = reconstruct_torsion_from_contorsion(result["contorsion_lower"])
    assert np.linalg.norm(reconstructed - torsion) < 1e-13
    assert result["connection_antisymmetry_residual"] < 1e-13


def test_quadratic_torsion_invariant_matches_axial_current_identity():
    j5 = np.array([0.2, -0.1, 0.3, 0.4])
    kappa = 2.3
    torsion = axial_torsion_from_current(j5, kappa)["torsion"]
    audit = audit_dirac_torsion_identities(torsion, j5, kappa)
    assert audit["contorsion_reconstruction_residual"] < 1e-13
    assert audit["torsion_invariant_identity_residual"] < 1e-12
    assert np.isclose(
        torsion_quadratic_invariant(torsion),
        -1.5 * kappa**2 * (j5 @ np.diag([-1., 1., 1., 1.]) @ j5),
    )


def test_reject_non_antisymmetric_torsion_for_dirac_shortcut():
    torsion = np.zeros((4, 4, 4))
    torsion[0, 1, 2] = 1.0
    with pytest.raises(ValueError, match="totally antisymmetric"):
        contorsion_from_dirac_torsion(torsion)

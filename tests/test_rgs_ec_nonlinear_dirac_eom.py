import numpy as np
import pytest

from src.rgs_ec_algebraic_torsion import (
    dirac_adjoint, dirac_axial_current,
)
from src.rgs_ec_nonlinear_dirac_eom import (
    contact_lagrangian_from_independent_bilinears,
    effective_dirac_contact_sources,
    effective_dirac_eom_residual,
)


def test_contact_sources_match_independent_barpsi_finite_difference():
    psi = np.array([0.7 + 0.2j, -0.3 + 0.5j, 0.4 - 0.1j, 0.2 + 0.6j])
    bar = dirac_adjoint(psi)
    current = dirac_axial_current(psi)
    kappa = 0.8
    result = effective_dirac_contact_sources(psi, current, kappa)
    delta_bar = np.array([0.2 - 0.1j, -0.3j, 0.1 + 0.4j, -0.2 + 0.1j])
    h = 1e-6
    plus = contact_lagrangian_from_independent_bilinears(
        bar + h * delta_bar, psi, kappa
    )
    minus = contact_lagrangian_from_independent_bilinears(
        bar - h * delta_bar, psi, kappa
    )
    fd = (plus - minus) / (2 * h)
    analytic = np.dot(delta_bar, result["dL_dbarpsi"])
    assert fd == pytest.approx(float(np.real(analytic)), rel=2e-7, abs=2e-9)


def test_contact_sources_match_independent_psi_finite_difference():
    psi = np.array([0.7 + 0.2j, -0.3 + 0.5j, 0.4 - 0.1j, 0.2 + 0.6j])
    bar = dirac_adjoint(psi)
    current = dirac_axial_current(psi)
    kappa = 0.8
    result = effective_dirac_contact_sources(psi, current, kappa)
    delta_psi = np.array([0.1 + 0.3j, -0.2 + 0.1j, 0.4j, -0.1 - 0.2j])
    h = 1e-6
    plus = contact_lagrangian_from_independent_bilinears(
        bar, psi + h * delta_psi, kappa
    )
    minus = contact_lagrangian_from_independent_bilinears(
        bar, psi - h * delta_psi, kappa
    )
    fd = (plus - minus) / (2 * h)
    analytic = result["dL_dpsi"] @ delta_psi
    assert fd == pytest.approx(float(np.real(analytic)), rel=2e-7, abs=2e-9)


def test_nonlinear_eom_residual_adds_contact_source_and_does_not_claim_conservation():
    psi = np.array([1.0 + 0.2j, 0.3j, -0.4, 0.7 + 0.1j])
    current = dirac_axial_current(psi)
    kinetic = np.array([0.1j, 0.2, -0.3j, 0.4])
    result = effective_dirac_eom_residual(kinetic, psi, current, 0.5)
    source = effective_dirac_contact_sources(psi, current, 0.5)
    assert np.allclose(result["residual"], kinetic + source["dL_dbarpsi"])
    assert result["kinetic_eom_supplied_externally"] is True
    assert result["conservation_claim"].startswith("NOT_ESTABLISHED")


def test_zero_current_has_zero_contact_euler_source():
    psi = np.zeros(4, dtype=complex)
    result = effective_dirac_contact_sources(psi, np.zeros(4), 1.0)
    assert np.allclose(result["dL_dbarpsi"], 0.0)
    assert np.allclose(result["dL_dpsi"], 0.0)

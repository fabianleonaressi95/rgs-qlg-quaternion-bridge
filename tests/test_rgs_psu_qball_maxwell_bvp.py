import numpy as np

from src.rgs_psu_qball_bvp import qball_stationary_bvp
from src.rgs_psu_qball_maxwell_bvp import qball_maxwell_stationary_bvp


def test_zero_gauge_charge_recovers_global_qball_limit():
    global_report = qball_stationary_bvp(
        central_amplitude=1.0, radius_max=20.0, points=120
    )
    gauged_report = qball_maxwell_stationary_bvp(
        gauge_charge=0.0, central_amplitude=1.0,
        radius_max=20.0, points=120,
    )
    assert abs(gauged_report["omega"] - global_report["omega"]) < 2e-7
    assert abs(gauged_report["energy_per_particle_charge"]
               - global_report["energy_per_charge"]) < 2e-7
    assert gauged_report["maxwell_electric_energy"] < 1e-18
    assert gauged_report["gauss_law_residual"] < 1e-10
    assert np.max(np.abs(gauged_report["electrostatic_potential"])) < 1e-12


def test_gauged_qball_solves_scalar_and_maxwell_equations_and_gauss_law():
    report = qball_maxwell_stationary_bvp(
        gauge_charge=0.05, central_amplitude=1.0,
        radius_max=20.0, points=160,
    )
    assert report["converged"]
    assert report["profile_is_stationary_solution"]
    assert report["maxwell_backreaction_included"]
    assert report["frequency_window"]["satisfied"]
    assert report["origin_scalar_regularity_residual"] < 1e-10
    assert report["origin_maxwell_regularity_residual"] < 1e-10
    assert report["outer_scalar_residual"] < 1e-10
    assert report["outer_gauge_residual"] < 1e-10
    assert report["scalar_eom_residual_max"] < 2e-4
    assert report["maxwell_eom_residual_max"] < 2e-5
    assert report["gauss_law_residual"] < 2e-5
    assert report["electric_charge"] > 0.0
    assert abs(report["gauss_flux_charge"] - report["electric_charge"]) < 2e-5
    assert report["energy"] > 0.0
    assert report["noether_particle_charge"] > 0.0
    assert report["electric_charge"] < 0.0
    assert report["maxwell_electric_energy"] > 0.0
    assert np.all(np.isfinite(report["local_frequency"]))
    assert not report["vierbein_backreaction_included"]
    assert not report["linear_stability_established"]


def test_gauged_qball_rejects_invalid_inputs():
    import pytest

    with pytest.raises(ValueError):
        qball_maxwell_stationary_bvp(gauge_charge=float("nan"))
    with pytest.raises(ValueError):
        qball_maxwell_stationary_bvp(gauge_charge=0.0, points=20)

import numpy as np
import pytest

from src.rgs_psu_qball_bvp import qball_stationary_bvp
from src.rgs_psu_antitwin_qball_bvp import psu_antitwin_qball_maxwell_bvp


def test_two_field_zero_gauge_limit_recovers_mass_shifted_global_qball():
    params = dict(
        central_amplitude=0.5, radius_max=20.0, points=140,
        mass_scale=1.0, quartic=1.5, sextic=1.0,
        imbalance_penalty=0.2, exchange_coupling=0.1,
    )
    effective_mass = np.sqrt(params["mass_scale"]**2 - params["exchange_coupling"])
    global_report = qball_stationary_bvp(
        central_amplitude=params["central_amplitude"],
        radius_max=params["radius_max"], points=params["points"],
        mass_scale=effective_mass, quartic=params["quartic"],
        sextic=params["sextic"],
    )
    report = psu_antitwin_qball_maxwell_bvp(gauge_charge=0.0, **params)
    assert report["converged"]
    assert abs(report["omega"] - global_report["omega"]) < 2e-6
    assert abs(report["energy_per_particle_charge"]
               - global_report["energy_per_charge"]) < 2e-6
    assert report["maxwell_electric_energy"] < 1e-18
    assert report["gauss_law_residual"] < 1e-10
    assert report["exchange_symmetry_residual"] < 1e-10


def test_two_field_gauged_branch_solves_both_scalar_equations_and_gauss_law():
    report = psu_antitwin_qball_maxwell_bvp(
        gauge_charge=0.05, central_amplitude=0.5,
        radius_max=20.0, points=160,
        mass_scale=1.0, quartic=1.5, sextic=1.0,
        imbalance_penalty=0.2, exchange_coupling=0.1,
    )
    assert report["converged"]
    assert report["frequency_window"]["satisfied"]
    assert report["plus_eom_residual_max"] < 2e-4
    assert report["minus_eom_residual_max"] < 2e-4
    assert report["maxwell_eom_residual_max"] < 2e-5
    assert report["gauss_law_residual"] < 2e-5
    assert report["exchange_symmetry_residual"] < 1e-10
    assert report["noether_particle_charge"] > 0.0
    assert report["electric_charge"] > 0.0
    assert report["energy"] > 0.0
    assert report["maxwell_electric_energy"] > 0.0
    assert report["relative_amplitude_curvature_min"] == pytest.approx(
        np.min(report["relative_amplitude_curvature_min"]), abs=0.0
    )
    assert not report["full_linear_stability_established"]
    assert not report["vierbein_backreaction_included"]


def test_imbalance_penalty_changes_relative_curvature_not_symmetric_background():
    common = dict(
        gauge_charge=0.0, central_amplitude=0.5,
        radius_max=20.0, points=120,
        mass_scale=1.0, quartic=1.5, sextic=1.0,
        exchange_coupling=0.1,
    )
    no_penalty = psu_antitwin_qball_maxwell_bvp(imbalance_penalty=0.0, **common)
    with_penalty = psu_antitwin_qball_maxwell_bvp(imbalance_penalty=0.2, **common)
    assert abs(no_penalty["omega"] - with_penalty["omega"]) < 2e-6
    assert with_penalty["relative_amplitude_curvature_min"] > (
        no_penalty["relative_amplitude_curvature_min"]
    )


def test_two_field_qball_rejects_invalid_frequency_window_and_parameters():
    with pytest.raises(ValueError):
        psu_antitwin_qball_maxwell_bvp(exchange_coupling=0.5)
    with pytest.raises(ValueError):
        psu_antitwin_qball_maxwell_bvp(imbalance_penalty=-0.1)
    with pytest.raises(ValueError):
        psu_antitwin_qball_maxwell_bvp(central_amplitude=float("nan"))

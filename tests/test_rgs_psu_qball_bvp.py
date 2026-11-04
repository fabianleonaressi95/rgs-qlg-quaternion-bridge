import numpy as np

from src.rgs_psu_qball_bvp import qball_stationary_bvp


def test_stationary_qball_bvp_converges_to_nontrivial_regular_profile():
    report = qball_stationary_bvp(
        central_amplitude=1.0,
        radius_max=20.0,
        points=180,
        mass_scale=1.0,
        quartic=1.0,
        sextic=0.5,
    )
    assert report["converged"]
    assert report["profile_is_stationary_solution"]
    assert report["frequency_window"]["satisfied"]
    assert 0.70710678 < report["omega"] < 1.0
    assert report["origin_regular_boundary_residual"] < 1e-12
    assert report["outer_dirichlet_residual"] < 1e-12
    assert report["profile_nonnegative_within_tolerance"]
    assert report["profile_max"] == 1.0
    assert report["collocation_max_rms_residual"] < 2e-7
    assert report["radial_eom_residual_max"] < 5e-5
    assert report["energy"] > 0.0
    assert report["noether_charge"] > 0.0
    assert report["below_free_particle_threshold"]
    assert not report["linear_stability_established"]
    assert not report["nonlinear_stability_established"]


def test_stationary_qball_frequency_and_energy_converge_under_mesh_refinement():
    coarse = qball_stationary_bvp(
        central_amplitude=1.0, radius_max=20.0, points=120
    )
    fine = qball_stationary_bvp(
        central_amplitude=1.0, radius_max=20.0, points=180
    )
    assert abs(coarse["omega"] - fine["omega"]) < 2e-7
    assert abs(coarse["energy_per_charge"] - fine["energy_per_charge"]) < 2e-7
    assert abs(coarse["noether_charge"] - fine["noether_charge"]) < 2e-5


def test_stationary_qball_rejects_potential_without_positive_frequency_window():
    import pytest

    with pytest.raises(ValueError, match="positive Q-ball frequency lower bound"):
        qball_stationary_bvp(
            central_amplitude=1.0,
            mass_scale=1.0,
            quartic=2.0,
            sextic=0.5,
        )

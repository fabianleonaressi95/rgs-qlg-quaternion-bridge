import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.rgs_multiscale_relational_audit import (
    acoustic_metric_audit, discrete_scale_toy_audit, gravity_scale_audit,
    markov_entropy_audit, q6_audit, relational_clock_audit,
)


def test_q6_exact_graph_and_spectrum():
    r = q6_audit()
    assert (r["vertices"], r["undirected_edges"], r["betti_1_graph"]) == (64, 192, 129)
    assert r["pass"]


def test_finite_relational_clock_constraint():
    r = relational_clock_audit()
    assert r["global_constraint_residual"] < 1e-11
    assert r["conditional_state_distinguishability_theta_0_vs_0p37"] > 0
    assert "not calibrated proper time" in r["interpretation_limit"]
    assert r["pass"]


def test_markov_equilibrium_and_driven_cycle():
    r = markov_entropy_audit(seed=20261104, trajectories=500, horizon=2.0)
    assert r["equilibrium"]["entropy_production_rate_nats_per_model_time"] == 0.0
    assert r["driven_cycle"]["entropy_production_rate_nats_per_model_time"] > 0
    assert r["pass"]


def test_acoustic_metric_and_pure_vortex():
    r = acoustic_metric_audit(cs=1, radial_speed=0, circulation=1, radius=2)
    assert r["metric_determinant"] == -1
    assert r["acoustic_horizon_radius"] is None
    assert r["ergoradius"] == 1
    assert r["pass"]


def test_discrete_scale_toy_is_not_misrepresented_as_paper_replication():
    r = discrete_scale_toy_audit()
    assert r["pass"]
    assert r["status"] == "TOY_IDENTITY_ONLY_NOT_A_REPRODUCTION_OF_ECKER_ET_AL"


def test_dimensionless_gravity_scale_diagnostic():
    r = gravity_scale_audit()
    assert r["species"]["electron"]["gravitational_coupling_alpha_g"] < 1e-40
    assert abs(r["species"]["planck_mass"]["gravitational_coupling_alpha_g"] - 1) < 1e-12
    assert r["pass"]

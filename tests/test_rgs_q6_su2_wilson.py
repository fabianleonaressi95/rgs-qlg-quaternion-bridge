import numpy as np

from src.q6_complex import build_q6_complex
from src.rgs_q6_su2_wilson import (
    audit_rgs_q6_su2_wilson,
    recovered_q6_symmetry_subgroups,
    random_q6_connection,
    plaquette_holonomies,
    wilson_action,
)


def test_recovered_symmetry_subgroups_have_expected_orders():
    subgroups = recovered_q6_symmetry_subgroups()
    assert subgroups["O_h_order"] == 48
    assert subgroups["Q6_automorphism_order"] == 16
    assert subgroups["proper_Q6_subgroup_order"] == 8


def test_connection_and_plaquettes_have_canonical_dimensions():
    q6 = build_q6_complex()
    lattice = random_q6_connection(seed=1234)
    assert len(q6.edges) == 192
    assert len(lattice.edges()) == 192
    assert len(plaquette_holonomies(lattice)) == 240


def test_wilson_bridge_audit_passes():
    report = audit_rgs_q6_su2_wilson(seed=20261103, beta=1.25)
    assert report["pass"]
    assert report["wilson"]["gauge_action_residual"] < 1e-9
    assert report["wilson"]["plaquette_gauge_covariance_max_frobenius"] < 1e-9
    assert report["wilson"]["proper_Q6_symmetry_action_max_residual"] < 1e-8
    assert report["wilson"]["all_Q6_base_symmetry_action_max_residual"] < 1e-8

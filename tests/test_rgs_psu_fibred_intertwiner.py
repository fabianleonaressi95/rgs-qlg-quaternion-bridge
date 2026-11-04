import numpy as np

from src.rgs_psu_fibred_intertwiner import (
    audit_fibred_intertwiner,
    build_fibred_intertwiner_data,
)


def test_q6_base_sign_fibre_reconstructs_archived_rgs_seed_exactly():
    data = build_fibred_intertwiner_data()
    report = audit_fibred_intertwiner(data)
    assert report["factorization"]["base_group"] == "Z2^3"
    assert report["factorization"]["sign_fibre_group"] == "Z2^3"
    assert report["factorization"]["base_states"] == 8
    assert report["factorization"]["fibre_states_per_base"] == 8
    assert report["factorization"]["bijective_labelled_reconstruction"]
    assert report["factorization"]["rgs_coordinate_reconstruction_max_abs"] == 0.0


def test_pullback_edge_metric_and_graph_fem_are_positive():
    report = audit_fibred_intertwiner()
    checks = report["checks"]
    assert report["edge_geometry"]["total_edges"] == 192
    assert checks["96_base_plus_96_fibre_edges"]
    assert checks["golden_projection_edges_uniform"]
    assert checks["rgs_pullback_metric_nonuniform"]
    assert checks["base_and_fibre_subgraphs_each_rank_56"]
    assert checks["combined_graph_connected"]
    assert checks["mass_positive_definite"]
    assert checks["stiffness_positive_semidefinite"]
    assert checks["single_constant_zero_mode"]
    assert np.isclose(report["metric_graph_fem"]["generalized_spectral_gap"],
                      0.266804737, atol=1e-8)

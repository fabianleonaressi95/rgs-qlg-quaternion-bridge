from src.rgs_psu_carrier_connection_audit import audit_rgs_psu_carrier_connection

def test_archived_psu_connection_is_carrier_connection_in_local_radial_frames():
    report = audit_rgs_psu_carrier_connection()
    c = report["connection"]
    assert c["max_SO3_gauge_transform_residual"] < 1e-12
    assert c["max_SO3_orthogonality_residual"] < 1e-12
    assert c["max_SO3_determinant_residual"] < 1e-12
    assert c["max_radial_axis_fixed_residual"] < 1e-12
    assert c["max_SU2_unitarity_residual"] < 1e-12
    assert c["max_SU2_adjoint_to_SO3_residual"] < 1e-12
    assert c["max_SU2_local_offdiagonal"] < 1e-12

def test_raw_spatial_torsion_reproduces_prior_report_and_carrier_coframe_is_distinct():
    report = audit_rgs_psu_carrier_connection()
    t = report["Q6_square_torsion"]
    assert t["square_count"] == 240
    assert t["raw_spatial_torsion_matches_previous_report"]
    assert t["radial_plain_cochain_torsion_rms"] < 1e-12
    assert report["coframe_comparison"]["carrier_is_same_as_archived_4D_feature"] is False

def test_carrier_coframe_orientation_and_edge_bound_hold_on_all_q6_links():
    report = audit_rgs_psu_carrier_connection()
    o = report["orientation"]
    assert o["max_carrier_reverse_edge_residual"] < 1e-12
    assert o["max_norm_bound_excess"] < 1e-12
    assert report["edge_alignment"]["edge_set_equal"] is True

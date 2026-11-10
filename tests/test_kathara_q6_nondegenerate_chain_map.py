import numpy as np
from src.kathara_q6_nondegenerate_chain_map import audit_nondegenerate_kathara_q6_chain_map


def test_nondegenerate_kathara_q6_chain_identities():
    audit = audit_nondegenerate_kathara_q6_chain_map()
    assert audit.chain_identity_0 and audit.chain_identity_1


def test_face_map_is_nonzero_and_injective_on_face_chain_space():
    audit = audit_nondegenerate_kathara_q6_chain_map()
    assert audit.R2.shape == (240, 10)
    assert audit.R2_rank_mod_prime == 10
    assert np.all(np.count_nonzero(audit.R2, axis=0) > 0)


def test_dimensions_and_vertex_chart_are_preserved():
    audit = audit_nondegenerate_kathara_q6_chain_map()
    assert audit.R0.shape == (64, 12)
    assert audit.R1.shape == (192, 21)
    assert len(set(audit.vertex_labels)) == 12
    assert audit.mapped_rgs_coordinates.shape == (12, 3)
    assert audit.minimum_source_edge_length > 0


def test_each_face_has_a_nonzero_cubical_fill():
    audit = audit_nondegenerate_kathara_q6_chain_map()
    assert len(audit.face_fill_l1) == 10
    assert all(value > 0 for value in audit.face_fill_l1)

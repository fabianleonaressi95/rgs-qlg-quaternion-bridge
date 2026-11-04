"""Tests for local-frame representations on the Möbius dual graph."""
import numpy as np
import pytest

from src.mobius_representations import (
    audit_local_frame_representations,
    build_frame_transitions,
    representation_matrix,
    transform_transition,
    transition_holonomy,
    transport_feature,
)


def test_six_shared_edge_transitions_and_nontrivial_mobius_holonomy():
    transitions = build_frame_transitions()
    assert len(transitions) == 6
    assert sorted(item.sign for item in transitions) == [-1, 1, 1, 1, 1, 1]
    assert transition_holonomy(transitions, (0, 1, 2, 5, 4, 3, 0)) == -1


@pytest.mark.parametrize("kind", ["scalar", "sign", "regular"])
def test_each_representation_obeys_z2_homomorphism(kind):
    for a in (-1, 1):
        for b in (-1, 1):
            assert np.allclose(
                representation_matrix(kind, a * b),
                representation_matrix(kind, a) @ representation_matrix(kind, b),
                atol=1e-14,
            )


def test_scalar_sign_and_regular_cycle_actions():
    assert np.allclose(representation_matrix("scalar", -1), [[1]])
    assert np.allclose(representation_matrix("sign", -1), [[-1]])
    assert np.allclose(representation_matrix("regular", -1), [[0, 1], [1, 0]])
    assert np.allclose(representation_matrix("regular", -1) @
                       representation_matrix("regular", -1), np.eye(2))


def test_regular_feature_transport_swaps_channels_on_orientation_reversing_transition():
    value = np.array([2.0, 5.0])
    assert np.array_equal(transport_feature(value, -1, "regular"), [5.0, 2.0])
    assert np.array_equal(transport_feature(value, 1, "regular"), value)


def test_local_frame_change_transforms_transition_and_preserves_covariance():
    # f_j = rho(t) f_i; f'_i = rho(h_i) f_i; hence t' = h_j t h_i.
    for kind in ("scalar", "sign", "regular"):
        for t in (-1, 1):
            for hi in (-1, 1):
                for hj in (-1, 1):
                    f = np.array([0.7]) if kind != "regular" else np.array([0.7, -1.3])
                    lhs = representation_matrix(kind, transform_transition(t, hi, hj)) @ (
                        representation_matrix(kind, hi) @ f
                    )
                    rhs = representation_matrix(kind, hj) @ (
                        representation_matrix(kind, t) @ f
                    )
                    assert np.allclose(lhs, rhs, atol=1e-14)


def test_invalid_transition_and_feature_shapes_are_rejected():
    with pytest.raises(ValueError):
        representation_matrix("regular", 0)
    with pytest.raises(ValueError):
        transport_feature(np.array([1.0]), -1, "regular")


def test_audit_reports_orientation_obstruction_without_claiming_torsion():
    result = audit_local_frame_representations()
    assert result["dual_faces"] == 6
    assert result["shared_edge_transitions"] == 6
    assert result["dual_cycle_holonomy"] == -1
    assert result["orientation_obstruction_detected"] is True
    assert result["representation_homomorphism_max_residual"] == 0.0

def test_reflection_transitions_are_explicit_od_matrices_not_pin_lifts():
    from src.mobius_representations import (
        audit_reflection_transition_representation,
        reflection_transition_matrix,
    )
    report = audit_reflection_transition_representation()
    reflection = reflection_transition_matrix(-1, dimension=4, reflected_axis=2)
    assert report["is_Od_transition"] is True
    assert report["det_reflection"] == -1.0
    assert report["reflection_involution_fro"] == 0.0
    assert report["pin_lift_defined"] is False
    assert np.linalg.det(reflection) == -1.0

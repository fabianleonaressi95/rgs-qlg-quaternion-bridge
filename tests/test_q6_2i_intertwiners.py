import numpy as np
import pytest

from src.binary_icosahedral import binary_icosahedral_elements
from src.rgs_icosahedral_q6_intertwiner import icosahedral_signed_axis_representation
from src.q6_2i_intertwiners import (
    audit_q6_2i_intertwiner_compatibility,
    exterior_signed_permutation_matrix,
    fundamental_to_twisted_module_intertwiner,
)


def test_exterior_signed_permutation_is_orthogonal_and_preserves_vacuum():
    s = icosahedral_signed_axis_representation()["signed_permutations"][7]
    d = exterior_signed_permutation_matrix(s)
    assert d.shape == (64, 64)
    assert np.allclose(d.conj().T @ d, np.eye(64), atol=1e-14)
    assert np.allclose(d[:, 0], np.eye(64)[:, 0], atol=1e-14)


def test_central_character_forbids_untwisted_doublet_intertwiner():
    report = audit_q6_2i_intertwiner_compatibility()
    assert report["untwisted_intertwiner"]["nullity"] == 0
    assert report["untwisted_intertwiner"]["dimension"] == 0
    assert report["central_minus_one"]["fundamental_action"] == "-I_2"
    assert report["central_minus_one"]["Q6_exterior_action"] == "+I_64"


def test_doublet_twist_has_explicit_equivariant_vacuum_embedding():
    report = audit_q6_2i_intertwiner_compatibility()
    assert report["twisted_intertwiner"]["nonzero"]
    assert report["twisted_intertwiner"]["max_group_equivariance_residual"] < 2e-8
    assert report["passed"]


def test_canonical_embedding_has_correct_shape_and_rank():
    t = fundamental_to_twisted_module_intertwiner()
    assert t.shape == (128, 2)
    assert np.linalg.matrix_rank(t) == 2


def test_exterior_action_rejects_non_signed_permutations():
    with pytest.raises(ValueError):
        exterior_signed_permutation_matrix(np.zeros((6, 6), dtype=int))


def test_exterior_action_respects_signed_permutation_composition():
    signed = icosahedral_signed_axis_representation()["signed_permutations"]
    for i, j in [(0, 1), (3, 11), (7, 23), (19, 41), (52, 59)]:
        left = exterior_signed_permutation_matrix(signed[i])
        right = exterior_signed_permutation_matrix(signed[j])
        composed = exterior_signed_permutation_matrix(signed[i] @ signed[j])
        assert np.allclose(left @ right, composed, atol=1e-14)

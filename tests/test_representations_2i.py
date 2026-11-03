"""Tests for explicit binary-icosahedral representations and intertwiners."""
import numpy as np
import pytest

from src.binary_icosahedral import IDENTITY, binary_icosahedral_elements, qmul
from src.representations_2i import (
    intertwiner_basis,
    intertwiner_dimension_audit,
    quaternion_to_fundamental_su2,
    representation_audit,
    representation_matrix,
)


def test_fundamental_representation_matches_su2_lift_convention():
    q = binary_icosahedral_elements()[17]
    w, x, y, z = q
    expected = np.array([[w-1j*z, -y-1j*x], [y-1j*x, w+1j*z]])
    assert np.allclose(quaternion_to_fundamental_su2(q), expected, atol=1e-14)


def test_three_explicit_representations_are_homomorphisms():
    result = representation_audit()
    assert result["passed"] is True
    assert result["group_order"] == 120
    assert max(result["homomorphism_max_residual"].values()) < 2e-10
    assert max(result["unitarity_max_residual"].values()) < 2e-10


def test_central_minus_one_distinguishes_spinorial_and_vector_representations():
    minus_one = -IDENTITY
    assert np.allclose(representation_matrix(minus_one, "trivial"), [[1]])
    assert np.allclose(representation_matrix(minus_one, "fundamental"), -np.eye(2))
    assert np.allclose(representation_matrix(minus_one, "vector"), np.eye(3))


def test_intertwiner_dimensions_match_selected_irreducible_representations():
    result = intertwiner_dimension_audit()
    assert result["dimension_match"] is True
    assert result["dimensions"]["fundamental->fundamental"] == 1
    assert result["dimensions"]["vector->vector"] == 1
    assert result["dimensions"]["trivial->fundamental"] == 0
    assert result["dimensions"]["fundamental->vector"] == 0
    assert result["max_intertwiner_residual"] < 2e-9


def test_fundamental_commutant_is_complex_scalar_identity():
    basis, dimension, residual = intertwiner_basis("fundamental", "fundamental")
    assert dimension == 1
    assert np.allclose(basis[0], np.eye(2) / np.sqrt(2), atol=2e-9)
    assert residual < 2e-9


def test_invalid_group_element_shape_is_rejected():
    with pytest.raises(ValueError):
        representation_matrix(np.array([1.0, 0.0, 0.0]), "fundamental")

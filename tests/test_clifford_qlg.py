"""Tests for Lorentzian Clifford algebra and 2I-equivariant spinor kernels."""
import numpy as np
import pytest

from src.binary_icosahedral import binary_icosahedral_elements
from src.clifford_qlg import (
    ETA,
    I4,
    clifford_blade_basis,
    clifford_grade_projections,
    clifford_relation_audit,
    gamma_matrices,
    spinor_group_action_audit,
    spinor_lift_2i,
    spinor_steerable_kernel_audit,
)


def test_clifford_anticommutator_signature_minus_plus_plus_plus():
    result = clifford_relation_audit()
    assert result["passed"] is True
    assert result["signature"] == [-1, 1, 1, 1]
    assert result["blade_count"] == 16
    assert result["blade_basis_rank"] == 16
    assert result["max_anticommutator_residual"] < 2e-12


def test_grade_projections_reconstruct_scalar_vector_bivector_and_pseudoscalar():
    gs = gamma_matrices()
    value = 2.0*I4 + 3.0*gs[1] - 0.7*(gs[1]@gs[3]) + 1.2*(gs[0]@gs[1]@gs[2]@gs[3])
    parts = clifford_grade_projections(value)
    assert set(parts) == {0, 1, 2, 3, 4}
    assert np.allclose(sum(parts.values()), value, atol=2e-12)
    assert np.allclose(parts[0], 2.0*I4, atol=2e-12)
    assert np.allclose(parts[1], 3.0*gs[1], atol=2e-12)
    assert np.allclose(parts[2], -0.7*(gs[1]@gs[3]), atol=2e-12)
    assert np.allclose(parts[4], 1.2*(gs[0]@gs[1]@gs[2]@gs[3]), atol=2e-12)
    assert np.linalg.norm(parts[3]) < 2e-12


def test_2i_spinor_action_is_unitary_and_covariant_on_spatial_gamma_matrices():
    result = spinor_group_action_audit()
    assert result["passed"] is True
    assert result["group_order"] == 120
    assert result["homomorphism_max_residual"] < 3e-10
    assert result["gamma_spatial_covariance_max_residual"] < 3e-10
    assert result["central_minus_one_residual_to_minus_I4"] < 1e-12


def test_spinor_steerable_kernel_basis_is_four_dimensional():
    result = spinor_steerable_kernel_audit()
    assert result["passed"] is True
    assert result["intertwiner_dimension"] == 4
    assert result["basis_shape"] == [4, 4, 4]
    assert result["max_intertwining_residual"] < 2e-9


def test_spinor_lift_is_block_diagonal_copy_of_fundamental():
    q = binary_icosahedral_elements()[11]
    s = spinor_lift_2i(q)
    assert np.allclose(s[:2, :2], s[2:, 2:])
    assert np.allclose(s[:2, 2:], 0)
    assert np.allclose(s[2:, :2], 0)


def test_invalid_grade_projection_input_is_rejected():
    with pytest.raises(ValueError):
        clifford_grade_projections(np.eye(3))


def test_clifford_blade_basis_has_six_bivectors():
    labels, matrix = clifford_blade_basis()
    assert sum(len(x) == 2 for x in labels) == 6
    assert matrix.shape == (16, 16)

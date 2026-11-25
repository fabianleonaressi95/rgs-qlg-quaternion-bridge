import numpy as np
import pytest

from src.pellis_d6_operator_gate import (
    audit_d6_radial_lift,
    lifted_radial_transfer,
    parity_boundary_operator,
    sheet_parity_projectors,
)


def test_lifted_radial_transfer_commutes_with_full_d6_channel_action():
    m = np.array([[1.2, 0.3], [-0.4, 0.9]], dtype=complex)
    report = audit_d6_radial_lift(m, b=2.0, d=3.0)
    assert report["passed"]
    assert report["max_residual"] < 1e-12
    assert lifted_radial_transfer(m).shape == (12, 12)


def test_parity_boundary_operator_has_three_dimensional_branch_kernels():
    # b=0 selects the antisymmetric sheet sector; d=0 selects the symmetric one.
    b_zero = parity_boundary_operator(0.0, 2.0)
    d_zero = parity_boundary_operator(2.0, 0.0)
    assert np.linalg.matrix_rank(b_zero) == 3
    assert np.linalg.matrix_rank(d_zero) == 3
    assert np.linalg.matrix_rank(b_zero @ b_zero) == 3
    assert np.linalg.matrix_rank(d_zero @ d_zero) == 3


def test_sheet_parity_projectors_are_orthogonal_rank_three_projectors():
    p_plus, p_minus = sheet_parity_projectors()
    assert np.allclose(p_plus @ p_plus, p_plus)
    assert np.allclose(p_minus @ p_minus, p_minus)
    assert np.allclose(p_plus @ p_minus, 0)
    assert np.linalg.matrix_rank(p_plus) == 3
    assert np.linalg.matrix_rank(p_minus) == 3


def test_invalid_transfer_shapes_and_nonfinite_data_are_rejected():
    with pytest.raises(ValueError):
        lifted_radial_transfer(np.eye(3))
    with pytest.raises(ValueError):
        parity_boundary_operator(np.nan, 1.0)

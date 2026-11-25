import numpy as np
import pytest
from src.pellis_d6_commutant import audit_d6_commutant, d6_covariant_radial_transfer, mode_isotypic_projectors


def test_mode_projectors_are_orthogonal_rank_one_and_rank_two():
    p1, p2 = mode_isotypic_projectors()
    assert np.allclose(p1 @ p1, p1)
    assert np.allclose(p2 @ p2, p2)
    assert np.allclose(p1 @ p2, 0)
    assert np.allclose(p1 + p2, np.eye(3))
    assert np.linalg.matrix_rank(p1) == 1
    assert np.linalg.matrix_rank(p2) == 2


def test_d6_covariance_allows_distinct_axis_and_plane_transfers():
    mp = np.array([[1.0, 0.2], [-0.1, 0.8]])
    mt = np.array([[1.4, 0.0], [0.3, 0.7]])
    report = audit_d6_commutant(mp, mt)
    assert report["passed"] and report["max_residual"] < 1e-12
    assert not report["triplet_degenerate"]
    assert d6_covariant_radial_transfer(mp, mt).shape == (12, 12)


def test_threefold_degeneracy_is_the_special_equal_transfer_case():
    m = np.array([[1.1, 0.4], [-0.2, 0.9]])
    report = audit_d6_commutant(m, m)
    assert report["passed"] and report["triplet_degenerate"]


def test_nonfinite_and_wrong_shape_matrices_are_rejected():
    with pytest.raises(ValueError):
        d6_covariant_radial_transfer(np.eye(3), np.eye(2))
    with pytest.raises(ValueError):
        d6_covariant_radial_transfer(np.eye(2), np.full((2, 2), np.nan))

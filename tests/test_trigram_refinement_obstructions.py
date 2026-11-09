import numpy as np

from src.trigram_refinement_obstructions import (
    algebraic_realization, q3_adjacency, q6_lower_channel_compression,
    symmetric_channel_form, translation_permutation, tribonacci_target,
    walsh_first_order,
)


def test_translation_equivariance_forces_diagonal_first_order_compression():
    E, A3 = walsh_first_order(), q3_adjacency()
    R = 2*np.eye(8) + A3
    for g in [(a,b,c) for a in (0,1) for b in (0,1) for c in (0,1)]:
        T = translation_permutation(g)
        assert np.linalg.norm(R@T-T@R) < 1e-12
    M = E.T@R@E
    assert np.allclose(M, np.diag(np.diag(M)), atol=1e-12)


def test_q6_adjacency_compresses_to_four_identity_on_lower_channels():
    assert np.allclose(q6_lower_channel_compression(), 4*np.eye(3), atol=1e-12)


def test_tribonacci_breaks_translation_and_full_coordinate_permutation_symmetry():
    M = tribonacci_target()
    assert not np.allclose(M, np.diag(np.diag(M)))
    assert not symmetric_channel_form(M)


def test_coordinate_permutation_commutant_has_constant_diagonal_and_offdiagonal():
    M = 2*np.eye(3) + (np.ones((3,3))-np.eye(3))
    assert symmetric_channel_form(M)


def test_algebraic_realization_is_not_geometric_evidence():
    E, M = walsh_first_order(), tribonacci_target()
    R = algebraic_realization(M, lam=0.7)
    assert np.allclose(E.T@R@E, M, atol=1e-12)
    assert np.linalg.norm((np.eye(8)-E@E.T)@R@E) < 1e-12

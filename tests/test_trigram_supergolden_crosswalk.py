import numpy as np


def test_phase21_and_trigram_matrices_are_permutation_similar():
    A = np.array([[1, 0, 1], [1, 0, 0], [0, 1, 0]], dtype=int)
    M = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 1]], dtype=int)
    P = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], dtype=int)
    assert np.array_equal(P.T @ P, np.eye(3, dtype=int))
    assert np.array_equal(P @ A @ P.T, M)
    assert np.array_equal(np.poly(A), np.poly(M))
    assert np.array_equal(np.poly(M), np.array([1, -1, 0, -1]))


def test_shared_polynomial_does_not_encode_a_geometric_intertwiner():
    # The algebraic test deliberately contains no cell embedding or subdivision map.
    A = np.array([[1, 0, 1], [1, 0, 0], [0, 1, 0]], dtype=int)
    M = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 1]], dtype=int)
    P = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], dtype=int)
    assert np.array_equal(P @ A @ P.T, M)
    # This verifies only finite-dimensional similarity, not a map of geometric cells.

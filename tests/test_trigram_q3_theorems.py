import numpy as np

from src.trigram_q3_theorems import (
    TRIGRAMS, affine_face_compatibility, q3_adjacency,
    q6_product_adjacency, refinement_reduction, tribonacci_target, walsh_matrix,
)


def test_walsh_channels_are_orthonormal_and_invariant_under_q3_adjacency():
    E, A = walsh_matrix(), q3_adjacency()
    assert np.allclose(E.T @ E, np.eye(3), atol=1e-14)
    assert np.allclose(A @ E, E, atol=1e-14)
    assert sorted(np.linalg.eigvalsh(A).round(10).tolist()) == [
        -3.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 3.0
    ]


def test_q3_vertex_edge_and_square_face_counts():
    A = q3_adjacency()
    assert len(TRIGRAMS) == 8
    assert int(A.sum() // 2) == 12
    faces = [[t for t in TRIGRAMS if t[j] == bit]
             for j in range(3) for bit in (0, 1)]
    assert len(faces) == 6 and all(len(face) == 4 for face in faces)


def test_affine_face_compatibility_criterion():
    A, b = np.eye(3), np.zeros(3)
    A_beta = np.diag([2.0, 1.0, 1.0])
    assert affine_face_compatibility(A, b, A_beta, b, [0, 0, 0],
                                     [[0, 1, 0], [0, 0, 1]])
    assert not affine_face_compatibility(A, b, A_beta, b, [1, 0, 0],
                                         [[0, 1, 0], [0, 0, 1]])


def test_refinement_reduction_and_leakage_criterion():
    E, M = walsh_matrix(), tribonacci_target().astype(float)
    R_closed = E @ M @ E.T
    reduced, leakage = refinement_reduction(R_closed)
    assert np.allclose(reduced, M, atol=1e-14)
    assert leakage < 1e-14


def test_q6_is_product_of_two_q3_graphs():
    A6 = q6_product_adjacency()
    assert A6.shape == (64, 64)
    assert np.all(A6.sum(axis=1) == 6)
    assert int(A6.sum() // 2) == 192


def test_tribonacci_target_polynomial_and_perron_root():
    M = tribonacci_target()
    assert np.allclose(np.poly(M), [1, -1, 0, -1], atol=1e-12)
    psi = max(np.linalg.eigvals(M.astype(float)).real)
    assert abs(psi - 1.465571231876768) < 1e-12

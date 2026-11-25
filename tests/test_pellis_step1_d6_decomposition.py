import numpy as np

from src.pellis_d6_commutant import mode_isotypic_projectors


def test_three_mode_representation_splits_as_trivial_plus_irreducible_doublet():
    # C3 and F3 on the standard three-coordinate permutation module.
    c3 = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], dtype=complex)
    f3 = np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=complex)
    p0, p2 = mode_isotypic_projectors()
    omega = np.exp(2j * np.pi / 3)
    v_plus = np.array([1, omega, omega**2], dtype=complex)
    v_minus = np.array([1, omega**2, omega], dtype=complex)

    # Orthogonal invariant decomposition C^3 = span(1,1,1) ⊕ sum-zero plane.
    assert np.allclose(c3 @ p0, p0 @ c3)
    assert np.allclose(f3 @ p0, p0 @ f3)
    assert np.allclose(c3 @ p2, p2 @ c3)
    assert np.allclose(f3 @ p2, p2 @ f3)
    assert np.linalg.matrix_rank(p0) == 1
    assert np.linalg.matrix_rank(p2) == 2

    # The two C3 eigenspaces in the plane are exchanged by reflection.
    # For this orientation of the cyclic permutation, v_plus has eigenvalue
    # omega**2 and v_minus has eigenvalue omega.
    assert np.allclose(c3 @ v_plus, omega**2 * v_plus)
    assert np.allclose(c3 @ v_minus, omega * v_minus)
    assert np.allclose(f3 @ v_plus, v_minus)
    assert np.allclose(f3 @ v_minus, v_plus)

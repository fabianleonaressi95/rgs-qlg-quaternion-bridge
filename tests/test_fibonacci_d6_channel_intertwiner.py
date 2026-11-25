import numpy as np

from src.fibonacci_d6_channel_intertwiner import (
    audit_d6_channel_intertwiner,
    d6_channel_intertwiner,
)


def test_fibonacci_hexagon_has_exact_two_sheet_three_mode_factorization():
    report = audit_d6_channel_intertwiner()
    assert report["passed"]
    assert report["max_residual"] == 0.0


def test_half_turn_is_sheet_exchange_and_second_iterate_is_c3_rotation():
    data = d6_channel_intertwiner()
    phi = data["Phi"]
    r = data["r_orbit"]
    s = data["s_orbit"]
    rc = data["r_channel"]
    sc = data["s_channel"]
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)

    assert np.allclose(phi @ np.linalg.matrix_power(r, 3) @ phi.conj().T,
                       np.kron(sigma_x, np.eye(3)))
    assert np.allclose(phi @ np.linalg.matrix_power(r, 2) @ phi.conj().T,
                       np.kron(np.eye(2), data["C3"]))
    assert np.allclose(phi @ s @ phi.conj().T,
                       np.kron(np.eye(2), data["F3"]))
    assert np.allclose(sc @ rc @ sc, np.linalg.matrix_power(rc, 5))


def test_sheet_parity_projectors_have_rank_three():
    data = d6_channel_intertwiner()
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
    S = np.kron(sigma_x, np.eye(3))
    P_plus = (np.eye(6) + S) / 2
    P_minus = (np.eye(6) - S) / 2

    assert np.linalg.matrix_rank(P_plus) == 3
    assert np.linalg.matrix_rank(P_minus) == 3
    assert np.allclose(P_plus @ P_minus, 0)
    assert np.allclose(P_plus + P_minus, np.eye(6))


def test_triangle_modes_are_equivalent_to_spherical_l1_dipole_modes():
    from src.fibonacci_d6_channel_intertwiner import (
        audit_triangle_to_dipole_intertwiner,
        triangle_to_dipole_intertwiner,
    )

    report = audit_triangle_to_dipole_intertwiner()
    assert report["passed"]
    assert report["max_residual"] < 1e-12

    data = triangle_to_dipole_intertwiner()
    W = data["W3"]
    assert np.allclose(W @ data["C3"] @ W.T, data["dipole_rotation"])
    assert np.allclose(W @ data["F3"] @ W.T, data["dipole_reflection"])

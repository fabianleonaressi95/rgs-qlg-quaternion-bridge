"""Tests for the SU(2) representation and SU(2) -> SO(3) map."""

import math

import numpy as np
import pytest

from src.quaternions import Quaternion
from src.rotations import (
    axis_angle_to_quaternion,
    quaternion_to_matrix,
    rotate_vector,
    rotate_vector_quaternion,
    rotation_distance,
    is_rotation_matrix,
)
from src.su2_lifts import (
    quaternion_to_su2,
    su2_to_quaternion,
    su2_to_so3,
    su2_diagnostics,
    is_su2,
    opposite_lift,
    su2_lifts_of_rotation,
    su2_multiplication_residual,
    lift_sign_residual,
    projection_sign_residual,
    su2_so3_roundtrip_residual,
    quaternion_inverse_representation_residual,
)


TOL = 1e-10


def test_identity_maps_to_identity_matrix():
    q = Quaternion.identity()

    u = quaternion_to_su2(q)
    r = quaternion_to_matrix(q)

    assert np.allclose(u, np.eye(2), atol=TOL)
    assert np.allclose(r, np.eye(3), atol=TOL)


def test_su2_unitarity_and_determinant():
    q = axis_angle_to_quaternion(
        axis=(1.0, 2.0, -1.0),
        angle=1.234,
    )

    u = quaternion_to_su2(q)
    diagnostics = su2_diagnostics(u)

    assert diagnostics["unitarity_error"] <= TOL
    assert diagnostics["determinant_error"] <= TOL
    assert is_su2(u, TOL)


def test_quaternion_to_su2_and_back():
    q = axis_angle_to_quaternion(
        axis=(1.0, -2.0, 3.0),
        angle=0.83,
    )

    recovered = su2_to_quaternion(quaternion_to_su2(q))

    assert recovered.almost_equal(q, TOL)


def test_opposite_lifts_are_opposite_su2_matrices():
    q = axis_angle_to_quaternion(
        axis=(0.0, 0.0, 1.0),
        angle=0.7,
    )

    u_plus, u_minus = su2_lifts_of_rotation(q)

    assert np.allclose(u_minus, -u_plus, atol=TOL)


def test_opposite_lifts_have_same_so3_projection():
    q = axis_angle_to_quaternion(
        axis=(1.0, 1.0, 0.0),
        angle=1.1,
    )

    assert projection_sign_residual(q) <= TOL

    r_plus = su2_to_so3(quaternion_to_su2(q))
    r_minus = su2_to_so3(quaternion_to_su2(-q))

    assert np.allclose(r_plus, r_minus, atol=TOL)


def test_double_cover_identity():
    q = axis_angle_to_quaternion(
        axis=(0.0, 1.0, 0.0),
        angle=math.pi,
    )

    assert np.allclose(
        quaternion_to_matrix(q),
        quaternion_to_matrix(-q),
        atol=TOL,
    )


def test_two_pi_rotation_is_minus_identity_in_su2():
    q = axis_angle_to_quaternion(
        axis=(0.0, 0.0, 1.0),
        angle=2.0 * math.pi,
    )

    u = quaternion_to_su2(q)

    assert np.allclose(u, -np.eye(2), atol=TOL)
    assert np.allclose(quaternion_to_matrix(q), np.eye(3), atol=TOL)


def test_four_pi_rotation_returns_to_su2_identity():
    q = axis_angle_to_quaternion(
        axis=(0.0, 0.0, 1.0),
        angle=4.0 * math.pi,
    )

    u = quaternion_to_su2(q)

    assert np.allclose(u, np.eye(2), atol=TOL)


def test_quaternion_product_matches_su2_matrix_product():
    q1 = axis_angle_to_quaternion((1.0, 0.0, 0.0), 0.4)
    q2 = axis_angle_to_quaternion((0.0, 1.0, 0.0), 0.8)

    assert su2_multiplication_residual(q1, q2) <= TOL


def test_su2_inverse_matches_quaternion_inverse():
    q = axis_angle_to_quaternion(
        axis=(1.0, 2.0, 3.0),
        angle=0.9,
    )

    assert quaternion_inverse_representation_residual(q) <= TOL


def test_rotation_matrix_is_in_so3():
    q = axis_angle_to_quaternion(
        axis=(2.0, -1.0, 4.0),
        angle=2.1,
    )

    r = quaternion_to_matrix(q)

    assert is_rotation_matrix(r, TOL)
    assert np.linalg.det(r) == pytest.approx(1.0, abs=TOL)


def test_matrix_action_matches_quaternion_action():
    q = axis_angle_to_quaternion(
        axis=(1.0, 2.0, 3.0),
        angle=1.7,
    )
    vector = np.array([0.3, -1.2, 2.4])

    from_matrix = rotate_vector(q, vector)
    from_quaternion = rotate_vector_quaternion(q, vector)

    assert np.allclose(from_matrix, from_quaternion, atol=TOL)


def test_rotation_about_z_axis():
    q = axis_angle_to_quaternion(
        axis=(0.0, 0.0, 1.0),
        angle=math.pi / 2.0,
    )

    result = rotate_vector(q, (1.0, 0.0, 0.0))

    assert np.allclose(result, (0.0, 1.0, 0.0), atol=TOL)


def test_rotation_composition_order():
    qx = axis_angle_to_quaternion((1.0, 0.0, 0.0), math.pi / 2.0)
    qy = axis_angle_to_quaternion((0.0, 1.0, 0.0), math.pi / 2.0)

    q_composed = qx * qy

    vector = np.array([0.0, 0.0, 1.0])

    sequential = rotate_vector(
        qx,
        rotate_vector(qy, vector),
    )
    composed = rotate_vector(q_composed, vector)

    assert np.allclose(sequential, composed, atol=TOL)


def test_rotation_distance_ignores_lift_sign():
    q = axis_angle_to_quaternion(
        axis=(1.0, 0.0, 0.0),
        angle=1.2,
    )

    assert rotation_distance(q, -q) == pytest.approx(0.0, abs=TOL)


def test_rotation_distance_is_bounded_by_pi():
    q1 = axis_angle_to_quaternion((1.0, 0.0, 0.0), 0.3)
    q2 = axis_angle_to_quaternion((0.0, 1.0, 0.0), 2.4)

    distance = rotation_distance(q1, q2)

    assert 0.0 <= distance <= math.pi + TOL


def test_lift_sign_residual():
    q = axis_angle_to_quaternion(
        axis=(1.0, 2.0, 0.0),
        angle=1.5,
    )

    assert lift_sign_residual(q) <= TOL


def test_so3_roundtrip_residual():
    q = axis_angle_to_quaternion(
        axis=(1.0, -3.0, 2.0),
        angle=2.3,
    )

    assert su2_so3_roundtrip_residual(q) <= TOL


def test_opposite_lift_function():
    q = Quaternion(0.5, 0.5, 0.5, 0.5)

    assert opposite_lift(q).almost_equal(-q, TOL)


def test_invalid_su2_matrix_is_rejected():
    invalid = np.array(
        [[2.0, 0.0], [0.0, 0.5]],
        dtype=np.complex128,
    )

    assert not is_su2(invalid, TOL)

    with pytest.raises(ValueError):
        su2_to_quaternion(invalid)

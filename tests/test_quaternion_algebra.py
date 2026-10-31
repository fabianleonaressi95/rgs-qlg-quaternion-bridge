"""Tests for the Hamilton quaternion algebra."""

import math

import pytest

from src.quaternions import (
    Quaternion,
    hamilton_product,
    quaternion_conjugate,
    quaternion_dot,
    quaternion_inverse,
    same_rotation,
)


TOL = 1e-12


def test_identity_element():
    q = Quaternion(0.5, 0.5, 0.5, 0.5)
    identity = Quaternion.identity()

    assert (q * identity).almost_equal(q, TOL)
    assert (identity * q).almost_equal(q, TOL)


def test_zero_element():
    q = Quaternion(1.0, -2.0, 3.0, 4.0)
    zero = Quaternion.zero()

    assert (q + zero).almost_equal(q, TOL)
    assert (q * zero).almost_equal(zero, TOL)


def test_conjugate():
    q = Quaternion(1.0, 2.0, 3.0, 4.0)
    expected = Quaternion(1.0, -2.0, -3.0, -4.0)

    assert quaternion_conjugate(q).almost_equal(expected, TOL)


def test_conjugate_product_reverses_order():
    q1 = Quaternion(1.0, 2.0, 3.0, 4.0)
    q2 = Quaternion(0.5, -1.0, 2.0, 0.25)

    left = (q1 * q2).conjugate()
    right = q2.conjugate() * q1.conjugate()

    assert left.almost_equal(right, TOL)


def test_inverse_left_and_right():
    q = Quaternion(1.0, 2.0, -3.0, 4.0)
    identity = Quaternion.identity()
    inverse = quaternion_inverse(q)

    assert (q * inverse).almost_equal(identity, TOL)
    assert (inverse * q).almost_equal(identity, TOL)


def test_inverse_of_unit_quaternion_is_conjugate():
    q = Quaternion(0.5, 0.5, 0.5, 0.5)

    assert q.inverse().almost_equal(q.conjugate(), TOL)


def test_zero_quaternion_has_no_inverse():
    with pytest.raises(ZeroDivisionError):
        Quaternion.zero().inverse()


def test_zero_quaternion_cannot_be_normalized():
    with pytest.raises(ZeroDivisionError):
        Quaternion.zero().normalized()


def test_norm_squared():
    q = Quaternion(1.0, 2.0, 3.0, 4.0)

    assert q.norm_squared() == pytest.approx(30.0)
    assert q.norm() == pytest.approx(math.sqrt(30.0))


def test_norm_multiplicativity():
    q1 = Quaternion(1.0, 2.0, 3.0, 4.0)
    q2 = Quaternion(-2.0, 0.5, 1.0, -3.0)

    assert (q1 * q2).norm() == pytest.approx(
        q1.norm() * q2.norm(),
        abs=TOL,
    )


def test_hamilton_product_is_associative():
    q1 = Quaternion(1.0, 2.0, 0.0, -1.0)
    q2 = Quaternion(0.5, -1.0, 3.0, 2.0)
    q3 = Quaternion(-2.0, 0.25, 1.0, 0.5)

    left = (q1 * q2) * q3
    right = q1 * (q2 * q3)

    assert left.almost_equal(right, TOL)


def test_hamilton_product_is_not_commutative():
    i = Quaternion(0.0, 1.0, 0.0, 0.0)
    j = Quaternion(0.0, 0.0, 1.0, 0.0)

    ij = i * j
    ji = j * i

    assert ij.almost_equal(
        Quaternion(0.0, 0.0, 0.0, 1.0),
        TOL,
    )
    assert ji.almost_equal(
        Quaternion(0.0, 0.0, 0.0, -1.0),
        TOL,
    )
    assert not ij.almost_equal(ji, TOL)


def test_dot_product():
    q1 = Quaternion(1.0, 2.0, 3.0, 4.0)
    q2 = Quaternion(2.0, 0.0, -1.0, 3.0)

    assert quaternion_dot(q1, q2) == pytest.approx(13.0)


def test_normalization():
    q = Quaternion(1.0, 2.0, 2.0, 4.0)
    unit = q.normalized()

    assert unit.norm() == pytest.approx(1.0, abs=TOL)


def test_scalar_multiplication():
    q = Quaternion(1.0, 2.0, 3.0, 4.0)

    assert (2.0 * q).almost_equal(
        Quaternion(2.0, 4.0, 6.0, 8.0),
        TOL,
    )


def test_from_iterable_requires_four_components():
    with pytest.raises(ValueError):
        Quaternion.from_iterable((1.0, 2.0, 3.0))


def test_non_finite_components_are_rejected():
    with pytest.raises(ValueError):
        Quaternion(float("nan"), 0.0, 0.0, 0.0)


def test_opposite_lifts_are_distinct_quaternions():
    q = Quaternion(0.5, 0.5, 0.5, 0.5)
    minus_q = -q

    assert not q.almost_equal(minus_q, TOL)
    assert same_rotation(q, minus_q, TOL)


def test_inverse_is_not_the_opposite_lift_of_conjugate():
    q = Quaternion(0.5, 0.5, 0.5, 0.5)

    inverse = q.inverse()
    opposite_lift_of_inverse = -q.conjugate()

    assert inverse.almost_equal(q.conjugate(), TOL)
    assert not inverse.almost_equal(opposite_lift_of_inverse, TOL)

    assert (q * inverse).almost_equal(
        Quaternion.identity(),
        TOL,
    )

    assert (q * opposite_lift_of_inverse).almost_equal(
        -Quaternion.identity(),
        TOL,
    )


@pytest.mark.parametrize(
    "q",
    [
        Quaternion.identity(),
        Quaternion(0.5, 0.5, 0.5, 0.5),
        Quaternion(0.0, 1.0, 0.0, 0.0),
        Quaternion(1.0, -2.0, 3.0, -4.0),
    ],
)
def test_functional_hamilton_product(q):
    identity = Quaternion.identity()

    assert hamilton_product(q, identity).almost_equal(q, TOL)

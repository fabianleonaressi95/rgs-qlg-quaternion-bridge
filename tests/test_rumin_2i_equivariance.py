"""Finite convention check for the commuting left 2I / right Hopf actions.

This test does not test the analytic Rumin operator. It only checks the
quaternionic action convention used by the geometric descent theorem.
"""
import numpy as np

from src.binary_icosahedral import binary_icosahedral_elements, qmul


def _unit_quaternion(axis, angle):
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    return np.concatenate(([np.cos(angle / 2.0)], np.sin(angle / 2.0) * axis))


def test_left_2i_action_commutes_with_right_hopf_circle():
    elements = binary_icosahedral_elements()
    samples = [
        np.array([1.0, 0.0, 0.0, 0.0]),
        _unit_quaternion([1.0, 2.0, -1.0], 0.731),
        _unit_quaternion([-2.0, 1.0, 3.0], 1.217),
    ]
    angles = (0.0, 0.17, 0.5 * np.pi, np.pi, 1.73)
    hopf_axis = np.array([0.0, 1.0, 0.0, 0.0])  # fixed imaginary unit i

    max_commutator = 0.0
    max_norm_error = 0.0
    for g in elements:
        for q in samples:
            for angle in angles:
                h = _unit_quaternion(hopf_axis[1:], angle)
                lhs = qmul(g, qmul(q, h))
                rhs = qmul(qmul(g, q), h)
                max_commutator = max(max_commutator, float(np.linalg.norm(lhs-rhs)))
                max_norm_error = max(
                    max_norm_error,
                    abs(float(np.linalg.norm(lhs))-1.0),
                    abs(float(np.linalg.norm(rhs))-1.0),
                )
    assert max_commutator < 2e-12
    assert max_norm_error < 2e-12

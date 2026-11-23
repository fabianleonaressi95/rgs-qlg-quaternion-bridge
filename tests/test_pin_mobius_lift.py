import numpy as np
import pytest
from src.discrete_spin_geometry import euclidean_gamma_matrices
from src.pin_mobius_lift import audit_mobius_pin_plus_lift, chirality_operator, clifford_transition, pin_plus_lift, pin_twisted_adjoint_vector, reflection_matrix

def test_pin_plus_generator_squares_to_positive_identity():
    p = pin_plus_lift(-1, 0)
    assert np.allclose(p @ p, np.eye(4), atol=1e-14)
    assert np.allclose(pin_plus_lift(1, 0), np.eye(4), atol=1e-14)

def test_pin_twisted_adjoint_projects_to_orientation_reflection():
    gamma, r = euclidean_gamma_matrices(), reflection_matrix(-1, 0)
    for j in range(4):
        expected = sum((r[k, j] * gamma[k] for k in range(4)), np.zeros((4, 4), complex))
        assert np.allclose(pin_twisted_adjoint_vector(-1, gamma[j]), expected, atol=1e-14)

def test_chirality_implements_grade_automorphism():
    gamma, gamma5 = euclidean_gamma_matrices(), chirality_operator()
    assert np.allclose(gamma5 @ gamma5, np.eye(4), atol=1e-14)
    for g in gamma:
        assert np.allclose(gamma5 @ g @ gamma5.conj().T, -g, atol=1e-14)

def test_clifford_transition_is_an_algebra_automorphism_and_z2_action():
    gamma = euclidean_gamma_matrices()
    for a in (-1, 1):
        for b in (-1, 1):
            for x in [np.eye(4), *gamma, gamma[0] @ gamma[1]]:
                assert np.allclose(clifford_transition(a, clifford_transition(b, x)), clifford_transition(a*b, x), atol=1e-14)
    x, y = gamma[0] + 2*gamma[1], gamma[2] @ gamma[3]
    assert np.allclose(clifford_transition(-1, x@y), clifford_transition(-1, x) @ clifford_transition(-1, y), atol=1e-14)

def test_mobius_dual_cycle_has_consistent_pin_plus_holonomy_and_cocycle():
    result = audit_mobius_pin_plus_lift()
    assert result["base_orientation_holonomy"] == -1
    assert result["pin_type"] == "Pin+(4)"
    assert result["pin_holonomy_equals_lift_of_base_holonomy"] is True
    assert result["passed"] is True
    assert result["max_residual"] < 2e-12

def test_invalid_pin_transition_is_rejected():
    with pytest.raises(ValueError):
        pin_plus_lift(0)

"""Pin+ lift of the finite Möbius orientation transitions.

The lift is defined on the six-face dual transition graph, not asserted to be
a Pin structure on an unspecified tangent bundle over the full 2-complex.

Convention: Cl(4,0), Hermitian generators gamma_a with gamma_a^2=+1.
The Pin+ element for the orientation-reversing transition is gamma_axis.
Its twisted adjoint on vectors is v -> -p v p^{-1}. The induced
Clifford-algebra transition is Phi_p(X)=p alpha(X) p^{-1}, where alpha is
the grade automorphism. In the 4D matrix representation alpha is implemented
by the chirality volume element Gamma=gamma_0 gamma_1 gamma_2 gamma_3.
"""
from __future__ import annotations
from typing import Iterable
import numpy as np
from src.discrete_spin_geometry import euclidean_gamma_matrices
from src.mobius_representations import build_frame_transitions, transition_holonomy, transform_transition

I4 = np.eye(4, dtype=complex)

def pin_plus_lift(sign: int, reflected_axis: int = 0) -> np.ndarray:
    """Lift a Z2 orientation sign to Pin+(4) in Cl(4,0)."""
    if sign not in (-1, 1):
        raise ValueError("A Z2 transition sign must be +1 or -1.")
    if not 0 <= reflected_axis < 4:
        raise ValueError("reflected_axis must be in range(4).")
    return I4.copy() if sign == 1 else euclidean_gamma_matrices()[reflected_axis].copy()

def reflection_matrix(sign: int, reflected_axis: int = 0) -> np.ndarray:
    if sign not in (-1, 1) or not 0 <= reflected_axis < 4:
        raise ValueError("Invalid Z2 sign or reflected axis.")
    result = np.eye(4)
    if sign == -1:
        result[reflected_axis, reflected_axis] = -1.0
    return result

def chirality_operator() -> np.ndarray:
    gamma = euclidean_gamma_matrices()
    return gamma[0] @ gamma[1] @ gamma[2] @ gamma[3]

def pin_twisted_adjoint_vector(sign: int, vector: np.ndarray, reflected_axis: int = 0) -> np.ndarray:
    """Apply the Pin twisted-adjoint O(4) action to a Clifford vector."""
    value = np.asarray(vector, dtype=complex)
    if value.shape != (4, 4):
        raise ValueError("vector must be a 4x4 Clifford-vector matrix.")
    if sign == 1:
        return value.copy()
    if sign != -1:
        raise ValueError("A Z2 transition sign must be +1 or -1.")
    p = pin_plus_lift(sign, reflected_axis)
    return -p @ value @ p.conj().T

def clifford_transition(sign: int, value: np.ndarray, reflected_axis: int = 0) -> np.ndarray:
    """Clifford-algebra transport p alpha(X) p^{-1}; on vectors it is reflection."""
    x = np.asarray(value, dtype=complex)
    if x.shape != (4, 4):
        raise ValueError("value must be a 4x4 Clifford-algebra matrix.")
    if sign == 1:
        return x.copy()
    if sign != -1:
        raise ValueError("A Z2 transition sign must be +1 or -1.")
    p = pin_plus_lift(sign, reflected_axis)
    gamma5 = chirality_operator()
    return p @ (gamma5 @ x @ gamma5.conj().T) @ p.conj().T

def _ordered_pin_holonomy(transitions: Iterable[object], cycle: tuple[int, ...], reflected_axis: int) -> np.ndarray:
    lookup = {}
    for item in transitions:
        lookup[(item.source_face, item.target_face)] = item.sign
        lookup[(item.target_face, item.source_face)] = item.sign
    result = I4.copy()
    for source, target in zip(cycle[:-1], cycle[1:]):
        result = pin_plus_lift(lookup[(source, target)], reflected_axis) @ result
    return result

def audit_mobius_pin_plus_lift(reflected_axis: int = 0, tol: float = 2e-12) -> dict[str, object]:
    """Audit Pin projection, Pin+ relations, lifted cocycle and Clifford covariance."""
    if not 0 <= reflected_axis < 4:
        raise ValueError("reflected_axis must be in range(4).")
    gamma = euclidean_gamma_matrices()
    transitions = build_frame_transitions()
    cycle = (0, 1, 2, 5, 4, 3, 0)
    base_holonomy = transition_holonomy(transitions, cycle)
    p = pin_plus_lift(-1, reflected_axis)
    gamma5 = chirality_operator()
    pin_square = float(np.linalg.norm(p @ p - I4))
    chirality_square = float(np.linalg.norm(gamma5 @ gamma5 - I4))
    grade_residual = max(float(np.linalg.norm(gamma5 @ g @ gamma5.conj().T + g)) for g in gamma)
    projection_residual = 0.0
    clifford_residual = 0.0
    r = reflection_matrix(-1, reflected_axis)
    for j, g in enumerate(gamma):
        expected = sum((r[k, j] * gamma[k] for k in range(4)), np.zeros((4, 4), complex))
        projection_residual = max(projection_residual, float(np.linalg.norm(pin_twisted_adjoint_vector(-1, g, reflected_axis) - expected)))
        clifford_residual = max(clifford_residual, float(np.linalg.norm(clifford_transition(-1, g, reflected_axis) - expected)))
    local_residual = 0.0
    for item in transitions:
        for hs in (-1, 1):
            for ht in (-1, 1):
                changed = transform_transition(item.sign, hs, ht)
                lhs = pin_plus_lift(ht, reflected_axis) @ pin_plus_lift(item.sign, reflected_axis) @ pin_plus_lift(hs, reflected_axis).conj().T
                local_residual = max(local_residual, float(np.linalg.norm(lhs - pin_plus_lift(changed, reflected_axis))))
    pin_holonomy = _ordered_pin_holonomy(transitions, cycle, reflected_axis)
    expected_holonomy = pin_plus_lift(base_holonomy, reflected_axis)
    algebra_law = 0.0
    probes = [I4, *gamma]
    for a in (-1, 1):
        for b in (-1, 1):
            for probe in probes:
                lhs = clifford_transition(a, clifford_transition(b, probe, reflected_axis), reflected_axis)
                rhs = clifford_transition(a * b, probe, reflected_axis)
                algebra_law = max(algebra_law, float(np.linalg.norm(lhs - rhs)))
    residuals = {
        "pin_plus_square": pin_square,
        "chirality_square": chirality_square,
        "grade_automorphism": grade_residual,
        "twisted_adjoint_projection": projection_residual,
        "clifford_vector_covariance": clifford_residual,
        "local_pin_cocycle": local_residual,
        "cycle_holonomy_lift": float(np.linalg.norm(pin_holonomy - expected_holonomy)),
        "clifford_transition_group_law": algebra_law,
    }
    return {
        "audit": "MOBIUS_PIN_PLUS_CLIFFORD_LIFT",
        "clifford_signature": [1, 1, 1, 1],
        "pin_type": "Pin+(4)",
        "reflected_axis": reflected_axis,
        "dual_cycle": list(cycle),
        "base_orientation_holonomy": base_holonomy,
        "pin_holonomy_equals_lift_of_base_holonomy": bool(np.allclose(pin_holonomy, expected_holonomy, atol=tol, rtol=tol)),
        "pin_holonomy_square_residual": float(np.linalg.norm(pin_holonomy @ pin_holonomy - I4)),
        "residuals": residuals,
        "max_residual": max(residuals.values()),
        "passed": max(residuals.values()) <= tol and base_holonomy == -1,
        "scope": "Pin+ lift and Clifford-algebra transport on the specified six-face dual transition graph; not a full tangent-bundle Pin structure",
    }

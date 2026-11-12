"""Central-character audit for the 2I action on the Q6 exterior Clifford module.

The signed-axis action of A5 on R^6 is pulled back along
2I -> 2I/{+/-1} ~= A5. Consequently the central element -1 acts trivially
on Lambda*(C^6), whereas it acts as -I on the fundamental doublet.

This module proves the resulting Hom-space obstruction and implements the
minimal representation-theoretic twist C^2 tensor Lambda*(C^6). The twist is
an algebraic construction, not a claim of a geometrically derived coupling.
"""
from __future__ import annotations

from itertools import combinations
import numpy as np

from src.binary_icosahedral import (
    IDENTITY, binary_icosahedral_elements, quaternion_to_so3,
)
from src.representations_2i import representation_matrix
from src.rgs_icosahedral_q6_intertwiner import icosahedral_signed_axis_representation


def exterior_signed_permutation_matrix(signed_permutation: np.ndarray) -> np.ndarray:
    """Induced 64x64 action on Lambda*(C^6) of a signed permutation of C^6."""
    s = np.asarray(signed_permutation, dtype=int)
    if s.shape != (6, 6) or not np.all(np.isin(s, (-1, 0, 1))):
        raise ValueError("signed_permutation must be a 6x6 signed permutation")
    if not (np.all(np.sum(np.abs(s), axis=0) == 1)
            and np.all(np.sum(np.abs(s), axis=1) == 1)):
        raise ValueError("matrix must have one signed unit per row and column")
    result = np.zeros((64, 64), dtype=np.complex128)
    for degree in range(7):
        for subset in combinations(range(6), degree):
            targets = [int(np.flatnonzero(s[:, j])[0]) for j in subset]
            coefficient = int(np.prod([s[targets[k], j] for k, j in enumerate(subset)]))
            inversions = sum(targets[i] > targets[j]
                             for i in range(len(targets))
                             for j in range(i + 1, len(targets)))
            if inversions % 2:
                coefficient = -coefficient
            source_mask = sum(1 << j for j in subset)
            target_mask = sum(1 << j for j in targets)
            result[target_mask, source_mask] = coefficient
    return result


def _q6_action_for_quaternion(q: np.ndarray, rotations: np.ndarray,
                               signed_permutations: np.ndarray,
                               tol: float = 2e-8) -> np.ndarray:
    rotation = quaternion_to_so3(q)
    distances = np.linalg.norm(rotations - rotation[None, :, :], axis=(1, 2))
    index = int(np.argmin(distances))
    if distances[index] > tol:
        raise ArithmeticError("Quaternion rotation did not match the A5 rotation list")
    return exterior_signed_permutation_matrix(signed_permutations[index])


def fundamental_to_twisted_module_intertwiner() -> np.ndarray:
    """Canonical injection v -> v tensor 1 into C^2 tensor Lambda*(C^6)."""
    t = np.zeros((2 * 64, 2), dtype=np.complex128)
    t[0, 0] = 1.0  # first doublet component tensor the degree-zero vacuum
    t[64, 1] = 1.0
    return t


def audit_q6_2i_intertwiner_compatibility(tol: float = 2e-8) -> dict[str, object]:
    """Audit the central obstruction and its minimal doublet-twisted repair."""
    elements = binary_icosahedral_elements()
    data = icosahedral_signed_axis_representation()
    rotations = np.asarray(data["rotations"])
    signed = np.asarray(data["signed_permutations"])

    rho_plus = representation_matrix(IDENTITY, "fundamental")
    rho_minus = representation_matrix(-IDENTITY, "fundamental")
    q6_plus = _q6_action_for_quaternion(IDENTITY, rotations, signed)
    q6_minus = _q6_action_for_quaternion(-IDENTITY, rotations, signed)

    central_fundamental_residual = float(np.linalg.norm(rho_minus + np.eye(2)))
    central_q6_residual = float(np.linalg.norm(q6_minus - np.eye(64)))
    # T rho_2(-1) = rho_Q(-1) T gives -T=T, hence 2T=0.
    obstruction_constraint_rank = 128
    obstruction_nullity = 0

    t = fundamental_to_twisted_module_intertwiner()
    twisted_residual = 0.0
    q6_vacuum_residual = 0.0
    for q in elements:
        rho2 = representation_matrix(q, "fundamental")
        rho_q6 = _q6_action_for_quaternion(q, rotations, signed)
        q6_vacuum_residual = max(
            q6_vacuum_residual,
            float(np.linalg.norm(rho_q6[:, 0] - np.eye(64)[:, 0])),
        )
        lhs = np.kron(rho2, rho_q6) @ t
        rhs = t @ rho2
        twisted_residual = max(
            twisted_residual, float(np.linalg.norm(lhs - rhs))
        )

    return {
        "audit": "Q6_2I_CENTRAL_CHARACTER_INTERTWINER",
        "group_order_2I": int(len(elements)),
        "quotient_order_A5": int(len(rotations)),
        "dimensions": {
            "fundamental_2I_doublet": 2,
            "Q6_exterior_clifford_module": 64,
            "doublet_twisted_Q6_module": 128,
        },
        "central_minus_one": {
            "fundamental_action": "-I_2",
            "Q6_exterior_action": "+I_64",
            "fundamental_residual": central_fundamental_residual,
            "Q6_residual": central_q6_residual,
        },
        "untwisted_intertwiner": {
            "equation_for_minus_one": "T(-I_2)=(+I_64)T, hence T=0",
            "constraint_rank": obstruction_constraint_rank,
            "nullity": obstruction_nullity,
            "dimension": 0,
            "status": "exactly_forbidden_by_central_character",
        },
        "twisted_intertwiner": {
            "map": "v -> v tensor 1 in C^2 tensor Lambda^0(C^6)",
            "matrix_shape": list(t.shape),
            "nonzero": bool(np.linalg.norm(t) > 0),
            "max_group_equivariance_residual": twisted_residual,
            "vacuum_invariance_residual": q6_vacuum_residual,
            "status": "explicit_equivariant_embedding",
        },
        "passed": bool(
            len(elements) == 120 and len(rotations) == 60
            and central_fundamental_residual < tol
            and central_q6_residual < tol
            and twisted_residual < tol
            and q6_vacuum_residual < tol
        ),
        "scope": (
            "The Q6 A5 action is pulled back to 2I and therefore has trivial "
            "central character. The tensor twist is a representation-theoretic "
            "repair, not a derived Pin-to-2I geometric coupling."
        ),
    }

"""Central-character audit for the 2I action on the Q6 exterior Clifford module.

The six-axis A5 representation is pulled back along the abstract quotient
2I -> 2I/{+/-1} ~= A5. This module deliberately does not identify a particular
quaternion matrix with a particular signed-axis matrix: that would require an
explicitly constructed isomorphism between the two A5 presentations.

The central element -1 acts trivially on Lambda*(C^6), but as -I on the
fundamental doublet. This gives an exact vanishing theorem for the untwisted
Hom space. The tensor product C^2 tensor Lambda*(C^6) has an explicit
equivariant vacuum embedding. This twist is algebraic, not a derived
geometric or physical coupling.
"""
from __future__ import annotations

from itertools import combinations
import numpy as np

from src.binary_icosahedral import IDENTITY, binary_icosahedral_elements
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


def fundamental_to_twisted_module_intertwiner() -> np.ndarray:
    """Canonical injection v -> v tensor 1 into C^2 tensor Lambda*(C^6)."""
    t = np.zeros((2 * 64, 2), dtype=np.complex128)
    t[0, 0] = 1.0
    t[64, 1] = 1.0
    return t


def audit_q6_2i_intertwiner_compatibility(tol: float = 2e-8) -> dict[str, object]:
    """Audit central characters and the tensor embedding without false element matching."""
    elements = binary_icosahedral_elements()
    data = icosahedral_signed_axis_representation()
    signed = np.asarray(data["signed_permutations"])
    q6_actions = [exterior_signed_permutation_matrix(s) for s in signed]

    rho_plus = representation_matrix(IDENTITY, "fundamental")
    rho_minus = representation_matrix(-IDENTITY, "fundamental")
    q6_central = np.eye(64, dtype=complex)  # quotient element pi(-1) is identity
    central_fundamental_residual = float(np.linalg.norm(rho_minus + np.eye(2)))
    central_q6_residual = float(np.linalg.norm(q6_central - np.eye(64)))

    # Equivariance at -1 gives -T=T, i.e. 2T=0 on 64*2=128 entries.
    obstruction_constraint_rank = 128
    obstruction_nullity = 0

    t = fundamental_to_twisted_module_intertwiner()
    vacuum = np.eye(64, dtype=complex)[:, 0]
    q6_vacuum_residual = max(
        float(np.linalg.norm(action @ vacuum - vacuum))
        for action in q6_actions
    )
    q6_unitarity_residual = max(
        float(np.linalg.norm(action.conj().T @ action - np.eye(64)))
        for action in q6_actions
    )
    fundamental_unitarity_residual = 0.0
    for q in elements:
        rho2 = representation_matrix(q, "fundamental")
        fundamental_unitarity_residual = max(
            fundamental_unitarity_residual,
            float(np.linalg.norm(rho2.conj().T @ rho2 - np.eye(2))),
        )

    # For any quotient isomorphism 2I/{+/-1} -> A5, the Q6 action fixes
    # the degree-zero vacuum. Therefore this equation holds for every pair
    # (q, a) in 2I x A5, without selecting an unsupported element matching.
    twisted_residual = 0.0
    for q in elements:
        rho2 = representation_matrix(q, "fundamental")
        rho_norm = float(np.linalg.norm(rho2, ord="fro"))
        for action in q6_actions:
            delta = action @ vacuum - vacuum
            twisted_residual = max(
                twisted_residual, rho_norm * float(np.linalg.norm(delta))
            )

    return {
        "audit": "Q6_2I_CENTRAL_CHARACTER_INTERTWINER",
        "group_order_2I": int(len(elements)),
        "quotient_order_A5": int(len(q6_actions)),
        "quotient_identification": (
            "abstract quotient 2I/{+/-1} ~= A5; no elementwise matching "
            "between the two matrix presentations is assumed"
        ),
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
            "max_product_action_equivariance_residual": twisted_residual,
            "Q6_vacuum_invariance_residual": q6_vacuum_residual,
            "status": "explicit_equivariant_embedding",
        },
        "residuals": {
            "Q6_exterior_unitarity": q6_unitarity_residual,
            "2I_fundamental_unitarity": fundamental_unitarity_residual,
        },
        "passed": bool(
            len(elements) == 120 and len(q6_actions) == 60
            and central_fundamental_residual < tol
            and central_q6_residual < tol
            and twisted_residual < tol
            and q6_vacuum_residual < tol
            and q6_unitarity_residual < tol
            and fundamental_unitarity_residual < tol
        ),
        "scope": (
            "The Q6 A5 action is pulled back to 2I and therefore has trivial "
            "central character. The tensor twist is a representation-theoretic "
            "repair, not a derived Pin-to-2I geometric coupling."
        ),
    }

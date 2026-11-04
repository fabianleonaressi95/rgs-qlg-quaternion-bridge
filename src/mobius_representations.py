"""Discrete local-frame representations on the six-face Möbius complex.

This is a Z2 frame-transition system on the dual graph of faces. It is
separate from the independent internal SU(2) lattice gauge field.

For adjacent faces i,j sharing an edge, the transition sign is derived from
their chosen cyclic boundary orientations: t_ij = -s_i s_j, where s is the
orientation of the shared edge relative to the canonical increasing-vertex
orientation. The six dual transitions have non-trivial product -1 around
the dual cycle, encoding the Möbius orientation obstruction.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from src.mobius_complex import build_mobius_complex

Representation = Literal["scalar", "sign", "regular"]
I2 = np.eye(2, dtype=float)
SWAP = np.array([[0.0, 1.0], [1.0, 0.0]])


@dataclass(frozen=True)
class FrameTransition:
    source_face: int
    target_face: int
    shared_edge: tuple[int, int]
    sign: int


def _edge_orientation(face: tuple[int, ...], edge: tuple[int, int]) -> int:
    for u, v in zip(face, face[1:] + face[:1]):
        if (min(u, v), max(u, v)) == edge:
            return 1 if u < v else -1
    raise ValueError("The requested edge is not in the face boundary.")


def build_frame_transitions() -> tuple[FrameTransition, ...]:
    """Derive the six Z2 transitions from shared-edge incidence data."""
    complex_ = build_mobius_complex()
    incidences: dict[tuple[int, int], list[tuple[int, int]]] = {}
    for face_id, face in enumerate(complex_.faces):
        for u, v in zip(face, face[1:] + face[:1]):
            edge = (min(u, v), max(u, v))
            incidences.setdefault(edge, []).append((face_id, 1 if u < v else -1))

    transitions = []
    for edge, incident in sorted(incidences.items()):
        if len(incident) == 2:
            (i, si), (j, sj) = sorted(incident)
            transitions.append(
                FrameTransition(i, j, edge, int(-si * sj))
            )
    return tuple(transitions)


def representation_matrix(kind: Representation, sign: int) -> np.ndarray:
    """Return rho(sign) for the scalar, sign, or regular Z2 representation."""
    if sign not in (-1, 1):
        raise ValueError("A Z2 transition sign must be +1 or -1.")
    if kind == "scalar":
        return np.ones((1, 1), dtype=float)
    if kind == "sign":
        return np.array([[float(sign)]])
    if kind == "regular":
        return I2.copy() if sign == 1 else SWAP.copy()
    raise ValueError(f"Unknown representation: {kind!r}")


def representation_dimension(kind: Representation) -> int:
    if kind in ("scalar", "sign"):
        return 1
    if kind == "regular":
        return 2
    raise ValueError(f"Unknown representation: {kind!r}")


def transition_holonomy(
    transitions: tuple[FrameTransition, ...],
    cycle: tuple[int, ...],
) -> int:
    """Multiply Z2 transitions along an ordered closed dual-face cycle."""
    if len(cycle) < 3 or cycle[0] != cycle[-1]:
        raise ValueError("A dual cycle must be closed and contain at least three steps.")
    lookup = {}
    for item in transitions:
        lookup[(item.source_face, item.target_face)] = item.sign
        lookup[(item.target_face, item.source_face)] = item.sign
    product = 1
    for i, j in zip(cycle[:-1], cycle[1:]):
        try:
            product *= lookup[(i, j)]
        except KeyError as exc:
            raise ValueError(f"Faces {i} and {j} do not share an interior edge.") from exc
    return int(product)


def transport_feature(
    feature: np.ndarray,
    sign: int,
    kind: Representation,
) -> np.ndarray:
    """Transport a feature between adjacent face frames."""
    value = np.asarray(feature, dtype=float)
    dim = representation_dimension(kind)
    if value.shape != (dim,):
        raise ValueError(f"{kind} features must have shape ({dim},).")
    return representation_matrix(kind, sign) @ value


def transform_transition(
    sign: int,
    source_frame_change: int,
    target_frame_change: int,
) -> int:
    """Change local Z2 frames at both endpoints: t' = h_target t h_source."""
    if sign not in (-1, 1) or source_frame_change not in (-1, 1) or target_frame_change not in (-1, 1):
        raise ValueError("All frame transition factors must be Z2 signs.")
    return int(target_frame_change * sign * source_frame_change)


def audit_local_frame_representations() -> dict[str, object]:
    """Audit transition counts, representation laws and Möbius cycle holonomy."""
    transitions = build_frame_transitions()
    cycle = (0, 1, 2, 5, 4, 3, 0)
    holonomy = transition_holonomy(transitions, cycle)
    rep_holonomies = {
        kind: representation_matrix(kind, holonomy).tolist()
        for kind in ("scalar", "sign", "regular")
    }
    representation_law_residual = 0.0
    for kind in ("scalar", "sign", "regular"):
        for a in (-1, 1):
            for b in (-1, 1):
                lhs = representation_matrix(kind, a * b)
                rhs = representation_matrix(kind, a) @ representation_matrix(kind, b)
                representation_law_residual = max(
                    representation_law_residual,
                    float(np.linalg.norm(lhs - rhs, ord="fro")),
                )

    return {
        "audit": "MOBIUS_LOCAL_FRAME_REPRESENTATIONS",
        "dual_faces": 6,
        "shared_edge_transitions": len(transitions),
        "transition_signs": [item.sign for item in transitions],
        "dual_cycle": list(cycle),
        "dual_cycle_holonomy": holonomy,
        "representations": rep_holonomies,
        "representation_homomorphism_max_residual": representation_law_residual,
        "transition_inverse_consistency": all(
            item.sign * item.sign == 1 for item in transitions
        ),
        "orientation_obstruction_detected": holonomy == -1,
    }

def reflection_transition_matrix(
    sign: int, dimension: int = 4, reflected_axis: int = 0
) -> np.ndarray:
    """Choose an explicit O(d) frame transition for the orientation sign.

    sign=-1 is represented by a single coordinate reflection. This is an
    O(d) transition matrix only; no Pin lift or fermionic reflection action
    is defined by this function.
    """
    if sign not in (-1, 1):
        raise ValueError("A frame transition sign must be +1 or -1.")
    if dimension < 2 or not 0 <= reflected_axis < dimension:
        raise ValueError("Invalid reflection dimension or axis.")
    result = np.eye(dimension)
    if sign == -1:
        result[reflected_axis, reflected_axis] = -1.0
    return result


def audit_reflection_transition_representation(
    dimension: int = 4, reflected_axis: int = 0
) -> dict[str, object]:
    """Check O(d) determinant, involution and Z2 composition laws."""
    identity = reflection_transition_matrix(1, dimension, reflected_axis)
    reflection = reflection_transition_matrix(-1, dimension, reflected_axis)
    return {
        "audit": "MOBIUS_ORTHOGONAL_REFLECTION_TRANSITIONS",
        "dimension": dimension,
        "det_identity": float(np.linalg.det(identity)),
        "det_reflection": float(np.linalg.det(reflection)),
        "reflection_involution_fro": float(np.linalg.norm(reflection @ reflection - identity)),
        "reflection_orthogonality_fro": float(np.linalg.norm(reflection.T @ reflection - identity)),
        "is_Od_transition": bool(np.linalg.det(reflection) < 0 and
                                 np.allclose(reflection.T @ reflection, identity)),
        "pin_lift_defined": False,
        "scope": "explicit O(d) frame transition only; no Pin structure inferred",
    }

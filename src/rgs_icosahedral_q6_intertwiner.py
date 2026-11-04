"""Exact finite-group intertwiner between icosahedral axes and Q6 directions.

The 60 orientation-preserving symmetries of the icosahedron act on its six
unoriented vertex axes by signed permutations. The same signed permutations
act on Q6 coordinate directions, with a negative sign represented by bit
complementation. The coordinate identification is explicit and reproducible;
it is an equivariant group-action identification, not a claim that the
icosahedron is metrically embedded in the 6-cube.
"""
from __future__ import annotations
from itertools import permutations
import numpy as np
from src.rgs_weyl_phase_closure import build_icosahedral_c6


def _icosahedral_rotation_group():
    data = build_icosahedral_c6()
    vertices = np.asarray(data["vertices"], dtype=float)
    axes = np.asarray(data["axes"], dtype=float)
    source_inverse = np.linalg.inv(vertices[:3].T)
    rotations = []
    for target_indices in permutations(range(12), 3):
        rotation = vertices[list(target_indices)].T @ source_inverse
        if np.linalg.norm(rotation.T @ rotation - np.eye(3)) > 1e-8:
            continue
        if abs(np.linalg.det(rotation) - 1.0) > 1e-8:
            continue
        image = vertices @ rotation.T
        if max(min(np.linalg.norm(point - vertex) for vertex in vertices)
               for point in image) > 1e-8:
            continue
        if not any(np.linalg.norm(rotation - old) < 1e-8 for old in rotations):
            rotations.append(rotation)
    rotations.sort(key=lambda matrix: tuple(np.round(matrix.ravel(), 12)))
    if len(rotations) != 60:
        raise ArithmeticError(f"expected 60 proper icosahedral rotations, got {len(rotations)}")
    return np.asarray(rotations), axes


def icosahedral_signed_axis_representation():
    rotations, axes = _icosahedral_rotation_group()
    signed = np.zeros((60, 6, 6), dtype=np.int8)
    for g, rotation in enumerate(rotations):
        for source_axis, axis in enumerate(axes):
            image = rotation @ axis
            positive = np.linalg.norm(axes - image, axis=1)
            negative = np.linalg.norm(axes + image, axis=1)
            if positive.min() < 1e-8:
                target, sign = int(positive.argmin()), 1
            elif negative.min() < 1e-8:
                target, sign = int(negative.argmin()), -1
            else:
                raise ArithmeticError("rotation failed to preserve antipodal axes")
            signed[g, target, source_axis] = sign
    return {"rotations": rotations, "axes": axes,
            "signed_permutations": signed, "group_order": 60}


def q6_vertex_action(vertex: int, signed_permutation: np.ndarray) -> int:
    """Apply signed coordinate permutation to big-endian Q6; minus means complement."""
    if not isinstance(vertex, (int, np.integer)) or not 0 <= int(vertex) < 64:
        raise ValueError("vertex must be an integer in 0..63")
    s = np.asarray(signed_permutation)
    if s.shape != (6, 6) or not np.all(np.isin(s, (-1, 0, 1))):
        raise ValueError("signed_permutation must be a 6x6 signed permutation")
    if not (np.all(np.sum(np.abs(s), axis=0) == 1)
            and np.all(np.sum(np.abs(s), axis=1) == 1)):
        raise ValueError("matrix must have one signed unit per row and column")
    source = np.array([(int(vertex) >> (5-j)) & 1 for j in range(6)], dtype=np.int8)
    target = np.zeros(6, dtype=np.int8)
    for k in range(6):
        j = int(np.flatnonzero(s[k])[0])
        target[k] = source[j] ^ int(s[k, j] < 0)
    result = 0
    for bit in target:
        result = (result << 1) | int(bit)
    return result


def audit_icosahedral_q6_intertwiner():
    data = icosahedral_signed_axis_representation()
    rotations = np.asarray(data["rotations"])
    signed = np.asarray(data["signed_permutations"])
    c6 = np.asarray(build_icosahedral_c6()["C6"], dtype=float)
    ortho = max(float(np.linalg.norm(r.T @ r - np.eye(3))) for r in rotations)
    det = max(abs(float(np.linalg.det(r))-1.0) for r in rotations)
    signed_ortho = max(float(np.linalg.norm(s.T @ s - np.eye(6))) for s in signed)
    c6_cov = max(float(np.linalg.norm(s.T @ c6 @ s - c6)) for s in signed)
    edge_violations = 0
    for s in signed:
        images = [q6_vertex_action(v, s) for v in range(64)]
        if len(set(images)) != 64:
            edge_violations += 1
            continue
        for v in range(64):
            for axis in range(6):
                neighbour = v ^ (1 << (5-axis))
                if (images[v] ^ images[neighbour]).bit_count() != 1:
                    edge_violations += 1
    hom = 0.0
    for r1, s1 in zip(rotations, signed):
        for r2, s2 in zip(rotations, signed):
            product = r1 @ r2
            k = int(np.argmin(np.linalg.norm(rotations-product, axis=(1,2))))
            hom = max(hom, float(np.linalg.norm(signed[k]-s1@s2)))
    return {
        "group_order": 60, "q6_vertex_count": 64, "q6_edge_count": 192,
        "orthogonal_residual": ortho, "determinant_residual": det,
        "signed_permutation_residual": signed_ortho,
        "C6_equivariance_residual": c6_cov,
        "group_homomorphism_residual": hom,
        "q6_edge_action_violations": int(edge_violations),
        "all_checks_passed": bool(ortho < 1e-8 and det < 1e-8
            and signed_ortho == 0.0 and c6_cov < 1e-8
            and hom < 1e-8 and edge_violations == 0),
    }

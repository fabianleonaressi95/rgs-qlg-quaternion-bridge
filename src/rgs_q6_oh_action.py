"""Exact O_h action on the recovered 64-point RGS seed and Q6 cells.

The seed is reconstructed by load_rgs_seed_data(), whose coordinate and
coframe hashes are checked against the archived provenance constants.
All 48 signed permutation matrices are then matched to the 64 coordinates.
No nearest-neighbour graph is inferred from distances in this module:
the Q6 graph and square cells are the canonical combinatorial objects.
"""
from __future__ import annotations

from itertools import permutations, product
from pathlib import Path
import argparse
import json

import numpy as np
from scipy.spatial import cKDTree

from src.q6_complex import build_q6_complex
from src.rgs_ellis_refinement import load_rgs_seed_data


def signed_permutation_matrices() -> tuple[np.ndarray, ...]:
    """Return all 48 signed coordinate-permutation matrices in O_h."""
    matrices = []
    for perm in permutations(range(3)):
        for signs in product((-1, 1), repeat=3):
            matrix = np.zeros((3, 3), dtype=np.int64)
            for row, col in enumerate(perm):
                matrix[row, col] = signs[row]
            matrices.append(matrix)
    matrices.sort(key=lambda m: tuple(int(x) for x in m.ravel()))
    return tuple(matrices)


def _matrix_key(matrix: np.ndarray) -> tuple[int, ...]:
    return tuple(int(x) for x in np.asarray(matrix).ravel())


def induced_vertex_permutation(
    coordinates: np.ndarray,
    matrix: np.ndarray,
    tolerance: float = 1e-10,
) -> tuple[np.ndarray, dict]:
    """Match x -> matrix @ x to the unique point in the recovered seed."""
    x = np.asarray(coordinates, dtype=float)
    s = np.asarray(matrix, dtype=float)
    if x.shape != (64, 3) or s.shape != (3, 3):
        raise ValueError("expected coordinates (64,3) and matrix (3,3)")
    if not np.allclose(s.T @ s, np.eye(3), atol=0, rtol=0):
        raise ValueError("matrix must be an exact orthogonal signed permutation")
    tree = cKDTree(x)
    transformed = x @ s.T
    distances, indices = tree.query(transformed, k=2)
    permutation = np.asarray(indices[:, 0], dtype=np.int64)
    max_residual = float(np.max(distances[:, 0]))
    min_match_gap = float(np.min(distances[:, 1] - distances[:, 0]))
    bijective = np.unique(permutation).size == 64
    unique_match = min_match_gap > tolerance
    if max_residual > tolerance or not bijective or not unique_match:
        raise ValueError(
            "O_h action does not induce a unique exact vertex permutation: "
            f"max_residual={max_residual:.3e}, bijective={bijective}, "
            f"minimum_match_gap={min_match_gap:.3e}"
        )
    return permutation, {
        "max_coordinate_residual": max_residual,
        "minimum_match_gap": min_match_gap,
        "bijective": bool(bijective),
    }


def _edge_set(edges: np.ndarray) -> set[tuple[int, int]]:
    return {tuple(sorted((int(u), int(v)))) for u, v in np.asarray(edges)[:, :2]}


def _face_set(faces) -> set[tuple[int, int, int, int]]:
    return {tuple(sorted(map(int, face))) for face in faces}


def _cell_action(permutation: np.ndarray, edges, faces) -> tuple[bool, bool]:
    edge_set = _edge_set(np.asarray(edges))
    face_set = _face_set(faces)
    edge_ok = all(
        tuple(sorted((int(permutation[u]), int(permutation[v])))) in edge_set
        for u, v in edge_set
    )
    face_ok = all(
        tuple(sorted(int(permutation[v]) for v in face)) in face_set
        for face in faces
    )
    return bool(edge_ok), bool(face_ok)


def _orbits(permutations_: tuple[np.ndarray, ...]) -> list[list[int]]:
    unseen = set(range(64))
    output = []
    while unseen:
        seed = min(unseen)
        orbit = sorted({int(p[seed]) for p in permutations_})
        output.append(orbit)
        unseen.difference_update(orbit)
    return sorted(output, key=lambda orbit: (len(orbit), orbit[0]))


def audit_rgs_q6_oh(tolerance: float = 1e-10) -> dict:
    """Audit O_h on X_Oh and test equivariance on Q6 vertices/edges/faces."""
    seed = load_rgs_seed_data()  # also verifies the archived SHA-256 hashes
    x = seed.coordinates
    group = signed_permutation_matrices()
    group_keys = {_matrix_key(g) for g in group}
    if len(group) != 48 or len(group_keys) != 48:
        raise RuntimeError("signed permutation group must have order 48")

    q6 = build_q6_complex()
    edge_array = np.asarray(q6.edges, dtype=np.int64)
    induced = []
    residuals = []
    for matrix in group:
        p, residual = induced_vertex_permutation(x, matrix, tolerance)
        induced.append(p)
        residuals.append(residual)

    # The vertex action must be faithful and preserve group multiplication.
    permutation_keys = {tuple(map(int, p)) for p in induced}
    faithful = len(permutation_keys) == 48
    matrix_to_perm = {
        _matrix_key(g): p for g, p in zip(group, induced)
    }
    homomorphism = True
    for a in group:
        for b in group:
            composed = a @ b
            lhs = matrix_to_perm[_matrix_key(composed)]
            pa = matrix_to_perm[_matrix_key(a)]
            pb = matrix_to_perm[_matrix_key(b)]
            if not np.array_equal(lhs, pa[pb]):
                homomorphism = False
                break
        if not homomorphism:
            break

    edge_face_checks = [
        _cell_action(p, edge_array, q6.faces) for p in induced
    ]
    edge_preserving_count = sum(e for e, _ in edge_face_checks)
    face_preserving_count = sum(f for _, f in edge_face_checks)

    # The central inversion is -I; its induced action is the coordinate
    # antipode and must be an involution with no fixed points.
    inversion = matrix_to_perm[_matrix_key(-np.eye(3, dtype=np.int64))]
    antipodal_involution = np.array_equal(inversion[inversion], np.arange(64))
    antipodal_fixed_points = int(np.count_nonzero(inversion == np.arange(64)))

    orbit_list = _orbits(tuple(induced))
    orbit_sizes = sorted(len(orbit) for orbit in orbit_list)

    return {
        "audit": "RGS_XOH_Q6_INDUCED_OH_ACTION",
        "coordinate_source": "load_rgs_seed_data; archived coordinate/coframe hashes checked",
        "coordinate_shape": list(x.shape),
        "group": {
            "name": "O_h signed permutation group",
            "order": len(group),
            "faithful_on_64_vertices": bool(faithful),
            "action_homomorphism_exact": bool(homomorphism),
            "coordinate_residual_max": max(r["max_coordinate_residual"] for r in residuals),
            "minimum_unique_match_gap": min(r["minimum_match_gap"] for r in residuals),
        },
        "vertex_orbits": orbit_list,
        "vertex_orbit_sizes": orbit_sizes,
        "antipodal": {
            "matrix": (-np.eye(3, dtype=np.int64)).tolist(),
            "involution_exact": bool(antipodal_involution),
            "fixed_point_count": antipodal_fixed_points,
            "coordinate_residual": float(np.max(np.linalg.norm(x[inversion] + x, axis=1))),
        },
        "q6": {
            "vertices": 64,
            "edges": len(q6.edges),
            "square_faces": len(q6.faces),
            "group_elements_preserving_edges": int(edge_preserving_count),
            "group_elements_preserving_faces": int(face_preserving_count),
            "full_Oh_graph_automorphism": bool(edge_preserving_count == 48),
            "full_Oh_cellular_automorphism": bool(
                edge_preserving_count == 48 and face_preserving_count == 48
            ),
        },
        "status": {
            "recovered_coordinate_Oh_action": "NUMERICALLY_VERIFIED",
            "induced_Q6_graph_action": (
                "VERIFIED" if edge_preserving_count == 48 else "NOT_EQUIVARIANT"
            ),
            "induced_Q6_square_complex_action": (
                "VERIFIED" if face_preserving_count == 48 else "NOT_EQUIVARIANT"
            ),
            "SU2_gauge_equivariance": "OPEN",
        },
        "pass": bool(
            len(group) == 48 and faithful and homomorphism
            and antipodal_involution and antipodal_fixed_points == 0
            and edge_preserving_count == 48 and face_preserving_count == 48
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path,
                        default=Path("results/json/rgs_q6_oh_equivariance.json"))
    args = parser.parse_args()
    report = audit_rgs_q6_oh()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if not report["pass"]:
        raise SystemExit("RGS/Q6 O_h equivariance audit failed")


if __name__ == "__main__":
    main()

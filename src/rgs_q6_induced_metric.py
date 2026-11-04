"""Induced RGS metric audit on the canonical Q6 one- and two-skeleton.

The recovered point map is X: {0,1}^6 -> R^3. Its vertexwise
finite-difference Jacobian is 3x6, so its pullback Gram matrix is positive
semidefinite with rank at most three. It is therefore a degenerate
six-dimensional pullback metric, not a nondegenerate 6D Riemannian metric.

This module builds only constructions justified by the data:
* edge lengths and square-face Gram/area diagnostics;
* the vertexwise 3x6 finite-difference map and its Gram matrix;
* 1D P1 edge-element mass/stiffness assembled on the Q6 graph;
* separate O_h-coordinate and Q6-cell equivariance gates.

The assembled graph P1 matrices are not a full 6D Whitney Hodge star,
not a curved continuum limit, and not a Yang-Mills/TEGR construction.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

import numpy as np

from src.q6_complex import build_q6_complex
from src.rgs_ellis_refinement import load_rgs_seed_data
from src.rgs_q6_oh_action import audit_rgs_q6_oh


@dataclass(frozen=True)
class RGSQ6InducedMetric:
    """Finite geometric data in canonical Q6 vertex/edge/face order."""
    coordinates: np.ndarray
    bit_coordinates: np.ndarray
    edges: np.ndarray
    faces: np.ndarray
    edge_lengths: np.ndarray
    edge_mass: np.ndarray
    edge_stiffness: np.ndarray
    local_jacobians: np.ndarray
    local_pullback_grams: np.ndarray
    face_grams: np.ndarray
    face_parallelogram_areas: np.ndarray
    face_triangle_areas: np.ndarray
    face_closure_defects: np.ndarray
    global_affine_fit: np.ndarray
    global_affine_residuals: np.ndarray


def _q6_bits() -> np.ndarray:
    """Canonical big-endian six-bit coordinates, rows indexed 0..63."""
    return np.asarray(
        [[(v >> (5 - axis)) & 1 for axis in range(6)] for v in range(64)],
        dtype=np.float64,
    )


def build_rgs_q6_induced_metric() -> RGSQ6InducedMetric:
    """Construct the finite metric data from hash-verified RGS coordinates."""
    seed = load_rgs_seed_data()
    q6 = build_q6_complex()
    x = np.asarray(seed.coordinates, dtype=np.float64)
    bits = _q6_bits()
    edges = np.asarray(q6.edges, dtype=np.int64)
    faces = np.asarray(q6.faces, dtype=np.int64)

    # The archived edge convention and canonical Q6 edge convention have
    # different axis enumeration order; compare unordered endpoint sets.
    archived = {tuple(sorted(map(int, e[:2]))) for e in seed.q6_edges}
    canonical = {tuple(map(int, e)) for e in edges}
    if archived != canonical:
        raise RuntimeError("recovered RGS edge set does not match canonical Q6")

    vectors = x[edges[:, 1]] - x[edges[:, 0]]
    lengths = np.linalg.norm(vectors, axis=1)
    if np.any(~np.isfinite(lengths)) or np.any(lengths <= 1e-13):
        raise RuntimeError("RGS/Q6 embedding has a zero or non-finite edge")

    # 1D affine P1 element matrices, assembled over all 192 graph edges.
    mass = np.zeros((64, 64), dtype=np.float64)
    stiffness = np.zeros((64, 64), dtype=np.float64)
    for (u, v), length in zip(edges, lengths):
        local_m = (length / 6.0) * np.asarray([[2.0, 1.0], [1.0, 2.0]])
        local_k = (1.0 / length) * np.asarray([[1.0, -1.0], [-1.0, 1.0]])
        ids = np.ix_([int(u), int(v)], [int(u), int(v)])
        mass[ids] += local_m
        stiffness[ids] += local_k

    # Six finite differences at every vertex. Each column is the positively
    # oriented coordinate-axis edge vector, evaluated at the lower endpoint.
    jac = np.empty((64, 3, 6), dtype=np.float64)
    for v in range(64):
        for axis in range(6):
            mask = 1 << (5 - axis)
            low = v & ~mask
            high = low | mask
            jac[v, :, axis] = x[high] - x[low]
    pullback = np.einsum("vki,vkj->vij", jac, jac)

    face_grams = np.empty((len(faces), 2, 2), dtype=np.float64)
    face_para_area = np.empty(len(faces), dtype=np.float64)
    face_tri_area = np.empty(len(faces), dtype=np.float64)
    face_closure = np.empty(len(faces), dtype=np.float64)
    for i, (v, va, vab, vb) in enumerate(faces):
        a = x[va] - x[v]
        b = x[vb] - x[v]
        closure = x[vab] - x[va] - x[vb] + x[v]
        gram = np.asarray([[a @ a, a @ b], [a @ b, b @ b]])
        det = max(float(np.linalg.det(gram)), 0.0)
        face_grams[i] = gram
        face_para_area[i] = np.sqrt(det)
        face_tri_area[i] = 0.5 * (
            np.linalg.norm(np.cross(x[va] - x[v], x[vab] - x[v]))
            + np.linalg.norm(np.cross(x[vab] - x[v], x[vb] - x[v]))
        )
        face_closure[i] = np.linalg.norm(closure)

    # Least-squares affine map X ~= [1, bits] A. Nonzero residual measures
    # failure of the 64-point embedding to be a single affine 6->3 map.
    design = np.column_stack((np.ones(64), bits))
    fit, *_ = np.linalg.lstsq(design, x, rcond=None)
    residuals = np.linalg.norm(design @ fit - x, axis=1)

    return RGSQ6InducedMetric(
        coordinates=x,
        bit_coordinates=bits,
        edges=edges,
        faces=faces,
        edge_lengths=lengths,
        edge_mass=mass,
        edge_stiffness=stiffness,
        local_jacobians=jac,
        local_pullback_grams=pullback,
        face_grams=face_grams,
        face_parallelogram_areas=face_para_area,
        face_triangle_areas=face_tri_area,
        face_closure_defects=face_closure,
        global_affine_fit=fit,
        global_affine_residuals=residuals,
    )


def _numerical_rank(values: np.ndarray, rtol: float = 1e-10) -> int:
    s = np.asarray(values, dtype=np.float64)
    if s.size == 0 or s[0] <= 0.0:
        return 0
    return int(np.count_nonzero(s > rtol * s[0]))


def audit_rgs_q6_induced_metric(
    model: RGSQ6InducedMetric | None = None,
    atol: float = 1e-10,
) -> dict:
    """Return numerical diagnostics with explicit degeneracy/equivariance gates."""
    data = build_rgs_q6_induced_metric() if model is None else model
    mass_eigs = np.linalg.eigvalsh(data.edge_mass)
    stiffness_eigs = np.linalg.eigvalsh(data.edge_stiffness)
    local_ranks = []
    local_min_eigs = []
    local_max_eigs = []
    for gram in data.local_pullback_grams:
        eig = np.linalg.eigvalsh(gram)
        local_min_eigs.append(float(eig[0]))
        local_max_eigs.append(float(eig[-1]))
        local_ranks.append(_numerical_rank(eig[::-1]))

    face_ranks = [
        _numerical_rank(np.linalg.eigvalsh(g)[::-1]) for g in data.face_grams
    ]
    oh = audit_rgs_q6_oh()
    constant = np.ones(64)

    checks = {
        "edge_count_192": len(data.edges) == 192,
        "face_count_240": len(data.faces) == 240,
        "all_edge_lengths_positive": bool(np.all(data.edge_lengths > 0.0)),
        "edge_mass_positive_definite": bool(mass_eigs[0] > 0.0),
        "edge_stiffness_constant_kernel": bool(
            np.linalg.norm(data.edge_stiffness @ constant) <= atol
        ),
        "edge_stiffness_psd": bool(stiffness_eigs[0] >= -atol),
        "local_pullback_rank_at_most_3": bool(max(local_ranks) <= 3),
        "face_gram_rank_at_most_2": bool(max(face_ranks) <= 2),
        "coordinate_Oh_action_verified": bool(oh["geometry_pass"]),
    }
    return {
        "audit": "RGS_Q6_INDUCED_METRIC_V1",
        "source": "load_rgs_seed_data() with archived coordinate, edge and coframe SHA-256 checks",
        "conventions": {
            "vertex_order": "integer labels 0..63; six-bit big-endian",
            "edge_order": "canonical lexicographic (u,v), u<v",
            "face_order": "axis pairs a<b, base vertex ascending",
            "edge_element": "1D affine P1, M=L/6[[2,1],[1,2]], K=1/L[[1,-1],[-1,1]]",
            "local_map": "six finite differences define J_v in R^(3x6); G_v=J_v.T@J_v",
            "face_metric": "Gram matrix of the two edges from the canonical face base vertex",
        },
        "counts": {"vertices": 64, "edges": len(data.edges), "square_faces": len(data.faces)},
        "edge_lengths": {
            "minimum": float(data.edge_lengths.min()),
            "maximum": float(data.edge_lengths.max()),
            "mean": float(data.edge_lengths.mean()),
            "standard_deviation": float(data.edge_lengths.std()),
        },
        "face_geometry": {
            "parallelogram_area_min": float(data.face_parallelogram_areas.min()),
            "parallelogram_area_max": float(data.face_parallelogram_areas.max()),
            "triangulated_quadrilateral_area_min": float(data.face_triangle_areas.min()),
            "triangulated_quadrilateral_area_max": float(data.face_triangle_areas.max()),
            "closure_defect_max": float(data.face_closure_defects.max()),
            "rank_min": int(min(face_ranks)),
            "rank_max": int(max(face_ranks)),
        },
        "pullback_6d": {
            "jacobian_shape": list(data.local_jacobians.shape),
            "gram_shape": list(data.local_pullback_grams.shape),
            "rank_min": int(min(local_ranks)),
            "rank_max": int(max(local_ranks)),
            "eigenvalue_min": float(min(local_min_eigs)),
            "eigenvalue_max": float(max(local_max_eigs)),
            "interpretation": "DEGENERATE_SEMIMETRIC: rank <= 3 because the target is R^3; not a positive-definite 6D metric",
        },
        "global_affine_fit": {
            "max_residual": float(data.global_affine_residuals.max()),
            "rms_residual": float(np.sqrt(np.mean(data.global_affine_residuals**2))),
            "interpretation": "single affine 6-to-3 map diagnostic; nonzero residual means the embedding is not globally affine",
        },
        "graph_p1": {
            "mass_min_eigenvalue": float(mass_eigs[0]),
            "mass_condition_number": float(np.linalg.cond(data.edge_mass)),
            "stiffness_min_eigenvalue": float(stiffness_eigs[0]),
            "stiffness_max_eigenvalue": float(stiffness_eigs[-1]),
            "stiffness_constant_mode_residual": float(np.linalg.norm(data.edge_stiffness @ constant)),
            "scope": "metric-weighted 1D P1 finite elements on the Q6 graph; not a 6D volume FEM or Whitney 1-form Hodge star",
        },
        "symmetry": {
            "coordinate_Oh_order": oh["group"]["order"],
            "coordinate_Oh_action": oh["status"]["recovered_coordinate_Oh_action"],
            "q6_edge_preserving_elements": oh["q6"]["group_elements_preserving_edges"],
            "q6_face_preserving_elements": oh["q6"]["group_elements_preserving_faces"],
            "q6_full_cellular_equivariance": oh["q6"]["full_Oh_cellular_automorphism"],
            "q6_edge_equivariance_status": oh["status"]["induced_Q6_graph_action"],
            "rule": "O_h invariance of the point cloud is not promoted to Q6 metric equivariance unless the induced vertex permutations preserve the Q6 cells",
        },
        "checks": checks,
        "pass": bool(all(checks.values())),
        "evidence": {
            "coordinate_seed": "NUMERICALLY_VERIFIED_BY_HASHED_RECONSTRUCTION",
            "edge_and_face_metric": "NUMERICALLY_VERIFIED_FINITE_COMPLEX",
            "six_dimensional_pullback": "DEGENERATE_BY_RANK_BOUND",
            "graph_p1": "NUMERICALLY_VERIFIED_FINITE_ELEMENT",
            "continuum_hodge_mosco": "OPEN",
            "curved_RGS_whitney_hodge": "OPEN",
            "TEGR_and_Yang_Mills": "OPEN",
        },
    }


def export_rgs_q6_induced_metric(path: str) -> dict:
    """Write the finite matrices and a compact audit report to NPZ."""
    import hashlib
    import json
    from pathlib import Path

    data = build_rgs_q6_induced_metric()
    report = audit_rgs_q6_induced_metric(data)
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output,
        coordinates=data.coordinates,
        bit_coordinates=data.bit_coordinates,
        edges=data.edges,
        faces=data.faces,
        edge_lengths=data.edge_lengths,
        edge_mass=data.edge_mass,
        edge_stiffness=data.edge_stiffness,
        local_jacobians=data.local_jacobians,
        local_pullback_grams=data.local_pullback_grams,
        face_grams=data.face_grams,
        face_parallelogram_areas=data.face_parallelogram_areas,
        face_triangle_areas=data.face_triangle_areas,
        face_closure_defects=data.face_closure_defects,
        global_affine_fit=data.global_affine_fit,
        global_affine_residuals=data.global_affine_residuals,
        report_json=np.asarray(json.dumps(report, sort_keys=True)),
    )
    digest = hashlib.sha256(output.read_bytes()).hexdigest()
    return {"path": str(output), "sha256": digest, "report": report}


def main() -> None:
    import argparse
    import json
    from pathlib import Path

    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("results/rgs_q6_induced_metric"))
    args = parser.parse_args()
    result = export_rgs_q6_induced_metric(str(args.output.with_suffix(".npz")))
    report_path = args.output.with_suffix(".json")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(result["report"], indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

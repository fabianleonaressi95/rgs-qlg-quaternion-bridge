"""Golden 6D->3D cut-and-project construction for RGS–PSU–Q6.

The projection convention is the recovered icosahedral six-axis convention.
The internal acceptance window is the projection of the centred unit 6-cube;
its convex hull is the rhombic triacontahedral zonotope (32 vertices,
30 rhombic facets). A finite physical-radius cut produces a finite patch of
the regular model set.

The module also implements the previously specified torsion-coupled internal
time-axis criterion. A real (64,4) torsion trace must be supplied; the Q6
edge coframe is not silently relabelled as torsion.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from pathlib import Path
import argparse
import hashlib
import json

import numpy as np
from scipy.spatial import ConvexHull
from scipy.spatial.distance import pdist


PHI = (1.0 + np.sqrt(5.0)) / 2.0
PHI_STAR = (1.0 - np.sqrt(5.0)) / 2.0


@dataclass(frozen=True)
class CutProjectData:
    p_parallel: np.ndarray
    p_internal: np.ndarray
    orthogonal_frame: np.ndarray
    q6_bits_centered: np.ndarray
    q6_parallel: np.ndarray
    q6_internal: np.ndarray
    window_vertices: np.ndarray
    window_hull_equations: np.ndarray
    window_vertex_indices: np.ndarray
    window_facet_count: int


@dataclass(frozen=True)
class ModelSetPatch:
    lattice_points: np.ndarray
    physical_points: np.ndarray
    internal_points: np.ndarray
    boundary_margins: np.ndarray
    physical_radius: float
    integer_search_bound: int


def icosahedral_axis_matrix(parameter: float) -> np.ndarray:
    """Recovered 3x6 six-axis icosahedral matrix A(parameter)."""
    p = float(parameter)
    return np.asarray([
        [0.0, 0.0, 1.0, -1.0, p, p],
        [1.0, -1.0, p, p, 0.0, 0.0],
        [p, p, 0.0, 0.0, 1.0, -1.0],
    ], dtype=np.float64)


def golden_cut_project_frame() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return orthonormal physical/internal 3x6 projections and 6x6 frame.

    Uses phi* = (1-sqrt(5))/2 = -1/phi. The row spaces are exactly
    orthogonal in algebra and orthonormal after the displayed normalization.
    """
    a_parallel = icosahedral_axis_matrix(PHI)
    a_internal = icosahedral_axis_matrix(PHI_STAR)
    p_parallel = a_parallel / np.sqrt(2.0 * (2.0 + PHI))
    p_internal = a_internal / np.sqrt(2.0 * (2.0 + PHI_STAR))
    frame = np.vstack((p_parallel, p_internal))
    return p_parallel, p_internal, frame


def _cube_vertices_centered() -> np.ndarray:
    return np.asarray([
        [((v >> (5 - j)) & 1) - 0.5 for j in range(6)]
        for v in range(64)
    ], dtype=np.float64)


def _unique_facet_planes(equations: np.ndarray, decimals: int = 9) -> np.ndarray:
    """Deduplicate triangulated coplanar hull facets by oriented plane."""
    unique: dict[tuple[float, ...], np.ndarray] = {}
    for equation in np.asarray(equations):
        key = tuple(np.round(equation, decimals))
        unique[key] = equation.copy()
    return np.asarray(list(unique.values()), dtype=np.float64)


def build_cut_project_data() -> CutProjectData:
    """Build the icosahedral projections, Q6 seed and canonical window."""
    p_parallel, p_internal, frame = golden_cut_project_frame()
    cube = _cube_vertices_centered()
    q6_parallel = cube @ p_parallel.T
    q6_internal = cube @ p_internal.T

    hull = ConvexHull(q6_internal)
    window_vertices = q6_internal[np.asarray(hull.vertices, dtype=np.int64)]
    facet_planes = _unique_facet_planes(hull.equations)
    return CutProjectData(
        p_parallel=p_parallel,
        p_internal=p_internal,
        orthogonal_frame=frame,
        q6_bits_centered=cube,
        q6_parallel=q6_parallel,
        q6_internal=q6_internal,
        window_vertices=window_vertices,
        window_hull_equations=facet_planes,
        window_vertex_indices=np.asarray(hull.vertices, dtype=np.int64),
        window_facet_count=int(len(facet_planes)),
    )


def audit_cut_project_frame(data: CutProjectData | None = None) -> dict:
    """Audit golden conjugacy, orthogonal decomposition and window geometry."""
    d = build_cut_project_data() if data is None else data
    identity = np.eye(6)
    parallel_gram = d.p_parallel @ d.p_parallel.T
    internal_gram = d.p_internal @ d.p_internal.T
    cross_gram = d.p_parallel @ d.p_internal.T
    frame_gram = d.orthogonal_frame @ d.orthogonal_frame.T
    hull = ConvexHull(d.q6_internal)
    unique_planes = _unique_facet_planes(hull.equations)
    distances = pdist(d.q6_parallel)
    return {
        "audit": "RGS_PSU_GOLDEN_CUT_PROJECT_V1",
        "constants": {
            "phi": float(PHI),
            "phi_star": float(PHI_STAR),
            "phi_phi_star_plus_1": float(PHI * PHI_STAR + 1.0),
        },
        "projection": {
            "P_parallel_shape": list(d.p_parallel.shape),
            "P_internal_shape": list(d.p_internal.shape),
            "parallel_orthonormal_residual": float(np.linalg.norm(parallel_gram - np.eye(3))),
            "internal_orthonormal_residual": float(np.linalg.norm(internal_gram - np.eye(3))),
            "cross_orthogonality_residual": float(np.linalg.norm(cross_gram)),
            "full_O6_frame_residual": float(np.linalg.norm(frame_gram - identity)),
            "rank_parallel": int(np.linalg.matrix_rank(d.p_parallel)),
            "rank_internal": int(np.linalg.matrix_rank(d.p_internal)),
            "exact_decomposition_status": "NUMERICALLY_VERIFIED_ORTHOGONAL_3_PLUS_3_SPLIT",
        },
        "q6_seed": {
            "vertices": int(len(d.q6_bits_centered)),
            "distinct_physical_points": int(len(np.unique(np.round(d.q6_parallel, 11), axis=0))),
            "physical_radius_min": float(np.linalg.norm(d.q6_parallel, axis=1).min()),
            "physical_radius_max": float(np.linalg.norm(d.q6_parallel, axis=1).max()),
            "pair_distance_min": float(distances.min()),
            "pair_distance_max": float(distances.max()),
            "all_internal_points_inside_window": bool(
                np.max(d.window_hull_equations[:, :3] @ d.q6_internal.T
                       + d.window_hull_equations[:, 3:4]) <= 1e-9
            ),
        },
        "acceptance_window": {
            "construction": "P_internal([-1/2,1/2]^6), the projected centred 6-cube zonotope",
            "shape": "rhombic_triacontahedron",
            "vertices": int(len(d.window_vertices)),
            "facets": int(len(unique_planes)),
            "triangulated_hull_facets": int(len(hull.simplices)),
            "boundary_measure": "zero (finite union of planar facets)",
            "status": "CANONICAL_CUBIC_WINDOW_FOR_THIS_PROJECTION_CONVENTION",
        },
        "checks": {
            "parallel_orthonormal": bool(np.linalg.norm(parallel_gram - np.eye(3)) < 1e-12),
            "internal_orthonormal": bool(np.linalg.norm(internal_gram - np.eye(3)) < 1e-12),
            "parallel_internal_orthogonal": bool(np.linalg.norm(cross_gram) < 1e-12),
            "full_frame_orthogonal": bool(np.linalg.norm(frame_gram - identity) < 1e-12),
            "window_is_32_vertex_30_facet_triacontahedron": bool(
                len(d.window_vertices) == 32 and len(unique_planes) == 30
            ),
            "q6_64_points_inside_window": bool(
                np.max(d.window_hull_equations[:, :3] @ d.q6_internal.T
                       + d.window_hull_equations[:, 3:4]) <= 1e-9
            ),
        },
    } | {"pass": bool(
        np.linalg.norm(parallel_gram - np.eye(3)) < 1e-12
        and np.linalg.norm(internal_gram - np.eye(3)) < 1e-12
        and np.linalg.norm(cross_gram) < 1e-12
        and np.linalg.norm(frame_gram - identity) < 1e-12
        and len(d.window_vertices) == 32
        and len(unique_planes) == 30
    )}


def _window_margins(points: np.ndarray, equations: np.ndarray) -> np.ndarray:
    values = np.asarray(points) @ equations[:, :3].T + equations[:, 3]
    return -np.max(values, axis=1)


def enumerate_model_set_patch(
    physical_radius: float = 2.5,
    data: CutProjectData | None = None,
    window_scale: float = 1.0,
    boundary_tolerance: float = 1e-10,
) -> ModelSetPatch:
    """Enumerate a finite physical-radius patch of the regular model set.

    The integer search bound is certified from orthogonality:
      ||n||^2 = ||pi_parallel(n)||^2 + ||pi_internal(n)||^2.
    Thus |pi_parallel(n)|<=R and pi_internal(n) in W imply
      |n_i| <= sqrt(R^2 + max_{y in W}||y||^2).
    """
    if physical_radius <= 0.0 or window_scale <= 0.0:
        raise ValueError("physical_radius and window_scale must be positive")
    d = build_cut_project_data() if data is None else data
    w_vertices = window_scale * d.window_vertices
    w_hull = ConvexHull(w_vertices)
    equations = _unique_facet_planes(w_hull.equations)
    rho = float(np.max(np.linalg.norm(w_vertices, axis=1)))
    bound = int(np.ceil(np.sqrt(physical_radius**2 + rho**2))) + 1

    lattice = np.asarray(list(product(range(-bound, bound + 1), repeat=6)), dtype=np.int64)
    physical = lattice @ d.p_parallel.T
    internal = lattice @ d.p_internal.T
    radial = np.linalg.norm(physical, axis=1)
    margins = _window_margins(internal, equations)
    keep = (radial <= physical_radius + boundary_tolerance) & (margins >= -boundary_tolerance)

    lattice = lattice[keep]
    physical = physical[keep]
    internal = internal[keep]
    margins = margins[keep]
    order = np.lexsort(tuple(lattice[:, k] for k in range(5, -1, -1)))
    return ModelSetPatch(
        lattice_points=lattice[order],
        physical_points=physical[order],
        internal_points=internal[order],
        boundary_margins=margins[order],
        physical_radius=float(physical_radius),
        integer_search_bound=bound,
    )


def audit_model_set_patch(patch: ModelSetPatch) -> dict:
    """Report finite patch counts, separation, window boundary hits and radii."""
    x = patch.physical_points
    distances = pdist(x) if len(x) > 1 else np.asarray([], dtype=float)
    duplicate_count = 0
    if len(x) > 1:
        duplicate_count = int(len(x) - len(np.unique(np.round(x, 10), axis=0)))
    return {
        "physical_radius": patch.physical_radius,
        "integer_search_bound": patch.integer_search_bound,
        "lattice_points": int(len(patch.lattice_points)),
        "physical_radius_max_observed": float(np.linalg.norm(x, axis=1).max()) if len(x) else 0.0,
        "minimum_physical_separation": float(distances.min()) if len(distances) else None,
        "maximum_physical_separation": float(distances.max()) if len(distances) else None,
        "rounded_projection_duplicates": duplicate_count,
        "window_boundary_hits": int(np.count_nonzero(
            np.abs(patch.boundary_margins) <= 1e-9
        )),
        "window_margin_min": float(patch.boundary_margins.min()) if len(x) else None,
        "window_margin_max": float(patch.boundary_margins.max()) if len(x) else None,
        "status": "FINITE_PATCH_OF_CUT_AND_PROJECT_MODEL_SET",
    }


def compare_rgs_seed_to_cut_project(
    rgs_coordinates: np.ndarray,
    data: CutProjectData | None = None,
) -> dict:
    """Compare pair-distance invariants; does not assume matching vertex labels."""
    d = build_cut_project_data() if data is None else data
    x = np.asarray(rgs_coordinates, dtype=np.float64)
    if x.shape != (64, 3):
        raise ValueError("rgs_coordinates must have shape (64,3)")
    dx = np.sort(pdist(d.q6_parallel))
    dr = np.sort(pdist(x))
    dx /= np.sqrt(np.mean(dx**2))
    dr /= np.sqrt(np.mean(dr**2))
    residual = float(np.linalg.norm(dx - dr) / np.sqrt(len(dx)))
    return {
        "comparison": "UNLABELLED_PAIR_DISTANCE_SPECTRUM",
        "rgs_points": 64,
        "cut_project_q6_points": 64,
        "normalized_pair_distance_rms_residual": residual,
        "similarity_equivalent_by_pair_distance_test": bool(residual < 1e-10),
        "interpretation": (
            "A nonzero residual rejects equality up to a common scale and rigid motion; "
            "a zero residual is necessary but not sufficient for point-cloud congruence."
        ),
    }


def _standardize_columns(array: np.ndarray) -> np.ndarray:
    z = np.asarray(array, dtype=np.float64)
    z = z - z.mean(axis=0, keepdims=True)
    norms = np.linalg.norm(z, axis=0)
    if np.any(norms <= 1e-14):
        raise ValueError("torsion/internal trace has a constant or zero-variance column")
    return z / norms


def _fix_svd_sign(u: np.ndarray, vh: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    pivot = int(np.argmax(np.abs(u)))
    if u[pivot] < 0.0:
        return -u, -vh
    return u, vh


def audit_torsion_coupled_internal_time_axis(
    internal_coordinates: np.ndarray,
    torsion_trace: np.ndarray,
    n_permutations: int = 5000,
    n_bootstrap: int = 2000,
    seed: int = 20261021,
) -> dict:
    """Select a candidate internal time axis by maximal torsion cross-coupling.

    Inputs are 64x3 internal coordinates and a separately derived 64xm
    torsion trace. Columns are centered and L2-normalized before the cross
    matrix C=Yz.T@Tz/64 is decomposed. The leading left singular vector is
    u_tau; its paired right singular vector is the torsion mode. A permutation
    null and row bootstrap quantify association and directional stability.
    """
    y = np.asarray(internal_coordinates, dtype=np.float64)
    t = np.asarray(torsion_trace, dtype=np.float64)
    if y.ndim != 2 or y.shape != (64, 3):
        raise ValueError("internal_coordinates must have shape (64,3)")
    if t.ndim != 2 or t.shape[0] != 64 or t.shape[1] < 1:
        raise ValueError("torsion_trace must have shape (64,m), m>=1")
    if not np.all(np.isfinite(y)) or not np.all(np.isfinite(t)):
        raise ValueError("internal coordinates and torsion trace must be finite")
    yz = _standardize_columns(y)
    tz = _standardize_columns(t)
    cross = yz.T @ tz / 64.0
    u, singular, vh = np.linalg.svd(cross, full_matrices=False)
    u0, v0 = _fix_svd_sign(u[:, 0].copy(), vh[0].copy())
    tau = (y - y.mean(axis=0, keepdims=True)) @ u0

    rng = np.random.default_rng(seed)
    null_max = np.empty(n_permutations, dtype=np.float64)
    for k in range(n_permutations):
        permuted = tz[rng.permutation(64)]
        null_max[k] = np.linalg.svd(yz.T @ permuted / 64.0, compute_uv=False)[0]
    p_value = float((1 + np.count_nonzero(null_max >= singular[0])) / (n_permutations + 1))

    angles = np.empty(n_bootstrap, dtype=np.float64)
    for k in range(n_bootstrap):
        indices = rng.integers(0, 64, size=64)
        yb = y[indices]
        tb = t[indices]
        try:
            ybz = _standardize_columns(yb)
            tbz = _standardize_columns(tb)
            ub = np.linalg.svd(ybz.T @ tbz / 64.0, full_matrices=False)[0][:, 0]
            cosine = min(1.0, max(-1.0, abs(float(np.dot(ub, u0)))))
            angles[k] = np.arccos(cosine)
        except (ValueError, np.linalg.LinAlgError):
            angles[k] = np.nan
    finite_angles = angles[np.isfinite(angles)]
    if finite_angles.size == 0:
        raise RuntimeError("all bootstrap samples were singular")
    return {
        "criterion": "leading singular mode of centered, column-normalized Y_internal.T @ torsion_trace / 64",
        "internal_time_axis_u_tau": u0.tolist(),
        "torsion_coupled_mode_v": v0.tolist(),
        "candidate_tau_coordinate_64": tau.tolist(),
        "cross_coupling_matrix": cross.tolist(),
        "singular_values": singular.tolist(),
        "leading_coupling_strength": float(singular[0]),
        "permutation_null": {
            "n": int(n_permutations),
            "seed": int(seed),
            "null_mean": float(null_max.mean()) if n_permutations else None,
            "null_95_percentile": float(np.quantile(null_max, 0.95)) if n_permutations else None,
            "p_value": p_value if n_permutations else None,
        },
        "bootstrap": {
            "n_requested": int(n_bootstrap),
            "n_valid": int(finite_angles.size),
            "median_axis_angle_radians": float(np.median(finite_angles)),
            "p95_axis_angle_radians": float(np.quantile(finite_angles, 0.95)),
            "median_axis_angle_degrees": float(np.degrees(np.median(finite_angles))),
            "p95_axis_angle_degrees": float(np.degrees(np.quantile(finite_angles, 0.95))),
        },
        "epistemic_status": (
            "MODEL_CANDIDATE_ONLY: a stable, significant internal direction coupled to the supplied "
            "torsion observable; it does not establish physical time or Lorentzian dynamics."
        ),
    }


def export_cut_project_audit(
    output_prefix: str | Path,
    physical_radius: float = 2.5,
    rgs_coordinates: np.ndarray | None = None,
    torsion_trace: np.ndarray | None = None,
    n_permutations: int = 5000,
    n_bootstrap: int = 2000,
) -> dict:
    """Export matrices, Q6 projections, model-set patch and JSON audit."""
    prefix = Path(output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    data = build_cut_project_data()
    frame_report = audit_cut_project_frame(data)
    patch = enumerate_model_set_patch(physical_radius, data)
    report = {
        "frame_and_window": frame_report,
        "finite_model_set_patch": audit_model_set_patch(patch),
        "RGS_seed_comparison": (
            compare_rgs_seed_to_cut_project(rgs_coordinates, data)
            if rgs_coordinates is not None else
            {"status": "NOT_RUN_NO_RGS_COORDINATES_SUPPLIED"}
        ),
        "temporal_axis": (
            audit_torsion_coupled_internal_time_axis(
                data.q6_internal, torsion_trace, n_permutations, n_bootstrap
            ) if torsion_trace is not None else {
                "status": "NOT_RUN",
                "required_input": "independently derived torsion_trace with shape (64,m); do not substitute the 192x4 edge coframe",
            }
        ),
        "evidence_ledger": {
            "golden_orthogonal_3_plus_3_split": "NUMERICALLY_VERIFIED_FROM_EXPLICIT_MATRICES",
            "rhombic_triacontahedral_window": "NUMERICALLY_VERIFIED_AS_PROJECTED_CUBE_ZONOTOPE",
            "finite_cut_and_project_patch": "NUMERICALLY_GENERATED_WITH_CERTIFIED_INTEGER_SEARCH_BOUND",
            "RGS_64_point_cloud_identification": "COMPARISON_ONLY_NOT_ASSUMED",
            "torsion_selected_internal_time_axis": (
                "CANDIDATE_IF_TORSION_TRACE_SUPPLIED_AND_NULL_BOOTSTRAP_AUDIT_PASSES"
                if torsion_trace is not None else "OPEN_MISSING_INDEPENDENT_64xm_TORSION_TRACE"
            ),
            "4D_Lorentzian_metric_and_TEGR": "OPEN",
            "continuum_quasicrystal_limit_and_physical_identification": "OPEN",
        },
    }
    npz_path = prefix.with_suffix(".npz")
    json_path = prefix.with_suffix(".json")
    np.savez_compressed(
        npz_path,
        P_parallel=data.p_parallel,
        P_internal=data.p_internal,
        orthogonal_frame=data.orthogonal_frame,
        Q6_bits_centered=data.q6_bits_centered,
        Q6_parallel=data.q6_parallel,
        Q6_internal=data.q6_internal,
        window_vertices=data.window_vertices,
        window_hull_equations=data.window_hull_equations,
        model_set_lattice_points=patch.lattice_points,
        model_set_physical_points=patch.physical_points,
        model_set_internal_points=patch.internal_points,
        model_set_boundary_margins=patch.boundary_margins,
    )
    digest = hashlib.sha256(npz_path.read_bytes()).hexdigest()
    report["archive"] = {"npz": str(npz_path), "sha256": digest}
    json_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("results/rgs_psu_cut_project"))
    parser.add_argument("--physical-radius", type=float, default=2.5)
    parser.add_argument("--rgs-coordinates", type=Path, default=None)
    parser.add_argument("--torsion-trace", type=Path, default=None)
    parser.add_argument("--permutations", type=int, default=5000)
    parser.add_argument("--bootstrap", type=int, default=2000)
    args = parser.parse_args()
    x = np.load(args.rgs_coordinates) if args.rgs_coordinates else None
    torsion = np.load(args.torsion_trace) if args.torsion_trace else None
    report = export_cut_project_audit(
        args.output, args.physical_radius, x, torsion,
        args.permutations, args.bootstrap,
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

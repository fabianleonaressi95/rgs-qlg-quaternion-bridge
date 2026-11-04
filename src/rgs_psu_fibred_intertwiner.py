"""Exact finite Q6 = base x sign-fibre map into the archived RGS 64B seed.

The 64 Q6 labels split canonically as (b,s) in Z2^3 x Z2^3:
the first three most-significant bits select one of eight RGS orbit
representatives p_b; the last three bits select coordinate signs.
The RGS map is Phi(b,s) = (2s-1) componentwise-multiplied by p_b.

This is a finite nonlinear/fibred correspondence, not an affine map and
not a continuum bundle construction. Pulling the RGS Euclidean edge lengths
back to Q6 gives a metric-weighted 1D P1 graph FEM action.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from .rgs_ellis_refinement import load_rgs_seed_data
from .rgs_psu_cut_project import build_cut_project_data


@dataclass(frozen=True)
class FibredIntertwinerData:
    base_labels: np.ndarray
    fibre_bits: np.ndarray
    orbit_representatives: np.ndarray
    signs: np.ndarray
    rgs_coordinates: np.ndarray
    q6_edges: np.ndarray
    edge_lengths_rgs: np.ndarray
    edge_lengths_physical_projection: np.ndarray
    edge_base_mask: np.ndarray
    edge_fibre_mask: np.ndarray
    incidence: np.ndarray
    stiffness_base: np.ndarray
    stiffness_fibre: np.ndarray
    stiffness: np.ndarray
    mass: np.ndarray
    generalized_eigenvalues: np.ndarray


def _edge_fem_matrices(
    edges: np.ndarray, lengths: np.ndarray, n_vertices: int = 64
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return incidence, P1 stiffness and consistent P1 mass on weighted edges."""
    if edges.ndim != 2 or edges.shape[1] != 3 or lengths.shape != (len(edges),):
        raise ValueError("expected an (n_edges,3) Q6 edge array and one length per edge")
    if np.any(~np.isfinite(lengths)) or np.any(lengths <= 0.0):
        raise ValueError("edge lengths must be finite and strictly positive")
    incidence = np.zeros((n_vertices, len(edges)), dtype=np.float64)
    stiffness = np.zeros((n_vertices, n_vertices), dtype=np.float64)
    mass = np.zeros((n_vertices, n_vertices), dtype=np.float64)
    local_stiffness = np.asarray([[1.0, -1.0], [-1.0, 1.0]])
    local_mass = np.asarray([[2.0, 1.0], [1.0, 2.0]])
    for k, ((u, v, _axis), length) in enumerate(zip(edges, lengths)):
        u, v = int(u), int(v)
        incidence[u, k] = -1.0
        incidence[v, k] = 1.0
        ids = np.ix_([u, v], [u, v])
        stiffness[ids] += local_stiffness / length
        mass[ids] += (length / 6.0) * local_mass
    return incidence, stiffness, mass


def build_fibred_intertwiner_data() -> FibredIntertwinerData:
    """Construct the exact labelled Q6-to-RGS map and induced graph FEM."""
    seed = load_rgs_seed_data()
    cut = build_cut_project_data()

    # The cut-project cube uses the same integer label g=0,...,63:
    # axes 0:3 are the base bits, axes 3:6 are the sign-fibre bits.
    bits = cut.q6_bits_centered + 0.5
    if not np.allclose(bits, np.rint(bits), atol=0.0):
        raise RuntimeError("cut-project cube coordinates are not binary")
    bits = np.rint(bits).astype(np.int64)
    base_bits = bits[:, :3]
    fibre_bits = bits[:, 3:]
    base_labels = base_bits @ np.asarray([4, 2, 1], dtype=np.int64)
    signs = 2 * fibre_bits - 1

    # load_rgs_seed_data orders 8 orbit representatives, each followed by
    # itertools.product((-1,+1), repeat=3). The first point in each block
    # therefore has sign (-1,-1,-1); absolute values recover p_b exactly.
    orbit_representatives = np.abs(seed.coordinates.reshape(8, 8, 3)[:, 0, :])
    mapped = signs * orbit_representatives[base_labels]
    if not np.array_equal(mapped, seed.coordinates):
        residual = float(np.max(np.abs(mapped - seed.coordinates)))
        raise RuntimeError(f"base/sign factorization failed to reproduce RGS seed: {residual}")

    # Invert the golden orthogonal 3+3 frame: [x_parallel,x_internal] -> cube.
    lifted = np.column_stack((cut.q6_parallel, cut.q6_internal))
    recovered_cube = lifted @ cut.orthogonal_frame
    lift_residual = float(np.max(np.abs(recovered_cube - cut.q6_bits_centered)))
    if lift_residual > 1e-12:
        raise RuntimeError("golden 3+3 frame failed to reconstruct Q6 cube labels")

    edges = seed.q6_edges.astype(np.int64)
    edge_lengths = np.linalg.norm(mapped[edges[:, 1]] - mapped[edges[:, 0]], axis=1)
    physical_lengths = np.linalg.norm(
        cut.q6_parallel[edges[:, 1]] - cut.q6_parallel[edges[:, 0]], axis=1
    )
    # q6_edges uses mask 1<<axis: axes 0,1,2 are the three low/sign bits;
    # axes 3,4,5 are the three high/base bits.
    edge_fibre_mask = edges[:, 2] < 3
    edge_base_mask = ~edge_fibre_mask

    incidence, stiffness, mass = _edge_fem_matrices(edges, edge_lengths)
    _, stiffness_base, _ = _edge_fem_matrices(
        edges[edge_base_mask], edge_lengths[edge_base_mask]
    )
    _, stiffness_fibre, _ = _edge_fem_matrices(
        edges[edge_fibre_mask], edge_lengths[edge_fibre_mask]
    )
    # Consistent mass is assembled on all mapped edges; generalized graph
    # eigenvalues solve K u = lambda M u.
    generalized = eigh(stiffness, mass, eigvals_only=True, check_finite=True)
    return FibredIntertwinerData(
        base_labels=base_labels,
        fibre_bits=fibre_bits,
        orbit_representatives=orbit_representatives,
        signs=signs,
        rgs_coordinates=mapped,
        q6_edges=edges,
        edge_lengths_rgs=edge_lengths,
        edge_lengths_physical_projection=physical_lengths,
        edge_base_mask=edge_base_mask,
        edge_fibre_mask=edge_fibre_mask,
        incidence=incidence,
        stiffness_base=stiffness_base,
        stiffness_fibre=stiffness_fibre,
        stiffness=stiffness,
        mass=mass,
        generalized_eigenvalues=generalized,
    )


def audit_fibred_intertwiner(
    data: FibredIntertwinerData | None = None,
    atol: float = 1e-10,
) -> dict:
    """Audit the finite base/fibre bijection and its induced 1D graph metric."""
    d = build_fibred_intertwiner_data() if data is None else data
    eig_k = np.linalg.eigvalsh(d.stiffness)
    eig_m = np.linalg.eigvalsh(d.mass)
    eig_base = np.linalg.eigvalsh(d.stiffness_base)
    eig_fibre = np.linalg.eigvalsh(d.stiffness_fibre)
    n_base = int(np.linalg.matrix_rank(d.stiffness_base, tol=atol))
    n_fibre = int(np.linalg.matrix_rank(d.stiffness_fibre, tol=atol))
    n_total = int(np.linalg.matrix_rank(d.stiffness, tol=atol))
    lengths = d.edge_lengths_rgs
    physical = d.edge_lengths_physical_projection
    base_lengths = lengths[d.edge_base_mask]
    fibre_lengths = lengths[d.edge_fibre_mask]
    positive_spectrum = d.generalized_eigenvalues[d.generalized_eigenvalues > atol]
    report = {
        "audit": "RGS_PSU_Q6_EXACT_FIBRED_INTERTWINER_V1",
        "factorization": {
            "base_group": "Z2^3",
            "sign_fibre_group": "Z2^3",
            "vertices": 64,
            "base_states": int(len(np.unique(d.base_labels))),
            "fibre_states_per_base": 8,
            "map": "Phi(b,s) = (2*s-1) componentwise-multiplied by p_b",
            "orbit_representatives": d.orbit_representatives.tolist(),
            "bijective_labelled_reconstruction": bool(
                np.array_equal(d.rgs_coordinates, load_rgs_seed_data().coordinates)
            ),
            "rgs_coordinate_reconstruction_max_abs": float(
                np.max(np.abs(d.rgs_coordinates - load_rgs_seed_data().coordinates))
            ),
        },
        "edge_geometry": {
            "total_edges": int(len(lengths)),
            "base_edges": int(np.count_nonzero(d.edge_base_mask)),
            "sign_fibre_edges": int(np.count_nonzero(d.edge_fibre_mask)),
            "rgs_edge_length_min": float(lengths.min()),
            "rgs_edge_length_max": float(lengths.max()),
            "rgs_edge_length_distinct_rounded_10": int(
                len(np.unique(np.round(lengths, 10)))
            ),
            "base_edge_length_min": float(base_lengths.min()),
            "base_edge_length_max": float(base_lengths.max()),
            "base_edge_length_mean": float(base_lengths.mean()),
            "sign_fibre_edge_length_min": float(fibre_lengths.min()),
            "sign_fibre_edge_length_max": float(fibre_lengths.max()),
            "sign_fibre_edge_length_mean": float(fibre_lengths.mean()),
            "golden_physical_projection_edge_length_min": float(physical.min()),
            "golden_physical_projection_edge_length_max": float(physical.max()),
            "golden_projection_expected_edge_length": float(1.0 / np.sqrt(2.0)),
        },
        "metric_graph_fem": {
            "action": "S(u)=1/2 sum_(edges) (u_v-u_w)^2 / ell_edge",
            "stiffness_base_rank": n_base,
            "stiffness_fibre_rank": n_fibre,
            "stiffness_total_rank": n_total,
            "base_component_kernel_dimension": int(np.count_nonzero(np.abs(eig_base) <= atol)),
            "fibre_component_kernel_dimension": int(np.count_nonzero(np.abs(eig_fibre) <= atol)),
            "stiffness_kernel_dimension": int(np.count_nonzero(np.abs(eig_k) <= atol)),
            "mass_min_eigenvalue": float(eig_m.min()),
            "mass_condition_number": float(np.linalg.cond(d.mass)),
            "stiffness_min_eigenvalue": float(eig_k.min()),
            "generalized_spectral_gap": float(positive_spectrum.min()),
            "generalized_max_eigenvalue": float(d.generalized_eigenvalues.max()),
            "generalized_zero_mode_count": int(
                np.count_nonzero(np.abs(d.generalized_eigenvalues) <= atol)
            ),
        },
        "checks": {
            "64_point_fibred_map_exact": bool(
                np.array_equal(d.rgs_coordinates, load_rgs_seed_data().coordinates)
            ),
            "96_base_plus_96_fibre_edges": bool(
                np.count_nonzero(d.edge_base_mask) == 96
                and np.count_nonzero(d.edge_fibre_mask) == 96
            ),
            "golden_projection_edges_uniform": bool(
                np.max(np.abs(physical - 1.0 / np.sqrt(2.0))) <= 1e-12
            ),
            "rgs_pullback_metric_nonuniform": bool(
                len(np.unique(np.round(lengths, 10))) > 1
            ),
            "base_and_fibre_subgraphs_each_rank_56": bool(
                n_base == 56 and n_fibre == 56
            ),
            "combined_graph_connected": bool(n_total == 63),
            "mass_positive_definite": bool(eig_m.min() > 0.0),
            "stiffness_positive_semidefinite": bool(eig_k.min() >= -atol),
            "single_constant_zero_mode": bool(
                np.count_nonzero(np.abs(eig_k) <= atol) == 1
            ),
        },
        "epistemic_status": {
            "finite_labelled_base_fibre_map": "EXACTLY_RECONSTRUCTED_FROM_ARCHIVED_SEED_ORDER",
            "induced_weighted_graph_fem": "NUMERICALLY_VERIFIED_FINITE_MODEL",
            "continuum_bundle_or_manifold": "OPEN",
            "curved_rgs_whitney_hodge_star": "OPEN",
            "TEGR_or_lorentzian_time_identification": "OPEN",
            "physical_interpretation_of_base_fibre_split": (
                "STRUCTURAL_CANDIDATE_NOT_YET_A_DYNAMICAL_DERIVATION"
            ),
        },
        "pass": bool(
            np.array_equal(d.rgs_coordinates, load_rgs_seed_data().coordinates)
            and np.count_nonzero(d.edge_base_mask) == 96
            and np.count_nonzero(d.edge_fibre_mask) == 96
            and np.max(np.abs(physical - 1.0 / np.sqrt(2.0))) <= 1e-12
            and n_base == 56 and n_fibre == 56 and n_total == 63
            and eig_m.min() > 0.0 and eig_k.min() >= -atol
            and np.count_nonzero(np.abs(eig_k) <= atol) == 1
        ),
    }
    return report


def export_fibred_intertwiner_audit(path: str | Path) -> dict:
    """Export deterministic JSON and NPZ audit artifacts."""
    d = build_fibred_intertwiner_data()
    report = audit_fibred_intertwiner(d)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    npz_path = target.with_suffix(".npz")
    json_path = target.with_suffix(".json")
    np.savez_compressed(
        npz_path,
        base_labels=d.base_labels,
        fibre_bits=d.fibre_bits,
        orbit_representatives=d.orbit_representatives,
        signs=d.signs,
        rgs_coordinates=d.rgs_coordinates,
        q6_edges=d.q6_edges,
        edge_lengths_rgs=d.edge_lengths_rgs,
        edge_lengths_physical_projection=d.edge_lengths_physical_projection,
        incidence=d.incidence,
        stiffness_base=d.stiffness_base,
        stiffness_fibre=d.stiffness_fibre,
        stiffness=d.stiffness,
        mass=d.mass,
        generalized_eigenvalues=d.generalized_eigenvalues,
    )
    import hashlib
    report["archive"] = {
        "npz": str(npz_path),
        "sha256": hashlib.sha256(npz_path.read_bytes()).hexdigest(),
    }
    json_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("results/rgs_psu_fibred_intertwiner"))
    args = parser.parse_args()
    print(json.dumps(export_fibred_intertwiner_audit(args.output), indent=2))

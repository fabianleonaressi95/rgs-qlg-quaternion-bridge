#!/usr/bin/env python3
"""Run the reproducible finite RGS–PSU / Möbius / Q6 geometry audit.

Writes JSON containing fixed seeds, tolerances, source/input SHA-256 hashes,
operator residuals, and explicit non-claims for continuum limits.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from src.binary_icosahedral import audit_binary_icosahedral, binary_icosahedral_elements
from src.clifford_qlg import (
    clifford_grade_action_audit, clifford_relation_audit,
    spinor_group_action_audit, spinor_steerable_kernel_audit,
)
from src.discrete_differential import audit_discrete_differential_operators
from src.discrete_cartan_geometry import (\n    audit_cartan_faces, audit_metric_compatibility,\n    build_metric_compatible_connection, edge_coframes_from_embedding,\n)
from src.mobius_complex import audit_mobius_complex, build_mobius_complex
from src.mobius_gauge_audit import audit_mobius_gauge_covariance, build_random_mobius_lattice
from src.mobius_representations import (
    audit_local_frame_representations, audit_reflection_transition_representation,
)
from src.q6_complex import audit_q6_complex, build_q6_complex, export_q6_incidence_bundle
from src.q6_prism_map import audit_q6_prism_map
from src.rgs_multiscale import equidistribution_moment_audit
from src.representations_2i import representation_audit, intertwiner_dimension_audit


SEEDS = (7, 17, 31, 20261103)
TOLERANCES = {
    "gauge_covariance_atol": 1e-10,
    "clifford_anticommutator": 2e-12,
    "binary_icosahedral": 2e-10,
    "intertwiner": 2e-9,
    "discrete_differential": 1e-12,
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_array(value: np.ndarray) -> str:
    array = np.ascontiguousarray(value)
    digest = hashlib.sha256()
    digest.update(str(array.dtype).encode("ascii"))
    digest.update(json.dumps(list(array.shape)).encode("ascii"))
    digest.update(array.tobytes())
    return digest.hexdigest()


def _matrix_inventory(root: Path) -> dict[str, Any]:
    folder = root / "data" / "rgs_188d"
    records: dict[str, Any] = {}
    matrices: dict[str, np.ndarray] = {}
    for path in sorted(folder.glob("*.npy")):
        try:
            array = np.load(path, allow_pickle=False)
            records[path.name] = {
                "shape": list(array.shape),
                "dtype": str(array.dtype),
                "sha256_file": sha256_file(path),
                "sha256_array_canonical": sha256_array(array),
            }
            if array.ndim == 2 and array.shape == (188, 188):
                matrices[path.name] = np.asarray(array, dtype=float)
        except Exception as exc:
            records[path.name] = {"load_error": f"{type(exc).__name__}: {exc}",
                                  "sha256_file": sha256_file(path)}

    hessian_candidates = [
        name for name, item in records.items()
        if item.get("shape") == [192, 192]
    ]
    operator_comparisons = {}
    names = sorted(matrices)
    for i, left in enumerate(names):
        a = matrices[left]
        for right in names[i + 1:]:
            b = matrices[right]
            denom = max(float(np.linalg.norm(a, ord="fro")),
                        float(np.linalg.norm(b, ord="fro")), 1e-300)
            operator_comparisons[f"{left} <> {right}"] = float(
                np.linalg.norm(a - b, ord="fro") / denom
            )

    inverse_checks = {}
    pairs = (
        ("K0.npy", "K0_inv.npy"),
        ("K0 (1).npy", "Kinv.npy"),
        ("K_tangent_rebuilt.npy", "Kinv_rebuilt.npy"),
    )
    for left, right in pairs:
        if left in records and right in records:
            if left in matrices and right in matrices:
                product = matrices[left] @ matrices[right]
                inverse_checks[f"{left} @ {right}"] = {
                    "identity_residual_fro": float(np.linalg.norm(
                        product - np.eye(product.shape[0]), ord="fro")),
                    "reverse_identity_residual_fro": float(np.linalg.norm(
                        matrices[right] @ matrices[left] - np.eye(product.shape[0]),
                        ord="fro")),
                }
            else:
                inverse_checks[f"{left} @ {right}"] = {
                    "status": "shape_not_188x188_or_load_failed"
                }

    t2u_files = (
        "all_T2u_arrays.npz", "K_tangent_rebuilt.npy",
        "K_T2u_rebuilt.npy", "B_T2U_TAN_rebuilt.npy", "J0_rebuilt.npy",
    )
    t2u_hashes = {
        name: records[name]["sha256_file"]
        for name in t2u_files if name in records and "sha256_file" in records[name]
    }
    t2u_reconstruction: dict[str, Any] = {"status": "not_recomputed"}
    archive = folder / "all_T2u_arrays.npz"
    if archive.exists():
        with np.load(archive, allow_pickle=False) as data:
            keys = list(data.files)
            t2u_reconstruction["archive_keys_and_shapes"] = {
                key: list(data[key].shape) for key in keys
            }
            eigen_key = next((k for k in keys if "eigenvectors_T2u_tan" in k), None)
            index_key = next((k for k in keys if "first_triplet_indices" in k), None)
            if eigen_key is not None and index_key is not None:
                rebuilt = np.asarray(data[eigen_key])[:, np.asarray(data[index_key], dtype=int)]
                target_path = folder / "B_T2U_TAN_rebuilt.npy"
                if target_path.exists():
                    target = np.load(target_path, allow_pickle=False)
                    t2u_reconstruction = {
                        "status": "recomputed_from_archive",
                        "formula": f"{eigen_key}[:, {index_key}]",
                        "shape": list(rebuilt.shape),
                        "relative_frobenius_difference": float(
                            np.linalg.norm(rebuilt-target) /
                            max(np.linalg.norm(target), 1e-300)
                        ),
                        "reconstruction_sha256": sha256_array(rebuilt),
                        "target_sha256": sha256_array(target),
                    }
            else:
                t2u_reconstruction["status"] = "required_formula_keys_not_found"
    return {
        "array_file_count": len(records),
        "arrays": records,
        "candidate_192x192_hessians": hessian_candidates,
        "original_192x192_hessian_recovered": bool(hessian_candidates),
        "same_shape_188x188_operator_relative_differences": operator_comparisons,
        "candidate_inverse_residuals": inverse_checks,
        "t2u_input_file_sha256": t2u_hashes,
        "t2u_formula_reconstruction": t2u_reconstruction,
        "interpretation": (
            "A file with shape 192x192 is only a candidate Hessian. Formula, "
            "constraint convention and provenance must still be established."
        ),
    }


def _random_gauge_input_hash(seed: int) -> str:
    lattice, gauge = build_random_mobius_lattice(seed)
    payload = []
    for source, target in sorted(lattice.edges()):
        payload.append(lattice.link(source, target).astype(np.complex128).tobytes())
    for vertex in sorted(gauge):
        payload.append(gauge[vertex].astype(np.complex128).tobytes())
    return hashlib.sha256(b"".join(payload)).hexdigest()


def build_integrated_audit(root: Path) -> dict[str, Any]:
    mobius = build_mobius_complex()
    q6 = build_q6_complex()
    mobius_face_sets = {frozenset(face) for face in mobius.faces}
    q6_face_lookup = {frozenset(face): i for i, face in enumerate(q6.faces)}
    q6_face_indices = [
        q6_face_lookup[face] for face in mobius_face_sets if face in q6_face_lookup
    ]
    mobius_q6 = {
        "mobius_vertices_edges_faces": [
            len(mobius.vertices), len(mobius.edges), len(mobius.faces)
        ],
        "q6_vertices_edges_faces": [
            len(q6.vertices), len(q6.edges), len(q6.faces)
        ],
        "mobius_faces_found_in_q6_240_faces": len(q6_face_indices),
        "mobius_q6_face_indices_sorted": sorted(q6_face_indices),
        "distinct_subcomplex_not_full_q6_face_set": (
            len(mobius.faces) == 6 and len(q6.faces) == 240
            and len(q6_face_indices) == 6
        ),
    }

    gauge = []
    for seed in SEEDS:
        result = audit_mobius_gauge_covariance(
            seed=seed, atol=TOLERANCES["gauge_covariance_atol"]
        )
        result["randomized_input_sha256"] = _random_gauge_input_hash(seed)
        gauge.append(result)

    group = audit_binary_icosahedral(tol=TOLERANCES["binary_icosahedral"])
    reps = representation_audit(tol=TOLERANCES["binary_icosahedral"])
    intertwiners = intertwiner_dimension_audit(rtol=1e-10)
    clifford = clifford_relation_audit(tol=TOLERANCES["clifford_anticommutator"])
    grade = clifford_grade_action_audit()
    spin = spinor_group_action_audit()
    kernels = spinor_steerable_kernel_audit()
    moment = equidistribution_moment_audit(
        binary_icosahedral_elements(),
        ((2, 0, 0, 0), (0, 2, 0, 0), (0, 0, 2, 0), (0, 0, 0, 2),
         (4, 0, 0, 0), (2, 2, 0, 0)),
    )

    source_paths = (
        "src/mobius_complex.py", "src/mobius_gauge_audit.py",
        "src/mobius_representations.py", "src/q6_complex.py",
        "src/q6_prism_map.py", "src/discrete_differential.py", "src/discrete_cartan_geometry.py",
        "src/binary_icosahedral.py", "src/representations_2i.py",
        "src/clifford_qlg.py", "src/rgs_multiscale.py",
        "data/q6/README.md",
        "data/rgs_188d/A14_jacobian_report.json",
        "data/rgs_188d/188D_T2U_RECONSTRUCTION_REPORT.json",
    )
    input_hashes = {
        name: sha256_file(root / name)
        for name in source_paths if (root / name).is_file()
    }
    q6_arrays = {}
    for name in ("Q6_edges.npy", "Q6_faces.npy"):
        path = root / "data" / "q6" / name
        if path.exists():
            arr = np.load(path, allow_pickle=False)
            q6_arrays[name] = {
                "shape": list(arr.shape),
                "sha256_file": sha256_file(path),
                "sha256_array_canonical": sha256_array(arr),
            }

    q6_bundle_path = root / "results" / "artifacts" / "Q6_canonical_incidence_v1.npz"
    q6_bundle_export = export_q6_incidence_bundle(str(q6_bundle_path))\n    q6_bundle_export["path"] = q6_bundle_path.relative_to(root).as_posix()\n    q6_bundle_manifest = q6_bundle_path.with_suffix(".manifest.json")\n    q6_bundle_manifest.parent.mkdir(parents=True, exist_ok=True)\n    q6_bundle_manifest.write_text(json.dumps(q6_bundle_export, indent=2, sort_keys=True) + "\\n", encoding="utf-8")\n\n    discrete_cartan = {"status": "64B_coordinate_seed_not_available"}\n    x_path = root / "data" / "rgs_188d" / "X_Oh_reconstructed.npy"\n    if x_path.exists():\n        x_oh = np.load(x_path, allow_pickle=False)\n        if x_oh.shape == (64, 3):\n            radial = np.linalg.norm(x_oh, axis=1)\n            local_coordinates = np.column_stack((x_oh, radial))\n            metrics = {\n                v: np.eye(4) + 0.01 * np.outer(local_coordinates[v], local_coordinates[v])\n                for v in q6.vertices\n            }\n            edge_theta = edge_coframes_from_embedding(\n                x_oh, q6.edges, append_radial_difference=True\n            )\n            edge_connection = build_metric_compatible_connection(metrics, q6.edges)\n            metric_report = audit_metric_compatibility(metrics, edge_connection, q6.edges)\n            cartan_report = audit_cartan_faces(q6.faces, edge_theta, edge_connection)\n            orthonormal_frame_residual = 0.0\n            frames = {v: np.linalg.cholesky(metrics[v]).T for v in q6.vertices}\n            for u, v in q6.edges:\n                omega = frames[v] @ edge_connection[(u, v)] @ np.linalg.inv(frames[u])\n                orthonormal_frame_residual = max(\n                    orthonormal_frame_residual,\n                    float(np.linalg.norm(omega.T @ omega - np.eye(4), ord="fro")),\n                )\n            discrete_cartan = {\n                "status": "coordinate_seed_induced_spd_metric_diagnostic",\n                "coordinate_seed_sha256": sha256_file(x_path),\n                "metric_compatibility": metric_report,\n                "cartan_faces": cartan_report,\n                "orthonormal_frame_transport_max_residual": orthonormal_frame_residual,\n                "metric_construction": "G_v = I + 0.01 y_v y_v^T, y_v=(X_v, ||X_v||); diagnostic metric, not an asserted recovered physical vierbein",\n            }\n\n    discrete_mobius = audit_discrete_differential_operators(
        mobius.boundary_1, mobius.boundary_2,
        atol=TOLERANCES["discrete_differential"],
    )
    discrete_q6 = audit_discrete_differential_operators(
        q6.boundary_1, q6.boundary_2,
        atol=TOLERANCES["discrete_differential"],
    )

    return {
        "schema_version": 1,
        "audit": "RGS_PSU_Q6_PRISM_MOBIUS_INTEGRATED_FINITE_AUDIT",
        "seeds": list(SEEDS),
        "tolerances": TOLERANCES,
        "input_sha256": input_hashes,
        "mobius": audit_mobius_complex(),
        "mobius_vs_q6": mobius_q6,
        "q6": audit_q6_complex(),
        "local_frame_representations": audit_local_frame_representations(),
        "reflection_frame_transitions": audit_reflection_transition_representation(),
        "randomized_gauge_covariance": gauge,
        "binary_icosahedral_2i": group,
        "2i_representations": reps,
        "2i_intertwiners": intertwiners,
        "clifford": clifford,
        "clifford_grade_action": grade,
        "spinor_group_action": spin,
        "steerable_kernels": kernels,
        "discrete_differential_mobius": discrete_mobius,
        "discrete_differential_q6": discrete_q6,
        "discrete_cartan_geometry": discrete_cartan,
        "q6_to_prism_map": audit_q6_prism_map(),
        "finite_2i_haar_moment_diagnostic": moment,
        "q6_source_arrays": q6_arrays,
        "q6_canonical_bundle_export": q6_bundle_export,
        "rgs_188d_operator_and_t2u_inventory": _matrix_inventory(root),
        "continuum_convergence": {
            "equidistribution": "finite 2I moments reported; a nested refinement sequence is required for weak convergence",
            "gromov_hausdorff": "not certified; no compatible sequence of metric correspondences/embeddings supplied",
            "mosco": "not certified; liminf and recovery-sequence proof obligations remain open",
            "multifractal": "partition-spectrum utility implemented; RGS nested partition masses not yet supplied",
            "navier_stokes": "no derivation from the present finite graph action to incompressible NS; conservation, scaling and hydrodynamic limit remain open",
            "yang_mills": "finite Wilson action exists; continuum YM requires gauge-covariant refinement and action/measure convergence",
            "teleparallel_tegr": "coframe-to-connection/torsion/curvature sequence remains a derivation gate; non-orientability alone gives neither torsion nor Pin structure",
        },
        "overall_status": "FINITE_AUDITS_ONLY_CONTINUUM_AND_PHYSICAL_LIMITS_OPEN",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path,
                        default=Path("results/json/rgs_psu_integrated_geometric_gauge_audit.json"))
    args = parser.parse_args()
    root = args.repo_root.resolve()
    output = args.output if args.output.is_absolute() else root / args.output
    report = build_integrated_audit(root)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    print(f"INTEGRATED_AUDIT_JSON={output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Run a finite Wilson spectral baseline on the recovered 64-vertex RGS/Q6 graph.

The local frame is intentionally set to identity: this isolates the actual
nonuniform embedded Q6 edge lengths without pretending that the archived
edge coframe is a vertexwise 3x3 tetrad or a continuum spin connection.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np

from src.curved_wilson_spectral_audit import curved_graph_wilson_spectral_audit
from src.rgs_ellis_refinement import load_rgs_seed_data


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--wilson-r", type=float, default=1.0)
    parser.add_argument("--cutoff", type=float, default=1e-8)
    parser.add_argument("--output", type=Path,
                        default=Path("results/json/rgs_q6_wilson_spectral_baseline.json"))
    args = parser.parse_args()

    seed = load_rgs_seed_data()
    coordinates = seed.coordinates
    edges = [(int(row[0]), int(row[1])) for row in seed.q6_edges]
    edge_lengths = {
        (i, j): float(np.linalg.norm(coordinates[j] - coordinates[i]))
        for i, j in edges
    }
    # Explicit reference choice: identity local coframe, actual RGS edge lengths.
    vierbein = np.broadcast_to(np.eye(3), (len(coordinates), 3, 3)).copy()
    links = {(i, j): np.eye(2, dtype=complex) for i, j in edges}
    report = curved_graph_wilson_spectral_audit(
        coordinates, vierbein, edges, spin_links=links,
        wilson_r=args.wilson_r, physical_cutoff=args.cutoff,
    )
    lengths = np.asarray(list(edge_lengths.values()), dtype=float)
    report["audit"] = "RECOVERED_RGS_Q6_EMBEDDED_GRAPH_WILSON_BASELINE"
    report["geometry"] = {
        "coordinate_source": "load_rgs_seed_data() with archived SHA-256 checks",
        "vertices": int(len(coordinates)),
        "edges": int(len(edges)),
        "edge_length_min": float(lengths.min()),
        "edge_length_max": float(lengths.max()),
        "edge_length_mean": float(lengths.mean()),
        "local_coframe": "identity at every vertex (flat-frame baseline)",
        "spin_transport": "identity SU(2) links",
        "edge_lengths_used_by_wilson_kernel": "Euclidean lengths of recovered RGS/Q6 embedded edges",
    }
    report["interpretation"] = (
        "Finite irregular-embedded-graph baseline only; not a curved tetrad, "
        "not a continuum spin connection, and not a refinement limit."
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    if not report["all_finite_checks_passed"]:
        raise SystemExit("RGS/Q6 finite Wilson baseline failed")


if __name__ == "__main__":
    main()

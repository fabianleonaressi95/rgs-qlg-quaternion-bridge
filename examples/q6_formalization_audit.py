"""Run and export the canonical Q6 formalization audit."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from src.q6_complex import build_q6_complex
from src.q6_spectral_audit import (
    antipodal_vertex_map, audit_q6_formalization, graph_laplacian,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=Path("results/q6"))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    q6 = build_q6_complex()
    report = audit_q6_formalization()
    np.save(args.out / "edges64_Q6_generated.npy", np.asarray(q6.edges, dtype=np.int64))
    np.save(args.out / "plaquettes_Q6_240.npy", np.asarray(q6.faces, dtype=np.int64))
    np.save(args.out / "B1_Q6.npy", q6.boundary_1)
    np.save(args.out / "B2_Q6.npy", q6.boundary_2)
    np.save(args.out / "L0_Q6.npy", graph_laplacian(q6))
    np.save(args.out / "antipodal_map_Q6.npy", antipodal_vertex_map())
    with (args.out / "q6_formalization_audit.json").open("w") as stream:
        json.dump(report, stream, indent=2)
    print(json.dumps({
        "full_pass": report["full_pass"],
        "vertices": report["vertex_count"],
        "edges": report["edge_count"],
        "faces": report["face_count"],
        "boundary_squared": report["boundary_squared_max_abs"],
        "graph_cycle_rank": report["homology"]["graph_cycle_rank_beta1"],
        "cellular_betti": [
            report["homology"]["cellular_betti_b0"],
            report["homology"]["cellular_betti_b1"],
            report["homology"]["cellular_betti_b2"],
        ],
        "output": str(args.out),
    }, indent=2))


if __name__ == "__main__":
    main()

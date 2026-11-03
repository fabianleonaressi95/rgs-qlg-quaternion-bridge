"""Run the A13 Q6-to-RGS tangent projection and Euclidean Killing audits.

Usage:
  PYTHONPATH=. python examples/rgs_projected_koopman_killing_audit.py \
      --matrices A13_all_matrices.npz --output results/json/rgs_projected_koopman_killing.json
"""
import argparse
import json
from pathlib import Path
import numpy as np

from src.koopman_rgs_psu import (
    discrete_euclidean_killing_audit,
    project_q6_laplacian_to_rgs_tangent,
)


def run(matrices_path, output_path):
    z = np.load(matrices_path)
    x_key = "X_Oh_reconstructed" if "X_Oh_reconstructed" in z else "X"
    q_key = "Q_tangent" if "Q_tangent" in z else "Q"
    X = np.asarray(z[x_key], dtype=float)
    Q = np.asarray(z[q_key], dtype=float)
    if "Q6_laplacian" in z:
        L6 = np.asarray(z["Q6_laplacian"], dtype=float)
    elif "Q6_adjacency" in z:
        A = np.asarray(z["Q6_adjacency"], dtype=float)
        L6 = np.diag(A.sum(axis=1)) - A
    else:
        raise KeyError("archive must contain Q6_laplacian or Q6_adjacency")
    if "Q6_adjacency" in z:
        A = np.asarray(z["Q6_adjacency"], dtype=float)
    else:
        A = np.zeros((64, 64), dtype=float)
        states = np.asarray(z["Q6_states"], dtype=int)
        for i in range(64):
            for j in range(i + 1, 64):
                if np.count_nonzero(states[i] != states[j]) == 1:
                    A[i, j] = A[j, i] = 1.0
    edges = np.argwhere(np.triu(A, 1))
    projected = project_q6_laplacian_to_rgs_tangent(Q, L6)
    killing = discrete_euclidean_killing_audit(X, edges)
    report = {
        "audit": "RGS_PROJECTED_KOOPMAN_AND_DISCRETE_KILLING",
        "source_archive": str(matrices_path),
        "source_keys": {"coordinates": x_key, "tangent_frame": q_key},
        "projection": {
            k: v for k, v in projected.items()
            if k not in {"generator", "stiffness", "ambient_laplacian"}
        },
        "killing": killing,
        "evidence_boundary": [
            "The projected operator is a finite Galerkin compression.",
            "A per-vertex curved vierbein and its connection are not inferred from X.",
            "The Euclidean Killing residual is not a curved-metric Killing test.",
        ],
    }
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(f"JSON report: {out}")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrices", required=True)
    parser.add_argument("--output", default="results/json/rgs_projected_koopman_killing.json")
    args = parser.parse_args()
    run(args.matrices, args.output)

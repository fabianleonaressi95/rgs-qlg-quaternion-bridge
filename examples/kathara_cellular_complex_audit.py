"""Export the exact Kathara cellular-complex audit as JSON."""
from __future__ import annotations
import argparse
import json
from pathlib import Path

from src.kathara_cellular_complex import audit_kathara_cellular_complex


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="results/kathara_cellular_complex/audit.json")
    args = parser.parse_args()
    a = audit_kathara_cellular_complex()
    payload = {
        "status": "exact finite combinatorial audit; 2-cell filling is a modeling convention",
        "dimensions_C0_C1_C2": [len(a.vertices), len(a.edges), len(a.faces)],
        "vertices": list(a.vertices),
        "oriented_edges": [list(e) for e in a.edges],
        "oriented_faces": [list(f) for f in a.faces],
        "B1_shape": [len(a.b1), len(a.b1[0])],
        "B2_shape": [len(a.b2), len(a.b2[0])],
        "B1": [list(row) for row in a.b1],
        "B2": [list(row) for row in a.b2],
        "B1B2_zero": a.b1b2_zero,
        "rank_B1": a.rank_b1,
        "rank_B2": a.rank_b2,
        "Betti_b0_b1_b2": list(a.betti),
        "Euler_characteristic": a.euler_characteristic,
        "Hodge_Laplacian_0": [list(row) for row in a.laplacian_0],
        "Hodge_Laplacian_1": [list(row) for row in a.laplacian_1],
        "Hodge_Laplacian_2": [list(row) for row in a.laplacian_2],
    }
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\\n", encoding="utf-8")
    print(json.dumps({k: payload[k] for k in (
        "status", "dimensions_C0_C1_C2", "B1B2_zero", "rank_B1",
        "rank_B2", "Betti_b0_b1_b2", "Euler_characteristic")}, indent=2))
    print(f"Wrote {path}")


if __name__ == "__main__":
    main()

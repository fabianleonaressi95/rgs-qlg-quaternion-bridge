from __future__ import annotations

import argparse
import json
from pathlib import Path

from src.kathara_square_similarity import audit_kathara_square_similarity


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--width", type=float, default=1.0)
    parser.add_argument("--output", default="results/kathara_square_similarity/audit.json")
    args = parser.parse_args()
    audit = audit_kathara_square_similarity(args.width)
    report = {
        "transform": audit.transform.tolist(),
        "transform_squared": audit.transform_squared.tolist(),
        "centers": {k: v.tolist() for k, v in audit.centers.items()},
        "residuals": audit.residuals,
        "lattice_image_index": audit.lattice_image_index,
        "lattice_parity_rule": audit.lattice_parity_rule,
        "scope": "Exact square similarity and center overlay constraints; not a derivation of the trigram tribonacci matrix.",
    }
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

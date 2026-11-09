from __future__ import annotations

import argparse
import json
from pathlib import Path

from src.q6_dyadic_refinement import audit_q6_dyadic_refinement


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="results/q6_dyadic_refinement/audit.json")
    args = parser.parse_args()
    report = audit_q6_dyadic_refinement()
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

"""Run the exact/numerical Q6 double-sheet entanglement audit."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from src.q6_bisheet_entanglement import audit_q6_bisheet_entanglement


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path,
                        default=Path("results/json/q6_bisheet_entanglement.json"))
    args = parser.parse_args()
    report = audit_q6_bisheet_entanglement()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if not report["pass"]:
        raise SystemExit("Q6 double-sheet entanglement audit failed")


if __name__ == "__main__":
    main()

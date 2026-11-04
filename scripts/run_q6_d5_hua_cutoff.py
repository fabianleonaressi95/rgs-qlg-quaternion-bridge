"""Run and persist the exact D5 Hua cutoff-normalization audit."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from src.rgs_q6_d5_hua_cutoff import d5_hua_cutoff_report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("reports/q6_d5_hua_cutoff_normalization.json"),
    )
    args = parser.parse_args()
    result = d5_hua_cutoff_report()
    result["execution"] = {
        "command": "python scripts/run_q6_d5_hua_cutoff.py",
        "normalization": "relative to normalized Euclidean volume on D5",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    args.output.write_text(payload, encoding="utf-8")
    print(payload, end="")
    if not result["all_checks_passed"]:
        raise SystemExit("D5 Hua cutoff audit checks failed")


if __name__ == "__main__":
    main()

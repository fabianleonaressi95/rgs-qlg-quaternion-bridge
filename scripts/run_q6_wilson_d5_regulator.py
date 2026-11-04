"""Run and persist the Q6 Wilson-to-D5 regulator identifiability audit."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.rgs_q6_d5_wilson_regulator import audit_q6_wilson_d5_regulator_obstruction


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path,
                        default=Path("reports/q6_wilson_d5_regulator.json"))
    parser.add_argument("--beta-time", type=float, default=80.0)
    parser.add_argument("--temporal-spacing", type=float, default=0.05)
    parser.add_argument("--minimum-transfer-weight", type=float, default=0.25)
    args = parser.parse_args()
    report = audit_q6_wilson_d5_regulator_obstruction(
        beta_time=args.beta_time,
        temporal_spacing=args.temporal_spacing,
        minimum_transfer_weight=args.minimum_transfer_weight,
    )
    report["execution"] = {
        "command": "python scripts/run_q6_wilson_d5_regulator.py",
        "parameters": {
            "beta_time": args.beta_time,
            "temporal_spacing": args.temporal_spacing,
            "minimum_transfer_weight": args.minimum_transfer_weight,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    args.output.write_text(payload, encoding="utf-8")
    print(payload, end="")
    if not report["all_checks_passed"]:
        raise SystemExit("Q6 Wilson-D5 regulator audit failed")


if __name__ == "__main__":
    main()

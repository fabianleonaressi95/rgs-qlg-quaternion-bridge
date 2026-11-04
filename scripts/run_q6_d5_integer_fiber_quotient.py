"""Run the exact Q6 integer-fiber quotient audit."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.rgs_q6_d5_integer_fiber_quotient import audit_q6_d5_integer_fiber_quotient


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path,
                        default=Path("reports/q6_d5_integer_fiber_quotient.json"))
    parser.add_argument("--beta-time", type=float, default=80.0)
    args = parser.parse_args()
    report = audit_q6_d5_integer_fiber_quotient(beta_time=args.beta_time)
    report["execution"] = {
        "command": "python scripts/run_q6_d5_integer_fiber_quotient.py",
        "beta_time": args.beta_time,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    args.output.write_text(payload, encoding="utf-8")
    print(payload, end="")
    if not report["all_checks_passed"]:
        raise SystemExit("Q6-D5 integer-fiber quotient audit failed")


if __name__ == "__main__":
    main()

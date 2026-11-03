"""Reproduce the RGS -> Q6 -> SU(2)^192 -> Wilson audit."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from src.rgs_q6_su2_wilson import audit_rgs_q6_su2_wilson


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=20261103)
    parser.add_argument("--beta", type=float, default=1.0)
    parser.add_argument("--output", type=Path,
                        default=Path("results/json/rgs_q6_su2_wilson_bridge.json"))
    args = parser.parse_args()
    report = audit_rgs_q6_su2_wilson(args.seed, args.beta)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if not report["pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

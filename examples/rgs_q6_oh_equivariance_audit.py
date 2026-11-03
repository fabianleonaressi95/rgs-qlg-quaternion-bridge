"""Export the recovered RGS/O_h/Q6 equivariance audit."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from src.rgs_q6_oh_action import audit_rgs_q6_oh


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path,
                        default=Path("results/json/rgs_q6_oh_equivariance.json"))
    args = parser.parse_args()
    report = audit_rgs_q6_oh()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if not report["pass"]:
        raise SystemExit("Recovered coordinate O_h action audit failed")


if __name__ == "__main__":
    main()

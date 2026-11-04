"""Run the finite RGS-Q6 induced metric audit and export NPZ + JSON."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from src.rgs_q6_induced_metric import export_rgs_q6_induced_metric


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("results/rgs_q6_induced_metric"))
    args = parser.parse_args()
    archive = args.output.with_suffix(".npz")
    report = args.output.with_suffix(".json")
    result = export_rgs_q6_induced_metric(str(archive))
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(json.dumps(result["report"], indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

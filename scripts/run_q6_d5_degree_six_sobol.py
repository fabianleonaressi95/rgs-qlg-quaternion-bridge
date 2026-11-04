"""Run and persist the reproducible Q6-D5 degree-six Sobol audit."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from src.rgs_q6_d5_degree_six_sobol import (
    build_q6_d5_degree_six_sobol_crosscheck,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--power", type=int, default=12)
    parser.add_argument("--scrambles", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20261105)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("reports/q6_d5_degree_six_sobol.json"),
    )
    args = parser.parse_args()

    result = build_q6_d5_degree_six_sobol_crosscheck(
        power=args.power, scrambles=args.scrambles, seed=args.seed
    )
    result["execution"] = {
        "command": (
            "python scripts/run_q6_d5_degree_six_sobol.py "
            f"--power {args.power} --scrambles {args.scrambles} --seed {args.seed}"
        ),
        "python_seed": args.seed,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\\n"
    args.output.write_text(payload, encoding="utf-8")
    print(payload, end="")
    if not result["all_checks_passed"]:
        raise SystemExit("Q6-D5 Sobol audit checks failed")


if __name__ == "__main__":
    main()

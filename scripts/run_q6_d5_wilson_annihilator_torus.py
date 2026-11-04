"""Run the exact Q6-D5 annihilator-torus Wilson marginal audit."""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.rgs_q6_d5_wilson_fiber_torus import q6_d5_wilson_annihilator_torus_integral


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--beta-time", type=float, default=2.0)
    parser.add_argument("--power", type=int, default=12)
    parser.add_argument("--scrambles", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20261104)
    parser.add_argument("--output", type=Path, default=Path("reports/q6_d5_wilson_annihilator_torus.json"))
    args = parser.parse_args()
    report = q6_d5_wilson_annihilator_torus_integral(
        beta_time=args.beta_time, power=args.power, scrambles=args.scrambles, seed=args.seed
    )
    report["execution"] = {
        "beta_time": args.beta_time, "power": args.power,
        "scrambles": args.scrambles, "seed": args.seed,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    args.output.write_text(payload, encoding="utf-8")
    print(payload, end="")


if __name__ == "__main__":
    main()

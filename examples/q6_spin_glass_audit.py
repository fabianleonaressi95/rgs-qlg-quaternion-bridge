"""Run a reproducible Edwards–Anderson audit on Q6."""
from __future__ import annotations

import argparse
import json

from src.rgs_q6_spin_glass import Q6SpinGlass, parallel_tempering


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=20261104)
    parser.add_argument("--distribution", choices=("gaussian", "bimodal"), default="gaussian")
    parser.add_argument("--sweeps", type=int, default=2000)
    parser.add_argument("--burn-in", type=int, default=1000)
    parser.add_argument("--measure-every", type=int, default=10)
    args = parser.parse_args()

    model = Q6SpinGlass.random(seed=args.seed, distribution=args.distribution)
    report = parallel_tempering(
        model, sweeps=args.sweeps, burn_in=args.burn_in,
        measure_every=args.measure_every, seed=args.seed + 1,
    )
    report["disorder_seed"] = args.seed
    report["bond_distribution"] = args.distribution
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

"""Generate the canonical golden cut-and-project RGS–PSU–Q6 audit."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from src.rgs_ellis_refinement import load_rgs_seed_data
from src.rgs_psu_cut_project import export_cut_project_audit


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("results/rgs_psu_cut_project"))
    parser.add_argument("--physical-radius", type=float, default=2.5)
    parser.add_argument("--torsion-trace", type=Path, default=None)
    parser.add_argument("--permutations", type=int, default=5000)
    parser.add_argument("--bootstrap", type=int, default=2000)
    args = parser.parse_args()
    seed = load_rgs_seed_data()
    torsion = np.load(args.torsion_trace) if args.torsion_trace else None
    report = export_cut_project_audit(
        args.output,
        args.physical_radius,
        rgs_coordinates=seed.coordinates,
        torsion_trace=torsion,
        n_permutations=args.permutations,
        n_bootstrap=args.bootstrap,
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

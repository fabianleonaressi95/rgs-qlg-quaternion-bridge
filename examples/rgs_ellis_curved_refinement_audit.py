#!/usr/bin/env python3
"""Run the curved RGS Ellis refinement G1 audit."""
from __future__ import annotations

import json
from pathlib import Path

from src.rgs_ellis_refinement import audit_ellis_rgs_refinement


def main() -> None:
    report = audit_ellis_rgs_refinement(
        levels=(0, 1, 2),
        r0=0.8292,
        half_width=3.0 * 0.8292,
        radial_base=8,
    )
    output = Path("results/json/rgs_ellis_curved_refinement_g1.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(f"Saved {output}")


if __name__ == "__main__":
    main()

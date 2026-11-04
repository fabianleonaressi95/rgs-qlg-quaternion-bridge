#!/usr/bin/env python3
"""Export the explicit source-defined Q6 incidence bundle and manifest."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from src.q6_complex import export_q6_incidence_bundle


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("results/artifacts/Q6_canonical_incidence_v1.npz"))
    args = parser.parse_args()
    report = export_q6_incidence_bundle(str(args.output))
    manifest = args.output.with_suffix(".manifest.json")
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text(json.dumps(report, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

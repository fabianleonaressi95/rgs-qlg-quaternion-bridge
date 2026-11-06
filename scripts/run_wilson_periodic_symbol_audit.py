"""Run the exact Brillouin-corner Wilson-symbol audit."""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.wilson_periodic_symbol import wilson_brillouin_corner_audit, wilson_scaling_window

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dimension", type=int, default=3)
    parser.add_argument("--spacing", type=float, default=1.0)
    parser.add_argument("--wilson-r", type=float, default=1.0)
    parser.add_argument("--alpha", type=float, default=0.0)
    parser.add_argument("--output", type=Path, default=Path("reports/wilson_periodic_symbol_audit.json"))
    args = parser.parse_args()
    report = {
        "audit": "PERIODIC_EUCLIDEAN_WILSON_SYMBOL",
        "scope": "exact Fourier-symbol corner audit; not a curved-Q6 continuum theorem",
        "corner_audit": wilson_brillouin_corner_audit(args.dimension, args.spacing, args.wilson_r),
        "scaling": wilson_scaling_window(args.alpha),
    }
    report["all_checks_passed"] = all(report["corner_audit"]["checks"].values())
    report["execution"] = {"command": "python scripts/run_wilson_periodic_symbol_audit.py",
                           "parameters": {"dimension": args.dimension, "spacing": args.spacing,
                                          "wilson_r": args.wilson_r, "alpha": args.alpha}}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    args.output.write_text(payload, encoding="utf-8")
    print(payload, end="")
    if not report["all_checks_passed"]:
        raise SystemExit("periodic Wilson-symbol audit failed")

if __name__ == "__main__":
    main()

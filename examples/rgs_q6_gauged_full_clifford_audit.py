"""Reproducible audit for the gauge-covariant full cubical Q6 operator."""
import json
from pathlib import Path

from src.rgs_q6_gauged_full_clifford import audit_q6_gauged_full_clifford


if __name__ == "__main__":
    report = audit_q6_gauged_full_clifford()
    out = Path("results/json/rgs_q6_gauged_full_clifford.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if not report["all_checks_passed"]:
        raise SystemExit("gauged full cubical Q6 audit failed")

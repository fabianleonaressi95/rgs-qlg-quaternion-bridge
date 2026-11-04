#!/usr/bin/env python3
"""Generate a seeded JSON audit for spinor Kirchhoff currents on Q6."""
import argparse
import json
from pathlib import Path
from src.spinor_kirchhoff import write_audit_json

p = argparse.ArgumentParser()
p.add_argument("--seed", type=int, default=20261104)
p.add_argument("--gauge-seed", type=int, default=20261105)
p.add_argument("--atol", type=float, default=1e-10)
p.add_argument("--output", type=Path, default=Path("results/json/spinor_kirchhoff_q6_audit.json"))
a = p.parse_args()
report = write_audit_json(a.output, a.seed, a.gauge_seed, a.atol)
print(json.dumps(report, indent=2, sort_keys=True))
if not report["passed"]:
    raise SystemExit(1)

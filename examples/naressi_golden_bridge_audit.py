"""Reproduce the Phase 20 Naressi golden-coupling bridge audit."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from src.naressi_golden_bridge import (
    audit_golden_naressi_bridge, audit_qat_unitarity, tfim_ground_state_audit,
)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--sizes",nargs="+",type=int,default=[4,8,12,16])
    parser.add_argument("--g",type=float,default=1.0)
    parser.add_argument("--out",type=Path,default=Path("results/json/naressi_golden_bridge.json"))
    args=parser.parse_args()
    report={
        "bridge":audit_golden_naressi_bridge(),
        "qat_unitarity":audit_qat_unitarity(),
        "tfim_finite_size":[tfim_ground_state_audit(n,g=args.g) for n in args.sizes],
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps(report,indent=2))
    print(f"JSON report: {args.out}")

if __name__=="__main__":
    main()

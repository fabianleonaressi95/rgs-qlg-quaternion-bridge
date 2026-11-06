"""Audit curved common-L2 mass against the lumped P1 vertex measure."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from src.rgs_ellis_refinement import (
    audit_ellis_lumped_mass_equivalence,
    build_ellis_rgs_refinement,
)


def run(levels=(0,1),radial_base=8,quadrature_order=4):
    rows=[]
    for level in levels:
        mesh=build_ellis_rgs_refinement(
            level,r0=.8292,radial_base=radial_base,include_chain=False
        )
        rows.append({
            "level":int(level),
            **audit_ellis_lumped_mass_equivalence(mesh,quadrature_order),
        })
    checks={
        "all_generalized_masses_positive":all(x["positive_mass"] for x in rows),
        "finite_equivalence_constants":all(
            0.0<x["generalized_mass_eigenvalue_min"]<
            x["generalized_mass_eigenvalue_max"]<2.0 for x in rows
        ),
        "curved_mass_volume_accurate":all(
            x["curved_mass_volume_relative_error"]<0.03 for x in rows
        ),
    }
    return {
        "audit":"ELLIS_RGS_CURVED_L2_VS_LUMPED_MASS",
        "scope":"FINITE_LEVEL_NUMERICAL_EQUIVALENCE_ONLY",
        "parameters":{"levels":list(levels),"radial_base":radial_base,
                      "quadrature_order":quadrature_order},
        "levels":rows,
        "checks":checks,
        "all_checks_passed":all(checks.values()),
        "open":"uniform-in-level spectral equivalence and a quadrature-stable analytic bound",
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--levels",nargs="+",type=int,default=[0,1])
    p.add_argument("--radial-base",type=int,default=8)
    p.add_argument("--quadrature-order",type=int,default=4)
    p.add_argument("--output",type=Path)
    a=p.parse_args()
    report=run(a.levels,a.radial_base,a.quadrature_order)
    payload=json.dumps(report,indent=2,sort_keys=True)
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(payload+"\n",encoding="utf-8")
    print(payload)
    if not report["all_checks_passed"]:
        raise SystemExit(1)

if __name__=="__main__":
    main()

"""Multi-level shape-regularity and geometric convergence audit for Ellis-RGS meshes."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from src.rgs_ellis_refinement import audit_ellis_rgs_refinement


def run(levels=(0,1,2), radial_base=8, r0=0.8292):
    report=audit_ellis_rgs_refinement(
        levels=tuple(levels),r0=r0,radial_base=radial_base
    )
    rows=report["refinements"]
    for row in rows:
        row["shape_regular_proxy"] = float(row["element_condition_max"])
        row["shape_regular_inverse_proxy"] = float(1.0/row["element_condition_max"])
    h_rates=[]
    volume_rates=[]
    for left,right in zip(rows,rows[1:]):
        h_rates.append(float(np.log(left["h_max_chord"]/right["h_max_chord"])/np.log(2.0)))
        e0,e1=left["volume_relative_error"],right["volume_relative_error"]
        volume_rates.append(float(np.log(e0/e1)/np.log(2.0)) if e0>0 and e1>0 else None)
    checks={
        "diameter_strictly_decreases":all(a["h_max_chord"]>b["h_max_chord"] for a,b in zip(rows,rows[1:])),
        "volume_error_strictly_decreases":all(a["volume_relative_error"]>b["volume_relative_error"] for a,b in zip(rows,rows[1:])),
        "shape_condition_bounded_on_tested_levels":max(r["element_condition_max"] for r in rows)<100.0,
        "angular_triangle_quality_positive":min(r["angular_triangle_quality_min"] for r in rows)>0.0,
        "angular_triangle_edge_ratio_positive":min(r["angular_triangle_edge_ratio_min"] for r in rows)>0.0,
        "radial_angular_scale_ratio_positive":min(r["radial_step_over_angular_h_max"] for r in rows)>0.0,
        "tetra_normalized_volume_positive":min(r["tetra_normalized_volume_min"] for r in rows)>0.0,
        "gram_minimum_positive":min(r["element_gram_min_eigenvalue"] for r in rows)>0.0,
        "embedding_identity_small":max(r["embedding_identity_max_residual"] for r in rows)<1e-13,
    }
    return {
        "audit":"ELLIS_RGS_MULTI_LEVEL_SHAPE_REGULARITY_AND_GEOMETRIC_CONVERGENCE",
        "scope":"FINITE_MESH_NUMERICAL_EVIDENCE_NOT_A_UNIFORM_PROOF",
        "parameters":{"levels":list(levels),"radial_base":int(radial_base),"r0":float(r0)},
        "refinements":rows,
        "observed_h_halving_orders":h_rates,
        "observed_volume_error_orders":volume_rates,
        "checks":checks,
        "all_checks_passed":bool(all(checks.values())),
        "proof_obligations":[
            "derive a level-independent analytic bound for the tetrahedron Gram condition numbers",
            "prove uniform equivalence of discrete and continuum L2 measures",
            "prove discrete Poincare and Rellich compactness",
            "construct recovery sequences for the curved Ellis energy",
            "establish Mosco and strong-resolvent convergence",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--levels",nargs="+",type=int,default=[0,1,2])
    p.add_argument("--radial-base",type=int,default=8)
    p.add_argument("--r0",type=float,default=.8292)
    p.add_argument("--output",type=Path)
    a=p.parse_args()
    result=run(a.levels,a.radial_base,a.r0)
    payload=json.dumps(result,indent=2,sort_keys=True)
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(payload+"\n",encoding="utf-8")
    print(payload)
    if not result["all_checks_passed"]:
        raise SystemExit(1)

if __name__=="__main__":
    main()

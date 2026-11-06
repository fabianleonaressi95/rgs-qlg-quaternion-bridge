"""Run the full level-0 Ellis-RGS intrinsic spin-Wilson audit."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np

from src.ellis_spin_connection import ellis_edge_spin_links
from src.ellis_rgs_spin_dirac import assemble_ellis_spin_wilson_dirac
from src.rgs_ellis_refinement import build_ellis_rgs_refinement


def run(level=0, radial_base=2, r0=0.8292, steps=12, wilson_r=0.7):
    mesh=build_ellis_rgs_refinement(
        level=level,r0=r0,radial_base=radial_base,include_chain=True
    )
    edges=[tuple(map(int,e)) for e in mesh.chain["edges"]]
    links,link_audit=ellis_edge_spin_links(
        mesh.vertices_r4,edges,r0,steps=steps
    )
    operator,audit=assemble_ellis_spin_wilson_dirac(
        mesh,edges=edges,spin_links=links,wilson_r=wilson_r
    )
    checks={
        "chain_complex_exact":bool(mesh.chain["chain_complex_exact"]),
        "all_operator_entries_finite":bool(np.all(np.isfinite(operator))),
        "dirac_hermitian":audit["dirac_hermiticity_residual"]<1e-8,
        "wilson_hermitian":audit["wilson_hermiticity_residual"]<1e-8,
        "gamma5_hermitian":audit["gamma5_hermiticity_residual"]<1e-7,
        "wilson_positive_semidefinite":audit["wilson_minimum_eigenvalue"]>-1e-7,
        "link_transport_residual_small":link_audit["max_transport_residual"]<1e-7,
    }
    return {
        "audit":"FULL_LEVEL_ELLIS_RGS_INTRINSIC_SPIN_WILSON_DIRAC",
        "scope":"FINITE_MESH_ONLY",
        "mesh":{
            "level":int(level),"radial_base":int(radial_base),
            "vertices":int(len(mesh.vertices_r4)),
            "tetrahedra":int(len(mesh.tetrahedra)),
            "edges":int(len(edges)),
            "max_element_diameter":float(mesh.max_element_diameter),
            "max_element_condition":float(np.max(mesh.element_condition)),
            "lumped_volume_total":float(np.sum(audit["lumped_volume_min"]*0+0)
                + np.sum(mesh.element_volumes)),
        },
        "connection":link_audit,
        "operator":{k:v for k,v in audit.items()
                    if k!="gauge_covariance_transform"},
        "matrix_dimension":int(operator.shape[0]),
        "checks":checks,
        "all_checks_passed":bool(all(checks.values())),
        "limitations":[
            "one finite mesh does not establish a refinement-uniform estimate",
            "the graph Wilson term is not proven to remove all curved-lattice doublers",
            "no Mosco or strong-resolvent convergence is claimed",
        ],
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--level",type=int,default=0)
    parser.add_argument("--radial-base",type=int,default=2)
    parser.add_argument("--steps",type=int,default=12)
    parser.add_argument("--wilson-r",type=float,default=.7)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    report=run(args.level,args.radial_base,steps=args.steps,wilson_r=args.wilson_r)
    text=json.dumps(report,indent=2,sort_keys=True)
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(text+"\n",encoding="utf-8")
    print(text)
    if not report["all_checks_passed"]:
        raise SystemExit(1)

if __name__=="__main__":
    main()

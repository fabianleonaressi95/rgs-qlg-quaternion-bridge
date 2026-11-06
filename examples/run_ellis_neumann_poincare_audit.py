"""Finite-level Neumann spectral-gap audit for the curved Ellis P1 energy."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import eigsh
from src.rgs_ellis_refinement import (
    assemble_ellis_common_l2_operators,
    build_ellis_rgs_refinement,
)


def run(levels=(0,1),radial_base=8,quadrature_order=4):
    rows=[]
    for level in levels:
        mesh=build_ellis_rgs_refinement(
            level,r0=.8292,radial_base=radial_base,include_chain=False
        )
        mass,stiffness=assemble_ellis_common_l2_operators(
            mesh,quadrature_order=quadrature_order
        )
        # Natural Neumann form: constants are the one-dimensional kernel.
        values,vectors=eigsh(
            stiffness,k=4,M=mass,sigma=1e-8,which="LM",tol=1e-9
        )
        order=np.argsort(values)
        values=np.asarray(values[order],dtype=float)
        vectors=np.asarray(vectors[:,order],dtype=float)
        constant=np.ones(mesh.vertices_r4.shape[0])
        zero_residual=float(np.linalg.norm(stiffness@constant))
        gram=np.asarray(vectors.T@(mass@vectors))
        rows.append({
            "level":int(level),
            "vertices":int(mesh.vertices_r4.shape[0]),
            "tetrahedra":int(mesh.tetrahedra.shape[0]),
            "neumann_eigenvalues_low":values.tolist(),
            "constant_mode_stiffness_residual":zero_residual,
            "first_positive_eigenvalue":float(values[1]),
            "discrete_poincare_constant":float(1.0/values[1]),
            "mass_orthonormality_residual":float(np.linalg.norm(gram-np.eye(4))),
            "positive_gap":bool(values[1]>0.0),
            "single_constant_kernel_numerically":bool(abs(values[0])<1e-7 and values[1]>1e-6),
        })
    checks={
        "all_levels_have_positive_neumann_gap":all(x["positive_gap"] for x in rows),
        "constant_is_kernel":all(x["constant_mode_stiffness_residual"]<1e-7 for x in rows),
        "eigenvectors_mass_orthonormal":all(x["mass_orthonormality_residual"]<1e-6 for x in rows),
        "one_constant_kernel_detected":all(x["single_constant_kernel_numerically"] for x in rows),
    }
    return {
        "audit":"ELLIS_RGS_FINITE_LEVEL_NEUMANN_POINCARE_GAP",
        "scope":"FINITE_LEVEL_SPECTRAL_EVIDENCE_NOT_A_UNIFORM_RELlich_PROOF",
        "parameters":{"levels":list(levels),"radial_base":radial_base,
                      "quadrature_order":quadrature_order},
        "levels":rows,
        "checks":checks,
        "all_checks_passed":all(checks.values()),
        "open":"uniform discrete Poincare inequality, compactness/Rellich, and spectral convergence",
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

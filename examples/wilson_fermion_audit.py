"""Reproducible finite-graph Wilson-term audit.

Run from the repository root:
    python examples/wilson_fermion_audit.py
"""
"""Reproducible audit distinguishing Euclidean Wilson from a scalar Hermitian sum."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from src.discrete_dirac import (
    assemble_spatial_dirac, assemble_wilson_laplacian, chirality_matrix,
)

    assemble_euclidean_wilson_dirac, assemble_spatial_dirac,
    assemble_wilson_laplacian, chirality_matrix,
)
from src.wilson_refinement import (
    euclidean_wilson_singular_value, hermitian_scalar_sum_branches,
    scalar_sum_accidental_zero_angle,
)


def main():
    n=16
    x=np.column_stack((np.arange(n,dtype=float),np.zeros(n),np.zeros(n)))
    e=np.broadcast_to(np.eye(3),(n,3,3)).copy()
    edges=[(i,i+1) for i in range(n-1)]+[(0,n-1)]
    periodic_lengths={(i,j):1.0 for i,j in edges}
    wilson,audit=assemble_wilson_laplacian(
        x,e,edges,wilson_r=1.0,edge_lengths=periodic_lengths
    )
    def mode_value(k):
        phase=np.exp(2j*np.pi*k*np.arange(n)/n)
        v=np.zeros(4*n,dtype=complex)
        v[0::4]=phase/np.sqrt(n)
        return float(np.real(np.vdot(v,wilson@v)))
    p0=mode_value(0)
    ppi=mode_value(n//2)

    # Independent two-site test: the naive operator anticommutes with gamma_5,
    # while the added scalar Wilson Laplacian breaks that exact anticommutation.
    x2=np.array([[0.,0.,0.],[1.,0.,0.]])
    e2=np.broadcast_to(np.eye(3),(2,3,3)).copy()
    d0,_=assemble_spatial_dirac(x2,e2,[(0,1)])
    dw,dw_audit=assemble_spatial_dirac(x2,e2,[(0,1)],wilson_r=1.0)
    g5=chirality_matrix()
    gamma5=np.kron(np.eye(2),g5)
    naive_chiral=float(np.linalg.norm(d0@gamma5+gamma5@d0,ord="fro"))
    wilson_chiral=float(np.linalg.norm(dw@gamma5+gamma5@dw,ord="fro"))
    report={
        "audit":"RGS_QLG_WILSON_FERMION_AUDIT",
        "periodic_ring":{"sites":n,"spacing":1.0,"wilson_r":1.0,
                         "p_zero_eigenvalue":p0,"p_pi_eigenvalue":ppi,
                         "expected_p_pi":2.0,
                         "minimum_operator_eigenvalue":audit["minimum_eigenvalue"],
                         "hermiticity_residual":audit["hermiticity_residual"]},
        "chiral_test":{"naive_anticommutator_norm":naive_chiral,
                       "wilson_anticommutator_norm":wilson_chiral,
                       "corrected_dirac_hermiticity_residual":dw_audit["hermiticity_residual"]},
        "passed":bool(abs(p0)<1e-11 and abs(ppi-2.0)<1e-11
                      and audit["minimum_eigenvalue"]>-1e-11
                      and naive_chiral<1e-11 and wilson_chiral>1e-8
                      and dw_audit["hermiticity_residual"]<1e-11),
        "scope":"finite graph algebra and controlled 1D periodic Fourier check; not a Q6 continuum no-doubling theorem"
    w,audit=assemble_wilson_laplacian(
        x,e,edges,wilson_r=1.0,edge_lengths=periodic_lengths
    )
    def wilson_mode(k):
        phase=np.exp(2j*np.pi*k*np.arange(n)/n)
        v=np.zeros(4*n,dtype=complex); v[0::4]=phase/np.sqrt(n)
        return float(np.real(np.vdot(v,w@v)))
    p0=wilson_mode(0); ppi=wilson_mode(n//2)

    # At this angle the Hermitian scalar sum D0+W has an accidental zero.
    a=1.0; r=1.0
    theta=scalar_sum_accidental_zero_angle(r)
    scalar_low,scalar_high=hermitian_scalar_sum_branches([theta/a],a,r)
    euclidean_at_accidental=euclidean_wilson_singular_value([theta/a],a,r)
    euclidean_at_pi=euclidean_wilson_singular_value([np.pi/a],a,r)

    # The Euclidean Wilson matrix is gamma5-Hermitian, not ordinarily Hermitian.
    x2=np.array([[0.,0.,0.],[1.,0.,0.]])
    e2=np.broadcast_to(np.eye(3),(2,3,3)).copy()
    d0,_=assemble_spatial_dirac(x2,e2,[(0,1)])
    de,de_audit=assemble_euclidean_wilson_dirac(x2,e2,[(0,1)],wilson_r=1.0)
    gamma5=np.kron(np.eye(2),chirality_matrix())
    naive_chiral=float(np.linalg.norm(d0@gamma5+gamma5@d0,ord="fro"))
    euclidean_gamma5_herm=float(np.linalg.norm(
        de.conj().T-gamma5@de@gamma5,ord="fro"
    ))
    ordinary_herm=float(np.linalg.norm(de-de.conj().T,ord="fro"))

    report={
        "audit":"RGS_QLG_EUCLIDEAN_WILSON_CONVENTION_AUDIT",
        "periodic_ring":{"sites":n,"spacing":1.0,"wilson_r":1.0,
            "W_p0":p0,"W_ppi":ppi,"expected_W_ppi":2.0,
            "wilson_laplacian_minimum_eigenvalue":audit["minimum_eigenvalue"],
            "wilson_laplacian_hermiticity_residual":audit["hermiticity_residual"]},
        "hermitian_scalar_sum_warning":{
            "accidental_zero_angle":theta,
            "expected_angle_for_r1":float(2*np.arctan(2.0)),
            "lower_branch_at_accidental_zero":scalar_low,
            "upper_branch_at_accidental_zero":scalar_high},
        "euclidean_wilson_dispersion":{
            "singular_value_at_accidental_angle":euclidean_at_accidental,
            "singular_value_at_pi":euclidean_at_pi,
            "expected_at_pi":2.0*r/a},
        "finite_matrix_gamma5_hermiticity":{
            "naive_chiral_anticommutator_residual":naive_chiral,
            "euclidean_wilson_gamma5_hermiticity_residual":euclidean_gamma5_herm,
            "euclidean_wilson_ordinary_hermiticity_residual":ordinary_herm,
            "assembler_report":de_audit},
        "passed":bool(abs(p0)<1e-11 and abs(ppi-2.0)<1e-11
            and audit["minimum_eigenvalue"]>-1e-11
            and abs(scalar_low)<1e-10 and scalar_high>0
            and euclidean_at_accidental>1e-6
            and abs(euclidean_at_pi-2.0)<1e-10
            and naive_chiral<1e-11 and euclidean_gamma5_herm<1e-11
            and ordinary_herm>1e-6),
        "scope":"finite graph algebra and exact 1D Fourier-symbol convention audit; no Q6 continuum no-doubling theorem"
    }
    out=Path("results/json/wilson_fermion_audit.json")
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))
    if not report["passed"]:
        raise SystemExit("Wilson audit failed")
        raise SystemExit("Euclidean Wilson convention audit failed")


if __name__=="__main__":
    main()

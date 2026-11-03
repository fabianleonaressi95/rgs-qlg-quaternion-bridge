"""Reproducible finite-graph Wilson-term audit.

Run from the repository root:
    python examples/wilson_fermion_audit.py
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from src.discrete_dirac import (
    assemble_spatial_dirac, assemble_wilson_laplacian, chirality_matrix,
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
    }
    out=Path("results/json/wilson_fermion_audit.json")
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))
    if not report["passed"]:
        raise SystemExit("Wilson audit failed")

if __name__=="__main__":
    main()

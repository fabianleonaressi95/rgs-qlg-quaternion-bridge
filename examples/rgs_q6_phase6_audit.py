"""Reproducible Phase 6 audit on the RGS/Q6 NPZ archive.

Usage:
  python examples/rgs_q6_phase6_audit.py /path/to/RGS_PSU_all_matrices.npz
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from src.discrete_dirac import curved_clifford_audit, polar_spin_links, assemble_spatial_dirac, discrete_dirac_spectrum

ELL2=7.035514834154695
C0=43.96495130838031

def reconstruct_vierbein(x):
    delta=x[None,:,:]-x[:,None,:]
    r2=np.einsum("ijk,ijk->ij",delta,delta)
    weights=np.exp(-r2/ELL2)
    np.fill_diagonal(weights,0.0)
    c=np.einsum("ij,ija,ijb->iab",weights,delta,delta)
    e=np.empty_like(c)
    for i,ci in enumerate(c):
        lam,u=np.linalg.eigh(ci/C0)
        if np.min(lam)<=0: raise ValueError(f"non-positive local covariance at vertex {i}")
        e[i]=(u*np.sqrt(lam))@u.T
    return e

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("archive",type=Path)
    parser.add_argument("--json",type=Path,default=Path("results/json/rgs_q6_phase6_audit.json"))
    parser.add_argument("--matrix",type=Path,default=None)
    args=parser.parse_args()
    with np.load(args.archive,allow_pickle=False) as z:
        x=np.asarray(z["X_Oh_reconstructed"],float)
        adjacency=np.asarray(z["Q6_adjacency"])
    if x.shape!=(64,3) or adjacency.shape!=(64,64): raise ValueError("archive does not contain the expected RGS/Q6 arrays")
    edges=[tuple(map(int,p)) for p in np.argwhere(np.triu(adjacency,1)>0)]
    if len(edges)!=192: raise ValueError(f"expected 192 Q6 edges, found {len(edges)}")
    e=reconstruct_vierbein(x)
    cliff=curved_clifford_audit(e)
    links,polar_res=polar_spin_links(e,edges)
    d,da=assemble_spatial_dirac(x,e,edges,spin_links=links)
    spec=discrete_dirac_spectrum(d)
    values=spec.pop("eigenvalues")
    unique,multiplicity=np.unique(np.round(values,10),return_counts=True)
    report={
      "audit":"RGS_Q6_WEILER_PHASE6",
      "input_archive":args.archive.name,
      "constants":{"ell_squared":ELL2,"c0":C0},
      "complex":{"vertices":64,"edges":len(edges),"degree_min":int(adjacency.sum(axis=1).min()),"degree_max":int(adjacency.sum(axis=1).max())},
      "vierbein":{"shape":list(e.shape),"det_min":float(np.linalg.det(e).min()),"det_max":float(np.linalg.det(e).max()),"condition_max":float(max(np.linalg.cond(v) for v in e))},
      "curved_clifford":cliff,
      "polar_spin_transport":{"links":len(links),"max_SO3_orthogonality_residual":polar_res},
      "dirac":{k:v for k,v in da.items()},
      "spectrum":{**spec,"distinct_rounded_levels":len(unique),"levels":[{"value":float(v),"multiplicity":int(m)} for v,m in zip(unique,multiplicity)]},
      "interpretation":"NUMERICALLY VERIFIED finite-graph candidate only; not a continuum Dirac or Mosco-convergence result"
    }
    args.json.parent.mkdir(parents=True,exist_ok=True)
    args.json.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    if args.matrix is not None:
        args.matrix.parent.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(args.matrix,vierbein=e,dirac=d,edges=np.asarray(edges,dtype=int),eigenvalues=values)
    print(json.dumps(report,indent=2))

if __name__=="__main__": main()

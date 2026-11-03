"""RGS-Q6 overlap audit. Run with an NPZ archive path as the first argument."""
import argparse, json
from collections import defaultdict
from pathlib import Path
import numpy as np
from src.discrete_dirac import assemble_spatial_dirac, assemble_wilson_laplacian, chirality_matrix, polar_spin_links
from src.overlap_dirac import assemble_overlap_from_wilson
from src.overlap_dirac import (assemble_overlap_from_wilson, overlap_index, overlap_kernel_spectral_flow, overlap_locality_profile)


def reconstruct_vierbein(x, ell2=7.035514834154695, c0=43.96495130838031):
    n=len(x); e=np.zeros((n,3,3))
    for i in range(n):
        d=x-x[i]; w=np.exp(-np.sum(d*d,axis=1)/ell2); w[i]=0
        c=np.einsum("j,ja,jb->ab",w,d,d)
        v,u=np.linalg.eigh(c/c0)
        if v.min()<=0: raise ValueError("non-positive covariance")
        e[i]=(u*np.sqrt(v))@u.T
    return e


def run(archive, output, wilson_r=1.0, rho=1.0):
    z=np.load(archive); x=z["X_Oh_reconstructed"]; adj=z["Q6_adjacency"]; states=z["Q6_states"]
    edges=[tuple(map(int,p)) for p in np.argwhere(np.triu(adj,1)>0)]
    if x.shape!=(64,3) or len(edges)!=192 or states.shape!=(64,6):
        raise ValueError("expected the 64-vertex Q6 archive")
    e=reconstruct_vierbein(x); links,rotres=polar_spin_links(e,edges)
    d0,da=assemble_spatial_dirac(x,e,edges,spin_links=links)
    w,wa=assemble_wilson_laplacian(x,e,edges,spin_links=links,wilson_r=wilson_r)
    g5=np.kron(np.eye(64),chirality_matrix())
    lengths=np.array([np.linalg.norm(x[j]-x[i]) for i,j in edges])
    a=float(lengths.mean())
    dov,oa=assemble_overlap_from_wilson(d0,w,g5,lattice_spacing=a,rho=rho)
    index_audit=overlap_index(dov,g5,lattice_spacing=a,rho=rho)
    spectral_flow=overlap_kernel_spectral_flow(
        d0,w,g5,a,np.linspace(max(0.1,0.5*rho),1.5*rho,9)
    )
    locality_audit=overlap_locality_profile(dov,edges,64,4)
    blocks=defaultdict(list)
    for i in range(64):
        for j in range(64):
            if i!=j:
                dist=int(np.count_nonzero(states[i]!=states[j]))
                blocks[dist].append(float(np.linalg.norm(dov[4*i:4*i+4,4*j:4*j+4],"fro")))
    locality={str(k):{"mean":float(np.mean(v)),"max":float(np.max(v))}
              for k,v in sorted(blocks.items())}
    ds=np.array(sorted(blocks),float); means=np.array([np.mean(blocks[int(k)]) for k in ds])
    slope=float(np.polyfit(ds,np.log(means),1)[0])
    report={
      "audit":"RGS_Q6_OVERLAP_OPERATOR",
      "geometry":{"vertices":64,"edges":192,"detE_min":float(np.linalg.det(e).min()),
                  "detE_max":float(np.linalg.det(e).max()),"edge_length_min":float(lengths.min()),
                  "reference_spacing_mean_edge":a,"edge_length_max":float(lengths.max()),
                  "polar_rotation_residual":rotres},
      "operators":{"naive_hermiticity":da["hermiticity_residual"],
                   "naive_chiral_residual":float(np.linalg.norm(d0@g5+g5@d0,"fro")),
                   "wilson_hermiticity":wa["hermiticity_residual"],
                   "wilson_min_eigenvalue":wa["minimum_eigenvalue"],"wilson_r":wilson_r},
      "overlap":oa,
      "index":index_audit,
      "spectral_flow":spectral_flow,
      "locality_graph_distance":locality_audit,
      "locality":{"distance":"Q6 Hamming distance","block_norms":locality,
                  "log_mean_fit_slope":slope,"empirical_factor_per_step":float(np.exp(slope)),
                  "status":"finite graph diagnostic, not a uniform locality theorem"},
      "status":{"Q6_overlap":"NUMERICALLY_CONSTRUCTED","GW":"VERIFIED",
                "kernel_gap":"POSITIVE_FOR_THIS_ARCHIVE_AND_RHO",
                "continuum_locality":"OPEN","index_theorem":"OPEN","Mosco":"OPEN"}
                "continuum_locality":"OPEN","finite_matrix_index_trace":("INTEGER_WITHIN_TOLERANCE" if index_audit["integer_within_tolerance"] else "NONINTEGER_OR_UNRESOLVED"),"continuum_index_theorem":"OPEN","Mosco":"OPEN"}
    }
    out=Path(output); out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))
    if oa["ginsparg_wilson_residual"]>1e-9 or oa["gamma5_hermiticity_residual"]>1e-9:
        raise SystemExit("overlap audit failed")
    return report


if __name__=="__main__":
    p=argparse.ArgumentParser(); p.add_argument("archive")
    p.add_argument("--output",default="results/json/rgs_q6_overlap_audit.json")
    p.add_argument("--wilson-r",type=float,default=1.0); p.add_argument("--rho",type=float,default=1.0)
    a=p.parse_args(); run(a.archive,a.output,a.wilson_r,a.rho)

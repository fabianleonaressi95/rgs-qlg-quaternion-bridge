"""Finite compact U(1) quotient and Wilson sector on Q6.

Gauge convention: theta' = theta - D0 lambda (mod 2pi). A deterministic
spanning tree fixes 63 links; the 129 chord angles carry normalized quotient
Haar measure. This finite construction is not a derivation of Wyler's alpha.
"""
from __future__ import annotations
import math
import numpy as np
from src.rgs_q6_equivariant_cochains import build_q6_cochain_complex, induced_cochain_actions
from src.rgs_icosahedral_q6_intertwiner import icosahedral_signed_axis_representation


def q6_u1_quotient_data(complex_data=None):
    c=build_q6_cochain_complex() if complex_data is None else complex_data
    parent=list(range(64))
    def find(x):
        while parent[x]!=x:
            parent[x]=parent[parent[x]]; x=parent[x]
        return x
    tree=[]; chords=[]
    for k,(a,b,_d) in enumerate(c["edges"]):
        x,y=find(a),find(b)
        if x!=y: parent[x]=y; tree.append(k)
        else: chords.append(k)
    d0=np.asarray(c["D0"],dtype=float); d1=np.asarray(c["D1"],dtype=float)
    reduced=d1[:,chords]; rank=int(np.linalg.matrix_rank(reduced))
    if len(tree)!=63 or len(chords)!=129 or rank!=129:
        raise ArithmeticError("unexpected Q6 gauge quotient topology")
    return {"complex":c,"tree_edges":tuple(tree),"chord_edges":tuple(chords),
            "D0":d0,"D1":d1,"D1_chords":reduced,"quotient_dimension":129,
            "plaquette_map_rank":rank}


def gauge_transform_u1(theta, gauge_angles, quotient_data=None):
    q=q6_u1_quotient_data() if quotient_data is None else quotient_data
    t=np.asarray(theta,dtype=float); lam=np.asarray(gauge_angles,dtype=float)
    if t.shape!=(192,) or lam.shape!=(64,) or not np.all(np.isfinite(t)) or not np.all(np.isfinite(lam)):
        raise ValueError("expected finite theta (192,) and gauge angles (64,)")
    return (t-q["D0"]@lam+math.pi)%(2*math.pi)-math.pi


def tree_gauge_fix(theta, quotient_data=None):
    q=q6_u1_quotient_data() if quotient_data is None else quotient_data
    t=np.asarray(theta,dtype=float); c=q["complex"]
    if t.shape!=(192,) or not np.all(np.isfinite(t)): raise ValueError("theta must be finite shape (192,)")
    adj=[[] for _ in range(64)]
    for e in q["tree_edges"]:
        a,b,_=c["edges"][e]
        adj[a].append((b,e,1.0)); adj[b].append((a,e,-1.0))
    lam=np.full(64,np.nan); lam[0]=0.; stack=[0]
    while stack:
        v=stack.pop()
        for w,e,sgn in adj[v]:
            if np.isfinite(lam[w]): continue
            lam[w]=lam[v]+sgn*t[e]; stack.append(w)
    fixed=(t-q["D0"]@lam+math.pi)%(2*math.pi)-math.pi
    fixed[np.abs(fixed)<5e-14]=0.
    if np.max(np.abs(fixed[list(q["tree_edges"])]))>1e-11: raise ArithmeticError("tree gauge fix failed")
    return fixed,lam


def plaquette_fluxes(theta, quotient_data=None):
    q=q6_u1_quotient_data() if quotient_data is None else quotient_data
    t=np.asarray(theta,dtype=float)
    if t.shape!=(192,) or not np.all(np.isfinite(t)): raise ValueError("theta must be finite shape (192,)")
    return q["D1"]@t


def u1_wilson_action(theta,beta=1.,quotient_data=None):
    if not np.isfinite(beta) or beta<=0: raise ValueError("beta must be finite and positive")
    return float(beta*np.sum(1.-np.cos(plaquette_fluxes(theta,quotient_data))))


def u1_quadratic_action(theta,beta=1.,quotient_data=None):
    if not np.isfinite(beta) or beta<=0: raise ValueError("beta must be finite and positive")
    f=plaquette_fluxes(theta,quotient_data)
    return float(.5*beta*np.dot(f,f))


def sample_quotient_haar(count=1,seed=0):
    if not isinstance(count,(int,np.integer)) or count<1: raise ValueError("count must be positive")
    return np.random.default_rng(seed).uniform(-math.pi,math.pi,size=(int(count),129))


def lift_chord_angles(chord_angles,quotient_data=None):
    q=q6_u1_quotient_data() if quotient_data is None else quotient_data
    x=np.asarray(chord_angles,dtype=float)
    if x.shape[-1:]!=(129,) or not np.all(np.isfinite(x)): raise ValueError("final dimension must be 129")
    out=np.zeros(x.shape[:-1]+(192,),dtype=float); out[...,list(q["chord_edges"])]=x
    return out


def audit_q6_u1_quotient(seed=20261104):
    q=q6_u1_quotient_data(); c=q["complex"]; rng=np.random.default_rng(seed)
    theta=rng.uniform(-.7,.7,192); lam=rng.uniform(-math.pi,math.pi,64); beta=2.3
    moved=gauge_transform_u1(theta,lam,q); fixed,_=tree_gauge_fix(theta,q); fixed2,_=tree_gauge_fix(moved,q)
    action=u1_wilson_action(theta,beta,q)
    flux_res=float(np.max(np.abs(plaquette_fluxes(theta,q)-plaquette_fluxes(moved,q))))
    orbit_res=float(np.max(np.abs(np.exp(1j*fixed[list(q["chord_edges"])])-np.exp(1j*fixed2[list(q["chord_edges"])]))))
    tree_res=float(np.max(np.abs(fixed[list(q["tree_edges"])])))
    singular=np.linalg.svd(q["D1_chords"],compute_uv=False)
    small=rng.normal(size=192); eps=1e-4
    exact=u1_wilson_action(eps*small,beta,q); quad=u1_quadratic_action(eps*small,beta,q)
    small_err=abs(exact-quad)/max(quad,np.finfo(float).tiny)
    sym_res=0.
    group=icosahedral_signed_axis_representation()
    for s in group["signed_permutations"]:
        _p0,p1,p2=induced_cochain_actions(s,c); moved_theta=p1@theta
        sym_res=max(sym_res,float(np.max(np.abs(q["D1"]@moved_theta-p2@(q["D1"]@theta)))),
                    abs(u1_wilson_action(moved_theta,beta,q)-action))
    sample=sample_quotient_haar(4,seed+1)
    checks={"cell_counts":len(c["edges"])==192 and len(c["faces"])==240,
      "tree_chords":len(q["tree_edges"])==63 and len(q["chord_edges"])==129,
      "D1D0_zero":bool(np.max(np.abs(q["D1"]@q["D0"]))==0),
      "reduced_rank":q["plaquette_map_rank"]==129,"tree_gauge":tree_res<1e-11,
      "orbit_representative":orbit_res<1e-10,"flux_gauge_invariance":flux_res<1e-10,
      "wilson_gauge_invariance":abs(action-u1_wilson_action(moved,beta,q))<1e-10,
      "icosahedral_invariance":sym_res<1e-10,"haar_shape":sample.shape==(4,129),
      "small_field_limit":small_err<1e-7}
    return {"audit":"Q6_U1_GAUGE_QUOTIENT_HAAR_WILSON","vertices":64,"edges":192,"plaquettes":240,
      "tree_edges_gauge_fixed":63,"quotient_dimension":129,"reduced_plaquette_map_rank":q["plaquette_map_rank"],
      "reduced_plaquette_map_min_singular_value":float(singular[-1]),
      "normalized_haar_measure":"product_{129 chords} dtheta/(2*pi)",
      "wilson_action":"beta * sum_f (1-cos((D1 theta)_f))",
      "quadratic_limit":"beta/2 * ||D1 theta||^2","gauge_residual":abs(action-u1_wilson_action(moved,beta,q)),
      "tree_gauge_residual":tree_res,"orbit_residual":orbit_res,"icosahedral_residual":sym_res,
      "small_field_quadratic_relative_error":small_err,"checks":checks,"all_checks_passed":bool(all(checks.values())),
      "not_established":["D5 identification","10-real-dimensional effective sector","physical charge normalization",
        "four-dimensional continuum/QED limit","fine-structure constant prediction"]}

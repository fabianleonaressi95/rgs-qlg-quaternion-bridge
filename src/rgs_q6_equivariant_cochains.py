"""Icosahedrally equivariant cochain complex on the cubical Q6 graph."""
from __future__ import annotations
import numpy as np
from src.rgs_icosahedral_q6_intertwiner import icosahedral_signed_axis_representation, q6_vertex_action

def _edge_key(a, b):
    return (a, b) if a < b else (b, a)

def build_q6_cochain_complex():
    masks = tuple(1 << (5-d) for d in range(6))
    edges, edge_index = [], {}
    for v in range(64):
        for d, mask in enumerate(masks):
            w = v ^ mask
            if v < w:
                edge_index[(v,w)] = len(edges); edges.append((v,w,d))
    faces, face_index = [], {}
    for i in range(6):
        for j in range(i+1,6):
            for v in range(64):
                if v & masks[i] or v & masks[j]: continue
                cyc=(v,v^masks[i],v^masks[i]^masks[j],v^masks[j])
                face_index[frozenset(cyc)] = len(faces); faces.append((v,i,j,cyc))
    d0=np.zeros((192,64),dtype=np.int8)
    for e,(a,b,_d) in enumerate(edges): d0[e,a]=-1; d0[e,b]=1
    d1=np.zeros((240,192),dtype=np.int8)
    for f,(_v,_i,_j,cyc) in enumerate(faces):
        for a,b in zip(cyc,cyc[1:]+cyc[:1]):
            d1[f,edge_index[_edge_key(a,b)]] += 1 if a<b else -1
    return {"edges":tuple(edges),"faces":tuple(faces),"edge_index":edge_index,
            "face_index":face_index,"D0":d0,"D1":d1}

def _boundary_vector(cycle, edge_index):
    out=np.zeros(192,dtype=np.int8)
    for a,b in zip(cycle,cycle[1:]+cycle[:1]):
        out[edge_index[_edge_key(a,b)]] += 1 if a<b else -1
    return out

def induced_cochain_actions(signed_permutation, complex_data=None):
    c=build_q6_cochain_complex() if complex_data is None else complex_data
    s=np.asarray(signed_permutation)
    p0=np.zeros((64,64),dtype=np.int8); p1=np.zeros((192,192),dtype=np.int8); p2=np.zeros((240,240),dtype=np.int8)
    vi=[q6_vertex_action(v,s) for v in range(64)]
    for v,w in enumerate(vi): p0[w,v]=1
    for e,(a,b,_d) in enumerate(c["edges"]):
        x,y=vi[a],vi[b]; p1[c["edge_index"][_edge_key(x,y)],e]=1 if x<y else -1
    for f,(_v,_i,_j,cyc) in enumerate(c["faces"]):
        image=tuple(vi[x] for x in cyc)
        target=c["face_index"].get(frozenset(image))
        if target is None: raise ArithmeticError("action does not preserve square faces")
        ib=_boundary_vector(image,c["edge_index"])
        tb=_boundary_vector(c["faces"][target][3],c["edge_index"])
        orient=int(round(float(ib@tb)/4.0))
        if orient not in (-1,1): raise ArithmeticError("invalid induced face orientation")
        p2[target,f]=orient
    return p0,p1,p2

def audit_q6_equivariant_cochains():
    c=build_q6_cochain_complex(); d0,d1=c["D0"],c["D1"]
    r0=int(np.linalg.matrix_rank(d0.astype(float))); r1=int(np.linalg.matrix_rank(d1.astype(float)))
    d2=float(np.max(np.abs(d1@d0))); e0=e1=0.0; failures=0
    group=icosahedral_signed_axis_representation()
    for s in group["signed_permutations"]:
        p0,p1,p2=induced_cochain_actions(s,c)
        e0=max(e0,float(np.max(np.abs(d0@p0-p1@d0))))
        e1=max(e1,float(np.max(np.abs(d1@p1-p2@d1))))
        if (not np.array_equal(p0.T@p0,np.eye(64,dtype=np.int8))
            or not np.array_equal(p1.T@p1,np.eye(192,dtype=np.int8))
            or not np.array_equal(p2.T@p2,np.eye(240,dtype=np.int8))): failures+=1
    return {"n_vertices":64,"n_edges":192,"n_square_faces":240,"rank_D0":r0,"rank_D1":r1,
      "d1_d0_max_abs":d2,"group_order":int(group["group_order"]),
      "D0_equivariance_max_abs":e0,"D1_equivariance_max_abs":e1,
      "cochain_orthogonality_failures":failures,
      "all_checks_passed":bool(r0==63 and r1==129 and d2==0 and e0==0 and e1==0 and failures==0)}

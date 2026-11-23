"""Exact elliptic Jitterbug vertex and rigid-face pose audit.

Uses scalar-first Hamilton quaternion conventions documented in
docs/quaternion-conventions.md. The symmetric vertex trajectory is parameterized
by ellipse angle theta; face rotation gamma is a distinct angle.
"""
from __future__ import annotations
from itertools import product
import math
import numpy as np
from scipy.spatial.transform import Rotation

def vertices0():
    V=[]
    for zero in range(3):
        axes=[i for i in range(3) if i!=zero]
        for s1,s2 in product((-1.0,1.0),repeat=2):
            p=np.zeros(3); p[axes]=[s1,s2]; V.append(p/math.sqrt(2))
    return np.asarray(V)

def faces0(V):
    out=[]
    for eps in product((-1,1),repeat=3):
        ids=[]
        for zero in range(3):
            p=np.zeros(3); axes=[i for i in range(3) if i!=zero]
            p[axes]=[eps[i] for i in axes]; p/=math.sqrt(2)
            ids.append(int(np.argmin(np.linalg.norm(V-p,axis=1))))
        out.append(ids)
    return out

def edges0(V):
    return [(i,j) for i in range(12) for j in range(i+1,12)
            if abs(np.linalg.norm(V[i]-V[j])-1)<1e-10]

def vertex_path(theta):
    V=vertices0(); c,s=math.cos(theta),math.sin(theta)
    den=math.sqrt(3-2*c*c); Z,Y=c/den,s/den
    out=np.empty_like(V)
    for i,p in enumerate(V):
        k=int(np.argmin(np.abs(p))); axes=[j for j in range(3) if j!=k]
        u=p.copy(); v=np.cross(np.eye(3)[k],u); v/=np.linalg.norm(v)
        sigma=-float(np.sign(p[axes[0]]*p[axes[1]]))
        out[i]=Z*u+sigma*Y*v
    return out

def face_pose(eps,theta):
    gamma=math.atan(math.sqrt(3)*math.tan(theta))
    u=np.asarray(eps,dtype=float)/math.sqrt(3)
    sign=-int(np.prod(eps))
    q=np.r_[math.cos(gamma/2),sign*u*math.sin(gamma/2)] # scalar-first
    R=Rotation.from_quat([q[1],q[2],q[3],q[0]]).as_matrix()
    t=math.sqrt(2/3)*(math.cos(gamma)-1)*u
    return q,R,t,gamma

def audit(samples=101,tol=1e-12):
    V=vertices0(); F=faces0(V); E=edges0(V)
    max_edge=max_face=max_joint=0.0
    for theta in np.linspace(-math.pi/4,math.pi/4,samples):
        P=vertex_path(float(theta))
        max_edge=max(max_edge,max(abs(np.linalg.norm(P[i]-P[j])-1) for i,j in E)
                    )
        estimates=[{} for _ in range(12)]
        for eps,face in zip(product((-1,1),repeat=3),F):
            q,R,t,gamma=face_pose(eps,float(theta))
            fit=V[face]@R.T+t
            max_face=max(max_face,float(np.max(np.linalg.norm(fit-P[face],axis=1))))
            for idx,point in zip(face,fit): estimates[idx].setdefault(len(estimates[idx]),point)
        for vals in estimates:
            pts=list(vals.values())
            if len(pts)>1:
                max_joint=max(max_joint,max(np.linalg.norm(x-pts[0]) for x in pts))
    if max(max_edge,max_face,max_joint)>tol: raise AssertionError((max_edge,max_face,max_joint))
    return dict(vertices=12,edges=24,triangular_faces=8,samples=samples,
                max_edge_residual=max_edge,max_rigid_face_residual=max_face,
                max_joint_coincidence_residual=max_joint)

if __name__=="__main__":
    print(audit())

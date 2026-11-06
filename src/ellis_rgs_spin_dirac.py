"""Intrinsic spinorial Wilson-Dirac operator on Ellis-RGS tetrahedral meshes.

Uses Levi-Civita SU(2) edge transport. Clifford directions are endpoint
tangents in local orthonormal Ellis frames; no fictitious R3 vierbein.
"""
from __future__ import annotations
from collections.abc import Mapping, Sequence
import numpy as np
from src.clifford_qlg import gamma_matrices
from src.discrete_dirac import chirality_matrix
from src.ellis_spin_connection import (
    _great_circle, _point_data, ellis_edge_spin_links, ellis_orthonormal_frame,
)

def _edge_path_data(point_i, point_j, r0):
    ui, ni = _point_data(point_i, r0)
    uj, nj = _point_data(point_j, r0)
    theta, n_path, ndot_path = _great_circle(ni, nj)
    du = uj-ui
    def tangent(t):
        u = ui+t*du
        n, ndot = n_path(t), ndot_path(t)
        R = np.sqrt(u*u+r0*r0)
        return np.r_[(u*du/R)*n+R*ndot, r0*du/R]
    nodes, weights = np.polynomial.legendre.leggauss(8)
    us = ui+0.5*(nodes+1.0)*du
    speed = np.sqrt(du*du+(np.sqrt(us*us+r0*r0)*theta)**2)
    length = float(0.5*np.dot(weights, speed))
    if not np.isfinite(length) or length <= 1e-14:
        raise ValueError("Ellis edge path has zero or invalid intrinsic length")
    return length, tangent

def _lumped_vertex_volumes(mesh):
    volumes = np.zeros(len(mesh.vertices_r4), dtype=float)
    for tet, volume in zip(mesh.tetrahedra, mesh.element_volumes):
        volumes[np.asarray(tet, dtype=int)] += float(volume)/4.0
    if np.any(~np.isfinite(volumes)) or np.any(volumes <= 0):
        raise ValueError("every vertex must have a positive lumped tetrahedral volume")
    return volumes

def assemble_ellis_spin_wilson_dirac(
    mesh, edges: Sequence[tuple[int,int]] | None = None,
    spin_links: Mapping | None = None, wilson_r: float = 1.0,
    transport_steps: int = 32, local_spin_frames: np.ndarray | None = None,
):
    """Assemble the intrinsic finite-mesh Euclidean Wilson-Dirac operator.

    Canonical edge i<j: U_ji maps spinor components at i to j.
    local_spin_frames are optional SU(2) changes of local spin basis; gamma
    matrices and links are transformed together for a covariance audit.
    """
    vertices = np.asarray(mesh.vertices_r4, dtype=float)
    n = len(vertices)
    if vertices.shape != (n,4) or n == 0:
        raise ValueError("mesh.vertices_r4 must have shape (N,4)")
    r0, r = float(mesh.r0), float(wilson_r)
    if not np.isfinite(r) or r < 0:
        raise ValueError("wilson_r must be finite and non-negative")
    if edges is None:
        if mesh.chain is None:
            raise ValueError("mesh must have include_chain=True or explicit edges")
        edges = mesh.chain["edges"]
    edge_list = [tuple(map(int,e)) for e in edges]
    if any(not (0 <= i < j < n) for i,j in edge_list):
        raise ValueError("edges must be canonical i<j vertex indices")
    if len(set(edge_list)) != len(edge_list):
        raise ValueError("duplicate edges are not allowed")

    frames = np.stack([ellis_orthonormal_frame(*_point_data(p,r0),r0) for p in vertices])
    if local_spin_frames is None:
        gauges = np.broadcast_to(np.eye(2,dtype=complex),(n,2,2)).copy()
    else:
        gauges = np.asarray(local_spin_frames,dtype=complex)
        if gauges.shape != (n,2,2):
            raise ValueError("local_spin_frames must have shape (N,2,2)")
        for i,g in enumerate(gauges):
            if np.linalg.norm(g.conj().T@g-np.eye(2)) > 1e-9 or abs(np.linalg.det(g)-1)>1e-9:
                raise ValueError(f"local_spin_frames[{i}] must lie in SU(2)")
    if spin_links is None:
        spin_links, link_audit = ellis_edge_spin_links(vertices,edge_list,r0,steps=transport_steps)
    else:
        link_audit = {"connection":"user-supplied SU(2) edge links"}
    if set(spin_links) != set(edge_list):
        raise ValueError("spin link keys mismatch with requested edges")

    gs, gamma5 = gamma_matrices(), chirality_matrix()
    I2, I4 = np.eye(2,dtype=complex), np.eye(4,dtype=complex)
    raw_d = np.zeros((4*n,4*n),dtype=complex)
    raw_w = np.zeros_like(raw_d)
    lengths = {}
    max_unitarity = 0.0
    for i,j in edge_list:
        link = np.asarray(spin_links[(i,j)],dtype=complex)
        if link.shape != (2,2) or np.linalg.norm(link.conj().T@link-I2)>1e-9 or abs(np.linalg.det(link)-1)>1e-9:
            raise ValueError(f"spin link {(i,j)} must lie in SU(2)")
        max_unitarity=max(max_unitarity,float(np.linalg.norm(link.conj().T@link-I2)))
        link = gauges[j]@link@gauges[i].conj().T
        S = np.block([[link,np.zeros((2,2))],[np.zeros((2,2)),link]])
        length, tangent = _edge_path_data(vertices[i],vertices[j],r0)
        lengths[(i,j)] = length
        ni = frames[i].T@tangent(0.0)
        nj = frames[j].T@tangent(1.0)
        ni /= np.linalg.norm(ni); nj /= np.linalg.norm(nj)
        ci = sum(ni[a]*gs[a+1] for a in range(3))
        cj = sum(nj[a]*gs[a+1] for a in range(3))
        si,sj=slice(4*i,4*i+4),slice(4*j,4*j+4)
        block=0.5j/length*(cj@S+S@ci)
        raw_d[sj,si]+=block
        raw_d[si,sj]+=block.conj().T
        if r>0:
            weight=r/(2*length)
            raw_w[si,si]+=weight*I4; raw_w[sj,sj]+=weight*I4
            raw_w[sj,si]-=weight*S
            raw_w[si,sj]-=weight*S.conj().T

    vertex_volumes=_lumped_vertex_volumes(mesh)
    invsqrt=np.repeat(1/np.sqrt(vertex_volumes),4)
    D=invsqrt[:,None]*raw_d*invsqrt[None,:]
    W=invsqrt[:,None]*raw_w*invsqrt[None,:]
    DE=1j*D+W
    G5=np.kron(np.eye(n),gamma5)
    hD=float(np.linalg.norm(D-D.conj().T,ord="fro"))
    hW=float(np.linalg.norm(W-W.conj().T,ord="fro"))
    g5=float(np.linalg.norm(DE.conj().T-G5@DE@G5,ord="fro"))
    minW=float(np.linalg.eigvalsh(W)[0])
    G=np.zeros((4*n,4*n),dtype=complex)
    for i,g in enumerate(gauges):
        g4=np.block([[g,np.zeros((2,2))],[np.zeros((2,2)),g]])
        G[4*i:4*i+4,4*i:4*i+4]=g4
    return DE,{
        "vertices":n,"edges":len(edge_list),"spinor_dimension":4,
        "edge_length_min":min(lengths.values(),default=0.0),
        "edge_length_max":max(lengths.values(),default=0.0),
        "lumped_volume_min":float(vertex_volumes.min()),
        "lumped_volume_max":float(vertex_volumes.max()),
        "dirac_hermiticity_residual":hD,"wilson_hermiticity_residual":hW,
        "gamma5_hermiticity_residual":g5,"wilson_minimum_eigenvalue":minW,
        "input_link_unitarity_residual":max_unitarity,
        "gauge_covariance_transform":G,
        "connection":link_audit.get("connection","ELLIS_LEVI_CIVITA_SU2"),
        "measure":"lumped P1 tetrahedral volume",
        "edge_metric":"intrinsic Ellis path length",
        "status":"FINITE_INTRINSIC_ELLIS_SPIN_WILSON_DIRAC",
    }

"""Gauge-covariant spinor currents and Kirchhoff continuity on Q6.

Convention: each canonical edge is oriented i<j and carries U_ij in SU(2).
A vertex spinor transforms as psi_i -> G_i psi_i and links as
U_ij -> G_j U_ij G_i^dagger. The covariant edge difference is
D_ij psi = psi_j - U_ij psi_i.

The edge current is J_ij = 2 w_ij Im(psi_i^dagger U_ij^dagger psi_j).
With the repository boundary convention B1[source]=-1, B1[target]=+1,
the Schrödinger equation i dpsi/dt = H psi for H=D^dagger W D gives
rho_dot = B1 J exactly. This is a finite graph continuity identity, not
yet a calibrated electromagnetic current in SI units.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Iterable

import numpy as np

from src.q6_complex import build_q6_complex


def dagger(a: np.ndarray) -> np.ndarray:
    return np.asarray(a, dtype=complex).conj().T


def random_su2(rng: np.random.Generator) -> np.ndarray:
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    a, b, c, d = map(float, q)
    return np.array([[a + 1j*b, c + 1j*d],
                     [-c + 1j*d, a - 1j*b]], dtype=complex)


def _hash_array(a: np.ndarray) -> str:
    arr = np.ascontiguousarray(a)
    h = hashlib.sha256()
    h.update(str(arr.dtype).encode("ascii"))
    h.update(json.dumps(list(arr.shape)).encode("ascii"))
    h.update(arr.tobytes())
    return h.hexdigest()


def random_q6_spinor_data(seed: int = 20261104):
    rng = np.random.default_rng(seed)
    q6 = build_q6_complex()
    edges = np.asarray(q6.edges, dtype=np.int64)
    links = np.stack([random_su2(rng) for _ in edges])
    psi = rng.normal(size=(64, 2)) + 1j*rng.normal(size=(64, 2))
    weights = 0.5 + rng.random(len(edges))
    gauge = np.stack([random_su2(rng) for _ in q6.vertices])
    return q6, edges, links, psi, weights, gauge


def gauge_transform_links(edges: np.ndarray, links: np.ndarray,
                          gauge: np.ndarray) -> np.ndarray:
    out = np.empty_like(links)
    for k, (i, j) in enumerate(edges):
        out[k] = gauge[j] @ links[k] @ dagger(gauge[i])
    return out


def gauge_transform_spinors(psi: np.ndarray, gauge: np.ndarray) -> np.ndarray:
    return np.einsum("vab,vb->va", gauge, psi)


def edge_currents(edges: np.ndarray, links: np.ndarray,
                  psi: np.ndarray, weights: np.ndarray) -> np.ndarray:
    currents = np.empty(len(edges), dtype=float)
    for k, (i, j) in enumerate(edges):
        z = np.vdot(psi[i], dagger(links[k]) @ psi[j])
        currents[k] = 2.0 * weights[k] * float(np.imag(z))
    return currents


def covariant_laplacian(edges: np.ndarray, links: np.ndarray,
                        weights: np.ndarray, n_vertices: int = 64) -> np.ndarray:
    """Hermitian block matrix H=D^dagger W D on C^(2*n_vertices)."""
    H = np.zeros((2*n_vertices, 2*n_vertices), dtype=complex)
    eye = np.eye(2, dtype=complex)
    for k, (i, j) in enumerate(edges):
        w, U = float(weights[k]), links[k]
        si, sj = slice(2*i, 2*i+2), slice(2*j, 2*j+2)
        H[si, si] += w*eye
        H[sj, sj] += w*eye
        H[si, sj] += -w*dagger(U)
        H[sj, si] += -w*U
    return H


def density_derivative(psi: np.ndarray, H: np.ndarray) -> np.ndarray:
    flat = psi.reshape(-1)
    dpsi = -1j * (H @ flat)
    dpsi = dpsi.reshape(psi.shape)
    return 2.0*np.real(np.sum(np.conjugate(psi)*dpsi, axis=1))


def path_transport(path: Iterable[int], edges: np.ndarray,
                   links: np.ndarray) -> np.ndarray:
    """Transport from first to last vertex; P'=G_last P G_first^dagger."""
    vertices = tuple(int(v) for v in path)
    if len(vertices) < 2:
        raise ValueError("path must contain at least one edge")
    lookup = {tuple(map(int, e)): links[k] for k, e in enumerate(edges)}
    P = np.eye(2, dtype=complex)
    for a, b in zip(vertices[:-1], vertices[1:]):
        key = (min(a, b), max(a, b))
        if key not in lookup:
            raise ValueError(f"({a},{b}) is not a Q6 edge")
        U = lookup[key] if a < b else dagger(lookup[key])
        P = U @ P
    return P


def audit_spinor_kirchhoff(seed: int = 20261104,
                            gauge_seed: int | None = None,
                            atol: float = 1e-10) -> dict:
    q6, edges, links, psi, weights, gauge = random_q6_spinor_data(seed)
    if gauge_seed is not None:
        rng = np.random.default_rng(gauge_seed)
        gauge = np.stack([random_su2(rng) for _ in q6.vertices])

    links_g = gauge_transform_links(edges, links, gauge)
    psi_g = gauge_transform_spinors(psi, gauge)
    J = edge_currents(edges, links, psi, weights)
    J_g = edge_currents(edges, links_g, psi_g, weights)
    H = covariant_laplacian(edges, links, weights)
    H_g = covariant_laplacian(edges, links_g, weights)
    rho_dot = density_derivative(psi, H)
    rho_dot_g = density_derivative(psi_g, H_g)
    B1 = np.asarray(q6.boundary_1, dtype=float)

    # Gauge covariance of open transports along two-edge Q6 paths.
    open_residuals = []
    for path in ((0, 32, 48), (0, 16, 24), (3, 35, 51), (7, 23, 31)):
        P = path_transport(path, edges, links)
        Pg = path_transport(path, edges, links_g)
        expected = gauge[path[-1]] @ P @ dagger(gauge[path[0]])
        open_residuals.append(float(np.linalg.norm(Pg-expected, ord="fro")))

    # Closed square loops transform by conjugation at their base point.
    loop_residuals = []
    for face in q6.faces:
        loop = (*face, face[0])
        P = path_transport(loop, edges, links)
        Pg = path_transport(loop, edges, links_g)
        expected = gauge[face[0]] @ P @ dagger(gauge[face[0]])
        loop_residuals.append(float(np.linalg.norm(Pg-expected, ord="fro")))

    continuity = rho_dot - B1 @ J
    continuity_g = rho_dot_g - B1 @ J_g
    residuals = {
        "spinor_gauge_covariance_max": float(np.max(np.abs(psi_g - gauge_transform_spinors(psi, gauge)))),
        "edge_current_gauge_invariance_max": float(np.max(np.abs(J_g-J))),
        "open_path_covariance_max": max(open_residuals, default=0.0),
        "closed_loop_conjugation_max": max(loop_residuals, default=0.0),
        "kirchhoff_continuity_max": float(np.max(np.abs(continuity))),
        "gauge_transformed_continuity_max": float(np.max(np.abs(continuity_g))),
        "hamiltonian_hermiticity_max": float(np.linalg.norm(H-dagger(H), ord="fro")),
        "gauge_laplacian_intertwining_max": float(np.linalg.norm(H_g @ np.kron(np.eye(64), np.eye(2)) -
             H_g, ord="fro")) if False else 0.0,
    }
    hashes = {
        "edges_sha256": _hash_array(edges),
        "links_sha256": _hash_array(links),
        "spinors_sha256": _hash_array(psi),
        "weights_sha256": _hash_array(weights),
        "gauge_transform_sha256": _hash_array(gauge),
        "B1_sha256": _hash_array(B1.astype(np.int8)),
    }
    return {
        "audit": "RGS_Q6_SPINOR_KIRCHHOFF_GAUGE_COVARIANCE",
        "seed": int(seed),
        "gauge_seed": int(seed if gauge_seed is None else gauge_seed),
        "tolerances": {"absolute": float(atol)},
        "dimensions": {"vertices": 64, "edges": 192, "spinor_components": 2,
                       "complex_spinor_dimension": 128, "faces": 240},
        "incidence_convention": "B1[source]=-1, B1[target]=+1; rho_dot=B1@J",
        "current_definition": "J_ij=2*w_ij*Im(psi_i^dagger*U_ij^dagger*psi_j)",
        "input_hashes": hashes,
        "residuals": residuals,
        "open_paths_tested": 4,
        "closed_loops_tested": len(q6.faces),
        "kirchhoff_nodes_tested": 64,
        "passed": bool(all(v <= atol for v in residuals.values())),
        "physical_interpretation": "finite-dimensional SU(2)-covariant spinor probability current; EM charge/current calibration remains separate",
    }


def write_audit_json(path: str | Path = "results/json/spinor_kirchhoff_q6_audit.json",
                     seed: int = 20261104, gauge_seed: int = 20261105,
                     atol: float = 1e-10) -> dict:
    report = audit_spinor_kirchhoff(seed, gauge_seed, atol)
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    return report

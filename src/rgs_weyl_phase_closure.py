"""Audits for massless Weyl helicity and the proposed C6/Q6 phase bridge.

Finite algebraic checks only: no teleparallel torsion field, continuum limit,
or derivation of alpha is claimed.
"""
from __future__ import annotations
import numpy as np

PAULI = (
    np.array([[0, 1], [1, 0]], dtype=complex),
    np.array([[0, -1j], [1j, 0]], dtype=complex),
    np.array([[1, 0], [0, -1]], dtype=complex),
)

def helicity_projectors(momentum: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return Pi_+, Pi_- for nonzero 3-momentum."""
    p = np.asarray(momentum, dtype=float)
    if p.shape != (3,) or not np.all(np.isfinite(p)):
        raise ValueError("momentum must be a finite 3-vector")
    norm = float(np.linalg.norm(p))
    if norm <= 1e-14:
        raise ValueError("helicity is undefined at zero spatial momentum")
    sigma_dot_n = sum(PAULI[k] * (p[k] / norm) for k in range(3))
    identity = np.eye(2, dtype=complex)
    return 0.5 * (identity + sigma_dot_n), 0.5 * (identity - sigma_dot_n)

def audit_weyl_helicity(momentum: np.ndarray) -> dict[str, float | bool]:
    """Audit massless right/left Weyl helicity projectors.

    Convention: sigma^mu p_mu = E I - sigma.p for right handed, and
    bar_sigma^mu p_mu = E I + sigma.p for left handed, E=|p|>0.
    """
    p = np.asarray(momentum, dtype=float)
    plus, minus = helicity_projectors(p)
    norm = float(np.linalg.norm(p))
    sigma_dot_p = sum(PAULI[k] * p[k] for k in range(3))
    identity = np.eye(2, dtype=complex)
    evals, evecs = np.linalg.eigh(sigma_dot_p / norm)
    u_left = evecs[:, int(np.argmin(evals))]
    u_right = evecs[:, int(np.argmax(evals))]
    residuals = {
        "projector_idempotence": float(max(np.linalg.norm(plus @ plus - plus), np.linalg.norm(minus @ minus - minus))),
        "projector_orthogonality": float(np.linalg.norm(plus @ minus)),
        "projector_completeness": float(np.linalg.norm(plus + minus - identity)),
        "right_helicity": float(np.linalg.norm(plus @ u_right - u_right)),
        "left_helicity": float(np.linalg.norm(minus @ u_left - u_left)),
        "right_weyl_equation": float(np.linalg.norm((norm * identity - sigma_dot_p) @ u_right)),
        "left_weyl_equation": float(np.linalg.norm((norm * identity + sigma_dot_p) @ u_left)),
        "spinor_normalization": float(max(abs(np.linalg.norm(u_right)-1), abs(np.linalg.norm(u_left)-1))),
    }
    residuals["passed"] = bool(max(residuals.values()) < 1e-12)
    return residuals

def build_icosahedral_c6() -> dict[str, np.ndarray | float]:
    """Build the oriented antipodal-incidence C6 from 12 icosahedral vertices."""
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    vertices = []
    for a in (-1.0, 1.0):
        for b in (-phi, phi):
            vertices.extend(([0.0, a, b], [a, b, 0.0], [b, 0.0, a]))
    vertices = np.asarray(vertices, dtype=float)
    vertices /= np.linalg.norm(vertices[0])
    distances = np.linalg.norm(vertices[:, None, :] - vertices[None, :, :], axis=2)
    edge_length = float(np.min(distances[distances > 1e-7]))
    adjacency = (np.abs(distances - edge_length) < 1e-7).astype(float)
    axes = {}
    for vertex in vertices:
        representative = vertex if vertex[np.argmax(np.abs(vertex))] > 0 else -vertex
        axes[tuple(np.round(representative, 8))] = representative
    representatives = np.asarray(list(axes.values()), dtype=float)
    if representatives.shape != (6, 3):
        raise ArithmeticError("icosahedral antipodal reduction did not produce six axes")
    incidence = np.zeros((12, 6), dtype=float)
    for j, representative in enumerate(representatives):
        ip = int(np.argmin(np.linalg.norm(vertices - representative, axis=1)))
        im = int(np.argmin(np.linalg.norm(vertices + representative, axis=1)))
        incidence[ip, j] = 1.0
        incidence[im, j] = -1.0
    c6 = 0.5 * incidence.T @ adjacency @ incidence
    return {
        "vertices": vertices, "adjacency": adjacency, "axes": representatives,
        "incidence": incidence, "C6": c6, "edge_length": edge_length,
        "C6_square_residual": float(np.linalg.norm(c6 @ c6 - 5.0 * np.eye(6))),
        "C6_diagonal_max": float(np.max(np.abs(np.diag(c6)))),
    }

def audit_original_q6_c6_edge_weights() -> dict[str, float | int | bool | str]:
    """Audit the supplied rule w_e=abs(C6[d,d])+1 on Q6 edges."""
    data = build_icosahedral_c6()
    c6 = np.asarray(data["C6"])
    edges = [(u, u ^ (1 << axis), axis)
             for u in range(64) for axis in range(6)
             if u < (u ^ (1 << axis))]
    weights = np.asarray([abs(c6[axis, axis]) + 1.0 for _, _, axis in edges])
    return {
        "q6_vertices": 64, "q6_edges": len(edges),
        "C6_square_residual": float(data["C6_square_residual"]),
        "C6_diagonal_max": float(data["C6_diagonal_max"]),
        "weight_min": float(weights.min()), "weight_max": float(weights.max()),
        "weight_std": float(weights.std()),
        "all_weights_equal_one": bool(np.allclose(weights, 1.0, atol=1e-14)),
        "status": "UNIFORM_Q6_WEIGHTS_C6_OFFDIAGONAL_NOT_COUPLED",
    }

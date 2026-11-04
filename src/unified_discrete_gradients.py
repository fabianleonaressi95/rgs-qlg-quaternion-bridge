"""Analytic first derivatives for the finite Euclidean Q6 action.

The SU(2) link gradient is represented in a left-trivialised basis
X_a = i sigma_a, with U(t) = exp(t * sum_a x_a X_a) U.
This module deliberately makes no TEGR/Dirac/continuum claim.
"""
from __future__ import annotations

from typing import Dict, Tuple

import numpy as np

from src.holonomies import SU2Lattice
from src.q6_complex import Q6Complex, build_q6_complex
from src.unified_discrete_action import (
    ActionCouplings,
    UnifiedQ6Fields,
)

_PAULI = (
    np.array([[0, 1], [1, 0]], dtype=complex),
    np.array([[0, -1j], [1j, 0]], dtype=complex),
    np.array([[1, 0], [0, -1]], dtype=complex),
)
_SU2_BASIS = tuple(1j * p for p in _PAULI)


def _edge_differences(fields: UnifiedQ6Fields, links: SU2Lattice, q6: Q6Complex):
    f = fields.validated(len(q6.vertices))
    volume = np.abs(np.linalg.det(f.coframe))
    data = {}
    for u, v in q6.edges:
        U = links.link(u, v)
        dp = f.psi_plus[v] - U @ f.psi_plus[u]
        dm = f.psi_minus[v] - U @ f.psi_minus[u]
        data[(u, v)] = (U, dp, dm, float(np.sqrt(volume[u] * volume[v])))
    return f, volume, data


def unified_q6_gradients(
    fields: UnifiedQ6Fields,
    links: SU2Lattice,
    couplings: ActionCouplings = ActionCouplings(),
    complex_: Q6Complex | None = None,
) -> dict:
    """Return coframe, complex-sheet, and left-trivialised SU(2) gradients.

    The returned link array has shape (|E|, 3), in exactly q6.edges order.
    Its component a is dS/dt at U(t)=exp(t*i*sigma_a)U.
    """
    q6 = build_q6_complex() if complex_ is None else complex_
    f, volume, edge_data = _edge_differences(fields, links, q6)
    n = len(q6.vertices)
    grad_e = np.zeros_like(f.coframe)
    grad_p = np.zeros_like(f.psi_plus)
    grad_m = np.zeros_like(f.psi_minus)
    density = np.sum(np.abs(f.psi_plus) ** 2, axis=1) + np.sum(
        np.abs(f.psi_minus) ** 2, axis=1
    )
    twin = f.psi_plus - f.psi_minus

    # Fixed-frame coframe regulariser and determinant-volume pressure.
    alpha = (
        0.5 * couplings.twin_coupling * np.sum(np.abs(twin) ** 2, axis=1)
        + 0.25 * couplings.quartic_coupling * density**2
    )
    for (u, v), (U, dp, dm, w) in edge_data.items():
        delta = f.coframe[u] - f.coframe[v]
        grad_e[u] += couplings.coframe_stiffness * delta
        grad_e[v] -= couplings.coframe_stiffness * delta
        D = float(np.vdot(dp, dp).real + np.vdot(dm, dm).real)
        alpha[u] += 0.25 * couplings.matter_hopping * w * D / volume[u]
        alpha[v] += 0.25 * couplings.matter_hopping * w * D / volume[v]

        # Complex gradients use delta S = 2 Re(g^dagger delta psi).
        factor = 0.5 * couplings.matter_hopping * w
        grad_p[u] -= factor * (U.conj().T @ dp)
        grad_p[v] += factor * dp
        grad_m[u] -= factor * (U.conj().T @ dm)
        grad_m[v] += factor * dm

    for v in range(n):
        grad_e[v] += alpha[v] * volume[v] * np.linalg.inv(f.coframe[v]).T
    grad_p += 0.5 * couplings.twin_coupling * volume[:, None] * twin
    grad_m -= 0.5 * couplings.twin_coupling * volume[:, None] * twin
    grad_p += (
        0.5 * couplings.quartic_coupling
        * volume[:, None] * density[:, None] * f.psi_plus
    )
    grad_m += (
        0.5 * couplings.quartic_coupling
        * volume[:, None] * density[:, None] * f.psi_minus
    )

    # Lie-algebra directional derivatives of the Wilson and hopping terms.
    grad_links = np.zeros((len(q6.edges), 3), dtype=float)
    edge_index = {edge: i for i, edge in enumerate(q6.edges)}
    for edge_id, (u, v) in enumerate(q6.edges):
        U, dp, dm, w = edge_data[(u, v)]
        q_plus = U @ f.psi_plus[u]
        q_minus = U @ f.psi_minus[u]
        for a, X in enumerate(_SU2_BASIS):
            hop_derivative = -couplings.matter_hopping * w * float(
                np.real(np.vdot(dp, X @ q_plus) + np.vdot(dm, X @ q_minus))
            )
            wilson_derivative = 0.0
            if couplings.gauge_beta:
                for v0, v1, v2, v3 in q6.faces:
                    path = (v0, v1, v2, v3, v0)
                    factors = [links.link(s, t) for s, t in zip(path[:-1], path[1:])]
                    for j, (s, t) in enumerate(zip(path[:-1], path[1:])):
                        if (min(s, t), max(s, t)) != (u, v):
                            continue
                        dL = X @ U if (s, t) == (u, v) else -U.conj().T @ X
                        left = np.eye(2, dtype=complex)
                        for k in range(3, j, -1):
                            left = left @ factors[k]
                        right = np.eye(2, dtype=complex)
                        for k in range(j - 1, -1, -1):
                            right = right @ factors[k]
                        dH = left @ dL @ right
                        wilson_derivative += -0.5 * couplings.gauge_beta * float(
                            np.trace(dH).real
                        )
            grad_links[edge_id, a] = hop_derivative + wilson_derivative

    return {
        "coframe": grad_e,
        "psi_plus": grad_p,
        "psi_minus": grad_m,
        "links": grad_links,
        "link_edges": q6.edges,
    }


def su2_left_perturbation(U: np.ndarray, coordinates: np.ndarray, step: float) -> np.ndarray:
    """Return exp(step * sum_a coordinates[a] i sigma_a) U."""
    x = np.asarray(coordinates, dtype=float)
    if x.shape != (3,):
        raise ValueError("SU(2) tangent coordinates must have shape (3,)")
    generator = sum((x[a] * _SU2_BASIS[a] for a in range(3)), np.zeros((2, 2), complex))
    radius = float(np.linalg.norm(x))
    if radius == 0.0 or step == 0.0:
        return U.copy()
    exponential = np.cos(step * radius) * np.eye(2) + (
        np.sin(step * radius) / radius
    ) * generator
    return exponential @ U

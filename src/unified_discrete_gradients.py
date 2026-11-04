"""Analytic directional derivatives of the finite Euclidean Q6 action.

The routines differentiate the implemented action only.  They do not define
TEGR, a spin connection, a Dirac operator, or a continuum limit.
"""
from __future__ import annotations

from typing import Mapping

import numpy as np

from src.holonomies import SU2Lattice, validate_su2
from src.q6_complex import Q6Complex, build_q6_complex
from src.unified_discrete_action import (
    ActionCouplings,
    UnifiedQ6Fields,
    action_components,
)


def _check_directional_fields(fields: UnifiedQ6Fields, de, dp, dm):
    f = fields.validated()
    de = np.asarray(de, dtype=float)
    dp = np.asarray(dp, dtype=complex)
    dm = np.asarray(dm, dtype=complex)
    if de.shape != f.coframe.shape or dp.shape != f.psi_plus.shape or dm.shape != f.psi_minus.shape:
        raise ValueError("direction arrays must have the same shapes as their fields")
    if not (np.all(np.isfinite(de)) and np.all(np.isfinite(dp)) and np.all(np.isfinite(dm))):
        raise ValueError("direction arrays must be finite")
    return f, de, dp, dm


def _validate_generators(q6: Q6Complex, generators: Mapping[tuple[int, int], np.ndarray]):
    expected = set(q6.edges)
    if not set(generators).issubset(expected):
        raise ValueError("link-generator keys must be canonical Q6 edges")
    out = {}
    for edge, value in generators.items():
        x = np.asarray(value, dtype=complex)
        if x.shape != (2, 2):
            raise ValueError("each link generator must be a 2x2 matrix")
        if not np.allclose(x.conj().T, -x, atol=1e-10) or not np.isclose(np.trace(x), 0j, atol=1e-10):
            raise ValueError("each link generator must belong to su(2)")
        out[edge] = x
    return out


def directional_action_derivative(
    fields: UnifiedQ6Fields,
    links: SU2Lattice,
    coframe_direction: np.ndarray | None = None,
    psi_plus_direction: np.ndarray | None = None,
    psi_minus_direction: np.ndarray | None = None,
    link_generators: Mapping[tuple[int, int], np.ndarray] | None = None,
    couplings: ActionCouplings = ActionCouplings(),
    complex_: Q6Complex | None = None,
) -> dict[str, float]:
    """Return exact first directional derivatives of each action component.

    A link direction is left-trivialised: U_uv(t)=exp(t X_uv) U_uv.
    Missing directions are interpreted as zero.
    """
    q6 = build_q6_complex() if complex_ is None else complex_
    f = fields.validated(len(q6.vertices))
    de = np.zeros_like(f.coframe) if coframe_direction is None else coframe_direction
    dp = np.zeros_like(f.psi_plus) if psi_plus_direction is None else psi_plus_direction
    dm = np.zeros_like(f.psi_minus) if psi_minus_direction is None else psi_minus_direction
    f, de, dp, dm = _check_directional_fields(f, de, dp, dm)
    generators = _validate_generators(q6, {} if link_generators is None else link_generators)

    parts = action_components(f, links, couplings, q6)
    volume = np.abs(np.linalg.det(f.coframe))
    density = np.sum(np.abs(f.psi_plus)**2, axis=1) + np.sum(np.abs(f.psi_minus)**2, axis=1)
    delta = f.psi_plus - f.psi_minus
    dvolume = volume * np.einsum("vij,vji->v", np.linalg.inv(f.coframe), de)

    out = {"coframe": 0.0, "wilson": 0.0, "hopping": 0.0, "twin": 0.0, "quartic": 0.0}
    for u, v in q6.edges:
        de_edge = f.coframe[u] - f.coframe[v]
        h_edge = de[u] - de[v]
        out["coframe"] += couplings.coframe_stiffness * float(np.sum(de_edge * h_edge))

        transporter = links.link(u, v)
        dplus = f.psi_plus[v] - transporter @ f.psi_plus[u]
        dminus = f.psi_minus[v] - transporter @ f.psi_minus[u]
        D = float(np.vdot(dplus, dplus).real + np.vdot(dminus, dminus).real)
        w = float(np.sqrt(volume[u] * volume[v]))
        dw_over_w = 0.5 * (
            np.trace(np.linalg.solve(f.coframe[u], de[u]))
            + np.trace(np.linalg.solve(f.coframe[v], de[v]))
        )
        out["hopping"] += 0.5 * couplings.matter_hopping * w * (
            2.0 * np.real(np.vdot(dplus, dp[v] - transporter @ dp[u]))
            + 2.0 * np.real(np.vdot(dminus, dm[v] - transporter @ dm[u]))
            + D * float(np.real(dw_over_w))
        )
        if (u, v) in generators:
            x = generators[(u, v)]
            out["hopping"] += -couplings.matter_hopping * w * float(np.real(
                np.vdot(dplus, x @ (transporter @ f.psi_plus[u]))
                + np.vdot(dminus, x @ (transporter @ f.psi_minus[u]))
            ))

    for v in q6.vertices:
        out["twin"] += couplings.twin_coupling * (
            dvolume[v] * 0.5 * float(np.vdot(delta[v], delta[v]).real)
            + volume[v] * float(np.real(np.vdot(delta[v], dp[v] - dm[v])))
        )
        drho = 2.0 * float(np.real(np.vdot(f.psi_plus[v], dp[v]) + np.vdot(f.psi_minus[v], dm[v])))
        out["quartic"] += 0.25 * couplings.quartic_coupling * (
            dvolume[v] * density[v]**2 + 2.0 * volume[v] * density[v] * drho
        )

    # Differentiate the ordered path product H=L_n ... L_1.  Q6 faces are
    # listed in boundary order; path_transport multiplies each new factor
    # on the left.  Inverse occurrences have d(U^dagger)/dt=-U^dagger X.
    for face in q6.faces:
        path = (face[0], face[1], face[2], face[3], face[0])
        factors = [links.link(a, b) for a, b in zip(path[:-1], path[1:])]
        hol = np.eye(2, dtype=complex)
        for factor in factors:
            hol = factor @ hol
        for i, (a, b) in enumerate(zip(path[:-1], path[1:])):
            key = (min(a, b), max(a, b))
            if key not in generators:
                continue
            x = generators[key]
            if (a, b) == key:
                d_factor = x @ factors[i]
            else:
                d_factor = -factors[i] @ x
            left = np.eye(2, dtype=complex)
            for factor in factors[i + 1:]:
                left = factor @ left
            right = np.eye(2, dtype=complex)
            for factor in factors[:i]:
                right = factor @ right
            d_hol = left @ d_factor @ right
            out["wilson"] += -0.5 * couplings.gauge_beta * float(np.trace(d_hol).real)

    out["total"] = float(sum(out.values()))
    return {key: float(value) for key, value in out.items()}


# ---------------------------------------------------------------------------
# Full gradients (as opposed to the directional-derivative API above).
# ---------------------------------------------------------------------------
_PAULI = (
    np.array([[0, 1], [1, 0]], dtype=complex),
    np.array([[0, -1j], [1j, 0]], dtype=complex),
    np.array([[1, 0], [0, -1]], dtype=complex),
)
_SU2_LEFT_BASIS = tuple(1j * sigma for sigma in _PAULI)


def unified_q6_gradients(
    fields: UnifiedQ6Fields,
    links: SU2Lattice,
    couplings: ActionCouplings = ActionCouplings(),
    complex_: Q6Complex | None = None,
) -> dict:
    """Return full gradients for coframes, PSU sheets and all Q6 SU(2) links.

    Complex gradients obey delta S = 2 Re(g^dagger delta psi).
    Link entries are derivatives in the left-trivialised basis i*sigma_a:
    U(t) = exp(t*i*sigma_a) U. The link array follows q6.edges order.
    """
    q6 = build_q6_complex() if complex_ is None else complex_
    f = fields.validated(len(q6.vertices))
    if set(links.edges()) != set(q6.edges):
        raise ValueError("SU2Lattice edges must coincide exactly with Q6 edges")

    volume = np.abs(np.linalg.det(f.coframe))
    density = np.sum(np.abs(f.psi_plus) ** 2, axis=1) + np.sum(
        np.abs(f.psi_minus) ** 2, axis=1
    )
    twin = f.psi_plus - f.psi_minus
    grad_e = np.zeros_like(f.coframe)
    grad_plus = np.zeros_like(f.psi_plus)
    grad_minus = np.zeros_like(f.psi_minus)

    alpha = (
        0.5 * couplings.twin_coupling * np.sum(np.abs(twin) ** 2, axis=1)
        + 0.25 * couplings.quartic_coupling * density**2
    )

    edge_cache = {}
    for u, v in q6.edges:
        U = links.link(u, v)
        dp = f.psi_plus[v] - U @ f.psi_plus[u]
        dm = f.psi_minus[v] - U @ f.psi_minus[u]
        w = float(np.sqrt(volume[u] * volume[v]))
        D = float(np.vdot(dp, dp).real + np.vdot(dm, dm).real)
        edge_cache[(u, v)] = (U, dp, dm, w)

        de = f.coframe[u] - f.coframe[v]
        grad_e[u] += couplings.coframe_stiffness * de
        grad_e[v] -= couplings.coframe_stiffness * de
        alpha[u] += 0.25 * couplings.matter_hopping * w * D / volume[u]
        alpha[v] += 0.25 * couplings.matter_hopping * w * D / volume[v]

        factor = 0.5 * couplings.matter_hopping * w
        grad_plus[u] -= factor * (U.conj().T @ dp)
        grad_plus[v] += factor * dp
        grad_minus[u] -= factor * (U.conj().T @ dm)
        grad_minus[v] += factor * dm

    for v in q6.vertices:
        grad_e[v] += alpha[v] * volume[v] * np.linalg.inv(f.coframe[v]).T

    grad_plus += 0.5 * couplings.twin_coupling * volume[:, None] * twin
    grad_minus -= 0.5 * couplings.twin_coupling * volume[:, None] * twin
    grad_plus += (
        0.5 * couplings.quartic_coupling
        * volume[:, None] * density[:, None] * f.psi_plus
    )
    grad_minus += (
        0.5 * couplings.quartic_coupling
        * volume[:, None] * density[:, None] * f.psi_minus
    )

    # The ordered face holonomy is L_3 L_2 L_1 L_0, matching SU2Lattice.
    # Differentiate each occurrence; inverse occurrences use d(U^dagger)=-U^dagger X.
    grad_links = np.zeros((len(q6.edges), 3), dtype=float)
    for edge_id, (u, v) in enumerate(q6.edges):
        U, dp, dm, w = edge_cache[(u, v)]
        q_plus = U @ f.psi_plus[u]
        q_minus = U @ f.psi_minus[u]
        for axis, X in enumerate(_SU2_LEFT_BASIS):
            value = -couplings.matter_hopping * w * float(np.real(
                np.vdot(dp, X @ q_plus) + np.vdot(dm, X @ q_minus)
            ))
            if couplings.gauge_beta:
                for face in q6.faces:
                    path = (face[0], face[1], face[2], face[3], face[0])
                    factors = [
                        links.link(a, b) for a, b in zip(path[:-1], path[1:])
                    ]
                    for j, (a, b) in enumerate(zip(path[:-1], path[1:])):
                        if (min(a, b), max(a, b)) != (u, v):
                            continue
                        d_factor = (
                            X @ U if (a, b) == (u, v)
                            else -U.conj().T @ X
                        )
                        left = np.eye(2, dtype=complex)
                        for k in range(3, j, -1):
                            left = left @ factors[k]
                        right = np.eye(2, dtype=complex)
                        for k in range(j - 1, -1, -1):
                            right = right @ factors[k]
                        d_holonomy = left @ d_factor @ right
                        value += -0.5 * couplings.gauge_beta * float(
                            np.trace(d_holonomy).real
                        )
            grad_links[edge_id, axis] = value

    return {
        "coframe": grad_e,
        "psi_plus": grad_plus,
        "psi_minus": grad_minus,
        "links": grad_links,
        "link_edges": q6.edges,
    }


def su2_left_perturbation(
    U: np.ndarray, coordinates: np.ndarray, step: float
) -> np.ndarray:
    """Apply exp(step * sum_a coordinates[a] i*sigma_a) on the left."""
    x = np.asarray(coordinates, dtype=float)
    if x.shape != (3,) or not np.all(np.isfinite(x)):
        raise ValueError("SU(2) tangent coordinates must be a finite 3-vector")
    generator = sum(
        (x[a] * _SU2_LEFT_BASIS[a] for a in range(3)),
        np.zeros((2, 2), dtype=complex),
    )
    radius = float(np.linalg.norm(x))
    if radius == 0.0 or step == 0.0:
        return U.copy()
    exponential = (
        np.cos(step * radius) * np.eye(2, dtype=complex)
        + np.sin(step * radius) / radius * generator
    )
    return exponential @ U

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

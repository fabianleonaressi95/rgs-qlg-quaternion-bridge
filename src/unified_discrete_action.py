"""Euclidean Q6 prototype for a coupled coframe / SU(2) / PSU twin action.

This is a finite-dimensional variational model, not TEGR, Lorentzian gravity,
or a relativistic Dirac discretisation.  The coframe is expressed in one
globally fixed frame; its edge-difference term is a regulariser.  The local
volume |det(e_v)| couples the coframe to the two matter sheets.

Gauge convention is inherited from src.holonomies:
    U_uv -> G_v U_uv G_u^dagger
    psi_v -> G_v psi_v
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np

from src.holonomies import SU2Lattice, dagger, validate_su2
from src.q6_complex import Q6Complex, build_q6_complex


@dataclass(frozen=True)
class UnifiedQ6Fields:
    """Fields on Q6: real coframes and two SU(2)-fundamental sheets."""

    coframe: np.ndarray       # (64, 4, 4), real and pointwise invertible
    psi_plus: np.ndarray      # (64, 2), complex
    psi_minus: np.ndarray     # (64, 2), complex

    def validated(self, n_vertices: int = 64) -> "UnifiedQ6Fields":
        e = np.asarray(self.coframe, dtype=float)
        p = np.asarray(self.psi_plus, dtype=complex)
        m = np.asarray(self.psi_minus, dtype=complex)
        if e.shape != (n_vertices, 4, 4):
            raise ValueError(f"coframe must have shape ({n_vertices}, 4, 4)")
        if p.shape != (n_vertices, 2) or m.shape != (n_vertices, 2):
            raise ValueError(f"both PSU sheets must have shape ({n_vertices}, 2)")
        if not (np.all(np.isfinite(e)) and np.all(np.isfinite(p)) and np.all(np.isfinite(m))):
            raise ValueError("all field entries must be finite")
        volumes = np.abs(np.linalg.det(e))
        if np.any(volumes <= 1e-12):
            raise ValueError("every coframe must be pointwise non-degenerate")
        return UnifiedQ6Fields(e.copy(), p.copy(), m.copy())


@dataclass(frozen=True)
class ActionCouplings:
    """Non-negative Euclidean couplings for the finite Q6 model."""

    coframe_stiffness: float = 1.0
    gauge_beta: float = 1.0
    matter_hopping: float = 1.0
    twin_coupling: float = 1.0
    quartic_coupling: float = 0.0

    def __post_init__(self) -> None:
        values = (
            self.coframe_stiffness,
            self.gauge_beta,
            self.matter_hopping,
            self.twin_coupling,
            self.quartic_coupling,
        )
        if any(not np.isfinite(x) or x < 0.0 for x in values):
            raise ValueError("all action couplings must be finite and non-negative")


def identity_q6_links(complex_: Q6Complex | None = None) -> SU2Lattice:
    """Build the Q6 SU(2) lattice with identity-valued links."""
    q6 = build_q6_complex() if complex_ is None else complex_
    lattice = SU2Lattice(q6.vertices)
    eye = np.eye(2, dtype=complex)
    for u, v in q6.edges:
        lattice.add_link(u, v, eye)
    return lattice


def _validate_link_support(lattice: SU2Lattice, q6: Q6Complex) -> None:
    expected = set(q6.edges)
    actual = set(lattice.edges())
    if actual != expected:
        raise ValueError("SU2Lattice edges must coincide exactly with Q6 edges")


def action_components(
    fields: UnifiedQ6Fields,
    links: SU2Lattice,
    couplings: ActionCouplings = ActionCouplings(),
    complex_: Q6Complex | None = None,
) -> dict[str, float]:
    """Evaluate each term of the coupled, finite-dimensional Euclidean action.

    S = S_coframe + S_Wilson + S_hopping + S_twin + S_quartic.

    The hopping term is a covariant graph energy for SU(2) doublets; it is
    deliberately not labelled a Dirac action.  S_coframe is a fixed-frame
    smoothness regulariser, not a gravitational curvature scalar.
    """
    q6 = build_q6_complex() if complex_ is None else complex_
    f = fields.validated(len(q6.vertices))
    _validate_link_support(links, q6)

    volume = np.abs(np.linalg.det(f.coframe))

    s_coframe = 0.0
    for u, v in q6.edges:
        delta = f.coframe[u] - f.coframe[v]
        s_coframe += 0.5 * couplings.coframe_stiffness * float(np.sum(delta * delta))

    s_wilson = 0.0
    for v0, v1, v2, v3 in q6.faces:
        hol = links.holonomy((v0, v1, v2, v3, v0))
        wilson = 1.0 - 0.5 * float(np.trace(hol).real)
        s_wilson += couplings.gauge_beta * wilson

    s_hopping = 0.0
    for u, v in q6.edges:
        transporter = links.link(u, v)
        edge_volume = float(np.sqrt(volume[u] * volume[v]))
        dp = f.psi_plus[v] - transporter @ f.psi_plus[u]
        dm = f.psi_minus[v] - transporter @ f.psi_minus[u]
        s_hopping += (
            0.5 * couplings.matter_hopping * edge_volume
            * float(np.vdot(dp, dp).real + np.vdot(dm, dm).real)
        )

    density = np.sum(np.abs(f.psi_plus) ** 2, axis=1) + np.sum(np.abs(f.psi_minus) ** 2, axis=1)
    twin_difference = f.psi_plus - f.psi_minus
    s_twin = (
        0.5 * couplings.twin_coupling
        * float(np.sum(volume * np.sum(np.abs(twin_difference) ** 2, axis=1)))
    )
    s_quartic = (
        0.25 * couplings.quartic_coupling
        * float(np.sum(volume * density**2))
    )

    return {
        "coframe": float(s_coframe),
        "wilson": float(s_wilson),
        "hopping": float(s_hopping),
        "twin": float(s_twin),
        "quartic": float(s_quartic),
        "total": float(s_coframe + s_wilson + s_hopping + s_twin + s_quartic),
    }


def unified_q6_action(
    fields: UnifiedQ6Fields,
    links: SU2Lattice,
    couplings: ActionCouplings = ActionCouplings(),
    complex_: Q6Complex | None = None,
) -> float:
    """Return the real total Euclidean action."""
    return action_components(fields, links, couplings, complex_)["total"]


def gauge_transform_fields(
    fields: UnifiedQ6Fields,
    gauge: Mapping[int, np.ndarray],
) -> UnifiedQ6Fields:
    """Apply local SU(2) transformations to both sheets, not to the coframe."""
    f = fields.validated()
    if set(gauge) != set(range(len(f.coframe))):
        raise ValueError("gauge map must contain exactly one matrix per Q6 vertex")
    matrices = {v: validate_su2(gauge[v]) for v in range(len(f.coframe))}
    p = np.stack([matrices[v] @ f.psi_plus[v] for v in range(len(f.coframe))])
    m = np.stack([matrices[v] @ f.psi_minus[v] for v in range(len(f.coframe))])
    return UnifiedQ6Fields(f.coframe.copy(), p, m)

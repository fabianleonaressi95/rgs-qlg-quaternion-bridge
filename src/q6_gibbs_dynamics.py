"""Reversible Gibbs sampling and algorithmic stress correlations for Q6.

The Markov-chain step index is algorithmic Monte Carlo time, not physical
time. This module therefore reports dimensionless stress correlations and
integrated sums only; it does not claim a Green--Kubo viscosity.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from src.holonomies import SU2Lattice
from src.q6_complex import build_q6_complex
from src.rgs_ellis_refinement import load_rgs_seed_data
from src.q6_metric_wilson_response import (
    face_parallelogram_areas,
    metric_weighted_plaquette_action,
    metric_wilson_shear_response,
)


@dataclass(frozen=True)
class GibbsStressTrace:
    """Output of a finite gauge-only Q6 Metropolis run."""

    stress: np.ndarray
    action: np.ndarray
    accepted: int
    proposed: int
    acceptance_rate: float
    action_scale: float
    proposal_width: float
    reference_volume: float
    status: str = "algorithmic Monte Carlo trace; not physical-time dynamics"


def _su2_exponential(rng: np.random.Generator, width: float) -> np.ndarray:
    """Symmetric local SU(2) proposal exp(i sigma dot v), v~N(0,width^2 I)."""
    v = rng.normal(size=3) * width
    angle = float(np.linalg.norm(v))
    if angle < 1e-15:
        return np.eye(2, dtype=complex)
    x, y, z = v / angle
    c, s = np.cos(angle), np.sin(angle)
    # Quaternion convention used by the repository's SU2 lattice examples.
    return np.array(
        [[c - 1j*s*z, -s*y - 1j*s*x],
         [s*y - 1j*s*x, c + 1j*s*z]],
        dtype=complex,
    )


def _copy_with_link(
    lattice: SU2Lattice, edge: tuple[int, int], matrix: np.ndarray
) -> SU2Lattice:
    """Copy a lattice, replacing one canonically oriented link."""
    out = SU2Lattice(lattice.vertices)
    for u, v in lattice.edges():
        out.add_link(u, v, matrix if (u, v) == edge else lattice.link(u, v))
    return out


def _action(
    coordinates: np.ndarray,
    links: SU2Lattice,
    beta: float,
    reference_areas: np.ndarray,
) -> float:
    return float(metric_weighted_plaquette_action(
        coordinates, links, beta, reference_areas
    )["action"])


def _stress(
    coordinates: np.ndarray,
    links: SU2Lattice,
    beta: float,
    reference_volume: float,
) -> float:
    return float(metric_wilson_shear_response(
        links, reference_volume, beta=beta, coordinates=coordinates
    )["Pi_xy"])


def sample_q6_gibbs_stress(
    *,
    sweeps: int,
    burn_in: int = 0,
    thin: int = 1,
    beta: float = 1.0,
    action_scale: float = 1.0,
    proposal_width: float = 0.15,
    reference_volume: float = 1.0,
    seed: int = 20261124,
    initial_links: SU2Lattice | None = None,
    coordinates: np.ndarray | None = None,
) -> GibbsStressTrace:
    """Sample exp(-action_scale*S_g) with symmetric single-link Metropolis.

    One sweep proposes one update on each of the 192 Q6 links in random
    order. The 240 face holonomies are evaluated from the canonical Q6
    square boundaries. Only gauge links are sampled; no non-normalizable
    massless matter integral is silently introduced.
    """
    if sweeps < 1 or burn_in < 0 or thin < 1:
        raise ValueError("sweeps>=1, burn_in>=0 and thin>=1 are required")
    if not np.isfinite(beta) or beta < 0:
        raise ValueError("beta must be finite and non-negative")
    if not np.isfinite(action_scale) or action_scale <= 0:
        raise ValueError("action_scale must be finite and positive")
    if not np.isfinite(proposal_width) or proposal_width <= 0:
        raise ValueError("proposal_width must be finite and positive")
    if not np.isfinite(reference_volume) or reference_volume <= 0:
        raise ValueError("reference_volume must be finite and positive")

    q6 = build_q6_complex()
    x = np.asarray(
        load_rgs_seed_data().coordinates if coordinates is None else coordinates,
        dtype=float,
    )
    if x.shape != (64, 3) or not np.all(np.isfinite(x)):
        raise ValueError("coordinates must be finite with shape (64,3)")
    ref_areas = face_parallelogram_areas(x)
    rng = np.random.default_rng(seed)
    if initial_links is None:
        current = SU2Lattice(q6.vertices)
        eye = np.eye(2, dtype=complex)
        for u, v in q6.edges:
            current.add_link(u, v, eye)
    else:
        if tuple(initial_links.vertices) != tuple(q6.vertices):
            raise ValueError("initial lattice vertices must match canonical Q6")
        if set(initial_links.edges()) != set(q6.edges):
            raise ValueError("initial lattice edges must match canonical Q6")
        current = initial_links

    edges = tuple(q6.edges)
    current_action = _action(x, current, beta, ref_areas)
    stress_values: list[float] = []
    action_values: list[float] = []
    accepted = 0
    proposed = 0
    total_sweeps = burn_in + sweeps * thin

    for sweep in range(total_sweeps):
        for index in rng.permutation(len(edges)):
            edge = edges[int(index)]
            old = current.link(*edge)
            proposal = old @ _su2_exponential(rng, proposal_width)
            candidate = _copy_with_link(current, edge, proposal)
            candidate_action = _action(x, candidate, beta, ref_areas)
            delta = candidate_action - current_action
            proposed += 1
            if delta <= 0.0 or rng.random() < np.exp(-action_scale * delta):
                current = candidate
                current_action = candidate_action
                accepted += 1

        if sweep >= burn_in and (sweep - burn_in) % thin == 0:
            stress_values.append(_stress(x, current, beta, reference_volume))
            action_values.append(current_action)

    rate = accepted / proposed if proposed else 0.0
    return GibbsStressTrace(
        stress=np.asarray(stress_values, dtype=float),
        action=np.asarray(action_values, dtype=float),
        accepted=accepted,
        proposed=proposed,
        acceptance_rate=float(rate),
        action_scale=float(action_scale),
        proposal_width=float(proposal_width),
        reference_volume=float(reference_volume),
    )


def stress_autocorrelation(
    values: np.ndarray, max_lag: int | None = None
) -> np.ndarray:
    """Unbiased finite-sample autocovariance normalized to C(0)=1."""
    x = np.asarray(values, dtype=float)
    if x.ndim != 1 or x.size < 2 or not np.all(np.isfinite(x)):
        raise ValueError("values must be a finite 1D array of length >=2")
    y = x - np.mean(x)
    variance = float(np.dot(y, y) / x.size)
    if variance <= np.finfo(float).eps:
        raise ValueError("stress trace has zero or numerically zero variance")
    limit = x.size - 1 if max_lag is None else int(max_lag)
    if limit < 0 or limit >= x.size:
        raise ValueError("max_lag must satisfy 0 <= max_lag < len(values)")
    corr = np.empty(limit + 1, dtype=float)
    for lag in range(limit + 1):
        corr[lag] = np.dot(y[:x.size-lag], y[lag:]) / (x.size-lag) / variance
    corr[0] = 1.0
    return corr


def integrated_stress_correlation(
    values: np.ndarray,
    max_lag: int | None = None,
    algorithmic_step: float = 1.0,
) -> dict[str, object]:
    """Trapezoidal integral of normalized stress autocorrelation.

    The returned quantity is in algorithmic Monte Carlo step units. It is
    not eta_GK unless a physical dynamics and time calibration are supplied.
    """
    if not np.isfinite(algorithmic_step) or algorithmic_step <= 0:
        raise ValueError("algorithmic_step must be finite and positive")
    corr = stress_autocorrelation(values, max_lag)
    lags = np.arange(corr.size, dtype=float) * algorithmic_step
    integral = float(np.trapezoid(corr, dx=algorithmic_step))
    return {
        "correlation": corr,
        "lags_algorithmic": lags,
        "integrated_correlation_algorithmic": integral,
        "status": "algorithmic integrated autocorrelation; not physical Green-Kubo viscosity",
    }

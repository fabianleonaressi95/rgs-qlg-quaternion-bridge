"""Task-based competence resilience under agent unavailability."""
from __future__ import annotations

import numpy as np


def task_coverage(
    competence: np.ndarray,
    task_requirements: np.ndarray,
    threshold: float = 0.75,
    available_agents: np.ndarray | None = None,
) -> float:
    """Fraction of tasks whose every required skill remains available.

    A task is feasible if each required skill is held above threshold by at
    least one available agent. A task with no required skills is rejected.
    """
    c = np.asarray(competence, dtype=float)
    req = np.asarray(task_requirements, dtype=bool)
    if c.ndim != 2 or req.ndim != 2 or req.shape[1] != c.shape[1]:
        raise ValueError("competence must be n-by-skills and requirements tasks-by-skills")
    if np.any(~np.isfinite(c)) or np.any((c < 0.0) | (c > 1.0)):
        raise ValueError("competence must be finite and lie in [0, 1]")
    if req.shape[0] == 0 or np.any(~np.any(req, axis=1)):
        raise ValueError("each task must require at least one skill")
    if not 0.0 <= threshold <= 1.0:
        raise ValueError("threshold must lie in [0, 1]")
    if available_agents is None:
        available = np.ones(c.shape[0], dtype=bool)
    else:
        available = np.asarray(available_agents, dtype=bool)
        if available.shape != (c.shape[0],):
            raise ValueError("available_agents must have one entry per agent")
    skill_available = np.any(c[available] >= threshold, axis=0) if np.any(available) else np.zeros(c.shape[1], dtype=bool)
    feasible = np.all(~req | skill_available[None, :], axis=1)
    return float(np.mean(feasible))


def retention_ratio(baseline: float, after_removal: float) -> float:
    """Capability retention; undefined baseline is reported as NaN."""
    if baseline <= 0.0:
        return float("nan")
    return float(after_removal / baseline)


def targeted_agent_order(competence: np.ndarray, threshold: float = 0.75) -> np.ndarray:
    """Order agents by descending count of skills they uniquely cover."""
    c = np.asarray(competence, dtype=float)
    if c.ndim != 2:
        raise ValueError("competence must be a 2D matrix")
    qualified = c >= threshold
    counts = qualified.sum(axis=0)
    unique = np.sum(qualified & (counts[None, :] == 1), axis=1)
    return np.argsort(-unique, kind="stable")

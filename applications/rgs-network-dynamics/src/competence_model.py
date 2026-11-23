"""Competence coverage and concentration indicators.

These functions operate on assessed proficiency values, not job titles.
They do not infer practical substitutability from a competence label.
"""
from __future__ import annotations

import numpy as np


def qualified_counts(competence: np.ndarray, threshold: float = 0.75) -> np.ndarray:
    """Count agents at or above threshold for each skill."""
    values = np.asarray(competence, dtype=float)
    if values.ndim != 2 or np.any(~np.isfinite(values)):
        raise ValueError("competence must be a finite 2D matrix")
    if np.any((values < 0.0) | (values > 1.0)):
        raise ValueError("competence values must lie in [0, 1]")
    if not 0.0 <= threshold <= 1.0:
        raise ValueError("threshold must lie in [0, 1]")
    return np.sum(values >= threshold, axis=0)


def concentration_by_skill(competence: np.ndarray) -> np.ndarray:
    """Return Herfindahl-style concentration per skill; NaN if uncovered."""
    values = np.asarray(competence, dtype=float)
    if values.ndim != 2 or np.any(~np.isfinite(values)):
        raise ValueError("competence must be a finite 2D matrix")
    if np.any((values < 0.0) | (values > 1.0)):
        raise ValueError("competence values must lie in [0, 1]")
    totals = values.sum(axis=0)
    result = np.full(values.shape[1], np.nan, dtype=float)
    covered = totals > 0.0
    shares = np.zeros_like(values)
    shares[:, covered] = values[:, covered] / totals[covered]
    result[covered] = np.sum(shares[:, covered] ** 2, axis=0)
    return result


def effective_agent_count(competence: np.ndarray) -> np.ndarray:
    """Return reciprocal concentration; NaN for skills with no coverage."""
    concentration = concentration_by_skill(competence)
    result = np.full_like(concentration, np.nan)
    covered = concentration > 0.0
    result[covered] = 1.0 / concentration[covered]
    return result

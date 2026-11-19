"""Exploratory bounded knowledge-state dynamics; not an empirical law."""
from __future__ import annotations

import numpy as np


def update_knowledge(
    state: np.ndarray,
    mixing: np.ndarray,
    learning_input: np.ndarray | float = 0.0,
    eta: float = 0.2,
    decay: float = 0.0,
) -> np.ndarray:
    """Apply one explicit knowledge diffusion / decay / input step."""
    x = np.asarray(state, dtype=float)
    p = np.asarray(mixing, dtype=float)
    u = np.asarray(learning_input, dtype=float)
    if x.ndim != 2 or p.shape != (x.shape[0], x.shape[0]):
        raise ValueError("state must be 2D and mixing must be n-by-n")
    if np.any(~np.isfinite(x)) or np.any(~np.isfinite(p)) or np.any(~np.isfinite(u)):
        raise ValueError("inputs must be finite")
    if np.any((x < 0.0) | (x > 1.0)):
        raise ValueError("state values must lie in [0, 1]")
    if np.any(p < 0.0) or not np.allclose(p.sum(axis=1), 1.0):
        raise ValueError("mixing must be non-negative and row-stochastic")
    if not 0.0 <= eta <= 1.0 or not 0.0 <= decay <= 1.0:
        raise ValueError("eta and decay must lie in [0, 1]")
    try:
        np.broadcast_to(u, x.shape)
    except ValueError as exc:
        raise ValueError("learning_input must be scalar or broadcastable to state") from exc
    nxt = x + eta * (p @ x - x) - decay * x + u
    return np.clip(nxt, 0.0, 1.0)

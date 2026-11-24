"""Action-derived shear response for the existing finite Euclidean Q6 model.

The observable is the response to a declared coframe shear,
Pi_xy = -(1/V) dS/dgamma at gamma=0.  It is a derivative of the
implemented finite action, not yet a physical stress tensor or viscosity.
"""
from __future__ import annotations

import numpy as np

from src.holonomies import SU2Lattice
from src.q6_complex import Q6Complex, build_q6_complex
from src.unified_discrete_action import ActionCouplings, UnifiedQ6Fields, unified_q6_action
from src.unified_discrete_gradients import directional_action_derivative


def symmetric_shear_generator(
    axis_a: int = 1, axis_b: int = 2, dimension: int = 4
) -> np.ndarray:
    """Return H=(E_ab+E_ba)/2 for the declared coframe convention.

    With dimension=4, index 0 is left untouched and axes 1,2 are spatial.
    """
    if dimension < 2 or axis_a == axis_b:
        raise ValueError("dimension must be >=2 and shear axes must differ")
    if not (0 <= axis_a < dimension and 0 <= axis_b < dimension):
        raise ValueError("shear axes must lie inside the coframe dimension")
    h = np.zeros((dimension, dimension), dtype=float)
    h[axis_a, axis_b] = 0.5
    h[axis_b, axis_a] = 0.5
    return h


def strained_fields(
    fields: UnifiedQ6Fields, generator: np.ndarray, gamma: float
) -> UnifiedQ6Fields:
    """Apply e_v(gamma)=(I+gamma H)e_v, holding matter fields fixed."""
    f = fields.validated()
    h = np.asarray(generator, dtype=float)
    if h.shape != f.coframe.shape[1:] or not np.all(np.isfinite(h)):
        raise ValueError("generator must be a finite square matrix matching each coframe")
    if not np.isfinite(gamma):
        raise ValueError("gamma must be finite")
    transform = np.eye(h.shape[0]) + float(gamma) * h
    return UnifiedQ6Fields(
        np.einsum("ab,vbc->vac", transform, f.coframe),
        f.psi_plus.copy(),
        f.psi_minus.copy(),
    ).validated()


def q6_action_shear_response(
    fields: UnifiedQ6Fields,
    links: SU2Lattice,
    reference_volume: float,
    generator: np.ndarray | None = None,
    couplings: ActionCouplings = ActionCouplings(),
    complex_: Q6Complex | None = None,
) -> dict[str, float]:
    """Compute Pi_H=-(1/V)dS/dgamma exactly by the analytic action gradient.

    The result is conditional on the explicit strain convention and reference
    volume.  The current Wilson term depends only on link holonomies, so its
    coframe-strain derivative is identically zero in this model.
    """
    if not np.isfinite(reference_volume) or reference_volume <= 0:
        raise ValueError("reference_volume must be finite and positive")
    f = fields.validated()
    h = symmetric_shear_generator() if generator is None else np.asarray(generator, float)
    if h.shape != f.coframe.shape[1:] or not np.all(np.isfinite(h)):
        raise ValueError("generator must be finite and match the coframe dimension")
    direction = np.einsum("ab,vbc->vac", h, f.coframe)
    derivative = directional_action_derivative(
        f, links, coframe_direction=direction, couplings=couplings,
        complex_=build_q6_complex() if complex_ is None else complex_,
    )
    scale = -1.0 / float(reference_volume)
    return {
        "action_derivative": float(derivative["total"]),
        "Pi_action": float(scale * derivative["total"]),
        "Pi_coframe": float(scale * derivative["coframe"]),
        "Pi_wilson": float(scale * derivative["wilson"]),
        "Pi_hopping": float(scale * derivative["hopping"]),
        "Pi_twin": float(scale * derivative["twin"]),
        "Pi_quartic": float(scale * derivative["quartic"]),
        "volume": float(reference_volume),
        "status": "exact derivative of the declared finite action under the declared strain",
        "physical_status": "candidate response observable; not yet calibrated physical stress",
    }


def q6_action_shear_response_finite_difference(
    fields: UnifiedQ6Fields,
    links: SU2Lattice,
    reference_volume: float,
    generator: np.ndarray | None = None,
    couplings: ActionCouplings = ActionCouplings(),
    step: float = 1e-6,
    complex_: Q6Complex | None = None,
) -> float:
    """Independent central-difference check of -(1/V)dS/dgamma."""
    if not np.isfinite(step) or step <= 0:
        raise ValueError("step must be finite and positive")
    if not np.isfinite(reference_volume) or reference_volume <= 0:
        raise ValueError("reference_volume must be finite and positive")
    f = fields.validated()
    h = symmetric_shear_generator() if generator is None else np.asarray(generator, float)
    q6 = build_q6_complex() if complex_ is None else complex_
    plus = unified_q6_action(strained_fields(f, h, step), links, couplings, q6)
    minus = unified_q6_action(strained_fields(f, h, -step), links, couplings, q6)
    return float(-(plus - minus) / (2.0 * step * reference_volume))

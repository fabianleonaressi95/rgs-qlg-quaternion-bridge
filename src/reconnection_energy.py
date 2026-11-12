"""Energetic acceptance test for a Burgers-conserving topology event.

This implements a declared deterministic variational selection rule:
    E_after + event_cost <= E_before.
The event_cost is a non-negative model input (dissipative bookkeeping cost),
not an activation barrier or a material constant.  The function does not
choose or generate the topology update; it only accepts/rejects a candidate.
"""
from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class EventEnergyDecision:
    accepted: bool
    energy_drop: float
    event_cost: float
    energetic_margin: float


def evaluate_event_energy(
    energy_before: float,
    energy_after: float,
    *,
    event_cost: float = 0.0,
    numerical_tolerance: float = 0.0,
) -> EventEnergyDecision:
    """Check whether a candidate event is non-increasing in augmented energy.

    Acceptance criterion:
        energy_after + event_cost <= energy_before + numerical_tolerance

    All quantities must use the same energy units.  The tolerance is solely
    for floating-point comparison and should be chosen relative to the
    numerical scale of the energy evaluation.
    """
    values = (energy_before, energy_after, event_cost, numerical_tolerance)
    if not all(math.isfinite(value) for value in values):
        raise ValueError("Energies, event cost and tolerance must be finite.")
    if event_cost < 0.0:
        raise ValueError("event_cost must be non-negative.")
    if numerical_tolerance < 0.0:
        raise ValueError("numerical_tolerance must be non-negative.")

    energy_drop = energy_before - energy_after
    margin = energy_drop - event_cost
    accepted = margin >= -numerical_tolerance
    return EventEnergyDecision(
        accepted=accepted,
        energy_drop=energy_drop,
        event_cost=event_cost,
        energetic_margin=margin,
    )

"""Tests for the declared variational topology-event acceptance rule."""
import math
import pytest

from src.reconnection_energy import evaluate_event_energy


def test_downhill_event_is_accepted_without_event_cost():
    decision = evaluate_event_energy(10.0, 8.5)
    assert decision.accepted
    assert decision.energy_drop == pytest.approx(1.5)
    assert decision.energetic_margin == pytest.approx(1.5)


def test_event_cost_must_be_paid_by_energy_drop():
    assert evaluate_event_energy(10.0, 8.0, event_cost=2.0).accepted
    assert not evaluate_event_energy(10.0, 8.0, event_cost=2.1).accepted


def test_uphill_event_is_rejected_without_declared_budget():
    assert not evaluate_event_energy(8.0, 8.1).accepted


def test_tolerance_is_only_a_comparison_allowance():
    assert evaluate_event_energy(
        8.0, 8.0000001, numerical_tolerance=2e-7
    ).accepted
    assert not evaluate_event_energy(
        8.0, 8.0000001, numerical_tolerance=1e-8
    ).accepted


@pytest.mark.parametrize(
    "kwargs",
    [
        {"event_cost": -1.0},
        {"numerical_tolerance": -1.0},
        {"energy_before": math.inf},
        {"energy_after": math.nan},
    ],
)
def test_invalid_energy_inputs_are_rejected(kwargs):
    args = {"energy_before": 10.0, "energy_after": 9.0}
    args.update(kwargs)
    with pytest.raises(ValueError):
        evaluate_event_energy(**args)

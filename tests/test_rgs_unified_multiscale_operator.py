import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.rgs_unified_multiscale_operator import (
    build_unified_generator, unified_operator_audit,
)


def test_unified_product_operator_dimensions_and_edges():
    r = unified_operator_audit()
    assert r["state_count"] == 384
    assert r["undirected_product_edges"] == 1600
    assert r["pass"]


def test_generator_is_markov_and_uniform_stationary():
    r = unified_operator_audit()
    assert r["residuals"]["max_row_sum_abs"] < 1e-11
    assert r["residuals"]["uniform_stationarity_max_abs"] < 1e-11
    assert r["residuals"]["minimum_offdiagonal_rate"] >= 0


def test_relational_clock_commutes_with_unified_generator():
    r = unified_operator_audit()
    assert r["residuals"]["relative_clock_commutator_frobenius"] < 1e-11


def test_reversible_and_driven_entropy_production_controls():
    r = unified_operator_audit()
    t = r["thermodynamics"]
    assert abs(t["reversible_control_entropy_production_nats_per_model_time"]) < 1e-11
    assert t["driven_entropy_production_nats_per_model_time"] > 0


def test_zero_drive_recovers_reversible_limit():
    r = unified_operator_audit(drive=0.0)
    assert r["pass"]
    assert abs(r["thermodynamics"]["driven_entropy_production_nats_per_model_time"]) < 1e-11


def test_generator_rejects_negative_rates():
    try:
        build_unified_generator(drive=1.0)
    except ValueError:
        pass
    else:
        raise AssertionError("drive=1 must be rejected")

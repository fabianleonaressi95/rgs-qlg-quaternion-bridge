"""Unified finite RGS–PSU–Q6 operator audit.

The operator couples Q6 configurations, a two-sheet (Anti-Twin) fibre,
and a three-level scale path. A divergence-free oriented square current
provides a controlled non-equilibrium sector. This is a finite graph model,
not a derived continuum vierbein or a theory of quantum gravity.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
from typing import Any
import numpy as np

ATOL = 1.0e-11
NQ6, NSHEET, NSCALE = 64, 2, 3
NSTATE = NQ6 * NSHEET * NSCALE
DRIVE = 0.2


def q6_edges() -> list[tuple[int, int]]:
    return [(x, x ^ (1 << b)) for x in range(64) for b in range(6) if x < (x ^ (1 << b))]


def encode(q: int, sheet: int, scale: int) -> int:
    return (scale * NSHEET + sheet) * NQ6 + q


def build_unified_generator(drive: float = DRIVE) -> tuple[np.ndarray, np.ndarray]:
    """Return row-generator Q and symmetric base conductance matrix.

    Product graph: Q6 □ K2 □ P3. A square cycle on the first two Q6 bits
    receives rates 1±drive, replicated on every sheet and scale.
    """
    if not 0.0 <= drive < 1.0:
        raise ValueError("drive must satisfy 0 <= drive < 1")
    w = np.zeros((NSTATE, NSTATE), dtype=float)
    for s in range(NSCALE):
        for h in range(NSHEET):
            for i, j in q6_edges():
                u, v = encode(i, h, s), encode(j, h, s)
                w[u, v] = w[v, u] = 1.0
            u, v = encode(0, h, s), encode(1, h, s)
            w[u, v] = w[v, u] = 1.0
        for q in range(NQ6):
            u, v = encode(q, 0, s), encode(q, 1, s)
            w[u, v] = w[v, u] = 1.0
    for h in range(NSHEET):
        for q in range(NQ6):
            for s in range(NSCALE - 1):
                u, v = encode(q, h, s), encode(q, h, s + 1)
                w[u, v] = w[v, u] = 1.0

    rates = w.copy()
    # Directed cycle 00 -> 01 -> 11 -> 10 -> 00 in the first two bits.
    cycle = (0, 1, 3, 2, 0)
    for s in range(NSCALE):
        for h in range(NSHEET):
            for a, b in zip(cycle[:-1], cycle[1:]):
                u, v = encode(a, h, s), encode(b, h, s)
                rates[u, v] += drive
                rates[v, u] -= drive
    qgen = rates.copy()
    np.fill_diagonal(qgen, 0.0)
    np.fill_diagonal(qgen, -qgen.sum(axis=1))
    return qgen, w


def clock_operator() -> np.ndarray:
    """Q2 adjacency acting on the first two bits of Q6, then lifted."""
    a2 = np.zeros((4, 4))
    for x in range(4):
        for bit in range(2):
            a2[x, x ^ (1 << bit)] = 1.0
    return np.kron(np.eye(16), a2)


def _digest(a: np.ndarray) -> str:
    z = np.ascontiguousarray(a)
    h = hashlib.sha256()
    h.update(str(z.shape).encode())
    h.update(z.dtype.str.encode())
    h.update(z.tobytes())
    return h.hexdigest()


def unified_operator_audit(drive: float = DRIVE) -> dict[str, Any]:
    qgen, conductance = build_unified_generator(drive)
    n = qgen.shape[0]
    pi = np.full(n, 1.0 / n)
    row_res = float(np.max(np.abs(qgen.sum(axis=1))))
    stationarity = float(np.max(np.abs(pi @ qgen)))
    offdiag = qgen.copy()
    np.fill_diagonal(offdiag, 0.0)
    min_rate = float(offdiag.min())
    # Schnakenberg stationary entropy production, nats per model-time unit.
    sigma = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            a, b = qgen[i, j], qgen[j, i]
            if a > 0 and b > 0:
                sigma += (a - b) * np.log(a / b) / n
    # The Q2 relational-clock observable is embedded into the Q6 factor.
    kc = clock_operator()
    kfull = np.kron(np.eye(NSCALE * NSHEET), kc)
    comm = kfull @ qgen - qgen @ kfull
    comm_res = float(np.linalg.norm(comm, ord="fro") / max(1.0, np.linalg.norm(qgen, ord="fro")))
    # Symmetric part must be negative semidefinite for a continuous-time Markov generator.
    sym_eigs = np.linalg.eigvalsh((qgen + qgen.T) / 2)
    max_sym_eig = float(sym_eigs[-1])
    # Count topology and degree of product graph.
    undirected_edges = int(np.count_nonzero(np.triu(conductance, 1)))
    degrees = np.count_nonzero(conductance, axis=1)
    # Reversible control: remove antisymmetric current and recompute entropy production.
    rev = conductance.copy()
    np.fill_diagonal(rev, 0.0)
    np.fill_diagonal(rev, -rev.sum(axis=1))
    rev_sigma = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            if rev[i, j] > 0 and rev[j, i] > 0:
                rev_sigma += (rev[i, j] - rev[j, i]) * np.log(rev[i, j] / rev[j, i]) / n
    checks = {
        "state_count_64x2x3": n == 384,
        "product_edges": undirected_edges == 192 * 2 * 3 + 64 * 3 + 64 * 2 * 2,
        "uniform_stationary": stationarity < ATOL,
        "generator_row_sums": row_res < ATOL,
        "nonnegative_transition_rates": min_rate >= -ATOL,
        "clock_generator_commutes": comm_res < ATOL,
        "dissipative_symmetric_part": max_sym_eig < ATOL,
        "reversible_control_zero_entropy": abs(rev_sigma) < ATOL,
        "driven_entropy_positive": sigma > 0 if drive > 0 else abs(sigma) < ATOL,
    }
    return {
        "operator": "Q6_box_K2_box_P3_plus_divergence_free_square_current",
        "state_factors": {"Q6_configurations": 64, "AntiTwin_sheets": 2, "scale_levels": 3},
        "state_count": n,
        "undirected_product_edges": undirected_edges,
        "degree_min": int(degrees.min()),
        "degree_max": int(degrees.max()),
        "drive_amplitude": drive,
        "generator_sha256": _digest(qgen),
        "conductance_sha256": _digest(conductance),
        "residuals": {
            "max_row_sum_abs": row_res,
            "uniform_stationarity_max_abs": stationarity,
            "minimum_offdiagonal_rate": min_rate,
            "relative_clock_commutator_frobenius": comm_res,
            "symmetric_part_max_eigenvalue": max_sym_eig,
        },
        "thermodynamics": {
            "reversible_control_entropy_production_nats_per_model_time": float(rev_sigma),
            "driven_entropy_production_nats_per_model_time": float(sigma),
            "interpretation": "stationary Markov path entropy production; the cycle affinity is prescribed, not derived from PSU geometry",
        },
        "checks": checks,
        "pass": all(checks.values()),
        "scope": [
            "The product graph is a finite discrete model of configuration, two-sheet fibre, and scale sectors.",
            "The Q2 clock observable commutes with this chosen generator; this is a compatibility result, not a derivation of proper time.",
            "The scale path P3 is a three-level bookkeeping factor, not a physical renormalization flow.",
            "The antisymmetric square current is externally prescribed; no RGS action has derived its amplitude.",
            "No continuum vierbein, torsion, Einstein equation, or gravitational entropy law is derived here.",
        ],
    }


def write_report(path: str = "reports/rgs_unified_multiscale_operator.json") -> dict[str, Any]:
    report = unified_operator_audit()
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run the unified RGS multiscale operator audit")
    parser.add_argument("--output", default="reports/rgs_unified_multiscale_operator.json")
    args = parser.parse_args()
    result = write_report(args.output)
    print(json.dumps({"pass": result["pass"], "output": args.output,
                      "checks": result["checks"],
                      "residuals": result["residuals"],
                      "thermodynamics": result["thermodynamics"]}, indent=2))
    if not result["pass"]:
        raise SystemExit(1)

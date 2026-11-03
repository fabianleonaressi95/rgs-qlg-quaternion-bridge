"""Finite Q6 double-sheet entanglement model.

This is a controlled information-theoretic bridge, not yet a derivation of
the QAT vierbein or a gravitational metric.

For Q6 with combinatorial Laplacian L, define K_tau = exp(-tau L/2) and
|Psi_tau> = Z_tau^{-1/2} sum_{x,y} K_tau[x,y] |x>_+ |y>_-.
Then rho_+ = K_tau K_tau^dagger / Tr(K_tau K_tau^dagger).

Since spec(L) = {2j, multiplicity C(6,j)}, the Schmidt probabilities are
q_j/Z with q_j=exp(-2 tau j), multiplicity C(6,j), and
Z=(1+exp(-2 tau))^6. Hence S=6 h(p), p=1/(1+exp(2 tau)).
tau is dimensionless; identifying it with a QAT ovalization, proper time,
or physical coupling remains an open calibration problem.
"""
from __future__ import annotations

from math import comb, exp, log
import numpy as np
from scipy.linalg import expm, eigh

from src.q6_complex import build_q6_complex


def q6_laplacian() -> np.ndarray:
    """Return the 64x64 combinatorial graph Laplacian of Q6."""
    q6 = build_q6_complex()
    adjacency = np.zeros((64, 64), dtype=float)
    for u, v in q6.edges:
        adjacency[u, v] = adjacency[v, u] = 1.0
    return 6.0 * np.eye(64) - adjacency


def normalized_bisheet_state(tau: float) -> tuple[np.ndarray, np.ndarray]:
    """Return coefficient matrix C and reduced density matrix rho_plus.

    The bipartite pure state is represented by
        |Psi> = sum_ij C[i,j] |i>_+ |j>_-.
    """
    tau = float(tau)
    if not np.isfinite(tau) or tau < 0.0:
        raise ValueError("tau must be finite and non-negative")
    L = q6_laplacian()
    K = expm(-0.5 * tau * L)
    norm2 = float(np.linalg.norm(K, ord="fro") ** 2)
    C = K / np.sqrt(norm2)
    rho = C @ C.conj().T
    rho = 0.5 * (rho + rho.conj().T)
    return C, rho


def von_neumann_entropy(rho: np.ndarray, atol: float = 1e-12) -> float:
    """Compute -Tr(rho log rho), rejecting non-density matrices."""
    rho = np.asarray(rho, dtype=complex)
    if rho.ndim != 2 or rho.shape[0] != rho.shape[1]:
        raise ValueError("rho must be square")
    if not np.allclose(rho, rho.conj().T, atol=atol):
        raise ValueError("rho must be Hermitian")
    values = np.linalg.eigvalsh(rho)
    if values.min() < -atol or not np.isclose(values.sum(), 1.0, atol=atol):
        raise ValueError("rho must be positive semidefinite with unit trace")
    positive = values[values > atol]
    return float(-np.sum(positive * np.log(positive)))


def analytic_q6_entropy(tau: float) -> float:
    """Exact spectral formula S(tau)=6 h(1/(1+exp(2 tau)))."""
    tau = float(tau)
    if not np.isfinite(tau) or tau < 0.0:
        raise ValueError("tau must be finite and non-negative")
    p = 1.0 / (1.0 + exp(2.0 * tau))
    if p == 0.0 or p == 1.0:
        return 0.0
    return float(-6.0 * (p * log(p) + (1.0 - p) * log(1.0 - p)))


def analytic_schmidt_spectrum(tau: float) -> np.ndarray:
    """Return all 64 Schmidt probabilities with Q6 multiplicities."""
    tau = float(tau)
    if not np.isfinite(tau) or tau < 0.0:
        raise ValueError("tau must be finite and non-negative")
    weights = []
    for j in range(7):
        value = exp(-2.0 * tau * j)
        weights.extend([value] * comb(6, j))
    weights = np.asarray(weights, dtype=float)
    return weights / weights.sum()


def audit_q6_bisheet_entanglement(tau_values=(0.0, 0.1, 0.25, 0.5, 1.0, 2.0)) -> dict:
    """Compare matrix calculation with the exact Q6 spectral entropy."""
    rows = []
    max_entropy_residual = 0.0
    max_spectrum_residual = 0.0
    for tau in tau_values:
        C, rho = normalized_bisheet_state(float(tau))
        numeric_entropy = von_neumann_entropy(rho)
        exact_entropy = analytic_q6_entropy(float(tau))
        numeric_spectrum = np.linalg.eigvalsh(rho)[::-1]
        exact_spectrum = np.sort(analytic_schmidt_spectrum(float(tau)))[::-1]
        entropy_residual = abs(numeric_entropy - exact_entropy)
        spectrum_residual = float(np.max(np.abs(numeric_spectrum - exact_spectrum)))
        max_entropy_residual = max(max_entropy_residual, entropy_residual)
        max_spectrum_residual = max(max_spectrum_residual, spectrum_residual)
        rows.append({
            "tau": float(tau),
            "entropy_numeric_nats": numeric_entropy,
            "entropy_exact_nats": exact_entropy,
            "purity": float(np.trace(rho @ rho).real),
            "mutual_information_pure_state_nats": float(2.0 * numeric_entropy),
            "entropy_residual": entropy_residual,
            "schmidt_spectrum_max_residual": spectrum_residual,
            "rho_trace_residual": abs(float(np.trace(rho).real) - 1.0),
            "rho_min_eigenvalue": float(np.linalg.eigvalsh(rho).min()),
        })
    return {
        "audit": "Q6_DOUBLE_SHEET_HEAT_KERNEL_ENTANGLEMENT",
        "dimensions": {"sheet_dimension": 64, "bipartite_dimension": 4096},
        "kernel": "K_tau = exp(-tau L_Q6 / 2)",
        "state": "vec(K_tau) / ||K_tau||_F",
        "exact_entropy": "6 h(1/(1+exp(2 tau)))",
        "parameter_status": "dimensionless model parameter; not yet calibrated to QAT geometry or physical time",
        "rows": rows,
        "max_entropy_residual": max_entropy_residual,
        "max_schmidt_spectrum_residual": max_spectrum_residual,
        "pass": bool(max_entropy_residual < 1e-10 and max_spectrum_residual < 1e-10),
    }

"""Exact U(1) character expansion and integer electric flux on Q6.

The temporal Wilson factor is expanded as exp(beta*cos(F)) =
sum_{n in Z} I_n(beta) exp(i*n*F). Integrating temporal link angles
projects onto the integer cycle lattice ker(D0.T) over Z. Q6 x C_Nt
remains a 6+1D Euclidean lattice; this is not a 3+1D reduction.
"""
from __future__ import annotations

import numpy as np
from scipy.special import ive

from src.rgs_q6_euclidean_time_wilson import build_q6_euclidean_time_data


def q6_u1_bessel_coefficients(beta: float, nmax: int) -> np.ndarray:
    """Return I_n(beta), n=0,...,nmax, using exponentially scaled Bessel I."""
    if not np.isfinite(beta) or beta <= 0:
        raise ValueError("beta must be finite and positive")
    if not isinstance(nmax, (int, np.integer)) or nmax < 0:
        raise ValueError("nmax must be a nonnegative integer")
    # ive(n,beta)=exp(-beta) I_n(beta); restoring the common factor
    # keeps the coefficients well-scaled for the moderate beta used here.
    return np.asarray(ive(np.arange(nmax + 1), beta) * np.exp(beta), dtype=float)


def q6_u1_character_kernel(beta: float, angles, nmax: int) -> np.ndarray:
    """Truncated symmetric Fourier-character sum for exp(beta*cos(angle))."""
    theta = np.asarray(angles, dtype=float)
    if not np.all(np.isfinite(theta)):
        raise ValueError("angles must be finite")
    coeff = q6_u1_bessel_coefficients(beta, nmax)
    out = coeff[0] * np.ones_like(theta, dtype=complex)
    for n in range(1, nmax + 1):
        out += coeff[n] * (np.exp(1j * n * theta) + np.exp(-1j * n * theta))
    return out


def q6_u1_integer_gauss_law(electric_flux, data=None) -> np.ndarray:
    """Return D0.T@n as exact integer vertex divergence."""
    d = build_q6_euclidean_time_data(2) if data is None else data
    n = np.asarray(electric_flux)
    if n.shape != (192,) or not np.all(np.isfinite(n)):
        raise ValueError("electric_flux must be a finite vector of shape (192,)")
    if not np.all(n == np.rint(n)):
        raise ValueError("electric_flux must be integer-valued")
    n = n.astype(np.int64)
    incidence = np.asarray(d["D0"], dtype=np.int64)
    return incidence.T @ n


def build_q6_integer_cycle_basis(data=None) -> np.ndarray:
    """Build 129 integral fundamental cycles from a spanning tree of Q6.

    Rows are edge-flux vectors in the stored edge orientation and satisfy
    D0.T @ basis.T == 0 exactly over the integers.
    """
    d = build_q6_euclidean_time_data(2) if data is None else data
    edges = d["complex"]["edges"]
    adjacency = [[] for _ in range(64)]
    for ei, (a, b, _axis) in enumerate(edges):
        adjacency[a].append((b, ei))
        adjacency[b].append((a, ei))

    parent = [-1] * 64
    parent_edge = [-1] * 64
    seen = {0}
    queue = [0]
    tree_edges = set()
    for v in queue:
        for w, ei in adjacency[v]:
            if w not in seen:
                seen.add(w)
                parent[w] = v
                parent_edge[w] = ei
                tree_edges.add(ei)
                queue.append(w)
    if len(seen) != 64 or len(tree_edges) != 63:
        raise ArithmeticError("Q6 graph is not connected")

    def tree_path(start, goal):
        # Unique tree path, returned as oriented vertex pairs.
        ancestors = {}
        v = start
        while v != -1:
            ancestors[v] = True
            v = parent[v]
        v = goal
        goal_chain = []
        while v not in ancestors:
            goal_chain.append(v)
            v = parent[v]
        lca = v
        start_chain = []
        v = start
        while v != lca:
            start_chain.append(v)
            v = parent[v]
        vertices = start_chain + [lca] + list(reversed(goal_chain))
        return list(zip(vertices[:-1], vertices[1:]))

    cycles = []
    for chord, (a, b, _axis) in enumerate(edges):
        if chord in tree_edges:
            continue
        # Chord carries unit flow a -> b; tree path returns b -> a.
        cycle = np.zeros(192, dtype=np.int64)
        cycle[chord] = 1
        for u, v in tree_path(b, a):
            key = (u, v) if u < v else (v, u)
            ei = d["complex"]["edge_index"][key]
            cycle[ei] += 1 if u < v else -1
        cycles.append(cycle)
    basis = np.asarray(cycles, dtype=np.int64)
    if basis.shape != (129, 192):
        raise ArithmeticError(f"unexpected cycle basis shape {basis.shape}")
    if np.any(np.asarray(d["D0"], dtype=np.int64).T @ basis.T):
        raise ArithmeticError("constructed cycle basis violates integer Gauss law")
    return basis


def q6_u1_electric_weights(beta_time: float, nmax: int) -> np.ndarray:
    """Single-link temporal character weights I_n(beta_time), n=0..nmax."""
    return q6_u1_bessel_coefficients(beta_time, nmax)


def audit_q6_u1_character_electric_flux(beta=2.3, nmax=18):
    d = build_q6_euclidean_time_data(2)
    basis = build_q6_integer_cycle_basis(d)
    divergence = np.asarray(d["D0"], dtype=np.int64).T @ basis.T
    angles = np.linspace(-np.pi, np.pi, 257)
    approx = q6_u1_character_kernel(beta, angles, nmax)
    exact = np.exp(beta * np.cos(angles))
    bessel_residual = float(np.max(np.abs(approx - exact)))
    relative_residual = float(np.max(np.abs(approx - exact) / exact))
    checks = {
        "cycle_basis_shape_129x192": basis.shape == (129, 192),
        "integer_cycle_basis": bool(np.issubdtype(basis.dtype, np.integer)),
        "integer_gauss_law_exact": bool(np.max(np.abs(divergence)) == 0),
        "character_expansion_converges": relative_residual < 1e-10,
        "positive_character_weights": bool(np.all(q6_u1_electric_weights(beta, nmax) > 0)),
    }
    return {
        "audit": "Q6_U1_CHARACTER_EXPANSION_INTEGER_ELECTRIC_FLUX",
        "beta_time": float(beta),
        "fourier_cutoff": int(nmax),
        "n_vertices": 64,
        "n_edges": 192,
        "rank_D0": 63,
        "integer_cycle_rank": int(basis.shape[0]),
        "bessel_absolute_residual": bessel_residual,
        "bessel_relative_residual": relative_residual,
        "integer_gauss_max_abs": int(np.max(np.abs(divergence))),
        "checks": checks,
        "all_checks_passed": bool(all(checks.values())),
        "exact_statement": (
            "Integrating each compact temporal link with normalized Haar measure "
            "imposes D0.T@n=0 for integer electric flux n in Z^192."
        ),
        "transfer_weight": "product_e I_{n_e}(beta_time)",
        "scope": (
            "Exact character expansion and integer Gauss projector on finite Q6; "
            "not a physical 3+1D reduction or a derived continuum Maxwell limit."
        ),
        "open": [
            "derive anisotropic continuum scaling including temporal lattice spacing",
            "derive beta_s and beta_t from the RGS-PSU action",
            "construct and validate the transfer operator with the spatial plaquette potential",
            "identify a physically justified three-dimensional spatial subcomplex",
        ],
    }

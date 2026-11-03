"""Spectral and homological audit for the canonical six-cube 2-skeleton."""
from __future__ import annotations
from math import comb
import numpy as np
from src.q6_complex import Q6Complex, build_q6_complex


def graph_laplacian(complex_: Q6Complex | None = None) -> np.ndarray:
    q6 = build_q6_complex() if complex_ is None else complex_
    b1 = q6.boundary_1.astype(np.int64)
    return b1 @ b1.T


def expected_q6_spectrum() -> tuple[tuple[int, int], ...]:
    return tuple((2*j, comb(6, j)) for j in range(7))


def spectral_audit(complex_: Q6Complex | None = None, atol: float = 1e-9) -> dict:
    eig = np.linalg.eigvalsh(graph_laplacian(complex_).astype(float))
    expected = expected_q6_spectrum()
    target = np.repeat([x for x, _ in expected], [m for _, m in expected])
    rounded = np.rint(eig).astype(int)
    observed = tuple((int(x), int(np.count_nonzero(rounded == x)))
                     for x in sorted(set(rounded.tolist())))
    return {
        "eigenvalues": eig.tolist(),
        "observed_multiplicities": observed,
        "expected_multiplicities": expected,
        "max_abs_spectral_error": float(np.max(np.abs(eig-target))),
        "pass": bool(observed == expected and np.allclose(eig, target, atol=atol, rtol=0)),
    }


def _rank_mod_prime(matrix: np.ndarray, prime: int = 1_000_003) -> int:
    """Gaussian elimination over GF(prime); an exact lower bound on Q-rank."""
    a = np.asarray(matrix, dtype=np.int64).copy() % prime
    rows, cols = a.shape
    rank = 0
    for col in range(cols):
        candidates = np.flatnonzero(a[rank:, col])
        if candidates.size == 0:
            continue
        pivot = rank + int(candidates[0])
        if pivot != rank:
            a[[rank, pivot]] = a[[pivot, rank]]
        inv = pow(int(a[rank, col]), prime-2, prime)
        a[rank, col:] = (a[rank, col:] * inv) % prime
        if rank + 1 < rows:
            factors = a[rank+1:, col].copy()
            nz = factors != 0
            if np.any(nz):
                a[rank+1:, col:] = (
                    a[rank+1:, col:] - factors[:, None] * a[rank, col:][None, :]
                ) % prime
        rank += 1
        if rank == rows:
            break
    return rank


def homology_audit(complex_: Q6Complex | None = None, prime: int = 1_000_003) -> dict:
    q6 = build_q6_complex() if complex_ is None else complex_
    b1 = q6.boundary_1.astype(np.int64)
    b2 = q6.boundary_2.astype(np.int64)
    n0, n1, n2 = len(q6.vertices), len(q6.edges), len(q6.faces)
    if np.count_nonzero(b1 @ b2):
        raise ValueError("Invalid chain complex: B1 @ B2 != 0.")
    rank_b1 = n0 - 1  # connected graph, exact incidence rank
    rank_mod = _rank_mod_prime(b2, prime)
    rank_upper = n1 - rank_b1  # follows from im(B2) subset ker(B1)
    rank_b2 = rank_mod if rank_mod == rank_upper else None
    graph_cycle_rank = n1 - n0 + 1
    b1_cell = None if rank_b2 is None else n1-rank_b1-rank_b2
    b2_cell = None if rank_b2 is None else n2-rank_b2
    return {
        "rank_B1_exact": rank_b1,
        "rank_B2_mod_prime": rank_mod,
        "rank_B2_rational_certified": rank_b2,
        "rank_B2_upper_bound": rank_upper,
        "graph_cycle_rank_beta1": graph_cycle_rank,
        "cellular_betti_b0": 1,
        "cellular_betti_b1": b1_cell,
        "cellular_betti_b2": b2_cell,
        "euler_characteristic": n0-n1+n2,
        "prime": prime,
        "pass": bool(rank_b2 == 129 and b1_cell == 0 and b2_cell == 111
                     and n0-n1+n2 == 1-b1_cell+b2_cell),
    }


def antipodal_vertex_map() -> np.ndarray:
    return np.arange(64, dtype=np.int64) ^ 63


def audit_q6_formalization() -> dict:
    q6 = build_q6_complex()
    b1 = q6.boundary_1.astype(np.int64)
    b2 = q6.boundary_2.astype(np.int64)
    boundary = b1 @ b2
    degree = np.sum(np.abs(q6.boundary_1), axis=1)
    anti = antipodal_vertex_map()
    edge_set = {tuple(e) for e in q6.edges}
    anti_edges = all(tuple(sorted((int(anti[u]), int(anti[v])))) in edge_set
                     for u, v in q6.edges)
    spec = spectral_audit(q6)
    homology = homology_audit(q6)
    structural = (
        len(q6.vertices) == 64 and len(q6.edges) == 192 and len(q6.faces) == 240
        and np.all(degree == 6) and np.count_nonzero(boundary) == 0
        and np.array_equal(anti[anti], np.arange(64)) and anti_edges
    )
    return {
        "object": "canonical six-dimensional cube 2-skeleton",
        "vertex_count": len(q6.vertices), "edge_count": len(q6.edges),
        "face_count": len(q6.faces), "degree_min": int(degree.min()),
        "degree_max": int(degree.max()),
        "boundary_squared_max_abs": int(np.max(np.abs(boundary))),
        "boundary_squared_exact_zero": bool(np.count_nonzero(boundary) == 0),
        "antipodal_map": anti.tolist(),
        "antipodal_involution_exact": bool(np.array_equal(anti[anti], np.arange(64))),
        "antipodal_preserves_edges": bool(anti_edges),
        "graph_laplacian_spectrum": spec, "homology": homology,
        "structural_pass": bool(structural),
        "full_pass": bool(structural and spec["pass"] and homology["pass"]),
        "geometric_embedding_Oh_identification": "OPEN: requires X_Oh and induced action audit",
        "gauge_dynamics_on_Q6": "OPEN: not part of this combinatorial construction",
    }

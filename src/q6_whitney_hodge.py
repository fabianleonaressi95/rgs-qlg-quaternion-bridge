"""Geometric cochain metrics from cubical Whitney forms on the Q6 2-skeleton.

The canonical Q6 complex is the 2-skeleton of the unit six-cube [0,1]^6.
For each p-cell, its cubical Whitney basis function is the tensor product
of endpoint hat functions in fixed directions and dx_i in active directions.
The L2 Gram matrices H_0, H_1, H_2 are assembled from the 1D endpoint mass
matrix [[1/3,1/6],[1/6,1/3]].

This is a metric-derived cochain inner product for the flat unit cube, not a
circumcentric dual-cell Hodge star and not a curved RGS or continuum-limit
claim. Dense matrices are intentional for this 64/192/240 finite audit.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from src.q6_complex import Q6Complex, build_q6_complex

_ONE_D_MASS = np.array([[1.0 / 3.0, 1.0 / 6.0],
                        [1.0 / 6.0, 1.0 / 3.0]], dtype=float)


@dataclass(frozen=True)
class Q6WhitneyMetrics:
    """L2 Gram matrices in the canonical Q6 cell order."""
    metric_0: np.ndarray
    metric_1: np.ndarray
    metric_2: np.ndarray
    orientation_blocks: dict[int, tuple[tuple[int, ...], ...]]
    cell_counts: tuple[int, int, int]


@dataclass(frozen=True)
class Q6WhitneyHodge:
    """Symmetric Hodge Laplacians in orthonormal cochain coordinates."""
    metrics: Q6WhitneyMetrics
    laplacian_0: np.ndarray
    laplacian_1: np.ndarray
    laplacian_2: np.ndarray


def _bits(vertex: int) -> tuple[int, ...]:
    return tuple((vertex >> (5 - axis)) & 1 for axis in range(6))


def _cell_metric(
    cells: list[tuple[tuple[int, ...], tuple[int, ...]]],
) -> tuple[np.ndarray, tuple[tuple[int, ...], ...]]:
    """Assemble the L2 Gram matrix, block diagonal by active-axis orientation."""
    n = len(cells)
    matrix = np.zeros((n, n), dtype=float)
    by_orientation: dict[tuple[int, ...], list[tuple[int, tuple[int, ...]]]] = {}
    for index, (active, fixed_bits) in enumerate(cells):
        by_orientation.setdefault(active, []).append((index, fixed_bits))

    for active, entries in by_orientation.items():
        fixed_axes = tuple(a for a in range(6) if a not in active)
        for i, (row, row_bits) in enumerate(entries):
            left = dict(zip(fixed_axes, row_bits))
            for col, col_bits in entries[: i + 1]:
                right = dict(zip(fixed_axes, col_bits))
                value = 1.0
                for axis in fixed_axes:
                    value *= _ONE_D_MASS[left[axis], right[axis]]
                matrix[row, col] = matrix[col, row] = value
    return matrix, tuple(sorted(by_orientation))


def build_q6_whitney_metrics(
    complex_: Q6Complex | None = None,
) -> Q6WhitneyMetrics:
    """Build Euclidean L2 Gram matrices for cubical Whitney forms."""
    q6 = build_q6_complex() if complex_ is None else complex_
    vertices = [((), _bits(v)) for v in q6.vertices]

    edges = []
    for u, v in q6.edges:
        ubits, vbits = _bits(u), _bits(v)
        active = tuple(a for a in range(6) if ubits[a] != vbits[a])
        if len(active) != 1:
            raise ValueError("Q6 edge is not an axis-aligned unit segment")
        edges.append((active, tuple(ubits[a] for a in range(6) if a not in active)))

    faces = []
    for v, va, vab, vb in q6.faces:
        b0, ba, bb = _bits(v), _bits(va), _bits(vb)
        active = tuple(a for a in range(6) if b0[a] != ba[a] or b0[a] != bb[a])
        if len(active) != 2 or active != tuple(sorted(active)):
            raise ValueError("Q6 face is not an axis-aligned square")
        faces.append((active, tuple(b0[a] for a in range(6) if a not in active)))

    h0, o0 = _cell_metric(vertices)
    h1, o1 = _cell_metric(edges)
    h2, o2 = _cell_metric(faces)
    return Q6WhitneyMetrics(h0, h1, h2, {0: o0, 1: o1, 2: o2},
                            (len(vertices), len(edges), len(faces)))


def build_q6_whitney_hodge(
    metrics: Q6WhitneyMetrics | None = None,
    complex_: Q6Complex | None = None,
) -> Q6WhitneyHodge:
    """Construct metric Hodge Laplacians using full SPD Gram matrices.

    For d0=B1.T and d1=B2.T, adjoints are H_k^{-1} d_k.T H_{k+1}.
    Cholesky coordinates yield symmetric representatives.
    """
    q6 = build_q6_complex() if complex_ is None else complex_
    model = build_q6_whitney_metrics(q6) if metrics is None else metrics
    h0, h1, h2 = model.metric_0, model.metric_1, model.metric_2
    sizes = (len(q6.vertices), len(q6.edges), len(q6.faces))
    if model.cell_counts != sizes:
        raise ValueError("Whitney metric cell counts do not match Q6 complex")
    b1 = np.asarray(q6.boundary_1, dtype=float)
    b2 = np.asarray(q6.boundary_2, dtype=float)
    if b1.shape != (sizes[0], sizes[1]) or b2.shape != (sizes[1], sizes[2]):
        raise ValueError("incidence matrix shapes are inconsistent")
    if np.count_nonzero(b1.astype(np.int64) @ b2.astype(np.int64)):
        raise ValueError("invalid chain complex: B1 @ B2 must vanish")

    l0, l1, l2 = (np.linalg.cholesky(h) for h in (h0, h1, h2))
    inv0, inv1 = np.linalg.inv(l0), np.linalg.inv(l1)
    d0, d1 = b1.T, b2.T
    a0 = l1.T @ d0 @ inv0.T
    a1 = l2.T @ d1 @ inv1.T
    lap0 = a0.T @ a0
    lap1 = a0 @ a0.T + a1.T @ a1
    lap2 = a1 @ a1.T
    return Q6WhitneyHodge(model, (lap0 + lap0.T) / 2.0,
                          (lap1 + lap1.T) / 2.0, (lap2 + lap2.T) / 2.0)


def audit_q6_whitney_hodge(
    hodge: Q6WhitneyHodge | None = None,
    complex_: Q6Complex | None = None,
    atol: float = 1e-9,
) -> dict[str, object]:
    """Audit SPD metrics, symmetry, PSD and harmonic dimensions."""
    q6 = build_q6_complex() if complex_ is None else complex_
    model = build_q6_whitney_hodge(complex_=q6) if hodge is None else hodge
    hs = (model.metrics.metric_0, model.metrics.metric_1, model.metrics.metric_2)
    ls = (model.laplacian_0, model.laplacian_1, model.laplacian_2)
    metric_eigs = tuple(np.linalg.eigvalsh(h) for h in hs)
    lap_eigs = tuple(np.linalg.eigvalsh(a) for a in ls)
    ranks = (int(np.linalg.matrix_rank(q6.boundary_1.astype(float))),
             int(np.linalg.matrix_rank(q6.boundary_2.astype(float))))
    sizes = (len(q6.vertices), len(q6.edges), len(q6.faces))
    betti = (sizes[0] - ranks[0], sizes[1] - ranks[0] - ranks[1],
             sizes[2] - ranks[1])
    harmonic = tuple(int(np.count_nonzero(np.abs(e) <= atol)) for e in lap_eigs)
    symmetry = tuple(float(np.linalg.norm(a - a.T)) for a in ls)
    min_lap = tuple(float(e[0]) for e in lap_eigs)
    min_metric = tuple(float(e[0]) for e in metric_eigs)
    chain_ok = not np.count_nonzero(
        q6.boundary_1.astype(np.int64) @ q6.boundary_2.astype(np.int64))
    passed = (all(v > 0.0 for v in min_metric)
              and all(v >= -atol for v in min_lap)
              and all(v <= atol for v in symmetry)
              and harmonic == betti and chain_ok)
    return {
        "cell_counts": {"C0": sizes[0], "C1": sizes[1], "C2": sizes[2]},
        "metric_min_eigenvalues": min_metric,
        "metric_condition_numbers": tuple(float(np.linalg.cond(h)) for h in hs),
        "laplacian_min_eigenvalues": min_lap,
        "laplacian_symmetry_residuals": symmetry,
        "cellular_betti_numbers": betti,
        "harmonic_dimensions": harmonic,
        "boundary_squared_exact_zero": bool(chain_ok),
        "pass": bool(passed),
        "scope": ("flat unit-cube cubical Whitney L2 metrics on the Q6 2-skeleton; "
                  "not a circumcentric dual, curved RGS Hodge star, continuum "
                  "limit, or Yang-Mills dynamics"),
    }

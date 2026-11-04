"""Metric-derived cubical Hodge operators on the full Euclidean 6-cube.

The old Q6 square complex is only the 2-skeleton. A Hodge star in ambient
R^6 needs complementary degrees, so this module constructs all cubical cells
(k=0,...,6) and the tensor-product lowest-order cubical Whitney mass matrices
for the unit cube [0,1]^6 with its Euclidean metric.

Cochain inner products are assembled L2 inner products of the corresponding
tensor-product differential-form basis. This is a metric-derived discrete
Hodge *inner product*, not a primal-to-dual-cell map. No curved RGS metric,
refinement limit, Yang-Mills dynamics, or physical spacetime is assumed.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from math import comb
import numpy as np


@dataclass(frozen=True)
class CubicalComplex:
    cells: tuple[tuple[int, ...], ...]  # -1=active axis, 0/1=fixed endpoint
    boundary: tuple[np.ndarray, ...]   # B_k: C_k -> C_{k-1}; B_0 is empty


@dataclass(frozen=True)
class CubicalWhitneyHodge:
    complex: CubicalComplex
    mass_matrices: tuple[np.ndarray, ...]
    laplacians: tuple[np.ndarray, ...]  # symmetric representatives in L2-orthonormal coordinates


def build_unit_cube_complex(dimension: int = 6) -> CubicalComplex:
    """Build the oriented full cubical cell complex of [0,1]^dimension."""
    if not isinstance(dimension, int) or dimension < 1 or dimension > 8:
        raise ValueError("dimension must be an integer in 1..8")
    cells_by_degree = []
    indices = []
    for k in range(dimension + 1):
        cells = tuple(
            state for state in product((-1, 0, 1), repeat=dimension)
            if state.count(-1) == k
        )
        cells_by_degree.append(cells)
        indices.append({cell: i for i, cell in enumerate(cells)})

    boundaries = [np.zeros((0, len(cells_by_degree[0])), dtype=np.int8)]
    for k in range(1, dimension + 1):
        lower_index = indices[k - 1]
        matrix = np.zeros((len(cells_by_degree[k - 1]), len(cells_by_degree[k])), dtype=np.int8)
        for column, cell in enumerate(cells_by_degree[k]):
            active_axes = [axis for axis, value in enumerate(cell) if value == -1]
            for position, axis in enumerate(active_axes):
                sign = -1 if position % 2 == 0 else 1
                for endpoint, coefficient in ((0, -sign), (1, sign)):
                    face = list(cell)
                    face[axis] = endpoint
                    matrix[lower_index[tuple(face)], column] += coefficient
        boundaries.append(matrix)
    return CubicalComplex(
        cells=tuple(cells_by_degree),
        boundary=tuple(boundaries),
    )


def cubical_whitney_mass_matrices(complex_: CubicalComplex) -> tuple[np.ndarray, ...]:
    """Assemble exact tensor-product L2 Gram matrices for unit-cube cubical forms.

    For a fixed active-axis set I, basis forms are
      prod_{j not in I} lambda_{endpoint_j}(x_j) dx_I.
    Distinct I are pointwise orthogonal. In each fixed coordinate,
    integral(lambda_a lambda_b) is 1/3 if a=b and 1/6 otherwise.
    """
    dimension = len(complex_.cells[0][0])
    one_dimensional_mass = np.array([[1.0 / 3.0, 1.0 / 6.0],
                                     [1.0 / 6.0, 1.0 / 3.0]])
    masses = []
    for k, cells in enumerate(complex_.cells):
        index = {cell: i for i, cell in enumerate(cells)}
        mass = np.zeros((len(cells), len(cells)), dtype=float)
        blocks: dict[tuple[int, ...], list[tuple[int, tuple[int, ...]]]] = {}
        for i, cell in enumerate(cells):
            active = tuple(axis for axis, value in enumerate(cell) if value == -1)
            fixed = tuple(cell[axis] for axis in range(dimension) if axis not in active)
            blocks.setdefault(active, []).append((i, fixed))
        for entries in blocks.values():
            for i, fixed_i in entries:
                for j, fixed_j in entries:
                    value = 1.0
                    for a, b in zip(fixed_i, fixed_j):
                        value *= one_dimensional_mass[a, b]
                    mass[i, j] = value
        # A structural guard: all cells of a fixed orientation form one SPD block.
        if not np.all(np.isfinite(mass)) or not np.allclose(mass, mass.T, atol=1e-14):
            raise ArithmeticError("invalid cubical Whitney mass matrix")
        masses.append(mass)
    return tuple(masses)


def build_cubical_whitney_hodge(dimension: int = 6) -> CubicalWhitneyHodge:
    """Build all degree Hodge Laplacians using Euclidean cubical Whitney masses."""
    complex_ = build_unit_cube_complex(dimension)
    masses = cubical_whitney_mass_matrices(complex_)
    laplacians = []
    for k, cells in enumerate(complex_.cells):
        hk = masses[k]
        # Similarity transform to an ordinary symmetric matrix:
        # H_k^(1/2) Delta_k H_k^(-1/2).
        evals, vectors = np.linalg.eigh(hk)
        if np.min(evals) <= 0:
            raise ArithmeticError(f"Whitney mass matrix H_{k} is not positive definite")
        sqrt_h = (vectors * np.sqrt(evals)) @ vectors.T
        inv_sqrt_h = (vectors * (1.0 / np.sqrt(evals))) @ vectors.T
        # Lower-degree contribution in L2-orthonormal coordinates:
        # H_k^(1/2) d_(k-1) H_(k-1)^(-1) d_(k-1)^T H_k^(1/2).
        symmetric = np.zeros_like(hk)
        if k > 0:
            d_prev = complex_.boundary[k].astype(float).T  # (n_k, n_(k-1))
            h_prev_inv_dT = np.linalg.solve(masses[k - 1], d_prev.T)
            symmetric += sqrt_h @ d_prev @ h_prev_inv_dT @ sqrt_h
        # Upper-degree contribution:
        # H_k^(-1/2) d_k^T H_(k+1) d_k H_k^(-1/2).
        if k < dimension:
            d_next = complex_.boundary[k + 1].astype(float).T  # (n_(k+1), n_k)
            symmetric += inv_sqrt_h @ d_next.T @ masses[k + 1] @ d_next @ inv_sqrt_h
        laplacians.append((symmetric + symmetric.T) / 2.0)
    return CubicalWhitneyHodge(complex_, masses, tuple(laplacians))


def audit_cubical_whitney_hodge(model: CubicalWhitneyHodge | None = None,
                                 dimension: int = 6,
                                 atol: float = 1e-9) -> dict[str, object]:
    """Audit cell counts, chain identities, mass positivity and harmonic kernels."""
    result = build_cubical_whitney_hodge(dimension) if model is None else model
    complex_ = result.complex
    dim = len(complex_.cells) - 1
    chain_residuals = []
    for k in range(2, dim + 1):
        residual = (complex_.boundary[k - 1].astype(np.int64)
                    @ complex_.boundary[k].astype(np.int64))
        chain_residuals.append(int(np.max(np.abs(residual))) if residual.size else 0)
    mass_min = [float(np.linalg.eigvalsh(h)[0]) for h in result.mass_matrices]
    spectra = [np.linalg.eigvalsh(lap) for lap in result.laplacians]
    kernels = [int(np.count_nonzero(np.abs(values) <= atol)) for values in spectra]
    ranks = [int(np.linalg.matrix_rank(b.astype(float))) for b in complex_.boundary[1:]]
    betti = []
    for k, cells in enumerate(complex_.cells):
        rank_lower = ranks[k - 1] if k > 0 else 0
        rank_upper = ranks[k] if k < dim else 0
        betti.append(len(cells) - rank_lower - rank_upper)
    expected_counts = [comb(dim, k) * (2 ** (dim - k)) for k in range(dim + 1)]
    chain_ok = all(value == 0 for value in chain_residuals)
    pass_ = (chain_ok and [len(x) for x in complex_.cells] == expected_counts
             and all(value > 0 for value in mass_min)
             and kernels == betti
             and all(float(values[0]) >= -atol for values in spectra))
    return {
        "dimension": dim,
        "cell_counts_by_degree": [len(x) for x in complex_.cells],
        "expected_cell_counts": expected_counts,
        "boundary_ranks": ranks,
        "betti_numbers": betti,
        "harmonic_dimensions": kernels,
        "mass_min_eigenvalues": mass_min,
        "laplacian_min_eigenvalues": [float(values[0]) for values in spectra],
        "boundary_squared_max_abs": chain_residuals,
        "pass": bool(pass_),
        "metric": "Euclidean metric on the unit cube [0,1]^d",
        "discretisation": "tensor-product lowest-order cubical Whitney differential forms",
        "scope": (
            "full finite cubical complex and its metric-derived L2 cochain inner products; "
            "not a primal-dual cell Hodge-star map, curved RGS geometry, refinement limit, "
            "Mosco convergence, Yang-Mills dynamics, or physical spacetime"
        ),
    }

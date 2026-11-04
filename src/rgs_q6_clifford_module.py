"""Full cubical Q6 cochains and the six-dimensional exterior Clifford module.

The 6-cube is treated as a genuine cubical cell complex, not as a periodic
lattice. Cells are oriented by increasing coordinate directions. The exterior
algebra Lambda*(R^6) supplies a 64-dimensional complex Clifford module.
"""
from __future__ import annotations

from itertools import combinations
import numpy as np


def _cells_by_degree(n=6):
    cells = []
    index = []
    for k in range(n + 1):
        layer = []
        lookup = {}
        for dirs in combinations(range(n), k):
            free_mask = sum(1 << d for d in dirs)
            for base in range(1 << n):
                if base & free_mask:
                    continue
                lookup[(base, dirs)] = len(layer)
                layer.append((base, dirs))
                cells.append((base, dirs))
        index.append(lookup)
    return cells, index


def build_full_q6_cubical_complex():
    """Return all Q6 cells, oriented boundaries and coboundaries."""
    cells, index = _cells_by_degree(6)
    counts = [sum(1 for _ in combinations(range(6), k)) * (2 ** (6-k))
              for k in range(7)]
    boundaries = [None]
    coboundaries = [None]
    for k in range(1, 7):
        b = np.zeros((counts[k-1], counts[k]), dtype=np.int8)
        for col, (base, dirs) in enumerate(
                cell for cell in cells if len(cell[1]) == k):
            for p, axis in enumerate(dirs):
                lower_dirs = dirs[:p] + dirs[p+1:]
                lower = index[k-1][(base, lower_dirs)]
                upper = index[k-1][(base | (1 << axis), lower_dirs)]
                sign = -1 if p % 2 else 1
                b[upper, col] += sign
                b[lower, col] -= sign
        boundaries.append(b)
        coboundaries.append(b.T.copy())
    offsets = np.cumsum([0] + counts)
    total = int(offsets[-1])
    d = np.zeros((total, total), dtype=np.int8)
    for k in range(6):
        block = coboundaries[k+1]
        # D_k maps C^k -> C^(k+1): higher-degree cells are rows,
        # lower-degree cells are columns. Keep the block orientation
        # consistent with coboundaries[k+1].shape == (n_{k+1}, n_k).
        d[offsets[k+1]:offsets[k+2], offsets[k]:offsets[k+1]] = block
    dk = d.astype(float) + d.T.astype(float)
    return {
        "cells": tuple(cells), "cell_counts": tuple(counts),
        "offsets": offsets, "boundaries": tuple(boundaries),
        "coboundaries": tuple(coboundaries), "d": d, "D_K": dk,
    }


def exterior_clifford_generators(n=6):
    """Hermitian gamma_j = wedge(e_j) + contraction(e_j) on Lambda*(R^n)."""
    if not isinstance(n, int) or n < 1 or n > 10:
        raise ValueError("n must be an integer in [1, 10]")
    dim = 1 << n
    gammas = []
    for axis in range(n):
        g = np.zeros((dim, dim), dtype=complex)
        bit = 1 << axis
        lower_mask = bit - 1
        for state in range(dim):
            occupied = bool(state & bit)
            sign = -1 if ((state & lower_mask).bit_count() % 2) else 1
            target = state ^ bit
            g[target, state] = sign
        gammas.append(g)
    return tuple(gammas)


def form_degree_chirality(n=6):
    """Z2 grading Gamma=(-1)^degree; anticommutes with each gamma_j."""
    dim = 1 << n
    return np.diag([1.0 if state.bit_count() % 2 == 0 else -1.0
                    for state in range(dim)]).astype(complex)


def audit_q6_clifford_module():
    c = build_full_q6_cubical_complex()
    gammas = exterior_clifford_generators(6)
    chi = form_degree_chirality(6)
    eye = np.eye(64, dtype=complex)
    clifford = hermitian = chirality = 0.0
    for i, gi in enumerate(gammas):
        hermitian = max(hermitian, float(np.linalg.norm(gi-gi.conj().T)))
        chirality = max(chirality, float(np.linalg.norm(chi@gi+gi@chi)))
        for j, gj in enumerate(gammas):
            target = 2.0*eye if i == j else np.zeros_like(eye)
            clifford = max(clifford, float(np.linalg.norm(gi@gj+gj@gi-target)))
    d2 = 0.0
    for k in range(1, 6):
        product = c["coboundaries"][k+1] @ c["coboundaries"][k]
        d2 = max(d2, float(np.max(np.abs(product))) if product.size else 0.0)
    counts = list(c["cell_counts"])
    ranks = [int(np.linalg.matrix_rank(c["boundaries"][k].astype(float)))
             for k in range(1, 7)]
    return {
        "cell_counts_by_degree": counts,
        "total_cells": int(sum(counts)),
        "clifford_module_dimension": 64,
        "clifford_anticommutator_residual": clifford,
        "gamma_hermiticity_residual": hermitian,
        "chirality_anticommutation_residual": chirality,
        "coboundary_squared_residual": d2,
        "boundary_ranks": ranks,
        "expected_boundary_ranks": [63, 129, 111, 49, 11, 1],
        "all_checks_passed": bool(
            counts == [64, 192, 240, 160, 60, 12, 1]
            and clifford < 1e-12 and hermitian < 1e-12
            and chirality < 1e-12 and d2 == 0
            and ranks == [63, 129, 111, 49, 11, 1]
        ),
    }

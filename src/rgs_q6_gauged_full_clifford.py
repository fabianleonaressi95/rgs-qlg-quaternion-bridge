"""Gauge-covariant full cubical Q6 Dirac--Kaehler operator.

This module couples the complete 6-cube cochain complex (729 cells in
degrees 0..6) to SU(2) parallel transport on its 192 oriented edges.
A cell fibre is anchored at the cell's lower vertex. In each cubical
coboundary, upper-face data are parallel transported back to that anchor.

The operator is a finite gauge-covariant Dirac--Kaehler operator. The
64-dimensional exterior Clifford module is kept as an explicit local
Clifford factor; this construction does NOT identify the SU(2) gauge
connection with a Spin(6) connection and does not claim a spin-Dirac
continuum limit.
"""
from __future__ import annotations

import numpy as np
from scipy import sparse

from src.rgs_q6_clifford_module import (
    build_full_q6_cubical_complex,
    exterior_clifford_generators,
    form_degree_chirality,
)
from src.rgs_q6_equivariant_cochains import build_q6_cochain_complex
from src.rgs_q6_su2_connection import (
    random_su2_matrices,
    validate_su2_field,
    gauge_transform_links,
    plaquette_holonomies,
)


def _cell_layers(complex_data):
    layers = [[] for _ in range(7)]
    for cell in complex_data["cells"]:
        layers[len(cell[1])].append(cell)
    return layers


def _cell_lookups(layers):
    return [{cell: i for i, cell in enumerate(layer)} for layer in layers]


def _edge_lookup(graph):
    return {
        tuple(sorted((int(edge[0]), int(edge[1])))): i
        for i, edge in enumerate(graph["edges"])
    }


def covariant_full_cubical_differential(links, cube=None, graph=None):
    """Return sparse D_U:C^k(Q6,C^2)->C^(k+1)(Q6,C^2), all degrees.

    For a (k+1)-cell c=(base, dirs), each lower face is anchored at base
    and each upper face at base+e_axis. The latter is transported to base
    by U_(base,base+e_axis). Thus D_U' = P D_U P^dagger exactly.
    """
    cube = build_full_q6_cubical_complex() if cube is None else cube
    graph = build_q6_cochain_complex() if graph is None else graph
    u = np.asarray(links, dtype=complex)
    validate_su2_field(u, len(graph["edges"]))
    layers = _cell_layers(cube)
    lookup = _cell_lookups(layers)
    edge_index = _edge_lookup(graph)
    offsets = np.asarray(cube["offsets"], dtype=int)
    n_cells = int(offsets[-1])
    dim = 2 * n_cells
    out = sparse.lil_matrix((dim, dim), dtype=complex)
    eye = np.eye(2, dtype=complex)

    for k in range(6):
        row_offset = int(offsets[k + 1])
        col_offset = int(offsets[k])
        for row_local, (base, dirs) in enumerate(layers[k + 1]):
            row = 2 * (row_offset + row_local)
            for p, axis in enumerate(dirs):
                lower_dirs = dirs[:p] + dirs[p + 1 :]
                lower = lookup[k][(base, lower_dirs)]
                upper_base = base | (1 << axis)
                upper = lookup[k][(upper_base, lower_dirs)]
                col_lower = 2 * (col_offset + lower)
                col_upper = 2 * (col_offset + upper)
                sign = -1.0 if p % 2 else 1.0
                edge_id = edge_index[(int(base), int(upper_base))]
                # U_(base,upper_base) maps the upper-face fibre back to base.
                out[row : row + 2, col_lower : col_lower + 2] += -sign * eye
                out[row : row + 2, col_upper : col_upper + 2] += sign * u[edge_id]

    return out.tocsr()


def cell_gauge_representation(gauge, cube=None):
    """Block-diagonal gauge action, one G_base per cell anchor."""
    cube = build_full_q6_cubical_complex() if cube is None else cube
    g = np.asarray(gauge, dtype=complex)
    validate_su2_field(g, 64)
    blocks = [g[int(base)] for base, _dirs in cube["cells"]]
    return sparse.block_diag(blocks, format="csr", dtype=complex)


def gauged_full_dirac_kaehler_operator(links, cube=None, graph=None):
    """Hermitian D_U + D_U^dagger on all cubical cochains with SU(2) fibre."""
    d = covariant_full_cubical_differential(links, cube=cube, graph=graph)
    return (d + d.getH()).tocsr()


def _max_abs_sparse(matrix):
    if matrix.nnz == 0:
        return 0.0
    return float(np.max(np.abs(matrix.data)))


def audit_q6_gauged_full_clifford(seed=20261104):
    """Audit all-degree transport, gauge covariance, curvature and Clifford factor."""
    cube = build_full_q6_cubical_complex()
    graph = build_q6_cochain_complex()
    links = random_su2_matrices(len(graph["edges"]), seed=seed)
    gauge = random_su2_matrices(64, seed=seed + 1)
    validate_su2_field(links, 192)
    validate_su2_field(gauge, 64)

    d = covariant_full_cubical_differential(links, cube, graph)
    dk = (d + d.getH()).tocsr()
    links_g = gauge_transform_links(links, gauge, graph)
    dg = covariant_full_cubical_differential(links_g, cube, graph)
    p = cell_gauge_representation(gauge, cube)
    covariance = _max_abs_sparse(dg - p @ d @ p.getH())
    hermiticity = _max_abs_sparse(dk - dk.getH())

    # Curvature is the obstruction to D_U^2=0. Flat identity links recover
    # the ordinary cubical cochain identity exactly.
    identity = np.repeat(np.eye(2, dtype=complex)[None, :, :], 192, axis=0)
    d_flat = covariant_full_cubical_differential(identity, cube, graph)
    flat_d2 = _max_abs_sparse(d_flat @ d_flat)
    curvature_d2 = _max_abs_sparse(d @ d)

    h = plaquette_holonomies(links, graph)
    hg = plaquette_holonomies(links_g, graph)
    trace_residual = float(np.max(np.abs(
        np.trace(h, axis1=1, axis2=2) - np.trace(hg, axis1=1, axis2=2)
    )))
    face_conjugation = 0.0
    for i, (_v, _a, _b, cycle) in enumerate(graph["faces"]):
        expected = gauge[cycle[0]] @ h[i] @ gauge[cycle[0]].conj().T
        face_conjugation = max(
            face_conjugation, float(np.linalg.norm(hg[i] - expected))
        )

    gammas = exterior_clifford_generators(6)
    chi = form_degree_chirality(6)
    eye2 = np.eye(2, dtype=complex)
    eye64 = np.eye(64, dtype=complex)
    # The internal SU(2) transporter acts trivially on the local exterior
    # Clifford factor; hence it commutes with every I_2 tensor gamma_j.
    test_link = links[0]
    lifted_link = np.kron(test_link, eye64)
    clifford_commutator = max(
        float(np.linalg.norm(
            lifted_link @ np.kron(eye2, gamma)
            - np.kron(eye2, gamma) @ lifted_link
        ))
        for gamma in gammas
    )
    chirality_residual = max(
        float(np.linalg.norm(chi @ gamma + gamma @ chi))
        for gamma in gammas
    )

    counts = list(cube["cell_counts"])
    passed = (
        counts == [64, 192, 240, 160, 60, 12, 1]
        and d.shape == (1458, 1458)
        and covariance < 1e-10
        and hermiticity < 1e-12
        and flat_d2 < 1e-12
        and curvature_d2 > 1e-8
        and trace_residual < 1e-10
        and face_conjugation < 1e-10
        and clifford_commutator < 1e-12
        and chirality_residual < 1e-12
    )
    return {
        "audit": "Q6_GAUGED_FULL_CUBICAL_DIRAC_KAEHLER_CLIFFORD",
        "cell_counts_by_degree": counts,
        "total_cells": int(sum(counts)),
        "edge_links": 192,
        "square_holonomies": 240,
        "su2_fiber_dimension": 2,
        "cochain_differential_shape": list(d.shape),
        "hermitian_dirac_kaehler_shape": list(dk.shape),
        "gauge_covariance_max_abs_residual": covariance,
        "hermiticity_max_abs_residual": hermiticity,
        "flat_connection_D_squared_max_abs": flat_d2,
        "curved_connection_D_squared_max_abs": curvature_d2,
        "wilson_trace_gauge_invariance_max_residual": trace_residual,
        "plaquette_conjugation_max_residual": face_conjugation,
        "su2_clifford_factor_commutator_residual": clifford_commutator,
        "clifford_chirality_anticommutator_residual": chirality_residual,
        "interpretation": {
            "proved_by_construction": [
                "cell-anchor gauge representation",
                "sparse all-degree gauge-covariant cubical coboundary",
                "Hermiticity of D_U + D_U^dagger",
            ],
            "numerically_verified": [
                "local SU(2) covariance",
                "240 plaquette conjugation and Wilson trace invariance",
                "flat D_U squared equals zero",
                "curved D_U squared is nonzero",
                "SU(2) transport commutes with the independent exterior Clifford factor",
            ],
            "not_claimed": [
                "Spin(6) or Spin(1,3) connection",
                "identification of internal SU(2) with the spin group",
                "continuum Dirac limit or fermion-doubling resolution",
                "physical Yang-Mills coupling derived from RGS data",
            ],
        },
        "all_checks_passed": bool(passed),
    }

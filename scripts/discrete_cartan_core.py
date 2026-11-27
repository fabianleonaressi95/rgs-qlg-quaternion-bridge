"""Minimal finite-complex Cartan core.

Convention:
- L[(v,w)] maps vectors in the frame at w to the frame at v.
- E[(v,w)] is the oriented edge coframe in the frame at v.
- Reverse edges are derived from L[(w,v)] = inv(L[(v,w)]) and
  E[(w,v)] = -inv(L[(v,w)]) @ E[(v,w)].

This module deliberately does not implement a continuum limit or a discrete
Bianchi theorem. It provides the local algebraic core for covariance tests.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence
import numpy as np
from numpy.typing import NDArray

Matrix = NDArray[np.float64]
Vector = NDArray[np.float64]
Edge = tuple[int, int]


def lorentz_boost_x(rapidity: float) -> Matrix:
    """Proper-orthochronous boost in the x direction, signature (+---)."""
    c = float(np.cosh(rapidity))
    s = float(np.sinh(rapidity))
    return np.array(
        [[c, s, 0.0, 0.0],
         [s, c, 0.0, 0.0],
         [0.0, 0.0, 1.0, 0.0],
         [0.0, 0.0, 0.0, 1.0]],
        dtype=float,
    )


def inverse_edge(
    edge: Edge,
    links: Mapping[Edge, Matrix],
    coframes: Mapping[Edge, Vector],
) -> tuple[Matrix, Vector]:
    """Return (L_wv, E_wv) for an edge (w,v), deriving from stored (v,w)."""
    v, w = edge
    if (v, w) in links and (v, w) in coframes:
        return np.asarray(links[(v, w)]), np.asarray(coframes[(v, w)])
    if (w, v) not in links or (w, v) not in coframes:
        raise KeyError(f"missing edge data for {edge} or its reverse")
    L = np.asarray(links[(w, v)], dtype=float)
    E = np.asarray(coframes[(w, v)], dtype=float)
    Linv = np.linalg.inv(L)
    return Linv, -Linv @ E


def edge_data(
    edge: Edge,
    links: Mapping[Edge, Matrix],
    coframes: Mapping[Edge, Vector],
) -> tuple[Matrix, Vector]:
    """Fetch oriented edge data, allowing only one stored orientation."""
    v, w = edge
    if (v, w) in links and (v, w) in coframes:
        return np.asarray(links[(v, w)], dtype=float), np.asarray(coframes[(v, w)], dtype=float)
    return inverse_edge(edge, links, coframes)


def face_holonomy(
    face: Sequence[int],
    links: Mapping[Edge, Matrix],
    coframes: Mapping[Edge, Vector] | None = None,
) -> Matrix:
    """Based face holonomy for an ordered triangular face (v0,v1,v2)."""
    if len(face) != 3:
        raise ValueError("face_holonomy expects an ordered triangle")
    v0, v1, v2 = map(int, face)
    L01, _ = edge_data((v0, v1), links, coframes or {})
    L12, _ = edge_data((v1, v2), links, coframes or {})
    L20, _ = edge_data((v2, v0), links, coframes or {})
    return L01 @ L12 @ L20


def face_torsion(
    face: Sequence[int],
    links: Mapping[Edge, Matrix],
    coframes: Mapping[Edge, Vector],
) -> Vector:
    """Transported coframe closure vector, based at the first face vertex."""
    if len(face) != 3:
        raise ValueError("face_torsion expects an ordered triangle")
    v0, v1, v2 = map(int, face)
    L01, E01 = edge_data((v0, v1), links, coframes)
    L12, E12 = edge_data((v1, v2), links, coframes)
    L20, E20 = edge_data((v2, v0), links, coframes)
    return E01 + L01 @ E12 + L01 @ L12 @ E20


def gauge_transform(
    links: Mapping[Edge, Matrix],
    coframes: Mapping[Edge, Vector],
    frames: Mapping[int, Matrix],
) -> tuple[dict[Edge, Matrix], dict[Edge, Vector]]:
    """Apply L_vw -> h_v L_vw h_w^-1 and E_vw -> h_v E_vw."""
    out_links: dict[Edge, Matrix] = {}
    out_coframes: dict[Edge, Vector] = {}
    for (v, w), L in links.items():
        if (v, w) not in coframes:
            raise KeyError(f"missing coframe for stored edge {(v, w)}")
        hv = np.asarray(frames[v], dtype=float)
        hw = np.asarray(frames[w], dtype=float)
        out_links[(v, w)] = hv @ np.asarray(L, dtype=float) @ np.linalg.inv(hw)
        out_coframes[(v, w)] = hv @ np.asarray(coframes[(v, w)], dtype=float)
    return out_links, out_coframes


def minkowski_metric() -> Matrix:
    """Return diag(+1,-1,-1,-1)."""
    return np.diag([1.0, -1.0, -1.0, -1.0])


def assert_lorentz(matrix: Matrix, atol: float = 1e-10) -> None:
    """Raise AssertionError unless matrix preserves eta and is proper/orthochronous."""
    eta = minkowski_metric()
    M = np.asarray(matrix, dtype=float)
    assert np.allclose(M.T @ eta @ M, eta, atol=atol)
    assert np.linalg.det(M) > 0.0
    assert M[0, 0] >= 1.0 - atol

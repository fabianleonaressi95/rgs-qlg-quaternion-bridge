"""Minimal finite-complex Cartan core.

Convention:
- L[(v,w)] maps vectors in the frame at w to the frame at v.
- E[(v,w)] is the oriented edge coframe in the frame at v.
- Reverse edges are derived from L[(w,v)] = inv(L[(v,w)]) and
  E[(w,v)] = -inv(L[(v,w)]) @ E[(v,w)].

The tetrahedral boundary-word identity is implemented exactly in the affine
Poincare group. It is a finite combinatorial identity, not by itself a proof
of continuum convergence or of the differential-form Bianchi equations.
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


def oriented_link(edge: Edge, links: Mapping[Edge, Matrix]) -> Matrix:
    """Fetch an oriented link, deriving the reverse solely from L_wv=L_vw^-1."""
    v, w = edge
    if (v, w) in links:
        return np.asarray(links[(v, w)], dtype=float)
    if (w, v) in links:
        return np.linalg.inv(np.asarray(links[(w, v)], dtype=float))
    raise KeyError(f"missing link data for {edge} or its reverse")


def face_holonomy(face: Sequence[int], links: Mapping[Edge, Matrix]) -> Matrix:
    """Based face holonomy for an ordered triangular face (v0,v1,v2)."""
    if len(face) != 3:
        raise ValueError("face_holonomy expects an ordered triangle")
    v0, v1, v2 = map(int, face)
    L01 = oriented_link((v0, v1), links)
    L12 = oriented_link((v1, v2), links)
    L20 = oriented_link((v2, v0), links)
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



def affine_multiply(
    left: tuple[Matrix, Vector],
    right: tuple[Matrix, Vector],
) -> tuple[Matrix, Vector]:
    """Poincare semidirect-product law: (H,T)(K,S)=(HK,T+H S)."""
    H, T = left
    K, S = right
    return np.asarray(H) @ np.asarray(K), np.asarray(T) + np.asarray(H) @ np.asarray(S)


def affine_inverse(element: tuple[Matrix, Vector]) -> tuple[Matrix, Vector]:
    """Inverse in the affine Lorentz/Poincare semidirect product."""
    H, T = element
    Hinv = np.linalg.inv(np.asarray(H))
    return Hinv, -Hinv @ np.asarray(T)


def affine_edge_data(
    edge: Edge,
    links: Mapping[Edge, Matrix],
    coframes: Mapping[Edge, Vector],
) -> tuple[Matrix, Vector]:
    """Fetch one oriented affine edge, deriving its reverse consistently."""
    return edge_data(edge, links, coframes)


def affine_face_holonomy(
    face: Sequence[int],
    links: Mapping[Edge, Matrix],
    coframes: Mapping[Edge, Vector],
) -> tuple[Matrix, Vector]:
    """Return (Lorentz holonomy, transported coframe closure) for a triangle."""
    if len(face) != 3:
        raise ValueError("affine_face_holonomy expects an ordered triangle")
    v0, v1, v2 = map(int, face)
    g01 = affine_edge_data((v0, v1), links, coframes)
    g12 = affine_edge_data((v1, v2), links, coframes)
    g20 = affine_edge_data((v2, v0), links, coframes)
    return affine_multiply(affine_multiply(g01, g12), g20)


def tetrahedron_boundary_bianchi(
    links: Mapping[Edge, Matrix],
    coframes: Mapping[Edge, Vector],
    *,
    atol: float = 1e-10,
) -> dict[str, Matrix | Vector | float]:
    """Audit the exact based boundary-word identity on tetrahedron [0,1,2,3].

    Define A=Hol(012), B=Hol(023), C=Hol(031), and
    D=G_01 Hol(123) G_01^-1, all based at vertex 0. Edge cancellation gives
    A B C D^-1 = identity in SO+(1,3) semidirect R^(1,3).

    The returned curvature residual is H_A H_B H_C H_D^-1-I.
    The returned translational residual is
        T_A + H_A T_B + H_A H_B T_C - T_D.
    This is an exact finite boundary identity, not a continuum-limit theorem.
    """
    A = affine_face_holonomy((0, 1, 2), links, coframes)
    B = affine_face_holonomy((0, 2, 3), links, coframes)
    C = affine_face_holonomy((0, 3, 1), links, coframes)
    face_123 = affine_face_holonomy((1, 2, 3), links, coframes)
    g01 = affine_edge_data((0, 1), links, coframes)
    D = affine_multiply(affine_multiply(g01, face_123), affine_inverse(g01))

    ABC = affine_multiply(affine_multiply(A, B), C)
    boundary_product = affine_multiply(ABC, affine_inverse(D))
    H_residual = boundary_product[0] - np.eye(4)
    T_residual = boundary_product[1]

    # Equivalent expanded translational identity, exposing curvature transport.
    H_A, T_A = A
    H_B, T_B = B
    H_C, T_C = C
    _, T_D = D
    expanded_T = T_A + H_A @ T_B + H_A @ H_B @ T_C - T_D

    return {
        "curvature_residual": H_residual,
        "torsion_residual": T_residual,
        "expanded_torsion_residual": expanded_T,
        "curvature_norm": float(np.linalg.norm(H_residual)),
        "torsion_norm": float(np.linalg.norm(T_residual)),
        "expanded_torsion_norm": float(np.linalg.norm(expanded_T)),
        "passed": float(
            np.linalg.norm(H_residual) <= atol
            and np.linalg.norm(T_residual) <= atol
            and np.linalg.norm(expanded_T) <= atol
        ),
    }

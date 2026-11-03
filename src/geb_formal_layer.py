"""Formal Godel--Escher--Bach layer for the RGS-PSU-Q6 program.

Bach: finite substitution dynamics.
Escher: triangulated torus and sheet involution.
Godel: canonical finite encoding/hash, explicitly not an incompleteness proof.
The strange-loop fixed point is an abstract research equation; existence and
uniqueness require a separately proved contraction or another fixed-point theorem.
"""
from __future__ import annotations
from hashlib import sha256
import json
import numpy as np

PSI = 1.465571231876768


def bach_substitution_matrix() -> np.ndarray:
    """Incidence matrix of the recorded three-symbol substitution."""
    return np.array([[0, 0, 1], [1, 0, 0], [0, 1, 1]], dtype=np.int64)


def bach_spectral_audit() -> dict:
    M = bach_substitution_matrix()
    eig = np.linalg.eigvals(M.astype(float))
    rho = float(np.max(np.abs(eig)))
    return {
        "matrix": M.tolist(),
        "characteristic_polynomial_coefficients_descending": [1, -1, 0, -1],
        "supergolden_root": PSI,
        "spectral_radius": rho,
        "spectral_radius_residual": abs(rho - PSI),
        "matrix_polynomial_residual": float(np.linalg.norm(
            np.linalg.matrix_power(M, 3) - np.linalg.matrix_power(M, 2) - np.eye(3))),
    }


def _torus_complex(n: int = 4):
    """Periodic n-by-n square grid, each square split into two oriented triangles."""
    if n < 3:
        raise ValueError("n must be >= 3")
    vid = lambda i, j: (i % n) * n + (j % n)
    coords = {vid(i, j): (i, j) for i in range(n) for j in range(n)}
    edges = {}
    faces = []

    def edge_id(u, v):
        key = (min(u, v), max(u, v))
        if key not in edges:
            edges[key] = len(edges)
        return edges[key]

    for i in range(n):
        for j in range(n):
            a, b = vid(i, j), vid(i + 1, j)
            c, d = vid(i + 1, j + 1), vid(i, j + 1)
            faces.append((a, b, c))
            faces.append((a, c, d))
            edge_id(a, b); edge_id(b, c); edge_id(c, a)
            edge_id(a, c); edge_id(c, d); edge_id(d, a)

    edge_list = [None] * len(edges)
    for key, k in edges.items():
        edge_list[k] = key
    B1 = np.zeros((n*n, len(edge_list)), dtype=np.int64)
    for k, (u, v) in enumerate(edge_list):
        B1[u, k] = -1
        B1[v, k] = 1
    B2 = np.zeros((len(edge_list), len(faces)), dtype=np.int64)
    for f, (a, b, c) in enumerate(faces):
        for u, v in ((a, b), (b, c), (c, a)):
            k = edges[(min(u, v), max(u, v))]
            B2[k, f] += 1 if (u, v) == edge_list[k] else -1
    return coords, edge_list, faces, B1, B2


def torus_audit(n: int = 4) -> dict:
    coords, edges, faces, B1, B2 = _torus_complex(n)
    r1, r2 = np.linalg.matrix_rank(B1), np.linalg.matrix_rank(B2)
    chi = len(coords) - len(edges) + len(faces)
    alpha = np.zeros(len(edges))
    beta = np.zeros(len(edges))
    for k, (u, v) in enumerate(edges):
        i0, j0 = coords[u]; i1, j1 = coords[v]
        # Lifted increments in the periodic grid, oriented along stored edge.
        di = ((i1 - i0 + n//2) % n) - n//2
        dj = ((j1 - j0 + n//2) % n) - n//2
        alpha[k] = di / n
        beta[k] = dj / n
    # Fundamental cycles, with orientation matching increasing coordinates.
    edge_lookup = {tuple(sorted(edge)): k for k, edge in enumerate(edges)}
    cycle_x, cycle_y = [], []
    for j in range(n):
        row = np.zeros(len(edges))
        for i in range(n):
            u, v = i*n+j, ((i+1)%n)*n+j
            k = edge_lookup[(min(u,v), max(u,v))]
            row[k] = 1 if (u,v)==edges[k] else -1
        cycle_x.append(row)
    for i in range(n):
        col = np.zeros(len(edges))
        for j in range(n):
            u, v = i*n+j, i*n+(j+1)%n
            k = edges[(min(u,v), max(u,v))]
            col[k] = 1 if (u,v)==edges[k] else -1
        cycle_y.append(col)
    periods = np.array([
        [cycle_x[0] @ alpha, cycle_x[0] @ beta],
        [cycle_y[0] @ alpha, cycle_y[0] @ beta],
    ])
    H0 = B1 @ B1.T
    H1 = B1.T @ B1 + B2 @ B2.T
    H2 = B2.T @ B2
    eig0 = np.linalg.eigvalsh(H0.astype(float))
    eig1 = np.linalg.eigvalsh(H1.astype(float))
    eig2 = np.linalg.eigvalsh(H2.astype(float))
    return {
        "n": n, "V": len(coords), "E": len(edges), "F": len(faces),
        "euler_characteristic": int(chi),
        "betti_numbers_from_ranks": [len(coords)-r1, len(edges)-r1-r2, len(faces)-r2],
        "rank_B1": int(r1), "rank_B2": int(r2),
        "boundary_squared_residual": int(np.max(np.abs(B1 @ B2))),
        "period_matrix": periods.tolist(),
        "closed_cochain_residuals": [
            float(np.linalg.norm(B2.T @ alpha)),
            float(np.linalg.norm(B2.T @ beta))],
        "hodge_zero_mode_counts": [
            int(np.sum(np.abs(eig0)<1e-9)),
            int(np.sum(np.abs(eig1)<1e-9)),
            int(np.sum(np.abs(eig2)<1e-9))],
    }


def escher_sheet_involution(i: int, j: int, sheet: int) -> tuple[int, int, int]:
    """Deck involution J(i,j,s)=(i,j,1-s)."""
    if sheet not in (0, 1):
        raise ValueError("sheet must be 0 or 1")
    return int(i), int(j), 1-int(sheet)


def connected_double_cover_criterion(rho: tuple[int, int]) -> dict:
    """A Z2 cover of T2 is connected iff its monodromy class is nonzero."""
    if len(rho) != 2 or any(x not in (0, 1) for x in rho):
        raise ValueError("rho must be a pair in Z2")
    connected = any(rho)
    return {
        "monodromy_class": list(rho),
        "connected": connected,
        "cover_cells_if_connected": {"V": 32, "E": 96, "F": 64} if connected else
                                    {"V": 32, "E": 96, "F": 64},
        "euler_characteristic": 0,
        "criterion": "connected iff [rho] != 0 in H^1(T^2; Z2)",
    }


def godel_finite_encoding(payload: dict) -> dict:
    """Canonical finite serialization and SHA-256 identifier, not Gödel incompleteness."""
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return {"canonical_serialization": canonical,
            "sha256": sha256(canonical.encode("utf-8")).hexdigest(),
            "interpretation": "finite syntactic fingerprint only; no incompleteness theorem"}


def strange_loop_residual(x, recursion_map) -> float:
    """Evaluate ||R(x)-x||; existence/uniqueness are not inferred."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(recursion_map(x), dtype=float)
    if y.shape != x.shape:
        raise ValueError("recursion map must preserve the state-space shape")
    return float(np.linalg.norm(y-x))


def audit_geb() -> dict:
    torus = torus_audit(4)
    bach = bach_spectral_audit()
    involution_ok = all(escher_sheet_involution(*escher_sheet_involution(2,3,s)) == (2,3,s)
                        for s in (0,1))
    cover = connected_double_cover_criterion((1,0))
    formal_payload = {
        "torus": {"V": torus["V"], "E": torus["E"], "F": torus["F"],
                  "betti": torus["betti_numbers_from_ranks"]},
        "bach_substitution": bach["matrix"],
        "escher_involution": "J(i,j,s)=(i,j,1-s)",
        "godel_layer": "canonical finite encoding",
    }
    return {
        "audit": "RGS_PSU_GODEL_ESCHER_BACH",
        "status": "finite combinatorial and algebraic identities; physical strange-loop interpretation open",
        "bach": bach,
        "escher": {**torus, "J_squared_identity": involution_ok,
                   "double_cover": cover},
        "godel": godel_finite_encoding(formal_payload),
        "strange_loop": {
            "fixed_point_equation": "G_* = R(G_*, Omega_*, Gamma_*)",
            "status": "research ansatz; fixed-point theorem requires hypotheses on R",
        },
        "open_obligations": [
            "define the exact self-reference language and provability predicate",
            "do not infer Gödel incompleteness from a hash or finite encoding",
            "derive the Bach substitution from the musical/recursive construction",
            "construct a nontrivial Z2 cocycle explicitly on the triangulated torus",
            "build a period-preserving intertwiner between refinement and Q6/PRISM",
            "derive the strange-loop map R and prove contraction or an alternative fixed-point theorem",
            "connect the formal system to the PSU action through an independently testable observable",
        ],
    }

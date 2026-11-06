"""Curved RGS refinement family from the recovered 64B/O_h seed.

The smooth target is the spatial Ellis-RGS metric
    ds^2 = du^2 + (u^2+r0^2) dOmega^2.
It has an exact isometric hypersurface embedding in R^4:
    F(u,n) = (sqrt(u^2+r0^2) n, r0 asinh(u/r0)).
The mesh is a radial product of a spherical triangulation built from the
actual 64-vertex RGS directions. Its P1 metric is the induced chordal
metric of each 3-simplex in R^4.

This is a finite-domain geometric refinement benchmark, not a claim that
the archived diagnostic profile solves the full RGS-EMS field equations.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, product
import hashlib
from typing import Any

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix, diags
from scipy.sparse.linalg import eigsh
from scipy.spatial import ConvexHull


# Exact absolute coordinate values recovered from
# A14_RGS_PSU_QLG_JACOBIAN_BUNDLE.zip:
# SHA256 80442238349cd6ddc720dabae93d083075dea88decba8acf0a78c3ec03e0abcc
_ORBIT_VALUES = {
    "cube_inner": 0.47863916086470326,
    "small_1": 0.48042984818593559,
    "small_2": 0.48230704367200045,
    "cube_outer": 0.95183963342906408,
    "large_1": 1.5073116806169751,
    "large_2": 1.9039106736116611,
}
SOURCE_ARCHIVES = {
    "A14": {
        "name": "A14_RGS_PSU_QLG_JACOBIAN_BUNDLE.zip",
        "sha256": "80442238349cd6ddc720dabae93d083075dea88decba8acf0a78c3ec03e0abcc",
        "array": "RGS_PSU_QlG_sources/RGS_PSU_all_matrices.npz::X_Oh_reconstructed",
    },
    "A13": {
        "name": "A13_RGS_PSU_QLG_KAGGLE_BUNDLE.zip",
        "sha256": "7ca0cdfeadbe1d693d4f776a5e22d957630b206b842f4b822f3afb578c488922",
        "array": "RGS_PSU_SEA_DERIVATION/sea_operators.npz::edges,edge_coframe",
    },
}
SOURCE_ARRAY_SHA256 = {
    "X_Oh_reconstructed_float64": "6a158cb30a671a5c994d8ab223e9283b57ac703808041d46e67f286126e419cf",
    "Q6_edges_int64": "da7854ca655577b629dc508208918326fde1d50563485f953ea0fdc1b2fb78e8",
    "edge_coframe_float64": "72980c38fdc652ece31895cf5c9f113a6206675fdf876c924bf6890d6e7fc263",
}


def _array_sha256(array: np.ndarray, dtype: str) -> str:
    canonical = np.ascontiguousarray(np.asarray(array, dtype=dtype))
    return hashlib.sha256(canonical.tobytes(order="C")).hexdigest()


@dataclass(frozen=True)
class RGSSeedData:
    coordinates: np.ndarray
    q6_edges: np.ndarray
    q6_raw_coframe: np.ndarray
    q6_unit_coframe: np.ndarray
    directions: np.ndarray
    vertex_to_direction: np.ndarray


@dataclass(frozen=True)
class EllisRGSRefinement:
    level: int
    r0: float
    half_width: float
    radial_segments: int
    radial_coordinates: np.ndarray
    angular_directions: np.ndarray
    angular_faces: np.ndarray
    vertices_r4: np.ndarray
    tetrahedra: np.ndarray
    tetrahedron_parameter_vertices: np.ndarray
    tetrahedron_angular_faces: np.ndarray
    tetrahedron_radial_layers: np.ndarray
    element_gram: np.ndarray
    element_volumes: np.ndarray
    element_condition: np.ndarray
    max_element_diameter: float
    chain: dict[str, Any] | None


def load_rgs_seed_data() -> RGSSeedData:
    """Reconstruct the archived 64B/O_h coordinate array exactly."""
    v = _ORBIT_VALUES
    orbit_representatives = [
        (v["cube_inner"], v["cube_inner"], v["cube_inner"]),
        (v["large_1"], v["large_1"], v["small_1"]),
        (v["large_1"], v["small_1"], v["large_1"]),
        (v["small_1"], v["large_1"], v["large_1"]),
        (v["large_2"], v["small_2"], v["small_2"]),
        (v["small_2"], v["large_2"], v["small_2"]),
        (v["small_2"], v["small_2"], v["large_2"]),
        (v["cube_outer"], v["cube_outer"], v["cube_outer"]),
    ]
    signs = tuple(product((-1.0, 1.0), repeat=3))
    X = np.asarray(
        [[s[0] * p[0], s[1] * p[1], s[2] * p[2]]
         for p in orbit_representatives for s in signs],
        dtype=float,
    )
    if X.shape != (64, 3):
        raise RuntimeError("reconstructed RGS seed has wrong shape")

    edges: list[tuple[int, int, int]] = []
    for g in range(64):
        for axis in range(6):
            mask = 1 << axis
            if (g & mask) == 0:
                edges.append((g, g | mask, axis))
    E = np.asarray(edges, dtype=np.int16)
    if E.shape != (192, 3):
        raise RuntimeError("reconstructed Q6 edge list has wrong shape")

    radii = np.linalg.norm(X, axis=1)
    raw = np.column_stack((
        X[E[:, 1]] - X[E[:, 0]],
        radii[E[:, 1]] - radii[E[:, 0]],
    ))
    lengths = np.linalg.norm(raw, axis=1)
    if np.any(lengths <= 0.0):
        raise RuntimeError("zero-length Q6 seed edge")
    unit = raw / lengths[:, None]
    hashes = {
        "X_Oh_reconstructed_float64": _array_sha256(X, "<f8"),
        "Q6_edges_int64": _array_sha256(E, "<i8"),
        "edge_coframe_float64": _array_sha256(unit, "<f8"),
    }
    for key, expected in SOURCE_ARRAY_SHA256.items():
        if hashes[key] != expected:
            raise RuntimeError(f"recovered RGS seed hash mismatch for {key}")

    # Deduplicate only coincident rays; keep the inverse map for all 64
    # original labels. This is an exact deterministic quotient of directions.
    directions, inverse = np.unique(
        X / radii[:, None], axis=0, return_inverse=True
    )
    if directions.shape != (56, 3):
        raise RuntimeError(
            f"expected 56 distinct RGS directions, got {directions.shape[0]}"
        )
    return RGSSeedData(
        coordinates=X,
        q6_edges=E,
        q6_raw_coframe=raw,
        q6_unit_coframe=unit,
        directions=directions,
        vertex_to_direction=inverse.astype(np.int64),
    )


def _orient_spherical_faces(
    directions: np.ndarray, faces: np.ndarray
) -> np.ndarray:
    out = np.asarray(faces, dtype=np.int64).copy()
    for k, tri in enumerate(out):
        a, b, c = directions[tri]
        if float(np.dot(np.cross(b - a, c - a), a + b + c)) < 0.0:
            out[k, 1], out[k, 2] = out[k, 2], out[k, 1]
    return out


def _initial_spherical_mesh(
    directions: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    hull = ConvexHull(directions)
    if hull.vertices.size != directions.shape[0]:
        raise RuntimeError("some distinct spherical seed directions are not hull vertices")
    faces = _orient_spherical_faces(directions, hull.simplices)
    if faces.shape != (108, 3):
        raise RuntimeError(f"expected 108 seed sphere faces, got {faces.shape[0]}")
    return directions.copy(), faces


def _refine_sphere(
    directions: np.ndarray, faces: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    vertices = [x.copy() for x in directions]
    midpoint_ids: dict[tuple[int, int], int] = {}

    def midpoint(i: int, j: int) -> int:
        key = (min(i, j), max(i, j))
        if key not in midpoint_ids:
            x = vertices[i] + vertices[j]
            norm = float(np.linalg.norm(x))
            if norm <= 1e-14:
                raise RuntimeError("antipodal edge encountered in spherical refinement")
            vertices.append(x / norm)
            midpoint_ids[key] = len(vertices) - 1
        return midpoint_ids[key]

    refined: list[tuple[int, int, int]] = []
    for a, b, c in faces:
        ab = midpoint(int(a), int(b))
        bc = midpoint(int(b), int(c))
        ca = midpoint(int(c), int(a))
        refined.extend((
            (int(a), ab, ca),
            (ab, int(b), bc),
            (ca, bc, int(c)),
            (ab, bc, ca),
        ))
    V = np.asarray(vertices, dtype=float)
    F = _orient_spherical_faces(V, np.asarray(refined, dtype=np.int64))
    return V, F


def ellis_embedding(
    u: np.ndarray, directions: np.ndarray, r0: float
) -> np.ndarray:
    """Exact R4 hypersurface embedding of du^2 + (u^2+r0^2)dOmega^2."""
    uu = np.asarray(u, dtype=float).reshape(-1)
    nn = np.asarray(directions, dtype=float)
    R = np.sqrt(uu * uu + r0 * r0)
    z = r0 * np.arcsinh(uu / r0)
    return np.concatenate(
        [np.column_stack((R[j] * nn, np.full(nn.shape[0], z[j])))
         for j in range(uu.size)],
        axis=0,
    )


def _product_tetrahedra(
    n_radial: int, n_angular: int, faces: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Conforming staircase prism split plus exact (t,alpha,beta) charts."""
    tets: list[tuple[int, int, int, int]] = []
    params: list[np.ndarray] = []
    angular_faces: list[tuple[int, int, int]] = []
    radial_layers: list[int] = []
    for j in range(n_radial):
        lo = j * n_angular
        hi = (j + 1) * n_angular
        for tri in faces:
            a, b, c = sorted(map(int, tri))
            A, B, C = lo + a, lo + b, lo + c
            Ap, Bp, Cp = hi + a, hi + b, hi + c
            # Parameter coordinates are (t,alpha,beta), where the
            # spherical face is n=normalize((1-alpha-beta)n_a
            # + alpha*n_b + beta*n_c).
            local = (
                ((A, B, C, Cp), np.asarray([
                    (0., 0., 0.), (0., 1., 0.),
                    (0., 0., 1.), (1., 0., 1.)
                ])),
                ((A, B, Bp, Cp), np.asarray([
                    (0., 0., 0.), (0., 1., 0.),
                    (1., 1., 0.), (1., 0., 1.)
                ])),
                ((A, Ap, Bp, Cp), np.asarray([
                    (0., 0., 0.), (1., 0., 0.),
                    (1., 1., 0.), (1., 0., 1.)
                ])),
            )
            for ids, coordinates in local:
                ordered = sorted(
                    zip(ids, coordinates), key=lambda item: item[0]
                )
                tets.append(tuple(int(item[0]) for item in ordered))
                params.append(np.asarray([item[1] for item in ordered]))
                angular_faces.append((a, b, c))
                radial_layers.append(j)
    return (
        np.asarray(tets, dtype=np.int64),
        np.asarray(params, dtype=float),
        np.asarray(angular_faces, dtype=np.int64),
        np.asarray(radial_layers, dtype=np.int64),
    )

def _element_geometry(
    vertices: np.ndarray, tetrahedra: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    points = vertices[tetrahedra]
    J = np.transpose(points[:, 1:, :] - points[:, :1, :], (0, 2, 1))
    gram = np.einsum("tki,tkj->tij", J, J)
    eigenvalues = np.linalg.eigvalsh(gram)
    if np.any(eigenvalues[:, 0] <= 1e-13):
        bad = int(np.count_nonzero(eigenvalues[:, 0] <= 1e-13))
        raise RuntimeError(f"{bad} degenerate or ill-conditioned curved tetrahedra")
    volumes = np.sqrt(np.linalg.det(gram)) / 6.0
    conditions = eigenvalues[:, -1] / eigenvalues[:, 0]
    edge_pairs = np.asarray(list(combinations(range(4), 2)), dtype=np.int64)
    edge_vectors = (
        points[:, edge_pairs[:, 0], :] - points[:, edge_pairs[:, 1], :]
    )
    diameters = np.linalg.norm(edge_vectors, axis=2)
    return gram, volumes, conditions, float(np.max(diameters))


def _build_chain_complex(
    n_vertices: int, tetrahedra: np.ndarray
) -> dict[str, Any]:
    """Build integer simplicial boundaries B1, B2, B3 and audit d^2=0."""
    edge_set: set[tuple[int, int]] = set()
    face_set: set[tuple[int, int, int]] = set()
    for tet in tetrahedra:
        for edge in combinations(map(int, tet), 2):
            edge_set.add(tuple(sorted(edge)))
        for face in combinations(map(int, tet), 3):
            face_set.add(tuple(sorted(face)))
    edges = sorted(edge_set)
    faces = sorted(face_set)
    edge_id = {e: k for k, e in enumerate(edges)}
    face_id = {f: k for k, f in enumerate(faces)}

    r1: list[int] = []
    c1: list[int] = []
    d1: list[int] = []
    for k, (a, b) in enumerate(edges):
        r1.extend((a, b)); c1.extend((k, k)); d1.extend((-1, 1))
    B1 = coo_matrix(
        (np.asarray(d1, dtype=np.int8), (r1, c1)),
        shape=(n_vertices, len(edges)), dtype=np.int8,
    ).tocsr()

    r2: list[int] = []
    c2: list[int] = []
    d2: list[int] = []
    for k, (a, b, c) in enumerate(faces):
        for edge, sign in (((b, c), 1), ((a, c), -1), ((a, b), 1)):
            r2.append(edge_id[edge]); c2.append(k); d2.append(sign)
    B2 = coo_matrix(
        (np.asarray(d2, dtype=np.int8), (r2, c2)),
        shape=(len(edges), len(faces)), dtype=np.int8,
    ).tocsr()

    r3: list[int] = []
    c3: list[int] = []
    d3: list[int] = []
    for k, (a, b, c, d) in enumerate(tetrahedra):
        boundary = (
            ((b, c, d), 1),
            ((a, c, d), -1),
            ((a, b, d), 1),
            ((a, b, c), -1),
        )
        for face, sign in boundary:
            r3.append(face_id[tuple(sorted(face))])
            c3.append(k)
            d3.append(sign)
    B3 = coo_matrix(
        (np.asarray(d3, dtype=np.int8), (r3, c3)),
        shape=(len(faces), len(tetrahedra)), dtype=np.int8,
    ).tocsr()

    chain12 = (B1.astype(np.int64) @ B2.astype(np.int64)).tocsr()
    chain23 = (B2.astype(np.int64) @ B3.astype(np.int64)).tocsr()
    return {
        "edges": edges,
        "faces": faces,
        "B1": B1,
        "B2": B2,
        "B3": B3,
        "B1B2_nnz": int(chain12.nnz),
        "B2B3_nnz": int(chain23.nnz),
        "chain_complex_exact": bool(chain12.nnz == 0 and chain23.nnz == 0),
    }


def build_ellis_rgs_refinement(
    level: int,
    r0: float = 0.8292,
    half_width: float | None = None,
    radial_base: int = 8,
    include_chain: bool = False,
) -> EllisRGSRefinement:
    """Build one level of the curved RGS Ellis-wormhole product family.

    At level l, angular edges are midpoint-refined l times and the radial
    interval has radial_base*2**l segments. The finite domain is [-L,L]xS2.
    """
    lev = int(level)
    if lev < 0 or lev > 3:
        raise ValueError("level must be between 0 and 3")
    r0 = float(r0)
    if not np.isfinite(r0) or r0 <= 0.0:
        raise ValueError("r0 must be finite and positive")
    L = 3.0 * r0 if half_width is None else float(half_width)
    if not np.isfinite(L) or L <= 0.0:
        raise ValueError("half_width must be finite and positive")
    rb = int(radial_base)
    if rb < 2:
        raise ValueError("radial_base must be at least 2")

    seed = load_rgs_seed_data()
    directions, faces = _initial_spherical_mesh(seed.directions)
    for _ in range(lev):
        directions, faces = _refine_sphere(directions, faces)

    n_radial = rb * (2**lev)
    u = np.linspace(-L, L, n_radial + 1)
    vertices = ellis_embedding(u, directions, r0)
    (tetrahedra, parameter_vertices, tetrahedron_faces,
     tetrahedron_layers) = _product_tetrahedra(
        n_radial, len(directions), faces
    )
    gram, volumes, conditions, max_diameter = _element_geometry(
        vertices, tetrahedra
    )
    chain = (
        _build_chain_complex(len(vertices), tetrahedra)
        if include_chain else None
    )
    return EllisRGSRefinement(
        level=lev,
        r0=r0,
        half_width=L,
        radial_segments=n_radial,
        radial_coordinates=u,
        angular_directions=directions,
        angular_faces=faces,
        vertices_r4=vertices,
        tetrahedra=tetrahedra,
        tetrahedron_parameter_vertices=parameter_vertices,
        tetrahedron_angular_faces=tetrahedron_faces,
        tetrahedron_radial_layers=tetrahedron_layers,
        element_gram=gram,
        element_volumes=volumes,
        element_condition=conditions,
        max_element_diameter=max_diameter,
        chain=chain,
    )


def assemble_ellis_p1_operators(
    mesh: EllisRGSRefinement,
) -> tuple[csr_matrix, csr_matrix]:
    """Assemble exact local P1 L2 mass and Dirichlet stiffness for chordal tets."""
    D = np.asarray(
        [[-1.0, -1.0, -1.0], [1.0, 0.0, 0.0],
         [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        dtype=float,
    )
    inverse_gram = np.linalg.inv(mesh.element_gram)
    local_K = mesh.element_volumes[:, None, None] * np.einsum(
        "ai,tij,bj->tab", D, inverse_gram, D
    )
    local_M = (
        mesh.element_volumes[:, None, None] / 20.0
        * (np.eye(4)[None, :, :] + np.ones((1, 4, 4)))
    )
    n = mesh.vertices_r4.shape[0]
    rows = np.repeat(mesh.tetrahedra, 4, axis=1).ravel()
    cols = np.tile(mesh.tetrahedra, (1, 4)).ravel()
    M = coo_matrix((local_M.ravel(), (rows, cols)), shape=(n, n)).tocsr()
    K = coo_matrix((local_K.ravel(), (rows, cols)), shape=(n, n)).tocsr()
    return M, K


def _duffy_tetra_quadrature(order: int) -> tuple[np.ndarray, np.ndarray]:
    """Tensor Gauss-Legendre quadrature mapped to the reference 3-simplex."""
    n = int(order)
    if n < 2 or n > 10:
        raise ValueError("quadrature_order must be between 2 and 10")
    nodes, weights = np.polynomial.legendre.leggauss(n)
    nodes = 0.5 * (nodes + 1.0)
    weights = 0.5 * weights
    points: list[tuple[float, float, float]] = []
    out_weights: list[float] = []
    for i, r in enumerate(nodes):
        for j, s in enumerate(nodes):
            for k, t in enumerate(nodes):
                points.append((r, (1.0 - r) * s,
                               (1.0 - r) * (1.0 - s) * t))
                out_weights.append(
                    weights[i] * weights[j] * weights[k]
                    * (1.0 - r) ** 2 * (1.0 - s)
                )
    return np.asarray(points), np.asarray(out_weights)


@dataclass(frozen=True)
class EllisP1Embedding:
    """P1 nodal data pulled back to the common L2 Ellis manifold."""

    mesh: EllisRGSRefinement
    nodal_values: np.ndarray

    def __post_init__(self) -> None:
        values = np.asarray(self.nodal_values)
        if values.ndim != 1 or values.size != self.mesh.vertices_r4.shape[0]:
            raise ValueError("nodal_values must have one value per curved-mesh vertex")
        if not np.all(np.isfinite(values)):
            raise ValueError("nodal_values must be finite")
        object.__setattr__(self, "nodal_values", values.copy())

    def evaluate_on_element(
        self, element_index: int, barycentric: np.ndarray
    ) -> np.ndarray | float:
        """Evaluate J_h u_h at barycentric points of one curved element."""
        k = int(element_index)
        if k < 0 or k >= self.mesh.tetrahedra.shape[0]:
            raise IndexError("element_index out of range")
        lam = np.asarray(barycentric, dtype=float)
        scalar = lam.ndim == 1
        if scalar:
            lam = lam[None, :]
        if lam.ndim != 2 or lam.shape[1] != 4:
            raise ValueError("barycentric must have shape (n_points,4)")
        if not np.all(np.isfinite(lam)):
            raise ValueError("barycentric coordinates must be finite")
        if np.max(np.abs(np.sum(lam, axis=1) - 1.0)) > 1e-10:
            raise ValueError("barycentric coordinates must sum to one")
        if np.min(lam) < -1e-10:
            raise ValueError("barycentric coordinates must be nonnegative")
        values = lam @ self.nodal_values[self.mesh.tetrahedra[k]]
        return float(values[0]) if scalar else values


def J_h(mesh: EllisRGSRefinement, nodal_values: np.ndarray) -> EllisP1Embedding:
    """Embed curved-mesh nodal data into H=L2(M_L,dmu_g) by P1 pullback."""
    return EllisP1Embedding(mesh=mesh, nodal_values=nodal_values)


def assemble_ellis_common_l2_operators(
    mesh: EllisRGSRefinement,
    quadrature_order: int = 6,
    chunk_size: int = 128,
) -> tuple[csr_matrix, csr_matrix]:
    """Integrate P1 mass and Dirichlet stiffness in the exact Ellis metric.

    Each reference tetrahedron maps into the smooth manifold by
    (t,alpha,beta) -> (u(t), normalize((1-alpha-beta)n_a+
    alpha*n_b+beta*n_c)) -> F(u,n). The integration is over the true
    metric du^2+(u^2+r0^2)dOmega^2, not the chordal PL surrogate.
    """
    qpoints, qweights = _duffy_tetra_quadrature(quadrature_order)
    D = np.asarray(
        [[-1.0, -1.0, -1.0], [1.0, 0.0, 0.0],
         [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        dtype=float,
    )
    barycentric = np.column_stack((
        1.0 - np.sum(qpoints, axis=1), qpoints
    ))
    n_tets = mesh.tetrahedra.shape[0]
    local_mass = np.empty((n_tets, 4, 4), dtype=float)
    local_stiffness = np.empty_like(local_mass)
    directions = mesh.angular_directions

    for first in range(0, n_tets, int(chunk_size)):
        last = min(first + int(chunk_size), n_tets)
        params = mesh.tetrahedron_parameter_vertices[first:last]
        affine = np.transpose(
            params[:, 1:, :] - params[:, :1, :], (0, 2, 1)
        )
        y = params[:, None, 0, :] + np.einsum(
            "tij,qj->tqi", affine, qpoints
        )
        tcoord, alpha, beta = y[..., 0], y[..., 1], y[..., 2]

        face = mesh.tetrahedron_angular_faces[first:last]
        na = directions[face[:, 0]]
        nb = directions[face[:, 1]]
        nc = directions[face[:, 2]]
        vab = nb - na
        vac = nc - na
        w = (
            na[:, None, :]
            + alpha[..., None] * vab[:, None, :]
            + beta[..., None] * vac[:, None, :]
        )
        wnorm = np.linalg.norm(w, axis=2)
        if np.any(wnorm <= 1e-12):
            raise RuntimeError("singular spherical chart in curved P1 element")
        nvec = w / wnorm[..., None]
        ndot_ab = np.einsum("tqi,ti->tq", nvec, vab)
        ndot_ac = np.einsum("tqi,ti->tq", nvec, vac)
        dn_alpha = (
            vab[:, None, :] - nvec * ndot_ab[..., None]
        ) / wnorm[..., None]
        dn_beta = (
            vac[:, None, :] - nvec * ndot_ac[..., None]
        ) / wnorm[..., None]

        layers = mesh.tetrahedron_radial_layers[first:last]
        u0 = mesh.radial_coordinates[layers]
        du = (
            mesh.radial_coordinates[layers + 1]
            - mesh.radial_coordinates[layers]
        )
        u = u0[:, None] + du[:, None] * tcoord
        R = np.sqrt(u * u + mesh.r0 * mesh.r0)
        Rp = u / R
        zp = mesh.r0 / R

        dF_dt = np.concatenate(
            (Rp[..., None] * nvec * du[:, None, None],
             (zp * du[:, None])[:, :, None]),
            axis=2,
        )
        dF_da = np.concatenate(
            (R[..., None] * dn_alpha, np.zeros((*R.shape, 1))), axis=2
        )
        dF_db = np.concatenate(
            (R[..., None] * dn_beta, np.zeros((*R.shape, 1))), axis=2
        )
        dF_dy = np.stack((dF_dt, dF_da, dF_db), axis=-1)
        dF_dxi = np.einsum("tqai,tij->tqaj", dF_dy, affine)
        gram = np.einsum("tqki,tqkj->tqij", dF_dxi, dF_dxi)
        determinant = np.linalg.det(gram)
        if np.any(determinant <= 1e-18):
            raise RuntimeError("nonpositive curved element metric determinant")
        inverse_gram = np.linalg.inv(gram)
        weighted_jacobian = np.sqrt(determinant) * qweights[None, :]
        local_mass[first:last] = np.einsum(
            "tq,qi,qj->tij",
            weighted_jacobian, barycentric, barycentric,
        )
        local_stiffness[first:last] = np.einsum(
            "tq,ai,tqij,bj->tab",
            weighted_jacobian, D, inverse_gram, D,
        )

    n_vertices = mesh.vertices_r4.shape[0]
    rows = np.repeat(mesh.tetrahedra, 4, axis=1).ravel()
    cols = np.tile(mesh.tetrahedra, (1, 4)).ravel()
    mass = coo_matrix(
        (local_mass.ravel(), (rows, cols)),
        shape=(n_vertices, n_vertices),
    ).tocsr()
    stiffness = coo_matrix(
        (local_stiffness.ravel(), (rows, cols)),
        shape=(n_vertices, n_vertices),
    ).tocsr()
    return mass, stiffness


def audit_ellis_lumped_mass_equivalence(
    mesh: EllisRGSRefinement, quadrature_order: int = 4
) -> dict[str, Any]:
    """Compare exact curved P1 mass with the chordal lumped vertex measure.

    The spectrum of D_lump^{-1/2} M_curved D_lump^{-1/2} gives finite-mesh
    equivalence constants without forming a dense global mass matrix.
    """
    mass, _ = assemble_ellis_common_l2_operators(
        mesh, quadrature_order=quadrature_order
    )
    lumped = np.zeros(mesh.vertices_r4.shape[0], dtype=float)
    np.add.at(
        lumped,
        mesh.tetrahedra.ravel(),
        np.repeat(mesh.element_volumes / 4.0, 4),
    )
    if np.any(lumped <= 0.0) or np.any(~np.isfinite(lumped)):
        raise RuntimeError("invalid chordal lumped vertex measure")
    invsqrt = 1.0 / np.sqrt(lumped)
    scaled = diags(invsqrt) @ mass @ diags(invsqrt)
    scaled = ((scaled + scaled.T) * 0.5).tocsr()
    if scaled.shape[0] <= 2:
        eigenvalues = np.linalg.eigvalsh(scaled.toarray())
        lower, upper = float(eigenvalues[0]), float(eigenvalues[-1])
    else:
        lower = float(eigsh(scaled, k=1, which="SA", return_eigenvectors=False)[0])
        upper = float(eigsh(scaled, k=1, which="LA", return_eigenvectors=False)[0])
    exact_volume = 4.0 * np.pi * (
        2.0 * mesh.half_width**3 / 3.0 + 2.0 * mesh.r0**2 * mesh.half_width
    )
    return {
        "gate": "CURVED_ELLIS_LUMPED_VS_COMMON_L2_MASS",
        "vertices": int(lumped.size),
        "tetrahedra": int(mesh.tetrahedra.shape[0]),
        "quadrature_order_per_duffy_axis": int(quadrature_order),
        "generalized_mass_eigenvalue_min": lower,
        "generalized_mass_eigenvalue_max": upper,
        "lumped_measure_total": float(np.sum(lumped)),
        "curved_mass_total": float(np.sum(mass)),
        "exact_ellis_volume": float(exact_volume),
        "curved_mass_volume_relative_error": float(
            abs(float(np.sum(mass)) - exact_volume) / exact_volume
        ),
        "lumped_measure_volume_relative_error": float(
            abs(float(np.sum(lumped)) - exact_volume) / exact_volume
        ),
        "positive_mass": bool(lower > 0.0),
        "scope": "finite-level spectral equivalence audit; uniform-in-level constants remain to be proved",
    }


def audit_ellis_g4_common_l2(
    mesh: EllisRGSRefinement, quadrature_order: int = 6
) -> dict[str, Any]:
    """Audit the common-L2 P1 embedding and true curved metric operators."""
    mass, stiffness = assemble_ellis_common_l2_operators(
        mesh, quadrature_order=quadrature_order
    )
    one = np.ones(mesh.vertices_r4.shape[0])
    exact_volume = 4.0 * np.pi * (
        2.0 * mesh.half_width**3 / 3.0 + 2.0 * mesh.r0**2 * mesh.half_width
    )
    volume = float(one @ (mass @ one))
    mass_dense = mass.toarray()
    mass_eigen_min = float(np.linalg.eigvalsh(
        (mass_dense + mass_dense.T) / 2.0
    )[0])
    return {
        "gate": "CURVED_ELLIS_G4_COMMON_L2_EMBEDDING",
        "hilbert_space": "L2([-L,L] x S2, R(u)^2 du dOmega)",
        "embedding": "continuous P1 on curved reference tetrahedra, pulled back by exact Ellis element maps",
        "mass_symmetry_residual": float(np.linalg.norm((mass - mass.T).data)),
        "stiffness_symmetry_residual": float(np.linalg.norm((stiffness - stiffness.T).data)),
        "mass_min_eigenvalue": mass_eigen_min,
        "mass_positive_definite": bool(mass_eigen_min > 0.0),
        "stiffness_constant_mode_residual": float(np.linalg.norm(stiffness @ one)),
        "volume_from_mass": volume,
        "volume_exact": exact_volume,
        "volume_relative_error": abs(volume - exact_volume) / exact_volume,
        "quadrature_order_per_duffy_axis": int(quadrature_order),
        "scope": "finite Ellis-RGS domain with boundary; common L2 embedding is defined for this curved benchmark",
    }


def audit_ellis_rgs_refinement(
    levels: tuple[int, ...] = (0, 1, 2),
    r0: float = 0.8292,
    half_width: float | None = None,
    radial_base: int = 8,
) -> dict[str, Any]:
    """Multi-level geometric, coframe, volume and chain-complex audit."""
    if len(levels) < 3 or any(b <= a for a, b in zip(levels, levels[1:])):
        raise ValueError("levels must contain at least three strictly increasing integers")
    seed = load_rgs_seed_data()
    L = 3.0 * r0 if half_width is None else float(half_width)
    reports: list[dict[str, Any]] = []
    for level in levels:
        mesh = build_ellis_rgs_refinement(
            level, r0=r0, half_width=L, radial_base=radial_base,
            include_chain=(int(level) == int(levels[0])),
        )
        exact_volume = 4.0 * np.pi * (
            2.0 * L**3 / 3.0 + 2.0 * r0**2 * L
        )
        volume = float(np.sum(mesh.element_volumes))
        Rprime = mesh.radial_coordinates / np.sqrt(
            mesh.radial_coordinates**2 + r0**2
        )
        zprime = r0 / np.sqrt(mesh.radial_coordinates**2 + r0**2)
        angular_edges = sorted({
            tuple(sorted((int(a), int(b))))
            for tri in mesh.angular_faces
            for a, b in combinations(map(int, tri), 2)
        })
        angular_edge_lengths = np.asarray([
            np.arccos(np.clip(
                np.dot(mesh.angular_directions[a], mesh.angular_directions[b]),
                -1.0, 1.0
            ))
            for a, b in angular_edges
        ])
        angular_face_areas = np.asarray([
            2.0 * np.arctan2(
                abs(float(np.dot(
                    mesh.angular_directions[int(a)],
                    np.cross(mesh.angular_directions[int(b)],
                             mesh.angular_directions[int(c)])
                ))),
                1.0 + float(np.dot(mesh.angular_directions[int(a)], mesh.angular_directions[int(b)]))
                  + float(np.dot(mesh.angular_directions[int(b)], mesh.angular_directions[int(c)]))
                  + float(np.dot(mesh.angular_directions[int(c)], mesh.angular_directions[int(a)]))
            )
            for a, b, c in mesh.angular_faces
        ])
        reports.append({
            "level": mesh.level,
            "angular_edge_min": float(np.min(angular_edge_lengths)),
            "angular_edge_max": float(np.max(angular_edge_lengths)),
            "angular_face_area_min": float(np.min(angular_face_areas)),
            "angular_face_area_max": float(np.max(angular_face_areas)),
            "angular_shape_proxy_min_area_over_h2": float(
                np.min(angular_face_areas) / np.max(angular_edge_lengths)**2
            ),
            "vertices": int(mesh.vertices_r4.shape[0]),
            "angular_vertices": int(mesh.angular_directions.shape[0]),
            "angular_faces": int(mesh.angular_faces.shape[0]),
            "radial_segments": mesh.radial_segments,
            "tetrahedra": int(mesh.tetrahedra.shape[0]),
            "h_max_chord": mesh.max_element_diameter,
            "element_condition_max": float(np.max(mesh.element_condition)),
            "element_condition_min": float(np.min(mesh.element_condition)),
            "element_gram_min_eigenvalue": float(
                np.min(np.linalg.eigvalsh(mesh.element_gram)[:, 0])
            ),
            "volume_pl": volume,
            "volume_exact_ellis": exact_volume,
            "volume_relative_error": abs(volume - exact_volume) / exact_volume,
            "embedding_identity_max_residual": float(
                np.max(np.abs(Rprime**2 + zprime**2 - 1.0))
            ),
            "chain_complex_exact": (
                mesh.chain["chain_complex_exact"] if mesh.chain is not None else None
            ),
            "B1B2_nnz": mesh.chain["B1B2_nnz"] if mesh.chain is not None else None,
            "B2B3_nnz": mesh.chain["B2B3_nnz"] if mesh.chain is not None else None,
        })

    radii = np.linalg.norm(seed.coordinates, axis=1)
    raw_reconstructed = np.column_stack((
        seed.coordinates[seed.q6_edges[:, 1]]
        - seed.coordinates[seed.q6_edges[:, 0]],
        radii[seed.q6_edges[:, 1]] - radii[seed.q6_edges[:, 0]],
    ))
    coframe_residual = float(
        np.linalg.norm(raw_reconstructed - seed.q6_raw_coframe)
    )
    unit_residual = float(np.max(np.abs(
        np.linalg.norm(seed.q6_unit_coframe, axis=1) - 1.0
    )))
    return {
        "audit": "RGS_ELLIS_CURVED_REFINEMENT_G1",
        "status": "GEOMETRIC_FAMILY_CONSTRUCTED_NUMERICALLY_AUDITED",
        "target_metric": "du^2 + (u^2+r0^2) dOmega^2",
        "exact_embedding": "F(u,n)=(sqrt(u^2+r0^2)*n, r0*asinh(u/r0)) in R4",
        "domain": [-L, L],
        "r0": float(r0),
        "radial_base": int(radial_base),
        "source_archives": SOURCE_ARCHIVES,
        "seed": {
            "original_vertices": int(seed.coordinates.shape[0]),
            "q6_edges": int(seed.q6_edges.shape[0]),
            "unique_angular_directions": int(seed.directions.shape[0]),
            "collapsed_duplicate_rays": int(
                seed.coordinates.shape[0] - seed.directions.shape[0]
            ),
            "q6_coframe_shape": list(seed.q6_unit_coframe.shape),
            "q6_coframe_reconstruction_residual": coframe_residual,
            "q6_coframe_unit_norm_residual": unit_residual,
            "archived_array_hashes_verified": True,
            "archived_array_sha256": SOURCE_ARRAY_SHA256,
            "q6_coframe_role": (
                "archived Q6 link data; retained and audited separately, "
                "not silently reinterpreted as tetrahedral finite-element edges"
            ),
        },
        "refinements": reports,
        "scope": {
            "finite_domain": "compact manifold with two boundary components",
            "mesh_metric": "induced chordal metric on R4 piecewise-linear tetrahedra",
            "continuum_metric": "exact Ellis-RGS warped-product metric",
            "Mosco_status": "OPEN for this curved family; G1 geometric construction only",
            "physical_status": (
                "r0 is taken from a diagnostic RGS Formulazione-4D profile; "
                "this mesh does not assert a full field-equation solution"
            ),
        },
    }

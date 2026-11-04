"""A5-equivariant finite relative-direction map from Q6 to S^2.

Each Q6 bit records the sign of one unoriented icosahedral axis. The signed
sum of the six axes is a nonzero vector in R^3; normalization gives a finite
set of directions on S^2. This is a 64-state map with 32 distinct directions,
each of multiplicity two. It is a discrete angular sampling, not a continuous
surjection Q6 -> S^2 and not a map from the five-dimensional Hodge eigenspace.
"""
from __future__ import annotations

import itertools
import numpy as np

from src.rgs_icosahedral_q6_intertwiner import (
    icosahedral_signed_axis_representation,
    q6_vertex_action,
)
from src.rgs_weyl_phase_closure import build_icosahedral_c6


def _sphere_moment_reference(order: int) -> np.ndarray:
    """Uniform S^2 moment tensor for even order 2,4,6; zero for odd order."""
    if order < 0:
        raise ValueError("order must be nonnegative")
    if order % 2:
        return np.zeros((3,) * order, dtype=float)
    if order == 0:
        return np.asarray(1.0)
    if order not in (2, 4, 6):
        raise ValueError("reference implemented for orders 0,2,4,6")
    tensor = np.zeros((3,) * order, dtype=float)
    indices = tuple(range(order))

    def pairings(items):
        if not items:
            yield ()
            return
        first = items[0]
        for j in range(1, len(items)):
            second = items[j]
            rest = items[1:j] + items[j + 1:]
            for tail in pairings(rest):
                yield ((first, second),) + tail

    denominator = {2: 3.0, 4: 15.0, 6: 105.0}[order]
    for pairing in pairings(indices):
        term = np.ones((3,) * order, dtype=float)
        for a, b in pairing:
            delta = np.eye(3)
            shape = [1] * order
            shape[a] = 3
            shape[b] = 3
            term *= delta.reshape(shape)
        tensor += term / denominator
    return tensor


def _unique_directions(directions: np.ndarray, tolerance: float = 1e-9):
    representatives = []
    multiplicities = []
    for direction in directions:
        match = next(
            (i for i, old in enumerate(representatives)
             if np.linalg.norm(direction - old) <= tolerance),
            None,
        )
        if match is None:
            representatives.append(direction.copy())
            multiplicities.append(1)
        else:
            multiplicities[match] += 1
    return np.asarray(representatives), np.asarray(multiplicities, dtype=int)


def _icosahedron_face_center_directions(vertices: np.ndarray) -> np.ndarray:
    """Return the 20 normalized face-centre directions of the icosahedron."""
    x = np.asarray(vertices, dtype=float)
    distances = np.linalg.norm(x[:, None, :] - x[None, :, :], axis=2)
    edge_length = float(np.min(distances[distances > 1e-7]))
    faces = []
    for i, j, k in itertools.combinations(range(len(x)), 3):
        if (
            abs(distances[i, j] - edge_length) < 1e-7
            and abs(distances[i, k] - edge_length) < 1e-7
            and abs(distances[j, k] - edge_length) < 1e-7
        ):
            faces.append((i, j, k))
    if len(faces) != 20:
        raise ArithmeticError(f"expected 20 icosahedron faces, got {len(faces)}")
    centers = np.asarray([x[list(face)].mean(axis=0) for face in faces])
    return centers / np.linalg.norm(centers, axis=1)[:, None]


def build_q6_a5_relative_sphere_map() -> dict:
    """Map all 64 Q6 vertices equivariantly to a 32-point subset of S^2."""
    group = icosahedral_signed_axis_representation()
    axes = np.asarray(group["axes"], dtype=float)  # (6,3), unit axes
    if axes.shape != (6, 3):
        raise ArithmeticError("expected six icosahedral axes in R^3")
    axis_matrix = axes.T
    if np.linalg.matrix_rank(axis_matrix) != 3:
        raise ArithmeticError("icosahedral axis frame must span R^3")

    bits = np.asarray(
        [[(v >> (5 - j)) & 1 for j in range(6)] for v in range(64)],
        dtype=np.int8,
    )
    signs = 1.0 - 2.0 * bits
    raw = signs @ axes
    radii = np.linalg.norm(raw, axis=1)
    if np.any(radii <= 1e-12):
        raise ArithmeticError("a Q6 sign state maps to the zero vector")
    directions = raw / radii[:, None]
    unique, multiplicities = _unique_directions(directions)
    icosa_data = build_icosahedral_c6()
    icosa_vertices = np.asarray(icosa_data["vertices"], dtype=float)
    icosa_vertices /= np.linalg.norm(icosa_vertices, axis=1)[:, None]
    dodeca_dual_vertices = _icosahedron_face_center_directions(icosa_vertices)
    vertex_matches = np.asarray([
        min(np.linalg.norm(direction - point) for point in icosa_vertices)
        for direction in unique
    ])
    dual_matches = np.asarray([
        min(np.linalg.norm(direction - point) for point in dodeca_dual_vertices)
        for direction in unique
    ])
    icosa_count = int(np.count_nonzero(vertex_matches < 1e-8))
    dual_count = int(np.count_nonzero(dual_matches < 1e-8))
    union_residual = float(max(
        np.max(np.minimum(vertex_matches, dual_matches)),
        np.max(np.linalg.norm(icosa_vertices[:, None, :] - unique[None, :, :], axis=2).min(axis=1)),
        np.max(np.linalg.norm(dodeca_dual_vertices[:, None, :] - unique[None, :, :], axis=2).min(axis=1)),
    ))

    # Group equivariance: x(g.v) = R_g x(v), before and after normalization.
    raw_residual = 0.0
    direction_residual = 0.0
    for rotation, signed in zip(group["rotations"], group["signed_permutations"]):
        for vertex in range(64):
            image = q6_vertex_action(vertex, signed)
            raw_residual = max(
                raw_residual,
                float(np.linalg.norm(raw[image] - rotation @ raw[vertex])),
            )
            direction_residual = max(
                direction_residual,
                float(np.linalg.norm(directions[image] - rotation @ directions[vertex])),
            )

    moments = {}
    moment_residuals = {}
    for order in range(1, 7):
        letters = "abcdef"
        subscripts = ",".join("n" + letters[j] for j in range(order))
        tensor = np.einsum(
            subscripts + "->" + letters[:order],
            *([unique] * order),
            optimize=True,
        ) / len(unique)
        reference = _sphere_moment_reference(order)
        moments[str(order)] = tensor
        moment_residuals[str(order)] = float(np.linalg.norm(tensor - reference))

    # The 32-point equal-weight set is antipodal. Count antipodal pairs.
    antipodal_residual = max(
        min(float(np.linalg.norm(direction + other)) for other in unique)
        for direction in unique
    )
    antipodal_pairs = len(unique) // 2 if antipodal_residual < 1e-9 else 0

    # Tight-frame identity for the six unoriented icosahedral axes.
    frame_residual = float(np.linalg.norm(axis_matrix @ axis_matrix.T - 2.0*np.eye(3)))
    checks = {
        "axis_frame_rank_three": int(np.linalg.matrix_rank(axis_matrix)) == 3,
        "axis_frame_tight": frame_residual < 1e-10,
        "no_zero_vertex_image": bool(np.all(radii > 1e-12)),
        "all_64_states_used": len(directions) == 64,
        "32_distinct_sphere_directions": len(unique) == 32,
        "uniform_two_to_one_fibres": bool(np.all(multiplicities == 2)),
        "antipodal_16_pair_structure": antipodal_pairs == 16,\n        "icosahedron_plus_dual_dodecahedron": icosa_count == 12 and dual_count == 20 and union_residual < 1e-8,
        "A5_raw_map_equivariant": raw_residual < 1e-9,
        "A5_sphere_map_equivariant": direction_residual < 1e-9,
        "second_moment_isotropic": moment_residuals["2"] < 1e-10,
        "fourth_moment_isotropic": moment_residuals["4"] < 1e-10,
        "sixth_moment_detects_finite_anisotropy": moment_residuals["6"] > 1e-5,
    }
    return {
        "audit": "Q6_A5_EQUIVARIANT_RELATIVE_SPHERE_MAP_V1",
        "q6_state_count": 64,
        "icosahedral_axis_count": 6,
        "linear_map_shape": [3, 6],
        "linear_map_rank": int(np.linalg.matrix_rank(axis_matrix)),
        "linear_map_kernel_dimension": 3,
        "group_order": int(group["group_order"]),
        "raw_vectors": raw,
        "radii": radii,
        "directions": directions,
        "unique_directions": unique,
        "fibre_multiplicities": multiplicities,
        "unique_direction_count": int(len(unique)),
        "antipodal_pair_count": int(antipodal_pairs),\n        "icosahedron_vertex_direction_count": icosa_count,\n        "dual_dodecahedron_vertex_direction_count": dual_count,\n        "icosahedron_dual_union_match_residual": union_residual,
        "radial_shells": {
            str(float(value)): int(np.count_nonzero(np.isclose(radii, value, atol=1e-9)))
            for value in np.unique(np.round(radii, 10))
        },
        "axis_frame_residual": frame_residual,
        "A5_raw_equivariance_residual": raw_residual,
        "A5_direction_equivariance_residual": direction_residual,
        "moment_residuals_to_uniform_S2": moment_residuals,
        "design_strength_at_least": 5 if (
            moment_residuals["1"] < 1e-10
            and moment_residuals["2"] < 1e-10
            and moment_residuals["3"] < 1e-10
            and moment_residuals["4"] < 1e-10
            and moment_residuals["5"] < 1e-10
        ) else 1,
        "sixth_moment_anisotropy": moment_residuals["6"],
        "checks": checks,
        "all_checks_passed": bool(all(checks.values())),
        "epistemic_status": {
            "numerically_verified": [
                "A5-equivariant map from 64 Q6 vertices to a finite subset of S2",
                "32 distinct directions with uniform multiplicity two",
                "16 antipodal direction pairs",
                "equal-weight spherical 5-design moment identities through degree five",
                "nonzero sixth-moment anisotropy",
            ],
            "not_proved": [
                "continuous or surjective Q6-to-S2 map",
                "identification of this direction map with the five-dimensional coexact Hodge eigenspace",
                "canonical coupling of the finite direction set to Smilga's continuous two-body phase-space measure",
                "global embedding of the resulting kinematic shell into the D5 Lie ball",
            ],
        },
    }


def q6_vertex_relative_direction(vertex: int) -> np.ndarray:
    """Return the normalized icosahedral-axis sum for one six-bit Q6 state."""
    if not isinstance(vertex, (int, np.integer)) or not 0 <= int(vertex) < 64:
        raise ValueError("vertex must be an integer in 0..63")
    axes = np.asarray(
        icosahedral_signed_axis_representation()["axes"], dtype=float
    )
    signs = np.asarray(
        [1.0 - 2.0*((int(vertex) >> (5-j)) & 1) for j in range(6)]
    )
    vector = signs @ axes
    norm = float(np.linalg.norm(vector))
    if norm <= 1e-12:
        raise ArithmeticError("Q6 state has zero relative-direction image")
    return vector / norm

"""Bridge the recovered RGS O_h action to SU(2) links and Q6 Wilson loops.

The code keeps three actions separate:
1. the full geometric O_h action on the 64 recovered points;
2. the subgroup H <= O_h preserving the canonical Q6 edges and squares;
3. the proper subgroup H+ = H intersect SO(3), which has an SU(2) lift.

SU(2) link variables are independent connection data in SU(2)^192.
The geometric seed determines their base-graph symmetry, not their values.
"""
from __future__ import annotations

from itertools import permutations, product
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

from src.holonomies import SU2Lattice, dagger
from src.q6_complex import build_q6_complex
from src.quaternions import Quaternion
from src.rotations import quaternion_to_matrix
from src.rgs_q6_oh_action import (
    induced_vertex_permutation,
    signed_permutation_matrices,
)
from src.rgs_ellis_refinement import load_rgs_seed_data
from src.su2_lifts import quaternion_to_su2


def _key(matrix: np.ndarray) -> tuple[int, ...]:
    return tuple(int(v) for v in np.asarray(matrix).ravel())


def _edge_set(edges):
    return {tuple(sorted((int(a), int(b)))) for a, b in np.asarray(edges)[:, :2]}


def _face_set(faces):
    return {tuple(sorted(map(int, f))) for f in faces}


def _preserves_q6(permutation: np.ndarray, edges, faces) -> tuple[bool, bool]:
    es, fs = _edge_set(edges), _face_set(faces)
    edge_ok = all(
        tuple(sorted((int(permutation[a]), int(permutation[b])))) in es
        for a, b in es
    )
    face_ok = all(
        tuple(sorted(int(permutation[v]) for v in face)) in fs
        for face in faces
    )
    return edge_ok, face_ok


def recovered_q6_symmetry_subgroups(tolerance: float = 1e-10) -> dict:
    """Recover H={S in O_h: S preserves Q6 edges and square faces}."""
    seed = load_rgs_seed_data()
    q6 = build_q6_complex()
    full, graph, proper = [], [], []
    for S in signed_permutation_matrices():
        p, _ = induced_vertex_permutation(seed.coordinates, S, tolerance)
        e_ok, f_ok = _preserves_q6(p, q6.edges, q6.faces)
        if e_ok and f_ok:
            item = {"matrix": S, "permutation": p}
            full.append(item)
            if round(float(np.linalg.det(S))) == 1:
                proper.append(item)
        if e_ok:
            graph.append(S)
    return {
        "O_h_order": 48,
        "Q6_automorphism_order": len(full),
        "Q6_graph_automorphism_order": len(graph),
        "proper_Q6_subgroup_order": len(proper),
        "H": full,
        "H_plus": proper,
    }


def _lift_proper_signed_permutation(S: np.ndarray) -> np.ndarray:
    """Choose one SU(2) lift of a proper signed-permutation rotation."""
    S = np.asarray(S, dtype=float)
    if S.shape != (3, 3) or not np.allclose(S.T @ S, np.eye(3), atol=0, rtol=0):
        raise ValueError("expected an exact orthogonal 3x3 matrix")
    if round(float(np.linalg.det(S))) != 1:
        raise ValueError("an improper O(3) transformation has no SU(2) lift")
    # scipy returns (x,y,z,w); the repository Quaternion is scalar-first.
    x, y, z, w = Rotation.from_matrix(S).as_quat()
    q = Quaternion(float(w), float(x), float(y), float(z)).normalized()
    if np.linalg.norm(quaternion_to_matrix(q) - S, ord="fro") > 1e-9:
        raise ArithmeticError("SU(2) lift does not project to the requested SO(3)")
    return quaternion_to_su2(q)


def _random_su2(rng: np.random.Generator) -> np.ndarray:
    qv = rng.normal(size=4)
    qv /= np.linalg.norm(qv)
    return quaternion_to_su2(Quaternion(*map(float, qv)))


def random_q6_connection(seed: int = 20261103) -> SU2Lattice:
    """Create a reproducible generic connection with one SU(2) matrix per edge."""
    rng = np.random.default_rng(seed)
    q6 = build_q6_complex()
    lattice = SU2Lattice(q6.vertices)
    for u, v in q6.edges:
        lattice.add_link(int(u), int(v), _random_su2(rng))
    return lattice


def plaquette_holonomies(lattice: SU2Lattice) -> tuple[np.ndarray, ...]:
    """Return the 240 ordered square holonomies in canonical Q6 orientation."""
    q6 = build_q6_complex()
    return tuple(
        lattice.holonomy((int(v0), int(v1), int(v2), int(v3), int(v0)))
        for v0, v1, v2, v3 in q6.faces
    )


def wilson_action(lattice: SU2Lattice, beta: float = 1.0) -> float:
    """S_W = beta sum_p [1 - (1/2) Re Tr(U_p)] for fundamental SU(2)."""
    beta = float(beta)
    if not np.isfinite(beta) or beta < 0:
        raise ValueError("beta must be finite and non-negative")
    total = 0.0
    for H in plaquette_holonomies(lattice):
        tr = np.trace(H)
        if abs(float(tr.imag)) > 1e-9:
            raise ArithmeticError("SU(2) plaquette trace is not real")
        total += 1.0 - 0.5 * float(tr.real)
    return beta * total


def spatial_symmetry_transform(
    lattice: SU2Lattice,
    permutation: np.ndarray,
    spin_lift: np.ndarray | None = None,
) -> SU2Lattice:
    """Push links through a Q6 automorphism; optionally conjugate by its spin lift.

    If an edge orientation reverses under the vertex permutation, the image
    link is inverted before insertion in the canonical lower-to-higher order.
    """
    q6 = build_q6_complex()
    p = np.asarray(permutation, dtype=np.int64)
    if p.shape != (64,) or set(map(int, p)) != set(range(64)):
        raise ValueError("permutation must be a bijection of 64 vertices")
    e_ok, f_ok = _preserves_q6(p, q6.edges, q6.faces)
    if not (e_ok and f_ok):
        raise ValueError("vertex permutation is not an automorphism of Q6's 2-complex")
    V = np.eye(2, dtype=complex) if spin_lift is None else np.asarray(spin_lift, dtype=complex)
    if V.shape != (2, 2):
        raise ValueError("spin_lift must be 2x2")
    out = SU2Lattice(q6.vertices)
    for u, v in q6.edges:
        a, b = int(p[u]), int(p[v])
        image = V @ lattice.link(int(u), int(v)) @ dagger(V)
        if a < b:
            out.add_link(a, b, image)
        else:
            out.add_link(b, a, dagger(image))
    return out


def audit_rgs_q6_su2_wilson(
    seed: int = 20261103, beta: float = 1.0
) -> dict:
    """Run geometry, group-lift, link, plaquette and Wilson invariance checks."""
    q6 = build_q6_complex()
    subgroups = recovered_q6_symmetry_subgroups()
    lattice = random_q6_connection(seed)
    action = wilson_action(lattice, beta)
    holonomies = plaquette_holonomies(lattice)

    rng = np.random.default_rng(seed + 1)
    gauge = {v: _random_su2(rng) for v in q6.vertices}
    gauge_lattice = lattice.gauge_transform(gauge)
    gauge_action = wilson_action(gauge_lattice, beta)
    gauge_residual = abs(gauge_action - action)

    gauge_covariance = 0.0
    for face, H0, H1 in zip(q6.faces, holonomies, plaquette_holonomies(gauge_lattice)):
        base = int(face[0])
        expected = gauge[base] @ H0 @ dagger(gauge[base])
        gauge_covariance = max(
            gauge_covariance, float(np.linalg.norm(H1 - expected, ord="fro"))
        )

    # Test every proper Q6 symmetry, not just a hand-picked generator.
    proper_action_residuals = []
    lift_projection_residuals = []
    for item in subgroups["H_plus"]:
        S = item["matrix"]
        p = item["permutation"]
        lift = _lift_proper_signed_permutation(S)
        lift_projection_residuals.append(
            float(np.linalg.norm(quaternion_to_matrix(
                Quaternion(*map(float, Rotation.from_matrix(S).as_quat()[[3,0,1,2]]))
            ) - S, ord="fro"))
        )
        moved = spatial_symmetry_transform(lattice, p, lift)
        proper_action_residuals.append(abs(wilson_action(moved, beta) - action))

    # All 16 base automorphisms preserve the Wilson sum without requiring
    # an SU(2) spin lift; this includes the improper elements.
    base_action_residuals = []
    for item in subgroups["H"]:
        moved = spatial_symmetry_transform(lattice, item["permutation"])
        base_action_residuals.append(abs(wilson_action(moved, beta) - action))

    return {
        "audit": "RGS_Q6_SU2_192_LINK_240_WILSON_BRIDGE",
        "dimensions": {
            "vertices": len(q6.vertices),
            "oriented_independent_links": len(q6.edges),
            "plaquettes": len(q6.faces),
            "configuration_space": "SU(2)^192",
            "gauge_group": "SU(2)^64",
        },
        "symmetry": {
            "O_h_order": subgroups["O_h_order"],
            "Q6_2_complex_automorphism_order": subgroups["Q6_automorphism_order"],
            "proper_Q6_subgroup_order": subgroups["proper_Q6_subgroup_order"],
            "improper_Q6_automorphisms": (
                subgroups["Q6_automorphism_order"]
                - subgroups["proper_Q6_subgroup_order"]
            ),
            "all_48_Oh_preserve_Q6": (
                subgroups["Q6_automorphism_order"] == 48
            ),
        },
        "connection": {
            "all_links_in_SU2": all(
                np.linalg.norm(dagger(lattice.link(u, v)) @ lattice.link(u, v)
                               - np.eye(2)) < 1e-10
                and abs(np.linalg.det(lattice.link(u, v)) - 1) < 1e-10
                for u, v in q6.edges
            ),
            "link_count": len(lattice.edges()),
        },
        "plaquettes": {
            "count": len(holonomies),
            "all_holonomies_in_SU2": all(
                np.linalg.norm(dagger(H) @ H - np.eye(2)) < 1e-9
                and abs(np.linalg.det(H) - 1) < 1e-9 for H in holonomies
            ),
        },
        "wilson": {
            "beta": float(beta),
            "action": float(action),
            "gauge_transformed_action": float(gauge_action),
            "gauge_action_residual": float(gauge_residual),
            "plaquette_gauge_covariance_max_frobenius": float(gauge_covariance),
            "proper_spin_lift_max_projection_residual": max(lift_projection_residuals, default=0.0),
            "proper_Q6_symmetry_action_max_residual": max(proper_action_residuals, default=0.0),
            "all_Q6_base_symmetry_action_max_residual": max(base_action_residuals, default=0.0),
        },
        "status": {
            "X_Oh_to_Q6_full_equivariance": "OBSTRUCTED: only 16 of 48 preserve the Q6 2-complex",
            "Q6_to_SU2_links": "CONSTRUCTED: independent SU(2) connection variables on 192 edges",
            "SU2_gauge_covariance": "NUMERICALLY_VERIFIED",
            "plaquette_holonomy": "NUMERICALLY_VERIFIED",
            "Wilson_gauge_invariance": "NUMERICALLY_VERIFIED",
            "proper_subgroup_SU2_lift": "NUMERICALLY_VERIFIED",
            "full_Oh_spinorial_lift": "NOT_AVAILABLE in SU(2); improper elements require Pin(3) or a separate base action",
            "continuum_Yang_Mills_limit": "OPEN",
        },
        "pass": bool(
            subgroups["Q6_automorphism_order"] == 16
            and subgroups["proper_Q6_subgroup_order"] == 8
            and len(lattice.edges()) == 192
            and len(holonomies) == 240
            and gauge_residual < 1e-9
            and gauge_covariance < 1e-9
            and max(proper_action_residuals, default=0.0) < 1e-8
            and max(base_action_residuals, default=0.0) < 1e-8
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=20261103)
    parser.add_argument("--beta", type=float, default=1.0)
    parser.add_argument("--output", type=Path,
                        default=Path("results/json/rgs_q6_su2_wilson_bridge.json"))
    args = parser.parse_args()
    report = audit_rgs_q6_su2_wilson(args.seed, args.beta)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if not report["pass"]:
        raise SystemExit("RGS-Q6-SU2-Wilson bridge audit failed")


if __name__ == "__main__":
    main()

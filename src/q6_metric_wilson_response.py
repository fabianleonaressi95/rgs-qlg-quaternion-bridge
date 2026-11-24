"""Metric-dependent Q6 plaquette response from the recovered RGS embedding.

This is a conditional finite face-Hodge model: S_g=beta sum_p
(A_ref/A_p) [1 - 1/2 Re Tr U_p], with A_p the parallelogram area of the
canonical embedded face and A_ref fixed at the reference configuration.
The missing circumcentric-dual measure is deliberately not inferred.
"""
from __future__ import annotations

import numpy as np

from src.holonomies import SU2Lattice
from src.q6_complex import build_q6_complex
from src.rgs_ellis_refinement import load_rgs_seed_data


def shear_generator_3d(axis_x: int = 0, axis_y: int = 1) -> np.ndarray:
    """Symmetric traceless ambient shear H=(E_xy+E_yx)/2 in R^3."""
    if axis_x == axis_y or not (0 <= axis_x < 3 and 0 <= axis_y < 3):
        raise ValueError("shear axes must be distinct members of {0,1,2}")
    h = np.zeros((3, 3), dtype=float)
    h[axis_x, axis_y] = h[axis_y, axis_x] = 0.5
    return h


def face_parallelogram_areas(coordinates: np.ndarray) -> np.ndarray:
    """Areas from the two canonical base edges of each embedded Q6 square."""
    x = np.asarray(coordinates, dtype=float)
    if x.shape != (64, 3) or not np.all(np.isfinite(x)):
        raise ValueError("coordinates must be a finite (64,3) array")
    q6 = build_q6_complex()
    areas = np.empty(len(q6.faces), dtype=float)
    for i, (v, va, vab, vb) in enumerate(q6.faces):
        a = x[va] - x[v]
        b = x[vb] - x[v]
        gram = np.array([[a @ a, a @ b], [a @ b, b @ b]], dtype=float)
        det = float(np.linalg.det(gram))
        if det <= 1e-24:
            raise ValueError(f"degenerate embedded face at index {i}")
        areas[i] = np.sqrt(det)
    return areas


def face_area_log_shear_derivative(
    coordinates: np.ndarray, generator: np.ndarray
) -> np.ndarray:
    """Exact d(log A_p)/dgamma at zero for X(gamma)=(I+gamma H)X."""
    x = np.asarray(coordinates, dtype=float)
    h = np.asarray(generator, dtype=float)
    if x.shape != (64, 3) or h.shape != (3, 3):
        raise ValueError("expected coordinates (64,3) and generator (3,3)")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(h)):
        raise ValueError("coordinates and generator must be finite")
    if not np.allclose(h, h.T, atol=1e-13) or abs(np.trace(h)) > 1e-12:
        raise ValueError("generator must be symmetric and traceless")
    q6 = build_q6_complex()
    out = np.empty(len(q6.faces), dtype=float)
    for i, (v, va, vab, vb) in enumerate(q6.faces):
        a = x[va] - x[v]
        b = x[vb] - x[v]
        gram = np.array([[a @ a, a @ b], [a @ b, b @ b]], dtype=float)
        da = h @ a
        db = h @ b
        dgram = np.array([
            [2.0 * (a @ da), da @ b + a @ db],
            [da @ b + a @ db, 2.0 * (b @ db)],
        ])
        if np.linalg.det(gram) <= 1e-24:
            raise ValueError(f"degenerate embedded face at index {i}")
        out[i] = 0.5 * float(np.trace(np.linalg.solve(gram, dgram)))
    return out


def metric_weighted_plaquette_action(
    coordinates: np.ndarray,
    links: SU2Lattice,
    beta: float = 1.0,
    reference_areas: np.ndarray | None = None,
) -> dict[str, np.ndarray | float]:
    """Evaluate conditional face-Hodge Wilson action with fixed reference areas."""
    q6 = build_q6_complex()
    x = np.asarray(coordinates, dtype=float)
    areas = face_parallelogram_areas(x)
    ref = areas.copy() if reference_areas is None else np.asarray(reference_areas, float)
    if ref.shape != areas.shape or np.any(~np.isfinite(ref)) or np.any(ref <= 0):
        raise ValueError("reference_areas must be 240 positive finite values")
    if not np.isfinite(beta) or beta < 0:
        raise ValueError("beta must be finite and non-negative")
    plaquette = np.empty(len(q6.faces), dtype=float)
    for i, face in enumerate(q6.faces):
        hol = links.holonomy(tuple(map(int, (*face, face[0]))))
        plaquette[i] = 1.0 - 0.5 * float(np.trace(hol).real)
    weights = ref / areas
    terms = float(beta) * weights * plaquette
    return {
        "action": float(np.sum(terms)),
        "areas": areas,
        "reference_areas": ref.copy(),
        "weights": weights,
        "plaquette_energies": plaquette,
        "terms": terms,
    }


def metric_wilson_shear_response(
    links: SU2Lattice,
    reference_volume: float,
    beta: float = 1.0,
    axis_x: int = 0,
    axis_y: int = 1,
    coordinates: np.ndarray | None = None,
) -> dict[str, object]:
    """Analytic -V^-1 dS_g/dgamma for the conditional metric-face action.

    RGS coordinates and reference face areas are fixed at gamma=0. Under
    shear, only current face areas vary; link variables and plaquette
    holonomies are held fixed. The result is not yet a full 3D DEC Hodge
    stress because no dual-face measure is specified by the seed data.
    """
    if not np.isfinite(reference_volume) or reference_volume <= 0:
        raise ValueError("reference_volume must be finite and positive")
    if coordinates is None:
        coordinates = load_rgs_seed_data().coordinates
    x = np.asarray(coordinates, dtype=float)
    h = shear_generator_3d(axis_x, axis_y)
    base = metric_weighted_plaquette_action(x, links, beta)
    dlog_area = face_area_log_shear_derivative(x, h)
    dS = -float(beta) * float(np.dot(base["weights"] * base["plaquette_energies"], dlog_area))
    return {
        "action": base["action"],
        "d_action_dgamma": dS,
        "Pi_xy": float(-dS / reference_volume),
        "face_areas": base["areas"],
        "face_weights": base["weights"],
        "plaquette_energies": base["plaquette_energies"],
        "dlog_face_area_dgamma": dlog_area,
        "reference_volume": float(reference_volume),
        "status": "exact derivative of the declared finite face-area-weighted plaquette action",
        "physical_status": "conditional face-Hodge candidate; dual measure and equilibrium dynamics remain open",
    }


def metric_wilson_shear_finite_difference(
    links: SU2Lattice,
    reference_volume: float,
    beta: float = 1.0,
    axis_x: int = 0,
    axis_y: int = 1,
    step: float = 1e-6,
    coordinates: np.ndarray | None = None,
) -> float:
    """Independent central-difference check, with reference areas held fixed."""
    if not np.isfinite(step) or step <= 0:
        raise ValueError("step must be finite and positive")
    x = np.asarray(load_rgs_seed_data().coordinates if coordinates is None else coordinates, float)
    h = shear_generator_3d(axis_x, axis_y)
    ref = face_parallelogram_areas(x)
    plus = metric_weighted_plaquette_action(((np.eye(3) + step*h) @ x.T).T, links, beta, ref)["action"]
    minus = metric_weighted_plaquette_action(((np.eye(3) - step*h) @ x.T).T, links, beta, ref)["action"]
    return float(-(plus - minus) / (2.0 * step * reference_volume))


def metric_edge_matter_action(
    coordinates: np.ndarray,
    fields,
    links: SU2Lattice,
    coupling: float = 1.0,
) -> float:
    """P1 graph Dirichlet matter action: kappa/2 sum_e ||D_e psi||^2/L_e.

    The inverse edge length is the exact 1D affine-P1 stiffness coefficient
    already used by the repository's RGS-Q6 induced-metric FEM. Both PSU
    sheets are included. This is a graph scalar-gradient action, not a Dirac
    action; fields are expected to expose psi_plus and psi_minus arrays.
    """
    x = np.asarray(coordinates, dtype=float)
    if x.shape != (64, 3) or not np.all(np.isfinite(x)):
        raise ValueError("coordinates must be finite with shape (64,3)")
    if not np.isfinite(coupling) or coupling < 0:
        raise ValueError("coupling must be finite and non-negative")
    q6 = build_q6_complex()
    plus = np.asarray(fields.psi_plus, dtype=complex)
    minus = np.asarray(fields.psi_minus, dtype=complex)
    if plus.shape != (64, 2) or minus.shape != (64, 2):
        raise ValueError("both matter sheets must have shape (64,2)")
    total = 0.0
    for u, v in q6.edges:
        length = float(np.linalg.norm(x[v] - x[u]))
        if length <= 1e-13:
            raise ValueError("zero-length embedded Q6 edge")
        U = links.link(int(u), int(v))
        dp = plus[v] - U @ plus[u]
        dm = minus[v] - U @ minus[u]
        total += (float(np.vdot(dp, dp).real) + float(np.vdot(dm, dm).real)) / length
    return 0.5 * float(coupling) * total


def metric_edge_matter_shear_response(
    fields,
    links: SU2Lattice,
    reference_volume: float,
    coupling: float = 1.0,
    axis_x: int = 0,
    axis_y: int = 1,
    coordinates: np.ndarray | None = None,
) -> dict[str, float]:
    """Exact -V^-1 shear derivative of the P1 inverse-edge-length action."""
    if not np.isfinite(reference_volume) or reference_volume <= 0:
        raise ValueError("reference_volume must be finite and positive")
    x = np.asarray(load_rgs_seed_data().coordinates if coordinates is None else coordinates, float)
    if x.shape != (64, 3) or not np.all(np.isfinite(x)):
        raise ValueError("coordinates must be finite with shape (64,3)")
    if not np.isfinite(coupling) or coupling < 0:
        raise ValueError("coupling must be finite and non-negative")
    h = shear_generator_3d(axis_x, axis_y)
    q6 = build_q6_complex()
    plus = np.asarray(fields.psi_plus, dtype=complex)
    minus = np.asarray(fields.psi_minus, dtype=complex)
    if plus.shape != (64, 2) or minus.shape != (64, 2):
        raise ValueError("both matter sheets must have shape (64,2)")
    dS = 0.0
    for u, v in q6.edges:
        d = x[v] - x[u]
        length2 = float(d @ d)
        if length2 <= 1e-26:
            raise ValueError("zero-length embedded Q6 edge")
        U = links.link(int(u), int(v))
        dp = plus[v] - U @ plus[u]
        dm = minus[v] - U @ minus[u]
        energy = float(np.vdot(dp, dp).real + np.vdot(dm, dm).real)
        dlog_length = float(d @ h @ d) / length2
        dS -= 0.5 * float(coupling) * energy / np.sqrt(length2) * dlog_length
    return {
        "d_action_dgamma": float(dS),
        "Pi_xy": float(-dS / reference_volume),
        "reference_volume": float(reference_volume),
        "status": "exact derivative of the 1D P1 graph Dirichlet matter action",
        "physical_status": "conditional graph scalar-gradient sector; continuum matter interpretation remains open",
    }


def metric_edge_matter_shear_finite_difference(
    fields,
    links: SU2Lattice,
    reference_volume: float,
    coupling: float = 1.0,
    axis_x: int = 0,
    axis_y: int = 1,
    step: float = 1e-6,
    coordinates: np.ndarray | None = None,
) -> float:
    """Central-difference audit with the link and matter fields held fixed."""
    if not np.isfinite(step) or step <= 0:
        raise ValueError("step must be finite and positive")
    x = np.asarray(load_rgs_seed_data().coordinates if coordinates is None else coordinates, float)
    h = shear_generator_3d(axis_x, axis_y)
    plus = metric_edge_matter_action(((np.eye(3)+step*h) @ x.T).T, fields, links, coupling)
    minus = metric_edge_matter_action(((np.eye(3)-step*h) @ x.T).T, fields, links, coupling)
    return float(-(plus-minus)/(2.0*step*reference_volume))


def q6_candidate_total_shear_response(
    fields,
    links: SU2Lattice,
    reference_volume: float,
    beta: float = 1.0,
    matter_coupling: float = 1.0,
    axis_x: int = 0,
    axis_y: int = 1,
    coordinates: np.ndarray | None = None,
) -> dict[str, object]:
    """Sum the declared conditional gauge-face and P1 matter shear sectors."""
    gauge = metric_wilson_shear_response(
        links, reference_volume, beta, axis_x, axis_y, coordinates,
    )
    matter = metric_edge_matter_shear_response(
        fields, links, reference_volume, matter_coupling,
        axis_x, axis_y, coordinates,
    )
    return {
        "Pi_gauge": float(gauge["Pi_xy"]),
        "Pi_matter": float(matter["Pi_xy"]),
        "Pi_total": float(gauge["Pi_xy"]) + float(matter["Pi_xy"]),
        "action_gauge": float(gauge["action"]),
        "action_matter": float(metric_edge_matter_action(
            np.asarray(load_rgs_seed_data().coordinates if coordinates is None else coordinates, float),
            fields, links, matter_coupling,
        )),
        "status": "sum of two exactly differentiated finite candidate sectors",
        "physical_status": "not yet a continuum-normalized stress; dual measure and equilibrium dynamics remain open",
    }

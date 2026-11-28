"""Manufactured-solution audit of local Cartan face consistency.

A smooth 4D Lorentzian Cartan pair is sampled on a shrinking coordinate
triangle. Edge affine transports solve dG/dt = G A_t, with
A_t = [[omega_t, e_t], [0, 0]]. The face product is compared with the
base-point curvature and torsion components.

This is a numerical consistency audit, not a proof of compactness or
convergence of arbitrary discrete solutions.
"""
from __future__ import annotations

import argparse
import json
import numpy as np


def boost_generator(spatial_axis: int) -> np.ndarray:
    """Generator of a (+---) Lorentz boost, spatial_axis in {1,2,3}."""
    if spatial_axis not in (1, 2, 3):
        raise ValueError("spatial_axis must be 1, 2, or 3")
    out = np.zeros((4, 4), dtype=float)
    out[0, spatial_axis] = 1.0
    out[spatial_axis, 0] = 1.0
    return out


KX = boost_generator(1)
KY = boost_generator(2)
A_CURV = 0.7 * KX
B_CURV = 0.4 * KY
TORSION_VECTOR = np.array([0.3, 0.2, -0.1, 0.4], dtype=float)
IDENTITY5 = np.eye(5)


def affine_generator(x: float, y: float, dx: float, dy: float) -> np.ndarray:
    """Pull back the smooth Cartan connection to a straight edge tangent."""
    omega_x = y * B_CURV
    omega_y = x * A_CURV
    e_x = np.array([1.0, 0.0, 0.0, 0.0])
    e_y = np.array([0.0, 1.0, 0.0, 0.0]) + x * TORSION_VECTOR

    out = np.zeros((5, 5), dtype=float)
    out[:4, :4] = dx * omega_x + dy * omega_y
    out[:4, 4] = dx * e_x + dy * e_y
    return out


def edge_transport(
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    steps: int = 64,
) -> np.ndarray:
    """Integrate dG/dt=G A_t along one oriented straight edge by RK4."""
    if steps < 1:
        raise ValueError("steps must be positive")
    p = np.asarray(start, dtype=float)
    delta = np.asarray(end, dtype=float) - p
    dt = 1.0 / steps
    G = IDENTITY5.copy()

    def rhs(t: float, value: np.ndarray) -> np.ndarray:
        point = p + t * delta
        return value @ affine_generator(
            float(point[0]), float(point[1]), float(delta[0]), float(delta[1])
        )

    for index in range(steps):
        t = index * dt
        k1 = rhs(t, G)
        k2 = rhs(t + 0.5 * dt, G + 0.5 * dt * k1)
        k3 = rhs(t + 0.5 * dt, G + 0.5 * dt * k2)
        k4 = rhs(t + dt, G + dt * k3)
        G = G + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    return G


def matrix_log_near_identity(matrix: np.ndarray, terms: int = 12) -> np.ndarray:
    """Evaluate log(I+X) by its convergent power series for small ||X||."""
    X = np.asarray(matrix, dtype=float) - np.eye(matrix.shape[0])
    result = np.zeros_like(X)
    power = np.eye(X.shape[0])
    for degree in range(1, terms + 1):
        power = power @ X
        result += ((-1.0) ** (degree + 1) / degree) * power
    return result


def face_observables(h: float, *, steps: int = 64) -> tuple[np.ndarray, np.ndarray]:
    """Return log Lorentz holonomy and affine translation for [0,(h,0),(0,h)]."""
    p0 = (0.0, 0.0)
    p1 = (h, 0.0)
    p2 = (0.0, h)
    G01 = edge_transport(p0, p1, steps=steps)
    G12 = edge_transport(p1, p2, steps=steps)
    G20 = edge_transport(p2, p0, steps=steps)
    face = G01 @ G12 @ G20
    return matrix_log_near_identity(face[:4, :4]), face[:4, 4]


def convergence_table(
    scales: tuple[float, ...] = (0.2, 0.1, 0.05, 0.025, 0.0125),
    *,
    steps: int = 64,
) -> list[dict[str, float]]:
    """Compare face observables with exact base-point R_xy and T_xy."""
    # At (0,0): R_xy = d_x omega_y - d_y omega_x + [omega_x,omega_y]
    #         = A_CURV - B_CURV; T_xy = d_x e_y - d_y e_x
    #         + omega_x e_y - omega_y e_x = TORSION_VECTOR.
    R_xy = A_CURV - B_CURV
    T_xy = TORSION_VECTOR
    rows: list[dict[str, float]] = []
    previous_R = None
    previous_T = None

    for h in scales:
        log_H, tau = face_observables(h, steps=steps)
        target_R = 0.5 * h * h * R_xy
        target_T = 0.5 * h * h * T_xy
        err_R = float(np.linalg.norm(log_H - target_R))
        err_T = float(np.linalg.norm(tau - target_T))
        density_err_R = 2.0 * err_R / (h * h)
        density_err_T = 2.0 * err_T / (h * h)

        order_R = (
            float(np.log(previous_R / err_R) / np.log(2.0))
            if previous_R is not None and err_R > 0.0 else float("nan")
        )
        order_T = (
            float(np.log(previous_T / err_T) / np.log(2.0))
            if previous_T is not None and err_T > 0.0 else float("nan")
        )
        rows.append({
            "h": float(h),
            "curvature_face_error": err_R,
            "torsion_face_error": err_T,
            "curvature_density_error": density_err_R,
            "torsion_density_error": density_err_T,
            "curvature_observed_order": order_R,
            "torsion_observed_order": order_T,
        })
        previous_R, previous_T = err_R, err_T
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--steps", type=int, default=64)
    args = parser.parse_args()
    rows = convergence_table(steps=args.steps)
    print(json.dumps({
        "model": "smooth noncommuting Lorentz connection with nonzero curvature and torsion",
        "edge_integrator": "RK4",
        "steps_per_edge": args.steps,
        "rows": rows,
        "interpretation": (
            "manufactured-solution local consistency only; no convergence theorem "
            "for arbitrary discrete fields or discrete solutions"
        ),
    }, indent=2))


if __name__ == "__main__":
    main()

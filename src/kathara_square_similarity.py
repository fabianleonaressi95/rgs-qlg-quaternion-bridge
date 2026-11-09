"""Exact square-similarity reconstruction of the Kathara/Krystal stage rule.

The module formalizes only the stated planar rule: clockwise 45-degree
rotation plus scale sqrt(2), with the 2:1 grid proportion and the reported
center-9 / center-12 / center-6 overlay conditions. It does not assert
physical validity or derive the three-channel tribonacci operator.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class KatharaSquareAudit:
    transform: np.ndarray
    transform_squared: np.ndarray
    centers: dict[str, np.ndarray]
    residuals: dict[str, float | int]
    lattice_image_index: int
    lattice_parity_rule: str


def kathara_stage_transform() -> np.ndarray:
    """Return T = sqrt(2) R(-pi/4) = [[1,1],[-1,1]]."""
    return np.asarray([[1, 1], [-1, 1]], dtype=np.int64)


def kathara_centers_from_overlay(width: float = 1.0) -> dict[str, np.ndarray]:
    """Reconstruct the constrained centers relative to AzurA=(0,0).

    width is the distance from center 6 to center 7. The source proportion
    sets the central 1-to-12 column length to 2*width, with AzurA at its
    midpoint. The two overlay equations then determine center 9.
    """
    if not np.isfinite(width) or width <= 0:
        raise ValueError("width must be finite and positive")
    T = kathara_stage_transform().astype(float)
    p6 = np.asarray([-width / 2.0, 0.0])
    p7 = np.asarray([width / 2.0, 0.0])
    p1 = np.asarray([0.0, -width])
    p12 = np.asarray([0.0, width])
    p9 = T @ p6
    return {
        "AzurA": np.zeros(2),
        "1": p1,
        "6": p6,
        "7": p7,
        "9": p9,
        "12": p12,
    }


def audit_kathara_square_similarity(width: float = 1.0) -> KatharaSquareAudit:
    """Audit the exact square similarity and the reported overlay constraints."""
    T = kathara_stage_transform()
    T2 = T @ T
    centers = kathara_centers_from_overlay(width)
    w = float(width)
    p6, p9, p12 = centers["6"], centers["9"], centers["12"]
    gram = T.T @ T
    expected_T2 = np.asarray([[0, 2], [-2, 0]], dtype=np.int64)
    residuals = {
        "TtT_minus_2I_max_abs": int(np.max(np.abs(gram - 2*np.eye(2, dtype=int)))),
        "T2_minus_2Rminus90_max_abs": int(np.max(np.abs(T2 - expected_T2))),
        "stage1_center9_to_center12_max_abs": float(np.max(np.abs(T @ p9 - p12))),
        "stage2_center6_to_center12_max_abs": float(np.max(np.abs(T2 @ p6 - p12))),
        "center9_equals_T_center6_max_abs": float(np.max(np.abs(p9 - T @ p6))),
        "central_column_ratio_residual": float(
            abs(np.linalg.norm(centers["12"] - centers["1"]) / np.linalg.norm(centers["7"] - centers["6"]) - 2.0)
        ),
        "center6_center7_midpoint_residual": float(
            np.linalg.norm((p6 + centers["7"]) / 2.0 - centers["AzurA"])
        ),
        "stage_scale_ratio_6_to_9": float(np.linalg.norm(p9) / np.linalg.norm(p6)),
        "stage_scale_ratio_9_to_12": float(np.linalg.norm(p12) / np.linalg.norm(p9)),
    }
    return KatharaSquareAudit(
        transform=T,
        transform_squared=T2,
        centers=centers,
        residuals=residuals,
        lattice_image_index=abs(int(round(np.linalg.det(T)))),
        lattice_parity_rule="T(Z^2) = {(u,v) in Z^2 : u == v (mod 2)}",
    )


def apply_stage(z: np.ndarray, stages: int = 1) -> np.ndarray:
    """Apply the integer stage map T repeatedly to a 2-vector."""
    if stages < 0:
        raise ValueError("stages must be nonnegative")
    z = np.asarray(z)
    if z.shape != (2,):
        raise ValueError("z must be a 2-vector")
    return np.linalg.matrix_power(kathara_stage_transform(), stages) @ z

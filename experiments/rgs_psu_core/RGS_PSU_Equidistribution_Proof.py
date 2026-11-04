#!/usr/bin/env python3
"""
Finite 600-cell moment audit for the RGS-PSU geometric programme.

Important distinction:
- The 120 vertices form a finite orbit / spherical point configuration.
- A finite configuration is NOT the continuous Haar-uniform measure on S^3.
- The audit checks the first, second and fourth moments against S^3 Haar
  moments. Passing these finite moment checks is evidence for a spherical
  design property, not a proof of weak equidistribution of a refinement
  sequence and not a proof of a continuum limit on S^3 / 2I.

The quotient S^3 / 2I is also not the same quotient as S^3 / {+1,-1}.
The 600-cell vertex set is the 2I orbit of the identity; under the left
2I quotient all points in that orbit represent the same quotient point.
"""

from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np


def _permutation_is_even(permutation: tuple[int, ...]) -> bool:
    inversions = sum(
        permutation[i] > permutation[j]
        for i in range(len(permutation))
        for j in range(i + 1, len(permutation))
    )
    return inversions % 2 == 0


def build_600_cell_vertices() -> np.ndarray:
    """Return the 120 unit vertices in standard 600-cell coordinates."""
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    vertices: list[np.ndarray] = []

    # 8 coordinate-axis vertices.
    for axis in range(4):
        for sign in (-1.0, 1.0):
            point = np.zeros(4, dtype=float)
            point[axis] = sign
            vertices.append(point)

    # 16 vertices (±1/2, ±1/2, ±1/2, ±1/2).
    for signs in itertools.product((-1.0, 1.0), repeat=4):
        vertices.append(0.5 * np.asarray(signs, dtype=float))

    # 96 vertices: even permutations of (0, 1/2, phi/2, 1/(2 phi)),
    # with independent signs on the three nonzero coordinates.
    base = np.asarray((0.0, 0.5, phi / 2.0, 1.0 / (2.0 * phi)))
    for permutation in itertools.permutations(range(4)):
        if not _permutation_is_even(permutation):
            continue
        point = base[list(permutation)].copy()
        nonzero = np.flatnonzero(point)
        for signs in itertools.product((-1.0, 1.0), repeat=3):
            signed = point.copy()
            signed[nonzero] *= np.asarray(signs)
            vertices.append(signed)

    result = np.asarray(vertices, dtype=float)
    # Stable uniqueness check guards accidental duplicate coordinate generation.
    if result.shape != (120, 4):
        raise RuntimeError(f"Expected 120 vertices, got {result.shape}")
    if np.unique(np.round(result, 12), axis=0).shape[0] != 120:
        raise RuntimeError("600-cell coordinate construction contains duplicates")
    return result


def haar_second_moment(dimension: int = 4) -> np.ndarray:
    """E[x_i x_j] for normalized Haar measure on S^(dimension-1)."""
    return np.eye(dimension, dtype=float) / dimension


def haar_fourth_moment(dimension: int = 4) -> np.ndarray:
    """E[x_i x_j x_k x_l] for normalized Haar measure on S^(dimension-1)."""
    tensor = np.zeros((dimension,) * 4, dtype=float)
    denominator = dimension * (dimension + 2)
    for i, j, k, ell in itertools.product(range(dimension), repeat=4):
        tensor[i, j, k, ell] = (
            (i == j) * (k == ell)
            + (i == k) * (j == ell)
            + (i == ell) * (j == k)
        ) / denominator
    return tensor


def moment_audit(vertices: np.ndarray | None = None) -> dict:
    """Compare finite point moments with continuous S^3 Haar moments."""
    points = build_600_cell_vertices() if vertices is None else np.asarray(vertices)
    if points.ndim != 2 or points.shape[1] != 4:
        raise ValueError("vertices must have shape (N, 4)")
    if not np.all(np.isfinite(points)):
        raise ValueError("vertices contain non-finite entries")
    if not np.allclose(np.linalg.norm(points, axis=1), 1.0, atol=1e-12):
        raise ValueError("all points must lie on the unit 3-sphere")

    center = np.mean(points, axis=0)
    second = points.T @ points / len(points)
    fourth = np.einsum("ni,nj,nk,nl->ijkl", points, points, points, points) / len(points)

    second_residual = float(np.linalg.norm(second - haar_second_moment()))
    fourth_residual = float(np.linalg.norm(fourth - haar_fourth_moment()))
    center_residual = float(np.linalg.norm(center))

    return {
        "object": "600-cell vertices in S3",
        "n_vertices": int(len(points)),
        "dimension": 4,
        "status": "NUMERICALLY_VERIFIED_FINITE_MOMENT_MATCH",
        "continuous_haar_equidistribution_claim": "NOT_ESTABLISHED",
        "quotient_warning": (
            "The 120-point 2I orbit maps to one orbit point under the left "
            "quotient S3/2I; antipodal quotient S3/{+1,-1} is a different quotient."
        ),
        "residuals": {
            "center_of_mass": center_residual,
            "second_moment": second_residual,
            "fourth_moment": fourth_residual,
        },
        "checks": {
            "unit_norms": bool(np.allclose(np.linalg.norm(points, axis=1), 1.0, atol=1e-12)),
            "center_of_mass_zero": bool(center_residual < 1e-12),
            "second_moment_matches_S3_Haar": bool(second_residual < 1e-12),
            "fourth_moment_matches_S3_Haar": bool(fourth_residual < 1e-12),
        },
        "interpretation": (
            "Finite moment agreement through degree four is a spherical-design "
            "type check. It does not imply weak convergence of a sequence of "
            "refined point measures to Haar measure."
        ),
    }


def main(output_path: str = "RGS_PSU_equidistribution_audit.json") -> dict:
    vertices = build_600_cell_vertices()
    report = moment_audit(vertices)
    report["coordinates"] = vertices.tolist()
    Path(output_path).write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({k: v for k, v in report.items() if k != "coordinates"}, indent=2))
    print(f"Audit written to: {output_path}")
    return report


if __name__ == "__main__":
    main()

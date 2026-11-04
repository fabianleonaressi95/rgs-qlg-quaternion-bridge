"""Affine/similarity and symmetry audit for the RGS 64B/O_h seed.

This module tests, rather than assumes, an identification with the 64-point
golden cut-and-project image of the centred 6-cube.

Key obstruction: both centred clouds have isotropic population covariance.
Consequently any invertible affine map carrying one equally weighted
64-point cloud onto the other must have linear part proportional to an
orthogonal matrix. Thus affine equivalence would imply similarity
equivalence, which is rejected by the normalized pair-distance spectrum.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.spatial.distance import pdist

from .rgs_psu_cut_project import build_cut_project_data
from .rgs_ellis_refinement import load_rgs_seed_data


def _radial_shells(points: np.ndarray, decimals: int = 10) -> dict:
    radii = np.linalg.norm(points, axis=1)
    values, counts = np.unique(np.round(radii, decimals), return_counts=True)
    return {
        "radii": values.tolist(),
        "multiplicities": counts.astype(int).tolist(),
        "number_of_shells": int(len(values)),
    }


def _isotropic_covariance(points: np.ndarray) -> tuple[np.ndarray, float, float]:
    centered = points - points.mean(axis=0, keepdims=True)
    covariance = centered.T @ centered / len(points)
    scalar = float(np.trace(covariance) / covariance.shape[0])
    residual = float(np.linalg.norm(covariance - scalar * np.eye(3)))
    relative = residual / max(float(np.linalg.norm(covariance)), 1e-15)
    return covariance, scalar, relative


def audit_rgs_golden_intertwiner(
    rgs_coordinates: np.ndarray | None = None,
    tolerance: float = 1e-10,
) -> dict:
    """Audit whether the RGS seed and projected Q6 cloud are affine-equivalent."""
    if rgs_coordinates is None:
        rgs_coordinates = load_rgs_seed_data().coordinates
    rgs = np.asarray(rgs_coordinates, dtype=np.float64)
    q6 = build_cut_project_data().q6_parallel
    if rgs.shape != (64, 3):
        raise ValueError("rgs_coordinates must have shape (64,3)")
    if not np.all(np.isfinite(rgs)):
        raise ValueError("rgs_coordinates must be finite")

    rgs0 = rgs - rgs.mean(axis=0, keepdims=True)
    q60 = q6 - q6.mean(axis=0, keepdims=True)
    cov_rgs, c_rgs, iso_rgs = _isotropic_covariance(rgs)
    cov_q6, c_q6, iso_q6 = _isotropic_covariance(q6)

    d_rgs = np.sort(pdist(rgs0))
    d_q6 = np.sort(pdist(q60))
    rms_rgs = float(np.sqrt(np.mean(d_rgs**2)))
    rms_q6 = float(np.sqrt(np.mean(d_q6**2)))
    dn_rgs = d_rgs / rms_rgs
    dn_q6 = d_q6 / rms_q6
    pair_residual = float(np.linalg.norm(dn_rgs - dn_q6) / np.sqrt(len(dn_rgs)))

    scale_squared = c_rgs / c_q6
    forced_similarity_scale = float(np.sqrt(scale_squared))
    affine_reduction_valid = iso_rgs < tolerance and iso_q6 < tolerance
    affine_equivalence_ruled_out = bool(
        affine_reduction_valid and pair_residual > tolerance
    )

    radial_rgs = np.linalg.norm(rgs0, axis=1) / np.sqrt(c_rgs)
    radial_q6 = np.linalg.norm(q60, axis=1) / np.sqrt(c_q6)
    radial_residual = float(
        np.linalg.norm(np.sort(radial_rgs) - np.sort(radial_q6)) / np.sqrt(64)
    )

    return {
        "audit": "RGS_GOLDEN_CUT_PROJECT_INTERTWINER_V1",
        "point_counts": {"rgs": 64, "golden_projected_q6": 64},
        "centroids": {
            "rgs_norm": float(np.linalg.norm(rgs.mean(axis=0))),
            "q6_norm": float(np.linalg.norm(q6.mean(axis=0))),
        },
        "population_covariance": {
            "rgs": cov_rgs.tolist(),
            "golden_projected_q6": cov_q6.tolist(),
            "rgs_isotropy_relative_residual": iso_rgs,
            "q6_isotropy_relative_residual": iso_q6,
            "rgs_scalar": c_rgs,
            "q6_scalar": c_q6,
            "forced_similarity_scale_q6_to_rgs": forced_similarity_scale,
        },
        "radial_shells": {
            "rgs": _radial_shells(rgs0),
            "golden_projected_q6": _radial_shells(q60),
            "normalized_sorted_radius_rms_residual": radial_residual,
        },
        "pair_distance_invariant": {
            "normalized_sorted_pair_distance_rms_residual": pair_residual,
            "pair_count": int(len(dn_rgs)),
        },
        "affine_no_go": {
            "covariance_argument_applies": affine_reduction_valid,
            "theorem": (
                "For centered equally weighted clouds with C_RGS=c_RGS I and "
                "C_Q6=c_Q6 I, any affine bijection has A A^T=(c_RGS/c_Q6)I; "
                "therefore it is a similarity up to translation."
            ),
            "affine_equivalence_ruled_out": affine_equivalence_ruled_out,
            "reason": (
                "The normalized pair-distance multisets differ, while isotropic "
                "covariances force every affine equivalence to be a similarity."
                if affine_equivalence_ruled_out else
                "The present invariant test does not establish affine equivalence."
            ),
        },
        "epistemic_status": (
            "NUMERICALLY_VERIFIED_NO_GO_FOR_AFFINE_IDENTIFICATION"
            if affine_equivalence_ruled_out else
            "INCONCLUSIVE_REQUIRES_ADDITIONAL_INVARIANTS"
        ),
        "interpretation": (
            "The golden 6D projection is not an affine re-embedding of the "
            "archived RGS 64B/O_h point cloud under the tested equal-weight "
            "point-set identification. A non-affine map, a different acceptance "
            "window/seed, or a larger fibred correspondence remains open."
        ),
    }


def export_intertwiner_audit(path: str | Path) -> dict:
    """Write the reproducible JSON audit and return its report."""
    report = audit_rgs_golden_intertwiner()
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("results/rgs_golden_intertwiner.json"))
    args = parser.parse_args()
    print(json.dumps(export_intertwiner_audit(args.output), indent=2))

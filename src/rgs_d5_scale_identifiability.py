"""RGS-to-D5 scale-identifiability and rank-obstruction audit.

The module tests what follows from the checked-in 64x3 RGS seed and the
archived Gaussian-covariance coframe rule. It does not construct a D5 domain.
"""
from __future__ import annotations

import numpy as np

from src.rgs_ellis_refinement import load_rgs_seed_data

ELL2 = 7.035514834154695
C0 = 43.96495130838031
D5_REAL_DIMENSION = 10
SMILGA_TWO_BODY_SHELL_DIMENSION = 5


def local_covariance_coframe(
    coordinates: np.ndarray,
    ell2: float = ELL2,
    c0: float = C0,
) -> tuple[np.ndarray, np.ndarray]:
    """Reproduce the documented local covariance and SPD square-root coframe."""
    x = np.asarray(coordinates, dtype=np.float64)
    if x.shape != (64, 3):
        raise ValueError("coordinates must have shape (64,3)")
    if not np.isfinite(ell2) or ell2 <= 0 or not np.isfinite(c0) or c0 <= 0:
        raise ValueError("ell2 and c0 must be finite and positive")
    covariance = np.empty((64, 3, 3), dtype=np.float64)
    coframe = np.empty_like(covariance)
    for i in range(64):
        delta = x - x[i]
        weights = np.exp(-np.einsum("ij,ij->i", delta, delta) / ell2)
        weights[i] = 0.0
        ci = np.einsum("i,ia,ib->ab", weights, delta, delta)
        covariance[i] = ci / c0
        eig, vec = np.linalg.eigh(covariance[i])
        if eig[0] <= 0.0:
            raise ValueError(f"covariance is not positive definite at site {i}")
        coframe[i] = (vec * np.sqrt(eig)) @ vec.T
    return covariance, coframe


def q6_pullback_jacobians(coordinates: np.ndarray) -> np.ndarray:
    """Six Q6 finite-difference columns in R^3 at each of 64 vertices."""
    x = np.asarray(coordinates, dtype=np.float64)
    if x.shape != (64, 3):
        raise ValueError("coordinates must have shape (64,3)")
    jac = np.empty((64, 3, 6), dtype=np.float64)
    for v in range(64):
        for axis in range(6):
            bit = 1 << (5 - axis)
            low = v & ~bit
            high = low | bit
            jac[v, :, axis] = x[high] - x[low]
    return jac


def audit_rgs_d5_scale_identifiability(scale: float = 1.7) -> dict:
    """Audit scale covariance and the rank obstruction to a D5 pullback metric."""
    if not np.isfinite(scale) or scale <= 0.0:
        raise ValueError("scale must be finite and positive")
    seed = load_rgs_seed_data()
    x = np.asarray(seed.coordinates, dtype=np.float64)

    c, e = local_covariance_coframe(x)
    c_scaled_fixed, e_scaled_fixed = local_covariance_coframe(scale * x)
    c_scaled_covariant, e_scaled_covariant = local_covariance_coframe(
        scale * x, ell2=scale**2 * ELL2, c0=scale**2 * C0
    )

    jac = q6_pullback_jacobians(x)
    ranks = [int(np.linalg.matrix_rank(j)) for j in jac]
    gram = np.einsum("vki,vkj->vij", jac, jac)
    gram_ranks = [int(np.linalg.matrix_rank(g, tol=1e-10)) for g in gram]

    return {
        "audit": "RGS_D5_SCALE_IDENTIFIABILITY_V1",
        "seed_provenance": "load_rgs_seed_data() with archived SHA-256 checks",
        "scale_probe": float(scale),
        "coframe_rule": "E_i = principal_sqrt(C_i/c0), C_i = sum_j exp(-|X_j-X_i|^2/ell2) deltaX deltaX^T",
        "fixed_parameters": {
            "ell2": ELL2,
            "c0": C0,
            "covariance_scaling_residual": float(
                np.linalg.norm(c_scaled_fixed - scale**2 * c)
            ),
            "coframe_scaling_residual": float(
                np.linalg.norm(e_scaled_fixed - scale * e)
            ),
            "interpretation": "with ell2 fixed, the dimensionless Gaussian weights change under X->sX; no homogeneous scaling law is expected",
        },
        "covariant_parameters": {
            "ell2_scaled": scale**2 * ELL2,
            "c0_scaled": scale**2 * C0,
            "covariance_invariance_residual": float(
                np.linalg.norm(c_scaled_covariant - c)
            ),
            "coframe_invariance_residual": float(
                np.linalg.norm(e_scaled_covariant - e)
            ),
            "interpretation": "if ell2 and c0 carry squared-length units and scale with X, normalized covariance and coframe are invariant",
        },
        "q6_to_r3_rank": {
            "jacobian_shape": list(jac.shape),
            "rank_min": min(ranks),
            "rank_max": max(ranks),
            "gram_shape": list(gram.shape),
            "gram_rank_min": min(gram_ranks),
            "gram_rank_max": max(gram_ranks),
            "rank_bound": 3,
        },
        "d5_target": {
            "real_tangent_dimension": D5_REAL_DIMENSION,
            "positive_definite_pullback_possible_from_current_q6_to_r3_map": False,
            "minimum_missing_rank": D5_REAL_DIMENSION - 3,
            "smilga_two_body_kinematic_slice_dimension": SMILGA_TWO_BODY_SHELL_DIMENSION,
            "smilga_slice_decomposition": "3 total-momentum coordinates + 2 relative-direction coordinates (S^2)",
            "positive_definite_pullback_possible_from_current_q6_r3_map_to_smilga_slice": False,
            "minimum_missing_rank_to_smilga_slice": SMILGA_TWO_BODY_SHELL_DIMENSION - 3,
        },
        "checks": {
            "fixed_kernel_length_breaks_naive_covariance_scaling": bool(
                not np.allclose(c_scaled_fixed, scale**2*c, rtol=2e-12, atol=2e-12)
            ),
            "fixed_kernel_length_breaks_naive_coframe_scaling": bool(
                not np.allclose(e_scaled_fixed, scale*e, rtol=2e-12, atol=2e-12)
            ),
            "covariant_normalized_covariance_invariance": bool(
                np.allclose(c_scaled_covariant, c, rtol=2e-12, atol=2e-12)
            ),
            "covariant_coframe_invariance": bool(
                np.allclose(e_scaled_covariant, e, rtol=2e-12, atol=2e-12)
            ),
            "q6_jacobian_rank_at_most_three": bool(max(ranks) <= 3),
            "q6_gram_rank_at_most_three": bool(max(gram_ranks) <= 3),
            "rank_gap_to_full_d5_is_seven": D5_REAL_DIMENSION - 3 == 7,
            "rank_gap_to_smilga_two_body_slice_is_two": SMILGA_TWO_BODY_SHELL_DIMENSION - 3 == 2,
        },
        "verdict": "NO_GO_FOR_DIRECT_FULL_D5_OR_SMILGA_5D_SLICE_PULLBACK_FROM_CURRENT_R3_SEED; ADD_TWO_RELATIVE_DIRECTION_MODES_FOR_THE_SMILGA_SLICE; SCALE_NOT_PHYSICALLY_IDENTIFIED_UNLESS_ell2_AND_c0_HAVE_INDEPENDENT_UNIT_CALIBRATION",
    }


if __name__ == "__main__":
    import json
    print(json.dumps(audit_rgs_d5_scale_identifiability(), indent=2))

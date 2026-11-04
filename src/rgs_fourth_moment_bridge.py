"""Fourth-moment bridge for the recovered 64-site RGS/O_h seed.

The code computes the rank-four moment of the actual hash-verified 64x3
coordinate seed, projects it orthogonally onto the SO(3)-isotropic tensor,
and retains the trace-free cubic anisotropy. It also computes the induced
linear map Sym_0(3)->Sym_0(3), whose O_h splitting is Eg + T2g.

This is a finite-seed numerical audit, not a derivation of torsion or a
claim of continuum convergence.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from src.rgs_ellis_refinement import load_rgs_seed_data
from src.rgs_q6_oh_action import signed_permutation_matrices


def symmetric_traceless_basis() -> tuple[np.ndarray, ...]:
    """Frobenius-orthonormal basis of Sym_0(3): Eg (2) + T2g (3)."""
    root2 = np.sqrt(2.0)
    root6 = np.sqrt(6.0)
    basis = [
        np.diag([1.0, -1.0, 0.0]) / root2,
        np.diag([1.0, 1.0, -2.0]) / root6,
    ]
    for i, j in ((0, 1), (0, 2), (1, 2)):
        q = np.zeros((3, 3), dtype=float)
        q[i, j] = q[j, i] = 1.0 / root2
        basis.append(q)
    return tuple(basis)


def isotropic_rank4_tensor(dimension: int, coefficient: float) -> np.ndarray:
    """c*(delta_ij delta_kl + delta_ik delta_jl + delta_il delta_jk)."""
    eye = np.eye(int(dimension), dtype=float)
    return float(coefficient) * (
        np.einsum("ij,kl->ijkl", eye, eye)
        + np.einsum("ik,jl->ijkl", eye, eye)
        + np.einsum("il,jk->ijkl", eye, eye)
    )


def fourth_moment_decomposition(points: np.ndarray) -> dict:
    """Return M2, M4, its SO(d)-isotropic projection and harmonic residual.

    M4 is the unnormalised sum over sites. For d dimensions the orthogonal
    SO(d)-invariant projection coefficient is
        c = M_ i i j j / (d*(d+2)).
    The residual has vanishing double trace and is orthogonal (Frobenius)
    to the isotropic rank-four tensor.
    """
    x = np.asarray(points, dtype=float)
    if x.ndim != 2 or x.shape[0] == 0 or not np.all(np.isfinite(x)):
        raise ValueError("points must be a nonempty finite (N,d) array")
    n, d = x.shape
    m2 = np.einsum("ni,nj->ij", x, x)
    m4 = np.einsum("ni,nj,nk,nl->ijkl", x, x, x, x)
    double_trace = float(np.einsum("iijj", m4))
    coefficient = double_trace / float(d * (d + 2))
    m4_iso = isotropic_rank4_tensor(d, coefficient)
    residual = m4 - m4_iso
    return {
        "n_points": int(n),
        "dimension": int(d),
        "M2": m2,
        "M4": m4,
        "isotropic_coefficient": float(coefficient),
        "M4_isotropic": m4_iso,
        "M4_anisotropic": residual,
        "double_trace_M4": double_trace,
        "anisotropic_double_trace_norm": float(
            np.linalg.norm(np.einsum("iikl->kl", residual))
        ),
        "isotropic_orthogonality": float(np.einsum("ijkl,ijkl", residual, m4_iso)),
        "anisotropic_frobenius_norm": float(np.linalg.norm(residual)),
    }


def anisotropic_response_operator(residual: np.ndarray) -> dict:
    """Construct Q_kl -> R_ijkl Q_kl on Sym_0(3), and its O_h sectors."""
    r = np.asarray(residual, dtype=float)
    if r.shape != (3, 3, 3, 3):
        raise ValueError("residual must have shape (3,3,3,3)")
    basis = symmetric_traceless_basis()
    matrix = np.asarray([
        [float(np.einsum("ij,ijkl,kl", basis[a], r, basis[b]))
         for b in range(5)]
        for a in range(5)
    ])
    eig = np.linalg.eigvalsh(matrix)
    eg = matrix[:2, :2]
    t2g = matrix[2:, 2:]
    return {
        "basis_order": ["Eg_1", "Eg_2", "T2g_xy", "T2g_xz", "T2g_yz"],
        "matrix": matrix,
        "eigenvalues": eig,
        "Eg_block_offdiagonal_norm": float(np.linalg.norm(eg - np.diag(np.diag(eg)))),
        "T2g_block_offdiagonal_norm": float(np.linalg.norm(t2g - np.diag(np.diag(t2g)))),
        "Eg_T2g_mixing_norm": float(np.linalg.norm(matrix[:2, 2:])),
        "Eg_eigenvalues": np.linalg.eigvalsh(eg),
        "T2g_eigenvalues": np.linalg.eigvalsh(t2g),
    }


def audit_rgs_fourth_moment(tolerance: float = 1e-9) -> dict:
    """Audit the 64-site RGS fourth moment and cubic harmonic bridge."""
    seed = load_rgs_seed_data()
    x = seed.coordinates
    parts = fourth_moment_decomposition(x)
    m2, m4, r = parts["M2"], parts["M4"], parts["M4_anisotropic"]
    group = signed_permutation_matrices()

    m2_residual = float(np.linalg.norm(m2 - 80.0 * np.eye(3)))
    group_residuals = []
    for g in group:
        transformed = np.einsum(
            "ia,jb,kc,ld,abcd->ijkl", g, g, g, g, m4
        )
        group_residuals.append(float(np.linalg.norm(transformed - m4)))

    a = float(m4[0, 0, 0, 0])
    b = float(m4[0, 0, 1, 1])
    c = float(parts["isotropic_coefficient"])
    response = anisotropic_response_operator(r)
    expected_eg = a - 3.0 * c - (b - c)
    expected_t2g = 2.0 * (b - c)

    report = {
        "audit": "RGS_64_SITE_FOURTH_MOMENT_CUBIC_BRIDGE",
        "coordinate_source": "load_rgs_seed_data; archived SHA-256 hashes checked",
        "point_shape": list(x.shape),
        "M2": m2.tolist(),
        "M2_equals_80I_residual": m2_residual,
        "M4_components": {
            "A_xxxx": a,
            "B_xxyy": b,
            "isotropic_projection_coefficient_c": c,
            "SO3_isotropic_xxxx": 3.0 * c,
            "SO3_isotropic_xxyy": c,
            "harmonic_xxxx": float(r[0, 0, 0, 0]),
            "harmonic_xxyy": float(r[0, 0, 1, 1]),
        },
        "M4_anisotropic_frobenius_norm": parts["anisotropic_frobenius_norm"],
        "anisotropic_double_trace_norm": parts["anisotropic_double_trace_norm"],
        "isotropic_orthogonality": parts["isotropic_orthogonality"],
        "O_h": {
            "group_order": len(group),
            "M4_max_invariance_residual": max(group_residuals),
        },
        "Sym0_response": {
            "interpretation": "Q_kl -> M4_anisotropic_ijkl Q_kl",
            "basis_order": response["basis_order"],
            "matrix": response["matrix"].tolist(),
            "eigenvalues": response["eigenvalues"].tolist(),
            "Eg_eigenvalues": response["Eg_eigenvalues"].tolist(),
            "T2g_eigenvalues": response["T2g_eigenvalues"].tolist(),
            "Eg_T2g_mixing_norm": response["Eg_T2g_mixing_norm"],
            "expected_Eg_eigenvalue": float(expected_eg),
            "expected_T2g_eigenvalue": float(expected_t2g),
        },
        "status": {
            "M2_isotropy": "NUMERICALLY_VERIFIED",
            "SO3_projection_and_trace_free_rank4_residual": "NUMERICALLY_VERIFIED",
            "O_h_invariance": "NUMERICALLY_VERIFIED",
            "Eg_T2g_response_splitting": "NUMERICALLY_VERIFIED",
            "vierbein_constitutive_coupling": "CANDIDATE",
            "torsion_derivation": "OPEN",
            "continuum_limit_and_Mosco_convergence": "OPEN",
        },
    }
    report["pass"] = bool(
        m2_residual < tolerance
        and report["anisotropic_double_trace_norm"] < tolerance
        and abs(report["isotropic_orthogonality"]) < tolerance
        and report["O_h"]["M4_max_invariance_residual"] < tolerance
        and response["Eg_T2g_mixing_norm"] < tolerance
        and np.allclose(
            response["Eg_eigenvalues"], [expected_eg, expected_eg],
            atol=tolerance, rtol=0.0
        )
        and np.allclose(
            response["T2g_eigenvalues"], [expected_t2g] * 3,
            atol=tolerance, rtol=0.0
        )
    )
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", type=Path,
        default=Path("results/json/rgs_64_fourth_moment_bridge.json")
    )
    args = parser.parse_args()
    report = audit_rgs_fourth_moment()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if not report["pass"]:
        raise SystemExit("RGS fourth-moment bridge audit failed")


if __name__ == "__main__":
    main()

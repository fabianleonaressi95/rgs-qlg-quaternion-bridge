"""Deterministic reconstruction of the RGS spatial coframe and 4D tetrad candidates.

This script reconstructs the 64x3 RGS seed from the checked-in provenance
manifest, applies the archived local-covariance prescription, and exports
separate 3+1 and analytic Lorentzian tetrads.

Important distinctions:
- RGS_PRISM_FINAL_vierbein_E.npy is the 64-site spatial coframe E_i.
- RGS_PRISM_tetrad_4d.npy is its static block-diagonal 3+1 lift.
- Ellis_tetrad_4d.npy is a coordinate coframe for the Ellis wormhole metric.
- Kerr_Schild_tetrad_4d.npy is a Kerr-Schild Cartesian coframe.

The outputs are deterministic reconstructions from archived formulas, not
proof that the tetrad is a solution of a coupled RGS-EMS field equation.
No TEGR torsion scalar is certified: that requires a separately fixed
inertial spin connection and convention audit.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.rgs_ellis_refinement import load_rgs_seed_data

ELL2 = 7.035514834154695
C0 = 43.96495130838031
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])


def reconstruct_spatial_vierbein(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Archived local Gaussian covariance and its principal SPD square root."""
    x = np.asarray(x, dtype=np.float64)
    if x.shape != (64, 3):
        raise ValueError(f"expected X shape (64,3), got {x.shape}")
    n = len(x)
    covariance = np.empty((n, 3, 3), dtype=np.float64)
    coframe = np.empty_like(covariance)
    for i in range(n):
        delta = x - x[i]
        weights = np.exp(-np.einsum("ij,ij->i", delta, delta) / ELL2)
        weights[i] = 0.0
        ci = np.einsum("i,ia,ib->ab", weights, delta, delta)
        covariance[i] = ci / C0
        lam, u = np.linalg.eigh(covariance[i])
        if lam[0] <= 0.0:
            raise ValueError(f"local covariance is not positive definite at site {i}")
        # Principal symmetric square root: E_i = U sqrt(Lambda) U^T.
        coframe[i] = (u * np.sqrt(lam)) @ u.T
    return coframe, covariance


def static_3plus1_tetrad(spatial_coframe: np.ndarray) -> np.ndarray:
    """Static lift e^A_mu = diag(1,E_i), with eta=(-,+,+,+)."""
    e = np.zeros((len(spatial_coframe), 4, 4), dtype=np.float64)
    e[:, 0, 0] = 1.0
    e[:, 1:, 1:] = spatial_coframe
    return e


def ellis_diagonal_tetrad(
    u: np.ndarray, theta: np.ndarray, r0: float = 0.8292,
    lapse: np.ndarray | float = 1.0,
) -> np.ndarray:
    """Coordinate coframe for ds²=-N(u)²dt²+du²+R(u)²dΩ²."""
    u = np.asarray(u, dtype=np.float64).reshape(-1)
    theta = np.asarray(theta, dtype=np.float64).reshape(-1)
    if u.shape != theta.shape:
        raise ValueError("u and theta must have the same shape")
    if r0 <= 0.0 or np.any((theta < 0.0) | (theta > np.pi)):
        raise ValueError("r0 must be positive and theta in [0,pi]")
    n = np.broadcast_to(np.asarray(lapse, dtype=np.float64), u.shape)
    if np.any(n <= 0.0):
        raise ValueError("lapse must be positive")
    radius = np.sqrt(u*u + r0*r0)
    e = np.zeros((len(u), 4, 4), dtype=np.float64)
    e[:, 0, 0] = n
    e[:, 1, 1] = 1.0
    e[:, 2, 2] = radius
    e[:, 3, 3] = radius * np.sin(theta)
    return e


def kerr_schild_tetrad(
    xyz: np.ndarray, mass: float = 1.0, spin: float = 0.6,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Kerr-Schild coframe e^a_mu=delta^a_mu+H l^a l_mu in Cartesian KS coords.

    Signature is (-,+,+,+), l_mu dx^mu = dt + l_i dx^i, and
    H=M*r/(r^2+a^2*cos(theta)^2). The oblate radius r is the positive
    solution of (x²+y²)/(r²+a²)+z²/r²=1.
    """
    xyz = np.asarray(xyz, dtype=np.float64)
    if xyz.ndim != 2 or xyz.shape[1] != 3:
        raise ValueError("xyz must have shape (N,3)")
    if mass <= 0.0 or abs(spin) > mass:
        raise ValueError("require mass>0 and |spin|<=mass")
    x, y, z = xyz.T
    rho2 = x*x + y*y + z*z
    a2 = spin*spin
    disc = (rho2 - a2)**2 + 4.0*a2*z*z
    r2 = 0.5 * (rho2 - a2 + np.sqrt(np.maximum(disc, 0.0)))
    r = np.sqrt(np.maximum(r2, 0.0))
    if np.any(r <= 1e-12):
        raise ValueError("Kerr-Schild chart is singular at the ring/coordinate singularity")
    cos_theta = z / r
    den = r2 + a2
    lx = (r*x + spin*y) / den
    ly = (r*y - spin*x) / den
    lz = z / r
    l_cov = np.column_stack((np.ones_like(r), lx, ly, lz))
    l_con = l_cov @ ETA
    H = mass*r / (r2 + a2*cos_theta*cos_theta)
    tetrad = np.broadcast_to(np.eye(4), (len(r), 4, 4)).copy()
    tetrad += H[:, None, None] * l_con[:, :, None] * l_cov[:, None, :]
    metric = ETA[None, :, :] + 2.0*H[:, None, None] * l_cov[:, :, None] * l_cov[:, None, :]
    return tetrad, metric, np.column_stack((r, H, cos_theta))


def metric_from_coframe(e: np.ndarray) -> np.ndarray:
    return np.einsum("...am,ab,...bn->...mn", e, ETA, e)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("results/vierbein_reconstruction"))
    args = parser.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)

    seed = load_rgs_seed_data()  # checks archived SHA-256 values internally
    X = seed.coordinates
    E, C = reconstruct_spatial_vierbein(X)
    e_static = static_3plus1_tetrad(E)

    # Ellis chart sample: 64 u-values, 64 polar angles; phi is a coordinate
    # label only because the diagonal coframe is independent of phi.
    u = np.linspace(-3.0*0.8292, 3.0*0.8292, 64)
    theta = np.linspace(0.15, np.pi-0.15, 64)
    e_ellis = ellis_diagonal_tetrad(u, theta, r0=0.8292, lapse=1.0)

    # Kerr sample away from the ring singularity; radii span both sides of r+.
    radii = np.linspace(1.95, 4.0, 64)
    ang = np.linspace(0.2, np.pi-0.2, 64)
    phi = np.linspace(0.0, 2.0*np.pi, 64, endpoint=False)
    xyz = np.column_stack((
        radii*np.sin(ang)*np.cos(phi),
        radii*np.sin(ang)*np.sin(phi),
        radii*np.cos(ang),
    ))
    e_kerr, g_kerr_exact, kerr_aux = kerr_schild_tetrad(xyz, mass=1.0, spin=0.6)

    g_static = metric_from_coframe(e_static)
    g_ellis = metric_from_coframe(e_ellis)
    g_kerr = metric_from_coframe(e_kerr)
    report = {
        "status": "DETERMINISTIC_ANALYTIC_RECONSTRUCTION",
        "provenance": {
            "seed": "src.rgs_ellis_refinement.load_rgs_seed_data",
            "seed_sha256_checked": True,
            "ell_squared": ELL2,
            "c0": C0,
            "spatial_rule": "E_i = principal_sqrt(C_i/c0)",
            "C_i_rule": "sum_j exp(-|X_j-X_i|^2/ell2) (X_j-X_i)(X_j-X_i)^T",
        },
        "outputs": {
            "X_Oh_reconstructed.npy": list(X.shape),
            "RGS_PRISM_FINAL_vierbein_E.npy": list(E.shape),
            "RGS_PRISM_tetrad_4d.npy": list(e_static.shape),
            "Ellis_tetrad_4d.npy": list(e_ellis.shape),
            "Kerr_Schild_tetrad_4d.npy": list(e_kerr.shape),
        },
        "spatial_vierbein": {
            "det_min": float(np.linalg.det(E).min()),
            "det_max": float(np.linalg.det(E).max()),
            "singular_min": float(min(np.linalg.svd(v, compute_uv=False).min() for v in E)),
            "condition_max": float(max(np.linalg.cond(v) for v in E)),
            "sqrt_residual": float(np.linalg.norm(E@np.swapaxes(E,1,2)-C)),
        },
        "static_3plus1": {
            "signature": "(-,+,+,+)",
            "metric_reconstruction_residual": float(np.linalg.norm(g_static[:,1:,1:]-np.einsum("...ki,...kj->...ij",E,E))),
            "time_leg": "e^0=dt; static adapter, not a derived dynamical lapse/shift",
        },
        "ellis": {
            "metric": "ds^2=-dt^2+du^2+(u^2+r0^2)(dtheta^2+sin(theta)^2 dphi^2)",
            "r0": 0.8292,
            "metric_reconstruction_residual": float(np.linalg.norm(g_ellis-np.einsum("...am,ab,...bn->...mn",e_ellis,ETA,e_ellis))),
            "det_min": float(np.linalg.det(e_ellis).min()),
        },
        "kerr_schild": {
            "mass": 1.0, "spin": 0.6,
            "horizons": [1.8, 0.2],
            "metric_reconstruction_residual": float(np.linalg.norm(g_kerr-g_kerr_exact)),
            "determinant_residual_max": float(np.max(np.abs(np.linalg.det(e_kerr)-1.0))),
            "det_min": float(np.linalg.det(e_kerr).min()),
            "aux_columns": ["oblate_r", "H", "cos_theta"],
            "aux": kerr_aux.tolist(),
        },
        "teleparallel": {
            "status": "NOT_CERTIFIED",
            "reason": "torsion and TEGR identity require an explicitly fixed inertial spin connection and convention audit",
        },
    }

    np.save(out/"X_Oh_reconstructed.npy", X)
    np.save(out/"RGS_PRISM_FINAL_vierbein_E.npy", E)
    np.save(out/"RGS_PRISM_tetrad_4d.npy", e_static)
    np.save(out/"Ellis_tetrad_4d.npy", e_ellis)
    np.save(out/"Kerr_Schild_tetrad_4d.npy", e_kerr)
    np.save(out/"Kerr_Schild_aux_r_H_costheta.npy", kerr_aux)
    (out/"vierbein_reconstruction_report.json").write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

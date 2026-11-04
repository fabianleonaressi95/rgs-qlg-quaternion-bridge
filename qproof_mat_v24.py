#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MCRG SU(2) 2D — Numba Optimized Pipeline
=========================================
"""

from __future__ import annotations

import json
import csv
import math
import time
import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from numba import jit, prange


# ============================================================
# CONSTANTS
# ============================================================

PHI = (1.0 + np.sqrt(5.0)) / 2.0


# ============================================================
# CONFIGURATION
# ============================================================

@dataclass
class SimulationConfig:
    L: int = 16
    beta: float = 2.3

    therm_sweeps: int = 300
    meas_sweeps: int = 600
    sweep_between_measurements: int = 2

    proposal_epsilon: float = 0.25
    blocking_factor: int = 2
    bootstrap_samples: int = 300
    seed: int = 12345
    output_dir: str = "mcrg_results"
    verbose: bool = True


# ============================================================
# SU(2) KERNELS (NUMBA JIT)
# ============================================================

@jit(nopython=True)
def su2_identity() -> np.ndarray:
    return np.array([1.0, 0.0, 0.0, 0.0], dtype=np.float64)


@jit(nopython=True)
def su2_dagger(q: np.ndarray) -> np.ndarray:
    return np.array([q[0], -q[1], -q[2], -q[3]], dtype=np.float64)


@jit(nopython=True)
def su2_mul(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    a0, a1, a2, a3 = a
    b0, b1, b2, b3 = b
    return np.array([
        a0*b0 - a1*b1 - a2*b2 - a3*b3,
        a0*b1 + a1*b0 + a2*b3 - a3*b2,
        a0*b2 - a1*b3 + a2*b0 + a3*b1,
        a0*b3 + a1*b2 - a2*b1 + a3*b0
    ], dtype=np.float64)


@jit(nopython=True)
def su2_normalize(q: np.ndarray) -> np.ndarray:
    n = np.sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3])
    if n < 1e-15:
        return su2_identity()
    return q / n


@jit(nopython=True)
def su2_retrace(q: np.ndarray) -> float:
    val = q[0]
    if val > 1.0:
        return 1.0
    if val < -1.0:
        return -1.0
    return float(val)


def su2_random_quaternion(rng: np.random.Generator) -> np.ndarray:
    x = rng.normal(size=4)
    x /= np.linalg.norm(x)
    return x


def su2_random_near_identity(rng: np.random.Generator, epsilon: float) -> np.ndarray:
    v = rng.normal(size=3)
    n = np.linalg.norm(v)
    if n < 1e-15:
        return su2_identity()
    v /= n
    angle = epsilon * rng.normal()
    c = np.cos(angle)
    s = np.sin(angle)
    return np.array([c, s * v[0], s * v[1], s * v[2]], dtype=np.float64)


# ============================================================
# LATTICE & LOW-LEVEL KERNELS
# ============================================================

@dataclass
class Lattice:
    L: int
    links: np.ndarray


def create_lattice(L: int) -> Lattice:
    links = np.zeros((L, L, 2, 4), dtype=np.float64)
    links[...] = su2_identity()
    return Lattice(L=L, links=links)


@jit(nopython=True)
def get_link_jit(links: np.ndarray, L: int, x: int, y: int, mu: int) -> np.ndarray:
    return links[x % L, y % L, mu]


@jit(nopython=True)
def staple_jit(links: np.ndarray, L: int, x: int, y: int, mu: int) -> np.ndarray:
    if mu == 0:
        # Link a
        l1 = links[x % L, y % L, 1]
        l2 = links[(x) % L, (y + 1) % L, 0]
        l3 = links[(x + 1) % L, y % L, 1]
        
        a = su2_mul(l1, l2)
        a = su2_mul(a, su2_dagger(l3))

        # Link b
        y_m1 = (y - 1) % L
        l4 = links[x % L, y_m1, 1]
        l5 = links[x % L, y_m1, 0]
        l6 = links[(x + 1) % L, y_m1, 1]

        b = su2_mul(su2_dagger(l4), l5)
        b = su2_mul(b, l6)
    else:
        # Link a
        l1 = links[x % L, y % L, 0]
        l2 = links[(x + 1) % L, y % L, 1]
        l3 = links[x % L, (y + 1) % L, 0]

        a = su2_mul(l1, l2)
        a = su2_mul(a, su2_dagger(l3))

        # Link b
        x_m1 = (x - 1) % L
        l4 = links[x_m1, y % L, 0]
        l5 = links[x_m1, y % L, 1]
        l6 = links[x_m1, (y + 1) % L, 0]

        b = su2_mul(su2_dagger(l4), l5)
        b = su2_mul(b, l6)

    return a + b


@jit(nopython=True, cache=True)
def metropolis_sweep_jit(links: np.ndarray, L: int, beta: float, epsilon: float, rand_vals: np.ndarray, prop_quats: np.ndarray) -> float:
    accepted = 0
    total = 0
    idx_r = 0
    idx_p = 0

    for x in range(L):
        for y in range(L):
            for mu in range(2):
                old = links[x, y, mu].copy()
                proposal = prop_quats[idx_p]
                idx_p += 1
                
                # New link = proposal * old
                a0, a1, a2, a3 = proposal
                b0, b1, b2, b3 = old
                new = np.array([
                    a0*b0 - a1*b1 - a2*b2 - a3*b3,
                    a0*b1 + a1*b0 + a2*b3 - a3*b2,
                    a0*b2 - a1*b3 + a2*b0 + a3*b1,
                    a0*b3 + a1*b2 - a2*b1 + a3*b0
                ], dtype=np.float64)
                
                # Normalize new
                n = np.sqrt(new[0]*new[0] + new[1]*new[1] + new[2]*new[2] + new[3]*new[3])
                if n >= 1e-15:
                    new /= n

                st = staple_jit(links, L, x, y, mu)

                # Local action change calculation
                old_val = old[0] * st[0] - old[1] * st[1] - old[2] * st[2] - old[3] * st[3]
                if old_val > 1.0: old_val = 1.0
                elif old_val < -1.0: old_val = -1.0

                new_val = new[0] * st[0] - new[1] * st[1] - new[2] * st[2] - new[3] * st[3]
                if new_val > 1.0: new_val = 1.0
                elif new_val < -1.0: new_val = -1.0

                dS = -beta * (new_val - old_val)

                if dS <= 0.0 or rand_vals[idx_r] < np.exp(-dS):
                    links[x, y, mu] = new
                    accepted += 1
                
                idx_r += 1
                total += 1

    return float(accepted) / float(total)


@jit(nopython=True, cache=True)
def loop_product_jit(links: np.ndarray, L: int, x: int, y: int, width: int, height: int) -> np.ndarray:
    U = su2_identity()

    for i in range(width):
        U = su2_mul(U, links[(x + i) % L, y % L, 0])
    for j in range(height):
        U = su2_mul(U, links[(x + width) % L, (y + j) % L, 1])
    for i in range(width):
        dag = su2_dagger(links[(x + width - 1 - i) % L, (y + height) % L, 0])
        U = su2_mul(U, dag)
    for j in range(height):
        dag = su2_dagger(links[x % L, (y + height - 1 - j) % L, 1])
        U = su2_mul(U, dag)

    return U


@jit(nopython=True, cache=True)
def wilson_loop_jit(links: np.ndarray, L: int, width: int, height: int) -> float:
    total = 0.0
    count = 0
    for x in range(L):
        for y in range(L):
            W = loop_product_jit(links, L, x, y, width, height)
            total += su2_retrace(W)
            count += 1
    return float(total / count)


@jit(nopython=True, cache=True)
def block_lattice_jit(links: np.ndarray, L: int, b: int) -> np.ndarray:
    Lc = L // b
    coarse_links = np.zeros((Lc, Lc, 2, 4), dtype=np.float64)

    for X in range(Lc):
        for Y in range(Lc):
            x = b * X
            y = b * Y

            # coarse x-link
            U_x = su2_identity()
            for k in range(b):
                U_x = su2_mul(U_x, links[(x + k) % L, y % L, 0])
            coarse_links[X, Y, 0] = su2_normalize(U_x)

            # coarse y-link
            U_y = su2_identity()
            for k in range(b):
                U_y = su2_mul(U_y, links[x % L, (y + k) % L, 1])
            coarse_links[X, Y, 1] = su2_normalize(U_y)

    return coarse_links


# ============================================================
# OPERATOR BASIS & MEASUREMENTS
# ============================================================

BASES = {
    "B1": [(1, 1), (1, 2), (2, 2)],
    "B2": [(1, 1), (1, 2), (1, 3), (2, 2), (2, 3)],
    "B3": [(1, 1), (1, 2), (1, 3), (1, 4), (2, 2), (2, 3), (2, 4), (3, 3)],
}


def measure_operator_basis(links: np.ndarray, L: int, basis) -> np.ndarray:
    return np.array([wilson_loop_jit(links, L, w, h) for w, h in basis], dtype=np.float64)


def generate_ensemble(cfg: SimulationConfig, basis, rng: np.random.Generator):
    lat = create_lattice(cfg.L)
    acceptance_history = []

    if cfg.verbose:
        print(f"[thermalization] L={cfg.L}, beta={cfg.beta}")

    total_links_per_sweep = cfg.L * cfg.L * 2

    for _ in range(cfg.therm_sweeps):
        rand_vals = rng.random(total_links_per_sweep)
        prop_quats = np.array([su2_random_near_identity(rng, cfg.proposal_epsilon) for _ in range(total_links_per_sweep)])
        acc = metropolis_sweep_jit(lat.links, cfg.L, cfg.beta, cfg.proposal_epsilon, rand_vals, prop_quats)
        acceptance_history.append(acc)

    fine_samples = []
    coarse_samples = []
    b = cfg.blocking_factor

    if cfg.verbose:
        print(f"[measurement] {cfg.meas_sweeps} measurements")

    for meas in range(cfg.meas_sweeps):
        for _ in range(cfg.sweep_between_measurements):
            rand_vals = rng.random(total_links_per_sweep)
            prop_quats = np.array([su2_random_near_identity(rng, cfg.proposal_epsilon) for _ in range(total_links_per_sweep)])
            acc = metropolis_sweep_jit(lat.links, cfg.L, cfg.beta, cfg.proposal_epsilon, rand_vals, prop_quats)
            acceptance_history.append(acc)

        fine_obs = measure_operator_basis(lat.links, cfg.L, basis)
        coarse_links = block_lattice_jit(lat.links, cfg.L, b=b)
        coarse_obs = measure_operator_basis(coarse_links, cfg.L // b, basis)

        fine_samples.append(fine_obs)
        coarse_samples.append(coarse_obs)

        if cfg.verbose and (meas + 1) % max(1, cfg.meas_sweeps // 10) == 0:
            print(f"  measurement {meas+1}/{cfg.meas_sweeps}")

    return np.asarray(fine_samples), np.asarray(coarse_samples), np.asarray(acceptance_history)


# ============================================================
# MCRG MAP & SPECTRAL ANALYSIS
# ============================================================

def estimate_rg_jacobian(fine_samples: np.ndarray, coarse_samples: np.ndarray) -> np.ndarray:
    Xc = fine_samples - fine_samples.mean(axis=0)
    Yc = coarse_samples - coarse_samples.mean(axis=0)
    A, _, _, _ = np.linalg.lstsq(Xc, Yc, rcond=None)
    return A.T


def spectral_analysis(J: np.ndarray):
    eigenvalues = np.linalg.eigvals(J)
    order = np.argsort(np.abs(eigenvalues))[::-1]
    eigenvalues = eigenvalues[order]

    rho1 = float(np.abs(eigenvalues[0]))
    rho2 = float(np.abs(eigenvalues[1])) if len(eigenvalues) > 1 else np.nan
    n_phi = -math.log(rho1) / math.log(PHI) if rho1 > 0 else np.nan

    return {
        "eigenvalues_real": [float(np.real(x)) for x in eigenvalues],
        "eigenvalues_imag": [float(np.imag(x)) for x in eigenvalues],
        "rho1": rho1,
        "rho2": rho2,
        "n_phi": n_phi,
    }


def bootstrap_jacobian(fine_samples: np.ndarray, coarse_samples: np.ndarray, n_bootstrap: int, rng: np.random.Generator):
    N = len(fine_samples)
    rho_values = []
    nphi_values = []

    for _ in range(n_bootstrap):
        indices = rng.integers(0, N, size=N)
        X = fine_samples[indices]
        Y = coarse_samples[indices]
        J = estimate_rg_jacobian(X, Y)
        spec = spectral_analysis(J)
        rho_values.append(spec["rho1"])
        nphi_values.append(spec["n_phi"])

    rho_values = np.asarray(rho_values)
    nphi_values = np.asarray(nphi_values)

    return {
        "rho_mean": float(np.mean(rho_values)),
        "rho_std": float(np.std(rho_values, ddof=1)),
        "rho_q025": float(np.quantile(rho_values, 0.025)),
        "rho_q975": float(np.quantile(rho_values, 0.975)),
        "nphi_mean": float(np.mean(nphi_values)),
        "nphi_std": float(np.std(nphi_values, ddof=1)),
        "nphi_q025": float(np.quantile(nphi_values, 0.025)),
        "nphi_q975": float(np.quantile(nphi_values, 0.975)),
    }


def phi_diagnostics(rho: float):
    targets = {
        "phi^-1": PHI ** (-1),
        "phi^-2": PHI ** (-2),
        "phi^-3": PHI ** (-3),
        "phi^-4": PHI ** (-4),
    }
    result = {}
    for name, target in targets.items():
        result[name] = {
            "target": target,
            "absolute_error": abs(rho - target),
            "relative_error": abs(rho - target) / target,
            "log_error": abs(math.log(rho / target)) if rho > 0 else np.inf
        }
    result["closest"] = min(result, key=lambda k: result[k]["log_error"])
    return result


# ============================================================
# WRITERS & RUNNER
# ============================================================

def save_json(path, data):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)


def save_csv(path, rows):
    if not rows:
        return
    keys = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def run_single(cfg: SimulationConfig, basis_name: str, rng: np.random.Generator):
    basis = BASES[basis_name]
    print("\n" + "=" * 72)
    print(f"MCRG SU(2) Numba | L={cfg.L} | beta={cfg.beta} | b={cfg.blocking_factor} | basis={basis_name}")
    print("=" * 72)

    t0 = time.time()
    fine, coarse, acceptance = generate_ensemble(cfg, basis, rng)
    J = estimate_rg_jacobian(fine, coarse)
    spec = spectral_analysis(J)
    boot = bootstrap_jacobian(fine, coarse, cfg.bootstrap_samples, rng)
    phi = phi_diagnostics(spec["rho1"])
    elapsed = time.time() - t0

    result = {
        "L": cfg.L,
        "beta": cfg.beta,
        "blocking_factor": cfg.blocking_factor,
        "basis": basis_name,
        "N_op": len(basis),
        "N_samples": len(fine),
        "acceptance": float(np.mean(acceptance)),
        "rho1": spec["rho1"],
        "rho2": spec["rho2"],
        "n_phi": spec["n_phi"],
        "bootstrap_rho_mean": boot["rho_mean"],
        "bootstrap_rho_std": boot["rho_std"],
        "bootstrap_rho_q025": boot["rho_q025"],
        "bootstrap_rho_q975": boot["rho_q975"],
        "bootstrap_nphi_mean": boot["nphi_mean"],
        "bootstrap_nphi_std": boot["nphi_std"],
        "bootstrap_nphi_q025": boot["nphi_q025"],
        "bootstrap_nphi_q975": boot["nphi_q975"],
        "phi_closest": phi["closest"],
        "phi_minus_1_error": phi["phi^-1"]["relative_error"],
        "phi_minus_2_error": phi["phi^-2"]["relative_error"],
        "phi_minus_3_error": phi["phi^-3"]["relative_error"],
        "phi_minus_4_error": phi["phi^-4"]["relative_error"],
        "elapsed_seconds": elapsed
    }
    return result, J


def main():
    parser = argparse.ArgumentParser(description="MCRG SU(2) 2D Numba pipeline")
    parser.add_argument("--L", type=int, nargs="+", default=[8, 16])
    parser.add_argument("--beta", type=float, nargs="+", default=[2.0, 2.3])
    parser.add_argument("--bases", type=str, nargs="+", default=["B1", "B2"])
    parser.add_argument("--blocking", type=int, nargs="+", default=[2])
    parser.add_argument("--therm", type=int, default=100)
    parser.add_argument("--measurements", type=int, default=200)
    parser.add_argument("--bootstrap", type=int, default=100)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--output", type=str, default="mcrg_results")
    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    all_results = []
    master_rng = np.random.default_rng(args.seed)

    for L in args.L:
        for beta in args.beta:
            for b in args.blocking:
                if L % b != 0:
                    continue
                Lc = L // b
                for basis_name in args.bases:
                    basis = BASES[basis_name]
                    max_loop = max(max(w, h) for w, h in basis)
                    if Lc < 2 * max_loop:
                        continue

                    seed = int(master_rng.integers(0, 2**32 - 1))
                    rng = np.random.default_rng(seed)

                    cfg = SimulationConfig(
                        L=L, beta=beta,
                        therm_sweeps=args.therm,
                        meas_sweeps=args.measurements,
                        blocking_factor=b,
                        bootstrap_samples=args.bootstrap,
                        seed=seed,
                        output_dir=str(output_dir)
                    )

                    result, J = run_single(cfg, basis_name, rng)
                    result["seed"] = seed
                    all_results.append(result)

                    J_name = f"J_L{L}_b{b}_beta{beta:.3f}_{basis_name}.npy"
                    np.save(output_dir / J_name, J)

    save_csv(output_dir / "mcrg_results.csv", all_results)
    save_json(output_dir / "mcrg_results.json", {"phi": float(PHI), "results": all_results})


if __name__ == "__main__":
    main()

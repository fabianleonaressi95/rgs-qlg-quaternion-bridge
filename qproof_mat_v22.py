#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MCRG SU(2) 2D — First Complete Pipeline
========================================

Obiettivo
---------
Prima pipeline numerica per studiare il flusso RG della teoria di gauge
SU(2) 2D con azione di Wilson.

La pipeline:

1. genera configurazioni SU(2) tramite Metropolis;
2. misura osservabili Wilson-loop;
3. applica blocking geometrico b=2;
4. misura nuovamente le osservabili coarse;
5. stima la mappa RG nello spazio degli operatori;
6. calcola il Jacobiano D R_b;
7. calcola gli autovalori del Jacobiano;
8. usa bootstrap per stimare le incertezze;
9. ripete per diverse dimensioni L;
10. ripete per diverse basi di operatori;
11. calcola n_Phi = -log(rho)/log(Phi) SOLO A POSTERIORI;
12. salva tutti i risultati in CSV e JSON.

IMPORTANTE
-----------
Questa è una pipeline MCRG numerica esplorativa.
NON costituisce una dimostrazione del mass gap di Yang-Mills 4D.
NON assume che Phi sia una costante RG.
NON usa rho_infinity come parametro imposto.

Dipendenze:
    numpy
    scipy

Esempio:
    python mcrg_su2_pipeline.py
"""

from __future__ import annotations

import json
import csv
import math
import time
import argparse
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


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
# SU(2) REPRESENTATION
# ============================================================

def su2_random_quaternion(rng: np.random.Generator) -> np.ndarray:
    """
    Haar-random SU(2) element represented as

        U = a0 I + i(a1 sigma_x + a2 sigma_y + a3 sigma_z)

    with a0^2 + a1^2 + a2^2 + a3^2 = 1.
    """

    x = rng.normal(size=4)
    x /= np.linalg.norm(x)
    return x


def su2_identity() -> np.ndarray:
    """Quaternion representation of identity."""
    return np.array([1.0, 0.0, 0.0, 0.0])


def su2_dagger(q: np.ndarray) -> np.ndarray:
    """Hermitian conjugate."""
    return np.array([q[0], -q[1], -q[2], -q[3]])


def su2_mul(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """
    Quaternion multiplication corresponding to SU(2) matrix multiplication.
    """

    a0, a1, a2, a3 = a
    b0, b1, b2, b3 = b

    return np.array([
        a0*b0 - a1*b1 - a2*b2 - a3*b3,

        a0*b1 + a1*b0 + a2*b3 - a3*b2,

        a0*b2 - a1*b3 + a2*b0 + a3*b1,

        a0*b3 + a1*b2 - a2*b1 + a3*b0
    ])


def su2_normalize(q: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(q)

    if n < 1e-15:
        return su2_identity()

    return q / n


def su2_retrace(q: np.ndarray) -> float:
    """
    For SU(2):

        Re Tr(U) / 2 = q0
    """

    return float(np.clip(q[0], -1.0, 1.0))


def su2_random_near_identity(
    rng: np.random.Generator,
    epsilon: float
) -> np.ndarray:
    """
    Small random SU(2) rotation around identity.
    """

    v = rng.normal(size=3)

    n = np.linalg.norm(v)

    if n < 1e-15:
        return su2_identity()

    v /= n

    angle = epsilon * rng.normal()

    return np.array([
        np.cos(angle),
        np.sin(angle) * v[0],
        np.sin(angle) * v[1],
        np.sin(angle) * v[2]
    ])


# ============================================================
# LATTICE
# ============================================================

@dataclass
class Lattice:
    L: int

    # links[x,y,mu,:]
    # mu=0 -> x direction
    # mu=1 -> y direction
    links: np.ndarray


def create_lattice(L: int) -> Lattice:

    links = np.zeros((L, L, 2, 4), dtype=np.float64)

    I = su2_identity()

    links[...] = I

    return Lattice(L=L, links=links)


# ============================================================
# LINK ACCESS
# ============================================================

def get_link(lat: Lattice, x: int, y: int, mu: int) -> np.ndarray:

    L = lat.L

    x %= L
    y %= L

    return lat.links[x, y, mu]


def set_link(
    lat: Lattice,
    x: int,
    y: int,
    mu: int,
    value: np.ndarray
):

    L = lat.L

    x %= L
    y %= L

    lat.links[x, y, mu] = su2_normalize(value)


def backward_link(
    lat: Lattice,
    x: int,
    y: int,
    mu: int
) -> np.ndarray:

    if mu == 0:
        return su2_dagger(get_link(lat, x - 1, y, 0))

    return su2_dagger(get_link(lat, x, y - 1, 1))


# ============================================================
# PLAQUETTE
# ============================================================

def plaquette(
    lat: Lattice,
    x: int,
    y: int
) -> np.ndarray:

    Ux = get_link(lat, x, y, 0)
    Uy = get_link(lat, x + 1, y, 1)

    Ux_d = su2_dagger(get_link(lat, x, y + 1, 0))
    Uy_d = su2_dagger(get_link(lat, x, y, 1))

    return su2_mul(
        su2_mul(
            su2_mul(Ux, Uy),
            Ux_d
        ),
        Uy_d
    )


def average_plaquette(lat: Lattice) -> float:

    L = lat.L

    values = []

    for x in range(L):
        for y in range(L):
            values.append(
                su2_retrace(
                    plaquette(lat, x, y)
                )
            )

    return float(np.mean(values))


# ============================================================
# STAPLES
# ============================================================

def staple(
    lat: Lattice,
    x: int,
    y: int,
    mu: int
) -> np.ndarray:

    """
    Sum of the two staples touching the link U_mu(x).
    """

    if mu == 0:

        # +y plaquette
        a = su2_mul(
            get_link(lat, x, y, 1),
            get_link(lat, x, y + 1, 0)
        )

        a = su2_mul(
            a,
            su2_dagger(
                get_link(lat, x + 1, y, 1)
            )
        )

        # -y plaquette
        b = su2_mul(
            su2_dagger(
                get_link(lat, x, y - 1, 1)
            ),
            get_link(lat, x, y - 1, 0)
        )

        b = su2_mul(
            b,
            get_link(lat, x + 1, y - 1, 1)
        )

    else:

        # +x plaquette
        a = su2_mul(
            get_link(lat, x, y, 0),
            get_link(lat, x + 1, y, 1)
        )

        a = su2_mul(
            a,
            su2_dagger(
                get_link(lat, x, y + 1, 0)
            )
        )

        # -x plaquette
        b = su2_mul(
            su2_dagger(
                get_link(lat, x - 1, y, 0)
            ),
            get_link(lat, x - 1, y, 1)
        )

        b = su2_mul(
            b,
            get_link(lat, x - 1, y + 1, 0)
        )

    return a + b


# ============================================================
# METROPOLIS
# ============================================================

def local_action_change(
    old_link: np.ndarray,
    new_link: np.ndarray,
    st: np.ndarray,
    beta: float
) -> float:

    old_val = su2_retrace(
        su2_mul(old_link, st)
    )

    new_val = su2_retrace(
        su2_mul(new_link, st)
    )

    return -beta * (new_val - old_val)


def metropolis_update(
    lat: Lattice,
    beta: float,
    rng: np.random.Generator,
    epsilon: float
) -> float:

    L = lat.L

    accepted = 0
    total = 0

    for x in range(L):
        for y in range(L):

            for mu in range(2):

                old = get_link(lat, x, y, mu).copy()

                proposal = su2_random_near_identity(
                    rng,
                    epsilon
                )

                new = su2_mul(proposal, old)

                st = staple(
                    lat,
                    x,
                    y,
                    mu
                )

                dS = local_action_change(
                    old,
                    new,
                    st,
                    beta
                )

                if dS <= 0.0 or rng.random() < np.exp(-dS):

                    set_link(
                        lat,
                        x,
                        y,
                        mu,
                        new
                    )

                    accepted += 1

                total += 1

    return accepted / total


# ============================================================
# WILSON LOOPS
# ============================================================

def loop_product(
    lat: Lattice,
    x: int,
    y: int,
    width: int,
    height: int
) -> np.ndarray:

    U = su2_identity()

    # +x
    for i in range(width):

        U = su2_mul(
            U,
            get_link(lat, x + i, y, 0)
        )

    # +y
    for j in range(height):

        U = su2_mul(
            U,
            get_link(
                lat,
                x + width,
                y + j,
                1
            )
        )

    # -x
    for i in range(width):

        U = su2_mul(
            U,
            su2_dagger(
                get_link(
                    lat,
                    x + width - 1 - i,
                    y + height,
                    0
                )
            )
        )

    # -y
    for j in range(height):

        U = su2_mul(
            U,
            su2_dagger(
                get_link(
                    lat,
                    x,
                    y + height - 1 - j,
                    1
                )
            )
        )

    return U


def wilson_loop(
    lat: Lattice,
    width: int,
    height: int
) -> float:

    L = lat.L

    values = []

    for x in range(L):
        for y in range(L):

            W = loop_product(
                lat,
                x,
                y,
                width,
                height
            )

            values.append(
                su2_retrace(W)
            )

    return float(np.mean(values))


# ============================================================
# OPERATOR BASIS
# ============================================================

BASES = {

    "B1": [
        (1, 1),   # plaquette
        (1, 2),   # rectangle
        (2, 2),   # square
    ],

    "B2": [
        (1, 1),
        (1, 2),
        (1, 3),
        (2, 2),
        (2, 3),
    ],

    "B3": [
        (1, 1),
        (1, 2),
        (1, 3),
        (1, 4),
        (2, 2),
        (2, 3),
        (2, 4),
        (3, 3),
    ],

    "B4": [
        (1, 1),
        (1, 2),
        (1, 3),
        (1, 4),
        (1, 5),
        (2, 2),
        (2, 3),
        (2, 4),
        (2, 5),
        (3, 3),
        (3, 4),
        (3, 5),
        (4, 4),
    ],
}


def measure_operator_basis(
    lat: Lattice,
    basis
) -> np.ndarray:

    return np.array([
        wilson_loop(lat, w, h)
        for w, h in basis
    ])


# ============================================================
# BLOCKING
# ============================================================

def block_lattice(
    lat: Lattice,
    b: int = 2
) -> Lattice:

    if lat.L % b != 0:
        raise ValueError(
            f"L={lat.L} must be divisible by b={b}"
        )

    Lc = lat.L // b

    coarse = create_lattice(Lc)

    for X in range(Lc):
        for Y in range(Lc):

            x = b * X
            y = b * Y

            # coarse x-link
            U = su2_identity()

            for k in range(b):

                U = su2_mul(
                    U,
                    get_link(
                        lat,
                        x + k,
                        y,
                        0
                    )
                )

            set_link(
                coarse,
                X,
                Y,
                0,
                U
            )

            # coarse y-link
            U = su2_identity()

            for k in range(b):

                U = su2_mul(
                    U,
                    get_link(
                        lat,
                        x,
                        y + k,
                        1
                    )
                )

            set_link(
                coarse,
                X,
                Y,
                1,
                U
            )

    return coarse


# ============================================================
# CONFIGURATION GENERATION
# ============================================================

def generate_ensemble(
    cfg: SimulationConfig,
    basis,
    rng: np.random.Generator
):

    lat = create_lattice(cfg.L)

    # --------------------------------------------------------
    # Thermalization
    # --------------------------------------------------------

    acceptance_history = []

    if cfg.verbose:
        print(
            f"[thermalization] L={cfg.L}, "
            f"beta={cfg.beta}"
        )

    for sweep in range(cfg.therm_sweeps):

        acc = metropolis_update(
            lat,
            cfg.beta,
            rng,
            cfg.proposal_epsilon
        )

        acceptance_history.append(acc)

    # --------------------------------------------------------
    # Measurements
    # --------------------------------------------------------

    fine_samples = []
    coarse_samples = []

    b = cfg.blocking_factor

    if cfg.verbose:
        print(
            f"[measurement] "
            f"{cfg.meas_sweeps} measurements"
        )

    for meas in range(cfg.meas_sweeps):

        for _ in range(
            cfg.sweep_between_measurements
        ):

            acc = metropolis_update(
                lat,
                cfg.beta,
                rng,
                cfg.proposal_epsilon
            )

            acceptance_history.append(acc)

        fine_obs = measure_operator_basis(
            lat,
            basis
        )

        coarse = block_lattice(
            lat,
            b=b
        )

        coarse_obs = measure_operator_basis(
            coarse,
            basis
        )

        fine_samples.append(fine_obs)
        coarse_samples.append(coarse_obs)

        if (
            cfg.verbose
            and (meas + 1) % max(
                1,
                cfg.meas_sweeps // 10
            ) == 0
        ):

            print(
                f"  measurement "
                f"{meas+1}/{cfg.meas_sweeps}"
            )

    return (
        np.asarray(fine_samples),
        np.asarray(coarse_samples),
        np.asarray(acceptance_history)
    )


# ============================================================
# MCRG MAP
# ============================================================

def covariance_matrix(X: np.ndarray) -> np.ndarray:

    Xc = X - np.mean(X, axis=0)

    return (
        Xc.T @ Xc
        /
        max(1, X.shape[0] - 1)
    )


def estimate_rg_jacobian(
    fine_samples: np.ndarray,
    coarse_samples: np.ndarray
) -> np.ndarray:

    """
    Linear MCRG estimator.

    We estimate the local relation

        coarse = A fine + c

    using least squares.

    This is a first numerical estimator of the local RG map.

    IMPORTANT:
    It is not yet the exact derivative of an analytically defined
    blocking transformation in coupling space.
    """

    X = fine_samples
    Y = coarse_samples

    Xc = X - X.mean(axis=0)
    Yc = Y - Y.mean(axis=0)

    A, _, _, _ = np.linalg.lstsq(
        Xc,
        Yc,
        rcond=None
    )

    return A.T


# ============================================================
# SPECTRAL ANALYSIS
# ============================================================

def spectral_analysis(J: np.ndarray):

    eigenvalues = np.linalg.eigvals(J)

    eigenvalues = np.asarray(eigenvalues)

    order = np.argsort(
        np.abs(eigenvalues)
    )[::-1]

    eigenvalues = eigenvalues[order]

    rho1 = float(np.abs(eigenvalues[0]))

    rho2 = (
        float(np.abs(eigenvalues[1]))
        if len(eigenvalues) > 1
        else np.nan
    )

    n_phi = (
        -math.log(rho1) / math.log(PHI)
        if rho1 > 0
        else np.nan
    )

    return {
        "eigenvalues_real": [
            float(np.real(x))
            for x in eigenvalues
        ],

        "eigenvalues_imag": [
            float(np.imag(x))
            for x in eigenvalues
        ],

        "rho1": rho1,
        "rho2": rho2,
        "n_phi": n_phi,
    }


# ============================================================
# BOOTSTRAP
# ============================================================

def bootstrap_jacobian(
    fine_samples: np.ndarray,
    coarse_samples: np.ndarray,
    n_bootstrap: int,
    rng: np.random.Generator
):

    N = len(fine_samples)

    rho_values = []

    nphi_values = []

    for _ in range(n_bootstrap):

        indices = rng.integers(
            0,
            N,
            size=N
        )

        X = fine_samples[indices]
        Y = coarse_samples[indices]

        J = estimate_rg_jacobian(
            X,
            Y
        )

        spec = spectral_analysis(J)

        rho_values.append(
            spec["rho1"]
        )

        nphi_values.append(
            spec["n_phi"]
        )

    rho_values = np.asarray(rho_values)
    nphi_values = np.asarray(nphi_values)

    return {
        "rho_mean": float(np.mean(rho_values)),
        "rho_std": float(np.std(rho_values, ddof=1)),

        "rho_q025": float(
            np.quantile(rho_values, 0.025)
        ),

        "rho_q975": float(
            np.quantile(rho_values, 0.975)
        ),

        "nphi_mean": float(np.mean(nphi_values)),
        "nphi_std": float(np.std(nphi_values, ddof=1)),

        "nphi_q025": float(
            np.quantile(nphi_values, 0.025)
        ),

        "nphi_q975": float(
            np.quantile(nphi_values, 0.975
            )
        ),
    }


# ============================================================
# PHI DIAGNOSTICS
# ============================================================

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

            "absolute_error": abs(
                rho - target
            ),

            "relative_error": abs(
                rho - target
            ) / target,

            "log_error": abs(
                math.log(rho / target)
            )
            if rho > 0
            else np.inf
        }

    best = min(
        result,
        key=lambda k: result[k]["log_error"]
    )

    result["closest"] = best

    return result


# ============================================================
# FINITE SIZE SCALING
# ============================================================

def fss_model(params, L):

    rho_inf, A, omega = params

    return (
        rho_inf
        +
        A * L ** (-omega)
    )


def fit_finite_size_scaling(
    L_values,
    rho_values
):

    L_values = np.asarray(
        L_values,
        dtype=float
    )

    rho_values = np.asarray(
        rho_values,
        dtype=float
    )

    if len(L_values) < 4:
        return None

    x0 = np.array([
        rho_values[-1],
        (
            rho_values[0]
            -
            rho_values[-1]
        )
        * L_values[0] ** 1.0,
        1.0
    ])

    def residuals(params):

        return (
            fss_model(params, L_values)
            -
            rho_values
        )

    result = least_squares(
        residuals,
        x0=x0,
        bounds=(
            [0.000001, -100.0, 0.05],
            [0.999999, 100.0, 10.0]
        )
    )

    rho_inf, A, omega = result.x

    return {
        "rho_inf": float(rho_inf),
        "A": float(A),
        "omega": float(omega),
        "cost": float(result.cost),
        "success": bool(result.success),
    }


# ============================================================
# RESULT WRITER
# ============================================================

def save_json(path, data):

    with open(
        path,
        "w",
        encoding="utf-8"
    ) as f:

        json.dump(
            data,
            f,
            indent=2
        )


def save_csv(path, rows):

    if not rows:
        return

    keys = list(rows[0].keys())

    with open(
        path,
        "w",
        newline="",
        encoding="utf-8"
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=keys
        )

        writer.writeheader()

        writer.writerows(rows)


# ============================================================
# SINGLE RUN
# ============================================================

def run_single(
    cfg: SimulationConfig,
    basis_name: str,
    rng: np.random.Generator
):

    basis = BASES[basis_name]

    print("\n" + "=" * 72)

    print(
        f"MCRG SU(2)"
        f" | L={cfg.L}"
        f" | beta={cfg.beta}"
        f" | b={cfg.blocking_factor}"
        f" | basis={basis_name}"
    )

    print("=" * 72)

    t0 = time.time()

    fine, coarse, acceptance = generate_ensemble(
        cfg,
        basis,
        rng
    )

    J = estimate_rg_jacobian(
        fine,
        coarse
    )

    spec = spectral_analysis(J)

    boot = bootstrap_jacobian(
        fine,
        coarse,
        cfg.bootstrap_samples,
        rng
    )

    phi = phi_diagnostics(
        spec["rho1"]
    )

    mean_acceptance = float(
        np.mean(acceptance)
    )

    elapsed = time.time() - t0

    result = {

        "L": cfg.L,
        "beta": cfg.beta,

        "blocking_factor":
            cfg.blocking_factor,

        "basis":
            basis_name,

        "N_op":
            len(basis),

        "N_samples":
            len(fine),

        "acceptance":
            mean_acceptance,

        "rho1":
            spec["rho1"],

        "rho2":
            spec["rho2"],

        "n_phi":
            spec["n_phi"],

        "bootstrap_rho_mean":
            boot["rho_mean"],

        "bootstrap_rho_std":
            boot["rho_std"],

        "bootstrap_rho_q025":
            boot["rho_q025"],

        "bootstrap_rho_q975":
            boot["rho_q975"],

        "bootstrap_nphi_mean":
            boot["nphi_mean"],

        "bootstrap_nphi_std":
            boot["nphi_std"],

        "bootstrap_nphi_q025":
            boot["nphi_q025"],

        "bootstrap_nphi_q975":
            boot["nphi_q975"],

        "phi_closest":
            phi["closest"],

        "phi_minus_1_error":
            phi["phi^-1"]["relative_error"],

        "phi_minus_2_error":
            phi["phi^-2"]["relative_error"],

        "phi_minus_3_error":
            phi["phi^-3"]["relative_error"],

        "phi_minus_4_error":
            phi["phi^-4"]["relative_error"],

        "elapsed_seconds":
            elapsed
    }

    print("\nRESULT")

    print(
        f"rho_1 = "
        f"{result['rho1']:.8f}"
    )

    print(
        f"rho_2 = "
        f"{result['rho2']:.8f}"
    )

    print(
        f"n_Phi = "
        f"{result['n_phi']:.8f}"
    )

    print(
        f"bootstrap rho = "
        f"{result['bootstrap_rho_mean']:.8f}"
        f" +/- "
        f"{result['bootstrap_rho_std']:.8f}"
    )

    print(
        f"bootstrap n_Phi = "
        f"{result['bootstrap_nphi_mean']:.8f}"
        f" +/- "
        f"{result['bootstrap_nphi_std']:.8f}"
    )

    print(
        f"closest golden target = "
        f"{result['phi_closest']}"
    )

    print(
        f"acceptance = "
        f"{result['acceptance']:.4f}"
    )

    return result, J


# ============================================================
# MAIN PIPELINE
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "MCRG SU(2) 2D first pipeline"
        )
    )

    parser.add_argument(
        "--L",
        type=int,
        nargs="+",
        default=[8, 16, 32]
    )

    parser.add_argument(
        "--beta",
        type=float,
        nargs="+",
        default=[2.0, 2.3, 2.6]
    )

    parser.add_argument(
        "--bases",
        type=str,
        nargs="+",
        default=["B1", "B2", "B3"]
    )

    parser.add_argument(
        "--blocking",
        type=int,
        nargs="+",
        default=[2]
    )

    parser.add_argument(
        "--therm",
        type=int,
        default=300
    )

    parser.add_argument(
        "--measurements",
        type=int,
        default=600
    )

    parser.add_argument(
        "--bootstrap",
        type=int,
        default=300
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=12345
    )

    parser.add_argument(
        "--output",
        type=str,
        default="mcrg_results"
    )

    args = parser.parse_args()

    output_dir = Path(
        args.output
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    all_results = []

    master_rng = np.random.default_rng(
        args.seed
    )

    for L in args.L:

        for beta in args.beta:

            for b in args.blocking:

                if L % b != 0:
                    continue

                # Coarse lattice must still be large
                # enough for the requested operators.
                Lc = L // b

                for basis_name in args.bases:

                    basis = BASES[basis_name]

                    max_loop = max(
                        max(w, h)
                        for w, h in basis
                    )

                    if Lc < 2 * max_loop:
                        print(
                            f"[SKIP] L={L}, b={b}, "
                            f"{basis_name}: "
                            f"coarse L={Lc} too small"
                        )
                        continue

                    # independent RNG for each run
                    seed = int(
                        master_rng.integers(
                            0,
                            2**32 - 1
                        )
                    )

                    rng = np.random.default_rng(
                        seed
                    )

                    cfg = SimulationConfig(
                        L=L,
                        beta=beta,

                        therm_sweeps=
                            args.therm,

                        meas_sweeps=
                            args.measurements,

                        blocking_factor=b,

                        bootstrap_samples=
                            args.bootstrap,

                        seed=seed,

                        output_dir=
                            str(output_dir)
                    )

                    result, J = run_single(
                        cfg,
                        basis_name,
                        rng
                    )

                    result["seed"] = seed

                    all_results.append(
                        result
                    )

                    # Save Jacobian
                    J_name = (
                        f"J_L{L}"
                        f"_b{b}"
                        f"_beta{beta:.3f}"
                        f"_{basis_name}.npy"
                    )

                    np.save(
                        output_dir / J_name,
                        J
                    )

    # --------------------------------------------------------
    # Save raw results
    # --------------------------------------------------------

    save_csv(
        output_dir / "mcrg_results.csv",
        all_results
    )

    save_json(
        output_dir / "mcrg_results.json",
        {
            "phi": float(PHI),
            "results": all_results
        }
    )

    # --------------------------------------------------------
    # FSS analysis grouped by beta/basis/b
    # --------------------------------------------------------

    fss_results = []

    groups = {}

    for row in all_results:

        key = (
            row["beta"],
            row["basis"],
            row["blocking_factor"]
        )

        groups.setdefault(
            key,
            []
        ).append(row)

    for key, rows in groups.items():

        rows = sorted(
            rows,
            key=lambda r: r["L"]
        )

        L_values = [
            r["L"]
            for r in rows
        ]

        rho_values = [
            r["rho1"]
            for r in rows
        ]

        fit = fit_finite_size_scaling(
            L_values,
            rho_values
        )

        if fit is None:
            continue

        beta, basis, b = key

        rho_inf = fit["rho_inf"]

        n_phi_inf = (
            -math.log(rho_inf)
            /
            math.log(PHI)
        )

        fss_results.append({

            "beta": beta,

            "basis": basis,

            "blocking_factor": b,

            "rho_inf": rho_inf,

            "A": fit["A"],

            "omega": fit["omega"],

            "n_phi_inf":
                n_phi_inf,

            "fit_cost":
                fit["cost"],

            "success":
                fit["success"]
        })

    save_csv(
        output_dir / "fss_results.csv",
        fss_results
    )

    save_json(
        output_dir / "fss_results.json",
        {
            "phi": float(PHI),
            "results": fss_results
        }
    )

    # --------------------------------------------------------
    # Final summary
    # --------------------------------------------------------

    print("\n")
    print("=" * 72)
    print("FINAL SUMMARY")
    print("=" * 72)

    for row in all_results:

        print(
            f"L={row['L']:3d} "
            f"beta={row['beta']:.3f} "
            f"b={row['blocking_factor']} "
            f"{row['basis']:>3s} "
            f"Nop={row['N_op']:2d} "
            f"rho={row['rho1']:.6f} "
            f"nPhi={row['n_phi']:.4f}"
        )

    print("\nFSS")

    for row in fss_results:

        print(
            f"beta={row['beta']:.3f} "
            f"{row['basis']:>3s} "
            f"b={row['blocking_factor']} "
            f"rho_inf={row['rho_inf']:.6f} "
            f"omega={row['omega']:.4f} "
            f"nPhi_inf={row['n_phi_inf']:.5f}"
        )

    print("\nOutput:")
    print(output_dir.resolve())


if __name__ == "__main__":
    main()

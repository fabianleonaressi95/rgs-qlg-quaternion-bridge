#!/usr/bin/env python3
"""
QVSM UNIVERSALITY / CONVERGENCE AUDIT

Scopo:
- generare dinamiche QVSM su sottoinsiemi casuali di microstati;
- costruire una densita radiale dinamica da occupazioni quantistiche;
- ottenere la massa cumulativa M_eff(r)/M;
- fittare la legge interna M_eff ~ r^p;
- fittare la legge esterna delta = 1-M_eff/M ~ r^p;
- confrontare stabilita degli esponenti al variare di:
    seed,
    numero di stati,
    width del kernel,
    coupling.

IMPORTANTE:
Questo è un modello efficace. La mappa:
    stato QVSM -> scala radiale
e il kernel radiale sono assunzioni dichiarate.

Il programma non dimostra PSU, gravita quantistica o bounce.
"""

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


# ============================================================
# Costanti
# ============================================================

PLANCK_LENGTH_M = 1.616255e-35


# ============================================================
# Funzioni QVSM
# ============================================================

def qvsm_index(G, L, C, T, B):
    return ((((G * 6 + L) * 6 + C) * 6 + T) * 5 + B)


def cyclic_distance(a, b, modulus):
    return min(abs(a - b), modulus - abs(a - b))


def qvsm_distance_squared(X, Y):
    moduli = [64, 6, 6, 6, 5]

    return sum(
        cyclic_distance(a, b, modulus) ** 2
        for a, b, modulus in zip(X, Y, moduli)
    )


def generate_states(n_states, seed):
    rng = np.random.default_rng(seed)

    states = []

    for _ in range(n_states):
        states.append((
            int(rng.integers(0, 64)),
            int(rng.integers(0, 6)),
            int(rng.integers(0, 6)),
            int(rng.integers(0, 6)),
            int(rng.integers(0, 5)),
        ))

    return states


def state_energy(X):
    G, L, C, T, B = X

    g = (G - 31.5) / 31.5
    l = (L - 2.5) / 2.5
    c = (C - 2.5) / 2.5
    t = (T - 2.5) / 2.5
    b = (B - 2.0) / 2.0

    return (
        0.20 * g**2
        + 0.10 * l**2
        + 0.10 * c**2
        + 0.10 * t**2
        + 0.10 * b**2
    )


# ============================================================
# Hamiltoniana QVSM locale
# ============================================================

def build_hamiltonian(states, coupling):
    """
    Hamiltoniana reale simmetrica:

        H_ii = E_X
        H_ij = -J/(1+d_X^2), solo per stati vicini.

    La condizione d_X^2 <= 2 rende l'accoppiamento sparso.
    """
    n = len(states)

    H = np.zeros((n, n), dtype=float)

    for i, X in enumerate(states):
        H[i, i] = state_energy(X)

    for i in range(n):
        for j in range(i):
            d2 = qvsm_distance_squared(
                states[i],
                states[j],
            )

            if d2 <= 2:
                J = -coupling / (1.0 + d2)

                H[i, j] = J
                H[j, i] = J

    return H


# ============================================================
# Mappa stato -> scala radiale
# ============================================================

def state_radial_scale(r_star_m, energy):
    """
    Ansatz esplicito:

        r_X = r_star exp(2 E_X)

    Non deriva da una teoria gravitazionale fondamentale.
    """
    return r_star_m * np.exp(2.0 * energy)


def lognormal_volume_density(
    radii_m,
    r_scale_m,
    kernel_width_log,
):
    """
    Kernel radiale di densita per volume.

    La massa cumulativa si ottiene poi tramite:

        dM ~ 4*pi*r^2*rho(r) dr

    Non viene normalizzato sul valore del bordo.
    """
    safe_r = np.maximum(radii_m, 1.0e-300)

    log_ratio = np.log(safe_r / r_scale_m)

    gaussian = np.exp(
        -0.5
        * (log_ratio / kernel_width_log) ** 2
    )

    denominator = (
        4.0
        * np.pi
        * safe_r**3
        * kernel_width_log
        * np.sqrt(2.0 * np.pi)
    )

    return gaussian / denominator


def cumulative_trapezoid(y, x):
    result = np.zeros_like(y, dtype=float)

    increments = (
        0.5
        * (y[1:] + y[:-1])
        * np.diff(x)
    )

    result[1:] = np.cumsum(increments)

    return result


def cumulative_mass_fraction(
    radii_m,
    density_volume,
):
    shell_density = (
        4.0
        * np.pi
        * radii_m**2
        * density_volume
    )

    cumulative = cumulative_trapezoid(
        shell_density,
        radii_m,
    )

    total = cumulative[-1]

    if (
        not np.isfinite(total)
        or total <= 0.0
    ):
        raise ValueError(
            "Massa totale non valida."
        )

    fraction = cumulative / total

    return np.clip(fraction, 0.0, 1.0)


# ============================================================
# Fit di legge di potenza
# ============================================================

def power_law_fit(
    x,
    y,
    x_min,
    x_max,
    y_min=None,
    y_max=None,
):
    """
    Fit log-log:

        y = A x^p

    Restituisce:
        p,
        A,
        n_points,
        R2,
        fit_status.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    mask = (
        np.isfinite(x)
        & np.isfinite(y)
        & (x > 0.0)
        & (y > 0.0)
        & (x >= x_min)
        & (x <= x_max)
    )

    if y_min is not None:
        mask &= y >= y_min

    if y_max is not None:
        mask &= y <= y_max

    n_points = int(np.count_nonzero(mask))

    if n_points < 8:
        return {
            "slope": np.nan,
            "prefactor": np.nan,
            "n_points": n_points,
            "R2": np.nan,
            "status": "INSUFFICIENT_POINTS",
        }

    log_x = np.log(x[mask])
    log_y = np.log(y[mask])

    slope, intercept = np.polyfit(
        log_x,
        log_y,
        1,
    )

    prediction = slope * log_x + intercept

    ss_res = np.sum(
        (log_y - prediction) ** 2
    )

    ss_tot = np.sum(
        (log_y - np.mean(log_y)) ** 2
    )

    R2 = (
        1.0 - ss_res / ss_tot
        if ss_tot > 0.0
        else np.nan
    )

    return {
        "slope": float(slope),
        "prefactor": float(np.exp(intercept)),
        "n_points": n_points,
        "R2": float(R2),
        "status": "FIT_OK",
    }


# ============================================================
# Singola simulazione
# ============================================================

def run_single_case(
    seed,
    n_states,
    coupling,
    kernel_width_log,
    steps,
    dt,
    r_star_planck,
    n_radii,
):
    """
    Esegue una simulazione e restituisce:

    - audit row;
    - profilo radiale;
    - spettro QVSM.
    """
    states = generate_states(
        n_states,
        seed,
    )

    H = build_hamiltonian(
        states,
        coupling,
    )

    eigenvalues, eigenvectors = np.linalg.eigh(H)

    r_star_m = (
        r_star_planck
        * PLANCK_LENGTH_M
    )

    radii_m = np.geomspace(
        1.0e-8 * r_star_m,
        1.0e8 * r_star_m,
        n_radii,
    )

    radial_scales = np.array([
        state_radial_scale(
            r_star_m,
            energy,
        )
        for energy in eigenvalues
    ])

    kernels = np.array([
        lognormal_volume_density(
            radii_m,
            r_scale,
            kernel_width_log,
        )
        for r_scale in radial_scales
    ])

    # Stato iniziale.
    psi = np.zeros(
        n_states,
        dtype=complex,
    )

    psi[0] = 1.0 + 0.0j

    coefficients = (
        eigenvectors.conj().T @ psi
    )

    phase_step = np.exp(
        -1j
        * eigenvalues
        * dt
    )

    for _ in range(steps):
        coefficients *= phase_step

    psi = eigenvectors @ coefficients
    psi /= np.linalg.norm(psi)

    probabilities = np.abs(psi)**2

    density_volume = probabilities @ kernels

    M_fraction = cumulative_mass_fraction(
        radii_m,
        density_volume,
    )

    x = radii_m / r_star_m

    # Regime interno:
    #
    # M_eff/M ~ x^p
    #
    # Fit su range interno definito in funzione di x.
    inner_fit = power_law_fit(
        x=x,
        y=M_fraction,
        x_min=1.0e-5,
        x_max=1.0e-2,
        y_min=1.0e-300,
        y_max=1.0e-2,
    )

    # Regime esterno:
    #
    # delta = 1 - M_eff/M ~ x^p
    #
    # Si escludono sia la transizione sia la saturazione
    # numerica a delta=0.
    delta = 1.0 - M_fraction

    outer_fit = power_law_fit(
        x=x,
        y=delta,
        x_min=2.0,
        x_max=1.0e6,
        y_min=1.0e-8,
        y_max=1.0e-2,
    )

    # Scala di recupero classico.
    recovery = {}

    for tolerance in [
        1.0e-1,
        1.0e-2,
        1.0e-3,
        1.0e-4,
    ]:
        mask = delta <= tolerance

        if np.any(mask):
            recovery[f"r_recovery_{tolerance:.0e}"] = float(
                x[np.where(mask)[0][0]]
            )
        else:
            recovery[f"r_recovery_{tolerance:.0e}"] = np.nan

    audit_row = {
        "seed": seed,
        "n_states": n_states,
        "coupling": coupling,
        "kernel_width_log": kernel_width_log,
        "steps": steps,
        "dt": dt,
        "r_star_planck": r_star_planck,
        "participation_ratio": float(
            1.0 / np.sum(probabilities**2)
        ),
        "max_probability": float(
            np.max(probabilities)
        ),
        "mean_energy": float(
            np.sum(probabilities * eigenvalues)
        ),
        "inner_p": inner_fit["slope"],
        "inner_prefactor": inner_fit["prefactor"],
        "inner_n_points": inner_fit["n_points"],
        "inner_R2": inner_fit["R2"],
        "inner_status": inner_fit["status"],
        "outer_p": outer_fit["slope"],
        "outer_prefactor": outer_fit["prefactor"],
        "outer_n_points": outer_fit["n_points"],
        "outer_R2": outer_fit["R2"],
        "outer_status": outer_fit["status"],
        "M_fraction_near": float(M_fraction[0]),
        "M_fraction_at_rstar": float(
            M_fraction[
                np.argmin(
                    np.abs(x - 1.0)
                )
            ]
        ),
        "M_fraction_far": float(M_fraction[-1]),
        "M_fraction_min": float(
            np.min(M_fraction)
        ),
        "M_fraction_max": float(
            np.max(M_fraction)
        ),
        **recovery,
    }

    radial_df = pd.DataFrame({
        "r_m": radii_m,
        "r_over_rstar": x,
        "rho_dynamic_volume": density_volume,
        "M_eff_over_M": M_fraction,
        "delta_outer": delta,
    })

    spectrum_rows = []

    for i, X in enumerate(states):
        G, L, C, T, B = X

        spectrum_rows.append({
            "state_order": i,
            "G": G,
            "L": L,
            "C": C,
            "T": T,
            "B": B,
            "qvsm_index": qvsm_index(
                G, L, C, T, B
            ),
            "eigenvalue": eigenvalues[i],
            "probability_final": probabilities[i],
            "radial_scale_over_rstar": (
                radial_scales[i] / r_star_m
            ),
        })

    spectrum_df = pd.DataFrame(
        spectrum_rows
    )

    return audit_row, radial_df, spectrum_df


# ============================================================
# Aggregazione
# ============================================================

def aggregate_results(audit_df):
    """
    Riassunto della stabilita degli esponenti.
    """
    rows = []

    for label, column in [
        ("inner", "inner_p"),
        ("outer", "outer_p"),
    ]:
        values = audit_df[column].dropna()

        if len(values) == 0:
            rows.append({
                "regime": label,
                "n_valid": 0,
                "mean_p": np.nan,
                "std_p": np.nan,
                "median_p": np.nan,
                "q025_p": np.nan,
                "q975_p": np.nan,
                "distance_from_3": np.nan,
            })

            continue

        target = 3.0

        if label == "outer":
            target = -3.0

        rows.append({
            "regime": label,
            "n_valid": len(values),
            "mean_p": float(values.mean()),
            "std_p": float(values.std(ddof=1)),
            "median_p": float(values.median()),
            "q025_p": float(values.quantile(0.025)),
            "q975_p": float(values.quantile(0.975)),
            "distance_from_3": float(
                abs(values.mean() - target)
            ),
        })

    return pd.DataFrame(rows)


def sha256_file(path):
    digest = hashlib.sha256()

    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)

    return digest.hexdigest()


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description=(
            "QVSM universality and convergence audit."
        )
    )

    parser.add_argument(
        "--outdir",
        default="qvsm_universality_audit_output",
    )

    parser.add_argument(
        "--seeds",
        nargs="+",
        type=int,
        default=[1, 2, 3, 4, 5],
    )

    parser.add_argument(
        "--n-states",
        nargs="+",
        type=int,
        default=[128, 256, 512],
    )

    parser.add_argument(
        "--kernel-widths",
        nargs="+",
        type=float,
        default=[0.4, 0.6, 0.8],
    )

    parser.add_argument(
        "--couplings",
        nargs="+",
        type=float,
        default=[0.02, 0.08, 0.20],
    )

    parser.add_argument(
        "--steps",
        type=int,
        default=500,
    )

    parser.add_argument(
        "--dt",
        type=float,
        default=0.01,
    )

    parser.add_argument(
        "--rstar-planck",
        type=float,
        default=1.0,
    )

    parser.add_argument(
        "--n-radii",
        type=int,
        default=1200,
    )

    parser.add_argument(
        "--save-all-profiles",
        action="store_true",
        help=(
            "Salva un CSV radiale e uno spettro per ogni "
            "combinazione. Senza questo flag salva solo "
            "il primo e l'ultimo profilo."
        ),
    )

    args = parser.parse_args()

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    profiles_dir = outdir / "profiles"
    spectra_dir = outdir / "spectra"

    profiles_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    spectra_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    parameter_grid = list(
        itertools.product(
            args.seeds,
            args.n_states,
            args.kernel_widths,
            args.couplings,
        )
    )

    print(
        "Numero combinazioni:",
        len(parameter_grid),
    )

    audit_rows = []

    for run_number, (
        seed,
        n_states,
        kernel_width,
        coupling,
    ) in enumerate(parameter_grid, start=1):
        print(
            f"[{run_number}/{len(parameter_grid)}] "
            f"seed={seed} "
            f"n_states={n_states} "
            f"width={kernel_width} "
            f"coupling={coupling}"
        )

        audit_row, radial_df, spectrum_df = run_single_case(
            seed=seed,
            n_states=n_states,
            coupling=coupling,
            kernel_width_log=kernel_width,
            steps=args.steps,
            dt=args.dt,
            r_star_planck=args.rstar_planck,
            n_radii=args.n_radii,
        )

        audit_row["run_number"] = run_number

        audit_rows.append(audit_row)

        tag = (
            f"seed{seed}_N{n_states}_"
            f"w{kernel_width:.3f}_"
            f"J{coupling:.3f}"
        )

        save_this = (
            args.save_all_profiles
            or run_number == 1
            or run_number == len(parameter_grid)
        )

        if save_this:
            radial_df.to_csv(
                profiles_dir / f"profile_{tag}.csv",
                index=False,
            )

            spectrum_df.to_csv(
                spectra_dir / f"spectrum_{tag}.csv",
                index=False,
            )

    audit_df = pd.DataFrame(audit_rows)

    aggregate_df = aggregate_results(
        audit_df
    )

    # Riassunto per gruppo di parametri.
    group_df = (
        audit_df
        .groupby(
            [
                "n_states",
                "kernel_width_log",
                "coupling",
            ],
            dropna=False,
        )
        .agg(
            n_runs=("seed", "count"),
            inner_p_mean=("inner_p", "mean"),
            inner_p_std=("inner_p", "std"),
            outer_p_mean=("outer_p", "mean"),
            outer_p_std=("outer_p", "std"),
            outer_R2_mean=("outer_R2", "mean"),
            inner_R2_mean=("inner_R2", "mean"),
        )
        .reset_index()
    )

    audit_path = outdir / "qvsm_universality_audit_runs.csv"
    aggregate_path = (
        outdir
        / "qvsm_universality_audit_summary.csv"
    )
    group_path = (
        outdir
        / "qvsm_universality_by_parameters.csv"
    )

    audit_df.to_csv(
        audit_path,
        index=False,
    )

    aggregate_df.to_csv(
        aggregate_path,
        index=False,
    )

    group_df.to_csv(
        group_path,
        index=False,
    )

    manifest = {
        "project": "QVSM_UNIVERSALITY_CONVERGENCE_AUDIT",
        "purpose": (
            "Test whether inner r^3 and outer r^-3 "
            "scalings emerge robustly from the sampled "
            "QVSM radial-density model."
        ),
        "parameter_grid": {
            "seeds": args.seeds,
            "n_states": args.n_states,
            "kernel_widths": args.kernel_widths,
            "couplings": args.couplings,
            "steps": args.steps,
            "dt": args.dt,
            "rstar_planck": args.rstar_planck,
            "n_radii": args.n_radii,
        },
        "expected_reference_scalings": {
            "inner": "M_eff/M proportional to r^3",
            "outer": (
                "1-M_eff/M proportional to r^-3"

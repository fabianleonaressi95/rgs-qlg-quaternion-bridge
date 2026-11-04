#!/usr/bin/env python3

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


G_NEWTON = 6.67430e-11
PLANCK_LENGTH_M = 1.616255e-35
PLANCK_MASS_KG = 2.176434e-8


def qvsm_index(G, L, C, T, B):
    return ((((G * 6 + L) * 6 + C) * 6 + T) * 5 + B)


def sha256_file(path):
    digest = hashlib.sha256()

    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)

    return digest.hexdigest()


def cyclic_distance(a, b, modulus):
    return min(abs(a - b), modulus - abs(a - b))


def qvsm_distance_squared(X, Y):
    moduli = [64, 6, 6, 6, 5]

    return sum(
        cyclic_distance(a, b, n) ** 2
        for a, b, n in zip(X, Y, moduli)
    )


def generate_states(n_states, seed):
    rng = np.random.default_rng(seed)

    states = []

    for _ in range(n_states):
        G = int(rng.integers(0, 64))
        L = int(rng.integers(0, 6))
        C = int(rng.integers(0, 6))
        T = int(rng.integers(0, 6))
        B = int(rng.integers(0, 5))

        states.append((G, L, C, T, B))

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


def build_hamiltonian(states, coupling):
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

            # Accoppiamento solo tra stati vicini.
            if d2 <= 2:
                value = -coupling / (1.0 + d2)
                H[i, j] = value
                H[j, i] = value

    return H


def state_radial_scale(r_star_m, energy):
    """
    Mappa esplicita stato -> scala radiale.

    Questa è una scelta efficace dichiarata. Non è una
    derivazione da un vincolo gravitazionale fondamentale.
    """
    return r_star_m * np.exp(2.0 * energy)


def radial_kernel(radii_m, r_scale_m, width_log=0.5):
    log_ratio = np.log(radii_m / r_scale_m)

    return np.exp(
        -(log_ratio**2) / width_log
    )


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--outdir",
        default="qvsm_dynamic_density_output",
    )

    parser.add_argument(
        "--n-states",
        type=int,
        default=256,
    )

    parser.add_argument(
        "--steps",
        type=int,
        default=1000,
    )

    parser.add_argument(
        "--dt",
        type=float,
        default=0.01,
    )

    parser.add_argument(
        "--coupling",
        type=float,
        default=0.08,
    )

    parser.add_argument(
        "--rstar-planck",
        type=float,
        default=1.0,
    )

    parser.add_argument(
        "--M-planck",
        type=float,
        default=1.0,
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=3,
    )

    args = parser.parse_args()

    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    print("Genero stati QVSM...")
    states = generate_states(
        args.n_states,
        args.seed,
    )

    print("Costruisco Hamiltoniana...")
    H = build_hamiltonian(
        states,
        args.coupling,
    )

    print("Diagonalizzo Hamiltoniana...")
    eigenvalues, eigenvectors = np.linalg.eigh(H)

    r_star_m = (
        args.rstar_planck
        * PLANCK_LENGTH_M
    )

    M_total_kg = (
        args.M_planck
        * PLANCK_MASS_KG
    )

    radii_m = np.geomspace(
        0.01 * r_star_m,
        1.0e6 * r_star_m,
        300,
    )

    # Scala radiale di ciascun autostato.
    radial_scales_m = np.array([
        state_radial_scale(
            r_star_m,
            energy,
        )
        for energy in eigenvalues
    ])

    kernels = np.array([
        radial_kernel(
            radii_m,
            scale,
        )
        for scale in radial_scales_m
    ])

    # Stato iniziale localizzato nel primo microstato.
    psi = np.zeros(
        args.n_states,
        dtype=complex,
    )

    psi[0] = 1.0 + 0.0j

    # Passaggio allo spazio degli autostati.
    coefficients = (
        eigenvectors.conj().T @ psi
    )

    phase_step = np.exp(
        -1j
        * eigenvalues
        * args.dt
    )

    trajectory_rows = []

    print("Evoluzione dinamica...")

    for step in range(args.steps):
        coefficients *= phase_step

        psi = eigenvectors @ coefficients
        psi /= np.linalg.norm(psi)

        probabilities = np.abs(psi)**2

        # Densita radiale dinamica.
        profile = probabilities @ kernels

        # Il limite lontano viene normalizzato a 1
        # per controllare il recupero classico.
        far_value = profile[-1]

        if far_value <= 1.0e-300:
            normalized_profile = np.zeros_like(profile)
        else:
            normalized_profile = profile / far_value

        M_eff_kg = (
            M_total_kg
            * normalized_profile
        )

        g_m_s2 = (
            -G_NEWTON
            * M_eff_kg
            / radii_m**2
        )

        flux = (
            4.0
            * np.pi
            * radii_m**2
            * g_m_s2
        )

        classical_flux = (
            -4.0
            * np.pi
            * G_NEWTON
            * M_total_kg
        )

        mean_energy = float(
            np.sum(
                probabilities
                * eigenvalues
            )
        )

        trajectory_rows.append({
            "step": step,
            "time_model_units": step * args.dt,
            "mean_state_energy": mean_energy,
            "F_dynamic_near": float(
                normalized_profile[0]
            ),
            "F_dynamic_far": float(
                normalized_profile[-1]
            ),
            "M_eff_near_over_M": float(
                M_eff_kg[0] / M_total_kg
            ),
            "M_eff_far_over_M": float(
                M_eff_kg[-1] / M_total_kg
            ),
            "min_M_eff_over_M": float(
                np.min(M_eff_kg / M_total_kg)
            ),
            "max_M_eff_over_M": float(
                np.max(M_eff_kg / M_total_kg)
            ),
            "flux_near_over_classical": float(
                flux[0] / classical_flux
            ),
            "flux_far_over_classical": float(
                flux[-1] / classical_flux
            ),
            "max_probability": float(
                np.max(probabilities)
            ),
            "participation_ratio": float(
                1.0 / np.sum(probabilities**2)
            ),
        })

    trajectory_df = pd.DataFrame(
        trajectory_rows
    )

    state_rows = []

    for i, X in enumerate(states):
        G, L, C, T, B = X

        state_rows.append({
            "state_order": i,
            "G": G,
            "L": L,
            "C": C,
            "T": T,
            "B": B,
            "qvsm_index": qvsm_index(
                G, L, C, T, B
            ),
            "eigenvalue_model_units": eigenvalues[i],
            "radial_scale_m": radial_scales_m[i],
            "radial_scale_over_rstar": (
                radial_scales_m[i]
                / r_star_m
            ),
        })

    states_df = pd.DataFrame(state_rows)

    # Salva il profilo finale completo.
    final_probabilities = np.abs(psi)**2
    final_profile = final_probabilities @ kernels
    final_profile /= final_profile[-1]

    final_mass = M_total_kg * final_profile

    final_g = (
        -G_NEWTON
        * final_mass
        / radii_m**2
    )

    final_flux = (
        4.0
        * np.pi
        * radii_m**2
        * final_g
    )

    radial_df = pd.DataFrame({
        "r_m": radii_m,
        "r_over_rstar": radii_m / r_star_m,
        "profile_dynamic": final_profile,
        "M_eff_kg": final_mass,
        "M_eff_over_M": final_mass / M_total_kg,
        "g_m_s2": final_g,
        "flux": final_flux,
        "flux_over_classical": (
            final_flux
            / (
                -4.0
                * np.pi
                * G_NEWTON
                * M_total_kg
            )
        ),
    })

    trajectory_path = (
        outdir
        / "dynamic_density_trajectory.csv"
    )

    spectrum_path = (
        outdir
        / "qvsm_dynamic_spectrum.csv"
    )

    radial_path = (
        outdir
        / "dynamic_radial_profile.csv"
    )

    trajectory_df.to_csv(
        trajectory_path,
        index=False,
    )

    states_df.to_csv(
        spectrum_path,
        index=False,
    )

    radial_df.to_csv(
        radial_path,
        index=False,
    )

    manifest = {
        "status": (
            "DYNAMIC_QVSM_DENSITY_EFFECTIVE_TEST"
        ),
        "n_states": args.n_states,
        "steps": args.steps,
        "dt": args.dt,
        "coupling": args.coupling,
        "r_star_m": r_star_m,
        "M_total_kg": M_total_kg,
        "model": {
            "state_evolution": (
                "unitary evolution under sampled QVSM Hamiltonian"
            ),
            "radial_profile": (
                "occupation-weighted logarithmic radial kernels"
            ),
            "critical_limit": (
                "profile normalized to one at largest sampled radius"
            ),
        },
        "limitations": [
            (
                "The state-to-radial-scale map is an explicit "
                "effective assumption."
            ),
            (
                "The radial kernel is not derived from a "
                "Hamiltonian gravitational constraint."
            ),
            (
                "This does not demonstrate PSU, LQG, "
                "or a physical black-to-white transition."
            ),
        ],
        "files_sha256": {},
    }

    for path in [
        trajectory_path,
        spectrum_path,
        radial_path,
    ]:
        manifest["files_sha256"][path.name] = (
            sha256_file(path)
        )

    manifest_path = outdir / "manifest.json"

    manifest_path.write_text(
        json.dumps(
            manifest,
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )

    print()
    print("=== COMPLETATO ===")
    print("Output:", outdir)
    print("File creati:")
    print(" ", trajectory_path.name)
    print(" ", spectrum_path.name)
    print(" ", radial_path.name)
    print(" ", manifest_path.name)
    print()
    print("Ultimo passo simulato:")
    print(
        trajectory_df.tail(1).to_string(
            index=False
        )
    )


if __name__ == "__main__":
    main()

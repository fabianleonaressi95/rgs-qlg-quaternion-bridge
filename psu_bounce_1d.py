#!/usr/bin/env python3

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


HBAR = 1.054571817e-34
G_NEWTON = 6.67430e-11
C_LIGHT = 299792458.0
PLANCK_LENGTH = 1.616255e-35
PLANCK_MASS = 2.176434e-8
PLANCK_TIME = 5.391247e-44


def laplacian_matrix(n, dx):
    matrix = np.zeros((n, n), dtype=complex)

    diagonal = -2.0 * np.ones(n)
    upper = np.ones(n - 1)
    lower = np.ones(n - 1)

    matrix += np.diag(diagonal)
    matrix += np.diag(upper, k=1)
    matrix += np.diag(lower, k=-1)

    return matrix / dx**2


def normalize(psi, dx):
    norm = np.sqrt(
        np.sum(np.abs(psi) ** 2) * dx
    )

    if norm == 0:
        raise ValueError("Funzione d'onda nulla")

    return psi / norm


def expectation(psi, operator, dx):
    return np.real(
        np.vdot(psi, operator @ psi) * dx
    )


def gaussian_wavepacket(q, center, width, momentum):
    envelope = np.exp(
        -0.5 * ((q - center) / width) ** 2
    )

    phase = np.exp(1j * momentum * q)

    return envelope * phase


def effective_potential(
    q,
    mass_kg,
    eta_l,
    alpha_qg,
    exponent_n,
):
    """
    Modello fenomenologico:

        V_GR = -G M mu / q
        V_QG = alpha * (ell_* / q)^n * E_scale

    Non è una derivazione LQG.
    """
    ell_star = eta_l * PLANCK_LENGTH

    # Massa gravitazionale di esempio del sistema radiale.
    # Deve essere sostituita da una relazione derivata.
    M_source = PLANCK_MASS

    # Massa efficace del grado radiale.
    mu_q = PLANCK_MASS

    q_safe = np.maximum(q, ell_star * 1.0e-3)

    V_gr = (
        -G_NEWTON
        * M_source
        * mu_q
        / q_safe
    )

    E_scale = (
        mu_q * C_LIGHT**2
    )

    V_qg = (
        alpha_qg
        * E_scale
        * (ell_star / q_safe) ** exponent_n
    )

    return V_gr + V_qg


def classify_probabilities(q, psi, p_operator, dx, q_min, p_threshold):
    probability_density = np.abs(psi) ** 2

    p_q_negative = expectation(
        psi,
        -p_operator,
        dx,
    )

    # Classificazione radiale tramite pressione locale.
    # Non è una misura covariante di orizzonte.
    inward = q < q_min
    near_bounce = np.abs(q - q_min) <= 0.05 * q_min
    outward = q > q_min

    P_bh = np.sum(
        probability_density[inward]
    ) * dx

    P_bounce = np.sum(
        probability_density[near_bounce]
    ) * dx

    P_wh = np.sum(
        probability_density[outward]
    ) * dx

    return {
        "P_BH_radial": float(P_bh),
        "P_bounce_radial": float(P_bounce),
        "P_WH_radial": float(P_wh),
        "mean_minus_p": float(p_q_negative),
    }


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--outdir",
        default="psu_bounce_1d_output",
    )

    parser.add_argument(
        "--grid-size",
        type=int,
        default=1200,
    )

    parser.add_argument(
        "--q-min-planck",
        type=float,
        default=0.5,
    )

    parser.add_argument(
        "--q-max-planck",
        type=float,
        default=20.0,
    )

    parser.add_argument(
        "--eta-l",
        type=float,
        default=1.0,
    )

    parser.add_argument(
        "--alpha-qg",
        type=float,
        default=1.0,
    )

    parser.add_argument(
        "--n-qg",
        type=int,
        default=4,
    )

    parser.add_argument(
        "--steps",
        type=int,
        default=2500,
    )

    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    ell_star = args.eta_l * PLANCK_LENGTH

    q = np.linspace(
        args.q_min_planck * ell_star,
        args.q_max_planck * ell_star,
        args.grid_size,
    )

    dx = q[1] - q[0]

    potential = effective_potential(
        q=q,
        mass_kg=PLANCK_MASS,
        eta_l=args.eta_l,
        alpha_qg=args.alpha_qg,
        exponent_n=args.n_qg,
    )

    hbar2_over_2m = (
        HBAR**2
        / (2.0 * PLANCK_MASS)
    )

    kinetic = (
        -hbar2_over_2m
        * laplacian_matrix(args.grid_size, dx)
    )

    hamiltonian = kinetic + np.diag(potential)

    eigenvalues, eigenvectors = np.linalg.eigh(hamiltonian)

    ground = eigenvectors[:, 0]
    ground = normalize(ground, dx)

    # Pacchetto iniziale in contrazione.
    center = 10.0 * ell_star
    width = 0.8 * ell_star
    momentum = -5.0 / ell_star

    psi = gaussian_wavepacket(
        q=q,
        center=center,
        width=width,
        momentum=momentum,
    )

    psi = normalize(psi, dx)

    # Propagazione split-step approssimata tramite autodecomposizione.
    # Per griglie grandi, sostituire con propagatore sparso.
    dt = 0.01 * PLANCK_TIME

    phases = np.exp(
        -1j * eigenvalues * dt / HBAR
    )

    q_min_index = int(
        np.argmin(potential)
    )

    q_min = q[q_min_index]

    rows = []

    eigenvector_matrix = eigenvectors
    coefficients = eigenvector_matrix.conj().T @ (
        psi * dx
    )

    p_operator = np.zeros(
        (args.grid_size, args.grid_size),
        dtype=complex,
    )

    first_derivative = np.zeros_like(p_operator)

    for i in range(1, args.grid_size - 1):
        first_derivative[i, i + 1] = 1.0
        first_derivative[i, i - 1] = -1.0

    first_derivative /= 2.0 * dx

    p_operator = -1j * HBAR * first_derivative

    for step in range(args.steps):
        coefficients *= phases
        psi = eigenvector_matrix @ coefficients
        psi = normalize(psi, dx)

        if step % 10 != 0:
            continue

        classification = classify_probabilities(
            q=q,
            psi=psi,
            p_operator=p_operator,
            dx=dx,
            q_min=q_min,
            p_threshold=0.0,
        )

        rows.append({
            "step": step,
            "time_planck_units": step * dt / PLANCK_TIME,
            "q_expectation_m": expectation(
                psi,
                np.diag(q),
                dx,
            ),
            "q_min_m": q_min,
            "potential_min_J": float(np.min(potential)),
            "energy_expectation_J": expectation(
                psi,
                hamiltonian,
                dx,
            ),
            **classification,
        })

    trajectory = pd.DataFrame(rows)
    trajectory.to_csv(
        outdir / "bounce_trajectory.csv",
        index=False,
    )

    pd.DataFrame({
        "q_m": q,
        "potential_J": potential,
        "ground_state_probability_density": np.abs(ground) ** 2,
    }).to_csv(
        outdir / "effective_potential.csv",
        index=False,
    )

    manifest = {
        "model": "1D effective radial PSU bounce",
        "status": "PHENOMENOLOGICAL_NOT_FUNDAMENTAL_QG",
        "parameters": {
            "eta_l": args.eta_l,
            "alpha_qg": args.alpha_qg,
            "n_qg": args.n_qg,
            "mass_source_kg": PLANCK_MASS,
            "mass_radial_kg": PLANCK_MASS,
        },
        "planck_length_m": PLANCK_LENGTH,
        "planck_mass_kg": PLANCK_MASS,
        "planck_time_s": PLANCK_TIME,
        "q_min_m": float(q_min),
        "q_min_over_ell_star": float(q_min / ell_star),
        "limitations": [
            "The potential is phenomenological.",
            "No covariant metric or horizon finder is implemented.",
            "BH/bounce/WH labels are radial wavepacket labels.",
            "No PSU ontology is established.",
            "No mass law is derived.",
            "No photon dispersion is derived from this Hamiltonian.",
        ],
    }

    (outdir / "bounce_manifest.json").write_text(
        json.dumps(
            manifest,
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )

    print("Output generato:", outdir)
    print("q_min / ell_star:", q_min / ell_star)
    print("File:")
    print("  bounce_trajectory.csv")
    print("  effective_potential.csv")
    print("  bounce_manifest.json")


if __name__ == "__main__":
    main()

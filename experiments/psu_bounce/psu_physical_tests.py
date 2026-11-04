#!/usr/bin/env python3

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


# Costanti SI
C_LIGHT = 299792458.0
H_PLANCK = 6.62607015e-34
HBAR = H_PLANCK / (2.0 * np.pi)
EV_J = 1.602176634e-19
PLANCK_ENERGY_EV = 1.220890e28
PLANCK_LENGTH_M = 1.616255e-35

# Cosmologia di riferimento semplice
H0_KM_S_MPC = 67.4
H0_SI = (
    H0_KM_S_MPC
    * 1000.0
    / 3.0856775814913673e22
)
OMEGA_M = 0.315
OMEGA_L = 0.685


def qvsm_index(G, L, C, T, B):
    if any(pd.isna(x) for x in [G, L, C, T, B]):
        return None

    G, L, C, T, B = map(int, [G, L, C, T, B])

    if not 0 <= G < 64:
        raise ValueError("G fuori range")

    if not 0 <= L < 6:
        raise ValueError("L fuori range")

    if not 0 <= C < 6:
        raise ValueError("C fuori range")

    if not 0 <= T < 6:
        raise ValueError("T fuori range")

    if not 0 <= B < 5:
        raise ValueError("B fuori range")

    return ((((G * 6 + L) * 6 + C) * 6 + T) * 5 + B)


def cosmological_integral(z, n, points=10000):
    z_values = np.linspace(0.0, z, points)

    denominator = np.sqrt(
        OMEGA_M * (1.0 + z_values) ** 3
        + OMEGA_L
    )

    integrand = (1.0 + z_values) ** n / denominator

    return float(np.trapezoid(integrand, z_values))


def photon_delay_seconds(
    redshift,
    E_high_eV,
    E_low_eV,
    E_qg_eV,
    n,
    xi=1.0,
):
    if E_qg_eV <= 0:
        raise ValueError("E_QG deve essere positivo")

    integral = cosmological_integral(redshift, n)

    delta_E_n = (
        E_high_eV ** n
        - E_low_eV ** n
    )

    return float(
        xi
        * (n + 1.0)
        / (2.0 * H0_SI)
        * delta_E_n
        / E_qg_eV ** n
        * integral
    )


def qg_bound_from_delay(
    redshift,
    E_high_eV,
    E_low_eV,
    delay_limit_seconds,
    n,
    xi=1.0,
):
    integral = cosmological_integral(redshift, n)

    delta_E_n = (
        E_high_eV ** n
        - E_low_eV ** n
    )

    numerator = (
        abs(xi)
        * (n + 1.0)
        * delta_E_n
        * integral
    )

    denominator = (
        2.0
        * H0_SI
        * delay_limit_seconds
    )

    return float(
        (numerator / denominator) ** (1.0 / n)
    )


def lorentz_anisotropy_tensor(Q):
    Q = np.asarray(Q, dtype=float)

    trace = np.trace(Q)

    if trace == 0:
        raise ValueError("Tensore con traccia nulla")

    Q_norm = Q / trace
    eigenvalues = np.linalg.eigvalsh(Q_norm)

    return {
        "lambda_min": float(np.min(eigenvalues)),
        "lambda_mid": float(np.sort(eigenvalues)[1]),
        "lambda_max": float(np.max(eigenvalues)),
        "anisotropy": float(
            np.max(eigenvalues) - np.min(eigenvalues)
        ),
    }


def lv_direction_factor(
    direction,
    Q,
    isotropic_part=1.0,
):
    """
    f(n) = isotropic_part + n^T Q n

    Q deve essere traceless per evitare una ridefinizione
    della parte isotropa.
    """
    direction = np.asarray(direction, dtype=float)
    direction /= np.linalg.norm(direction)

    Q = np.asarray(Q, dtype=float)

    return float(
        isotropic_part
        + direction @ Q @ direction
    )


def psu_vacuum_density(
    probabilities,
    energies_J,
    volume_m3,
):
    probabilities = np.asarray(probabilities, dtype=float)
    energies_J = np.asarray(energies_J, dtype=float)

    if volume_m3 <= 0:
        raise ValueError("Il volume deve essere positivo")

    if len(probabilities) != len(energies_J):
        raise ValueError("Dimensioni incompatibili")

    if not np.isclose(np.sum(probabilities), 1.0):
        probabilities = probabilities / np.sum(probabilities)

    return float(
        np.sum(probabilities * energies_J) / volume_m3
    )


def newtonian_mass_from_psu(
    m_star_kg,
    occupancies,
    weights,
):
    occupancies = np.asarray(occupancies, dtype=float)
    weights = np.asarray(weights, dtype=float)

    if len(occupancies) != len(weights):
        raise ValueError("Dimensioni incompatibili")

    return float(
        m_star_kg
        * np.sum(occupancies * weights)
    )


def newtonian_gravity(mass_kg, radius_m):
    G_NEWTON = 6.67430e-11

    if radius_m <= 0:
        raise ValueError("Il raggio deve essere positivo")

    return float(G_NEWTON * mass_kg / radius_m**2)


def psu_modulation(
    distance_m,
    amplitude,
    ell_obs_m,
    phase_rad=0.0,
):
    if ell_obs_m <= 0:
        raise ValueError("ell_obs deve essere positivo")

    return float(
        1.0
        + amplitude
        * np.cos(
            2.0 * np.pi * distance_m / ell_obs_m
            + phase_rad
        )
    )


def power_spectrum(frequencies_hz, amplitude, beta, white_noise):
    frequencies_hz = np.asarray(frequencies_hz, dtype=float)

    if np.any(frequencies_hz <= 0):
        raise ValueError("Le frequenze devono essere positive")

    return (
        white_noise
        + amplitude * frequencies_hz ** (-beta)
    )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Cinque test fenomenologici per un modello PSU. "
            "Nessun parametro viene presentato come misurato."
        )
    )

    parser.add_argument(
        "--outdir",
        default="psu_physical_tests",
    )

    parser.add_argument(
        "--redshift",
        type=float,
        default=0.903,
        help="Redshift di esempio, GRB 090510.",
    )

    parser.add_argument(
        "--E-high-GeV",
        type=float,
        default=31.0,
    )

    parser.add_argument(
        "--E-low-GeV",
        type=float,
        default=0.1,
    )

    parser.add_argument(
        "--delay-limit-s",
        type=float,
        default=1.0,
        help="Limite conservativo di ritardo in secondi.",
    )

    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    E_high_eV = args.E_high_GeV * 1.0e9
    E_low_eV = args.E_low_GeV * 1.0e9

    # ============================================================
    # TEST 1: dispersione fotonica
    # ============================================================

    dispersion_rows = []

    for n in [1, 2]:
        E_qg = PLANCK_ENERGY_EV

        predicted_delay = photon_delay_seconds(
            redshift=args.redshift,
            E_high_eV=E_high_eV,
            E_low_eV=E_low_eV,
            E_qg_eV=E_qg,
            n=n,
            xi=1.0,
        )

        bound = qg_bound_from_delay(
            redshift=args.redshift,
            E_high_eV=E_high_eV,
            E_low_eV=E_low_eV,
            delay_limit_seconds=args.delay_limit_s,
            n=n,
            xi=1.0,
        )

        dispersion_rows.append({
            "test": "photon_dispersion",
            "n": n,
            "redshift": args.redshift,
            "E_high_GeV": args.E_high_GeV,
            "E_low_GeV": args.E_low_GeV,
            "E_QG_assumed_eV": E_qg,
            "E_QG_assumed_over_EPlanck": (
                E_qg / PLANCK_ENERGY_EV
            ),
            "xi": 1.0,
            "predicted_delay_s": predicted_delay,
            "input_delay_limit_s": args.delay_limit_s,
            "derived_E_QG_lower_bound_eV": bound,
            "derived_E_QG_lower_bound_over_EPlanck": (
                bound / PLANCK_ENERGY_EV
            ),
            "status": "PHENOMENOLOGICAL_TEST",
        })

    dispersion_df = pd.DataFrame(dispersion_rows)
    dispersion_df.to_csv(
        outdir / "test_1_photon_dispersion.csv",
        index=False,
    )

    # ============================================================
    # TEST 2: anisotropia/Lorentz violation
    # ============================================================

    Q_isotropic = np.eye(3) / 3.0

    Q_anisotropic = np.diag([
        0.50,
        0.30,
        0.20,
    ])

    Q_anisotropic_traceless = (
        Q_anisotropic
        - np.eye(3) * np.trace(Q_anisotropic) / 3.0
    )

    directions = {
        "x": np.array([1.0, 0.0, 0.0]),
        "y": np.array([0.0, 1.0, 0.0]),
        "z": np.array([0.0, 0.0, 1.0]),
        "diagonal": np.array([1.0, 1.0, 1.0]),
    }

    lv_rows = []

    tensor_cases = {
        "isotropic": Q_isotropic,
        "anisotropic": Q_anisotropic,
    }

    for tensor_name, Q in tensor_cases.items():
        metrics = lorentz_anisotropy_tensor(Q)

        for direction_name, direction in directions.items():
            factor = lv_direction_factor(
                direction,
                Q_anisotropic_traceless
                if tensor_name == "anisotropic"
                else np.zeros((3, 3)),
            )

            lv_rows.append({
                "test": "lorentz_anisotropy",
                "tensor_case": tensor_name,
                "direction": direction_name,
                "lambda_min": metrics["lambda_min"],
                "lambda_mid": metrics["lambda_mid"],
                "lambda_max": metrics["lambda_max"],
                "anisotropy": metrics["anisotropy"],
                "direction_factor": factor,
                "xi": np.nan,
                "status": "MODEL_GEOMETRY_ONLY",
            })

    pd.DataFrame(lv_rows).to_csv(
        outdir / "test_2_lorentz_anisotropy.csv",
        index=False,
    )

    # ============================================================
    # TEST 3: densità di energia QVSM/PSU
    # ============================================================

    # Esempio esplicito: non è un valore PSU misurato.
    probabilities = np.array([
        0.25,
        0.25,
        0.25,
        0.25,
    ])

    energies_J = np.array([
        0.0,
        1.0e-20,
        2.0e-20,
        3.0e-20,
    ])

    volume_m3 = 1.0e-30

    rho = psu_vacuum_density(
        probabilities,
        energies_J,
        volume_m3,
    )

    vacuum_row = pd.DataFrame([{
        "test": "QVSM_vacuum_density",
        "volume_m3": volume_m3,
        "rho_QVSM_J_m3": rho,
        "rho_QFT_ren_J_m3": np.nan,
        "delta_rho_J_m3": np.nan,
        "Lambda_shift": np.nan,
        "status": "SYNTHETIC_PARAMETERIZATION_ONLY",
        "notes": (
            "Requires a physically defined PSU Hamiltonian and "
            "gravitational coupling before interpretation."
        ),
    }])

    vacuum_row.to_csv(
        outdir / "test_3_vacuum_density.csv",
        index=False,
    )

    # ============================================================
    # TEST 4: scala discreta/modulazione
    # ============================================================

    distances_m = np.linspace(
        0.0,
        100.0 * PLANCK_LENGTH_M,
        1001,
    )

    ell_obs_values = [
        PLANCK_LENGTH_M,
        10.0 * PLANCK_LENGTH_M,
        100.0 * PLANCK_LENGTH_M,
    ]

    modulation_rows = []

    for ell_obs in ell_obs_values:
        modulation = np.array([
            psu_modulation(
                distance_m=d,
                amplitude=1.0e-6,
                ell_obs_m=ell_obs,
                phase_rad=0.0,
            )
            for d in distances_m
        ])

        for distance, value in zip(distances_m, modulation):
            modulation_rows.append({
                "test": "discrete_scale_modulation",
                "distance_m": distance,
                "ell_obs_m": ell_obs,
                "amplitude": 1.0e-6,
                "modulation": value,
                "status": "HYPOTHESIS_PARAMETERIZATION",
            })

    pd.DataFrame(modulation_rows).to_csv(
        outdir / "test_4_discrete_scale_modulation.csv",
        index=False,
    )

    # ============================================================
    # TEST 5: massa/gravità da occupazioni PSU
    # ============================================================

    # Esempio minimo a tre microstati.
    # m_star non è la massa Planck automaticamente:
    # è un parametro libero dichiarato del modello.
    m_star_kg = 1.0e-30

    occupancies = np.array([
        10.0,
        20.0,
        30.0,
    ])

    weights = np.array([
        1.0,
        2.0,
        0.5,
    ])

    mass = newtonian_mass_from_psu(
        m_star_kg,
        occupancies,
        weights,
    )

    radii_m = np.array([
        1.0e-10,
        1.0e-9,
        1.0e-8,
    ])

    gravity_rows = []

    for radius in radii_m:
        gravity_rows.append({
            "test": "PSU_mass_gravity",
            "m_star_kg": m_star_kg,
            "weighted_occupancy": np.sum(
                occupancies * weights
            ),
            "mass_QVSM_kg": mass,
            "radius_m": radius,
            "g_newtonian_m_s2": newtonian_gravity(
                mass,
                radius,
            ),
            "eta_correction": np.nan,
            "status": "SYNTHETIC_PARAMETERIZATION_ONLY",
            "notes": (
                "Requires an independently derived PSU mass rule "
                "and measured inertial/gravitational comparison."
            ),
        })

    pd.DataFrame(gravity_rows).to_csv(
        outdir / "test_5_mass_gravity.csv",
        index=False,
    )

    # ============================================================
    # TEST 6: firma EM spettrale
    # ============================================================

    frequencies_hz = np.logspace(-3, 6, 1000)

    spectrum_rows = []

    for beta in [0.0, 1.0, 2.0]:
        spectrum = power_spectrum(
            frequencies_hz,
            amplitude=1.0e-12,
            beta=beta,
            white_noise=1.0e-15,
        )

        for frequency, value in zip(
            frequencies_hz,
            spectrum,
        ):
            spectrum_rows.append({
                "test": "EM_spectral_signature",
                "frequency_Hz": frequency,
                "beta": beta,
                "amplitude": 1.0e-12,
                "white_noise": 1.0e-15,
                "S_E": value,
                "status": "HYPOTHESIS_PARAMETERIZATION",
            })

    pd.DataFrame(spectrum_rows).to_csv(
        outdir / "test_6_em_spectral_signature.csv",
        index=False,
    )

    # ============================================================
    # MANIFEST
    # ============================================================

    manifest = {
        "project": "PSU-physical-tests",
        "status": (
            "phenomenological framework; no PSU parameter is "
            "treated as experimentally established"
        ),
        "tests": [
            "photon dispersion",
            "Lorentz anisotropy",
            "QVSM vacuum density",
            "discrete scale modulation",
            "PSU mass/gravity",
            "EM spectral signature",
        ],
        "experimental_constraint": {
            "GRB090510_linear": "E_QG_1 > 7.6 E_Planck",
            "GRB090510_quadratic": (
                "E_QG_2 > 1.3e11 GeV"
            ),
            "confidence": "95 percent in cited Fermi analysis",
        },
        "fixed_constants": {
            "Planck_energy_eV": PLANCK_ENERGY_EV,
            "Planck_length_m": PLANCK_LENGTH_M,
            "H0_km_s_Mpc": H0_KM_S_MPC,
            "Omega_m": OMEGA_M,
            "Omega_lambda": OMEGA_L,
        },
        "epistemic_rules": [
            "Planck units are dimensional combinations, not evidence of PSU ontology.",
            "ZPE is not a unique PSU signature.",
            "A PSU claim requires a distinct field/equation and a measurable prediction.",
            "Synthetic parameterizations are not experimental measurements.",
            "No fit to data is performed by this script.",
        ],
        "references": [
            "Fermi GRB LIV constraints: 1305.3463 / 10.1016/j.nima.2013.10.088",
            "SME framework: hep-ph/9809521",
        ],
    }

    manifest_path = outdir / "psu_test_manifest.json"

    manifest_path.write_text(
        json.dumps(
            manifest,
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )

    print("Creati i test PSU in:", outdir)
    print("File principali:")
    print("  test_1_photon_dispersion.csv")
    print("  test_2_lorentz_anisotropy.csv")
    print("  test_3_vacuum_density.csv")
    print("  test_4_discrete_scale_modulation.csv")
    print("  test_5_mass_gravity.csv")
    print("  test_6_em_spectral_signature.csv")
    print("  psu_test_manifest.json")


if __name__ == "__main__":
    main()

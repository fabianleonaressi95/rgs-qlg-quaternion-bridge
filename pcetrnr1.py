#!/usr/bin/env python3

import csv
import math
import sys
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser


DEFAULT_STRUCTURE_FILE = Path("6w4x.cif")
DEFAULT_EDGE_FILE = Path(
    "rnr_hybrid_results/rnr_hybrid_edges.csv"
)
OUTPUT_DIR = Path("rnr_pcet_results")

TEMPERATURE_K = 298.15
KB_EV_PER_K = 8.617333262e-5
HBAR_EV_S = 6.582119569e-16
AMU_KG = 1.66053906660e-27
EV_TO_J = 1.602176634e-19
ANGSTROM_M = 1.0e-10

# Parametri iniziali espliciti.
# Sono proxy sostituibili, non valori sperimentali universali.
BETA_ET_PER_A = 1.0
HAB0_EV = 0.010
LAMBDA_EV = 0.70
DELTA_G_ET_EV = 0.0

# Parametri del modello WKB protonico.
PROTON_BARRIER_EV = 0.40
PROTON_BARRIER_WIDTH_A = 0.60
ATTEMPT_FREQUENCY_S = 1.0e13

HBOND_MAX_A = 3.6
HETERO_ELEMENTS = {"N", "O", "S"}


def safe_float(value):
    if value is None:
        return None

    text = str(value).strip()

    if text == "":
        return None

    if text.upper() in {"NA", "N/A", ".", "NONE"}:
        return None

    try:
        return float(text)
    except ValueError:
        return None


def fmt(value):
    if value is None:
        return ""

    return f"{float(value):.8g}"


def atom_element(atom):
    element = getattr(atom, "element", "")

    if element:
        return element.strip().upper()

    name = atom.get_name().strip().upper()

    if not name:
        return ""

    return name[0]


def is_hydrogen(atom):
    return atom_element(atom) == "H"


def atom_distance(atom_a, atom_b):
    return float(np.linalg.norm(atom_a.coord - atom_b.coord))


def residue_key(residue):
    chain = residue.get_parent().id
    _, number, _ = residue.id
    name = str(residue.resname).strip().upper()

    return chain, int(number), name


def load_residues(structure):
    residues = {}

    for model in structure:
        for chain in model:
            for residue in chain:
                key = residue_key(residue)
                residues[key] = residue

    return residues


def parse_residue_label(label):
    """
    Converte:
        A:TYR:122
        C:TRP:48
        B:CYS:439

    in:
        chain, numero, resname
    """
    if not label:
        return None

    fields = str(label).split(":")

    if len(fields) < 3:
        return None

    chain = fields[0].strip()
    resname = fields[1].strip().upper()

    number_text = "".join(
        character
        for character in fields[2]
        if character.isdigit()
    )

    if not number_text:
        return None

    return chain, int(number_text), resname


def resolve_residue(residues, label):
    parsed = parse_residue_label(label)

    if parsed is None:
        return None

    chain, number, resname = parsed

    exact_key = (
        chain,
        number,
        resname,
    )

    if exact_key in residues:
        return residues[exact_key]

    # Fallback per residui modificati o nomenclature diverse.
    for key, residue in residues.items():
        candidate_chain, candidate_number, _ = key

        if (
            candidate_chain == chain
            and candidate_number == number
        ):
            return residue

    return None


def heavy_atoms(residue):
    return [
        atom
        for atom in residue.get_atoms()
        if not is_hydrogen(atom)
    ]


def count_hbond_like_contacts(residue_a, residue_b):
    """
    Proxy geometrico per contatti N/O/S entro 3.6 Å.

    Non è un conteggio definitivo dei legami H:
    mancano idrogeni espliciti e criteri angolari.
    """
    count = 0

    for atom_a in heavy_atoms(residue_a):
        element_a = atom_element(atom_a)

        if element_a not in HETERO_ELEMENTS:
            continue

        for atom_b in heavy_atoms(residue_b):
            element_b = atom_element(atom_b)

            if element_b not in HETERO_ELEMENTS:
                continue

            if atom_distance(atom_a, atom_b) <= HBOND_MAX_A:
                count += 1

    return count


def proton_distance_proxy(residue_a, residue_b):
    """
    Minima distanza tra eteroatomi N/O/S.
    """
    best = None

    for atom_a in heavy_atoms(residue_a):
        if atom_element(atom_a) not in HETERO_ELEMENTS:
            continue

        for atom_b in heavy_atoms(residue_b):
            if atom_element(atom_b) not in HETERO_ELEMENTS:
                continue

            distance = atom_distance(atom_a, atom_b)

            if best is None or distance < best:
                best = distance

    return best


def get_d_et_proxy(row):
    centroid_distance = safe_float(
        row.get("centroid_distance_A")
    )

    if centroid_distance is not None:
        return centroid_distance

    return safe_float(
        row.get("min_distance_A")
    )


def h_ab_proxy(distance_A):
    """
    Proxy:
        H_AB = H0 exp[-beta (d - d0)]
    """
    if distance_A is None:
        return None

    reference_distance_A = 3.5

    return HAB0_EV * math.exp(
        -BETA_ET_PER_A
        * max(0.0, distance_A - reference_distance_A)
    )


def marcus_rate(
    h_ab_ev,
    lambda_ev,
    delta_g_ev,
):
    """
    Tasso di Marcus nonadiabatico semplificato.
    """
    if h_ab_ev is None:
        return None

    if lambda_ev <= 0:
        return None

    thermal_energy = (
        KB_EV_PER_K
        * TEMPERATURE_K
    )

    prefactor = (
        2.0
        * math.pi
        / HBAR_EV_S
    )

    gaussian_prefactor = math.sqrt(
        math.pi
        / (
            lambda_ev
            * thermal_energy
        )
    )

    exponent = -(
        (delta_g_ev + lambda_ev) ** 2
        / (
            4.0
            * lambda_ev
            * thermal_energy
        )
    )

    return (
        prefactor
        * h_ab_ev ** 2
        * gaussian_prefactor
        * math.exp(exponent)
    )


def wkb_exponent(
    mass_amu,
    barrier_ev,
    width_A,
):
    """
    Esponente WKB:
        2 a sqrt(2 m V) / hbar
    """
    mass_kg = mass_amu * AMU_KG
    barrier_joule = barrier_ev * EV_TO_J
    width_m = width_A * ANGSTROM_M
    hbar_joule_s = HBAR_EV_S * EV_TO_J

    return (
        2.0
        * width_m
        * math.sqrt(
            2.0
            * mass_kg
            * barrier_joule
        )
        / hbar_joule_s
    )


def tunnelling_rate(mass_amu):
    exponent = wkb_exponent(
        mass_amu,
        PROTON_BARRIER_EV,
        PROTON_BARRIER_WIDTH_A,
    )

    return (
        ATTEMPT_FREQUENCY_S
        * math.exp(-exponent)
    )


def calculate_row(row, residues):
    source_label = row.get(
        "source_label",
        "",
    )

    target_label = row.get(
        "target_label",
        "",
    )

    source_residue = resolve_residue(
        residues,
        source_label,
    )

    target_residue = resolve_residue(
        residues,
        target_label,
    )

    d_et = get_d_et_proxy(row)
    h_ab = h_ab_proxy(d_et)

    k_et = marcus_rate(
        h_ab,
        LAMBDA_EV,
        DELTA_G_ET_EV,
    )

    k_h = tunnelling_rate(1.00784)
    k_d = tunnelling_rate(2.01410)
    kie = k_h / k_d

    if (
        source_residue is not None
        and target_residue is not None
    ):
        d_h = proton_distance_proxy(
            source_residue,
            target_residue,
        )

        n_hbond = count_hbond_like_contacts(
            source_residue,
            target_residue,
        )
    else:
        d_h = None
        n_hbond = None

    result = dict(row)

    result.update({
        "d_ET_proxy_A": d_et,
        "H_AB_proxy_eV": h_ab,
        "d_H_proxy_A": d_h,
        "N_Hbond_proxy": n_hbond,
        "lambda_proxy_eV": LAMBDA_EV,
        "deltaG_ET_proxy_eV": DELTA_G_ET_EV,
        "k_ET_proxy_s-1": k_et,
        "k_H_proxy_s-1": k_h,
        "k_D_proxy_s-1": k_d,
        "KIE_proxy": kie,
        "status": (
            "geometry plus Marcus/WKB proxy; "
            "not experimental"
        ),
    })

    return result


def read_edge_rows(edge_file):
    with edge_file.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_results(results):
    OUTPUT_DIR.mkdir(exist_ok=True)

    output_file = (
        OUTPUT_DIR
        / "rnr_pcet_proxy_results.csv"
    )

    columns = [
        "source",
        "target",
        "source_chain",
        "target_chain",
        "geometry_source",
        "source_label",
        "target_label",
        "ca_distance_A",
        "min_distance_A",
        "centroid_distance_A",
        "d_ET_proxy_A",
        "H_AB_proxy_eV",
        "d_H_proxy_A",
        "N_Hbond_proxy",
        "lambda_proxy_eV",
        "deltaG_ET_proxy_eV",
        "k_ET_proxy_s-1",
        "k_H_proxy_s-1",
        "k_D_proxy_s-1",
        "KIE_proxy",
        "status",
    ]

    with output_file.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            extrasaction="ignore",
        )

        writer.writeheader()

        for result in results:
            row = {}

            for column in columns:
                value = result.get(column)

                if column in {
                    "source",
                    "target",
                    "source_chain",
                    "target_chain",
                    "geometry_source",
                    "source_label",
                    "target_label",
                    "status",
                }:
                    row[column] = (
                        ""
                        if value is None
                        else str(value)
                    )
                else:
                    row[column] = fmt(value)

            writer.writerow(row)

    return output_file


def print_results(results):
    print()
    print("Risultati PCET/tunnelling proxy:")
    print()

    for result in results:
        print(
            f"{result.get('source', '')} -> "
            f"{result.get('target', '')} "
            f"("
            f"{result.get('geometry_source', '')}"
            f")"
        )

        print(
            "  d_ET proxy: "
            f"{fmt(result.get('d_ET_proxy_A'))} Å"
        )

        print(
            "  H_AB proxy: "
            f"{fmt(result.get('H_AB_proxy_eV'))} eV"
        )

        print(
            "  d_H proxy: "
            f"{fmt(result.get('d_H_proxy_A'))} Å"
        )

        print(
            "  N_Hbond proxy: "
            f"{result.get('N_Hbond_proxy', '')}"
        )

        print(
            "  k_ET proxy: "
            f"{fmt(result.get('k_ET_proxy_s-1'))} s^-1"
        )

        print(
            "  k_H proxy: "
            f"{fmt(result.get('k_H_proxy_s-1'))} s^-1"
        )

        print(
            "  k_D proxy: "
            f"{fmt(result.get('k_D_proxy_s-1'))} s^-1"
        )

        print(
            "  KIE proxy: "
            f"{fmt(result.get('KIE_proxy'))}"
        )

        print()


def main():
    edge_file = (
        Path(sys.argv[1])
        if len(sys.argv) > 1
        else DEFAULT_EDGE_FILE
    )

    structure_file = (
        Path(sys.argv[2])
        if len(sys.argv) > 2
        else DEFAULT_STRUCTURE_FILE
    )

    if not edge_file.exists():
        raise FileNotFoundError(
            f"File degli archi non trovato: "
            f"{edge_file}"
        )

    if not structure_file.exists():
        raise FileNotFoundError(
            f"File della struttura non trovato: "
            f"{structure_file}"
        )

    parser = MMCIFParser(QUIET=True)

    structure = parser.get_structure(
        "6W4X",
        str(structure_file),
    )

    residues = load_residues(structure)
    edge_rows = read_edge_rows(edge_file)

    results = [
        calculate_row(row, residues)
        for row in edge_rows
    ]

    output_file = write_results(results)

    print_results(results)

    print(
        f"File prodotto: {output_file}"
    )


if __name__ == "__main__":
    main()

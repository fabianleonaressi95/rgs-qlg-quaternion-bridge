#!/usr/bin/env python3

import csv
import math
import sys
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser


STRUCTURE_FILE = Path("6w4x.cif")
EDGE_FILE = Path(
    "rnr_hybrid_results/rnr_hybrid_edges.csv"
)
OUTPUT_DIR = Path("rnr_pcet_results")

TEMPERATURE_K = 298.15
KB_EV_PER_K = 8.617333262e-5
H_EV_S = 4.135667696e-15
HBAR_EV_S = 6.582119569e-16
AMU_KG = 1.66053906660e-27
EV_TO_J = 1.602176634e-19
ANGSTROM_M = 1.0e-10

# Parametri geometrici/prototipali.
# Sono espliciti e sostituibili, non dati sperimentali automatici.
BETA_ET_PER_A = 1.0
HAB0_EV = 0.010
LAMBDA_EV = 0.70
DELTA_G_ET_EV = 0.0

# Barriera protonica proxy.
PROTON_BARRIER_EV = 0.40
PROTON_BARRIER_WIDTH_A = 0.60

# Legame H: distanza donatore-accettore pesante.
HBOND_MAX_A = 3.6

# Solo contatti con questi elementi possono fungere
# da candidati donatore/accettore protonico.
H_ACCEPTOR_ELEMENTS = {"N", "O", "S"}


def safe_float(value):
    if value is None:
        return None

    text = str(value).strip()

    if not text or text.upper() in {"NA", "N/A", "."}:
        return None

    return float(text)


def fmt(value):
    if value is None:
        return ""

    return f"{float(value):.8g}"


def read_edges():
    if not EDGE_FILE.exists():
        raise FileNotFoundError(
            f"File non trovato: {EDGE_FILE}"
        )

    with EDGE_FILE.open(newline="") as handle:
        return list(csv.DictReader(handle))


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


def residue_number(residue):
    _, number, _ = residue.id
    return int(number)


def residue_chain(residue):
    return residue.get_parent().id


def residue_key(residue):
    return (
        residue_chain(residue),
        residue_number(residue),
        str(residue.resname).strip().upper(),
    )


def load_residues(structure):
    residues = {}

    for model in structure:
        for chain in model:
            for residue in chain:
                key = residue_key(residue)
                residues[key] = residue

    return residues


def parse_label(label):
    """
    Legge label come:
        A:TYR:122
        C:TRP:48
        B:CYS:439
    """
    if not label:
        return None

    fields = str(label).split(":")

    if len(fields) < 3:
        return None

    chain = fields[0]
    resname = fields[1]

    number_text = "".join(
        character
        for character in fields[2]
        if character.isdigit()
    )

    if not number_text:
        return None

    return chain, int(number_text), resname


def resolve_residue(residues, label):
    parsed = parse_label(label)

    if parsed is None:
        return None

    chain, number, resname = parsed

    key = (
        chain,
        number,
        resname.upper(),
    )

    if key in residues:
        return residues[key]

    # Fallback: stessa catena e numerazione, ignorando
    # la nomenclatura del residuo modificato.
    for candidate_key, residue in residues.items():
        candidate_chain, candidate_number, _ = candidate_key

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


def is_backbone(atom):
    return atom.get_name().strip().upper() in {
        "N",
        "CA",
        "C",
        "O",
    }


def count_hbond_contacts(residue_a, residue_b):
    """
    Conta contatti geometrici N/O/S entro 3.6 Å.
    È un proxy: senza idrogeni espliciti e angoli,
    non è una definizione completa di legame H.
    """
    count = 0

    for atom_a in heavy_atoms(residue_a):
        element_a = atom_element(atom_a)

        if element_a not in H_ACCEPTOR_ELEMENTS:
            continue

        for atom_b in heavy_atoms(residue_b):
            element_b = atom_element(atom_b)

            if element_b not in H_ACCEPTOR_ELEMENTS:
                continue

            distance = atom_distance(atom_a, atom_b)

            if distance <= HBOND_MAX_A:
                count += 1

    return count


def proton_distance_proxy(residue_a, residue_b):
    """
    Cerca la minima distanza tra possibili eteroatomi
    proton-donor/proton-acceptor: N, O, S.
    """
    best = None

    for atom_a in heavy_atoms(residue_a):
        if atom_element(atom_a) not in H_ACCEPTOR_ELEMENTS:
            continue

        for atom_b in heavy_atoms(residue_b):
            if atom_element(atom_b) not in H_ACCEPTOR_ELEMENTS:
                continue

            distance = atom_distance(atom_a, atom_b)

            if best is None or distance < best:
                best = distance

    return best


def inferred_centroid_distance(row):
    value = safe_float(
        row.get("centroid_distance_A")
    )

    return value


def d_et_proxy(row):
    """
    Usa la distanza tra centroidi aromatici se presente.
    """
    value = inferred_centroid_distance(row)

    if value is not None:
        return value

    return safe_float(
        row.get("min_distance_A")
    )


def hab_proxy(distance_A):
    """
    Proxy esponenziale per l'accoppiamento elettronico:
        H_AB = H0 exp[-beta (d - d0)]
    """
    if distance_A is None:
        return None

    d0 = 3.5

    return HAB0_EV * math.exp(
        -BETA_ET_PER_A * max(0.0, distance_A - d0)
    )


def marcus_rate(hab_ev, lambda_ev, delta_g_ev):
    """
    Formula di Marcus nonadiabatica semplificata.
    """
    if hab_ev is None:
        return None

    if lambda_ev <= 0:
        return None

    thermal_factor = math.sqrt(
        math.pi
        / (
            lambda_ev
            * KB_EV_PER_K
            * TEMPERATURE_K
        )
    )

    exponent = -(
        (delta_g_ev + lambda_ev) ** 2
        / (
            4.0
            * lambda_ev
            * KB_EV_PER_K
            * TEMPERATURE_K
        )
    )

    return (
        2.0
        * math.pi
        / HBAR_EV_S
        * hab_ev ** 2
        * thermal_factor
        * math.exp(exponent)
    )


def wkb_exponent(mass_amu, barrier_ev, width_A):
    """
    Calcola l'esponente WKB:
        2 a sqrt(2 m V) / hbar
    """
    mass_kg = mass_amu * AMU_KG
    barrier_joule = barrier_ev * EV_TO_J
    width_m = width_A * ANGSTROM_M

    return (
        2.0
        * width_m
        * math.sqrt(
            2.0
            * mass_kg
            * barrier_joule
        )
        / (HBAR_EV_S * EV_TO_J)
    )


def tunneling_rate(mass_amu):
    exponent = wkb_exponent(
        mass_amu,
        PROTON_BARRIER_EV,
        PROTON_BARRIER_WIDTH_A,
    )

    # Frequenza di tentativo proxy.
    attempt_frequency = 1.0e13

    return attempt_frequency * math.exp(-exponent)


def calculate_edge(row, residues):
    source_label = row.get("source_label", "")
    target_label = row.get("target_label", "")

    source_residue = resolve_residue(
        residues,
        source_label,
    )

    target_residue = resolve_residue(
        residues,
        target_label,
    )

    d_et = d_et_proxy(row)
    hab = hab_proxy(d_et)

    ket = marcus_rate(
        hab,
        LAMBDA_EV,
        DELTA_G_ET_EV,
    )

    kh = tunneling_rate(1.00784)
    kd = tunneling_rate(2.01410)
    kie = kh / kd

    if (
        source_residue is not None
        and target_residue is not None
    ):
        d_h = proton_distance_proxy(
            source_residue,
            target_residue,
        )

        hbond_count = count_hbond_contacts(
            source_residue,
            target_residue,
        )
    else:
        d_h = None
        hbond_count = None

    result = dict(row)

    result.update({
        "d_ET_proxy_A": d_et,
        "H_AB_proxy_eV": hab,
        "d_H_proxy_A": d_h,
        "N_Hbond_proxy": hbond_count,
        "lambda_proxy_eV": LAMBDA_EV,
        "deltaG_ET_proxy_eV": DELTA_G_ET_EV,
        "k_ET_proxy_s-1": ket,
        "k_H_proxy_s-1": kh,
        "k_D_proxy_s-1": kd,
        "KIE_proxy": kie,
        "status": (
            "proxy geometry + Marcus/WKB; "
            "not experimental"
        ),
    })

    return result


def write_output(rows):
    OUTPUT_DIR.mkdir(exist_ok=True)

    output = OUTPUT_DIR / "rnr_pcet_proxy_results.csv"

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

    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            extrasaction="ignore",
        )

        writer.writeheader()

        for row in rows:
            clean = {}

            for column in columns:
                clean[column] = fmt(
                    row.get(column)
                )

            for text_column in [
                "source",
                "target",
                "source_chain",
                "target_chain",
                "geometry_source",
                "source_label",
                "target_label",
                "status",
            ]:
                clean[text_column] = row.get(
                    text_column,
                    "",
                )

            writer.writerow(clean)

    return output


def print_results(rows):
    print()
    print("Risultati PCET/tunnelling proxy:")

    for row in rows:
        print(
            f"{row['source']} -> {row['target']} "
            f"({row.get('geometry_source', '')})"
        )

        print(
            "  d_ET proxy: "
            f"{fmt(row.get('d_ET_proxy_A'))} Å"
        )

        print(
            "  H_AB proxy: "
            f"{fmt(row.get('H_AB_proxy_eV'))} eV"
        )

        print(
            "  d_H proxy: "
            f"{fmt(row.get('d_H_proxy_A'))} Å"
        )

        print(
            "  N_Hbond proxy: "
            f"{row.get('N_Hbond_proxy', '')}"
        )

        print(
            "  k_ET proxy: "
            f"{fmt(row.get('k_ET_proxy_s-1'))} s^-1"
        )

        print(
            "  k_H proxy: "
            f"{fmt(row.get('k_H_proxy_s-1'))} s^-1"
        )

        print(
            "  k_D proxy: "
            f"{fmt(row.get('k_D_proxy_s-1'))} s^-1"
        )

        print(
            "  KIE proxy: "
            f"{fmt(row.get('KIE_proxy'))}"
        )

        print()


def main():
    edge_file = (
        Path(sys.argv[1])
        if len(sys.argv) > 1
        else EDGE_FILE
    )

    structure_file = (
        Path(sys.argv[2])
        if len(sys.argv) > 2
        else STRUCTURE_FILE
    )

    global EDGE_FILE
    EDGE_FILE = edge_file

    if not edge_file.exists():
        raise FileNotFoundError(
            f"File archi non trovato: {edge_file}"
        )

    if not structure_file.exists():
        raise FileNotFoundError(
            f"Struttura non trovata: {structure_file}"
        )

    parser = MMCIFParser(QUIET=True)

    structure = parser.get_structure(
        "6W4X",
        str(structure_file),
    )

    residues = load_residues(structure)
    edges = read_edges()

    results = [
        calculate_edge(row, residues)
        for row in edges
    ]

    output = write_output(results)

    print_results(results)
    print(f"File prodotto: {output}")


if __name__ == "__main__":
    main()

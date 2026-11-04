#!/usr/bin/env python3

import csv
import sys
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser


INPUT_DEFAULT = "6w4x.cif"
OUTPUT_DIR = Path("6w4x_results")

MAX_ET_DISTANCE = 8.0
MAX_HBOND_DISTANCE = 3.6


# Mapping strutturale usato nel file 6w4x.cif:
#
# FY122 biologico -> TYR 123 nel file atomico
#
PATHWAY = [
    ("FY122", 123, {"TYR", "FY3"}),
    ("W48", 48, {"TRP"}),
    ("Y356", 356, {"TYR"}),
    ("Y731", 731, {"TYR"}),
    ("Y730", 730, {"TYR"}),
    ("C439", 439, {"CYS"}),
]


AROMATIC_ATOMS = {
    "TYR": [
        "CG",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CZ",
    ],
    "FY3": [
        "CG",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CZ",
    ],
    "TRP": [
        "CG",
        "CD1",
        "CD2",
        "NE1",
        "CE2",
        "CE3",
        "CZ2",
        "CZ3",
        "CH2",
    ],
    "CYS": [
        "CB",
        "SG",
    ],
}


# Residui della catena principale.
BACKBONE_ATOMS = {
    "N",
    "CA",
    "C",
    "O",
}


def normalize(name):
    return str(name).strip().upper()


def is_hydrogen(atom):
    element = getattr(atom, "element", "")

    if element:
        return element.upper() == "H"

    return atom.get_name().upper().startswith("H")


def is_backbone_atom(atom):
    return atom.get_name().upper() in BACKBONE_ATOMS


def residue_label(residue):
    chain = residue.get_parent().id
    _, number, insertion = residue.id
    insertion = str(insertion).strip()

    if insertion:
        return f"{chain}:{residue.resname}:{number}{insertion}"

    return f"{chain}:{residue.resname}:{number}"


def atom_distance(atom_a, atom_b):
    return float(np.linalg.norm(atom_a.coord - atom_b.coord))


def heavy_atoms(residue):
    return [
        atom
        for atom in residue.get_atoms()
        if not is_hydrogen(atom)
    ]


def aromatic_centroid(residue):
    name = normalize(residue.resname)
    atom_names = AROMATIC_ATOMS.get(name, [])

    coords = [
        residue[atom_name].coord
        for atom_name in atom_names
        if atom_name in residue
    ]

    if not coords:
        return None

    return np.mean(np.asarray(coords), axis=0)


def centroid_distance(residue_a, residue_b):
    centroid_a = aromatic_centroid(residue_a)
    centroid_b = aromatic_centroid(residue_b)

    if centroid_a is None or centroid_b is None:
        return None

    return float(
        np.linalg.norm(centroid_a - centroid_b)
    )


def ca_distance(residue_a, residue_b):
    if "CA" not in residue_a:
        return None

    if "CA" not in residue_b:
        return None

    return atom_distance(
        residue_a["CA"],
        residue_b["CA"],
    )


def residue_sequence_distance(residue_a, residue_b):
    chain_a = residue_a.get_parent().id
    chain_b = residue_b.get_parent().id

    _, number_a, _ = residue_a.id
    _, number_b, _ = residue_b.id

    if chain_a != chain_b:
        return None

    return abs(int(number_a) - int(number_b))


def is_peptide_contact(atom_a, atom_b, residue_a, residue_b):
    """
    Elimina i contatti ovvi della catena peptidica:
    - residui consecutivi;
    - atomi backbone-backbone;
    - distanze covalenti C-N.
    """
    sequence_distance = residue_sequence_distance(
        residue_a,
        residue_b,
    )

    if sequence_distance != 1:
        return False

    if not is_backbone_atom(atom_a):
        return False

    if not is_backbone_atom(atom_b):
        return False

    atom_names = {
        atom_a.get_name().upper(),
        atom_b.get_name().upper(),
    }

    if atom_names == {"C", "N"}:
        return True

    return False


def minimum_relevant_distance(residue_a, residue_b):
    """
    Cerca la minima distanza tra atomi pesanti,
    escludendo contatti peptidici.
    """
    best_distance = None
    best_pair = None

    for atom_a in heavy_atoms(residue_a):
        for atom_b in heavy_atoms(residue_b):
            if is_peptide_contact(
                atom_a,
                atom_b,
                residue_a,
                residue_b,
            ):
                continue

            distance = atom_distance(atom_a, atom_b)

            if best_distance is None or distance < best_distance:
                best_distance = distance
                best_pair = (
                    atom_a.get_name(),
                    atom_b.get_name(),
                )

    return best_distance, best_pair


def find_residues(structure):
    """
    Restituisce tutte le copie dei residui del percorso.
    """
    found = []

    for model in structure:
        for chain in model:
            for residue in chain:
                _, number, insertion = residue.id
                name = normalize(residue.resname)

                for node, expected_number, allowed_names in PATHWAY:
                    if number != expected_number:
                        continue

                    if name not in allowed_names:
                        continue

                    found.append({
                        "model": model.id,
                        "chain": chain.id,
                        "node": node,
                        "resname": name,
                        "number": int(number),
                        "icode": str(insertion).strip(),
                        "residue": residue,
                    })

    return found


def node_matches(found, node):
    return [
        item
        for item in found
        if item["node"] == node
    ]


def candidate_score(candidate):
    """
    Punteggio per selezionare l'arco più plausibile.

    Priorità:
    1. distanza tra centroidi aromatici;
    2. distanza minima tra atomi pesanti;
    3. distanza CA-CA.

    Un valore None viene penalizzato.
    """
    d_centroid = candidate["centroid_distance_A"]
    d_min = candidate["min_distance_A"]
    d_ca = candidate["ca_distance_A"]

    centroid_score = (
        d_centroid
        if d_centroid is not None
        else 999.0
    )

    min_score = (
        d_min
        if d_min is not None
        else 999.0
    )

    ca_score = (
        d_ca
        if d_ca is not None
        else 999.0
    )

    return (
        int(candidate["within_8A"] == 0),
        centroid_score,
        min_score,
        ca_score,
    )


def build_candidates(found):
    """
    Costruisce tutte le combinazioni per ogni arco consecutivo.
    """
    candidates = []

    for left, right in zip(PATHWAY, PATHWAY[1:]):
        source_node = left[0]
        target_node = right[0]

        source_items = node_matches(
            found,
            source_node,
        )

        target_items = node_matches(
            found,
            target_node,
        )

        for source in source_items:
            for target in target_items:
                source_residue = source["residue"]
                target_residue = target["residue"]

                d_min, atom_pair = (
                    minimum_relevant_distance(
                        source_residue,
                        target_residue,
                    )
                )

                d_centroid = centroid_distance(
                    source_residue,
                    target_residue,
                )

                d_ca = ca_distance(
                    source_residue,
                    target_residue,
                )

                if d_min is None:
                    within_8 = 0
                else:
                    within_8 = int(
                        d_min <= MAX_ET_DISTANCE
                    )

                if d_centroid is None:
                    within_centroid = 0
                else:
                    within_centroid = int(
                        d_centroid <= MAX_ET_DISTANCE
                    )

                candidates.append({
                    "source": source_node,
                    "target": target_node,
                    "source_chain": source["chain"],
                    "target_chain": target["chain"],
                    "source_label": residue_label(
                        source_residue
                    ),
                    "target_label": residue_label(
                        target_residue
                    ),
                    "ca_distance_A": d_ca,
                    "min_distance_A": d_min,
                    "atom_pair": atom_pair,
                    "centroid_distance_A": d_centroid,
                    "within_8A": within_8,
                    "within_centroid_8A": within_centroid,
                    "same_chain": int(
                        source["chain"] == target["chain"]
                    ),
                })

    return candidates


def select_best_edges(candidates):
    """
    Seleziona una sola combinazione per ogni arco.
    """
    selected = []

    edge_names = [
        ("FY122", "W48"),
        ("W48", "Y356"),
        ("Y356", "Y731"),
        ("Y731", "Y730"),
        ("Y730", "C439"),
    ]

    for source_node, target_node in edge_names:
        edge_candidates = [
            item
            for item in candidates
            if (
                item["source"] == source_node
                and item["target"] == target_node
            )
        ]

        if not edge_candidates:
            continue

        best = sorted(
            edge_candidates,
            key=candidate_score,
        )[0]

        best["selected"] = 1

        selected.append(best)

    return selected


def write_all_candidates(candidates):
    output = OUTPUT_DIR / "rnr_all_edge_candidates.csv"

    with output.open("w", newline="") as handle:
        writer = csv.writer(handle)

        writer.writerow([
            "source",
            "target",
            "source_chain",
            "target_chain",
            "source_label",
            "target_label",
            "same_chain",
            "ca_distance_A",
            "min_distance_A",
            "atom_pair",
            "centroid_distance_A",
            "within_8A",
            "within_centroid_8A",
        ])

        for item in candidates:
            writer.writerow([
                item["source"],
                item["target"],
                item["source_chain"],
                item["target_chain"],
                item["source_label"],
                item["target_label"],
                item["same_chain"],
                format_value(
                    item["ca_distance_A"]
                ),
                format_value(
                    item["min_distance_A"]
                ),
                format_pair(
                    item["atom_pair"]
                ),
                format_value(
                    item["centroid_distance_A"]
                ),
                item["within_8A"],
                item["within_centroid_8A"],
            ])

    return output


def write_selected_edges(selected):
    output = OUTPUT_DIR / "rnr_selected_edges.csv"

    with output.open("w", newline="") as handle:
        writer = csv.writer(handle)

        writer.writerow([
            "source",
            "target",
            "source_chain",
            "target_chain",
            "source_label",
            "target_label",
            "same_chain",
            "ca_distance_A",
            "min_distance_A",
            "atom_pair",
            "centroid_distance_A",
            "within_8A",
            "within_centroid_8A",
            "d_ET_A",
            "d_H_A",
            "N_Hbond",
            "lambda_eV",
            "deltaG_ET_eV",
            "k_H_s-1",
            "k_D_s-1",
            "status",
        ])

        for item in selected:
            writer.writerow([
                item["source"],
                item["target"],
                item["source_chain"],
                item["target_chain"],
                item["source_label"],
                item["target_label"],
                item["same_chain"],
                format_value(
                    item["ca_distance_A"]
                ),
                format_value(
                    item["min_distance_A"]
                ),
                format_pair(
                    item["atom_pair"]
                ),
                format_value(
                    item["centroid_distance_A"]
                ),
                item["within_8A"],
                item["within_centroid_8A"],
                format_value(
                    item["centroid_distance_A"]
                ),
                "",
                "",
                "",
                "",
                "",
                "",
                "geometry extracted; kinetic/QM values missing",
            ])

    return output


def format_value(value):
    if value is None:
        return ""

    return f"{float(value):.4f}"


def format_pair(pair):
    if pair is None:
        return ""

    return f"{pair[0]}-{pair[1]}"


def print_selected_edges(selected):
    print()
    print("Archi selezionati:")

    for item in selected:
        print(
            f"  {item['source']} -> {item['target']} "
            f"({item['source_chain']}"
            f"-{item['target_chain']})"
        )

        print(
            f"    distanza minima: "
            f"{format_value(item['min_distance_A'])} Å"
        )

        print(
            f"    centroide: "
            f"{format_value(item['centroid_distance_A'])} Å"
        )

        print(
            f"    CA-CA: "
            f"{format_value(item['ca_distance_A'])} Å"
        )

        print(
            f"    coppia atomica: "
            f"{format_pair(item['atom_pair'])}"
        )


def main():
    input_file = (
        Path(sys.argv[1])
        if len(sys.argv) > 1
        else Path(INPUT_DEFAULT)
    )

    if not input_file.exists():
        raise FileNotFoundError(
            f"File non trovato: {input_file}"
        )

    OUTPUT_DIR.mkdir(exist_ok=True)

    parser = MMCIFParser(QUIET=True)

    structure = parser.get_structure(
        "6W4X",
        str(input_file),
    )

    found = find_residues(structure)
    candidates = build_candidates(found)
    selected = select_best_edges(candidates)

    all_candidates_file = write_all_candidates(
        candidates
    )

    selected_file = write_selected_edges(
        selected
    )

    detected = sorted({
        item["node"]
        for item in found
    })

    expected = [
        item[0]
        for item in PATHWAY
    ]

    missing = [
        node
        for node in expected
        if node not in detected
    ]

    print(f"Struttura analizzata: {input_file}")
    print(
        f"Residui del percorso trovati: "
        f"{len(found)}"
    )

    print(
        "Nodi rilevati: "
        + (
            ", ".join(detected)
            if detected
            else "nessuno"
        )
    )

    if missing:
        print()
        print("Nodi mancanti:")
        print("  " + ", ".join(missing))

    print_selected_edges(selected)

    print()
    print("File prodotti:")
    print(f"  {all_candidates_file}")
    print(f"  {selected_file}")


if __name__ == "__main__":
    main()

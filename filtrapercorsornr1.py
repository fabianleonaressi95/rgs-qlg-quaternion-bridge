#!/usr/bin/env python3

import csv
import sys
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser


INPUT_DEFAULT = "6w4x.cif"
OUTPUT_DIR = Path("6w4x_results")

CONTACT_CUTOFF = 4.0
MAX_ET_DISTANCE = 12.0
MAX_CENTROID_DISTANCE = 14.0


# Percorso radicalico.
#
# Nel file 6w4x.cif:
#
# FY122 biologico = FY3 strutturale 123
# FY3 è presente nelle catene C e D.
#
PATHWAY = [
    ("FY122", 123),
    ("W48", 48),
    ("Y356", 356),
    ("Y731", 731),
    ("Y730", 730),
    ("C439", 439),
]


EXPECTED_NAMES = {
    "W48": {"TRP"},
    "Y356": {"TYR"},
    "Y731": {"TYR"},
    "Y730": {"TYR"},
    "C439": {"CYS"},
}


# Compatibilità biologica tra copie del complesso.
#
# Il percorso è:
#
# FY122(C/D) -> W48(C/D)
# W48(C/D)   -> Y356(C/D)
# Y356(C/D)  -> Y731(A/B)
# Y731(A/B)  -> Y730(A/B)
# Y730(A/B)  -> C439(A/B)
#
ALLOWED_CHAIN_PAIRS = {
    ("FY122", "W48"): {
        ("C", "C"),
        ("D", "D"),
    },
    ("W48", "Y356"): {
        ("C", "C"),
        ("D", "D"),
    },
    ("Y356", "Y731"): {
        ("C", "A"),
        ("C", "B"),
        ("D", "A"),
        ("D", "B"),
    },
    ("Y731", "Y730"): {
        ("A", "A"),
        ("B", "B"),
    },
    ("Y730", "C439"): {
        ("A", "A"),
        ("B", "B"),
    },
}


AROMATIC_ATOMS = {
    "FY3": [
        "CG",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CZ",
    ],
    "TYR": [
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

    if chain_a != chain_b:
        return None

    _, number_a, _ = residue_a.id
    _, number_b, _ = residue_b.id

    return abs(int(number_a) - int(number_b))


def is_backbone_contact(atom_a, atom_b, residue_a, residue_b):
    """
    Identifica contatti backbone-backbone tra residui consecutivi.

    Questi contatti non vengono usati come candidati ET/PCET.
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

    return atom_names == {"C", "N"}


def minimum_relevant_distance(residue_a, residue_b):
    """
    Minima distanza tra atomi pesanti, escludendo
    i legami peptidici backbone C-N.
    """
    best_distance = None
    best_pair = None
    best_backbone = 0

    for atom_a in heavy_atoms(residue_a):
        for atom_b in heavy_atoms(residue_b):
            backbone = is_backbone_contact(
                atom_a,
                atom_b,
                residue_a,
                residue_b,
            )

            if backbone:
                continue

            distance = atom_distance(atom_a, atom_b)

            if best_distance is None or distance < best_distance:
                best_distance = distance
                best_pair = (
                    atom_a.get_name(),
                    atom_b.get_name(),
                )
                best_backbone = int(backbone)

    return best_distance, best_pair, best_backbone


def aromatic_centroid(residue):
    residue_name = normalize(residue.resname)
    atom_names = AROMATIC_ATOMS.get(residue_name, [])

    coordinates = []

    for atom_name in atom_names:
        if atom_name in residue:
            coordinates.append(residue[atom_name].coord)

    if not coordinates:
        return None

    return np.mean(
        np.asarray(coordinates),
        axis=0,
    )


def centroid_distance(residue_a, residue_b):
    centroid_a = aromatic_centroid(residue_a)
    centroid_b = aromatic_centroid(residue_b)

    if centroid_a is None or centroid_b is None:
        return None

    return float(
        np.linalg.norm(centroid_a - centroid_b)
    )


def add_residue(found, model, chain, node, residue):
    _, number, insertion = residue.id

    found.append({
        "model": model.id,
        "chain": chain.id,
        "node": node,
        "resname": normalize(residue.resname),
        "number": int(number),
        "icode": str(insertion).strip(),
        "residue": residue,
    })


def find_pathway_residues(structure):
    """
    Cerca i residui del percorso.

    FY122:
        accetta SOLO FY3 123 nelle catene C e D.

    Questo evita di confondere:
        TYR A/B 123
    con:
        FY3 C/D 123.
    """
    found = []

    for model in structure:
        for chain in model:
            chain_id = chain.id

            for residue in chain:
                _, number, _ = residue.id
                name = normalize(residue.resname)

                # FY122 biologico:
                # nel percorso PCET è FY3 123 C/D.
                if (
                    number == 123
                    and name == "FY3"
                    and chain_id in {"C", "D"}
                ):
                    add_residue(
                        found,
                        model,
                        chain,
                        "FY122",
                        residue,
                    )
                    continue

                # Nodi standard.
                for node, expected_names in EXPECTED_NAMES.items():
                    expected_number = dict(PATHWAY)[node]

                    if (
                        number == expected_number
                        and name in expected_names
                    ):
                        add_residue(
                            found,
                            model,
                            chain,
                            node,
                            residue,
                        )

    return found


def node_items(found, node):
    return [
        item
        for item in found
        if item["node"] == node
    ]


def compatible_chain_pair(
    source_node,
    target_node,
    source_chain,
    target_chain,
):
    allowed = ALLOWED_CHAIN_PAIRS.get(
        (source_node, target_node),
        set(),
    )

    return (
        source_chain,
        target_chain,
    ) in allowed


def build_candidates(found):
    """
    Costruisce i candidati rispettando le compatibilità
    di catena definite sopra.
    """
    candidates = []

    for left, right in zip(PATHWAY, PATHWAY[1:]):
        source_node = left[0]
        target_node = right[0]

        source_items = node_items(
            found,
            source_node,
        )

        target_items = node_items(
            found,
            target_node,
        )

        for source in source_items:
            for target in target_items:
                if not compatible_chain_pair(
                    source_node,
                    target_node,
                    source["chain"],
                    target["chain"],
                ):
                    continue

                source_residue = source["residue"]
                target_residue = target["residue"]

                d_min, atom_pair, backbone = (
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

                within_min = int(
                    d_min is not None
                    and d_min <= MAX_ET_DISTANCE
                )

                within_centroid = int(
                    d_centroid is not None
                    and d_centroid <= MAX_CENTROID_DISTANCE
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
                    "within_8A": within_min,
                    "within_centroid_14A": within_centroid,
                    "backbone_contact": backbone,
                })

    return candidates


def candidate_score(item):
    """
    Punteggio di selezione.

    Priorità:
    1. distanza minima entro soglia;
    2. centroide entro soglia;
    3. distanza tra centroidi;
    4. distanza minima tra atomi pesanti;
    5. CA-CA.

    I contatti backbone sono già rimossi dal calcolo
    della distanza minima.
    """
    d_centroid = item["centroid_distance_A"]
    d_min = item["min_distance_A"]
    d_ca = item["ca_distance_A"]

    centroid_value = (
        d_centroid
        if d_centroid is not None
        else 999.0
    )

    min_value = (
        d_min
        if d_min is not None
        else 999.0
    )

    ca_value = (
        d_ca
        if d_ca is not None
        else 999.0
    )

    return (
        0 if item["within_8A"] else 1,
        0 if item["within_centroid_14A"] else 1,
        centroid_value,
        min_value,
        ca_value,
    )


def select_best_edges(candidates):
    """
    Seleziona il candidato migliore per ogni arco.
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


def value_to_csv(value):
    if value is None:
        return ""

    return f"{float(value):.4f}"


def pair_to_csv(pair):
    if pair is None:
        return ""

    return f"{pair[0]}-{pair[1]}"


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
            "ca_distance_A",
            "min_distance_A",
            "atom_pair",
            "centroid_distance_A",
            "within_8A",
            "within_centroid_14A",
            "backbone_contact",
        ])

        for item in candidates:
            writer.writerow([
                item["source"],
                item["target"],
                item["source_chain"],
                item["target_chain"],
                item["source_label"],
                item["target_label"],
                value_to_csv(
                    item["ca_distance_A"]
                ),
                value_to_csv(
                    item["min_distance_A"]
                ),
                pair_to_csv(
                    item["atom_pair"]
                ),
                value_to_csv(
                    item["centroid_distance_A"]
                ),
                item["within_8A"],
                item["within_centroid_14A"],
                item["backbone_contact"],
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
            "ca_distance_A",
            "min_distance_A",
            "atom_pair",
            "centroid_distance_A",
            "within_8A",
            "within_centroid_14A",
            "backbone_contact",
            "d_ET_proxy_A",
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
                value_to_csv(
                    item["ca_distance_A"]
                ),
                value_to_csv(
                    item["min_distance_A"]
                ),
                pair_to_csv(
                    item["atom_pair"]
                ),
                value_to_csv(
                    item["centroid_distance_A"]
                ),
                item["within_8A"],
                item["within_centroid_14A"],
                item["backbone_contact"],
                value_to_csv(
                    item["centroid_distance_A"]
                ),
                "",
                "",
                "",
                "",
                "",
                "",
                "geometry only; QM and kinetic data missing",
            ])

    return output


def print_selected(selected):
    print()
    print("Archi selezionati:")

    for item in selected:
        print(
            f"  {item['source']} -> {item['target']} "
            f"({item['source_chain']}-"
            f"{item['target_chain']})"
        )

        print(
            "    distanza minima non-backbone: "
            f"{value_to_csv(item['min_distance_A'])} Å"
        )

        print(
            "    distanza centroidi: "
            f"{value_to_csv(item['centroid_distance_A'])} Å"
        )

        print(
            "    distanza CA-CA: "
            f"{value_to_csv(item['ca_distance_A'])} Å"
        )

        print(
            "    coppia atomica: "
            f"{pair_to_csv(item['atom_pair'])}"
        )

        print(
            "    entro 8 Å: "
            f"{item['within_8A']}"
        )

        print(
            "    entro 14 Å tra centroidi: "
            f"{item['within_centroid_14A']}"
        )


def main():
    input_file = (
        Path(sys.argv[1])
        if len(sys.argv) > 1
        else Path(INPUT_DEFAULT)
    )

    if not input_file.exists():
        raise FileNotFoundError(
            f"File non trovato: {input_file}\n"
            "Metti 6w4x.cif nella stessa cartella."
        )

    OUTPUT_DIR.mkdir(exist_ok=True)

    parser = MMCIFParser(QUIET=True)

    structure = parser.get_structure(
        "6W4X",
        str(input_file),
    )

    found = find_pathway_residues(structure)
    candidates = build_candidates(found)
    selected = select_best_edges(candidates)

    all_candidates_file = write_all_candidates(
        candidates
    )

    selected_file = write_selected_edges(
        selected
    )

    detected_nodes = sorted({
        item["node"]
        for item in found
    })

    expected_nodes = [
        item[0]
        for item in PATHWAY
    ]

    missing_nodes = [
        node
        for node in expected_nodes
        if node not in detected_nodes
    ]

    print(f"Struttura analizzata: {input_file}")
    print(
        f"Residui del percorso trovati: "
        f"{len(found)}"
    )

    print(
        "Nodi rilevati: "
        + (
            ", ".join(detected_nodes)
            if detected_nodes
            else "nessuno"
        )
    )

    if missing_nodes:
        print()
        print("ATTENZIONE: nodi non trovati:")
        print("  " + ", ".join(missing_nodes))

    print_selected(selected)

    print()
    print("File prodotti:")
    print(f"  {all_candidates_file}")
    print(f"  {selected_file}")


if __name__ == "__main__":
    main()

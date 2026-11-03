#!/usr/bin/env python3
import sys
import numpy as np
from Bio.PDB import MMCIFParser

CIF = sys.argv[1] if len(sys.argv) > 1 else "6W4X.cif"

parser = MMCIFParser(QUIET=True)
structure = parser.get_structure("6W4X", CIF)
model = next(structure.get_models())


def residue(chain_id, number, name):
    hits = [
        r for r in model[chain_id].get_residues()
        if r.id[0] == " "
        and r.id[1] == number
        and r.get_resname().strip() == name
    ]

    if len(hits) != 1:
        raise ValueError(
            f"Residuo non trovato univocamente: {chain_id}:{name}{number}; "
            f"hits={len(hits)}"
        )

    return hits[0]


def atom_xyz(atom):
    return atom.get_coord().astype(float)


def heavy(atom):
    return atom.element.strip().upper() not in {"H", "D"}


def sidechain_atoms(residue):
    backbone = {"N", "CA", "C", "O", "OXT"}
    return [
        atom for atom in residue.get_atoms()
        if atom.get_name().strip() not in backbone and heavy(atom)
    ]


def ring_atoms(residue):
    ring_names = {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"}
    return [
        atom for atom in residue.get_atoms()
        if atom.get_name().strip() in ring_names and heavy(atom)
    ]


def all_heavy_atoms(residue):
    return [atom for atom in residue.get_atoms() if heavy(atom)]


def min_pair(atoms_a, atoms_b):
    rows = []

    for atom_a in atoms_a:
        for atom_b in atoms_b:
            distance = float(np.linalg.norm(
                atom_xyz(atom_a) - atom_xyz(atom_b)
            ))

            rows.append({
                "atom1": atom_a.get_name().strip(),
                "altloc1": atom_a.get_altloc(),
                "occupancy1": atom_a.get_occupancy(),
                "atom2": atom_b.get_name().strip(),
                "altloc2": atom_b.get_altloc(),
                "occupancy2": atom_b.get_occupancy(),
                "distance_A": distance,
            })

    return sorted(rows, key=lambda row: row["distance_A"])


y731 = residue("B", 731, "TYR")
y730 = residue("B", 730, "TYR")

categories = {
    "ALL_HEAVY_ATOMS": (
        all_heavy_atoms(y731),
        all_heavy_atoms(y730),
    ),
    "SIDECHAIN_HEAVY_ATOMS": (
        sidechain_atoms(y731),
        sidechain_atoms(y730),
    ),
    "AROMATIC_RING_ATOMS": (
        ring_atoms(y731),
        ring_atoms(y730),
    ),
}

for label, (atoms_a, atoms_b) in categories.items():
    rows = min_pair(atoms_a, atoms_b)

    print(f"\n=== {label} ===")
    print(f"n_atoms_Y731 = {len(atoms_a)}")
    print(f"n_atoms_Y730 = {len(atoms_b)}")

    for row in rows[:10]:
        print(row)

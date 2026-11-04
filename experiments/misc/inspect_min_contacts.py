#!/usr/bin/env python3
import sys
import numpy as np
from Bio.PDB import MMCIFParser

cif = sys.argv[1] if len(sys.argv) > 1 else "6W4X.cif"
parser = MMCIFParser(QUIET=True)
structure = parser.get_structure("6W4X", cif)
model = next(structure.get_models())

def residue(chain, number, name):
    hits = [
        r for r in model[chain].get_residues()
        if r.id[0] == " "
        and r.id[1] == number
        and r.get_resname().strip() == name
    ]
    if len(hits) != 1:
        raise ValueError((chain, number, name, len(hits)))
    return hits[0]

def heavy_atoms(r):
    return [
        a for a in r.get_atoms()
        if a.element.strip().upper() not in {"H", "D"}
    ]

def contacts(r1, r2, cutoff=3.0):
    rows = []
    for a in heavy_atoms(r1):
        for b in heavy_atoms(r2):
            d = float(np.linalg.norm(a.get_coord() - b.get_coord()))
            if d <= cutoff:
                rows.append({
                    "atom1": a.get_full_id(),
                    "atom1_name": a.get_name(),
                    "atom1_altloc": a.get_altloc(),
                    "atom1_occ": a.get_occupancy(),
                    "atom2": b.get_full_id(),
                    "atom2_name": b.get_name(),
                    "atom2_altloc": b.get_altloc(),
                    "atom2_occ": b.get_occupancy(),
                    "distance_A": d,
                })
    return sorted(rows, key=lambda x: x["distance_A"])

pairs = [
    ("Y731_B-Y730_B", residue("B", 731, "TYR"), residue("B", 730, "TYR")),
    ("Y730_B-C439_B", residue("B", 730, "TYR"), residue("B", 439, "CYS")),
    ("Y356_C-Y731_B", residue("C", 356, "TYR"), residue("B", 731, "TYR")),
]

for label, r1, r2 in pairs:
    print(f"\n{label}")
    for row in contacts(r1, r2):
        print(row)

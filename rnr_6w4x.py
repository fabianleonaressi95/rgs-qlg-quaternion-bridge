#!/usr/bin/env python3
import sys
from collections import defaultdict
from Bio.PDB import MMCIFParser

CIF = sys.argv[1] if len(sys.argv) > 1 else "6W4X.cif"

targets = [
    ("Y122", "TYR", 122),
    ("W48",  "TRP", 48),
    ("Y356", "TYR", 356),
    ("Y731", "TYR", 731),
    ("Y730", "TYR", 730),
    ("C439", "CYS", 439),
]

parser = MMCIFParser(QUIET=True)
structure = parser.get_structure("6W4X", CIF)
model = next(structure.get_models())

print(f"\nFILE: {CIF}")
print(f"MODEL: {model.id}\n")

all_hits = defaultdict(list)

for chain in model:
    residues = [
        r for r in chain.get_residues()
        if r.id[0] == " "
    ]

    if not residues:
        print(f"CHAIN {chain.id}: no polymer residues")
        continue

    first = f"{residues[0].get_resname()}{residues[0].id[1]}"
    last = f"{residues[-1].get_resname()}{residues[-1].id[1]}"

    print(
        f"CHAIN {chain.id}: "
        f"n_polymer={len(residues)}, "
        f"range={first}..{last}"
    )

    for label, expected_resname, expected_number in targets:
        matches = [
            r for r in residues
            if r.id[1] == expected_number
        ]

        for residue in matches:
            found = residue.get_resname().strip()
            all_hits[label].append(
                (chain.id, found, residue.id[1], residue.id[2])
            )
            print(
                f"  TARGET_NUMBER {label}: "
                f"found {found}{residue.id[1]}"
                f"{str(residue.id[2]).strip()}"
            )

print("\n=== RISULTATI PER LE SEI POSIZIONI ===")
for label, expected_resname, expected_number in targets:
    hits = all_hits[label]
    if not hits:
        print(
            f"{label}: NOT FOUND as coordinate residue "
            f"with numeric ID {expected_number}"
        )
    else:
        formatted = ", ".join(
            f"chain={c}, residue={r}{n}{str(i).strip()}"
            for c, r, n, i in hits
        )
        print(f"{label}: {formatted}")

print("\n=== TYR / TRP / CYS NELLE REGIONI DI INTERESSE ===")
for chain in model:
    residues = [
        r for r in chain.get_residues()
        if r.id[0] == " "
    ]

    near = []
    for residue in residues:
        name = residue.get_resname().strip()
        num = residue.id[1]

        if (
            (name == "TYR" and num in range(110, 135))
            or (name == "TRP" and num in range(35, 60))
            or (name == "TYR" and num in range(340, 370))
            or (name == "TYR" and num in range(715, 745))
            or (name == "CYS" and num in range(425, 455))
        ):
            near.append(f"{name}{num}{str(residue.id[2]).strip()}")

    if near:
        print(f"CHAIN {chain.id}: {', '.join(near)}")

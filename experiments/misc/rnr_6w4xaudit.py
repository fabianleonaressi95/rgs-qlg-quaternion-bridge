#!/usr/bin/env python3
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser


TYR_RING = ["CG", "CD1", "CE1", "CZ", "CE2", "CD2"]
TRP_INDOLE = ["CG", "CD1", "NE1", "CE2", "CD2", "CE3", "CZ3", "CH2", "CZ2"]

BETA_CHAIN = "C"
ALPHA_CANDIDATES = ["A", "B"]

NODES = {
    "W48":  {"chain": "C", "resid": 48,  "resname": "TRP", "subunit": "beta"},
    "Y356": {"chain": "C", "resid": 356, "resname": "TYR", "subunit": "beta"},
}

ALPHA_NODES = {
    "Y731": {"resid": 731, "resname": "TYR", "subunit": "alpha"},
    "Y730": {"resid": 730, "resname": "TYR", "subunit": "alpha"},
    "C439": {"resid": 439, "resname": "CYS", "subunit": "alpha"},
}


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def preferred_atom(residue, atom_name):
    atoms = [
        atom for atom in residue.get_atoms()
        if atom.get_name().strip() == atom_name
    ]

    if not atoms:
        raise KeyError(
            f"Atom {atom_name} assente in {residue.get_resname()} {residue.id}"
        )

    def sort_key(atom):
        occupancy = atom.get_occupancy()
        occupancy = -1.0 if occupancy is None else float(occupancy)
        altloc = atom.get_altloc()
        altloc_preference = 1 if altloc in (" ", "A", "") else 0
        return occupancy, altloc_preference

    atoms.sort(key=sort_key, reverse=True)
    return atoms[0]


def xyz(residue, atom_name):
    return preferred_atom(residue, atom_name).get_coord().astype(float)


def atom_set_xyz(residue, atom_names):
    return np.vstack([xyz(residue, atom_name) for atom_name in atom_names])


def unit(v):
    norm = np.linalg.norm(v)
    if norm == 0:
        raise ValueError("Vettore nullo")
    return v / norm


def centroid(points):
    return points.mean(axis=0)


def plane_normal(points):
    centered = points - centroid(points)
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    return unit(vt[-1])


def planarity_rmsd(points, normal):
    centered = points - centroid(points)
    return float(np.sqrt(np.mean((centered @ normal) ** 2)))


def epsilon_hex(points):
    c = centroid(points)
    n = plane_normal(points)

    u = points[0] - c
    u = unit(u - n * np.dot(u, n))
    v = np.cross(n, u)

    phi = np.arctan2((points - c) @ v, (points - c) @ u)
    delta = np.array([
        (phi[(i + 1) % 6] - phi[i]) % (2 * np.pi)
        for i in range(6)
    ])

    epsilon = float(np.sqrt(np.mean((delta - np.pi / 3) ** 2)))
    return epsilon, delta


def heavy_coordinates(residue):
    coords = []
    for atom in residue.get_atoms():
        element = atom.element.strip().upper()
        if element not in {"H", "D"}:
            coords.append(atom.get_coord().astype(float))
    return np.vstack(coords)


def min_heavy_distance(residue_a, residue_b):
    a = heavy_coordinates(residue_a)
    b = heavy_coordinates(residue_b)
    distances = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
    return float(np.min(distances))


def get_residue(model, chain_id, resid, expected_name):
    if chain_id not in model:
        raise KeyError(f"Chain {chain_id} non presente nel modello.")

    chain = model[chain_id]
    matches = [
        residue for residue in chain.get_residues()
        if residue.id[0] == " "
        and residue.id[1] == resid
        and residue.get_resname().strip() == expected_name
    ]

    if len(matches) != 1:
        found = [
            f"{residue.get_resname()}{residue.id[1]}"
            for residue in chain.get_residues()
            if residue.id[0] == " " and residue.id[1] == resid
        ]
        raise ValueError(
            f"Atteso {expected_name}{resid} in chain {chain_id}; trovato: {found}"
        )

    return matches[0]


def aromatic_descriptor(residue):
    name = residue.get_resname().strip()

    if name == "TYR":
        points = atom_set_xyz(residue, TYR_RING)
        n = plane_normal(points)
        eps, increments = epsilon_hex(points)
        edges = np.linalg.norm(np.roll(points, -1, axis=0) - points, axis=1)

        return {
            "type": "TYR_hexagon",
            "centroid": centroid(points),
            "normal": n,
            "epsilon_hex_rad": eps,
            "azimuth_increments_rad": ";".join(f"{x:.8f}" for x in increments),
            "planarity_rmsd_A": planarity_rmsd(points, n),
            "mean_ring_edge_A": float(np.mean(edges)),
            "ring_edge_sd_A": float(np.std(edges)),
        }

    if name == "TRP":
        points = atom_set_xyz(residue, TRP_INDOLE)
        n = plane_normal(points)

        return {
            "type": "TRP_indole",
            "centroid": centroid(points),
            "normal": n,
            "epsilon_hex_rad": np.nan,
            "azimuth_increments_rad": np.nan,
            "planarity_rmsd_A": planarity_rmsd(points, n),
            "mean_ring_edge_A": np.nan,
            "ring_edge_sd_A": np.nan,
        }

    if name == "CYS":
        return {
            "type": "CYS_non_aromatic",
            "centroid": None,
            "normal": None,
            "epsilon_hex_rad": np.nan,
            "azimuth_increments_rad": np.nan,
            "planarity_rmsd_A": np.nan,
            "mean_ring_edge_A": np.nan,
            "ring_edge_sd_A": np.nan,
        }

    raise ValueError(f"Residuo non previsto: {name}")


def ring_angle(desc_a, desc_b):
    if desc_a["normal"] is None or desc_b["normal"] is None:
        return np.nan
    dot = abs(float(np.dot(desc_a["normal"], desc_b["normal"])))
    return float(np.degrees(np.arccos(np.clip(dot, -1.0, 1.0))))


def centroid_distance(desc_a, desc_b):
    if desc_a["centroid"] is None or desc_b["centroid"] is None:
        return np.nan
    return float(np.linalg.norm(desc_a["centroid"] - desc_b["centroid"]))


def select_alpha_partner(model, beta_y356):
    beta_desc = aromatic_descriptor(beta_y356)
    results = []

    for chain_id in ALPHA_CANDIDATES:
        y731 = get_residue(model, chain_id, 731, "TYR")
        y731_desc = aromatic_descriptor(y731)
        d = centroid_distance(beta_desc, y731_desc)
        results.append((chain_id, d))

    results.sort(key=lambda item: item[1])
    return results[0][0], results


def main():
    argp = argparse.ArgumentParser()
    argp.add_argument("--cif", default="6W4X.cif")
    argp.add_argument("--outdir", default="rnr_6W4X_actual_audit")
    args = argp.parse_args()

    cif = Path(args.cif)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not cif.exists():
        raise FileNotFoundError(f"File non trovato: {cif}")

    structure = MMCIFParser(QUIET=True).get_structure("6W4X", str(cif))
    model = next(structure.get_models())

    w48 = get_residue(model, "C", 48, "TRP")
    y356 = get_residue(model, "C", 356, "TYR")

    alpha_chain, alpha_candidates = select_alpha_partner(model, y356)
    y731 = get_residue(model, alpha_chain, 731, "TYR")
    y730 = get_residue(model, alpha_chain, 730, "TYR")
    c439 = get_residue(model, alpha_chain, 439, "CYS")

    residues = {
        "W48": w48,
        "Y356": y356,
        "Y731": y731,
        "Y730": y730,
        "C439": c439,
    }

    descriptors = {
        key: aromatic_descriptor(residue)
        for key, residue in residues.items()
    }

    edges = [
        ("e2", "W48", "Y356"),
        ("e3", "Y356", "Y731"),
        ("e4", "Y731", "Y730"),
        ("e5", "Y730", "C439"),
    ]

    rows = []

    rows.append({
        "edge": "e1",
        "source": "Y122",
        "target": "W48",
        "source_chain": "UNAVAILABLE",
        "target_chain": "C",
        "PDB_ID": "6W4X",
        "geometry_source": "6W4X_model_0",
        "d_centroid_A": np.nan,
        "d_CA_A": np.nan,
        "d_min_heavy_A": np.nan,
        "ring_normal_angle_deg": np.nan,
        "source_epsilon_hex_rad": np.nan,
        "target_epsilon_hex_rad": descriptors["W48"]["epsilon_hex_rad"],
        "d_Tyr_centroid_to_Cys_SG_A": np.nan,
        "angle_Tyr_normal_to_Cys_SG_deg": np.nan,
        "geometry_status": "UNAVAILABLE",
        "provenance": "6W4X",
        "notes": (
            "Y122_beta is not modelled in the deposited beta chain C. "
            "Do not substitute alpha-chain Asp122."
        ),
    })

    for edge, source, target in edges:
        source_res = residues[source]
        target_res = residues[target]
        source_desc = descriptors[source]
        target_desc = descriptors[target]

        d_tyr_sg = np.nan
        tyr_sg_angle = np.nan
        notes = []

        if edge == "e5":
            sg = xyz(target_res, "SG")
            delta = sg - source_desc["centroid"]
            d_tyr_sg = float(np.linalg.norm(delta))
            direction = unit(delta)
            dot = abs(float(np.dot(source_desc["normal"], direction)))
            tyr_sg_angle = float(
                np.degrees(np.arccos(np.clip(dot, -1.0, 1.0)))
            )
            notes.append(
                "C439 is non-aromatic; ring_normal_angle_deg is NA by design."
            )

        if edge == "e3":
            notes.append(
                "Selected alpha partner by minimum Y356_C-to-Y731_alpha centroid distance."
            )

        rows.append({
            "edge": edge,
            "source": source,
            "target": target,
            "source_chain": source_res.get_parent().id,
            "target_chain": target_res.get_parent().id,
            "PDB_ID": "6W4X",
            "geometry_source": f"6W4X_model_0_chains_{source_res.get_parent().id}_{target_res.get_parent().id}",
            "d_centroid_A": centroid_distance(source_desc, target_desc),
            "d_CA_A": float(np.linalg.norm(
                xyz(source_res, "CA") - xyz(target_res, "CA")
            )),
            "d_min_heavy_A": min_heavy_distance(source_res, target_res),
            "ring_normal_angle_deg": (
                np.nan if edge == "e5" else ring_angle(source_desc, target_desc)
            ),
            "source_epsilon_hex_rad": source_desc["epsilon_hex_rad"],
            "target_epsilon_hex_rad": target_desc["epsilon_hex_rad"],
            "d_Tyr_centroid_to_Cys_SG_A": d_tyr_sg,
            "angle_Tyr_normal_to_Cys_SG_deg": tyr_sg_angle,
            "geometry_status": "GEOMETRY_DERIVED",
            "provenance": f"6W4X.cif:sha256={sha256_file(cif)}",
            "notes": " ".join(notes),
        })

    df = pd.DataFrame(rows)
    edge_csv = outdir / "rnr_6W4X_observed_geometry.csv"
    df.to_csv(edge_csv, index=False)

    manifest = {
        "PDB_ID": "6W4X",
        "model_id": 0,
        "beta_chain_used": "C",
        "alpha_chain_selected": alpha_chain,
        "alpha_candidate_Y356_to_Y731_centroid_distances_A": {
            c: d for c, d in alpha_candidates
        },
        "unavailable_node": {
            "node": "Y122",
            "reason": (
                "Y122_beta not modelled in beta chain C of deposited coordinate model; "
                "6W4X contains a trapped F3Y122_beta system, but that coordinate is absent "
                "from the inspected coordinate segment."
            ),
        },
        "input_sha256": sha256_file(cif),
        "output_csv": str(edge_csv),
        "output_csv_sha256": sha256_file(edge_csv),
    }

    manifest_file = outdir / "rnr_6W4X_observed_geometry_manifest.json"
    manifest_file.write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    print("\n=== 6W4X: SELEZIONE FISICA ===")
    print(f"beta chain: C")
    print(f"alpha candidates (Y356_C -> Y731_alpha centroid distance, Å):")
    for chain_id, distance in alpha_candidates:
        print(f"  {chain_id}: {distance:.4f} Å")
    print(f"selected alpha chain: {alpha_chain}")

    print("\n=== RETE BIOLOGICA ===")
    print("Y122_beta -> W48_beta -> Y356_beta -> Y731_alpha -> Y730_alpha -> C439_alpha")

    print("\n=== COPERTURA GEOMETRICA 6W4X ===")
    print("e1 Y122->W48: UNAVAILABLE (Y122_beta non coordinata nella beta-C selezionata)")
    print("e2-e5: GEOMETRY_DERIVED")

    print("\n=== ARCHI ===")
    print(df.to_string(index=False))

    print("\n=== FILE ===")
    print(edge_csv)
    print(manifest_file)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERRORE AUDIT: {type(error).__name__}: {error}", file=sys.stderr)
        sys.exit(1)

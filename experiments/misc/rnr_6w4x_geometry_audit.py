#!/usr/bin/env python3
import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser


TYR_RING_ORDER = ["CG", "CD1", "CE1", "CZ", "CE2", "CD2"]
TRP_INDOLE_ATOMS = ["CG", "CD1", "NE1", "CE2", "CD2", "CE3", "CZ3", "CH2", "CZ2"]

DEFAULT_NODES = {
    "Y122": {"chain": "C", "resid": 122, "resname": "TYR", "subunit": "beta"},
    "W48": {"chain": "C", "resid": 48, "resname": "TRP", "subunit": "beta"},
    "Y356": {"chain": "C", "resid": 356, "resname": "TYR", "subunit": "beta"},
    "Y731": {"chain": "B", "resid": 731, "resname": "TYR", "subunit": "alpha"},
    "Y730": {"chain": "B", "resid": 730, "resname": "TYR", "subunit": "alpha"},
    "C439": {"chain": "B", "resid": 439, "resname": "CYS", "subunit": "alpha"},
}

EDGES = [
    ("e1", "Y122", "W48"),
    ("e2", "W48", "Y356"),
    ("e3", "Y356", "Y731"),
    ("e4", "Y731", "Y730"),
    ("e5", "Y730", "C439"),
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def finite_or_na(value):
    if value is None:
        return np.nan
    try:
        value = float(value)
    except (TypeError, ValueError):
        return np.nan
    return value if np.isfinite(value) else np.nan


def atom_candidates(residue, atom_name):
    atoms = []
    for atom in residue.get_atoms():
        if atom.get_name().strip() == atom_name:
            atoms.append(atom)
    return atoms


def preferred_atom(residue, atom_name):
    atoms = atom_candidates(residue, atom_name)
    if not atoms:
        raise KeyError(
            f"Atom {atom_name} non trovato in "
            f"{residue.get_resname()} {residue.id}"
        )

    def key(atom):
        occ = finite_or_na(atom.get_occupancy())
        occ = -1.0 if np.isnan(occ) else occ
        alt = atom.get_altloc()
        alt_rank = 1 if alt in (" ", "A", "") else 0
        return (occ, alt_rank)

    atoms.sort(key=key, reverse=True)
    return atoms[0]


def atom_coord(residue, atom_name):
    return preferred_atom(residue, atom_name).get_coord().astype(float)


def ordered_coordinates(residue, atom_names):
    return np.vstack([atom_coord(residue, atom_name) for atom_name in atom_names])


def centroid(xyz):
    return xyz.mean(axis=0)


def unit(vector):
    n = np.linalg.norm(vector)
    if n == 0:
        raise ValueError("Vettore nullo: impossibile normalizzare.")
    return vector / n


def plane_normal(xyz):
    centered = xyz - centroid(xyz)
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    return unit(vt[-1])


def planarity_rmsd(xyz, normal):
    centered = xyz - centroid(xyz)
    perpendicular = centered @ normal
    return float(np.sqrt(np.mean(perpendicular ** 2)))


def ring_edges(xyz):
    return np.linalg.norm(np.roll(xyz, -1, axis=0) - xyz, axis=1)


def oriented_ring_azimuths(xyz, normal):
    center = centroid(xyz)
    u = xyz[0] - center
    u = u - normal * np.dot(u, normal)
    u = unit(u)
    v = np.cross(normal, u)
    rel = xyz - center
    return np.arctan2(rel @ v, rel @ u)


def wrap_positive(angle):
    return angle % (2.0 * np.pi)


def epsilon_hexagon(xyz):
    normal = plane_normal(xyz)
    phi = oriented_ring_azimuths(xyz, normal)
    dphi = np.array(
        [wrap_positive(phi[(k + 1) % 6] - phi[k]) for k in range(6)]
    )
    return float(np.sqrt(np.mean((dphi - np.pi / 3.0) ** 2))), dphi


def heavy_atom_coordinates(residue):
    xyz = []
    for atom in residue.get_atoms():
        element = atom.element.strip().upper()
        if element not in {"H", "D"}:
            xyz.append(atom.get_coord().astype(float))
    if not xyz:
        raise ValueError(f"Nessun heavy atom per {residue.get_resname()} {residue.id}")
    return np.vstack(xyz)


def min_heavy_distance(residue_a, residue_b):
    a = heavy_atom_coordinates(residue_a)
    b = heavy_atom_coordinates(residue_b)
    differences = a[:, None, :] - b[None, :, :]
    distances = np.linalg.norm(differences, axis=2)
    return float(np.min(distances))


def residue_summary(residue):
    atoms = list(residue.get_atoms())
    altlocs = sorted(set(str(a.get_altloc()).strip() for a in atoms if str(a.get_altloc()).strip()))
    occupancies = [finite_or_na(a.get_occupancy()) for a in atoms]
    b_factors = [finite_or_na(a.get_bfactor()) for a in atoms]
    occupancies = [x for x in occupancies if not np.isnan(x)]
    b_factors = [x for x in b_factors if not np.isnan(x)]

    return {
        "resname": residue.get_resname(),
        "hetero_flag": residue.id[0],
        "auth_seq_id": int(residue.id[1]),
        "insertion_code": str(residue.id[2]).strip() or None,
        "altlocs_present": altlocs,
        "mean_occupancy": float(np.mean(occupancies)) if occupancies else np.nan,
        "mean_B_factor": float(np.mean(b_factors)) if b_factors else np.nan,
    }


def find_residue(chain, auth_resid, expected_resname):
    hits = [
        residue
        for residue in chain.get_residues()
        if residue.id[0] == " "
        and residue.id[1] == auth_resid
        and residue.get_resname().strip() == expected_resname
    ]

    if len(hits) != 1:
        seen = [
            f"{r.get_resname()}:{r.id[1]}{str(r.id[2]).strip()}"
            for r in chain.get_residues()
            if r.id[0] == " " and r.id[1] == auth_resid
        ]
        raise ValueError(
            f"Atteso esattamente un residuo {expected_resname} {auth_resid} "
            f"nella chain {chain.id}; trovati {len(hits)}. "
            f"Residui con lo stesso numero: {seen}"
        )
    return hits[0]


def aromatic_descriptor(residue):
    resname = residue.get_resname().strip()

    if resname == "TYR":
        xyz = ordered_coordinates(residue, TYR_RING_ORDER)
        normal = plane_normal(xyz)
        eps_hex, dphi = epsilon_hexagon(xyz)
        edges = ring_edges(xyz)

        return {
            "ring_type": "tyr_phenyl_hexagon",
            "ring_atoms": ";".join(TYR_RING_ORDER),
            "centroid": centroid(xyz),
            "normal": normal,
            "epsilon_hex_rad": eps_hex,
            "azimuth_increments_rad": ";".join(f"{x:.8f}" for x in dphi),
            "planarity_rmsd_A": planarity_rmsd(xyz, normal),
            "mean_ring_edge_A": float(np.mean(edges)),
            "ring_edge_sd_A": float(np.std(edges, ddof=0)),
        }

    if resname == "TRP":
        xyz = ordered_coordinates(residue, TRP_INDOLE_ATOMS)
        normal = plane_normal(xyz)

        return {
            "ring_type": "trp_indole_fused",
            "ring_atoms": ";".join(TRP_INDOLE_ATOMS),
            "centroid": centroid(xyz),
            "normal": normal,
            "epsilon_hex_rad": np.nan,
            "azimuth_increments_rad": np.nan,
            "planarity_rmsd_A": planarity_rmsd(xyz, normal),
            "mean_ring_edge_A": np.nan,
            "ring_edge_sd_A": np.nan,
        }

    if resname == "CYS":
        return {
            "ring_type": "none_cysteine",
            "ring_atoms": "",
            "centroid": np.nan,
            "normal": np.nan,
            "epsilon_hex_rad": np.nan,
            "azimuth_increments_rad": np.nan,
            "planarity_rmsd_A": np.nan,
            "mean_ring_edge_A": np.nan,
            "ring_edge_sd_A": np.nan,
        }

    raise ValueError(f"Residuo non supportato: {resname}")


def is_vector(value):
    return isinstance(value, np.ndarray) and value.shape == (3,)


def centroid_distance(desc_a, desc_b):
    if is_vector(desc_a["centroid"]) and is_vector(desc_b["centroid"]):
        return float(np.linalg.norm(desc_a["centroid"] - desc_b["centroid"]))
    return np.nan


def ring_normal_angle(desc_a, desc_b):
    if is_vector(desc_a["normal"]) and is_vector(desc_b["normal"]):
        cosine = abs(float(np.dot(desc_a["normal"], desc_b["normal"])))
        cosine = np.clip(cosine, -1.0, 1.0)
        return float(np.degrees(np.arccos(cosine)))
    return np.nan


def tyr_centroid_to_cys_sg(tyr_desc, cys_residue):
    if not is_vector(tyr_desc["centroid"]) or not is_vector(tyr_desc["normal"]):
        return np.nan, np.nan

    sg = atom_coord(cys_residue, "SG")
    delta = sg - tyr_desc["centroid"]
    distance = float(np.linalg.norm(delta))
    direction = unit(delta)
    cosine = abs(float(np.dot(tyr_desc["normal"], direction)))
    cosine = np.clip(cosine, -1.0, 1.0)
    angle = float(np.degrees(np.arccos(cosine)))
    return distance, angle


def find_structural_waters(model, residue, cutoff_A=4.0):
    anchor_atoms = []
    for atom_name in ("OH", "NE1", "SG"):
        atoms = atom_candidates(residue, atom_name)
        if atoms:
            anchor_atoms.append(preferred_atom(residue, atom_name).get_coord().astype(float))

    if not anchor_atoms:
        return []

    waters = []
    for chain in model:
        for candidate in chain.get_residues():
            if candidate.get_resname().strip() not in {"HOH", "WAT", "DOD"}:
                continue
            oxygen_atoms = atom_candidates(candidate, "O")
            if not oxygen_atoms:
                continue
            oxygen = preferred_atom(candidate, "O").get_coord().astype(float)
            d = min(float(np.linalg.norm(oxygen - anchor)) for anchor in anchor_atoms)
            if d <= cutoff_A:
                waters.append({
                    "water_chain": candidate.get_parent().id,
                    "water_resid": int(candidate.id[1]),
                    "distance_to_anchor_A": d,
                })

    waters.sort(key=lambda item: item["distance_to_anchor_A"])
    return waters


def main():
    parser = argparse.ArgumentParser(
        description="Audit geometrico atomico 6W4X per la rete RNR a sei nodi."
    )
    parser.add_argument(
        "--cif",
        default="6W4X.cif",
        help="Percorso al file mmCIF 6W4X. Default: 6W4X.cif",
    )
    parser.add_argument(
        "--outdir",
        default="rnr_6W4X_geometry_audit",
        help="Directory di output.",
    )
    parser.add_argument(
        "--model",
        type=int,
        default=1,
        help="Numero modello PDB, 1-based. Default: 1.",
    )
    parser.add_argument(
        "--beta-chain",
        default="C",
        help="Chain candidata beta. Default: C.",
    )
    parser.add_argument(
        "--alpha-chain",
        default="B",
        help="Chain candidata alpha. Default: B.",
    )
    parser.add_argument(
        "--water-cutoff",
        type=float,
        default=4.0,
        help="Cutoff per acqua strutturale in Å. Default: 4.0.",
    )
    args = parser.parse_args()

    cif_path = Path(args.cif).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    if not cif_path.exists():
        raise FileNotFoundError(f"File mmCIF non trovato: {cif_path}")

    parser_mmcif = MMCIFParser(QUIET=True)
    structure = parser_mmcif.get_structure("6W4X", str(cif_path))

    models = list(structure.get_models())
    if args.model < 1 or args.model > len(models):
        raise ValueError(
            f"Model {args.model} non disponibile. Modelli trovati: {len(models)}"
        )
    model = models[args.model - 1]

    nodes = {key: dict(value) for key, value in DEFAULT_NODES.items()}
    for key in ("Y122", "W48", "Y356"):
        nodes[key]["chain"] = args.beta_chain
    for key in ("Y731", "Y730", "C439"):
        nodes[key]["chain"] = args.alpha_chain

    available_chains = [chain.id for chain in model]
    missing_chains = sorted(set(v["chain"] for v in nodes.values()) - set(available_chains))
    if missing_chains:
        raise KeyError(
            f"Chain mancanti nel modello {args.model}: {missing_chains}. "
            f"Chain disponibili: {available_chains}"
        )

    residues = {}
    node_rows = []

    for label, spec in nodes.items():
        chain = model[spec["chain"]]
        residue = find_residue(chain, spec["resid"], spec["resname"])
        residues[label] = residue
        desc = aromatic_descriptor(residue)
        summary = residue_summary(residue)
        waters = find_structural_waters(model, residue, cutoff_A=args.water_cutoff)

        node_rows.append({
            "node": label,
            "subunit": spec["subunit"],
            "amino_acid": spec["resname"],
            "PDB_ID": "6W4X",
            "model_id": args.model,
            "auth_asym_id": chain.id,
            "auth_seq_id": summary["auth_seq_id"],
            "insertion_code": summary["insertion_code"],
            "ring_type": desc["ring_type"],
            "ring_atoms": desc["ring_atoms"],
            "centroid_x_A": desc["centroid"][0] if is_vector(desc["centroid"]) else np.nan,
            "centroid_y_A": desc["centroid"][1] if is_vector(desc["centroid"]) else np.nan,
            "centroid_z_A": desc["centroid"][2] if is_vector(desc["centroid"]) else np.nan,
            "normal_x": desc["normal"][0] if is_vector(desc["normal"]) else np.nan,
            "normal_y": desc["normal"][1] if is_vector(desc["normal"]) else np.nan,
            "normal_z": desc["normal"][2] if is_vector(desc["normal"]) else np.nan,
            "epsilon_hex_rad": desc["epsilon_hex_rad"],
            "azimuth_increments_rad": desc["azimuth_increments_rad"],
            "planarity_rmsd_A": desc["planarity_rmsd_A"],
            "mean_ring_edge_A": desc["mean_ring_edge_A"],
            "ring_edge_sd_A": desc["ring_edge_sd_A"],
            "altlocs_present": ";".join(summary["altlocs_present"]),
            "mean_occupancy": summary["mean_occupancy"],
            "mean_B_factor": summary["mean_B_factor"],
            "waters_within_cutoff": json.dumps(waters, sort_keys=True),
            "geometry_status": "GEOMETRY_DERIVED",
            "provenance": f"6W4X.cif:sha256={sha256_file(cif_path)}",
        })

    descriptors = {
        label: aromatic_descriptor(residue)
        for label, residue in residues.items()
    }

    edge_rows = []
    for edge_id, source, target in EDGES:
        source_res = residues[source]
        target_res = residues[target]
        source_desc = descriptors[source]
        target_desc = descriptors[target]

        d_tyr_to_sg = np.nan
        tyr_normal_to_sg_angle = np.nan
        notes = []

        if edge_id == "e5":
            d_tyr_to_sg, tyr_normal_to_sg_angle = tyr_centroid_to_cys_sg(
                source_desc, target_res
            )
            notes.append(
                "Cys non aromatico: ring_normal_angle_deg intenzionalmente NA; "
                "usare d_Tyr_centroid_to_Cys_SG_A e angle_Tyr_normal_to_Cys_SG_deg."
            )

        if edge_id in {"e1", "e2"}:
            notes.append(
                "Trp indolico: ring_normal_angle_deg è un descrittore piano-piano, "
                "non un confronto tra due esagoni Tyr."
            )

        if edge_id == "e3":
            notes.append(
                "Arco interfaciale beta-alpha; geometria estratta esclusivamente da 6W4X."
            )

        if edge_id == "e4":
            notes.append(
                "Coppia Tyr-Tyr alpha; valutare stacking tramite d_centroid_A, "
                "ring_normal_angle_deg e d_min_heavy_A."
            )

        edge_rows.append({
            "edge": edge_id,
            "source": source,
            "target": target,
            "source_amino_acid": source_res.get_resname().strip(),
            "target_amino_acid": target_res.get_resname().strip(),
            "source_chain": source_res.get_parent().id,
            "target_chain": target_res.get_parent().id,
            "PDB_ID": "6W4X",
            "model_id": args.model,
            "geometry_source": (
                f"6W4X_model_{args.model}_"
                f"chains_{source_res.get_parent().id}_{target_res.get_parent().id}"
            ),
            "source_auth_seq_id": int(source_res.id[1]),
            "target_auth_seq_id": int(target_res.id[1]),
            "d_centroid_A": centroid_distance(source_desc, target_desc),
            "d_CA_A": float(np.linalg.norm(
                atom_coord(source_res, "CA") - atom_coord(target_res, "CA")
            )),
            "d_min_heavy_A": min_heavy_distance(source_res, target_res),
            "ring_normal_angle_deg": (
                np.nan if edge_id == "e5"
                else ring_normal_angle(source_desc, target_desc)
            ),
            "source_epsilon_hex_rad": source_desc["epsilon_hex_rad"],
            "target_epsilon_hex_rad": target_desc["epsilon_hex_rad"],
            "source_planarity_rmsd_A": source_desc["planarity_rmsd_A"],
            "target_planarity_rmsd_A": target_desc["planarity_rmsd_A"],
            "source_mean_ring_edge_A": source_desc["mean_ring_edge_A"],
            "target_mean_ring_edge_A": target_desc["mean_ring_edge_A"],
            "d_Tyr_centroid_to_Cys_SG_A": d_tyr_to_sg,
            "angle_Tyr_normal_to_Cys_SG_deg": tyr_normal_to_sg_angle,
            "d_H_A": np.nan,
            "Hbond_angle_deg": np.nan,
            "geometry_status": "GEOMETRY_DERIVED",
            "provenance": f"6W4X.cif:sha256={sha256_file(cif_path)}",
            "notes": " ".join(notes),
        })

    node_df = pd.DataFrame(node_rows)
    edge_df = pd.DataFrame(edge_rows)

    node_csv = outdir / "rnr_6W4X_nodes_geometry.csv"
    edge_csv = outdir / "rnr_6W4X_edges_geometry.csv"
    node_df.to_csv(node_csv, index=False)
    edge_df.to_csv(edge_csv, index=False)

    manifest = {
        "PDB_ID": "6W4X",
        "input_mmcif": str(cif_path),
        "input_mmcif_sha256": sha256_file(cif_path),
        "model_id": args.model,
        "candidate_beta_chain": args.beta_chain,
        "candidate_alpha_chain": args.alpha_chain,
        "water_cutoff_A": args.water_cutoff,
        "nodes": nodes,
        "edges": EDGES,
        "output_files": {
            "nodes_csv": str(node_csv),
            "edges_csv": str(edge_csv),
            "nodes_csv_sha256": sha256_file(node_csv),
            "edges_csv_sha256": sha256_file(edge_csv),
        },
        "interpretation_constraints": [
            "d_min_heavy_A is geometry-derived and is not an ET or PT distance.",
            "No H-bond angle is reported because hydrogen coordinates/protonation are not reconstructed.",
            "No kinetic, KIE, H_AB, lambda, or deltaG parameter is included in this geometry audit.",
            "Tyr epsilon_hex is an atomic-geometry descriptor and does not imply QVSM state multiplicity.",
        ],
    }

    manifest_json = outdir / "rnr_6W4X_geometry_manifest.json"
    manifest_json.write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    print("\n=== AUDIT COMPLETATO ===")
    print(f"Input:          {cif_path}")
    print(f"Input SHA-256:  {manifest['input_mmcif_sha256']}")
    print(f"Node CSV:       {node_csv}")
    print(f"Edge CSV:       {edge_csv}")
    print(f"Manifest:       {manifest_json}")
    print("\n=== NODI ===")
    print(node_df.to_string(index=False))
    print("\n=== ARCHI ===")
    print(edge_df.to_string(index=False))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"\nERRORE AUDIT: {type(exc).__name__}: {exc}", file=sys.stderr)
        sys.exit(1)

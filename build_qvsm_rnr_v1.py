#!/usr/bin/env python3

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


HBAR_J_S = 1.054571817e-34
AMU_KG = 1.66053906660e-27
EV_J = 1.602176634e-19

MASS_H_AMU = 1.00782503223
MASS_D_AMU = 2.01410177812

NETWORK_NODES = [
    "Y122",
    "W48",
    "Y356",
    "Y731",
    "Y730",
    "C439",
]

NETWORK_EDGES = [
    ("e1", "Y122", "W48"),
    ("e2", "W48", "Y356"),
    ("e3", "Y356", "Y731"),
    ("e4", "Y731", "Y730"),
    ("e5", "Y730", "C439"),
]

CODON_TO_G = {
    "TAT": 8,
    "TAC": 9,
    "TGG": 15,
    "TGT": 12,
    "TGC": 13,
}

NODE_CODONS = {
    "Y122": ["TAT", "TAC"],
    "W48": ["TGG"],
    "Y356": ["TAT", "TAC"],
    "Y731": ["TAT", "TAC"],
    "Y730": ["TAT", "TAC"],
    "C439": ["TGT", "TGC"],
}

AMINO_ACID = {
    "Y122": "Tyr",
    "W48": "Trp",
    "Y356": "Tyr",
    "Y731": "Tyr",
    "Y730": "Tyr",
    "C439": "Cys",
}

SUBUNIT = {
    "Y122": "beta",
    "W48": "beta",
    "Y356": "beta",
    "Y731": "alpha",
    "Y730": "alpha",
    "C439": "alpha",
}

# Coordinate QVSM pilota.
# L: localizzazione radicalica lungo il percorso.
L_INDEX = {
    "Y122": 0,
    "W48": np.nan,
    "Y356": 1,
    "Y731": 2,
    "Y730": 3,
    "C439": 4,
}

L_LABEL = {
    "Y122": "Y122 radicale stabile",
    "W48": "nessuna assegnazione radicalica discreta diretta",
    "Y356": "Y356 radicale/intermedio",
    "Y731": "Y731 radicale/intermedio",
    "Y730": "Y730 radicale/intermedio",
    "C439": "C439 stato tiile/catalitico",
}

L_STATUS = {
    "Y122": "MEASURED",
    "W48": "UNAVAILABLE",
    "Y356": "MEASURED_OR_INFERRED",
    "Y731": "MEASURED_OR_INFERRED",
    "Y730": "MEASURED_OR_INFERRED",
    "C439": "PROCESS_LEVEL",
}

# C: conformazione.
# Per il primo dataset usiamo una configurazione primaria conservativa.
# Le configurazioni alternative saranno mantenute in un file separato.
C_INDEX_PRIMARY = {
    "Y122": np.nan,
    "W48": np.nan,
    "Y356": 1,
    "Y731": 2,   # stacked, coerente con e4 6W4X
    "Y730": 2,   # stacked
    "C439": 4,   # substrate-competent, inferenziale
}

C_LABEL_PRIMARY = {
    "Y122": "non assegnata",
    "W48": "non assegnata",
    "Y356": "interfaciale/water-assisted",
    "Y731": "stacked",
    "Y730": "stacked",
    "C439": "substrate-competent",
}

C_STATUS = {
    "Y122": "UNAVAILABLE",
    "W48": "UNAVAILABLE",
    "Y356": "QM_MM_SUPPORTED",
    "Y731": "STRUCTURE_QM_MM_SUPPORTED",
    "Y730": "STRUCTURE_QM_MM_SUPPORTED",
    "C439": "INFERRED",
}

# T: regime PCET per arco.
# Questo è un attributo dell'arco, non una misura edge-rate.
T_INDEX_EDGE = {
    "e1": np.nan,
    "e2": np.nan,
    "e3": 1,  # nonadiabatic
    "e4": 3,  # adiabatic
    "e5": 2,  # intermediate
}

T_LABEL_EDGE = {
    "e1": "non assegnato",
    "e2": "non assegnato",
    "e3": "nonadiabatic",
    "e4": "adiabatic",
    "e5": "intermediate",
}

T_STATUS_EDGE = {
    "e1": "UNAVAILABLE",
    "e2": "UNAVAILABLE",
    "e3": "QM_DERIVED",
    "e4": "QM_DERIVED",
    "e5": "QM_DERIVED",
}

# B: rete acqua/H-bond.
# Sono classi pilota. Non rappresentano energia.
B_INDEX_PRIMARY = {
    "Y122": np.nan,
    "W48": np.nan,
    "Y356": 3,
    "Y731": 3,
    "Y730": 4,
    "C439": 4,
}

B_LABEL_PRIMARY = {
    "Y122": "non assegnata",
    "W48": "non assegnata",
    "Y356": "bridging-water-capable interface",
    "Y731": "interfacial water/H-bond network",
    "Y730": "multiple H-bond/water network",
    "C439": "multiple H-bond/water network",
}

B_STATUS = {
    "Y122": "UNAVAILABLE",
    "W48": "UNAVAILABLE",
    "Y356": "MD_SUPPORTED",
    "Y731": "ENDOR_MD_SUPPORTED",
    "Y730": "ENDOR_MD_SUPPORTED",
    "C439": "ENDOR_MD_SUPPORTED",
}


def sha256_file(path):
    digest = hashlib.sha256()

    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)

    return digest.hexdigest()


def qvsm_index(G, L, C, T, B):
    values = [G, L, C, T, B]

    if any(pd.isna(value) for value in values):
        return None

    G = int(G)
    L = int(L)
    C = int(C)
    T = int(T)
    B = int(B)

    if not 0 <= G < 64:
        raise ValueError(f"G fuori range: {G}")

    if not 0 <= L < 6:
        raise ValueError(f"L fuori range: {L}")

    if not 0 <= C < 6:
        raise ValueError(f"C fuori range: {C}")

    if not 0 <= T < 6:
        raise ValueError(f"T fuori range: {T}")

    if not 0 <= B < 5:
        raise ValueError(f"B fuori range: {B}")

    return ((((G * 6 + L) * 6 + C) * 6 + T) * 5 + B)


def parse_float(value):
    if value is None:
        return np.nan

    if isinstance(value, float) and np.isnan(value):
        return np.nan

    text = str(value).strip()

    if text in {"", "NA", "nan", "None"}:
        return np.nan

    try:
        return float(text)
    except ValueError:
        return np.nan


def value_from_dataframe(row, column):
    if column not in row.index:
        return np.nan

    return parse_float(row[column])


def get_edge_row(edge_df, edge_id):
    rows = edge_df.loc[edge_df["edge"].eq(edge_id)]

    if len(rows) != 1:
        raise ValueError(
            f"Attesa una sola riga per {edge_id}, trovate {len(rows)}"
        )

    return rows.iloc[0]


def get_first_present(row, columns):
    for column in columns:
        if column in row.index:
            value = parse_float(row[column])
            if not np.isnan(value):
                return value
    return np.nan


def edge_water_defaults(edge_id):
    """
    Valori descrittivi per la prima versione.
    Non sono numeri QM/MM e non entrano in kappa/NQE.
    """
    values = {
        "e1": {
            "water_bridge_status": "UNAVAILABLE",
            "water_network_status": "UNAVAILABLE",
            "B_edge": np.nan,
        },
        "e2": {
            "water_bridge_status": "UNAVAILABLE",
            "water_network_status": "UNAVAILABLE",
            "B_edge": np.nan,
        },
        "e3": {
            "water_bridge_status": "CONFIGURATION_DEPENDENT",
            "water_network_status": "MD_SUPPORTED",
            "B_edge": 3,
        },
        "e4": {
            "water_bridge_status": "NOT_REQUIRED_FOR_PRIMARY_STACKED_STATE",
            "water_network_status": "ENDOR_MD_SUPPORTED",
            "B_edge": 4,
        },
        "e5": {
            "water_bridge_status": "UNAVAILABLE",
            "water_network_status": "ENDOR_MD_SUPPORTED",
            "B_edge": 4,
        },
    }

    return values[edge_id]


def edge_kclass_defaults(edge_id):
    """
    Kinetic process constraints, non edge-resolved.
    """
    values = {
        "e1": {
            "k_edge_exp_s-1": np.nan,
            "k_process_exp_s-1": np.nan,
            "k_process_uncertainty_s-1": np.nan,
            "evidence_scope": "UNAVAILABLE",
        },
        "e2": {
            "k_edge_exp_s-1": np.nan,
            "k_process_exp_s-1": np.nan,
            "k_process_uncertainty_s-1": np.nan,
            "evidence_scope": "UNAVAILABLE",
        },
        "e3": {
            "k_edge_exp_s-1": np.nan,
            "k_process_exp_s-1": 1.0e4,
            "k_process_uncertainty_s-1": 0.5e4,
            "evidence_scope": "INTERFACE",
        },
        "e4": {
            "k_edge_exp_s-1": np.nan,
            "k_process_exp_s-1": 1.0e5,
            "k_process_uncertainty_s-1": np.nan,
            "evidence_scope": "PROCESS",
        },
        "e5": {
            "k_edge_exp_s-1": np.nan,
            "k_process_exp_s-1": 1.4e4,
            "k_process_uncertainty_s-1": 0.4e4,
            "evidence_scope": "PROCESS",
        },
    }

    return values[edge_id]


def wkb_exponent(a_eff_A, barrier_eV, mass_amu):
    """
    WKB exponent:
      S = 2 a / hbar * sqrt(2 m DeltaV)

    a_eff_A is interpreted as half-width/effective width in Å.
    barrier_eV is in eV.
    """
    if any(pd.isna(x) for x in [a_eff_A, barrier_eV, mass_amu]):
        return np.nan

    if a_eff_A <= 0 or barrier_eV <= 0 or mass_amu <= 0:
        return np.nan

    a_m = a_eff_A * 1.0e-10
    barrier_joule = barrier_eV * EV_J
    mass_kg = mass_amu * AMU_KG

    exponent = (
        2.0
        * a_m
        / HBAR_J_S
        * math.sqrt(2.0 * mass_kg * barrier_joule)
    )

    return float(exponent)


def wkb_rate(k0_s, a_eff_A, barrier_eV, mass_amu):
    if pd.isna(k0_s):
        return np.nan

    exponent = wkb_exponent(a_eff_A, barrier_eV, mass_amu)

    if pd.isna(exponent):
        return np.nan

    return float(k0_s * math.exp(-exponent))


def reduced_kie(a_eff_A, barrier_eV):
    """
    KIE = exp(S_D - S_H) for equal prefactors.
    """
    s_h = wkb_exponent(a_eff_A, barrier_eV, MASS_H_AMU)
    s_d = wkb_exponent(a_eff_A, barrier_eV, MASS_D_AMU)

    if pd.isna(s_h) or pd.isna(s_d):
        return np.nan

    value = math.exp(s_d - s_h)

    if not np.isfinite(value):
        return np.nan

    return float(value)


def build_edge_table(edge_df, geometry_path):
    rows = []

    for edge_id, source, target in NETWORK_EDGES:
        edge = get_edge_row(edge_df, edge_id)
        water = edge_water_defaults(edge_id)
        kinetic = edge_kclass_defaults(edge_id)

        d_centroid = get_first_present(
            edge,
            ["d_centroid_A"],
        )

        d_ring = get_first_present(
            edge,
            ["d_min_ring_A"],
        )

        d_sidechain = get_first_present(
            edge,
            ["d_min_sidechain_heavy_A"],
        )

        d_oh_sg = get_first_present(
            edge,
            ["d_Tyr_OH_to_Cys_SG_A"],
        )

        d_oo = d_oh_sg

        if edge_id == "e4":
            # Per Y731-Y730 non esiste un campo proton donor-acceptor
            # già calcolato nel PDB audit. Il valore resta NA.
            d_oo = np.nan

        # T è una classificazione di regime PCET.
        T = T_INDEX_EDGE[edge_id]

        # C-edge è una coppia rappresentativa primaria.
        if edge_id == "e3":
            C_source = 1
            C_target = 1
        elif edge_id == "e4":
            C_source = 2
            C_target = 2
        elif edge_id == "e5":
            C_source = 2
            C_target = 4
        else:
            C_source = np.nan
            C_target = np.nan

        row = {
            "edge": edge_id,
            "source": source,
            "target": target,
            "source_chain": edge.get("source_chain", "NA"),
            "target_chain": edge.get("target_chain", "NA"),
            "PDB_ID": edge.get("PDB_ID", "6W4X"),
            "geometry_source": edge.get("geometry_source", ""),
            "d_centroid_A": d_centroid,
            "d_OO_A": d_oo,
            "d_min_ring_A": d_ring,
            "d_min_sidechain_heavy_A": d_sidechain,
            "ring_normal_angle_deg": get_first_present(
                edge,
                ["ring_normal_angle_deg"],
            ),
            "source_epsilon_hex_rad": get_first_present(
                edge,
                ["source_epsilon_hex_rad"],
            ),
            "target_epsilon_hex_rad": get_first_present(
                edge,
                ["target_epsilon_hex_rad"],
            ),
            "C_source_index": C_source,
            "C_target_index": C_target,
            "C_source_label": (
                "flipped-out"
                if C_source == 1
                else "stacked"
                if C_source == 2
                else "substrate-competent"
                if C_source == 4
                else "NA"
            ),
            "C_target_label": (
                "flipped-out"
                if C_target == 1
                else "stacked"
                if C_target == 2
                else "substrate-competent"
                if C_target == 4
                else "NA"
            ),
            "C_status": (
                "QM_MM_SUPPORTED"
                if edge_id in {"e3", "e4"}
                else "INFERRED"
                if edge_id == "e5"
                else "UNAVAILABLE"
            ),
            "T_index": T,
            "T_label": T_LABEL_EDGE[edge_id],
            "T_status": T_STATUS_EDGE[edge_id],
            "B_edge_index": water["B_edge"],
            "B_edge_label": water["water_network_status"],
            "B_status": (
                "MD_SUPPORTED"
                if edge_id in {"e3", "e4", "e5"}
                else "UNAVAILABLE"
            ),
            "water_bridge_status": water["water_bridge_status"],
            "water_network_status": water["water_network_status"],
            "N_HB": np.nan,
            "N_bridge": np.nan,
            "water_occupancy": np.nan,
            "water_Hbond_lifetime_ps": np.nan,
            "chi_water_classical": np.nan,
            "chi_water_status": "UNAVAILABLE",
            "H_AB_QM_eV": np.nan,
            "lambda_QM_eV": np.nan,
            "deltaG_ET_exp_eV": np.nan,
            "DeltaV_barrier_eV": np.nan,
            "a_eff_A": np.nan,
            "k0_s-1": np.nan,
            "k_edge_exp_s-1": kinetic["k_edge_exp_s-1"],
            "k_process_exp_s-1": kinetic["k_process_exp_s-1"],
            "k_process_uncertainty_s-1": kinetic[
                "k_process_uncertainty_s-1"
            ],
            "evidence_scope": kinetic["evidence_scope"],
            "evidence_status": (
                "GEOMETRY_DERIVED_PLUS_PROCESS_CONSTRAINT"
                if edge_id in {"e3", "e4", "e5"}
                else "GEOMETRY_DERIVED_OR_UNAVAILABLE"
            ),
            "geometry_provenance": edge.get("provenance", ""),
            "source_reference": (
                "10.1021/jacs.0c04325;"
                "10.1021/jacs.2c13615;"
                "10.1021/jacs.4c13955"
            ),
            "notes": (
                "T is a PCET-regime class, not an edge rate. "
                "NQE/WKB fields remain NA until QM/MM surface samples "
                "are supplied."
            ),
        }

        rows.append(row)

    return pd.DataFrame(rows)


def build_node_table():
    rows = []

    for node in NETWORK_NODES:
        for codon in NODE_CODONS[node]:
            G = CODON_TO_G[codon]

            rows.append({
                "node": node,
                "amino_acid": AMINO_ACID[node],
                "subunit": SUBUNIT[node],
                "codon_DNA": codon,
                "G_index": G,
                "G_status": "CANONICAL_CODON",
                "L_index": L_INDEX[node],
                "L_label": L_LABEL[node],
                "L_status": L_STATUS[node],
                "C_index_primary": C_INDEX_PRIMARY[node],
                "C_label_primary": C_LABEL_PRIMARY[node],
                "C_status": C_STATUS[node],
                "T_index_node": np.nan,
                "T_label_node": "EDGE_LEVEL_ONLY",
                "T_status_node": "EDGE_LEVEL_ONLY",
                "B_index_primary": B_INDEX_PRIMARY[node],
                "B_label_primary": B_LABEL_PRIMARY[node],
                "B_status": B_STATUS[node],
                "PDB_ID": "6W4X",
                "geometry_status": (
                    "UNAVAILABLE"
                    if node == "Y122"
                    else "GEOMETRY_DERIVED"
                ),
                "notes": (
                    "Y122 beta coordinate absent in selected 6W4X beta-C."
                    if node == "Y122"
                    else
                    "Node state is provisional; edge geometry is stored separately."
                ),
            })

    return pd.DataFrame(rows)


def build_microstate_examples(node_df, edge_df):
    """
    Genera solo gli esempi completi realmente definiti.
    Non crea 69120 righe fittizie.
    """
    rows = []

    for _, node in node_df.iterrows():
        G = node["G_index"]
        L = node["L_index"]
        C = node["C_index_primary"]
        T = node["T_index_node"]
        B = node["B_index_primary"]

        index = qvsm_index(G, L, C, T, B)

        rows.append({
            "node": node["node"],
            "codon_DNA": node["codon_DNA"],
            "G": G,
            "L": L,
            "C": C,
            "T": T,
            "B": B,
            "qvsm_index": index,
            "status": (
                "COMPLETE"
                if index is not None
                else "PARTIAL_UNAVAILABLE"
            ),
        })

    return pd.DataFrame(rows)


def build_synthetic_wkb_demo():
    """
    Produce una demo numerica volutamente separata dai dati reali.
    Non viene mescolata al dataset sperimentale.
    """
    rows = []

    parameter_sets = [
        {
            "label": "DEMO_NOT_EXPERIMENTAL",
            "edge": "e4",
            "a_eff_A": 0.25,
            "DeltaV_barrier_eV": 0.10,
            "k0_s-1": 1.0e13,
        },
    ]

    for params in parameter_sets:
        a_eff = params["a_eff_A"]
        barrier = params["DeltaV_barrier_eV"]
        k0 = params["k0_s-1"]

        s_h = wkb_exponent(a_eff, barrier, MASS_H_AMU)
        s_d = wkb_exponent(a_eff, barrier, MASS_D_AMU)
        k_h = wkb_rate(k0, a_eff, barrier, MASS_H_AMU)
        k_d = wkb_rate(k0, a_eff, barrier, MASS_D_AMU)
        kie = reduced_kie(a_eff, barrier)

        rows.append({
            "label": params["label"],
            "edge": params["edge"],
            "a_eff_A": a_eff,
            "DeltaV_barrier_eV": barrier,
            "k0_s-1": k0,
            "S_H": s_h,
            "S_D": s_d,
            "k_H_s-1": k_h,
            "k_D_s-1": k_d,
            "KIE_H_over_D": kie,
            "status": "SYNTHETIC_DEMO_ONLY",
            "notes": (
                "Not an RNR measurement and not a QM/MM result. "
                "Replace with extracted NEO-QM/MM surface data."
            ),
        })

    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Build QVSM-RNR v1 from 6W4X geometry and explicit state metadata."
    )

    parser.add_argument(
        "--geometry-csv",
        default=(
            "rnr_6W4X_canonical_geometry/"
            "rnr_6W4X_canonical_edges.csv"
        ),
    )

    parser.add_argument(
        "--outdir",
        default="qvsm_rnr_v1_output",
    )

    parser.add_argument(
        "--include-wkb-demo",
        action="store_true",
        help="Include a clearly labeled synthetic WKB demonstration.",
    )

    args = parser.parse_args()

    geometry_path = Path(args.geometry_csv).resolve()
    output_directory = Path(args.outdir).resolve()
    output_directory.mkdir(parents=True, exist_ok=True)

    if not geometry_path.exists():
        raise FileNotFoundError(
            f"CSV geometrico non trovato: {geometry_path}"
        )

    geometry_sha256 = sha256_file(geometry_path)
    edge_df = pd.read_csv(geometry_path)

    expected_edges = {edge[0] for edge in NETWORK_EDGES}
    observed_edges = set(edge_df["edge"].astype(str))

    missing_edges = expected_edges - observed_edges

    if missing_edges:
        raise ValueError(
            f"Archi mancanti dal CSV geometrico: {sorted(missing_edges)}"
        )

    node_df = build_node_table()
    qvsm_edge_df = build_edge_table(edge_df, geometry_path)
    examples_df = build_microstate_examples(node_df, qvsm_edge_df)

    node_output = output_directory / "qvsm_rnr_v1_nodes.csv"
    edge_output = output_directory / "qvsm_rnr_v1_edges.csv"
    examples_output = output_directory / "qvsm_rnr_v1_microstate_examples.csv"

    node_df.to_csv(node_output, index=False)
    qvsm_edge_df.to_csv(edge_output, index=False)
    examples_df.to_csv(examples_output, index=False)

    outputs = {
        "nodes_csv": str(node_output),
        "edges_csv": str(edge_output),
        "microstate_examples_csv": str(examples_output),
    }

    if args.include_wkb_demo:
        demo_df = build_synthetic_wkb_demo()
        demo_output = output_directory / "qvsm_rnr_v1_wkb_demo_NOT_EXPERIMENTAL.csv"
        demo_df.to_csv(demo_output, index=False)
        outputs["wkb_demo_csv"] = str(demo_output)

    manifest = {
        "project": "QVSM-RNR-v1",
        "encoding": "(G,L,C,T,B)",
        "ranges": {
            "G": [0, 63],
            "L": [0, 5],
            "C": [0, 5],
            "T": [0, 5],
            "B": [0, 4],
        },
        "total_microstates": 64 * 6 * 6 * 6 * 5,
        "index_formula": "((((G*6+L)*6+C)*6+T)*5+B)",
        "network_nodes": NETWORK_NODES,
        "network_edges": NETWORK_EDGES,
        "geometry_input": str(geometry_path),
        "geometry_input_sha256": geometry_sha256,
        "selected_structure": {
            "PDB_ID": "6W4X",
            "model_index": 0,
            "beta_chain": "C",
            "alpha_chain": "B",
        },
        "epistemic_constraints": [
            "Geometry-derived fields are not H_AB, lambda, DeltaG, or rate constants.",
            "Process-level rates are never written as edge-resolved rates.",
            "Y122 remains a biological node but is geometry-unavailable in selected 6W4X beta-C coordinates.",
            "T is assigned at edge level in this pilot implementation.",
            "NQE/WKB fields are NA unless supplied by explicit QM/MM surface samples.",
            "QED water coherence is not included as a measured parameter.",
        ],
        "references": [
            "10.1021/jacs.0c04325",
            "10.1021/jacs.2c13615",
            "10.1021/jacs.4c13955",
            "10.1021/ja510513z",
        ],
        "outputs": outputs,
    }

    for output_path in outputs.values():
        output_path = Path(output_path)
        manifest.setdefault("output_sha256", {})[
            output_path.name
        ] = sha256_file(output_path)

    manifest_path = output_directory / "qvsm_rnr_v1_manifest.json"
    manifest["manifest_sha256_before_write"] = None

    with open(manifest_path, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)

    print("\n=== QVSM-RNR v1 COMPLETATO ===")
    print(f"Input geometrico: {geometry_path}")
    print(f"SHA-256 input:    {geometry_sha256}")
    print(f"Output directory: {output_directory}")

    print("\n=== COVERAGE NODI ===")
    print(
        node_df[
            [
                "node",
                "codon_DNA",
                "G_index",
                "L_index",
                "C_index_primary",
                "T_index_node",
                "B_index_primary",
            ]
        ].to_string(index=False)
    )

    print("\n=== EDGE FINGERPRINT ===")
    print(
        qvsm_edge_df[
            [
                "edge",
                "source",
                "target",
                "d_centroid_A",
                "d_OO_A",
                "d_min_ring_A",
                "d_min_sidechain_heavy_A",
                "T_index",
                "B_edge_index",
                "k_edge_exp_s-1",
                "k_process_exp_s-1",
                "evidence_scope",
            ]
        ].to_string(index=False)
    )

    print("\n=== MICROSTATES ===")
    print(
        examples_df[
            [
                "node",
                "codon_DNA",
                "G",
                "L",
                "C",
                "T",
                "B",
                "qvsm_index",
                "status",
            ]
        ].to_string(index=False)
    )

    print("\n=== OUTPUT FILES ===")
    for output_path in outputs.values():
        print(output_path)

    print(manifest_path)


if __name__ == "__main__":
    main()

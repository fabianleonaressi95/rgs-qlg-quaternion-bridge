#!/usr/bin/env python3

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


R_KCAL_MOL_K = 0.00198720425864083
T_DEFAULT_K = 298.15
KCAL_TO_J = 4184.0
AVOGADRO = 6.02214076e23
H_PLANCK = 6.62607015e-34
K_BOLTZMANN = 1.380649e-23

NETWORK_EDGES = [
    ("e1", "Y122", "W48"),
    ("e2", "W48", "Y356"),
    ("e3", "Y356", "Y731"),
    ("e4", "Y731", "Y730"),
    ("e5", "Y730", "C439"),
]

# QVSM fixed schema
QVSM_RANGES = {
    "G": (0, 63),
    "L": (0, 5),
    "C": (0, 5),
    "T": (0, 5),
    "B": (0, 4),
}

NEO_REFERENCE = {
    "doi": "10.1021/jacs.4c13955",
    "pdb": "6W4X",
    "system": "E. coli RNR Y731-Y730",
    "method": "NEO-DFT/QM/MM finite-temperature string plus WHAM reweighting",
}

NONADIABATIC_REFERENCE = {
    "doi": "10.1021/acs.jpclett.3c03552",
    "method": "CDFT-CI/Fourier-grid vibronic coupling analysis",
}

DIRECT_PCET_REFERENCE = {
    "doi": "10.1021/jacs.4c15756",
    "method": "vibronically nonadiabatic PCET theory with MD-averaged donor-acceptor distance",
}


def sha256_file(path):
    digest = hashlib.sha256()

    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)

    return digest.hexdigest()


def parse_float(value):
    if value is None:
        return np.nan

    if isinstance(value, float) and np.isnan(value):
        return np.nan

    text = str(value).strip()

    if text in {"", "NA", "None", "nan", "NaN"}:
        return np.nan

    try:
        return float(text)
    except ValueError:
        return np.nan


def qvsm_index(G, L, C, T, B):
    values = [G, L, C, T, B]

    if any(pd.isna(value) for value in values):
        return None

    G, L, C, T, B = map(int, values)

    for name, value in zip(["G", "L", "C", "T", "B"], [G, L, C, T, B]):
        lower, upper = QVSM_RANGES[name]

        if not lower <= value <= upper:
            raise ValueError(f"{name}={value} fuori range {lower}..{upper}")

    return ((((G * 6 + L) * 6 + C) * 6 + T) * 5 + B)


def transition_state_rate(deltaG_dagger_kcal_mol, temperature_K):
    """
    Eyring TST:
    k = k_B T / h * exp(-DeltaG^‡ / RT)
    """
    if pd.isna(deltaG_dagger_kcal_mol):
        return np.nan

    if temperature_K <= 0:
        raise ValueError("Temperature deve essere positiva.")

    prefactor = K_BOLTZMANN * temperature_K / H_PLANCK

    exponent = -deltaG_dagger_kcal_mol / (
        R_KCAL_MOL_K * temperature_K
    )

    return float(prefactor * math.exp(exponent))


def log_eyring_rate(deltaG_dagger_kcal_mol, temperature_K):
    """
    Versione logaritmica per evitare overflow/underflow.
    """
    if pd.isna(deltaG_dagger_kcal_mol):
        return np.nan

    log_prefactor = math.log(
        K_BOLTZMANN * temperature_K / H_PLANCK
    )

    return float(
        log_prefactor
        - deltaG_dagger_kcal_mol
        / (R_KCAL_MOL_K * temperature_K)
    )


def kie_from_barriers(deltaG_H, deltaG_D, temperature_K):
    if pd.isna(deltaG_H) or pd.isna(deltaG_D):
        return np.nan

    return float(
        math.exp(
            (deltaG_D - deltaG_H)
            / (R_KCAL_MOL_K * temperature_K)
        )
    )


def kde_weighted_mean(values, weights):
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)

    mask = np.isfinite(values) & np.isfinite(weights) & (weights >= 0)

    if not np.any(mask):
        return np.nan

    values = values[mask]
    weights = weights[mask]

    if np.sum(weights) == 0:
        return np.nan

    return float(np.sum(values * weights) / np.sum(weights))


def weighted_logsumexp(log_values, weights):
    """
    Restituisce log(sum_i weights_i * exp(log_values_i)).
    """
    log_values = np.asarray(log_values, dtype=float)
    weights = np.asarray(weights, dtype=float)

    mask = (
        np.isfinite(log_values)
        & np.isfinite(weights)
        & (weights > 0)
    )

    if not np.any(mask):
        return np.nan

    x = log_values[mask] + np.log(weights[mask])
    x_max = np.max(x)

    return float(x_max + np.log(np.sum(np.exp(x - x_max))))


def rate_from_configurations(
    deltaG_classical,
    deltaE_nqe,
    weights,
    temperature_K,
):
    """
    k_NEO = < k_TST(config) >_ensemble

    k_config =
      kBT/h * exp[-(DeltaG_classical + DeltaE_NQE)/(RT)]
    """
    deltaG_classical = np.asarray(deltaG_classical, dtype=float)
    deltaE_nqe = np.asarray(deltaE_nqe, dtype=float)
    weights = np.asarray(weights, dtype=float)

    deltaG_quantum = deltaG_classical + deltaE_nqe

    log_rates = np.array([
        log_eyring_rate(value, temperature_K)
        for value in deltaG_quantum
    ])

    log_rate_mean = weighted_logsumexp(log_rates, weights)

    if pd.isna(log_rate_mean):
        return np.nan

    # Evita overflow in output.
    if log_rate_mean > math.log(np.finfo(float).max):
        return np.inf

    return float(math.exp(log_rate_mean))


def load_geometry(path):
    df = pd.read_csv(path)

    expected = {edge[0] for edge in NETWORK_EDGES}
    observed = set(df["edge"].astype(str))

    missing = expected - observed

    if missing:
        raise ValueError(
            f"Archi mancanti nel CSV geometrico: {sorted(missing)}"
        )

    return df


def geometry_lookup(geometry_df, edge_id):
    rows = geometry_df.loc[geometry_df["edge"].eq(edge_id)]

    if len(rows) != 1:
        raise ValueError(
            f"Attesa una riga geometrica per {edge_id}; trovate {len(rows)}"
        )

    return rows.iloc[0]


def build_edge_qvsm_table(geometry_df, temperature_K):
    rows = []

    # Classi T provenienti dall'analisi di nonadiabaticità.
    T_EDGE = {
        "e1": np.nan,
        "e2": np.nan,
        "e3": 1,
        "e4": 3,
        "e5": 2,
    }

    T_LABEL = {
        "e1": "unavailable",
        "e2": "unavailable",
        "e3": "nonadiabatic",
        "e4": "adiabatic",
        "e5": "intermediate",
    }

    for edge_id, source, target in NETWORK_EDGES:
        geom = geometry_lookup(geometry_df, edge_id)

        T = T_EDGE[edge_id]

        row = {
            "edge": edge_id,
            "source": source,
            "target": target,
            "PDB_ID": geom.get("PDB_ID", "6W4X"),
            "geometry_source": geom.get("geometry_source", ""),
            "d_centroid_A": parse_float(
                geom.get("d_centroid_A", np.nan)
            ),
            "d_CA_A": parse_float(
                geom.get("d_CA_A", np.nan)
            ),
            "d_min_ring_A": parse_float(
                geom.get("d_min_ring_A", np.nan)
            ),
            "d_min_sidechain_heavy_A": parse_float(
                geom.get("d_min_sidechain_heavy_A", np.nan)
            ),
            "ring_normal_angle_deg": parse_float(
                geom.get("ring_normal_angle_deg", np.nan)
            ),
            "T_index": T,
            "T_label": T_LABEL[edge_id],
            "T_status": (
                "QM_DERIVED"
                if edge_id in {"e3", "e4", "e5"}
                else "UNAVAILABLE"
            ),
            "temperature_K": temperature_K,
            "deltaG_classical_kcal_mol": np.nan,
            "deltaG_nqe_H_kcal_mol": np.nan,
            "deltaG_nqe_D_kcal_mol": np.nan,
            "Q_H_kcal_mol": np.nan,
            "Q_D_kcal_mol": np.nan,
            "deltaG_quantum_H_kcal_mol": np.nan,
            "deltaG_quantum_D_kcal_mol": np.nan,
            "k_H_s-1": np.nan,
            "k_D_s-1": np.nan,
            "KIE_H_over_D": np.nan,
            "N_NEO_configurations": np.nan,
            "N_DFT_configurations": np.nan,
            "NEO_status": "NO_NEO_SAMPLES_ATTACHED",
            "evidence_scope": (
                "MODEL"
                if edge_id in {"e3", "e4", "e5"}
                else "UNAVAILABLE"
            ),
            "source_reference": (
                "10.1021/jacs.4c13955;"
                "10.1021/acs.jpclett.3c03552"
            ),
            "notes": (
                "No NEO configuration-level data supplied. "
                "All NQE/rate fields remain NA."
            ),
        }

        # Published process/model calibration for e4.
        # These are not inserted as measured edge rates.
        if edge_id == "e4":
            row["deltaG_classical_kcal_mol"] = 4.2
            row["Q_H_kcal_mol"] = np.nan
            row["Q_D_kcal_mol"] = np.nan
            row["NEO_status"] = (
                "PUBLISHED_PROCESS_CALIBRATION_AVAILABLE_BUT_NOT_SPLIT_BY_ISOTOPE"
            )
            row["N_NEO_configurations_alpha"] = 50000
            row["N_NEO_configurations_alpha_prime"] = 85800
            row["notes"] = (
                "Classical QM/MM barrier approximately 4.2 kcal/mol. "
                "Published NQE shifts: +0.8 kcal/mol for E623-mediated alpha "
                "and -4.4 kcal/mol for direct alpha-prime pathway. "
                "These are pathway/model shifts, not native WT edge measurements."
            )

        rows.append(row)

    return pd.DataFrame(rows)


def apply_published_e4_scenarios(edge_df, temperature_K):
    """
    Adds explicit e4 scenario rows:
      alpha E623-mediated: Q = +0.8 kcal/mol
      alpha-prime direct:  Q = -4.4 kcal/mol

    These rows are model/reference scenarios, not WT edge-resolved experimental data.
    """
    base = edge_df.loc[edge_df["edge"].eq("e4")]

    if len(base) != 1:
        raise ValueError("e4 non trovato univocamente.")

    base = base.iloc[0].to_dict()
    scenarios = []

    published = [
        {
            "scenario_id": "e4_alpha_E623_mediated",
            "system": "alpha_pre_turnover",
            "mechanism": "E623-mediated",
            "Q_H_kcal_mol": 0.8,
            "Q_D_kcal_mol": np.nan,
            "Q_status": "PUBLISHED_NQE_SHIFT_NOT_ISOTOPE_RESOLVED",
            "notes": (
                "Published NQE barrier increase for H in E623-mediated "
                "pathway; typical chemical-accuracy caution."
            ),
        },
        {
            "scenario_id": "e4_alpha_prime_direct",
            "system": "alpha_prime_post_turnover",
            "mechanism": "direct",
            "Q_H_kcal_mol": -4.4,
            "Q_D_kcal_mol": np.nan,
            "Q_status": "PUBLISHED_NQE_SHIFT_NOT_ISOTOPE_RESOLVED",
            "notes": (
                "Published NQE barrier decrease for H in direct alpha-prime "
                "pathway; not the ordered pre-turnover alpha pair."
            ),
    ]

    for scenario in published:
        row = dict(base)

        row["scenario_id"] = scenario["scenario_id"]
        row["system"] = scenario["system"]
        row["mechanism"] = scenario["mechanism"]
        row["Q_H_kcal_mol"] = scenario["Q_H_kcal_mol"]
        row["Q_D_kcal_mol"] = scenario["Q_D_kcal_mol"]
        row["Q_status"] = scenario["Q_status"]
        row["deltaG_quantum_H_kcal_mol"] = (
            row["deltaG_classical_kcal_mol"]
            + row["Q_H_kcal_mol"]
        )

        row["k_H_s-1"] = transition_state_rate(
            row["deltaG_quantum_H_kcal_mol"],
            temperature_K,
        )

        row["KIE_H_over_D"] = np.nan
        row["k_D_s-1"] = np.nan
        row["evidence_scope"] = "MODEL"
        row["evidence_status"] = "PUBLISHED_NEO_PATHWAY_SCENARIO"
        row["source_reference"] = "10.1021/jacs.4c13955"
        row["notes"] = scenario["notes"]

        scenarios.append(row)

    return pd.DataFrame(scenarios)


def add_neo_samples(
    edge_df,
    neo_samples_path,
    temperature_K,
):
    """
    Input schema:
      sample_id,edge,configuration_id,delta,R_OO_A,
      deltaG_classical_kcal_mol,deltaE_NQE_H_kcal_mol,
      deltaE_NQE_D_kcal_mol,weight,
      source_reference,evidence_status

    For each configuration:
      G_NEO = G_DFT + DeltaE_NQE

    Then:
      k_H = <k_TST(H)>_weighted
      k_D = <k_TST(D)>_weighted
      KIE = k_H/k_D
    """
    if neo_samples_path is None:
        return pd.DataFrame()

    neo_path = Path(neo_samples_path)

    if not neo_path.exists():
        raise FileNotFoundError(
            f"File NEO samples non trovato: {neo_path}"
        )

    samples = pd.read_csv(neo_path)

    required = {
        "sample_id",
        "edge",
        "configuration_id",
        "delta",
        "R_OO_A",
        "deltaG_classical_kcal_mol",
        "deltaE_NQE_H_kcal_mol",
        "deltaE_NQE_D_kcal_mol",
        "weight",
    }

    missing = required - set(samples.columns)

    if missing:
        raise ValueError(
            f"Colonne NEO mancanti: {sorted(missing)}"
        )

    results = []

    for edge_id in samples["edge"].astype(str).unique():
        subset = samples.loc[samples["edge"].eq(edge_id)].copy()

        classical = subset["deltaG_classical_kcal_mol"].map(parse_float).to_numpy()
        nqe_h = subset["deltaE_NQE_H_kcal_mol"].map(parse_float).to_numpy()
        nqe_d = subset["deltaE_NQE_D_kcal_mol"].map(parse_float).to_numpy()
        weights = subset["weight"].map(parse_float).to_numpy()

        quantum_h = classical + nqe_h
        quantum_d = classical + nqe_d

        k_h = rate_from_configurations(
            classical,
            nqe_h,
            weights,
            temperature_K,
        )

        k_d = rate_from_configurations(
            classical,
            nqe_d,
            weights,
            temperature_K,
        )

        kie = np.nan
        if not pd.isna(k_h) and not pd.isna(k_d) and k_d > 0:
            kie = float(k_h / k_d)

        results.append({
            "edge": edge_id,
            "temperature_K": temperature_K,
            "N_configurations": len(subset),
            "weight_sum": float(np.nansum(weights)),
            "mean_deltaG_classical_kcal_mol": kde_weighted_mean(
                classical,
                weights,
            ),
            "mean_deltaE_NQE_H_kcal_mol": kde_weighted_mean(
                nqe_h,
                weights,
            ),
            "mean_deltaE_NQE_D_kcal_mol": kde_weighted_mean(
                nqe_d,
                weights,
            ),
            "mean_deltaG_quantum_H_kcal_mol": kde_weighted_mean(
                quantum_h,
                weights,
            ),
            "mean_deltaG_quantum_D_kcal_mol": kde_weighted_mean(
                quantum_d,
                weights,
            ),
            "k_H_s-1": k_h,
            "k_D_s-1": k_d,
            "KIE_H_over_D": kie,
            "NEO_status": "CONFIGURATION_LEVEL_REWEIGHTING",
            "source_reference": (
                "10.1021/jacs.4c13955"
                if "source_reference" not in subset.columns
                else ";".join(
                    sorted(
                        set(
                            subset["source_reference"]
                            .dropna()
                            .astype(str)
                        )
                    )
                )
            ),
        })

    return pd.DataFrame(results)


def add_qvsm_mapping(edge_df):
    """
    Adds QVSM edge-level mapping.

    Important:
      T belongs to the edge/process class here.
      G and L are node-level quantities and are not collapsed into
      a single edge index unless an explicit source-state convention is supplied.
    """
    rows = []

    G_NODE = {
        "Y122": [8, 9],
        "W48": [15],
        "Y356": [8, 9],
        "Y731": [8, 9],
        "Y730": [8, 9],
        "C439": [12, 13],
    }

    L_NODE = {
        "Y122": 0,
        "W48": np.nan,
        "Y356": 1,
        "Y731": 2,
        "Y730": 3,
        "C439": 4,
    }

    # Primary C/B values are intentionally provisional.
    C_NODE = {
        "Y122": np.nan,
        "W48": np.nan,
        "Y356": 1,
        "Y731": 2,
        "Y730": 2,
        "C439": 4,
    }

    B_NODE = {
        "Y122": np.nan,
        "W48": np.nan,
        "Y356": 3,
        "Y731": 3,
        "Y730": 4,
        "C439": 4,
    }

    for _, row in edge_df.iterrows():
        edge_id = row["edge"]
        source = row["source"]
        target = row["target"]

        T = parse_float(row.get("T_index", np.nan))
        B = parse_float(row.get("B_edge_index", np.nan))

        for G in G_NODE[source]:
            L = parse_float(L_NODE[source])
            C = parse_float(C_NODE[source])

            index = qvsm_index(G, L, C, T, B)

            rows.append({
                "edge": edge_id,
                "source": source,
                "target": target,
                "source_G": G,
                "source_L": L,
                "source_C": C,
                "edge_T": T,
                "edge_B": B,
                "qvsm_index_source_edge": index,
                "qvsm_index_status": (
                    "COMPUTED"
                    if index is not None
                    else "PARTIAL_UNAVAILABLE"
                ),
            })

    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "QVSM-RNR NEO reduction: "
            "NEO-DFT/DFT energy shifts, barriers, rates, KIE."
        )
    )

    parser.add_argument(
        "--geometry-csv",
        default=(
            "rnr_6W4X_canonical_geometry/"
            "rnr_6W4X_canonical_edges.csv"
        ),
    )

    parser.add_argument(
        "--neo-samples-csv",
        default=None,
        help=(
            "Optional CSV with configuration-level NEO/DFT data."
        ),
    )

    parser.add_argument(
        "--temperature-K",
        type=float,
        default=T_DEFAULT_K,
    )

    parser.add_argument(
        "--outdir",
        default="qvsm_rnr_neo_reduction_output",
    )

    args = parser.parse_args()

    geometry_path = Path(args.geometry_csv).resolve()
    output_directory = Path(args.outdir).resolve()
    output_directory.mkdir(parents=True, exist_ok=True)

    if not geometry_path.exists():
        raise FileNotFoundError(
            f"CSV geometrico non trovato: {geometry_path}"
        )

    geometry_df = load_geometry(geometry_path)

    edge_df = build_edge_qvsm_table(
        geometry_df,
        args.temperature_K,
    )

    e4_scenarios_df = apply_published_e4_scenarios(
        edge_df,
        args.temperature_K,
    )

    neo_results_df = add_neo_samples(
        edge_df,
        args.neo_samples_csv,
        args.temperature_K,
    )

    qvsm_mapping_df = add_qvsm_mapping(edge_df)

    edge_output = output_directory / "qvsm_rnr_neo_edges.csv"
    scenario_output = output_directory / "qvsm_rnr_e4_published_scenarios.csv"
    mapping_output = output_directory / "qvsm_rnr_neo_qvsm_mapping.csv"

    edge_df.to_csv(edge_output, index=False)
    e4_scenarios_df.to_csv(scenario_output, index=False)
    qvsm_mapping_df.to_csv(mapping_output, index=False)

    outputs = {
        "edge_table": str(edge_output),
        "e4_scenarios": str(scenario_output),
        "qvsm_mapping": str(mapping_output),
    }

    if not neo_results_df.empty:
        neo_output = output_directory / "qvsm_rnr_neo_configuration_results.csv"
        neo_results_df.to_csv(neo_output, index=False)
        outputs["neo_configuration_results"] = str(neo_output)

    manifest = {
        "project": "QVSM-RNR-NEO-v1",
        "temperature_K": args.temperature_K,
        "qvsm_encoding": "(G,L,C,T,B)",
        "qvsm_ranges": QVSM_RANGES,
        "total_microstates": 64 * 6 * 6 * 6 * 5,
        "index_formula": "((((G*6+L)*6+C)*6+T)*5+B)",
        "geometry_input": str(geometry_path),
        "geometry_input_sha256": sha256_file(geometry_path),
        "NEO_reference": NEO_REFERENCE,
        "nonadiabatic_reference": NONADIABATIC_REFERENCE,
        "direct_PCET_reference": DIRECT_PCET_REFERENCE,
        "published_NQE_calibrations": {
            "e4_alpha_E623_mediated_Q_H_kcal_mol": 0.8,
            "e4_alpha_prime_direct_Q_H_kcal_mol": -4.4,
            "published_NEO_configuration_counts": {
                "alpha": 50000,
                "alpha_prime": 85800,
            },
        },
        "epistemic_rules": [
            "NEO energy shifts are QM-derived, not experimental edge rates.",
            "Published +0.8 and -4.4 kcal/mol are pathway/scenario values.",
            "The +0.8 value is below typical chemical-accuracy scale and must not be overinterpreted.",
            "KIE requires isotope-resolved Q_H and Q_D; one Q_H shift alone cannot produce a valid KIE.",
            "T is a regime class, not a rate constant.",
            "Geometry fields are not H_AB or lambda.",
            "No QED water-coherence parameter is inserted.",
        ],
        "outputs": outputs,
    }

    manifest["output_sha256"] = {}

    for path in outputs.values():
        path = Path(path)
        manifest["output_sha256"][path.name] = sha256_file(path)

    manifest_path = output_directory / "qvsm_rnr_neo_manifest.json"

    with open(manifest_path, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)

    print("\n=== QVSM-RNR NEO REDUCTION COMPLETATA ===")
    print(f"Input geometrico: {geometry_path}")
    print(f"Temperatura:      {args.temperature_K} K")
    print(f"Output directory: {output_directory}")

    print("\n=== EDGE TABLE ===")
    print(
        edge_df[
            [
                "edge",
                "source",
                "target",
                "T_index",
                "T_label",
                "deltaG_classical_kcal_mol",
                "Q_H_kcal_mol",
                "k_H_s-1",
                "KIE_H_over_D",
                "NEO_status",
            ]
        ].to_string(index=False)
    )

    print("\n=== PUBLISED e4 SCENARIOS ===")
    print(
        e4_scenarios_df[
            [
                "scenario_id",
                "system",
                "mechanism",
                "deltaG_classical_kcal_mol",
                "Q_H_kcal_mol",
                "deltaG_quantum_H_kcal_mol",
                "k_H_s-1",
                "Q_status",
            ]
        ].to_string(index=False)
    )

    if not neo_results_df.empty:
        print("\n=== NEO CONFIGURATION RESULTS ===")
        print(neo_results_df.to_string(index=False))

    print("\n=== OUTPUT ===")
    for path in outputs.values():
        print(path)
    print(manifest_path)


if __name__ == "__main__":
    main()

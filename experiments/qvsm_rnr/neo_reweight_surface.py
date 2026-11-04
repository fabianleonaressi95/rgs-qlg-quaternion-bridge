#!/usr/bin/env python3

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


R_KCAL_MOL_K = 0.00198720425864083


def sha256_file(path):
    import hashlib

    h = hashlib.sha256()

    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)

    return h.hexdigest()


def logsumexp(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]

    if len(values) == 0:
        return np.nan

    m = np.max(values)
    return float(m + np.log(np.sum(np.exp(values - m))))


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--samples", required=True)
    parser.add_argument("--outdir", default="neo_surface_output")
    parser.add_argument("--temperature-K", type=float, default=298.15)
    parser.add_argument("--R-bin-A", type=float, default=0.05)
    parser.add_argument("--delta-bin-A", type=float, default=0.05)

    args = parser.parse_args()

    samples_path = Path(args.samples).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    required = {
        "sample_id",
        "edge",
        "R_OO_A",
        "delta_A",
        "weight_classical",
        "E_DFT_kcal_mol",
        "E_NEO_H_kcal_mol",
        "E_NEO_D_kcal_mol",
    }

    df = pd.read_csv(samples_path)
    missing = required - set(df.columns)

    if missing:
        raise ValueError(f"Colonne mancanti: {sorted(missing)}")

    for column in required - {"sample_id", "edge"}:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    df = df.dropna(
        subset=[
            "R_OO_A",
            "delta_A",
            "weight_classical",
            "E_DFT_kcal_mol",
            "E_NEO_H_kcal_mol",
            "E_NEO_D_kcal_mol",
        ]
    ).copy()

    df = df.loc[df["weight_classical"] > 0].copy()

    if df.empty:
        raise ValueError("Nessuna configurazione valida dopo i filtri.")

    beta = 1.0 / (R_KCAL_MOL_K * args.temperature_K)

    df["deltaE_NQE_H_kcal_mol"] = (
        df["E_NEO_H_kcal_mol"] - df["E_DFT_kcal_mol"]
    )

    df["deltaE_NQE_D_kcal_mol"] = (
        df["E_NEO_D_kcal_mol"] - df["E_DFT_kcal_mol"]
    )

    df["R_bin"] = (
        np.floor(df["R_OO_A"] / args.R_bin_A)
        * args.R_bin_A
    )

    df["delta_bin"] = (
        np.floor(df["delta_A"] / args.delta_bin_A)
        * args.delta_bin_A
    )

    rows = []

    for (edge, R_bin, delta_bin), group in df.groupby(
        ["edge", "R_bin", "delta_bin"],
        sort=True,
    ):
        w_cl = group["weight_classical"].to_numpy(dtype=float)
        dE_H = group["deltaE_NQE_H_kcal_mol"].to_numpy(dtype=float)
        dE_D = group["deltaE_NQE_D_kcal_mol"].to_numpy(dtype=float)

        P_cl = np.sum(w_cl)
        P_H = np.sum(w_cl * np.exp(-beta * dE_H))
        P_D = np.sum(w_cl * np.exp(-beta * dE_D))

        rows.append({
            "edge": edge,
            "R_bin_A": R_bin,
            "delta_bin_A": delta_bin,
            "N_samples": len(group),
            "P_classical_unnormalized": P_cl,
            "P_H_unnormalized": P_H,
            "P_D_unnormalized": P_D,
            "mean_deltaE_NQE_H_kcal_mol": np.average(
                dE_H,
                weights=w_cl,
            ),
            "mean_deltaE_NQE_D_kcal_mol": np.average(
                dE_D,
                weights=w_cl,
            ),
        })

    surface = pd.DataFrame(rows)

    for edge in surface["edge"].unique():
        mask = surface["edge"].eq(edge)

        for label, probability_column in [
            ("classical", "P_classical_unnormalized"),
            ("H", "P_H_unnormalized"),
            ("D", "P_D_unnormalized"),
        ]:
            p = surface.loc[mask, probability_column].to_numpy(dtype=float)
            p = p / np.sum(p)

            G = -R_KCAL_MOL_K * args.temperature_K * np.log(p)
            G = G - np.nanmin(G)

            surface.loc[mask, f"G_{label}_kcal_mol"] = G

    surface_csv = outdir / "neo_reweighted_surface.csv"
    samples_csv = outdir / "neo_samples_with_corrections.csv"

    surface.to_csv(surface_csv, index=False)
    df.to_csv(samples_csv, index=False)

    manifest = {
        "input_samples": str(samples_path),
        "input_sha256": sha256_file(samples_path),
        "temperature_K": args.temperature_K,
        "R_bin_A": args.R_bin_A,
        "delta_bin_A": args.delta_bin_A,
        "formula": {
            "deltaE_NQE_H": "E_NEO_H - E_DFT",
            "deltaE_NQE_D": "E_NEO_D - E_DFT",
            "P_H": "sum(w_classical * exp(-beta*deltaE_NQE_H))",
            "P_D": "sum(w_classical * exp(-beta*deltaE_NQE_D))",
            "G": "-RT ln(P) + constant",
        },
        "outputs": {
            "surface_csv": str(surface_csv),
            "sample_corrections_csv": str(samples_csv),
            "surface_sha256": sha256_file(surface_csv),
            "samples_sha256": sha256_file(samples_csv),
        },
    }

    manifest_path = outdir / "neo_reweight_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    print("Creati:")
    print(surface_csv)
    print(samples_csv)
    print(manifest_path)


if __name__ == "__main__":
    main()

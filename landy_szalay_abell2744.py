#!/usr/bin/env python3
"""
Landy-Szalay 2-point correlation function for an Abell 2744 galaxy catalog.

Supports:
  1) projected coordinates x_kpc,y_kpc (recommended for simulated catalogues);
  2) sky coordinates ra_deg,dec_deg (or common aliases), using angular separations.

Outputs:
  - CSV with bin centers, xi_LS, and delete-one spatial jackknife errors
  - PNG plot with error bars

Dependencies: numpy, pandas, scipy, matplotlib
Example:
  python landy_szalay_abell2744.py abell2744_galaxies_catalog.json \
      --xcol x_kpc --ycol y_kpc --bins 12 --random-factor 20 --jackknife 8

For sky coordinates:
  python landy_szalay_abell2744.py catalogue.csv --xcol ra_deg --ycol dec_deg \
      --bins 12 --random-factor 20 --jackknife 8

Important:
- The random catalogue must represent the SAME survey footprint/selection function.
- If the footprint is not a rectangle, supply a pre-generated random catalogue via
  --random-file with the same coordinate columns.
- Jackknife errors quantify spatial resampling variance; they do not replace a
  full covariance/selection-function analysis.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree


ALIASES_X = ["x_kpc", "x", "x_mpc", "ra_deg", "ra", "RA", "RAJ2000", "ra_j2000"]
ALIASES_Y = ["y_kpc", "y", "y_mpc", "dec_deg", "dec", "DEC", "DEJ2000", "dec_j2000"]


def load_catalog(path):
    path = Path(path)
    if path.suffix.lower() == ".json":
        obj = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(obj, dict):
            # Common container keys
            for key in ("galaxies", "objects", "catalog", "data", "sources"):
                if key in obj and isinstance(obj[key], list):
                    obj = obj[key]
                    break
        if not isinstance(obj, list):
            raise ValueError("JSON deve contenere una lista di oggetti o una chiave galaxies/objects/catalog/data/sources.")
        df = pd.DataFrame(obj)
    elif path.suffix.lower() in (".csv", ".txt"):
        df = pd.read_csv(path)
    else:
        raise ValueError("Formato input non supportato: usare JSON o CSV.")
    if df.empty:
        raise ValueError("Catalogo vuoto.")
    return df


def choose_column(df, requested, aliases, label):
    if requested:
        if requested not in df.columns:
            raise ValueError(f"Colonna '{requested}' assente. Colonne disponibili: {list(df.columns)}")
        return requested
    for c in aliases:
        if c in df.columns:
            return c
    raise ValueError(f"Impossibile identificare la colonna {label}. Specificarla con --xcol/--ycol.")


def to_xy(df, xcol, ycol, sky):
    xy = df[[xcol, ycol]].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    good = np.isfinite(xy).all(axis=1)
    xy = xy[good]
    if len(xy) < 10:
        raise ValueError(f"Solo {len(xy)} posizioni valide: servono almeno 10.")
    if sky:
        ra = np.deg2rad(xy[:, 0])
        dec = np.deg2rad(xy[:, 1])
        # Unit-sphere Cartesian coordinates: Euclidean chord distances give
        # exact angular separations via 2 asin(chord/2).
        return np.column_stack((np.cos(dec)*np.cos(ra),
                               np.cos(dec)*np.sin(ra),
                               np.sin(dec)))
    return xy


def pair_counts(data_a, data_b, edges, auto=False, sky=False):
    """Count pairs in separation bins using cKDTree; normalize to possible pairs."""
    if sky:
        # Chord-distance edges corresponding to angular edges in radians.
        edges_metric = 2.0 * np.sin(np.asarray(edges) / 2.0)
    else:
        edges_metric = np.asarray(edges)
    tree_a = cKDTree(data_a)
    tree_b = tree_a if auto else cKDTree(data_b)
    counts = np.zeros(len(edges_metric) - 1, dtype=float)
    for i in range(len(counts)):
        # count_neighbors is cumulative, so difference adjacent cumulative counts
        pass
    cumulative = tree_a.count_neighbors(tree_b, edges_metric, cumulative=True)
    cumulative = np.asarray(cumulative, dtype=float)
    if auto:
        # count_neighbors(auto) counts ordered pairs and self-pairs; convert
        # to unordered distinct pairs.
        cumulative = (cumulative - len(data_a)) / 2.0
        norm = len(data_a) * (len(data_a) - 1) / 2.0
    else:
        norm = len(data_a) * len(data_b)
    hist = np.diff(cumulative)
    if norm <= 0:
        raise ValueError("Numero di coppie insufficiente.")
    return hist / norm


def landy_szalay(data, random, edges, sky=False):
    dd = pair_counts(data, data, edges, auto=True, sky=sky)
    dr = pair_counts(data, random, edges, auto=False, sky=sky)
    rr = pair_counts(random, random, edges, auto=True, sky=sky)
    xi = np.full_like(rr, np.nan)
    ok = rr > 0
    xi[ok] = (dd[ok] - 2.0 * dr[ok] + rr[ok]) / rr[ok]
    return xi, dd, dr, rr


def random_in_rectangle(data, n, rng):
    lo = np.min(data, axis=0)
    hi = np.max(data, axis=0)
    if np.any(hi <= lo):
        raise ValueError("Estensione nulla in una o più coordinate.")
    return rng.uniform(lo, hi, size=(n, data.shape[1]))


def make_regions(data, n_regions, sky=False):
    """Create spatial jackknife regions by quantile-binning x (or RA-like axis)."""
    axis = data[:, 0]
    # For sky vectors, recover RA from xyz.
    if sky:
        axis = np.arctan2(data[:, 1], data[:, 0])
    cuts = np.quantile(axis, np.linspace(0, 1, n_regions + 1))
    labels = np.digitize(axis, cuts[1:-1], right=False)
    labels = np.clip(labels, 0, n_regions - 1)
    return labels


def main():
    p = argparse.ArgumentParser()
    p.add_argument("catalog")
    p.add_argument("--xcol", default=None)
    p.add_argument("--ycol", default=None)
    p.add_argument("--bins", type=int, default=12)
    p.add_argument("--rmin", type=float, default=None)
    p.add_argument("--rmax", type=float, default=None)
    p.add_argument("--random-factor", type=int, default=20)
    p.add_argument("--random-file", default=None)
    p.add_argument("--random-xcol", default=None)
    p.add_argument("--random-ycol", default=None)
    p.add_argument("--jackknife", type=int, default=8)
    p.add_argument("--seed", type=int, default=20261021)
    p.add_argument("--out", default="abell2744_landy_szalay")
    p.add_argument("--sky", action="store_true",
                   help="Interpret coordinates as RA/Dec in degrees; angular bins are in arcmin.")
    args = p.parse_args()

    df = load_catalog(args.catalog)
    xcol = choose_column(df, args.xcol, ALIASES_X, "x/RA")
    ycol = choose_column(df, args.ycol, ALIASES_Y, "y/Dec")
    data = to_xy(df, xcol, ycol, args.sky)
    ndata = len(data)
    rng = np.random.default_rng(args.seed)

    if args.random_file:
        rdf = load_catalog(args.random_file)
        rx = choose_column(rdf, args.random_xcol or args.xcol, ALIASES_X, "random x/RA")
        ry = choose_column(rdf, args.random_ycol or args.ycol, ALIASES_Y, "random y/Dec")
        random = to_xy(rdf, rx, ry, args.sky)
    else:
        if args.sky:
            # Rectangular RA/Dec sampling is only a small-field approximation;
            # explicit random catalogue is preferred for a real survey footprint.
            raw = df[[xcol, ycol]].apply(pd.to_numeric, errors="coerce").dropna().to_numpy(float)
            rraw = random_in_rectangle(raw, args.random_factor * len(raw), rng)
            random = to_xy(pd.DataFrame(rraw, columns=[xcol, ycol]), xcol, ycol, True)
        else:
            random = random_in_rectangle(data, args.random_factor * ndata, rng)

    if args.sky:
        # arcmin -> radians; cap at <180 degrees
        valid = np.clip(data[:, 2], -1, 1)
        # all angular pair separations in arcmin; use data's angular extent
        from scipy.spatial.distance import pdist
        max_chord = np.max(pdist(data)) if len(data) > 1 else 0
        max_ang_arcmin = np.rad2deg(2*np.arcsin(min(1.0, max_chord/2))) * 60
        rmin = args.rmin if args.rmin is not None else max(max_ang_arcmin / 1000, 0.05)
        rmax = args.rmax if args.rmax is not None else max_ang_arcmin * 0.95
        edges = np.deg2rad(np.geomspace(rmin, rmax, args.bins + 1) / 60.0)
        unit = "arcmin"
        centers = np.sqrt(edges[:-1] * edges[1:]) * 180/np.pi * 60
    else:
        from scipy.spatial.distance import pdist
        maxdist = float(np.max(pdist(data)))
        rmin = args.rmin if args.rmin is not None else maxdist / 1000
        rmax = args.rmax if args.rmax is not None else maxdist * 0.95
        if rmin <= 0 or rmax <= rmin:
            raise ValueError("Servono 0 < rmin < rmax.")
        edges = np.geomspace(rmin, rmax, args.bins + 1)
        unit = "coordinate_units"
        centers = np.sqrt(edges[:-1] * edges[1:])

    xi, dd, dr, rr = landy_szalay(data, random, edges, sky=args.sky)

    # Delete-one jackknife: divide catalogue spatially, remove one region at a time.
    jk = []
    if args.jackknife >= 3:
        labels = make_regions(data, args.jackknife, sky=args.sky)
        # Assign random points to the same spatial regions using the same coordinate
        # convention. For 2D projected data, quantile boundaries are estimated from data.
        if args.sky:
            d_axis = np.arctan2(data[:, 1], data[:, 0])
            r_axis = np.arctan2(random[:, 1], random[:, 0])
        else:
            d_axis, r_axis = data[:, 0], random[:, 0]
        cuts = np.quantile(d_axis, np.linspace(0, 1, args.jackknife + 1))
        rlabels = np.clip(np.digitize(r_axis, cuts[1:-1]), 0, args.jackknife - 1)
        for j in range(args.jackknife):
            dsub = data[labels != j]
            rsub = random[rlabels != j]
            if len(dsub) < 3 or len(rsub) < 3:
                continue
            xij, _, _, _ = landy_szalay(dsub, rsub, edges, sky=args.sky)
            jk.append(xij)
    jk = np.asarray(jk)
    if len(jk) >= 3:
        jkmean = np.nanmean(jk, axis=0)
        err = np.sqrt((len(jk)-1)/len(jk) * np.nansum((jk-jkmean)**2, axis=0))
    else:
        err = np.full_like(xi, np.nan)

    out_csv = Path(args.out + ".csv")
    out_png = Path(args.out + ".png")
    result = pd.DataFrame({
        "r_center": centers,
        "xi_landy_szalay": xi,
        "jackknife_sigma": err,
        "DD_norm": dd,
        "DR_norm": dr,
        "RR_norm": rr,
        "n_galaxies_used": ndata,
        "n_random_used": len(random),
        "coordinate_unit": unit,
    })
    result.to_csv(out_csv, index=False)

    plt.figure(figsize=(8, 5.5))
    plt.errorbar(centers, xi, yerr=err, fmt="o-", capsize=3, label="Landy–Szalay")
    plt.axhline(0, linewidth=1)
    plt.xscale("log")
    plt.xlabel(f"Separazione r [{unit}]")
    plt.ylabel(r"$\hat{\xi}(r)$")
    plt.title("Abell 2744 — funzione di correlazione a due punti")
    plt.grid(True, which="both", alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=180)
    print(f"Galassie valide: {ndata}")
    print(f"Random: {len(random)}")
    print(f"Coordinate: {'RA/Dec (angolare)' if args.sky else 'proiettate'}")
    print(f"CSV: {out_csv.resolve()}")
    print(f"Plot: {out_png.resolve()}")
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()

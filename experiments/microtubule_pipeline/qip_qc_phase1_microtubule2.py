"""
QIP-QC PHASE 1B
============================================================

REAL MICROTUBULE -> PROTOFILAMENT / DIMER GEOMETRY

Input:
    PDB 5IJ0

Goal:
    Extract alpha/beta tubulin units
    reconstruct geometric ordering
    build structural signals:
        d_n       = longitudinal spacing
        theta_n   = local orientation/twist
        r_n       = radial coordinate

Compare:
    1. Real microtubule
    2. Periodic surrogate
    3. Fibonacci surrogate
    4. Random surrogate
    5. Fourier phase-randomized surrogate

Important:
    This script does NOT assume the microtubule is
    quasicrystalline.

It tests whether its geometry contains
aperiodic/quasiperiodic signatures.
"""

import os
import urllib.request
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from Bio.PDB import PDBParser


# ============================================================
# CONFIG
# ============================================================

PDB_ID = "5IJ0"

DATA_DIR = "qip_qc_data"
FIG_DIR = "qip_qc_figures"

PDB_FILE = os.path.join(
    DATA_DIR,
    f"{PDB_ID}.pdb"
)

os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

RANDOM_SEED = 42


# ============================================================
# DOWNLOAD
# ============================================================

def download_pdb():

    if os.path.exists(PDB_FILE):
        print("[OK] PDB already exists")
        return

    url = (
        f"https://files.rcsb.org/download/"
        f"{PDB_ID}.pdb"
    )

    print("[DOWNLOAD]", url)

    urllib.request.urlretrieve(
        url,
        PDB_FILE
    )

    print("[OK]", PDB_FILE)


# ============================================================
# LOAD
# ============================================================

def load_structure():

    parser = PDBParser(QUIET=True)

    structure = parser.get_structure(
        PDB_ID,
        PDB_FILE
    )

    return structure


# ============================================================
# IDENTIFY TUBULIN CHAINS
# ============================================================

def extract_chains(structure):

    model = next(
        structure.get_models()
    )

    chains = []

    for chain in model:

        residues = []

        for residue in chain:

            if residue.id[0] != " ":
                continue

            if "CA" not in residue:
                continue

            residues.append(
                residue
            )

        if len(residues) < 50:
            continue

        coords = np.array(
            [
                r["CA"].coord
                for r in residues
            ]
        )

        chains.append({
            "id": chain.id,
            "residues": residues,
            "coords": coords
        })

    print(
        f"[CHAINS] usable chains = {len(chains)}"
    )

    for c in chains:

        print(
            f"  Chain {c['id']}: "
            f"{len(c['residues'])} residues"
        )

    return chains


# ============================================================
# CHAIN CENTROID
# ============================================================

def centroid(coords):

    return np.mean(
        coords,
        axis=0
    )


# ============================================================
# CHAIN AXIS
# ============================================================

def principal_axis(coords):

    """
    Principal component of the chain coordinates.
    """

    x = coords - np.mean(
        coords,
        axis=0
    )

    cov = np.cov(
        x.T
    )

    values, vectors = np.linalg.eigh(
        cov
    )

    axis = vectors[
        :,
        np.argmax(values)
    ]

    axis = axis / np.linalg.norm(axis)

    return axis


# ============================================================
# PROJECT RESIDUES ONTO AXIS
# ============================================================

def order_chain(coords):

    axis = principal_axis(
        coords
    )

    center = np.mean(
        coords,
        axis=0
    )

    projection = (
        coords - center
    ) @ axis

    order = np.argsort(
        projection
    )

    return order, axis, projection


# ============================================================
# EXTRACT RESIDUE SIGNAL
# ============================================================

def extract_ordered_chain_signal(chain):

    coords = chain["coords"]

    order, axis, projection = (
        order_chain(coords)
    )

    ordered = coords[
        order
    ]

    return ordered, axis, projection[order]


# ============================================================
# LONGITUDINAL SPACING
# ============================================================

def longitudinal_spacing(
    ordered_coords,
    axis
):

    projected = (
        ordered_coords
        @ axis
    )

    d = np.diff(
        projected
    )

    return np.abs(d)


# ============================================================
# RADIAL DISTANCE
# ============================================================

def radial_signal(
    ordered_coords
):

    center = np.mean(
        ordered_coords,
        axis=0
    )

    radial = np.linalg.norm(
        ordered_coords - center,
        axis=1
    )

    return radial


# ============================================================
# LOCAL ORIENTATION
# ============================================================

def local_orientation(
    ordered_coords
):

    vectors = np.diff(
        ordered_coords,
        axis=0
    )

    norms = np.linalg.norm(
        vectors,
        axis=1,
        keepdims=True
    )

    norms[norms == 0] = 1

    vectors = (
        vectors / norms
    )

    angles = []

    for i in range(
        len(vectors) - 1
    ):

        a = vectors[i]
        b = vectors[i + 1]

        c = np.clip(
            np.dot(a, b),
            -1.0,
            1.0
        )

        angles.append(
            np.arccos(c)
        )

    return np.array(
        angles
    )


# ============================================================
# FIBONACCI WORD
# ============================================================

def fibonacci_word(n):

    word = "L"

    while len(word) < n:

        new_word = ""

        for c in word:

            if c == "L":
                new_word += "LS"
            else:
                new_word += "L"

        word = new_word

    return word[:n]


# ============================================================
# FIBONACCI NUMERIC
# ============================================================

def fibonacci_numeric(
    n,
    long_value=1.0,
    short_value=0.61803398875
):

    word = fibonacci_word(
        n
    )

    return np.array(
        [
            long_value
            if c == "L"
            else short_value
            for c in word
        ]
    )


# ============================================================
# PERIODIC CONTROL
# ============================================================

def periodic_signal(
    n,
    mean,
    amplitude
):

    x = np.arange(n)

    return (
        mean
        +
        amplitude
        *
        np.cos(
            2 * np.pi * x / 2
        )
    )


# ============================================================
# RANDOM CONTROL
# ============================================================

def random_signal(
    n,
    mean,
    std,
    seed=42
):

    rng = np.random.default_rng(
        seed
    )

    return rng.normal(
        mean,
        std,
        size=n
    )


# ============================================================
# FOURIER PHASE SURROGATE
# ============================================================

def fourier_surrogate(
    x,
    seed=42
):

    rng = np.random.default_rng(
        seed
    )

    x = np.asarray(x)

    x0 = x - np.mean(x)

    fft = np.fft.rfft(
        x0
    )

    amplitude = np.abs(
        fft
    )

    phase = np.angle(
        fft
    )

    random_phase = rng.uniform(
        0,
        2 * np.pi,
        size=len(fft)
    )

    random_phase[0] = phase[0]

    if len(x) % 2 == 0:

        random_phase[-1] = phase[-1]

    new_fft = (
        amplitude
        *
        np.exp(
            1j * random_phase
        )
    )

    surrogate = np.fft.irfft(
        new_fft,
        n=len(x)
    )

    return surrogate + np.mean(x)


# ============================================================
# AUTOCORRELATION
# ============================================================

def autocorrelation(x):

    x = np.asarray(x)

    x = x - np.mean(x)

    corr = np.correlate(
        x,
        x,
        mode="full"
    )

    corr = corr[
        len(corr)//2:
    ]

    if abs(corr[0]) > 1e-14:

        corr /= corr[0]

    return corr


# ============================================================
# NORMALIZED PSD
# ============================================================

def power_spectrum(x):

    x = np.asarray(x)

    x = x - np.mean(x)

    fft = np.fft.rfft(
        x
    )

    power = np.abs(
        fft
    ) ** 2

    if np.max(power) > 0:

        power /= np.max(
            power
        )

    freq = np.fft.rfftfreq(
        len(x)
    )

    return freq, power


# ============================================================
# SPECTRAL ENTROPY
# ============================================================

def spectral_entropy(x):

    x = np.asarray(x)

    x = x - np.mean(x)

    power = np.abs(
        np.fft.rfft(x)
    ) ** 2

    power = power[
        1:
    ]

    if np.sum(power) <= 0:

        return np.nan

    p = power / np.sum(
        power
    )

    entropy = -np.sum(
        p * np.log(
            p + 1e-15
        )
    )

    entropy /= np.log(
        len(p)
    )

    return entropy


# ============================================================
# FIBONACCI MATCHING
# ============================================================

def fibonacci_distance(
    x
):

    """
    Fit the real signal to a two-level
    Fibonacci structural template.

    We optimize the two levels directly.
    """

    x = np.asarray(x)

    n = len(x)

    fib = fibonacci_word(
        n
    )

    mask_L = np.array(
        [
            c == "L"
            for c in fib
        ]
    )

    mask_S = ~mask_L

    if (
        np.sum(mask_L) == 0
        or np.sum(mask_S) == 0
    ):
        return np.nan

    L = np.mean(
        x[mask_L]
    )

    S = np.mean(
        x[mask_S]
    )

    model = np.where(
        mask_L,
        L,
        S
    )

    rmse = np.sqrt(
        np.mean(
            (x - model) ** 2
        )
    )

    normalized = (
        rmse
        /
        (
            np.std(x)
            + 1e-12
        )
    )

    return normalized


# ============================================================
# PERIODIC FIT
# ============================================================

def periodic_distance(
    x
):

    x = np.asarray(x)

    n = len(x)

    t = np.arange(n)

    best = np.inf

    for period in np.linspace(
        1.5,
        max(2.0, n / 2),
        200
    ):

        omega = (
            2 * np.pi
            / period
        )

        X = np.column_stack(
            [
                np.ones(n),
                np.cos(omega * t),
                np.sin(omega * t)
            ]
        )

        beta, *_ = np.linalg.lstsq(
            X,
            x,
            rcond=None
        )

        pred = X @ beta

        rmse = np.sqrt(
            np.mean(
                (x - pred) ** 2
            )
        )

        if rmse < best:

            best = rmse

    return (
        best
        /
        (
            np.std(x)
            + 1e-12
        )
    )


# ============================================================
# PROCESS ONE CHAIN
# ============================================================

def process_chain(chain):

    ordered, axis, projection = (
        extract_ordered_chain_signal(
            chain
        )
    )

    spacing = longitudinal_spacing(
        ordered,
        axis
    )

    radial = radial_signal(
        ordered
    )

    angles = local_orientation(
        ordered
    )

    return {
        "ordered": ordered,
        "axis": axis,
        "projection": projection,
        "spacing": spacing,
        "radial": radial,
        "angles": angles
    }


# ============================================================
# PLOT SIGNAL
# ============================================================

def plot_signal(
    signal,
    name,
    filename
):

    plt.figure(
        figsize=(12, 5)
    )

    plt.plot(
        signal,
        linewidth=1
    )

    plt.xlabel(
        "Ordered structural index"
    )

    plt.ylabel(
        name
    )

    plt.title(
        f"Microtubule structural signal: {name}"
    )

    plt.tight_layout()

    path = os.path.join(
        FIG_DIR,
        filename
    )

    plt.savefig(
        path,
        dpi=200
    )

    plt.close()

    print(
        "[FIGURE]",
        path
    )


# ============================================================
# COMPARISON PLOT
# ============================================================

def plot_comparison(
    x,
    name
):

    n = len(x)

    mean = np.mean(x)
    std = np.std(x)

    fib = fibonacci_numeric(
        n,
        mean + std,
        mean - std
    )

    periodic = periodic_signal(
        n,
        mean,
        std
    )

    random = random_signal(
        n,
        mean,
        std,
        RANDOM_SEED
    )

    surrogate = fourier_surrogate(
        x,
        RANDOM_SEED
    )

    plt.figure(
        figsize=(14, 8)
    )

    plt.plot(
        x,
        label="Microtubule",
        linewidth=1.5
    )

    plt.plot(
        fib,
        label="Fibonacci",
        alpha=0.8
    )

    plt.plot(
        periodic,
        label="Periodic",
        alpha=0.8
    )

    plt.plot(
        random,
        label="Random",
        alpha=0.6
    )

    plt.plot(
        surrogate,
        label="Fourier surrogate",
        alpha=0.8
    )

    plt.xlabel(
        "Structural index"
    )

    plt.ylabel(
        name
    )

    plt.title(
        f"Structural comparison — {name}"
    )

    plt.legend()

    plt.tight_layout()

    path = os.path.join(
        FIG_DIR,
        f"comparison_{name}.png"
    )

    plt.savefig(
        path,
        dpi=200
    )

    plt.close()

    print(
        "[FIGURE]",
        path
    )


# ============================================================
# REPORT
# ============================================================

def analyze_signal(
    x,
    name
):

    ac = autocorrelation(
        x
    )

    entropy = spectral_entropy(
        x
    )

    fib_d = fibonacci_distance(
        x
    )

    periodic_d = periodic_distance(
        x
    )

    print("\n")
    print("-" * 70)
    print(name)
    print("-" * 70)

    print(
        f"N = {len(x)}"
    )

    print(
        f"mean = {np.mean(x):.8f}"
    )

    print(
        f"std  = {np.std(x):.8f}"
    )

    if len(ac) > 20:

        short_corr = np.mean(
            np.abs(
                ac[10:20]
            )
        )

    else:

        short_corr = np.nan

    print(
        f"short-range correlation = "
        f"{short_corr:.8f}"
    )

    print(
        f"spectral entropy = "
        f"{entropy:.8f}"
    )

    print(
        f"Fibonacci normalized error = "
        f"{fib_d:.8f}"
    )

    print(
        f"Periodic normalized error = "
        f"{periodic_d:.8f}"
    )

    return {
        "N": len(x),
        "mean": np.mean(x),
        "std": np.std(x),
        "short_corr": short_corr,
        "spectral_entropy": entropy,
        "fibonacci_error": fib_d,
        "periodic_error": periodic_d
    }


# ============================================================
# MAIN
# ============================================================

def main():

    print("=" * 70)

    print(
        "QIP-QC PHASE 1B"
    )

    print(
        "MICROTUBULE PROTOFILAMENT GEOMETRY"
    )

    print("=" * 70)

    # --------------------------------------------------------
    # Download
    # --------------------------------------------------------

    download_pdb()

    # --------------------------------------------------------
    # Load
    # --------------------------------------------------------

    structure = load_structure()

    # --------------------------------------------------------
    # Chains
    # --------------------------------------------------------

    chains = extract_chains(
        structure
    )

    if len(chains) == 0:

        raise RuntimeError(
            "No suitable chains found."
        )

    # --------------------------------------------------------
    # Process all chains
    # --------------------------------------------------------

    results = []

    for chain in chains:

        print("\n")
        print(
            "=" * 50
        )

        print(
            f"PROCESSING CHAIN {chain['id']}"
        )

        print(
            "=" * 50
        )

        result = process_chain(
            chain
        )

        results.append(
            (
                chain["id"],
                result
            )
        )

        spacing = result[
            "spacing"
        ]

        radial = result[
            "radial"
        ]

        angles = result[
            "angles"
        ]

        # ----------------------------------------------------
        # Analyze
        # ----------------------------------------------------

        spacing_report = analyze_signal(
            spacing,
            f"Chain {chain['id']} longitudinal spacing"
        )

        radial_report = analyze_signal(
            radial,
            f"Chain {chain['id']} radial coordinate"
        )

        angle_report = analyze_signal(
            angles,
            f"Chain {chain['id']} local angle"
        )

        # ----------------------------------------------------
        # Plots
        # ----------------------------------------------------

        plot_signal(
            spacing,
            "Spacing (Å)",
            f"chain_{chain['id']}_spacing.png"
        )

        plot_signal(
            radial,
            "Radial coordinate (Å)",
            f"chain_{chain['id']}_radial.png"
        )

        plot_signal(
            angles,
            "Local orientation (rad)",
            f"chain_{chain['id']}_angles.png"
        )

        plot_comparison(
            spacing,
            "spacing"
        )

        plot_comparison(
            radial,
            "radial"
        )

        plot_comparison(
            angles,
            "angles"
        )

    # ========================================================
    # COMBINE CHAINS
    # ========================================================

    print("\n")
    print("=" * 70)
    print("COMBINED STRUCTURAL ANALYSIS")
    print("=" * 70)

    all_spacing = np.concatenate(
        [
            r["spacing"]
            for _, r in results
        ]
    )

    all_radial = np.concatenate(
        [
            r["radial"]
            for _, r in results
        ]
    )

    all_angles = np.concatenate(
        [
            r["angles"]
            for _, r in results
        ]
    )

    reports = {}

    reports["spacing"] = analyze_signal(
        all_spacing,
        "ALL CHAINS — spacing"
    )

    reports["radial"] = analyze_signal(
        all_radial,
        "ALL CHAINS — radial"
    )

    reports["angles"] = analyze_signal(
        all_angles,
        "ALL CHAINS — angles"
    )

    # ========================================================
    # SAVE DATA
    # ========================================================

    np.save(
        os.path.join(
            DATA_DIR,
            "phase1b_spacing.npy"
        ),
        all_spacing
    )

    np.save(
        os.path.join(
            DATA_DIR,
            "phase1b_radial.npy"
        ),
        all_radial
    )

    np.save(
        os.path.join(
            DATA_DIR,
            "phase1b_angles.npy"
        ),
        all_angles
    )

    # ========================================================
    # FINAL SUMMARY
    # ========================================================

    print("\n")
    print("=" * 70)
    print("QIP-QC PHASE 1B — FINAL SUMMARY")
    print("=" * 70)

    for name, report in reports.items():

        print(
            f"\n{name}"
        )

        print(
            f"  N = {report['N']}"
        )

        print(
            f"  mean = {report['mean']:.8f}"
        )

        print(
            f"  std = {report['std']:.8f}"
        )

        print(
            f"  short corr = "
            f"{report['short_corr']:.8f}"
        )

        print(
            f"  spectral entropy = "
            f"{report['spectral_entropy']:.8f}"
        )

        print(
            f"  Fibonacci error = "
            f"{report['fibonacci_error']:.8f}"
        )

        print(
            f"  Periodic error = "
            f"{report['periodic_error']:.8f}"
        )

    print("\n")
    print("=" * 70)
    print("IMPORTANT")
    print("=" * 70)

    print(
        """
This analysis does NOT establish that the microtubule
is a quasicrystal.

A positive Fibonacci similarity is only evidence of
quasiperiodic compatibility.

A rigorous claim requires:
    - correct protofilament assignment
    - multiple independent structures
    - long-range correlations
    - diffraction / structure-factor analysis
    - surrogate statistics
    - comparison with periodic and random controls
"""
    )

    print("=" * 70)


if __name__ == "__main__":

    main()

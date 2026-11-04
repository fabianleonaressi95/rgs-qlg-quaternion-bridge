"""
QIP-QC PHASE 1C
=============================================================
3D STRUCTURE FACTOR / DIFFRACTION-LIKE ANALYSIS

Input:
    PDB 5IJ0

Tests:
    REAL MICROTUBULE
    PERIODIC CONTROL
    FIBONACCI CONTROL
    RANDOM CONTROL

Main observable:

        S(k) = | sum_j exp(i k.r_j) |^2 / N

The script does NOT assume that the microtubule is
quasicrystalline.

It asks whether the real structure contains:

    - long-range order
    - discrete reciprocal-space peaks
    - non-periodic but deterministic structure
    - structure-factor characteristics compatible
      with quasiperiodicity

IMPORTANT:
    This is a structural screening test, not a proof
    of consciousness or quantum effects.
"""

import os
import urllib.request
import numpy as np
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from Bio.PDB import PDBParser


# ============================================================
# CONFIGURATION
# ============================================================

PDB_ID = "5IJ0"

DATA_DIR = "qip_qc_phase1c"
FIG_DIR = os.path.join(
    DATA_DIR,
    "figures"
)

os.makedirs(
    DATA_DIR,
    exist_ok=True
)

os.makedirs(
    FIG_DIR,
    exist_ok=True
)

SEED = 42

N_POINTS = 856

N_K = 1200

K_MIN = 0.02

K_MAX = 2.50

N_SURROGATES = 100


# ============================================================
# DOWNLOAD PDB
# ============================================================

def download_pdb():

    filename = os.path.join(
        DATA_DIR,
        f"{PDB_ID}.pdb"
    )

    if os.path.exists(filename):

        print(
            "[OK] PDB already available:",
            filename
        )

        return filename

    url = (
        "https://files.rcsb.org/download/"
        f"{PDB_ID}.pdb"
    )

    print(
        "[DOWNLOAD]",
        url
    )

    urllib.request.urlretrieve(
        url,
        filename
    )

    return filename


# ============================================================
# LOAD PDB
# ============================================================

def load_ca_coordinates(
    pdb_file
):

    parser = PDBParser(
        QUIET=True
    )

    structure = parser.get_structure(
        PDB_ID,
        pdb_file
    )

    coords = []

    model = next(
        structure.get_models()
    )

    for chain in model:

        for residue in chain:

            if residue.id[0] != " ":
                continue

            if "CA" not in residue:
                continue

            coords.append(
                residue["CA"].coord
            )

    coords = np.asarray(
        coords,
        dtype=float
    )

    print(
        "[STRUCTURE] points:",
        len(coords)
    )

    return coords


# ============================================================
# CENTER + NORMALIZE
# ============================================================

def normalize_structure(
    coords
):

    coords = np.asarray(
        coords,
        dtype=float
    )

    center = np.mean(
        coords,
        axis=0
    )

    x = coords - center

    scale = np.sqrt(
        np.mean(
            np.sum(
                x**2,
                axis=1
            )
        )
    )

    x /= scale

    return x


# ============================================================
# PRINCIPAL AXIS
# ============================================================

def principal_axis(
    coords
):

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

    axis /= np.linalg.norm(
        axis
    )

    return axis


# ============================================================
# ALIGN STRUCTURE TO PRINCIPAL AXIS
# ============================================================

def align_to_z(
    coords
):

    axis = principal_axis(
        coords
    )

    z = np.array(
        [0.0, 0.0, 1.0]
    )

    v = np.cross(
        axis,
        z
    )

    s = np.linalg.norm(v)

    c = np.dot(
        axis,
        z
    )

    if s < 1e-12:

        return coords

    vx = np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

    R = (
        np.eye(3)
        + vx
        + vx @ vx
        * ((1 - c) / (s**2))
    )

    return coords @ R.T


# ============================================================
# RESAMPLE POINTS
# ============================================================

def resample_points(
    coords,
    n=N_POINTS,
    seed=42
):

    if len(coords) <= n:

        return coords.copy()

    rng = np.random.default_rng(
        seed
    )

    indices = np.sort(
        rng.choice(
            len(coords),
            size=n,
            replace=False
        )
    )

    return coords[
        indices
    ]


# ============================================================
# STRUCTURE FACTOR — 1D
# ============================================================

def structure_factor_1d(
    coordinates,
    axis,
    k_values
):

    r = coordinates @ axis

    N = len(r)

    S = np.zeros(
        len(k_values)
    )

    for i, k in enumerate(
        k_values
    ):

        phase = np.exp(
            1j * k * r
        )

        amplitude = np.sum(
            phase
        )

        S[i] = (
            np.abs(amplitude)**2
            / N
        )

    return S


# ============================================================
# STRUCTURE FACTOR — 3D POWDER
# ============================================================

def structure_factor_3d_powder(
    coords,
    k_values,
    n_directions=64,
    seed=42
):

    """
    Directionally averaged 3D structure factor.

    S(k) = average_dir
           |sum_j exp(i k u.r_j)|² / N
    """

    rng = np.random.default_rng(
        seed
    )

    N = len(coords)

    directions = rng.normal(
        size=(
            n_directions,
            3
        )
    )

    directions /= np.linalg.norm(
        directions,
        axis=1,
        keepdims=True
    )

    S = np.zeros(
        len(k_values)
    )

    for u in directions:

        projection = (
            coords @ u
        )

        for i, k in enumerate(
            k_values
        ):

            amplitude = np.sum(
                np.exp(
                    1j * k * projection
                )
            )

            S[i] += (
                np.abs(amplitude)**2
                / N
            )

    S /= n_directions

    return S


# ============================================================
# PERIODIC CONTROL
# ============================================================

def generate_periodic(
    n,
    radius=1.0,
    pitch=0.35
):

    theta = np.arange(n) * (
        2 * np.pi / 13
    )

    z = (
        np.arange(n)
        * pitch
    )

    x = radius * np.cos(
        theta
    )

    y = radius * np.sin(
        theta
    )

    return np.column_stack(
        [
            x,
            y,
            z
        ]
    )


# ============================================================
# FIBONACCI CONTROL
# ============================================================

def fibonacci_word(
    n
):

    word = "L"

    while len(word) < n:

        new = ""

        for c in word:

            if c == "L":

                new += "LS"

            else:

                new += "L"

        word = new

    return word[:n]


def generate_fibonacci(
    n,
    radius=1.0,
    long_step=0.35,
    short_step=0.35 / 1.61803398875
):

    word = fibonacci_word(
        n
    )

    z = np.zeros(n)

    for i in range(
        1,
        n
    ):

        if word[i - 1] == "L":

            step = long_step

        else:

            step = short_step

        z[i] = z[i - 1] + step

    # quasiperiodic angular modulation
    phi = (
        (1 + np.sqrt(5))
        / 2
    )

    theta = (
        2 * np.pi
        * np.arange(n)
        / phi
    )

    x = radius * np.cos(
        theta
    )

    y = radius * np.sin(
        theta
    )

    return np.column_stack(
        [
            x,
            y,
            z
        ]
    )


# ============================================================
# RANDOM CONTROL
# ============================================================

def generate_random(
    n,
    radius=1.0,
    pitch=0.35,
    seed=42
):

    rng = np.random.default_rng(
        seed
    )

    theta = rng.uniform(
        0,
        2 * np.pi,
        size=n
    )

    z = np.arange(n) * pitch

    z += rng.normal(
        0,
        pitch * 0.35,
        size=n
    )

    r = radius + rng.normal(
        0,
        radius * 0.15,
        size=n
    )

    x = r * np.cos(
        theta
    )

    y = r * np.sin(
        theta
    )

    return np.column_stack(
        [
            x,
            y,
            z
        ]
    )


# ============================================================
# RANDOM PHASE SURROGATE
# ============================================================

def fourier_phase_surrogate(
    x,
    seed=42
):

    rng = np.random.default_rng(
        seed
    )

    x = np.asarray(
        x
    )

    X = np.fft.rfft(
        x - np.mean(x)
    )

    amp = np.abs(X)

    phase = rng.uniform(
        0,
        2 * np.pi,
        len(X)
    )

    phase[0] = 0

    if len(x) % 2 == 0:

        phase[-1] = 0

    Y = (
        amp
        *
        np.exp(
            1j * phase
        )
    )

    return np.fft.irfft(
        Y,
        n=len(x)
    )


# ============================================================
# PEAK DETECTION
# ============================================================

def detect_peaks(
    k,
    S,
    threshold_sigma=3.0
):

    baseline = np.median(
        S
    )

    mad = np.median(
        np.abs(
            S - baseline
        )
    )

    sigma = (
        1.4826
        * mad
        + 1e-12
    )

    threshold = (
        baseline
        + threshold_sigma * sigma
    )

    peaks = []

    for i in range(
        1,
        len(S) - 1
    ):

        if (
            S[i] > S[i - 1]
            and S[i] > S[i + 1]
            and S[i] > threshold
        ):

            peaks.append(i)

    return np.asarray(
        peaks,
        dtype=int
    )


# ============================================================
# PEAK STATISTICS
# ============================================================

def peak_statistics(
    k,
    S
):

    peaks = detect_peaks(
        k,
        S
    )

    if len(peaks) < 2:

        return {
            "n_peaks": len(peaks),
            "peak_density": 0.0,
            "mean_spacing": np.nan,
            "spacing_cv": np.nan
        }

    kp = k[
        peaks
    ]

    spacings = np.diff(
        kp
    )

    mean_spacing = np.mean(
        spacings
    )

    cv = (
        np.std(spacings)
        /
        (
            mean_spacing
            + 1e-12
        )
    )

    density = (
        len(peaks)
        /
        (
            k[-1] - k[0]
        )
    )

    return {
        "n_peaks": len(peaks),
        "peak_density": density,
        "mean_spacing": mean_spacing,
        "spacing_cv": cv
    }


# ============================================================
# ORDER SCORE
# ============================================================

def order_score(
    S
):

    """
    Ratio between structured peak power and
    diffuse background.

    Higher = stronger long-range order.
    """

    median = np.median(
        S
    )

    p95 = np.percentile(
        S,
        95
    )

    p99 = np.percentile(
        S,
        99
    )

    return (
        p99 - median
    ) / (
        p95 - median
        + 1e-12
    )


# ============================================================
# SPECTRAL ENTROPY
# ============================================================

def spectral_entropy(
    S
):

    p = np.asarray(
        S,
        dtype=float
    )

    p = p / (
        np.sum(p)
        + 1e-15
    )

    H = -np.sum(
        p * np.log(
            p + 1e-15
        )
    )

    return H / np.log(
        len(p)
    )


# ============================================================
# COMPARE STRUCTURES
# ============================================================

def compare_spectra(
    S_real,
    S_control
):

    a = S_real / (
        np.sum(S_real)
        + 1e-15
    )

    b = S_control / (
        np.sum(S_control)
        + 1e-15
    )

    return np.sqrt(
        np.mean(
            (a - b)**2
        )
    )


# ============================================================
# RADIAL DISTRIBUTION
# ============================================================

def radial_distribution(
    coords,
    bins=120
):

    r = np.linalg.norm(
        coords[:, :2],
        axis=1
    )

    hist, edges = np.histogram(
        r,
        bins=bins,
        density=True
    )

    centers = (
        edges[:-1]
        + edges[1:]
    ) / 2

    return centers, hist


# ============================================================
# PLOT STRUCTURE FACTORS
# ============================================================

def plot_structure_factors(
    k,
    spectra,
    filename
):

    plt.figure(
        figsize=(13, 7)
    )

    for name, S in spectra.items():

        S_norm = S / (
            np.max(S)
            + 1e-12
        )

        plt.plot(
            k,
            S_norm,
            label=name,
            linewidth=1
        )

    plt.xlabel(
        "Wave number k"
    )

    plt.ylabel(
        "Normalized S(k)"
    )

    plt.title(
        "QIP-QC Phase 1C — Structure Factor"
    )

    plt.legend()

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
# PLOT LOG STRUCTURE FACTOR
# ============================================================

def plot_log_structure_factors(
    k,
    spectra
):

    plt.figure(
        figsize=(13, 7)
    )

    for name, S in spectra.items():

        plt.semilogy(
            k,
            S + 1e-12,
            label=name
        )

    plt.xlabel(
        "Wave number k"
    )

    plt.ylabel(
        "S(k)"
    )

    plt.title(
        "QIP-QC Phase 1C — Log Structure Factor"
    )

    plt.legend()

    plt.tight_layout()

    path = os.path.join(
        FIG_DIR,
        "structure_factor_log.png"
    )

    plt.savefig(
        path,
        dpi=200
    )

    plt.close()


# ============================================================
# MAIN
# ============================================================

def main():

    print("=" * 75)

    print(
        "QIP-QC PHASE 1C"
    )

    print(
        "3D STRUCTURE FACTOR / DIFFRACTION ANALYSIS"
    )

    print("=" * 75)

    # --------------------------------------------------------
    # Load real structure
    # --------------------------------------------------------

    pdb_file = download_pdb()

    coords = load_ca_coordinates(
        pdb_file
    )

    coords = normalize_structure(
        coords
    )

    coords = align_to_z(
        coords
    )

    coords = resample_points(
        coords,
        N_POINTS,
        SEED
    )

    N = len(coords)

    print(
        "[REAL STRUCTURE]",
        N,
        "points"
    )

    # --------------------------------------------------------
    # Controls
    # --------------------------------------------------------

    periodic = generate_periodic(
        N
    )

    fibonacci = generate_fibonacci(
        N
    )

    random_structure = generate_random(
        N,
        seed=SEED
    )

    periodic = normalize_structure(
        periodic
    )

    fibonacci = normalize_structure(
        fibonacci
    )

    random_structure = normalize_structure(
        random_structure
    )

    # --------------------------------------------------------
    # k grid
    # --------------------------------------------------------

    k = np.linspace(
        K_MIN,
        K_MAX,
        N_K
    )

    # --------------------------------------------------------
    # 3D powder spectra
    # --------------------------------------------------------

    print(
        "\nCalculating REAL structure factor..."
    )

    S_real = structure_factor_3d_powder(
        coords,
        k,
        n_directions=96,
        seed=SEED
    )

    print(
        "Calculating PERIODIC control..."
    )

    S_periodic = structure_factor_3d_powder(
        periodic,
        k,
        n_directions=96,
        seed=SEED
    )

    print(
        "Calculating FIBONACCI control..."
    )

    S_fibonacci = structure_factor_3d_powder(
        fibonacci,
        k,
        n_directions=96,
        seed=SEED
    )

    print(
        "Calculating RANDOM control..."
    )

    S_random = structure_factor_3d_powder(
        random_structure,
        k,
        n_directions=96,
        seed=SEED
    )

    spectra = {

        "Microtubule": S_real,

        "Periodic": S_periodic,

        "Fibonacci": S_fibonacci,

        "Random": S_random
    }

    # --------------------------------------------------------
    # Statistics
    # --------------------------------------------------------

    print("\n")
    print("=" * 75)
    print("STRUCTURE FACTOR STATISTICS")
    print("=" * 75)

    stats = {}

    for name, S in spectra.items():

        peaks = peak_statistics(
            k,
            S
        )

        entropy = spectral_entropy(
            S
        )

        order = order_score(
            S
        )

        stats[name] = {
            **peaks,
            "entropy": entropy,
            "order_score": order
        }

        print(
            f"\n{name}"
        )

        print(
            f"  peaks          = "
            f"{peaks['n_peaks']}"
        )

        print(
            f"  peak density   = "
            f"{peaks['peak_density']:.8f}"
        )

        print(
            f"  mean spacing   = "
            f"{peaks['mean_spacing']:.8f}"
        )

        print(
            f"  spacing CV     = "
            f"{peaks['spacing_cv']:.8f}"
        )

        print(
            f"  entropy        = "
            f"{entropy:.8f}"
        )

        print(
            f"  order score    = "
            f"{order:.8f}"
        )

    # --------------------------------------------------------
    # Spectral distances
    # --------------------------------------------------------

    print("\n")
    print("=" * 75)
    print("SPECTRAL DISTANCES")
    print("=" * 75)

    d_periodic = compare_spectra(
        S_real,
        S_periodic
    )

    d_fibonacci = compare_spectra(
        S_real,
        S_fibonacci
    )

    d_random = compare_spectra(
        S_real,
        S_random
    )

    print(
        f"Microtubule vs Periodic  = "
        f"{d_periodic:.10f}"
    )

    print(
        f"Microtubule vs Fibonacci = "
        f"{d_fibonacci:.10f}"
    )

    print(
        f"Microtubule vs Random    = "
        f"{d_random:.10f}"
    )

    # --------------------------------------------------------
    # Plot
    # --------------------------------------------------------

    plot_structure_factors(
        k,
        spectra,
        "structure_factor.png"
    )

    plot_log_structure_factors(
        k,
        spectra
    )

    # --------------------------------------------------------
    # Save
    # --------------------------------------------------------

    np.savez(
        os.path.join(
            DATA_DIR,
            "phase1c_structure_factors.npz"
        ),
        k=k,
        S_real=S_real,
        S_periodic=S_periodic,
        S_fibonacci=S_fibonacci,
        S_random=S_random
    )

    # --------------------------------------------------------
    # 1D axial structure factor
    # --------------------------------------------------------

    axis = np.array(
        [0.0, 0.0, 1.0]
    )

    k1d = np.linspace(
        K_MIN,
        K_MAX,
        N_K
    )

    print(
        "\nCalculating axial structure factor..."
    )

    S1d_real = structure_factor_1d(
        coords,
        axis,
        k1d
    )

    S1d_periodic = structure_factor_1d(
        periodic,
        axis,
        k1d
    )

    S1d_fibonacci = structure_factor_1d(
        fibonacci,
        axis,
        k1d
    )

    S1d_random = structure_factor_1d(
        random_structure,
        axis,
        k1d
    )

    plt.figure(
        figsize=(13, 7)
    )

    plt.plot(
        k1d,
        S1d_real /
        np.max(S1d_real),
        label="Microtubule"
    )

    plt.plot(
        k1d,
        S1d_periodic /
        np.max(S1d_periodic),
        label="Periodic"
    )

    plt.plot(
        k1d,
        S1d_fibonacci /
        np.max(S1d_fibonacci),
        label="Fibonacci"
    )

    plt.plot(
        k1d,
        S1d_random /
        np.max(S1d_random),
        label="Random"
    )

    plt.xlabel(
        "Axial wave number k"
    )

    plt.ylabel(
        "Normalized S(k)"
    )

    plt.title(
        "Axial Structure Factor"
    )

    plt.legend()

    plt.tight_layout()

    plt.savefig(
        os.path.join(
            FIG_DIR,
            "axial_structure_factor.png"
        ),
        dpi=200
    )

    plt.close()

    # --------------------------------------------------------
    # Final interpretation
    # --------------------------------------------------------

    print("\n")
    print("=" * 75)
    print("PHASE 1C INTERPRETATION")
    print("=" * 75)

    print(
        """
IMPORTANT:

A large order_score does NOT automatically imply
quasicrystallinity.

For a quasicrystalline interpretation we need:

    1. Strong non-random long-range order
    2. Discrete reciprocal-space structure
    3. Absence of a single conventional period
    4. Reproducibility across independent structures
    5. Robustness against resolution/noise
    6. Comparison with known quasiperiodic controls

The strongest possible result is:

    Random << Microtubule ~ Fibonacci
                       and
    Microtubule != Periodic

But even that is only evidence of
quasiperiodic compatibility.

It is NOT evidence of consciousness.
"""
    )

    print("=" * 75)

    print(
        "\nData saved in:",
        DATA_DIR
    )


# ============================================================
# EXECUTE
# ============================================================

if __name__ == "__main__":

    main()

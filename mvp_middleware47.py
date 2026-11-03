"""
QIP-QC PHASE 1D — VECTORIZED
============================================================

Complete cryo-EM microtubule
EMDB-25156
PDB 7SJ7

Pipeline:

    cryo-EM map
        |
        v
    density coordinates
        |
        v
    PCA axis alignment
        |
        v
    cylindrical coordinates
        |
        +----------------------+
        |                      |
        v                      v
    angular spectrum       axial spectrum
        |                      |
        +----------+-----------+
                   |
                   v
             S(m, kz)
                   |
                   v
       periodic / quasiperiodic
             structure test

Optimizations:
    - vectorized NumPy
    - matrix multiplication
    - no Python loop over voxels
    - block processing for S(m,kz)
    - float32 where appropriate
"""

import os
import urllib.request

import numpy as np
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt


# ============================================================
# CONFIGURATION
# ============================================================

EMDB_ID = "EMD-25156"
PDB_ID = "7SJ7"

OUTDIR = "QIP_QC_EMD25156"

FIGDIR = os.path.join(
    OUTDIR,
    "figures"
)

os.makedirs(
    FIGDIR,
    exist_ok=True
)


MAP_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/emdb/"
    "structures/EMD-25156/map/"
    "emd_25156.map.gz"
)

PDB_URL = (
    "https://files.rcsb.org/download/"
    f"{PDB_ID}.pdb"
)

MAP_FILE = os.path.join(
    OUTDIR,
    "emd_25156.map.gz"
)

PDB_FILE = os.path.join(
    OUTDIR,
    f"{PDB_ID}.pdb"
)


# ============================================================
# DOWNLOAD
# ============================================================

def download_file(url, filename):

    if os.path.exists(filename):

        print("[OK] Existing:", filename)

        return

    print("[DOWNLOAD]", url)

    urllib.request.urlretrieve(
        url,
        filename
    )

    print("[OK]", filename)


def download_dataset():

    print("=" * 80)
    print("DOWNLOADING DATASET")
    print("=" * 80)

    download_file(
        MAP_URL,
        MAP_FILE
    )

    download_file(
        PDB_URL,
        PDB_FILE
    )


# ============================================================
# LOAD MAP
# ============================================================

def load_em_map():

    try:
        import mrcfile
    except ImportError:

        raise RuntimeError(
            "Install mrcfile:\n"
            "pip install mrcfile"
        )

    print("\nLoading cryo-EM map...")

    with mrcfile.open(
        MAP_FILE,
        permissive=True
    ) as mrc:

        density = np.asarray(
            mrc.data,
            dtype=np.float32
        )

        voxel = np.asarray(
            mrc.voxel_size.tolist(),
            dtype=np.float64
        )

        origin = np.asarray(
            mrc.header.origin.tolist(),
            dtype=np.float64
        )

    print("Map shape:", density.shape)
    print("Voxel size:", voxel)
    print("Density range:",
          float(density.min()),
          float(density.max()))

    return density, voxel, origin


# ============================================================
# THRESHOLD
# ============================================================

def threshold_density(
    density,
    percentile=95
):

    threshold = np.percentile(
        density,
        percentile
    )

    mask = density >= threshold

    print(
        "\nThreshold:",
        threshold
    )

    print(
        "Occupied voxels:",
        int(mask.sum())
    )

    return mask, threshold


# ============================================================
# VECTORIZED VOXEL COORDINATES
# ============================================================

def density_coordinates(
    mask,
    voxel
):

    print(
        "\nGenerating coordinates..."
    )

    indices = np.argwhere(
        mask
    )

    # MRC convention:
    # indices = z, y, x

    xyz = (
        indices[:, ::-1]
        *
        voxel
    )

    xyz = (
        xyz
        -
        xyz.mean(
            axis=0,
            keepdims=True
        )
    )

    return xyz.astype(
        np.float32
    )


# ============================================================
# PCA AXIS
# ============================================================

def principal_axis(
    coords
):

    centered = (
        coords
        -
        coords.mean(
            axis=0,
            keepdims=True
        )
    )

    covariance = (
        centered.T
        @
        centered
    )

    covariance /= (
        len(coords) - 1
    )

    eigenvalues, eigenvectors = (
        np.linalg.eigh(
            covariance
        )
    )

    axis = eigenvectors[
        :,
        np.argmax(eigenvalues)
    ]

    axis /= np.linalg.norm(
        axis
    )

    return axis.astype(
        np.float32
    )


# ============================================================
# ROTATION MATRIX
# ============================================================

def rotation_matrix(
    a,
    b
):

    a = (
        a /
        np.linalg.norm(a)
    )

    b = (
        b /
        np.linalg.norm(b)
    )

    v = np.cross(
        a,
        b
    )

    c = np.dot(
        a,
        b
    )

    s = np.linalg.norm(v)

    if s < 1e-12:

        if c > 0:
            return np.eye(
                3,
                dtype=np.float32
            )

        return np.diag(
            [
                1.,
                -1.,
                -1.
            ]
        ).astype(
            np.float32
        )

    K = np.array(
        [
            [0., -v[2], v[1]],
            [v[2], 0., -v[0]],
            [-v[1], v[0], 0.]
        ],
        dtype=np.float32
    )

    return (
        np.eye(
            3,
            dtype=np.float32
        )
        +
        K
        +
        (
            K @ K
            *
            (
                (1. - c)
                /
                (s * s)
            )
        )
    )


# ============================================================
# ALIGN AXIS
# ============================================================

def align_axis(
    coords
):

    axis = principal_axis(
        coords
    )

    target = np.array(
        [
            0.,
            0.,
            1.
        ],
        dtype=np.float32
    )

    R = rotation_matrix(
        axis,
        target
    )

    # IMPORTANT:
    # one matrix multiplication instead
    # of point-by-point rotation

    aligned = (
        coords
        @
        R.T
    )

    return aligned.astype(
        np.float32
    )


# ============================================================
# CYLINDRICAL COORDINATES
# ============================================================

def cylindrical_coordinates(
    coords
):

    x = coords[:, 0]
    y = coords[:, 1]
    z = coords[:, 2]

    radius = np.hypot(
        x,
        y
    ).astype(
        np.float32
    )

    theta = np.arctan2(
        y,
        x
    ).astype(
        np.float32
    )

    return theta, radius, z


# ============================================================
# VECTORIZED ANGULAR STRUCTURE FACTOR
# ============================================================

def angular_structure_factor(
    theta,
    m_values
):

    """
    Original:

        for m:
            for theta:
                exp(i*m*theta)

    New:

        exp(i * outer(m, theta))

    Then sum over theta.

    Complexity remains O(MN), but the
    inner operation is executed by optimized
    NumPy kernels instead of Python.
    """

    theta64 = theta.astype(
        np.float64,
        copy=False
    )

    m64 = m_values.astype(
        np.float64,
        copy=False
    )

    phase = np.exp(
        1j *
        np.outer(
            m64,
            theta64
        )
    )

    amplitude = phase.sum(
        axis=1
    )

    N = theta.size

    S = (
        np.abs(amplitude) ** 2
        /
        N
    )

    return S.astype(
        np.float64
    )


# ============================================================
# VECTORIZED AXIAL STRUCTURE FACTOR
# ============================================================

def axial_structure_factor(
    z,
    kz_values,
    block_size=64
):

    """
    Block-vectorized implementation.

    We avoid creating a gigantic
    len(kz) x len(z) matrix.
    """

    z = z.astype(
        np.float64,
        copy=False
    )

    kz_values = kz_values.astype(
        np.float64,
        copy=False
    )

    N = z.size

    S = np.empty(
        len(kz_values),
        dtype=np.float64
    )

    for start in range(
        0,
        len(kz_values),
        block_size
    ):

        stop = min(
            start + block_size,
            len(kz_values)
        )

        kz_block = (
            kz_values[start:stop]
        )

        phase = np.exp(
            1j
            *
            (
                kz_block[:, None]
                *
                z[None, :]
            )
        )

        amplitude = phase.sum(
            axis=1
        )

        S[start:stop] = (
            np.abs(amplitude) ** 2
            /
            N
        )

    return S


# ============================================================
# VECTORIZED 2D CYLINDRICAL SPECTRUM
# ============================================================

def cylindrical_structure_factor(
    theta,
    z,
    m_values,
    kz_values,
    block_m=8,
    block_kz=32
):

    """
    Calculates

        S(m,kz)
        =
        | Σ_j exp(i[m theta_j + kz z_j])|² / N

    without Python loops over particles.

    Blocks prevent RAM explosion.
    """

    theta = theta.astype(
        np.float64,
        copy=False
    )

    z = z.astype(
        np.float64,
        copy=False
    )

    m_values = m_values.astype(
        np.float64,
        copy=False
    )

    kz_values = kz_values.astype(
        np.float64,
        copy=False
    )

    N = theta.size

    S = np.empty(
        (
            len(m_values),
            len(kz_values)
        ),
        dtype=np.float64
    )

    # --------------------------------------------------------
    # Block over angular harmonics
    # and axial frequencies
    # --------------------------------------------------------

    for m_start in range(
        0,
        len(m_values),
        block_m
    ):

        m_stop = min(
            m_start + block_m,
            len(m_values)
        )

        m_block = (
            m_values[m_start:m_stop]
        )

        angular_phase = np.exp(
            1j
            *
            (
                m_block[:, None]
                *
                theta[None, :]
            )
        )

        for k_start in range(
            0,
            len(kz_values),
            block_kz
        ):

            k_stop = min(
                k_start + block_kz,
                len(kz_values)
            )

            k_block = (
                kz_values[k_start:k_stop]
            )

            axial_phase = np.exp(
                1j
                *
                (
                    k_block[:, None]
                    *
                    z[None, :]
                )
            )

            # ------------------------------------------------
            # Matrix multiplication
            #
            # angular_phase:
            #       M x N
            #
            # conjugate(axial_phase):
            #       K x N
            #
            # product:
            #       M x K
            # ------------------------------------------------

            amplitude = (
                angular_phase
                @
                axial_phase.T
            )

            S[
                m_start:m_stop,
                k_start:k_stop
            ] = (
                np.abs(amplitude) ** 2
                /
                N
            )

    return S


# ============================================================
# SPECTRAL ENTROPY
# ============================================================

def spectral_entropy(
    spectrum
):

    spectrum = np.asarray(
        spectrum,
        dtype=np.float64
    )

    total = spectrum.sum()

    if total <= 0:
        return np.nan

    p = spectrum / total

    p = p[p > 0]

    H = -np.sum(
        p * np.log(p)
    )

    Hmax = np.log(
        len(spectrum)
    )

    return H / Hmax


# ============================================================
# SHORT-RANGE CORRELATION
# ============================================================

def short_range_correlation(
    x,
    lag=1
):

    x = np.asarray(
        x,
        dtype=np.float64
    )

    if len(x) <= lag:
        return np.nan

    a = x[:-lag]
    b = x[lag:]

    a -= a.mean()
    b -= b.mean()

    denom = (
        np.linalg.norm(a)
        *
        np.linalg.norm(b)
    )

    if denom == 0:
        return np.nan

    return np.dot(
        a,
        b
    ) / denom


# ============================================================
# RADIAL PROFILE
# ============================================================

def radial_profile(
    radius,
    z,
    bins=856
):

    zmin = z.min()
    zmax = z.max()

    edges = np.linspace(
        zmin,
        zmax,
        bins + 1
    )

    idx = np.digitize(
        z,
        edges
    ) - 1

    valid = (
        (idx >= 0)
        &
        (idx < bins)
    )

    sums = np.bincount(
        idx[valid],
        weights=radius[valid],
        minlength=bins
    )

    counts = np.bincount(
        idx[valid],
        minlength=bins
    )

    profile = (
        sums
        /
        np.maximum(
            counts,
            1
        )
    )

    centers = (
        edges[:-1]
        +
        edges[1:]
    ) / 2

    return centers, profile


# ============================================================
# ANGULAR PROFILE
# ============================================================

def angular_profile(
    theta,
    bins=360
):

    edges = np.linspace(
        -np.pi,
        np.pi,
        bins + 1
    )

    idx = np.digitize(
        theta,
        edges
    ) - 1

    valid = (
        (idx >= 0)
        &
        (idx < bins)
    )

    counts = np.bincount(
        idx[valid],
        minlength=bins
    )

    centers = (
        edges[:-1]
        +
        edges[1:]
    ) / 2

    return centers, counts


# ============================================================
# PLOTS
# ============================================================

def plot_1d(
    x,
    y,
    xlabel,
    ylabel,
    title,
    filename
):

    plt.figure(
        figsize=(12, 6)
    )

    plt.plot(
        x,
        y
    )

    plt.xlabel(
        xlabel
    )

    plt.ylabel(
        ylabel
    )

    plt.title(
        title
    )

    plt.tight_layout()

    plt.savefig(
        os.path.join(
            FIGDIR,
            filename
        ),
        dpi=200
    )

    plt.close()


def plot_2d(
    m,
    kz,
    S
):

    plt.figure(
        figsize=(14, 8)
    )

    plt.imshow(
        np.log10(
            S + 1e-12
        ),
        origin="lower",
        aspect="auto",
        extent=[
            kz[0],
            kz[-1],
            m[0],
            m[-1]
        ]
    )

    plt.xlabel(
        "kz (1/Å)"
    )

    plt.ylabel(
        "Angular harmonic m"
    )

    plt.title(
        "QIP-QC — S(m,kz)"
    )

    plt.colorbar(
        label="log10 S"
    )

    plt.tight_layout()

    plt.savefig(
        os.path.join(
            FIGDIR,
            "S_m_kz.png"
        ),
        dpi=200
    )

    plt.close()


# ============================================================
# MAIN
# ============================================================

def main():

    print(
        "\n"
        + "=" * 80
    )

    print(
        "QIP-QC PHASE 1D — VECTORIZED"
    )

    print(
        "COMPLETE CRYO-EM MICROTUBULE"
    )

    print(
        f"{EMDB_ID} / {PDB_ID}"
    )

    print(
        "=" * 80
    )

    # --------------------------------------------------------
    # DOWNLOAD
    # --------------------------------------------------------

    download_dataset()

    # --------------------------------------------------------
    # MAP
    # --------------------------------------------------------

    density, voxel, origin = (
        load_em_map()
    )

    # --------------------------------------------------------
    # MASK
    # --------------------------------------------------------

    mask, threshold = (
        threshold_density(
            density,
            percentile=95
        )
    )

    # --------------------------------------------------------
    # COORDINATES
    # --------------------------------------------------------

    coords = density_coordinates(
        mask,
        voxel
    )

    print(
        "Coordinates:",
        coords.shape
    )

    # --------------------------------------------------------
    # ALIGN
    # --------------------------------------------------------

    print(
        "\nAligning principal axis..."
    )

    coords = align_axis(
        coords
    )

    # --------------------------------------------------------
    # CYLINDRICAL
    # --------------------------------------------------------

    theta, radius, z = (
        cylindrical_coordinates(
            coords
        )
    )

    print(
        "\n"
        + "=" * 70
    )

    print(
        "STRUCTURAL GEOMETRY"
    )

    print(
        "=" * 70
    )

    print(
        "N:",
        len(z)
    )

    print(
        "Radius mean:",
        np.mean(radius)
    )

    print(
        "Radius std:",
        np.std(radius)
    )

    print(
        "Z range:",
        z.min(),
        z.max()
    )

    print(
        "Angular short correlation:",
        short_range_correlation(
            angular_profile(theta)[1]
        )
    )

    # --------------------------------------------------------
    # SPECTRA
    # --------------------------------------------------------

    m = np.arange(
        0,
        101
    )

    kz = np.linspace(
        0.01,
        1.5,
        400
    )

    print(
        "\nCalculating angular spectrum..."
    )

    S_m = angular_structure_factor(
        theta,
        m
    )

    print(
        "Calculating axial spectrum..."
    )

    S_kz = axial_structure_factor(
        z,
        kz,
        block_size=64
    )

    # --------------------------------------------------------
    # 2D
    # --------------------------------------------------------

    m2 = np.arange(
        0,
        51
    )

    kz2 = np.linspace(
        0.01,
        1.5,
        250
    )

    print(
        "Calculating S(m,kz)..."
    )

    S2 = cylindrical_structure_factor(
        theta,
        z,
        m2,
        kz2,
        block_m=8,
        block_kz=32
    )

    # --------------------------------------------------------
    # ENTROPY
    # --------------------------------------------------------

    Hm = spectral_entropy(
        S_m
    )

    Hk = spectral_entropy(
        S_kz
    )

    H2 = spectral_entropy(
        S2.ravel()
    )

    # --------------------------------------------------------
    # DOMINANT HARMONICS
    # --------------------------------------------------------

    order = np.argsort(
        S_m[1:]
    )[::-1]

    dominant = (
        order[:20]
        + 1
    )

    # --------------------------------------------------------
    # RADIAL PROFILE
    # --------------------------------------------------------

    zr, radial = radial_profile(
        radius,
        z
    )

    angular_x, angular = (
        angular_profile(theta)
    )

    # --------------------------------------------------------
    # RESULTS
    # --------------------------------------------------------

    print(
        "\n"
        + "=" * 80
    )

    print(
        "QIP-QC PHASE 1 — STRUCTURAL SUMMARY"
    )

    print(
        "=" * 80
    )

    print(
        "Dataset:",
        EMDB_ID,
        "/",
        PDB_ID
    )

    print(
        "Real structural signal length:",
        len(z)
    )

    print(
        "Radius mean:",
        f"{np.mean(radius):.8f}"
    )

    print(
        "Radius std:",
        f"{np.std(radius):.8f}"
    )

    print(
        "Angular short-range corr:",
        f"{short_range_correlation(angular):.8f}"
    )

    print(
        "Angular spectral entropy:",
        f"{Hm:.8f}"
    )

    print(
        "Axial spectral entropy:",
        f"{Hk:.8f}"
    )

    print(
        "2D spectral entropy:",
        f"{H2:.8f}"
    )

    print(
        "\nDominant angular harmonics:"
    )

    for mm in dominant:

        print(
            f"m={mm:3d} "
            f"S={S_m[mm]:.8e}"
        )

    # --------------------------------------------------------
    # SAVE
    # --------------------------------------------------------

    np.savez_compressed(
        os.path.join(
            OUTDIR,
            "QIP_QC_EMD25156_vectorized.npz"
        ),
        coordinates=coords,
        theta=theta,
        radius=radius,
        z=z,
        m=m,
        S_m=S_m,
        kz=kz,
        S_kz=S_kz,
        m2=m2,
        kz2=kz2,
        S2=S2,
        radial_z=zr,
        radial_profile=radial,
        angular_theta=angular_x,
        angular_profile=angular,
        threshold=threshold,
        voxel=voxel
    )

    # --------------------------------------------------------
    # PLOTS
    # --------------------------------------------------------

    plot_1d(
        kz,
        S_kz,
        "kz (1/Å)",
        "S(kz)",
        "Microtubule Axial Structure Factor",
        "axial_structure.png"
    )

    plot_1d(
        m,
        S_m,
        "Angular harmonic m",
        "S(m)",
        "Microtubule Angular Structure Factor",
        "angular_structure.png"
    )

    plot_1d(
        zr,
        radial,
        "z (Å)",
        "Mean radius (Å)",
        "Microtubule Radial Profile",
        "radial_profile.png"
    )

    plot_1d(
        angular_x,
        angular,
        "theta",
        "Density count",
        "Microtubule Angular Density",
        "angular_density.png"
    )

    plot_2d(
        m2,
        kz2,
        S2
    )

    print(
        "\n"
        + "=" * 80
    )

    print(
        "COMPLETE"
    )

    print(
        "Results:",
        OUTDIR
    )

    print(
        "=" * 80
    )


if __name__ == "__main__":
    main()

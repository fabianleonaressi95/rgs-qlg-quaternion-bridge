"""
QIP-QC PHASE 1D — ULTRA-FAST VECTORIZED (FFT Optimized)
============================================================

Complete cryo-EM microtubule
EMDB-25156
PDB 7SJ7
"""

import os
import urllib.request
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Prova a usare scipy.fft per performance multithread superiori
try:
    from scipy.fft import rfft, fft2, fftshift
except ImportError:
    from numpy.fft import rfft, fft2, fftshift


# ============================================================
# CONFIGURATION
# ============================================================

EMDB_ID = "EMD-25156"
PDB_ID = "7SJ7"

OUTDIR = "QIP_QC_EMD25156"
FIGDIR = os.path.join(OUTDIR, "figures")

os.makedirs(FIGDIR, exist_ok=True)

MAP_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/emdb/"
    "structures/EMD-25156/map/"
    "emd_25156.map.gz"
)
PDB_URL = f"https://files.rcsb.org/download/{PDB_ID}.pdb"

MAP_FILE = os.path.join(OUTDIR, "emd_25156.map.gz")
PDB_FILE = os.path.join(OUTDIR, f"{PDB_ID}.pdb")


# ============================================================
# DOWNLOAD & MAP LOADING
# ============================================================

def download_file(url, filename):
    if os.path.exists(filename):
        print("[OK] Existing:", filename)
        return
    print("[DOWNLOAD]", url)
    urllib.request.urlretrieve(url, filename)
    print("[OK]", filename)

def download_dataset():
    print("=" * 80)
    print("DOWNLOADING DATASET")
    print("=" * 80)
    download_file(MAP_URL, MAP_FILE)
    download_file(PDB_URL, PDB_FILE)

def load_em_map():
    try:
        import mrcfile
    except ImportError:
        raise RuntimeError("Install mrcfile:\npip install mrcfile")

    print("\nLoading cryo-EM map...")
    with mrcfile.open(MAP_FILE, permissive=True) as mrc:
        density = np.asarray(mrc.data, dtype=np.float32)
        voxel = np.asarray(mrc.voxel_size.tolist(), dtype=np.float64)
        origin = np.asarray(mrc.header.origin.tolist(), dtype=np.float64)

    print("Map shape:", density.shape)
    print("Voxel size:", voxel)
    print("Density range:", float(density.min()), float(density.max()))
    return density, voxel, origin

def threshold_density(density, percentile=95):
    threshold = np.percentile(density, percentile)
    mask = density >= threshold
    print("\nThreshold:", threshold)
    print("Occupied voxels:", int(mask.sum()))
    return mask, threshold


# ============================================================
# GEOMETRY & ALIGNMENT
# ============================================================

def density_coordinates(mask, voxel):
    print("\nGenerating coordinates...")
    indices = np.argwhere(mask)
    # MRC convention: indices = z, y, x
    xyz = indices[:, ::-1] * voxel
    xyz = xyz - xyz.mean(axis=0, keepdims=True)
    return xyz.astype(np.float32)

def principal_axis(coords):
    centered = coords - coords.mean(axis=0, keepdims=True)
    covariance = (centered.T @ centered) / (len(coords) - 1)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    axis = eigenvectors[:, np.argmax(eigenvalues)]
    return (axis / np.linalg.norm(axis)).astype(np.float32)

def rotation_matrix(a, b):
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)

    if s < 1e-12:
        if c > 0:
            return np.eye(3, dtype=np.float32)
        return np.diag([1.0, -1.0, -1.0]).astype(np.float32)

    K = np.array([
        [0.0, -v[2], v[1]],
        [v[2], 0.0, -v[0]],
        [-v[1], v[0], 0.0]
    ], dtype=np.float32)

    return np.eye(3, dtype=np.float32) + K + (K @ K) * ((1.0 - c) / (s * s))

def align_axis(coords):
    axis = principal_axis(coords)
    target = np.array([0.0, 0.0, 1.0], dtype=np.float32)
    R = rotation_matrix(axis, target)
    return (coords @ R.T).astype(np.float32)

def cylindrical_coordinates(coords):
    x, y, z = coords[:, 0], coords[:, 1], coords[:, 2]
    radius = np.hypot(x, y).astype(np.float32)
    theta = np.arctan2(y, x).astype(np.float32)
    return theta, radius, z


# ============================================================
# HIGH-SPEED SPECTRAL CALCULATIONS (GRID + FFT)
# ============================================================

def angular_structure_factor_fft(theta, m_max=100, n_bins=2048):
    """
    S(m) tramite FFT 1D su un istogramma angolare.
    Velocità: ~100x rispetto all'exp outer product.
    """
    counts, _ = np.histogram(theta, bins=n_bins, range=(-np.pi, np.pi))
    
    # FFT Reale
    fft_vals = rfft(counts)
    N = theta.size
    S = (np.abs(fft_vals) ** 2) / N
    
    m_values = np.arange(0, min(m_max + 1, len(S)))
    return m_values, S[:len(m_values)]


def axial_structure_factor_fft(z, kz_target, n_bins=4096):
    """
    S(kz) interpolato tramite FFT 1D uniforme lungo l'asse Z.
    """
    z_min, z_max = z.min(), z.max()
    counts, _ = np.histogram(z, bins=n_bins, range=(z_min, z_max))
    
    fft_vals = rfft(counts)
    N = z.size
    S_full = (np.abs(fft_vals) ** 2) / N
    
    # Frequenze spatiali associate alla FFT
    L = z_max - z_min
    kz_fft = 2 * np.pi * np.fft.rfftfreq(n_bins, d=L / n_bins)
    
    # Interpolazione sui punti kz_target richiesti
    S_kz = np.interp(kz_target, kz_fft, S_full)
    return S_kz


def cylindrical_structure_factor_fft(theta, z, m_max=50, n_kz_bins=512, n_theta_bins=512):
    """
    S(m, kz) tramite FFT 2D su griglia (theta, z).
    Riduce il calcolo da minuti a frazioni di secondo.
    """
    z_min, z_max = z.min(), z.max()
    
    # Istogramma 2D per mappare la densità voxel su griglia cilindrica
    grid, _, _ = np.histogram2d(
        theta, z,
        bins=[n_theta_bins, n_kz_bins],
        range=[[-np.pi, np.pi], [z_min, z_max]]
    )
    
    # FFT 2D
    fft2d = fft2(grid)
    N = theta.size
    S2_full = (np.abs(fft2d) ** 2) / N
    
    # Estrazione delle frequenze m e kz valide
    m_arr = np.arange(0, m_max + 1)
    
    L = z_max - z_min
    kz_fft = 2 * np.pi * np.fft.fftfreq(n_kz_bins, d=L / n_kz_bins)
    
    # Prendiamo le sole frequenze kz positive
    pos_kz_idx = np.where((kz_fft >= 0.01) & (kz_fft <= 1.5))[0]
    sorted_idx = pos_kz_idx[np.argsort(kz_fft[pos_kz_idx])]
    
    kz_arr = kz_fft[sorted_idx]
    S2 = S2_full[m_arr[:, None], sorted_idx[None, :]]
    
    return m_arr, kz_arr, S2


# ============================================================
# METRICS & PROFILES
# ============================================================

def spectral_entropy(spectrum):
    spectrum = np.asarray(spectrum, dtype=np.float64)
    total = spectrum.sum()
    if total <= 0:
        return np.nan
    p = spectrum / total
    p = p[p > 0]
    return -np.sum(p * np.log(p)) / np.log(len(spectrum))

def short_range_correlation(x, lag=1):
    x = np.asarray(x, dtype=np.float64)
    if len(x) <= lag:
        return np.nan
    a, b = x[:-lag] - x[:-lag].mean(), x[lag:] - x[lag:].mean()
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    return np.dot(a, b) / denom if denom != 0 else np.nan

def radial_profile(radius, z, bins=856):
    edges = np.linspace(z.min(), z.max(), bins + 1)
    idx = np.digitize(z, edges) - 1
    valid = (idx >= 0) & (idx < bins)
    
    sums = np.bincount(idx[valid], weights=radius[valid], minlength=bins)
    counts = np.bincount(idx[valid], minlength=bins)
    
    profile = sums / np.maximum(counts, 1)
    centers = (edges[:-1] + edges[1:]) / 2
    return centers, profile

def angular_profile(theta, bins=360):
    edges = np.linspace(-np.pi, np.pi, bins + 1)
    idx = np.digitize(theta, edges) - 1
    valid = (idx >= 0) & (idx < bins)
    counts = np.bincount(idx[valid], minlength=bins)
    centers = (edges[:-1] + edges[1:]) / 2
    return centers, counts


# ============================================================
# PLOTS
# ============================================================

def plot_1d(x, y, xlabel, ylabel, title, filename):
    plt.figure(figsize=(12, 6))
    plt.plot(x, y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, filename), dpi=200)
    plt.close()

def plot_2d(m, kz, S):
    plt.figure(figsize=(14, 8))
    plt.imshow(
        np.log10(S + 1e-12),
        origin="lower",
        aspect="auto",
        extent=[kz[0], kz[-1], m[0], m[-1]]
    )
    plt.xlabel("kz (1/Å)")
    plt.ylabel("Angular harmonic m")
    plt.title("QIP-QC — S(m,kz)")
    plt.colorbar(label="log10 S")
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, "S_m_kz.png"), dpi=200)
    plt.close()


# ============================================================
# MAIN EXECUTION
# ============================================================

def main():
    print("\n" + "=" * 80)
    print("QIP-QC PHASE 1D — ULTRA-FAST VECTORIZED (FFT Optimized)")
    print("=" * 80)

    download_dataset()
    density, voxel, origin = load_em_map()
    mask, threshold = threshold_density(density, percentile=95)
    
    coords = density_coordinates(mask, voxel)
    print("\nAligning principal axis...")
    coords = align_axis(coords)
    
    theta, radius, z = cylindrical_coordinates(coords)

    print("\n" + "=" * 70)
    print("STRUCTURAL GEOMETRY")
    print("=" * 70)
    print("N:", len(z))
    print("Radius mean:", np.mean(radius))
    print("Radius std:", np.std(radius))
    print("Z range:", z.min(), z.max())
    
    angular_x, angular = angular_profile(theta)
    print("Angular short correlation:", short_range_correlation(angular))

    # --------------------------------------------------------
    # FAST SPECTRA VIA FFT
    # --------------------------------------------------------
    print("\nCalculating angular spectrum (FFT)...")
    m, S_m = angular_structure_factor_fft(theta, m_max=100)

    print("Calculating axial spectrum (FFT)...")
    kz_target = np.linspace(0.01, 1.5, 400)
    S_kz = axial_structure_factor_fft(z, kz_target)

    print("Calculating 2D S(m,kz) (FFT2)...")
    m2, kz2, S2 = cylindrical_structure_factor_fft(theta, z, m_max=50)

    # Entropy
    Hm = spectral_entropy(S_m)
    Hk = spectral_entropy(S_kz)
    H2 = spectral_entropy(S2.ravel())

    # Dominant Harmonics
    order = np.argsort(S_m[1:])\n    dominant = order[::-1][:20] + 1

    zr, radial = radial_profile(radius, z)

    # --------------------------------------------------------
    # RESULTS
    # --------------------------------------------------------
    print("\n" + "=" * 80)
    print("QIP-QC PHASE 1 — STRUCTURAL SUMMARY")
    print("=" * 80)
    print("Dataset:", EMDB_ID, "/", PDB_ID)
    print("Real structural signal length:", len(z))
    print(f"Radius mean: {np.mean(radius):.8f}")
    print(f"Radius std: {np.std(radius):.8f}")
    print(f"Angular short-range corr: {short_range_correlation(angular):.8f}")
    print(f"Angular spectral entropy: {Hm:.8f}")
    print(f"Axial spectral entropy: {Hk:.8f}")
    print(f"2D spectral entropy: {H2:.8f}")

    print("\nDominant angular harmonics:")
    for mm in dominant:
        if mm < len(S_m):
            print(f"m={mm:3d} S={S_m[mm]:.8e}")

    # --------------------------------------------------------
    # SAVE & PLOT
    # --------------------------------------------------------
    np.savez_compressed(
        os.path.join(OUTDIR, "QIP_QC_EMD25156_vectorized.npz"),
        coordinates=coords, theta=theta, radius=radius, z=z,
        m=m, S_m=S_m, kz=kz_target, S_kz=S_kz,
        m2=m2, kz2=kz2, S2=S2,
        radial_z=zr, radial_profile=radial,
        angular_theta=angular_x, angular_profile=angular,
        threshold=threshold, voxel=voxel
    )

    plot_1d(kz_target, S_kz, "kz (1/Å)", "S(kz)", "Microtubule Axial Structure Factor", "axial_structure.png")
    plot_1d(m, S_m, "Angular harmonic m", "S(m)", "Microtubule Angular Structure Factor", "angular_structure.png")
    plot_1d(zr, radial, "z (Å)", "Mean radius (Å)", "Microtubule Radial Profile", "radial_profile.png")
    plot_1d(angular_x, angular, "theta", "Density count", "Microtubule Angular Density", "angular_density.png")
    plot_2d(m2, kz2, S2)

    print("\n" + "=" * 80)
    print("COMPLETE")
    print("Results:", OUTDIR)
    print("=" * 80)


if __name__ == "__main__":
    main()

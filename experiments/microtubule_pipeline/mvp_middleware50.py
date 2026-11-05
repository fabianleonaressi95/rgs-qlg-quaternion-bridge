"""
QIP-QC PHASE 1E — ULTRA-FAST VECTORIZED + LOCAL PATH & DENSITY WEIGHTS (MEMORY FIX)
===================================================================================

Complete cryo-EM microtubule spectral analysis
EMDB-25156 | PDB 7SJ7
"""

import os
import gzip
import urllib.request
import numpy as np
from pipeline_io import download_to_path as _download_to_path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from scipy.fft import rfft, fft2
except ImportError:
    from numpy.fft import rfft, fft2


# ============================================================
# CONFIGURATION & LOCAL PATHS
# ============================================================

EMDB_ID = "EMD-25156"
PDB_ID = "7SJ7"

OUTDIR = os.environ.get(
    "RGS_QIP_OUTDIR",
    os.path.join(os.path.expanduser("~"), "Scaricati", "QIP_QC_EMD25156"),
)
FIGDIR = os.path.join(OUTDIR, "figures")

os.makedirs(FIGDIR, exist_ok=True)

MAP_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/emdb/"
    "structures/EMD-25156/map/"
    "emd_25156.map.gz"
)
PDB_URL = f"https://files.rcsb.org/download/{PDB_ID}.pdb"

MAP_FILE_GZ = os.path.join(OUTDIR, "emd_25156.map.gz")
MAP_FILE_RAW = os.path.join(OUTDIR, "emd_25156.map")

if os.path.exists(MAP_FILE_RAW):
    MAP_FILE = MAP_FILE_RAW
elif os.path.exists(MAP_FILE_GZ):
    MAP_FILE = MAP_FILE_GZ
else:
    MAP_FILE = MAP_FILE_GZ

PDB_FILE = os.path.join(OUTDIR, f"{PDB_ID}.pdb")


# ============================================================
# FILE INTEGRITY & DATA LOADING
# ============================================================

def is_valid_gzip(filename):
    if not os.path.exists(filename):
        return False
    try:
        with gzip.open(filename, "rb") as f:
            while f.read(1024 * 1024):
                pass
        return True
    except Exception:
        return False

def download_file(url, filename):
    """Download atomically with a finite timeout and optional pinned SHA-256."""
    checksum_name = (
        "RGS_QIP_MAP_SHA256"
        if os.path.basename(filename).startswith("emd_25156")
        else "RGS_QIP_PDB_SHA256"
    )
    expected_sha256 = os.environ.get(checksum_name)

    if os.path.exists(filename):
        if filename.endswith(".gz") and not is_valid_gzip(filename):
            print(f"[WARNING] File gzip corrotto ({filename}). Riavvio download...")
            os.remove(filename)
        else:
            print(f"[OK] Trovato in locale: {filename}")
            if expected_sha256 is None:
                print(f"[WARN] {checksum_name} non impostato: integrità crittografica non verificata.")
            return

    if expected_sha256 is None:
        print(f"[WARN] {checksum_name} non impostato: il download sarà verificato per formato, non contro un digest atteso.")
    print(f"[DOWNLOAD] {url}")
    digest = _download_to_path(
        url,
        filename,
        timeout=120.0,
        expected_sha256=expected_sha256,
    )
    print(f"[OK] Completato: {filename} (SHA-256 {digest})")

def download_dataset():
    print("=" * 80)
    print(f"VERIFICA DATASET LOCALE IN: {OUTDIR}")
    print("=" * 80)

    if os.path.exists(MAP_FILE):
        print(f"[OK] Mappa Cryo-EM presente: {MAP_FILE}")
    else:
        print("[INFO] Mappa non trovata in locale. Scaricamento...")
        download_file(MAP_URL, MAP_FILE_GZ)

    if os.path.exists(PDB_FILE):
        print(f"[OK] Struttura PDB presente: {PDB_FILE}")
    else:
        print("[INFO] PDB non trovato in locale. Scaricamento...")
        download_file(PDB_URL, PDB_FILE)

def load_em_map():
    try:
        import mrcfile
    except ImportError:
        raise RuntimeError("Modulo 'mrcfile' mancante. Installa con:\npip install mrcfile")

    print(f"\nCaricamento mappa Cryo-EM da: {MAP_FILE}")
    with mrcfile.open(MAP_FILE, permissive=True) as mrc:
        density = np.asarray(mrc.data, dtype=np.float32)
        voxel = np.asarray(mrc.voxel_size.tolist(), dtype=np.float64)
        origin = np.asarray(mrc.header.origin.tolist(), dtype=np.float64)

    print("Dimensioni griglia (shape):", density.shape)
    print("Dimensione voxel (Å):", voxel)
    print(f"Range densità: min={float(density.min()):.4f}, max={float(density.max()):.4f}")
    return density, voxel, origin

def threshold_density(density, percentile=95):
    threshold = np.percentile(density, percentile)
    mask = density >= threshold
    print(f"\nSoglia di densità (percentile {percentile}%): {threshold:.4f}")
    print(f"Voxel occupati sopra soglia: {int(mask.sum())}")
    return mask, threshold


# ============================================================
# GEOMETRY & ALIGNMENT (MEMORY EFFICIENT)
# ============================================================

def density_coordinates(mask, voxel, density):
    print("\nEstrazione coordinate e pesi di densità (rho_j)...")
    indices = np.argwhere(mask)
    xyz = indices[:, ::-1] * voxel
    xyz = xyz - xyz.mean(axis=0, keepdims=True)
    
    weights = density[mask].astype(np.float32)
    return xyz.astype(np.float32), weights

def principal_axis(coords, weights):
    centered = coords - coords.mean(axis=0, keepdims=True)
    # FIX MEMORIA: O(N) memoria invece di O(N^2)
    weighted_centered = centered * weights[:, np.newaxis]
    covariance = (centered.T @ weighted_centered) / np.sum(weights)
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

def align_axis(coords, weights):
    axis = principal_axis(coords, weights)
    target = np.array([0.0, 0.0, 1.0], dtype=np.float32)
    R = rotation_matrix(axis, target)
    return (coords @ R.T).astype(np.float32)

def cylindrical_coordinates(coords):
    x, y, z = coords[:, 0], coords[:, 1], coords[:, 2]
    radius = np.hypot(x, y).astype(np.float32)
    theta = np.arctan2(y, x).astype(np.float32)
    return theta, radius, z


# ============================================================
# SPECTRAL CALCULATIONS (DENSITY-WEIGHTED FFT)
# ============================================================

def angular_structure_factor_fft(theta, weights, m_max=100, n_bins=2048):
    counts, _ = np.histogram(theta, bins=n_bins, range=(-np.pi, np.pi), weights=weights)
    fft_vals = rfft(counts)
    W = np.sum(weights)
    S = (np.abs(fft_vals) ** 2) / W 
    m_values = np.arange(0, min(m_max + 1, len(S)))
    return m_values, S[:len(m_values)]

def axial_structure_factor_fft(z, weights, kz_target, n_bins=4096):
    z_min, z_max = z.min(), z.max()
    counts, _ = np.histogram(z, bins=n_bins, range=(z_min, z_max), weights=weights)
    fft_vals = rfft(counts)
    W = np.sum(weights)
    S_full = (np.abs(fft_vals) ** 2) / W
    L = z_max - z_min
    kz_fft = 2 * np.pi * np.fft.rfftfreq(n_bins, d=L / n_bins)
    S_kz = np.interp(kz_target, kz_fft, S_full)
    return S_kz

def cylindrical_structure_factor_fft(theta, z, weights, m_max=50, n_kz_bins=512, n_theta_bins=512):
    z_min, z_max = z.min(), z.max()
    grid, _, _ = np.histogram2d(
        theta, z,
        bins=[n_theta_bins, n_kz_bins],
        range=[[-np.pi, np.pi], [z_min, z_max]],
        weights=weights
    )
    fft2d = fft2(grid)
    W = np.sum(weights)
    S2_full = (np.abs(fft2d) ** 2) / W
    m_arr = np.arange(0, m_max + 1)
    L = z_max - z_min
    kz_fft = 2 * np.pi * np.fft.fftfreq(n_kz_bins, d=L / n_kz_bins)
    pos_kz_idx = np.where((kz_fft >= 0.01) & (kz_fft <= 1.5))[0]
    sorted_idx = pos_kz_idx[np.argsort(kz_fft[pos_kz_idx])]
    kz_arr = kz_fft[sorted_idx]
    S2 = S2_full[m_arr[:, None], sorted_idx[None, :]]
    return m_arr, kz_arr, S2


# ============================================================
# METRICS & STATISTICAL PROFILES
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

def radial_profile(radius, z, weights, bins=856):
    edges = np.linspace(z.min(), z.max(), bins + 1)
    idx = np.digitize(z, edges) - 1
    valid = (idx >= 0) & (idx < bins)
    sums = np.bincount(idx[valid], weights=(radius[valid] * weights[valid]), minlength=bins)
    counts = np.bincount(idx[valid], weights=weights[valid], minlength=bins)
    profile = sums / np.maximum(counts, 1e-9)
    centers = (edges[:-1] + edges[1:]) / 2
    return centers, profile

def angular_profile(theta, weights, bins=360):
    edges = np.linspace(-np.pi, np.pi, bins + 1)
    idx = np.digitize(theta, edges) - 1
    valid = (idx >= 0) & (idx < bins)
    counts = np.bincount(idx[valid], weights=weights[valid], minlength=bins)
    centers = (edges[:-1] + edges[1:]) / 2
    return centers, counts


# ============================================================
# PLOTTING FUNCTIONS
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
    plt.ylabel("Armonica angolare m")
    plt.title("QIP-QC — Density-Weighted Spettro S(m,kz)")
    plt.colorbar(label="log10 S")
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, "S_m_kz.png"), dpi=200)
    plt.close()


# ============================================================
# MAIN EXECUTION
# ============================================================

def main():
    print("\n" + "=" * 80)
    print("QIP-QC PHASE 1E — SPECCRO PESATO DALLA DENSITÀ")
    print("=" * 80)

    download_dataset()
    density, voxel, origin = load_em_map()
    mask, threshold = threshold_density(density, percentile=95)

    coords, weights = density_coordinates(mask, voxel, density)
    print("\nAllineamento asse principale (pesato sulla massa/densità)...")
    coords = align_axis(coords, weights)

    theta, radius, z = cylindrical_coordinates(coords)

    print("\n" + "=" * 70)
    print("GEOMETRIA STRUTTURALE PESATA")
    print("=" * 70)
    print("Numero di voxel analizzati (N):", len(z))
    mean_radius = np.average(radius, weights=weights)
    std_radius = np.sqrt(np.average((radius - mean_radius)**2, weights=weights))
    print(f"Raggio medio pesato: {mean_radius:.5f} Å")
    print(f"Deviazione standard raggio pesato: {std_radius:.5f} Å")
    print(f"Estensione Z: da {z.min():.2f} a {z.max():.2f} Å")

    angular_x, angular = angular_profile(theta, weights)
    print(f"Correlazione angolare a corto raggio: {short_range_correlation(angular):.8f}")

    # --------------------------------------------------------
    # CALCOLO SPETTRI VIA FFT (PESATI)
    # --------------------------------------------------------
    print("\nCalcolo fattore di struttura angolare pesato (FFT)...")
    m, S_m = angular_structure_factor_fft(theta, weights, m_max=100)

    print("Calcolo fattore di struttura assiale pesato (FFT)...")
    kz_target = np.linspace(0.01, 1.5, 400)
    S_kz = axial_structure_factor_fft(z, weights, kz_target)

    print("Calcolo fattore di struttura 2D S(m,kz) pesato (FFT2)...")
    m2, kz2, S2 = cylindrical_structure_factor_fft(theta, z, weights, m_max=50)

    # Entropia spettrale
    Hm = spectral_entropy(S_m)
    Hk = spectral_entropy(S_kz)
    H2 = spectral_entropy(S2.ravel())

    # Armoniche dominanti
    order = np.argsort(S_m[1:])
    dominant = order[::-1][:20] + 1

    zr, radial = radial_profile(radius, z, weights)

    # --------------------------------------------------------
    # OUTPUT RISULTATI
    # --------------------------------------------------------
    print("\n" + "=" * 80)
    print("QIP-QC PHASE 1 — RIASSUNTO SPETTRALE PESATO")
    print("=" * 80)
    print("Dataset:", EMDB_ID, "/", PDB_ID)
    print("Percorso elaborato:", MAP_FILE)
    print(f"Raggio medio: {mean_radius:.8f} Å")
    print(f"Dev. Std raggio: {std_radius:.8f} Å")
    print(f"Entropia spettrale angolare (Hm): {Hm:.8f}")
    print(f"Entropia spettrale assiale (Hk): {Hk:.8f}")
    print(f"Entropia spettrale 2D (H2): {H2:.8f}")

    print("\nTop 20 armoniche angolari dominanti (pesate):")
    for mm in dominant:
        if mm < len(S_m):
            print(f"m={mm:3d}  -->  S(m) = {S_m[mm]:.8e}")

    # --------------------------------------------------------
    # SALVATAGGIO DATI E GRAFICI
    # --------------------------------------------------------
    np.savez_compressed(
        os.path.join(OUTDIR, "QIP_QC_EMD25156_weighted.npz"),
        coordinates=coords, weights=weights, theta=theta, radius=radius, z=z,
        m=m, S_m=S_m, kz=kz_target, S_kz=S_kz,
        m2=m2, kz2=kz2, S2=S2,
        radial_z=zr, radial_profile=radial,
        angular_theta=angular_x, angular_profile=angular,
        threshold=threshold, voxel=voxel
    )

    plot_1d(kz_target, S_kz, "kz (1/Å)", "S(kz)", "Fattore di Struttura Assiale (Pesato)", "axial_structure.png")
    plot_1d(m, S_m, "Armonica angolare m", "S(m)", "Fattore di Struttura Angolare (Pesato)", "angular_structure.png")
    plot_1d(zr, radial, "z (Å)", "Raggio medio (Å)", "Profilo Radiale (Pesato)", "radial_profile.png")
    plot_1d(angular_x, angular, "theta (rad)", "Densità integrata", "Profilo Angolare di Densità (Pesato)", "angular_density.png")
    plot_2d(m2, kz2, S2)

    print("\n" + "=" * 80)
    print("ELABORAZIONE COMPLETATA CON SUCCESSO")
    print("Risultati e figure salvati in:", OUTDIR)
    print("=" * 80)


if __name__ == "__main__":
    main()

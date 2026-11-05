"""
QIP-QC PHASE 1F — SYNTHETIC LATTICE GENERATOR & SPECTRAL DISTANCE
==================================================================

Generates Periodico, Fibonacci, and Random synthetic microtubule lattices,
computes their spectral profiles S(m,kz), calculates spectral distances relative 
to cryo-EM data (EMD-25156), and extracts the residual spectrum S_R(m,kz).
"""

import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from scipy.fft import fft2
except ImportError:
    from numpy.fft import fft2

# ============================================================
# CONFIGURATION & PATHS
# ============================================================

OUTDIR = os.environ.get(
    "RGS_QIP_OUTDIR",
    os.path.join(os.path.expanduser("~"), "Scaricati", "QIP_QC_EMD25156"),
)
FIGDIR = os.path.join(OUTDIR, "figures")
EM_DATA_FILE = os.path.join(OUTDIR, "QIP_QC_EMD25156_weighted.npz")

os.makedirs(FIGDIR, exist_ok=True)

# parametri fisici derivati dai dati reali
N_PROTOFILAMENTS = 13
MONOMER_HEIGHT = 40.0  # Å (subunità di tubulina alpha/beta)
HELICAL_SHIFT = 3.0    # Offset di 3 monomeri per giro (reticolo B canonico)
PHI = (1.0 + np.sqrt(5.0)) / 2.0  # Sezione aurea


# ============================================================
# SYNTHETIC LATTICE GENERATORS
# ============================================================

def generate_periodic_lattice(radius, z_range, num_proto=13, dz=40.0, shift=3.0):
    """Genera un reticolo cilindrico canonico B-lattice (13-3)."""
    z_min, z_max = z_range
    num_rings = int((z_max - z_min) / dz)
    z_steps = np.linspace(z_min, z_max, num_rings)
    
    theta_list, z_list, weights_list = [], [], []
    d_theta = 2 * np.pi / num_proto

    for ring_idx, z_val in enumerate(z_steps):
        for p in range(num_proto):
            # Sfalsamento elicoidale canonico
            th = p * d_theta + (ring_idx * shift / num_proto) * d_theta
            th = (th + np.pi) % (2 * np.pi) - np.pi
            
            theta_list.append(th)
            z_list.append(z_val)
            weights_list.append(1.0)

    return np.array(theta_list, dtype=np.float32), np.array(z_list, dtype=np.float32), np.array(weights_list, dtype=np.float32)


def generate_fibonacci_lattice(radius, z_range, num_proto=13, base_dz=40.0, shift=3.0, alpha=0.15):
    """Genera un reticolo cilindrico con modulazione assiale di Fibonacci (A/B = PHI)."""
    z_min, z_max = z_range
    total_length = z_max - z_min
    
    # Sequenza di Fibonacci per gli intervalli dz
    # L1 = base_dz, L2 = base_dz * PHI
    fib_sequence = []
    curr_z = z_min
    idx = 0
    while curr_z < z_max:
        # Alternanza basata sulla sequenza di Sturmian/Fibonacci
        step = base_dz * (1.0 + alpha * (1.0 if ((idx * PHI) % 1.0) < 0.5 else -1.0 / PHI))
        fib_sequence.append(curr_z)
        curr_z += step
        idx += 1
        
    z_steps = np.array(fib_sequence, dtype=np.float32)
    
    theta_list, z_list, weights_list = [], [], []
    d_theta = 2 * np.pi / num_proto

    for ring_idx, z_val in enumerate(z_steps):
        for p in range(num_proto):
            th = p * d_theta + (ring_idx * shift / num_proto) * d_theta
            th = (th + np.pi) % (2 * np.pi) - np.pi
            
            theta_list.append(th)
            z_list.append(z_val)
            # Pesi modulati aperiodicamente
            w = 1.0 + 0.2 * np.sin(2 * np.pi * ring_idx / PHI)
            weights_list.append(w)

    return np.array(theta_list, dtype=np.float32), np.array(z_list, dtype=np.float32), np.array(weights_list, dtype=np.float32)


def generate_random_lattice(radius, z_range, num_proto=13, base_dz=40.0, shift=3.0, noise_scale=8.0):
    """Genera un reticolo con perturbazione gaussiana stocastica."""
    z_min, z_max = z_range
    num_rings = int((z_max - z_min) / base_dz)
    z_steps = np.linspace(z_min, z_max, num_rings)
    
    theta_list, z_list, weights_list = [], [], []
    d_theta = 2 * np.pi / num_proto

    np.random.seed(42)  # Riproducibilità

    for ring_idx, z_val in enumerate(z_steps):
        for p in range(num_proto):
            z_noisy = z_val + np.random.normal(0, noise_scale)
            th_noisy = p * d_theta + (ring_idx * shift / num_proto) * d_theta + np.random.normal(0, 0.05)
            th_noisy = (th_noisy + np.pi) % (2 * np.pi) - np.pi
            
            theta_list.append(th_noisy)
            z_list.append(z_noisy)
            weights_list.append(1.0 + np.random.uniform(-0.2, 0.2))

    return np.array(theta_list, dtype=np.float32), np.array(z_list, dtype=np.float32), np.array(weights_list, dtype=np.float32)


# ============================================================
# SPECTRAL COMPUTATION
# ============================================================

def compute_2d_spectrum(theta, z, weights, m_max=50, n_kz_bins=512, n_theta_bins=512, z_range=None):
    if z_range is None:
        z_min, z_max = z.min(), z.max()
    else:
        z_min, z_max = z_range

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
# DISTANCE METRICS
# ============================================================

def relative_spectral_distance(S_target, S_model, eps=1e-8):
    """Calcola la distanza relativa in scala logaritmica (Norma di Frobenius)."""
    log_target = np.log10(S_target + eps)
    log_model = np.log10(S_model + eps)
    diff = log_target - log_model
    return np.linalg.norm(diff) / np.linalg.norm(log_target)


# ============================================================
# MAIN EXECUTION
# ============================================================

def main():
    print("=" * 80)
    print("QIP-QC PHASE 1F — SYNTHETIC BENCHMARK & RESIDUAL EXTRACTION")
    print("=" * 80)

    if not os.path.exists(EM_DATA_FILE):
        raise FileNotFoundError(f"File dati non trovato: {EM_DATA_FILE}\nEsegui prima lo script della Phase 1E.")

    em_data = np.load(EM_DATA_FILE)
    S2_em = em_data['S2']
    m_arr = em_data['m2']
    kz_arr = em_data['kz2']
    z_em = em_data['z']
    r_mean = float(np.average(em_data['radius'], weights=em_data['weights']))
    z_range = (float(z_em.min()), float(z_em.max()))

    print(f"[OK] Caricato spettro Cryo-EM reale: S2 shape = {S2_em.shape}")
    print(f"Raggio medio target: {r_mean:.2f} Å | Range Z: {z_range[0]:.2f} a {z_range[1]:.2f} Å")

    # 1. Generazione dei 3 Reticoli
    print("\nGenerazione Reticolo 1: Canonico Periodico (13-3 B-Lattice)...")
    th_p, z_p, w_p = generate_periodic_lattice(r_mean, z_range)
    _, _, S2_periodic = compute_2d_spectrum(th_p, z_p, w_p, m_max=50, n_kz_bins=512, n_theta_bins=512, z_range=z_range)

    print("Generazione Reticolo 2: Aperiodico Fibonacci (Gold Modulation)...")
    th_f, z_f, w_f = generate_fibonacci_lattice(r_mean, z_range)
    _, _, S2_fibonacci = compute_2d_spectrum(th_f, z_f, w_f, m_max=50, n_kz_bins=512, n_theta_bins=512, z_range=z_range)

    print("Generazione Reticolo 3: Stocastico Random (Gaussian Noise)...")
    th_r, z_r, w_r = generate_random_lattice(r_mean, z_range)
    _, _, S2_random = compute_2d_spectrum(th_r, z_r, w_r, m_max=50, n_kz_bins=512, n_theta_bins=512, z_range=z_range)

    # 2. Calcolo Matrice delle Distanze Spettrali Relative
    d_periodic = relative_spectral_distance(S2_em, S2_periodic)
    d_fibonacci = relative_spectral_distance(S2_em, S2_fibonacci)
    d_random = relative_spectral_distance(S2_em, S2_random)

    print("\n" + "=" * 80)
    print("MATRICE DI DISTANZA SPETTRALE RELATIVA D(S_EM, S_model)")
    print("=" * 80)
    print(f"1. D(EM, Periodic)  : {d_periodic:.6f}")
    print(f"2. D(EM, Fibonacci) : {d_fibonacci:.6f}")
    print(f"3. D(EM, Random)    : {d_random:.6f}")

    # 3. Estrazione Spettro Residuo S_R(m, kz)
    # S_R = | S_EM - S_Periodic |
    S_residual = np.abs(S2_em - S2_periodic)

    # Top armoniche del segnale residuo
    S_res_m = np.sum(S_residual, axis=1)
    top_res_m = np.argsort(S_res_m[1:])[::-1][:10] + 1

    print("\n" + "=" * 80)
    print("TOP HARMONICS IN RESIDUAL SPECTRUM S_R(m, kz)")
    print("=" * 80)
    for mm in top_res_m:
        print(f"m={mm:2d}  -->  Residuo integrato = {S_res_m[mm]:.8e}")

    # 4. Grafici Comparativi
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    im0 = axes[0, 0].imshow(np.log10(S2_em + 1e-12), origin="lower", aspect="auto", extent=[kz_arr[0], kz_arr[-1], m_arr[0], m_arr[-1]])
    axes[0, 0].set_title("Cryo-EM S(m,kz) (Sperimentale)")
    axes[0, 0].set_ylabel("Armonica m")
    plt.colorbar(im0, ax=axes[0, 0])

    im1 = axes[0, 1].imshow(np.log10(S2_periodic + 1e-12), origin="lower", aspect="auto", extent=[kz_arr[0], kz_arr[-1], m_arr[0], m_arr[-1]])
    axes[0, 1].set_title(f"Periodico Canonico (D={d_periodic:.4f})")
    plt.colorbar(im1, ax=axes[0, 1])

    im2 = axes[1, 0].imshow(np.log10(S2_fibonacci + 1e-12), origin="lower", aspect="auto", extent=[kz_arr[0], kz_arr[-1], m_arr[0], m_arr[-1]])
    axes[1, 0].set_title(f"Fibonacci Aperiodico (D={d_fibonacci:.4f})")
    axes[1, 0].set_xlabel("kz (1/Å)")
    axes[1, 0].set_ylabel("Armonica m")
    plt.colorbar(im2, ax=axes[1, 0])

    im3 = axes[1, 1].imshow(np.log10(S_residual + 1e-12), origin="lower", aspect="auto", extent=[kz_arr[0], kz_arr[-1], m_arr[0], m_arr[-1]])
    axes[1, 1].set_title("Spettro Residuo S_R(m,kz) = |EM - Periodico|")
    axes[1, 1].set_xlabel("kz (1/Å)")
    plt.colorbar(im3, ax=axes[1, 1])

    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, "spectral_benchmark_comparison.png"), dpi=200)
    plt.close()

    # Salvataggio Dati
    np.savez_compressed(
        os.path.join(OUTDIR, "QIP_QC_EMD25156_residual_bench.npz"),
        S_residual=S_residual,
        S_periodic=S2_periodic,
        S_fibonacci=S2_fibonacci,
        S_random=S2_random,
        distances={"periodic": d_periodic, "fibonacci": d_fibonacci, "random": d_random}
    )

    print("\n[OK] Analisi comparativa completata. Grafico salvato in:")
    print(os.path.join(FIGDIR, "spectral_benchmark_comparison.png"))


if __name__ == "__main__":
    main()

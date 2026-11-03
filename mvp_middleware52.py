"""
QIP-QC PHASE 2 — TIGHT-BINDING HAMILTONIAN & EXACT DIAGONALIZATION
==================================================================

Reconstructs real-space potential V(theta, z) from spectral residual S_R(m, kz),
constructs the 2D helical microtubule Tight-Binding Hamiltonian (13-3 B-lattice),
performs exact diagonalization, computes Band Structure, DOS, and Inverse 
Participation Ratio (IPR) to evaluate quantum localization.
"""

import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.linalg import eigh

# ============================================================
# CONFIGURATION & PATHS
# ============================================================

OUTDIR = "/home/fabianleo/Scaricati/QIP_QC_EMD25156"
FIGDIR = os.path.join(OUTDIR, "figures")
RESIDUAL_FILE = os.path.join(OUTDIR, "QIP_QC_EMD25156_residual_bench.npz")
EM_DATA_FILE = os.path.join(OUTDIR, "QIP_QC_EMD25156_weighted.npz")

os.makedirs(FIGDIR, exist_ok=True)

# Parametri Reticolo e Tight-Binding
N_PROTO = 13          # Numero di protofilamenti
N_RINGS = 64          # Numero di anelli lungo Z per la diagonalizzazione
N_SITES = N_PROTO * N_RINGS

E0 = 0.0              # Energia di sito base (eV)
T_PARALLEL = -1.0     # Hopping intra-protofilamento (longitudinale, Z) [eV]
T_PERP = -0.3         # Hopping inter-protofilamento (trasversale, theta) [eV]
T_SEAM = -0.15        # Hopping difettivo sulla giunzione A-type/B-type seam [eV]
HELICAL_SHIFT = 3     # Sfalsamento elicoidale B-lattice (3 monomeri per giro)
GAMMA_POTENTIAL = 0.25 # Accoppiamento con il potenziale residuo S_R [eV]


# ============================================================
# POTENTIAL RECONSTRUCTION FROM S_R(m, kz)
# ============================================================

def reconstruct_realspace_potential(s_residual, n_proto=13, n_rings=64):
    """
    Ricostruisce il potenziale reale V(theta_p, z_z) a partire dal residuo S_R(m, kz).
    Utilizza la radice dell'ampiezza spettrale con fase casuale deterministica/simmetrica.
    """
    m_max, n_kz = s_residual.shape
    
    # Griglia reale
    theta_grid = np.linspace(-np.pi, np.pi, n_proto, endpoint=False)
    z_grid = np.linspace(0, n_rings - 1, n_rings)
    
    V_real = np.zeros((n_proto, n_rings), dtype=np.float64)
    
    # Proiezione armonica inversa
    np.random.seed(1337)
    random_phases = np.random.uniform(0, 2 * np.pi, size=(m_max, n_kz))
    
    for m in range(m_max):
        amplitude_m = np.sqrt(s_residual[m, :])
        if np.sum(amplitude_m) < 1e-12:
            continue
        for ikz, kz_val in enumerate(np.linspace(0, np.pi, n_kz)):
            A = amplitude_m[ikz]
            phi = random_phases[m, ikz]
            
            # Somma delle componenti armoniche
            term = A * np.cos(m * theta_grid[:, None] + kz_val * z_grid[None, :] + phi)
            V_real += term
            
    # Normalizzazione std-dev
    if np.std(V_real) > 0:
        V_real = (V_real - np.mean(V_real)) / np.std(V_real)
        
    return V_real


# ============================================================
# HAMILTONIAN BUILDER
# ============================================================

def build_tight_binding_hamiltonian(V_potential, n_proto=13, n_rings=64, 
                                     t_par=-1.0, t_perp=-0.3, t_seam=-0.15, 
                                     shift=3, gamma=0.25):
    """
    Costruisce la matrice Hamiltoniana H_TB su reticolo cilindrico 2D con seam elicoidale.
    Indice del sito: i(p, z) = p * n_rings + z
    """
    n_sites = n_proto * n_rings
    H = np.zeros((n_sites, n_sites), dtype=np.float64)
    
    def site_index(p, z):
        return p * n_rings + z

    for p in range(n_proto):
        for z in range(n_rings):
            i = site_index(p, z)
            
            # 1. Energie di sito modulate dal potenziale residuo V(p, z)
            H[i, i] = E0 + gamma * V_potential[p, z]
            
            # 2. Hopping Longitudinale Intra-Protofilamento (Z ± 1)
            if z + 1 < n_rings:
                j_up = site_index(p, z + 1)
                H[i, j_up] += t_par
                H[j_up, i] += t_par
                
            # 3. Hopping Trasversale Inter-Protofilamento (p ± 1)
            # Condizioni al contorno azimutali + giunzione Seam (p = 0 <-> p = 12)
            if p < n_proto - 1:
                # Legame inter-protofilamento regolare
                j_right = site_index(p + 1, z)
                H[i, j_right] += t_perp
                H[j_right, i] += t_perp
            else:
                # Giunzione SEAM (p = 13-1 -> p = 0) con sfalsamento elicoidale Z + shift
                z_seam = z + shift
                if z_seam < n_rings:
                    j_seam = site_index(0, z_seam)
                    H[i, j_seam] += t_seam
                    H[j_seam, i] += t_seam

    return H


# ============================================================
# QUANTUM OBSERVABLES (IPR & DOS)
# ============================================================

def compute_ipr(eigenvectors):
    """
    Calcola l'Inverse Participation Ratio (IPR) per ogni autostatopsi_n:
    IPR_n = sum_i |psi_n(i)|^4
    IPR ~ 1/N (stato delocalizzato/esteso)
    IPR ~ 1   (stato fortemente localizzato)
    """
    return np.sum(np.abs(eigenvectors) ** 4, axis=0)


def compute_dos(eigenvalues, energy_grid, sigma=0.05):
    """Calcola la Densità degli Stati (DOS) con broadening gaussiano."""
    dos = np.zeros_like(energy_grid)
    for E_n in eigenvalues:
        dos += np.exp(-((energy_grid - E_n) ** 2) / (2 * sigma ** 2))
    dos /= (np.sqrt(2 * np.pi) * sigma)
    return dos


def compute_angular_projection(eigenvectors, n_proto=13, n_rings=64):
    """Calcola il valore atteso del momento angolare <m> per ciascun autostato."""
    m_values = np.zeros(eigenvectors.shape[1])
    theta = np.linspace(-np.pi, np.pi, n_proto, endpoint=False)
    
    for n in range(eigenvectors.shape[1]):
        psi = eigenvectors[:, n].reshape(n_proto, n_rings)
        prob_p = np.sum(np.abs(psi) ** 2, axis=1) # densità azimutale
        
        # Trasformata di Fourier azimutale
        fft_p = np.abs(np.fft.fft(prob_p)) ** 2
        m_eff = np.argmax(fft_p[:n_proto // 2 + 1])
        m_values[n] = m_eff
        
    return m_values


# ============================================================
# MAIN EXECUTION
# ============================================================

def main():
    print("=" * 80)
    print("QIP-QC PHASE 2 — TIGHT-BINDING MODEL & EXACT DIAGONALIZATION")
    print("=" * 80)

    # 1. Caricamento Residuo o Generazione di Fallback
    if os.path.exists(RESIDUAL_FILE):
        data = np.load(RESIDUAL_FILE)
        S_residual = data['S_residual']
        print(f"[OK] Caricato spettro residuo S_R: shape = {S_residual.shape}")
    else:
        print("[WARN] File residuo non trovato. Generazione sintetica di test S_R(m, kz)...")
        S_residual = np.zeros((51, 141))
        S_residual[13, 20] = 4.27e5
        S_residual[1, 10]  = 3.82e5
        S_residual[2, 5]   = 3.11e5

    # 2. Ricostruzione Potenziale Reale
    print("\nRicostruzione potenziale spaziale V(theta, z) da S_R(m, kz)...")
    V_potential = reconstruct_realspace_potential(S_residual, n_proto=N_PROTO, n_rings=N_RINGS)
    print(f" Potenziale V(p, z) calcolato: min={V_potential.min():.3f}, max={V_potential.max():.3f}, std={V_potential.std():.3f}")

    # 3. Costruzione Hamiltoniana
    print(f"\nAssemblaggio Hamiltoniana Tight-Binding (Dimensione: {N_SITES}x{N_SITES})...")
    print(f" Parametri: t_par={T_PARALLEL} eV | t_perp={T_PERP} eV | t_seam={T_SEAM} eV | gamma={GAMMA_POTENTIAL} eV")
    H_tb = build_tight_binding_hamiltonian(
        V_potential, n_proto=N_PROTO, n_rings=N_RINGS,
        t_par=T_PARALLEL, t_perp=T_PERP, t_seam=T_SEAM,
        shift=HELICAL_SHIFT, gamma=GAMMA_POTENTIAL
    )

    # 4. Diagonalizzazione Esatta
    print("\nEsecuzione Diagonalizzazione Esatta (eigh)...")
    evals, evecs = eigh(H_tb)
    print(f"[OK] Spettro energetico calcolato.")
    print(f" Band Gap / Estensione Energetica: E_min = {evals.min():.4f} eV | E_max = {evals.max():.4f} eV")

    # 5. Calcolo IPR, DOS e Momento Angolare
    print("\nCalcolo osservabili quantistiche (IPR, DOS, <m>)...")
    ipr = compute_ipr(evecs)
    energy_grid = np.linspace(evals.min() - 0.5, evals.max() + 0.5, 1000)
    dos = compute_dos(evals, energy_grid, sigma=0.04)
    m_proj = compute_angular_projection(evecs, n_proto=N_PROTO, n_rings=N_RINGS)

    ipr_mean = np.mean(ipr)
    ipr_max = np.max(ipr)
    print(f" IPR medio: {ipr_mean:.6f} (Soglia delocalizzata = {1.0/N_SITES:.6f})")
    print(f" IPR massimo: {ipr_max:.6f} (Picco di localizzazione)")

    # 6. Generazione Grafici
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # A. Potenziale Reale V(p, z)
    im0 = axes[0, 0].imshow(V_potential, aspect="auto", cmap="twilight", origin="lower")
    axes[0, 0].set_title("Potenziale Modulato Reale V(θ, z) da S_R")
    axes[0, 0].set_xlabel("Anelli lungo Z (z)")
    axes[0, 0].set_ylabel("Protofilamento (p)")
    plt.colorbar(im0, ax=axes[0, 0], label="V(p,z) [norm.]")

    # B. Spettro Energetico vs IPR (Localizzazione)
    sc = axes[0, 1].scatter(range(N_SITES), evals, c=ipr, cmap="plasma", s=18)
    axes[0, 1].set_title("Spettro Energetico & Localizzazione (IPR)")
    axes[0, 1].set_xlabel("Indice dell'Autostato n")
    axes[0, 1].set_ylabel("Energia E_n [eV]")
    plt.colorbar(sc, ax=axes[0, 1], label="IPR (Localizzazione)")

    # C. Densità degli Stati (DOS)
    axes[1, 0].plot(energy_grid, dos, color="navy", lw=1.8)
    axes[1, 0].fill_between(energy_grid, 0, dos, color="skyblue", alpha=0.5)
    axes[1, 0].set_title("Densità degli Stati (DOS)")
    axes[1, 0].set_xlabel("Energia E [eV]")
    axes[1, 0].set_ylabel("DOS (Stati / eV)")
    axes[1, 0].grid(True, linestyle="--", alpha=0.4)

    # D. Proiezione Armonica <m> vs Energia
    axes[1, 1].scatter(evals, m_proj, c=ipr, cmap="viridis", s=20, alpha=0.8)
    axes[1, 1].set_title("Proiezione Armonica Modale <m> vs Energia")
    axes[1, 1].set_xlabel("Energia E_n [eV]")
    axes[1, 1].set_ylabel("Armonica Azimutale m")
    axes[1, 1].grid(True, linestyle="--", alpha=0.4)

    plt.tight_layout()
    chart_path = os.path.join(FIGDIR, "tight_binding_spectrum_dos.png")
    plt.savefig(chart_path, dpi=200)
    plt.close()

    # Salvataggio Risultati
    np.savez_compressed(
        os.path.join(OUTDIR, "QIP_QC_PHASE2_TB_results.npz"),
        evals=evals,
        evecs=evecs,
        ipr=ipr,
        dos=dos,
        energy_grid=energy_grid,
        V_potential=V_potential,
        H_tb=H_tb
    )

    print("\n" + "=" * 80)
    print("FASE 2 COMPLETATA CON SUCCESSO")
    print("=" * 80)
    print(f"Risultati e figure salvati in:\n - {chart_path}\n - {os.path.join(OUTDIR, 'QIP_QC_PHASE2_TB_results.npz')}")


if __name__ == "__main__":
    main()

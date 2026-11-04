"""
QIP-QC PHASE 8 — ADVANCED CONTROL MODELS (FIBONACCI, DEFORMED, CORRELATED RANDOM)
=============================================================================

Builds and solves advanced control Hamiltonians to benchmark against the 13-3 model:
1. Fibonacci quasiperiodic chain modulation
2. Deterministic deformed periodic potential
3. Colored noise / PSD-matched correlated random potential
"""

import os
import numpy as np
from scipy.linalg import eigh
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTDIR = "/home/fabianleo/Scaricati/QIP_QC_EMD25156"
FIGDIR = os.path.join(OUTDIR, "figures")
TB_RESULTS_FILE = os.path.join(OUTDIR, "QIP_QC_PHASE2_TB_results.npz")

os.makedirs(FIGDIR, exist_ok=True)

N_PROTO = 13
N_RINGS = 64
N_SITES = N_PROTO * N_RINGS


def build_fibonacci_potential():
    """Genera un potenziale basato su sequenza di Fibonacci/sezione aurea."""
    phi_golden = (1.0 + np.sqrt(5.0)) / 2.0
    V_fib = np.zeros((N_PROTO, N_RINGS))
    for p in range(N_PROTO):
        for z in range(N_RINGS):
            # Modulazione quasiperiodica standard incommensurabile
            val = np.cos(2.0 * np.pi * p / phi_golden + 2.0 * np.pi * z / (phi_golden**2))
            V_fib[p, z] = val
    return V_fib * 0.5  # Riscalato per mantenere l'ordine di grandezza energetico


def build_deformed_periodic_potential():
    """Genera un potenziale periodico con deformazione armonica esplicita."""
    V_def = np.zeros((N_PROTO, N_RINGS))
    for p in range(N_PROTO):
        for z in range(N_RINGS):
            V_def[p, z] = 0.4 * np.sin(2.0 * np.pi * p / 3.0) * np.cos(2.0 * np.pi * z / 16.0)
    return V_def


def build_colored_random_potential(V_reference):
    """Genera un potenziale stocastico con lo stesso spettro di potenza (PSD) del 13-3."""
    # Trasformata di Fourier 2D del potenziale originale
    V_fft = np.fft.fft2(V_reference)
    amplitudes = np.abs(V_fft)
    
    # Fasi casuali uniformi preservando le ampiezze (PSD matching)
    phases = np.random.uniform(0, 2 * np.pi, size=V_reference.shape)
    V_surrogate_fft = amplitudes * np.exp(1j * phases)
    
    # Antitrasformata per ottenere il potenziale correlato spazialmente
    V_colored = np.real(np.fft.ifft2(V_surrogate_fft))
    # Normalizzazione della varianza
    V_colored = V_colored / np.std(V_colored) * np.std(V_reference)
    return V_colored


def build_hamiltonian(V_pot):
    H = np.zeros((N_SITES, N_SITES), dtype=np.complex128)
    t_parallel = -1.0
    t_perp = -0.3
    
    for p in range(N_PROTO):
        for z in range(N_RINGS):
            i = p * N_RINGS + z
            H[i, i] = V_pot[p, z]
            
            if z + 1 < N_RINGS:
                j = p * N_RINGS + (z + 1)
                H[i, j] = t_parallel
                H[j, i] = t_parallel
                
            p_next = (p + 1) % N_PROTO
            j = p_next * N_RINGS + z
            H[i, j] = t_perp
            H[j, i] = t_perp
            
    return H


def compute_phi_array(evals, evecs):
    n_states = len(evals)
    phi_metric = np.zeros(n_states)
    for n in range(n_states):
        psi_2d = evecs[:, n].reshape(N_PROTO, N_RINGS)
        rho_A = psi_2d @ psi_2d.conj().T
        tr_rho = np.trace(rho_A).real
        if tr_rho > 1e-12:
            rho_A /= tr_rho
        eig_rho = np.linalg.eigvalsh(rho_A)
        eig_rho = np.clip(eig_rho, 1e-15, 1.0)
        s_vn = -np.sum(eig_rho * np.log2(eig_rho))
        
        marginal_entropies = 0.0
        for p in range(N_PROTO):
            p_val = np.abs(rho_A[p, p])
            if p_val > 1e-15:
                marginal_entropies -= p_val * np.log2(p_val)
        phi_metric[n] = max(0.0, marginal_entropies - s_vn)
    return phi_metric


def main():
    print("=" * 80)
    print("QIP-QC PHASE 8 — ADVANCED CONTROL MODELS BENCHMARK")
    print("=" * 80)

    if not os.path.exists(TB_RESULTS_FILE):
        raise FileNotFoundError(f"File risultati non trovato: {TB_RESULTS_FILE}")

    data = np.load(TB_RESULTS_FILE)
    V_orig = data['V_potential']
    
    # 1. Risoluzione controlli avanzati
    print("\n[1/3] Risoluzione modello Fibonacci...")
    H_fib = build_hamiltonian(build_fibonacci_potential())
    evals_fib, evecs_fib = eigh(H_fib)
    phi_fib = compute_phi_array(evals_fib, evecs_fib)

    print("[2/3] Risoluzione modello Periodico Deformato...")
    H_def = build_hamiltonian(build_deformed_periodic_potential())
    evals_def, evecs_def = eigh(H_def)
    phi_def = compute_phi_array(evals_def, evecs_def)

    print("[3/3] Risoluzione modello Disordine Correlato (PSD-Matched)...")
    np.random.seed(101)
    H_col = build_hamiltonian(build_colored_random_potential(V_orig))
    evals_col, evecs_col = eigh(H_col)
    phi_col = compute_phi_array(evals_col, evecs_col)

    # Salvataggio dati avanzati
    np.savez_compressed(
        os.path.join(OUTDIR, "QIP_QC_PHASE8_controls_advanced.npz"),
        evals_fib=evals_fib, phi_fib=phi_fib,
        evals_def=evals_def, phi_def=phi_def,
        evals_col=evals_col, phi_col=phi_col
    )

    # Plot comparativo avanzato
    plt.figure(figsize=(12, 7))
    plt.scatter(evals_fib, phi_fib, c="crimson", s=10, alpha=0.5, label="Controllo Fibonacci (Quasicristallino)")
    plt.scatter(evals_def, phi_def, c="dodgerblue", s=10, alpha=0.5, label="Controllo Periodico Deformato")
    plt.scatter(evals_col, phi_col, c="forestgreen", s=10, alpha=0.4, label="Controllo Disordine Correlato (PSD)")
    
    plt.title("Confronto Avanzato: Modelli di Controllo vs Profilo Informazionale $\Phi(E)$")
    plt.xlabel("Energia E [eV]")
    plt.ylabel("$\Phi$ (Informazione Integrata)")
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.legend(loc="upper right")
    plt.tight_layout()
    
    chart_path = os.path.join(FIGDIR, "phase8_advanced_controls_phi.png")
    plt.savefig(chart_path, dpi=200)
    plt.close()

    print("\n" + "=" * 80)
    print("FASE 8 (CONTROLLI AVANZATI) COMPLETATA CON SUCCESSO")
    print("=" * 80)
    print(f" Grafico comparativo avanzato salvato in : {chart_path}")


if __name__ == "__main__":
    main()

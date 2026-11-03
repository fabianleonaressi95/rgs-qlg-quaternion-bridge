"""
QIP-QC PHASE 7 — STATISTICAL BOOTSTRAP & ENSEMBLE ROBUSTNESS OF PHI(E)
=============================================================================

Implements a randomized ensemble/bootstrap analysis to compute confidence 
intervals for the integrated information metric Phi(E), isolating robust 
structural effects from local numerical fluctuations.
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


def build_hamiltonian_from_potential(V_pot):
    """Costruisce l'Hamiltoniana Tight-Binding a partire da un potenziale V_pot dato."""
    H = np.zeros((N_SITES, N_SITES), dtype=np.complex128)
    t_parallel = -1.0
    t_perp = -0.3
    
    for p in range(N_PROTO):
        for z in range(N_RINGS):
            i = p * N_RINGS + z
            H[i, i] = V_pot[p, z]
            
            # Hopping longitudinale
            if z + 1 < N_RINGS:
                j = p * N_RINGS + (z + 1)
                H[i, j] = t_parallel
                H[j, i] = t_parallel
                
            # Hopping trasversale
            p_next = (p + 1) % N_PROTO
            j = p_next * N_RINGS + z
            H[i, j] = t_perp
            H[j, i] = t_perp
            
    return H


def compute_phi_array(evals, evecs):
    """Calcola l'indice Phi operativo per tutti gli autostati."""
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
    print("QIP-QC PHASE 7 — STATISTICAL BOOTSTRAP & ENSEMBLE ROBUSTNESS")
    print("=" * 80)

    if not os.path.exists(TB_RESULTS_FILE):
        raise FileNotFoundError(f"File risultati non trovato: {TB_RESULTS_FILE}")

    data = np.load(TB_RESULTS_FILE)
    V_original = data['V_potential']
    
    # Risoluzione del modello originale di riferimento
    H_orig = build_hamiltonian_from_potential(V_original)
    evals_orig, evecs_orig = eigh(H_orig)
    phi_orig = compute_phi_array(evals_orig, evecs_orig)

    # Configurazione del ciclo di ensemble / bootstrap
    n_ensembles = 20  # Numero di realizzazioni surrogate stocastiche
    print(f"\nEsecuzione ciclo di ensemble (N = {n_ensembles}) con surrogate randomizzati...")
    
    all_evals = []
    all_phis = []
    
    np.random.seed(42)
    for ens in range(n_ensembles):
        # Aggiunta di rumore gaussiano controllato (5% della deviazione standard)
        noise = np.random.normal(0.0, 0.05 * np.std(V_original), size=V_original.shape)
        V_surrogate = V_original + noise
        
        H_surr = build_hamiltonian_from_potential(V_surrogate)
        evals_s, evecs_s = eigh(H_surr)
        phi_s = compute_phi_array(evals_s, evecs_s)
        
        all_evals.append(evals_s)
        all_phis.append(phi_s)
        print(f"   [Ensemble {ens+1}/{n_ensembles}] Completato.")

    # Salvataggio dei risultati dell'ensemble
    np.savez_compressed(
        os.path.join(OUTDIR, "QIP_QC_PHASE7_bootstrap_results.npz"),
        evals_orig=evals_orig, phi_orig=phi_orig,
        all_evals=np.array(all_evals), all_phis=np.array(all_phis)
    )

    # Elaborazione statistica (Intervallo di confidenza al 95%) tramite binning energetico
    flatten_evals = np.concatenate(all_evals)
    flatten_phis = np.concatenate(all_phis)
    
    bins = np.linspace(-3.0, 3.0, 60)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    mean_phi = []
    ci_lower = []
    ci_upper = []
    
    for i in range(len(bins) - 1):
        mask = (flatten_evals >= bins[i]) & (flatten_evals < bins[i+1])
        if np.sum(mask) > 5:
            vals = flatten_phis[mask]
            mean_phi.append(np.mean(vals))
            ci_lower.append(np.percentile(vals, 2.5))
            ci_upper.append(np.percentile(vals, 97.5))
        else:
            mean_phi.append(np.nan)
            ci_lower.append(np.nan)
            ci_upper.append(np.nan)
            
    # Generazione grafico della banda di confidenza statistica
    plt.figure(figsize=(10, 6))
    plt.scatter(evals_orig, phi_orig, c="purple", s=12, alpha=0.6, label="Modello 13-3 Base", zorder=3)
    plt.plot(bin_centers, mean_phi, color="darkorange", lw=2, label="Media Ensemble (Bootstrap)", zorder=2)
    plt.fill_between(bin_centers, ci_lower, ci_upper, color="orange", alpha=0.25, label="Intervallo di Confidenza 95%", zorder=1)
    
    plt.title("Bootstrap Statistico e Banda di Confidenza di $\Phi(E)$")
    plt.xlabel("Energia E [eV]")
    plt.ylabel("$\Phi$ (Informazione Integrata)")
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.legend(loc="upper center")
    plt.tight_layout()
    
    chart_path = os.path.join(FIGDIR, "phase7_bootstrap_phi.png")
    plt.savefig(chart_path, dpi=200)
    plt.close()

    print("\n" + "=" * 80)
    print("FASE 7 (BOOTSTRAP STATISTICO) COMPLETATA CON SUCCESSO")
    print("=" * 80)
    print(f" Grafico bootstrap salvato in : {chart_path}")


if __name__ == "__main__":
    main()


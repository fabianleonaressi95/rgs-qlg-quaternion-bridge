"""
QIP-QC PHASE 4 — QUANTUM INFORMATIONAL METRICS (KAPPA & ENTROPY)
=============================================================================

Loads Tight-Binding eigenstates and calculates:
1. Spatial Quantum Correlation Matrix (kappa)
2. State Participation Entropy (informational spread)
3. Preliminary metrics for Quantum Integrated Information (QII / Phi analysis)
"""

import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================================================
# CONFIGURATION & PATHS
# ============================================================

OUTDIR = os.environ.get(
    "RGS_QIP_OUTDIR",
    os.path.join(os.path.expanduser("~"), "Scaricati", "QIP_QC_EMD25156"),
)
FIGDIR = os.path.join(OUTDIR, "figures")
TB_RESULTS_FILE = os.path.join(OUTDIR, "QIP_QC_PHASE2_TB_results.npz")

os.makedirs(FIGDIR, exist_ok=True)

N_PROTO = 13
N_RINGS = 64
N_SITES = N_PROTO * N_RINGS


def compute_informational_metrics(evals, evecs):
    """
    Calcola l'entropia di partecipazione e le correlazioni spaziali per ogni autostato.
    """
    n_states = len(evals)
    participation_entropy = np.zeros(n_states)
    
    # 1. Entropia di Shannon degli autostati (misura di localizzazione/informazione)
    # S_n = - sum_i |psi_i|^2 log(|psi_i|^2)
    for n in range(n_states):
        prob_density = np.abs(evecs[:, n])**2
        # Evitiamo log(0)
        prob_density = np.clip(prob_density, 1e-15, 1.0)
        participation_entropy[n] = -np.sum(prob_density * np.log(prob_density))
        
    return participation_entropy


def compute_spatial_correlation_matrix(evecs, selected_state_idx=416):
    """
    Calcola la matrice di correlazione spaziale kappa_ij per un dato autostato n:
    kappa_ij = |psi_i|^2 * |psi_j|^2 (mappa di accoppiamento informazionale locale)
    """
    psi_n = evecs[:, selected_state_idx]
    prob_density = np.abs(psi_n)**2
    prob_2d = prob_density.reshape(N_PROTO, N_RINGS)
    
    # Matrice di correlazione di densità spaziale (13 x 13 tra protofilamenti mediati)
    kappa_proto = np.zeros((N_PROTO, N_PROTO))
    for p1 in range(N_PROTO):
        for p2 in range(N_PROTO):
            kappa_proto[p1, p2] = np.sum(prob_2d[p1, :] * prob_2d[p2, :])
            
    return prob_2d, kappa_proto


def main():
    print("=" * 80)
    print("QIP-QC PHASE 4 — INFORMATIONAL METRICS & CORRELATION ANALYSIS")
    print("=" * 80)

    if not os.path.exists(TB_RESULTS_FILE):
        raise FileNotFoundError(f"File risultati non trovato: {TB_RESULTS_FILE}")

    data = np.load(TB_RESULTS_FILE)
    evals = data['evals']
    evecs = data['evecs']
    print(f"[OK] Caricati {len(evals)} autostati da H_TB.")

    # 1. Calcolo Entropia di Partecipazione
    print("\nCalcolo Entropia di Partecipazione degli autostati...")
    s_entropy = compute_informational_metrics(evals, evecs)
    
    # 2. Correlazione spaziale sul bulk (stato centrale di energia E ~ 0)
    mid_state = len(evals) // 2
    print(f"Estrazione matrice di correlazione kappa per lo stato centrale (indice {mid_state}, E = {evals[mid_state]:.3f} eV)...")
    prob_2d, kappa_proto = compute_spatial_correlation_matrix(evecs, selected_state_idx=mid_state)

    # 3. Salvataggio Risultati
    np.savez_compressed(
        os.path.join(OUTDIR, "QIP_QC_PHASE4_informational_results.npz"),
        participation_entropy=s_entropy,
        kappa_proto=kappa_proto,
        prob_2d=prob_2d
    )

    # 4. Generazione Grafici Analitici
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # A. Entropia di Partecipazione vs Energia
    axes[0].scatter(evals, s_entropy, c=s_entropy, cmap="viridis", s=12, alpha=0.8)
    axes[0].set_title("Entropia di Partecipazione S(E)")
    axes[0].set_xlabel("Energia E [eV]")
    axes[0].set_ylabel("Entropia di Shannon S")
    axes[0].grid(True, linestyle="--", alpha=0.4)

    # B. Matrice di Correlazione Protofilamento-Protofilamento kappa
    im1 = axes[1].imshow(kappa_proto, cmap="inferno", origin="lower")
    axes[1].set_title(f"Matrice Correlazione Protofilamenti kappa (Stato {mid_state})")
    axes[1].set_xlabel("Protofilamento p2")
    axes[1].set_ylabel("Protofilamento p1")
    plt.colorbar(im1, ax=axes[1], label="Correlazione di Densità")

    plt.tight_layout()
    chart_path = os.path.join(FIGDIR, "informational_metrics_phase4.png")
    plt.savefig(chart_path, dpi=200)
    plt.close()

    print("\n" + "=" * 80)
    print("FASE 4 (METRICHE INFORMAZIONALI) COMPLETATA CON SUCCESSO")
    print("=" * 80)
    print(f" Grafico salvato in : {chart_path}")
    print(f" Dati salvati in   : {os.path.join(OUTDIR, 'QIP_QC_PHASE4_informational_results.npz')}")


if __name__ == "__main__":
    main()

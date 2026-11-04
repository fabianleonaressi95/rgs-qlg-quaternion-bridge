"""
QIP-QC PHASE 5 — QUANTUM INTEGRATED INFORMATION (QII) & REDUCED DENSITY MATRICES
=============================================================================

Computes subsystem reduced density matrices via partial trace over the ring coordinates,
evaluating entanglement entropy and integrated information (Phi proxy) across the 13 protofilaments.
"""

import os
import numpy as np
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


def compute_reduced_density_matrix_and_phi(evals, evecs):
    """
    Per ogni autostato, calcola la matrice di densità ridotta sui protofilamenti 
    tramite traccia parziale sulle coordinate dei ring, ricavando l'entropia 
    di entanglement e la metrica di integrazione Phi.
    """
    n_states = len(evals)
    entanglement_entropy = np.zeros(n_states)
    phi_integration_metric = np.zeros(n_states)

    for n in range(n_states):
        # Autostato psi come matrice (N_PROTO, N_RINGS)
        psi_2d = evecs[:, n].reshape(N_PROTO, N_RINGS)
        
        # Matrice di densità ridotta sui protofilamenti: rho_A = psi_2d @ psi_2d^H
        rho_A = psi_2d @ psi_2d.conj().T
        
        # Normalizzazione della traccia
        tr_rho = np.trace(rho_A).real
        if tr_rho > 1e-12:
            rho_A /= tr_rho
            
        # Autovalori di rho_A per il calcolo dell'entropia di von Neumann
        # Utilizziamo eigh poiché rho_A è Hermitiana
        eig_rho = np.linalg.eigvalsh(rho_A)
        eig_rho = np.clip(eig_rho, 1e-15, 1.0) # Evita log(0)
        
        # Entropia di von Neumann S(rho_A) = - sum lambda_i log2(lambda_i)
        s_vn = -np.sum(eig_rho * np.log2(eig_rho))
        entanglement_entropy[n] = s_vn
        
        # Metrica Phi proxy: misura l'irriducibilità/integrità rispetto alla bipartizione ottimale
        # (confronto tra l'entropia globale e la somma delle entropie marginali dei singoli protofilamenti)
        marginal_entropies = 0.0
        for p in range(N_PROTO):
            p_val = np.abs(rho_A[p, p])
            if p_val > 1e-15:
                marginal_entropies -= p_val * np.log2(p_val)
                
        # Phi stimato come informazione eccedente (mutua informazione strutturale bipartita)
        phi_val = max(0.0, marginal_entropies - s_vn)
        phi_integration_metric[n] = phi_val

    return entanglement_entropy, phi_integration_metric


def main():
    print("=" * 80)
    print("QIP-QC PHASE 5 — QUANTUM INTEGRATED INFORMATION & RDM ANALYSIS")
    print("=" * 80)

    if not os.path.exists(TB_RESULTS_FILE):
        raise FileNotFoundError(f"File risultati non trovato: {TB_RESULTS_FILE}")

    data = np.load(TB_RESULTS_FILE)
    evals = data['evals']
    evecs = data['evecs']
    print(f"[OK] Caricati {len(evals)} autostati per l'analisi QII.")

    print("\nCalcolo delle Matrici di Densità Ridotte e metriche Phi in corso...")
    ent_entropy, phi_metric = compute_reduced_density_matrix_and_phi(evals, evecs)

    # Salvataggio dei risultati
    np.savez_compressed(
        os.path.join(OUTDIR, "QIP_QC_PHASE5_qii_results.npz"),
        entanglement_entropy=ent_entropy,
        phi_metric=phi_metric
    )

    # Generazione dei grafici analitici
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # A. Entropia di Entanglement vs Energia
    axes[0].scatter(evals, ent_entropy, c=ent_entropy, cmap="plasma", s=12, alpha=0.8)
    axes[0].set_title("Entropia di Entanglement S_vN(E)")
    axes[0].set_xlabel("Energia E [eV]")
    axes[0].set_ylabel("Entropia di von Neumann (bit)")
    axes[0].grid(True, linestyle="--", alpha=0.4)

    # B. Metrica di Informazione Integrata Phi vs Energia
    axes[1].scatter(evals, phi_metric, c=phi_metric, cmap="inferno", s=12, alpha=0.8)
    axes[1].set_title("Indice di Integrazione Quantistica Phi(E)")
    axes[1].set_xlabel("Energia E [eV]")
    axes[1].set_ylabel("Phi (Informazione Integrata)")
    axes[1].grid(True, linestyle="--", alpha=0.4)

    plt.tight_layout()
    chart_path = os.path.join(FIGDIR, "qii_integrated_information_phase5.png")
    plt.savefig(chart_path, dpi=200)
    plt.close()

    print("\n" + "=" * 80)
    print("FASE 5 (ANALISI QII & RDM) COMPLETATA CON SUCCESSO")
    print("=" * 80)
    print(f" Grafico salvato in : {chart_path}")
    print(f" Dati salvati in   : {os.path.join(OUTDIR, 'QIP_QC_PHASE5_qii_results.npz')}")


if __name__ == "__main__":
    main()

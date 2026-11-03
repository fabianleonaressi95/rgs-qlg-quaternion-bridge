"""
QIP-QC PHASE 3C — PEIERLS PHASE INJECTION & TOPOLOGICAL PHASE TRANSITION
=============================================================================

Injects a gauge phase phi_P into the transverse hopping t_perp (Aharonov-Bohm flux)
to break Time-Reversal Symmetry and scans the Bott Index B(E_F, phi_P).
"""

import os
import numpy as np
from scipy.linalg import eigh, logm

OUTDIR = "/home/fabianleo/Scaricati/QIP_QC_EMD25156"
TB_RESULTS_FILE = os.path.join(OUTDIR, "QIP_QC_PHASE2_TB_results.npz")

N_PROTO = 13
N_RINGS = 64
N_SITES = N_PROTO * N_RINGS

E0 = 0.0
T_PARALLEL = -1.0
T_PERP = -0.3
T_SEAM = -0.15
HELICAL_SHIFT = 3
GAMMA_POTENTIAL = 0.25


def build_complex_hamiltonian(V_potential, peierls_phase=0.0, gamma=0.25):
    """
    Costruisce l'Hamiltoniana complessa H(phi_P) integrando la fase di Peierls.
    """
    H = np.zeros((N_SITES, N_SITES), dtype=np.complex128)
    
    def site_index(p, z):
        return p * N_RINGS + z

    t_perp_complex = T_PERP * np.exp(1j * peierls_phase)

    for p in range(N_PROTO):
        for z in range(N_RINGS):
            i = site_index(p, z)
            H[i, i] = E0 + gamma * V_potential[p, z]
            
            # Hopping Longitudinale (Z)
            if z + 1 < N_RINGS:
                j_up = site_index(p, z + 1)
                H[i, j_up] += T_PARALLEL
                H[j_up, i] += np.conj(T_PARALLEL)
                
            # Hopping Trasversale con Fase di Peierls (theta)
            if p < N_PROTO - 1:
                j_right = site_index(p + 1, z)
                H[i, j_right] += t_perp_complex
                H[j_right, i] += np.conj(t_perp_complex)
            else:
                # Giunzione Seam (p=12 -> p=0)
                z_seam = z + HELICAL_SHIFT
                if z_seam < N_RINGS:
                    j_seam = site_index(0, z_seam)
                    H[i, j_seam] += T_SEAM * np.exp(1j * peierls_phase)
                    H[j_seam, i] += np.conj(T_SEAM * np.exp(1j * peierls_phase))

    return H


def compute_bott_index(evals, evecs, e_fermi=0.0):
    """Calcola l'invariante topologico Bott Index per matrici complesse."""
    occ_idx = np.where(evals < e_fermi)[0]
    if len(occ_idx) == 0 or len(occ_idx) == N_SITES:
        return 0.0
        
    P = evecs[:, occ_idx] @ evecs[:, occ_idx].T.conj()
    
    theta_coords = np.tile(np.linspace(0, 2 * np.pi, N_PROTO, endpoint=False), N_RINGS)
    z_coords = np.repeat(np.linspace(0, 2 * np.pi, N_RINGS, endpoint=False), N_PROTO)
    
    X = np.diag(np.exp(1j * theta_coords))
    Y = np.diag(np.exp(1j * z_coords))
    
    I = np.eye(N_SITES)
    U = P @ X @ P + (I - P)
    V = P @ Y @ P + (I - P)
    
    W = V @ U @ V.T.conj() @ U.T.conj()
    log_W = logm(W)
    bott_index = np.real(np.trace(log_W) / (2j * np.pi))
    
    return bott_index


def main():
    print("=" * 80)
    print("QIP-QC PHASE 3C — PEIERLS PHASE & BOTT INDEX SCAN")
    print("=" * 80)

    data = np.load(TB_RESULTS_FILE)
    V_potential = data['V_potential']

    # Fasi di Peierls da testare (in unita di rad): phi = [0, pi/4, pi/2, pi]
    phases = [0.0, np.pi / 4, np.pi / 2, np.pi]
    e_fermi_levels = [-0.5, 0.0, 0.5]

    for phi in phases:
        print(f"\n--- Analisi con Fase di Peierls phi_P = {phi/np.pi:.2f} * pi rad ---")
        H_complex = build_complex_hamiltonian(V_potential, peierls_phase=phi, gamma=GAMMA_POTENTIAL)
        evals, evecs = eigh(H_complex)
        
        for ef in e_fermi_levels:
            b_idx = compute_bott_index(evals, evecs, e_fermi=ef)
            print(f"  E_Fermi = {ef:+.2f} eV  |  Bott Index B = {b_idx:+.4f}  --> Round: {int(np.round(b_idx))}")

if __name__ == "__main__":
    main()

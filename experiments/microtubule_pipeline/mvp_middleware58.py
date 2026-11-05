"""
QIP-QC PHASE 3D — SPINFUL HAMILTONIAN (1664x1664) & Z2 TOPOLOGICAL INDEX (FIXED)
=============================================================================

Extends the Hilbert space to Spin-1/2 (Dimension: 2 * 832 = 1664).
Includes intrinsic Kane-Mele Spin-Orbit Coupling (SOC) and Rashba SOC.
Calculates the Z2 Topological Index with custom Pfaffian solver.
"""

import os
import numpy as np
from scipy.linalg import eigh

OUTDIR = os.environ.get(
    "RGS_QIP_OUTDIR",
    os.path.join(os.path.expanduser("~"), "Scaricati", "QIP_QC_EMD25156"),
)
TB_RESULTS_FILE = os.path.join(OUTDIR, "QIP_QC_PHASE2_TB_results.npz")

N_PROTO = 13
N_RINGS = 64
N_SITES = N_PROTO * N_RINGS
N_SPIN_SITES = 2 * N_SITES  # 1664 x 1664

E0 = 0.0
T_PARALLEL = -1.0
T_PERP = -0.3
T_SEAM = -0.15
HELICAL_SHIFT = 3
GAMMA_POTENTIAL = 0.25

# Parametri Spin-Orbita (eV)
LAMBDA_KM = 0.08     # Kane-Mele intrinsic SOC
LAMBDA_RASHBA = 0.05 # Rashba SOC


def pfaffian(A):
    """
    Calcola il Pfaffiano di una matrice complessa antisimmetrica A (A^T = -A).
    Utilizza la decomposizione per eliminazione gaussiana orientata alle coppie.
    """
    A = np.array(A, dtype=np.complex128, copy=True)
    n = A.shape[0]
    
    if n % 2 != 0:
        return 0.0 + 0.0j
        
    pf = 1.0 + 0.0j
    for k in range(0, n - 1, 2):
        # Pivotazione per garantire stabilità numerica
        max_idx = np.argmax(np.abs(A[k:, k])) + k
        if max_idx != k + 1:
            # Scambio righe e colonne per mantenere l'antisimmetria
            A[[k + 1, max_idx], :] = A[[max_idx, k + 1], :]
            A[:, [k + 1, max_idx]] = A[:, [max_idx, k + 1]]
            pf = -pf
            
        pivot = A[k, k + 1]
        if np.abs(pivot) < 1e-12:
            return 0.0 + 0.0j
            
        pf *= pivot
        
        # Eliminazione sui blocchi rimanenti
        if k + 2 < n:
            tau = A[k + 2:, k] / pivot
            A[k + 2:, k + 2:] += np.outer(tau, A[k + 2:, k + 1]) - np.outer(A[k + 2:, k + 1], tau)
            
    return pf


def build_spinful_hamiltonian(V_potential, lambda_km=0.08, lambda_r=0.05, gamma=0.25):
    """
    Costruisce l'Hamiltoniana con Spin 1/2 (Dimensione 1664x1664).
    """
    H = np.zeros((N_SPIN_SITES, N_SPIN_SITES), dtype=np.complex128)
    
    sigma_0 = np.eye(2, dtype=np.complex128)
    sigma_x = np.array([[0, 1], [1, 0]], dtype=np.complex128)
    sigma_y = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
    sigma_z = np.array([[1, 0], [0, -1]], dtype=np.complex128)

    for p in range(N_PROTO):
        for z in range(N_RINGS):
            site_base = (p * N_RINGS + z) * 2
            
            # 1. Energia di sito + Potenziale Residuo S_R
            H[site_base:site_base+2, site_base:site_base+2] = (E0 + gamma * V_potential[p, z]) * sigma_0
            
            # 2. Hopping Longitudinale Intra-Protofilamento (Z) + Rashba Z
            if z + 1 < N_RINGS:
                target_base = (p * N_RINGS + (z + 1)) * 2
                h_hop = T_PARALLEL * sigma_0 - 1j * lambda_r * sigma_y
                H[site_base:site_base+2, target_base:target_base+2] += h_hop
                H[target_base:target_base+2, site_base:site_base+2] += h_hop.T.conj()
                
            # 3. Hopping Trasversale Inter-Protofilamento (theta) + Kane-Mele SOC
            if p < N_PROTO - 1:
                target_base = ((p + 1) * N_RINGS + z) * 2
                h_hop = T_PERP * sigma_0 + 1j * lambda_km * sigma_z + 1j * lambda_r * sigma_x
                H[site_base:site_base+2, target_base:target_base+2] += h_hop
                H[target_base:target_base+2, site_base:site_base+2] += h_hop.T.conj()
            else:
                # Giunzione SEAM (p = 12 -> p = 0)
                z_seam = z + HELICAL_SHIFT
                if z_seam < N_RINGS:
                    target_base = (0 * N_RINGS + z_seam) * 2
                    h_hop = T_SEAM * sigma_0 + 1j * lambda_km * sigma_z
                    H[site_base:site_base+2, target_base:target_base+2] += h_hop
                    H[target_base:target_base+2, site_base:site_base+2] += h_hop.T.conj()

    return H


def compute_z2_index_fukui(evals, evecs, e_fermi=0.0):
    """
    Calcola l'invariante Z2 valutando la matrice di Time-Reversal W_ij = <psi_i | Theta | psi_j>.
    """
    occ_idx = np.where(evals < e_fermi)[0]
    n_occ = len(occ_idx)
    
    if n_occ == 0 or n_occ % 2 != 0:
        return 0, 0.0
        
    psi_occ = evecs[:, occ_idx]
    
    # Operatore di Inversione Temporale Theta = i * sigma_y (x) I_N * K
    theta_psi = np.zeros_like(psi_occ)
    for i in range(N_SITES):
        idx_up = 2 * i
        idx_dn = 2 * i + 1
        theta_psi[idx_up, :] = -np.conj(psi_occ[idx_dn, :])
        theta_psi[idx_dn, :] =  np.conj(psi_occ[idx_up, :])
        
    # Matrice antisimmetrica W_ij = <psi_i | Theta | psi_j>
    W = psi_occ.T.conj() @ theta_psi
    
    # Calcolo del Pfaffiano e della fase topologica
    pf_W = pfaffian(W)
    det_W = np.abs(np.linalg.det(W))
    
    # L'invariante Z2 nu = 0 (banale) o nu = 1 (topologico)
    phase = np.angle(pf_W)
    z2_nu = 1 if np.abs(np.abs(phase) - np.pi) < 0.5 else 0
    
    return z2_nu, det_W


def main():
    print("=" * 80)
    print("QIP-QC PHASE 3D — SPINFUL HAMILTONIAN (1664x1664) & Z2 TOPOLOGICAL INDEX")
    print("=" * 80)

    if not os.path.exists(TB_RESULTS_FILE):
        raise FileNotFoundError(f"File risultati Fase 2 non trovato: {TB_RESULTS_FILE}")

    data = np.load(TB_RESULTS_FILE)
    V_potential = data['V_potential']

    print(f"\nAssemblaggio Hamiltoniana 1664x1664 (Spin 1/2)...")
    print(f" Parametri SOC: Kane-Mele lambda_KM = {LAMBDA_KM} eV | Rashba lambda_R = {LAMBDA_RASHBA} eV")
    
    H_spin = build_spinful_hamiltonian(V_potential, lambda_km=LAMBDA_KM, lambda_r=LAMBDA_RASHBA, gamma=GAMMA_POTENTIAL)
    
    print("\nDiagonalizzazione Esatta in corso...")
    evals, evecs = eigh(H_spin)
    print(f"[OK] Diagonalizzazione completata. Range Spettro: [{evals.min():.4f}, {evals.max():.4f}] eV")

    e_fermi_levels = [-1.0, -0.5, 0.0, 0.5, 1.0]
    
    print("\nCalcolo Invariante Topologico Z2 (Kane-Mele Index):")
    print("-" * 65)
    for ef in e_fermi_levels:
        z2, det_val = compute_z2_index_fukui(evals, evecs, e_fermi=ef)
        status = "ISOLANTE TOPOLOGICO (nu = 1)" if z2 == 1 else "BANALE (nu = 0)"
        print(f" E_Fermi = {ef:+.2f} eV | Z2 Index nu = {z2} | Det[W] = {det_val:.6f} --> {status}")
    print("-" * 65)


if __name__ == "__main__":
    main()

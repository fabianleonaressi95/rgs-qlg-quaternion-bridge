"""
QIP-QC PHASE 3B — TOPOLOGICAL BOTT INDEX CALCULATION
=============================================================================

Calculates the Bott Index B for the projector matrix P = \sum_{E_n < E_F} |n><n|
to verify non-trivial topological invariance and protection against local defects.
"""

import os
import numpy as np
from scipy.linalg import logm

OUTDIR = os.environ.get(
    "RGS_QIP_OUTDIR",
    os.path.join(os.path.expanduser("~"), "Scaricati", "QIP_QC_EMD25156"),
)
TB_RESULTS_FILE = os.path.join(OUTDIR, "QIP_QC_PHASE2_TB_results.npz")

N_PROTO = 13
N_RINGS = 64
N_SITES = N_PROTO * N_RINGS

def compute_bott_index(evals, evecs, e_fermi=0.0):
    """
    Calcola il Bott Index B = (1 / 2pi i) Tr(log(V U V^\dagger U^\dagger))
    dove U e V sono gli operatori di posizione proiettati.
    """
    # 1. Costruzione Matrice di Proiezione P sugli stati occupati (E < E_F)
    occ_idx = np.where(evals < e_fermi)[0]
    P = evecs[:, occ_idx] @ evecs[:, occ_idx].T.conj()
    
    # 2. Operatori Posizione Angolare (Theta) e Assiale (Z) resi periodici su [0, 2pi)
    theta_coords = np.tile(np.linspace(0, 2 * np.pi, N_PROTO, endpoint=False), N_RINGS)
    z_coords = np.repeat(np.linspace(0, 2 * np.pi, N_RINGS, endpoint=False), N_PROTO)
    
    X = np.diag(np.exp(1j * theta_coords))
    Y = np.diag(np.exp(1j * z_coords))
    
    # 3. Proiezione nello spazio occupato: U = P X P + (I - P), V = P Y P + (I - P)
    I = np.eye(N_SITES)
    U = P @ X @ P + (I - P)
    V = P @ Y @ P + (I - P)
    
    # 4. Calcolo Commutatore Modale W = V U V^\dagger U^\dagger
    W = V @ U @ V.T.conj() @ U.T.conj()
    
    # 5. Bott Index
    log_W = logm(W)
    bott_index = np.real(np.trace(log_W) / (2j * np.pi))
    
    return bott_index

def main():
    print("=" * 80)
    print("QIP-QC PHASE 3B — TOPOLOGICAL BOTT INDEX")
    print("=" * 80)

    data = np.load(TB_RESULTS_FILE)
    evals = data['evals']
    evecs = data['evecs']

    # Scansione del Bott Index al variare del Livello di Fermi E_F
    e_fermi_range = np.linspace(-1.5, 1.5, 7)
    
    print("\nCalcolo Bott Index in funzione del Livello di Fermi E_F:")
    print("-" * 50)
    for ef in e_fermi_range:
        b_idx = compute_bott_index(evals, evecs, e_fermi=ef)
        print(f" E_Fermi = {ef:+.2f} eV  |  Bott Index B = {b_idx:+.4f} (Arrotondato: {int(np.round(b_idx))})")
    print("-" * 50)

if __name__ == "__main__":
    main()

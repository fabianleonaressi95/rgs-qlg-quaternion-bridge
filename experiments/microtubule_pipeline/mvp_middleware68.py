"""
QIP-QC PHASE 10 — MULTI-UNIT NETWORK TRANSPORT SIMULATION (COMPLETE)
=============================================================================
Simulates external energy wavepacket propagation across a 2D matrix 
of interconnected 13-3 tubular units using sparse matrix evolution.
"""

import os
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTDIR = os.environ.get(
    "RGS_QIP_OUTDIR",
    os.path.join(os.path.expanduser("~"), "Scaricati", "QIP_QC_EMD25156"),
)
FIGDIR = os.path.join(OUTDIR, "figures")
os.makedirs(FIGDIR, exist_ok=True)

N_PROTO = 13
N_RINGS = 32  # Ridotto per efficienza della matrice di rete
N_SITES_UNIT = N_PROTO * N_RINGS

# Configurazione della matrice di rete (griglia 2x2 = 4 unità)
GRID_X = 2
GRID_Y = 2
N_UNITS = GRID_X * GRID_Y
TOTAL_SITES = N_UNITS * N_SITES_UNIT


def build_single_unit_sparse(V_pot):
    """Costruisce l'Hamiltoniana sparsa per una singola unità 13-3."""
    rows, cols, data = [], [], []
    t_parallel = -1.0
    t_perp = -0.3

    for p in range(N_PROTO):
        for z in range(N_RINGS):
            i = p * N_RINGS + z
            rows.append(i)
            cols.append(i)
            data.append(V_pot[p, z])

            if z + 1 < N_RINGS:
                j = p * N_RINGS + (z + 1)
                rows.extend([i, j])
                cols.extend([j, i])
                data.extend([t_parallel, t_parallel])

            p_next = (p + 1) % N_PROTO
            j = p_next * N_RINGS + z
            rows.extend([i, j])
            cols.extend([j, i])
            data.extend([t_perp, t_perp])

    return sp.coo_matrix((data, (rows, cols)), shape=(N_SITES_UNIT, N_SITES_UNIT)).tocsr()


def build_network_hamiltonian():
    """Assembla la supra-Hamiltoniana per la matrice di unità 13-3."""
    np.random.seed(42)
    p_indices = np.linspace(0, 2 * np.pi, N_PROTO)
    z_indices = np.linspace(0, 2 * np.pi, N_RINGS)
    P, Z = np.meshgrid(p_indices, z_indices, indexing='ij')
    V_base = 0.5 * np.cos(13 * P - 2 * Z)

    unit_H = build_single_unit_sparse(V_base)
    
    # Assemblaggio dei blocchi diagonali usando lil_matrix per evitare avvisi di efficienza
    blocks = [unit_H for _ in range(N_UNITS)]
    H_net = sp.block_diag(blocks, format='lil')

    # Aggiunta accoppiamenti inter-tubo
    t_inter = -0.15
    for uy in range(GRID_Y):
        for ux in range(GRID_X):
            u_idx = uy * GRID_X + ux
            if ux + 1 < GRID_X:
                u_right = uy * GRID_X + (ux + 1)
                for p in range(N_PROTO):
                    s_left = u_idx * N_SITES_UNIT + p * N_RINGS + (N_RINGS - 1)
                    s_right = u_right * N_SITES_UNIT + p * N_RINGS + 0
                    H_net[s_left, s_right] = t_inter
                    H_net[s_right, s_left] = t_inter

    return H_net.tocsr()


def main():
    print("=" * 80)
    print("QIP-QC PHASE 10 — MULTI-UNIT NETWORK TRANSPORT SIMULATION")
    print("=" * 80)

    print("\n[1/3] Assemblaggio della supra-Hamiltoniana di rete...")
    H_net = build_network_hamiltonian()
    print(f"   Dimensione totale matrice di rete: {H_net.shape[0]} siti.")

    print("\n[2/3] Simulazione propagazione pacchetto d'onda temporale...")
    psi_0 = np.zeros(TOTAL_SITES, dtype=np.complex128)
    psi_0[0] = 1.0  # Eccitazione localizzata all'ingresso della matrice

    times = np.linspace(0.0, 30.0, 60)
    probabilities_per_unit = np.zeros((len(times), N_UNITS))

    for idx, t in enumerate(times):
        psi_t = sla.expm_multiply(-1j * H_net * t, psi_0)
        
        for u in range(N_UNITS):
            start_idx = u * N_SITES_UNIT
            end_idx = (u + 1) * N_SITES_UNIT
            probabilities_per_unit[idx, u] = np.sum(np.abs(psi_t[start_idx:end_idx])**2)

    print("[3/3] Generazione grafico di transito energetico nella matrice...")
    plt.figure(figsize=(10, 6))
    for u in range(N_UNITS):
        plt.plot(times, probabilities_per_unit[:, u], lw=2, label=f"Unità Tubolare #{u+1}")

    plt.title("Flusso di Energia Esterna attraverso la Matrice di Unità $13\text{-}3$")
    plt.xlabel("Tempo Evolutivo [unità arbitrarie]")
    plt.ylabel("Probabilità / Densità di Energia Integrata")
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.legend()
    plt.tight_layout()
    
    chart_path = os.path.join(FIGDIR, "phase10_network_transport.png")
    plt.savefig(chart_path, dpi=200)
    plt.close()

    print("\n" + "=" * 80)
    print("SIMULAZIONE DI RETE MULTI-UNITÀ COMPLETATA CON SUCCESSO")
    print("=" * 80)
    print(f" Grafico di transito salvato in : {chart_path}")


if __name__ == "__main__":
    main()

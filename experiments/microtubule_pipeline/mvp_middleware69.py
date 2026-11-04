"""
QIP-QC PHASE 11 — THERMAL DECOHERENCE & NOISE ROBUSTNESS SIMULATION
=============================================================================
Evaluates the persistence of energy routing across the 13-3 multi-unit network 
under varying amplitudes of thermal/stochastic site disorder (sigma_th).
"""

import os
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTDIR = "/home/fabianleo/Scaricati/QIP_QC_EMD25156"
FIGDIR = os.path.join(OUTDIR, "figures")
os.makedirs(FIGDIR, exist_ok=True)

N_PROTO = 13
N_RINGS = 32
N_SITES_UNIT = N_PROTO * N_RINGS

GRID_X = 2
GRID_Y = 2
N_UNITS = GRID_X * GRID_Y
TOTAL_SITES = N_UNITS * N_SITES_UNIT


def build_base_network():
    """Costruisce la matrice di rete di base."""
    np.random.seed(42)
    p_indices = np.linspace(0, 2 * np.pi, N_PROTO)
    z_indices = np.linspace(0, 2 * np.pi, N_RINGS)
    P, Z = np.meshgrid(p_indices, z_indices, indexing='ij')
    V_base = 0.5 * np.cos(13 * P - 2 * Z)

    rows, cols, data = [], [], []
    t_parallel = -1.0
    t_perp = -0.3

    for p in range(N_PROTO):
        for z in range(N_RINGS):
            i = p * N_RINGS + z
            rows.append(i); cols.append(i); data.append(V_base[p, z])
            if z + 1 < N_RINGS:
                j = p * N_RINGS + (z + 1)
                rows.extend([i, j]); cols.extend([j, i]); data.extend([t_parallel, t_parallel])
            p_next = (p + 1) % N_PROTO
            j = p_next * N_RINGS + z
            rows.extend([i, j]); cols.extend([j, i]); data.extend([t_perp, t_perp])

    unit_H = sp.coo_matrix((data, (rows, cols)), shape=(N_SITES_UNIT, N_SITES_UNIT)).tocsr()
    blocks = [unit_H for _ in range(N_UNITS)]
    H_net = sp.block_diag(blocks, format='lil')

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
    print("QIP-QC PHASE 11 — THERMAL DECOHERENCE & NOISE ROBUSTNESS")
    print("=" * 80)

    H_base = build_base_network()
    noise_levels = [0.0, 0.1, 0.3, 0.6, 1.0]
    target_time = 25.0  # Istante temporale di osservazione del transito

    psi_0 = np.zeros(TOTAL_SITES, dtype=np.complex128)
    psi_0[0] = 1.0

    plt.figure(figsize=(10, 6))
    
    np.random.seed(100)
    for sigma in noise_levels:
        print(f"\n[+] Simulazione con ampiezza rumore termico sigma_th = {sigma:.1f}...")
        
        # Aggiunta di disordine stocastico diagonale (rumore termico)
        noise_vector = np.random.normal(0.0, sigma, size=TOTAL_SITES)
        H_noisy = H_base + sp.diags(noise_vector, format='csr')

        times = np.linspace(0.0, target_time, 40)
        prob_unit_2 = []

        for t in times:
            psi_t = sla.expm_multiply(-1j * H_noisy * t, psi_0)
            # Calcola la probabilità accumulata nell'Unità Tubolare #2
            start_u2 = N_SITES_UNIT
            end_u2 = 2 * N_SITES_UNIT
            p2 = np.sum(np.abs(psi_t[start_u2:end_u2])**2)
            prob_unit_2.append(p2)

        plt.plot(times, prob_unit_2, lw=2, label=f"Rumore $\\sigma_{{th}} = {sigma:.1f}$")

    plt.title("Robustezza del Trasporto Energetico sotto Rumore Termico (Unità #2)")
    plt.xlabel("Tempo Evolutivo [unità arbitrarie]")
    plt.ylabel("Probabilità di Trasferimento (Unità #2)")
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.legend()
    plt.tight_layout()

    chart_path = os.path.join(FIGDIR, "phase11_thermal_robustness.png")
    plt.savefig(chart_path, dpi=200)
    plt.close()

    print("\n" + "=" * 80)
    print("FASE 11 (ANALISI ROBUSTEZZA TERMICA) COMPLETATA CON SUCCESSO")
    print("=" * 80)
    print(f" Grafico di robustezza termica salvato in : {chart_path}")


if __name__ == "__main__":
    main()

"""
QIP-QC PHASE 3A — TIMECUT WAVEPACKET DYNAMICS & TRANSPORT REGIME
================================================================

Simulates time evolution |psi(t)> = exp(-i H t / hbar) |psi(0)> for a 
localized Gaussian wavepacket on the microtubule lattice.
Computes Mean Position <Z(t)>, Mean Squared Displacement (MSD), 
Survival Probability P_surv(t), and identifies transport regime (Ballistic vs Diffusive).
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

# Parametri temporali (Unità: hbar = 1, t in femtosecondi / unita adimensionali)
HBAR = 0.65821195  # eV * fs
T_MAX = 40.0       # fs
N_TSTEPS = 200     # Passi temporali

# Parametri Pacchetto d'Onde Iniziale
Z_0 = 10           # Anello iniziale (vicino all'estremita per osservare la propagazione)
P_0 = 6            # Protofilamento centrale
SIGMA_Z = 2.0      # Dispersione spaziale iniziale lungo Z
K_Z0 = 1.0         # Impulso iniziale lungo Z


# ============================================================
# INITIAL WAVEPACKET CREATION
# ============================================================

def create_gaussian_wavepacket(z0=10, p0=6, sigma_z=2.0, kz0=1.0, n_proto=13, n_rings=64):
    """
    Crea un pacchetto d'onde Gaussiano localizzato attorno a (p0, z0)
    con impulso iniziale kz0 lungo l'asse assiale Z.
    """
    psi0 = np.zeros((n_proto, n_rings), dtype=np.complex128)
    
    for p in range(n_proto):
        for z in range(n_rings):
            # Inviluppo Gaussiano azimutale e assiale
            dist_p = min(abs(p - p0), n_proto - abs(p - p0)) # distanza ciclica
            dist_z = z - z0
            
            spatial_envelope = np.exp(- (dist_z**2) / (4 * sigma_z**2)) * np.exp(- (dist_p**2) / 2.0)
            phase = np.exp(1j * kz0 * z)
            
            psi0[p, z] = spatial_envelope * phase

    # Normalizzazione della norma di stato
    psi0_flat = psi0.flatten()
    norm = np.linalg.norm(psi0_flat)
    if norm > 0:
        psi0_flat /= norm
        
    return psi0_flat


# ============================================================
# MAIN EXECUTION
# ============================================================

def main():
    print("=" * 80)
    print("QIP-QC PHASE 3A — QUANTUM WAVEPACKET DYNAMICS SIMULATION")
    print("=" * 80)

    if not os.path.exists(TB_RESULTS_FILE):
        raise FileNotFoundError(f"File risultati Tight-Binding non trovato: {TB_RESULTS_FILE}\nEsegui prima la Fase 2.")

    tb_data = np.load(TB_RESULTS_FILE)
    evals = tb_data['evals']    # Autovalori E_n
    evecs = tb_data['evecs']    # Autovettori |phi_n>
    print(f"[OK] Caricati autostati H_TB: {evecs.shape[1]} autostati, range E = [{evals.min():.3f}, {evals.max():.3f}] eV")

    # 1. Inizializzazione dello Stato |psi(0)>
    print(f"\nGenerazione pacchetto d'onde iniziale (z0={Z_0}, p0={P_0}, sigma={SIGMA_Z}, kz0={K_Z0})...")
    psi0 = create_gaussian_wavepacket(z0=Z_0, p0=P_0, sigma_z=SIGMA_Z, kz0=K_Z0, n_proto=N_PROTO, n_rings=N_RINGS)

    # Proiezione di |psi(0)> sulla base degli autostati: c_n = <phi_n | psi(0)>
    c_n = evecs.T.conj() @ psi0
    print(f" Proiezione completata. Norma iniziale conservata: {np.sum(np.abs(c_n)**2):.6f}")

    # 2. Propagazione Temporale
    time_grid = np.linspace(0, T_MAX, N_TSTEPS)
    z_coords = np.arange(N_RINGS)
    
    # Matrici per registrare l'evoluzione
    prob_z_t = np.zeros((N_RINGS, N_TSTEPS))       # P(z, t)
    mean_z_t = np.zeros(N_TSTEPS)                  # <Z(t)>
    msd_t = np.zeros(N_TSTEPS)                     # <Delta Z^2(t)>
    survival_prob = np.zeros(N_TSTEPS)             # P_surv(t) = |<psi(0)|psi(t)>|^2

    print("\nSimulazione evoluzione temporale |psi(t)> = exp(-i H t / hbar) |psi(0)>...")

    for it, t in enumerate(time_grid):
        # Fasi temporali exp(-i E_n t / hbar)
        phase_factors = np.exp(-1j * (evals / HBAR) * t)
        c_n_t = c_n * phase_factors
        
        # Ricostruzione nello spazio reale: |psi(t)> = U * c(t)
        psi_t = evecs @ c_n_t
        psi_t_2d = psi_t.reshape(N_PROTO, N_RINGS)
        
        # Densita di probabilita per anello Z: P(z, t) = sum_p |psi(p, z, t)|^2
        p_z = np.sum(np.abs(psi_t_2d)**2, axis=0)
        prob_z_t[:, it] = p_z
        
        # Valori attesi
        z_avg = np.sum(z_coords * p_z)
        z2_avg = np.sum((z_coords**2) * p_z)
        
        mean_z_t[it] = z_avg
        msd_t[it] = z2_avg - (z_avg**2)
        survival_prob[it] = np.abs(np.vdot(psi0, psi_t))**2

    # 3. Analisi del Regime di Trasporto (Fit esponente alpha: MSD ~ t^alpha)
    # alpha = 2 (balistico), alpha = 1 (diffusivo), alpha = 0 (sub-diffusivo/localizzato)
    valid_mask = (time_grid > 2.0) & (time_grid < T_MAX * 0.7)
    log_t = np.log(time_grid[valid_mask])
    log_msd = np.log(msd_t[valid_mask] + 1e-12)
    poly = np.polyfit(log_t, log_msd, 1)
    alpha_exponent = poly[0]

    print("\n" + "=" * 80)
    print("RISULTATI DINAMICA QUANTISTICA E TRASPORTO")
    print("=" * 80)
    print(f" Esponente di Trasporto alpha (MSD ~ t^alpha) : {alpha_exponent:.4f}")
    if alpha_exponent > 1.7:
        print(" Regime identificato: TRASPORTO QUASI-BALISTICO (Coerenza conservata)")
    elif alpha_exponent > 0.8:
        print(" Regime identificato: TRASPORTO DIFFUSIVO")
    else:
        print(" Regime identificato: LOCALIZZAZIONE DI ANDERSON / SUB-DIFFUSIVO")
    print(f" Posizione finale <Z(t_max)>                   : {mean_z_t[-1]:.2f} / {N_RINGS} anelli")

    # 4. Generazione Grafici
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # A. Mappa Spazio-Tempo P(z, t)
    im0 = axes[0, 0].imshow(prob_z_t, aspect="auto", origin="lower", cmap="inferno",
                            extent=[time_grid[0], time_grid[-1], 0, N_RINGS - 1])
    axes[0, 0].plot(time_grid, mean_z_t, color="cyan", lw=2, linestyle="--", label="<Z(t)> Traiettoria")
    axes[0, 0].set_title("Propagazione del Pacchetto d'Onde P(z, t)")
    axes[0, 0].set_xlabel("Tempo t [fs]")
    axes[0, 0].set_ylabel("Posizione Assiale Z [Anelli]")
    axes[0, 0].legend(loc="upper left")
    plt.colorbar(im0, ax=axes[0, 0], label="Densità di Probabilità")

    # B. Mean Squared Displacement MSD(t)
    axes[0, 1].plot(time_grid, msd_t, color="crimson", lw=2, label=f"MSD(t) [alpha={alpha_exponent:.2f}]")
    axes[0, 1].set_title("Mean Squared Displacement ΔZ²(t)")
    axes[0, 1].set_xlabel("Tempo t [fs]")
    axes[0, 1].set_ylabel("ΔZ² [Anelli²]")
    axes[0, 1].grid(True, linestyle="--", alpha=0.4)
    axes[0, 1].legend()

    # C. Probabilità di Sopravvivenza P_surv(t)
    axes[1, 0].plot(time_grid, survival_prob, color="purple", lw=2)
    axes[1, 0].set_title("Probabilità di Sopravvivenza P_surv(t) = |<ψ(0)|ψ(t)>|²")
    axes[1, 0].set_xlabel("Tempo t [fs]")
    axes[1, 0].set_ylabel("Autosovrapposizione P_surv")
    axes[1, 0].grid(True, linestyle="--", alpha=0.4)

    # D. Istantanee del Profilo Spaziale a tempi t_i
    t_snapshots = [0, N_TSTEPS // 4, N_TSTEPS // 2, N_TSTEPS - 1]
    colors = ["black", "teal", "orange", "red"]
    for idx, st in enumerate(t_snapshots):
        t_val = time_grid[st]
        axes[1, 1].plot(z_coords, prob_z_t[:, st], label=f"t = {t_val:.1f} fs", color=colors[idx], lw=1.8)
    axes[1, 1].set_title("Istantanee del Profilo di Densità P(z)")
    axes[1, 1].set_xlabel("Posizione Assiale Z [Anelli]")
    axes[1, 1].set_ylabel("Densità P(z)")
    axes[1, 1].grid(True, linestyle="--", alpha=0.4)
    axes[1, 1].legend()

    plt.tight_layout()
    chart_path = os.path.join(FIGDIR, "time_dynamics_wavepacket.png")
    plt.savefig(chart_path, dpi=200)
    plt.close()

    # Salvataggio Dati
    np.savez_compressed(
        os.path.join(OUTDIR, "QIP_QC_PHASE3_dynamics_results.npz"),
        time_grid=time_grid,
        prob_z_t=prob_z_t,
        mean_z_t=mean_z_t,
        msd_t=msd_t,
        survival_prob=survival_prob,
        alpha_exponent=alpha_exponent
    )

    print("\n" + "=" * 80)
    print("FASE 3A COMPLETATA CON SUCCESSO")
    print("=" * 80)
    print(f"Grafici e dati salvati in:\n - {chart_path}\n - {os.path.join(OUTDIR, 'QIP_QC_PHASE3_dynamics_results.npz')}")


if __name__ == "__main__":
    main()

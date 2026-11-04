"""
QIP-QC PHASE 9 — HELICAL HARMONIC DECOMPOSITION & ANALYTICAL EXTRACTION
=============================================================================
Extracts the analytical helical harmonics Am,kz and phases from V(theta, z)
to formulate the closed-form potential equation governing the 13-3 lattice.
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


def main():
    print("=" * 80)
    print("QIP-QC PHASE 9 — HELICAL HARMONIC DECOMPOSITION")
    print("=" * 80)

    if not os.path.exists(TB_RESULTS_FILE):
        raise FileNotFoundError(f"File risultati non trovato: {TB_RESULTS_FILE}")

    data = np.load(TB_RESULTS_FILE)
    V = data['V_potential']  # Forma: (N_PROTO=13, N_RINGS=64)

    # 2D Fast Fourier Transform per estrarre i modi (m, k_z)
    V_fft = np.fft.fft2(V)
    amplitudes = np.abs(V_fft) / V.size
    phases = np.angle(V_fft)

    n_proto, n_rings = V.shape
    m_indices = np.fft.fftfreq(n_proto, d=1.0 / n_proto) * n_proto
    z_indices = np.fft.fftfreq(n_rings, d=1.0 / n_rings)

    # Ordinamento dei modi per ampiezza decrescente
    flat_indices = np.argsort(amplitudes.ravel())[::-1]

    print("\n[+] Top 5 Armoniche Elicoidali Dominanti del Potenziale:")
    top_harmonics = []
    for idx in flat_indices[:5]:
        p_idx, z_idx = np.unravel_index(idx, amplitudes.shape)
        m = int(round(m_indices[p_idx]))
        kz = z_indices[z_idx] * 2.0 * np.pi
        amp = amplitudes[p_idx, z_idx]
        phase = phases[p_idx, z_idx]
        top_harmonics.append((m, kz, amp, phase))
        print(f"    m = {m:2d}, k_z = {kz:6.3f} rad/ring, Ampiezza = {amp:.4f}, Fase = {phase:.3f}")

    # Salvataggio dei risultati analitici
    np.savez_compressed(
        os.path.join(OUTDIR, "QIP_QC_PHASE9_harmonics.npz"),
        amplitudes=amplitudes, phases=phases, top_harmonics=top_harmonics
    )

    # Generazione mappa spettrale delle armoniche
    plt.figure(figsize=(10, 6))
    plt.imshow(np.fft.fftshift(amplitudes, axes=(0,)), aspect='auto', cmap='plasma', origin='lower')
    plt.colorbar(label="Ampiezza del Modo $A_{m, k_z}$")
    plt.title("Spettro delle Armoniche Elicoidali del Potenziale $V(\\theta, z)$")
    plt.xlabel("Indice Frequenza Assiale ($k_z$)")
    plt.ylabel("Indice Modulo Azimutale ($m$ - shiftato)")
    plt.tight_layout()

    chart_path = os.path.join(FIGDIR, "phase9_helical_harmonics.png")
    plt.savefig(chart_path, dpi=200)
    plt.close()

    print("\n" + "=" * 80)
    print("FASE 9 (DECOMPOSIZIONE ELICOIDALE) COMPLETATA CON SUCCESSO")
    print("=" * 80)
    print(f" Mappa delle armoniche salvata in : {chart_path}")


if __name__ == "__main__":
    main()

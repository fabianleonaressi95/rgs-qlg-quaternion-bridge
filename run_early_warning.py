"""
run_early_warning.py

Il test che conta davvero per il caso d'uso del middleware: non "quanto e'
critico il sistema a un certo p fissato", ma "le metriche di criticita' si
muovono PRIMA che M_Q crolli, dando un margine di preavviso utile?".

Si simula una singola run in cui fault_prob cresce lentamente nel tempo
(rete che degrada in modo realistico, es. per un attacco progressivo o un
guasto hardware che si diffonde), poi si applica il Criticality Engine con
finestra scorrevole e si controlla se chi/autocorrelazione anticipano il
crollo di M_Q oppure si muovono in contemporanea (nel qual caso non danno
nessun vantaggio predittivo rispetto a guardare M_Q direttamente).
"""

import numpy as np
import matplotlib.pyplot as plt
from simulator import ToyConsensusNetwork
from criticality import CriticalityEngine

N_NODES = 200
N_ROUNDS = 1000
WINDOW = 80
STEP = 4
COLLAPSE_THRESHOLD = 0.65  # sotto questa soglia consideriamo il consenso "a rischio"


def ramping_fault_prob(t, n_rounds, p_start=0.0, p_end=0.35):
    """Fault probability che cresce linearmente nel tempo."""
    return p_start + (p_end - p_start) * (t / n_rounds)


def run():
    net = ToyConsensusNetwork(n_nodes=N_NODES, coupling=0.15, seed=42)
    mq_series = np.zeros(N_ROUNDS)
    p_series = np.zeros(N_ROUNDS)

    for t in range(N_ROUNDS):
        p = ramping_fault_prob(t, N_ROUNDS)
        net.step(p)
        mq_series[t] = net.order_parameter()
        p_series[t] = p

    engine = CriticalityEngine(n_nodes=N_NODES, window_size=WINDOW)
    sliding = engine.sliding_metrics(mq_series, step=STEP)

    rounds = [m["round"] for m in sliding]
    chi = [m["chi"] for m in sliding]
    autocorr = [m["autocorr_lag1"] for m in sliding]
    mq_win = [m["M_Q"] for m in sliding]

    # trova il primo round in cui M_Q scende sotto la soglia di rischio
    collapse_round = next((t for t, v in enumerate(mq_series) if v < COLLAPSE_THRESHOLD), None)

    # trova il primo round in cui chi supera una soglia di allerta
    # (definita qui come 1.5x il valore medio di chi nei primi 20 punti = baseline "sano")
    baseline_chi = np.mean(chi[:20])
    chi_alert_threshold = 1.5 * baseline_chi
    chi_alert_round = next((r for r, c in zip(rounds, chi) if c > chi_alert_threshold), None)

    print(f"Baseline chi (sistema sano): {baseline_chi:.4f}")
    print(f"Soglia di allerta chi (1.5x baseline): {chi_alert_threshold:.4f}")
    print(f"Round in cui M_Q scende sotto {COLLAPSE_THRESHOLD}: {collapse_round}")
    print(f"Round in cui chi supera la soglia di allerta: {chi_alert_round}")
    if collapse_round is not None and chi_alert_round is not None:
        margin = collapse_round - chi_alert_round
        print(f"Margine di preavviso: {margin} round "
              f"({'utile' if margin > 0 else 'NESSUN vantaggio predittivo'})")

    # plot
    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)

    axes[0].plot(mq_series, color="#2b6cb0", linewidth=1)
    axes[0].axhline(COLLAPSE_THRESHOLD, color="red", linestyle="--", label="soglia di rischio")
    if collapse_round:
        axes[0].axvline(collapse_round, color="red", linestyle=":")
    axes[0].set_ylabel("M_Q")
    axes[0].set_title("Order parameter M_Q nel tempo (rete che degrada)")
    axes[0].legend()

    axes[1].plot(rounds, chi, color="#c53030", linewidth=1.5)
    axes[1].axhline(chi_alert_threshold, color="black", linestyle="--", label="soglia di allerta chi")
    if chi_alert_round:
        axes[1].axvline(chi_alert_round, color="black", linestyle=":")
    axes[1].set_ylabel("chi")
    axes[1].set_title("Susceptibility chi (finestra scorrevole)")
    axes[1].legend()

    axes[2].plot(rounds, autocorr, color="#805ad5", linewidth=1.5)
    axes[2].set_ylabel("autocorr lag-1")
    axes[2].set_xlabel("round")
    axes[2].set_title("Autocorrelazione (critical slowing down)")

    plt.tight_layout()
    plt.savefig("early_warning_plot.png", dpi=150)
    print("\nSalvato early_warning_plot.png")


if __name__ == "__main__":
    run()

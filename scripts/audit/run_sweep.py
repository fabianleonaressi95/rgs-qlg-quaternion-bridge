"""
run_sweep.py

Esperimento principale dell'MVP: fa variare il parametro di fault p in un
range (0 = rete perfetta, 1 = rumore totale), e per ogni valore misura le
metriche di criticita' del Criticality Engine in condizioni stazionarie.

Obiettivo: capire se esiste una regione di p in cui chi (susceptibility)
ha un picco netto (comportamento critico, "transizione di fase" nel
consenso) oppure se il sistema degrada in modo graduale e lineare (nel
qual caso le metriche di criticita' non aggiungono molto rispetto a
guardare direttamente M_Q, e l'ipotesi del documento originale va rivista).

Uso:
    python run_sweep.py
Produce:
    sweep_results.csv
    sweep_plot.png
"""

import numpy as np
import csv
from simulator import run_simulation
from criticality import CriticalityEngine

N_NODES = 200
N_ROUNDS = 400          # round totali per ogni valore di p
BURN_IN = 100           # round scartati per far stabilizzare il sistema
WINDOW = 200            # finestra usata dal Criticality Engine (post burn-in)
P_VALUES = np.linspace(0.0, 0.6, 25)  # range del parametro di fault
N_SEEDS = 5             # ripetizioni per ridurre il rumore statistico


def run_sweep():
    engine = CriticalityEngine(n_nodes=N_NODES, window_size=WINDOW)
    rows = []

    for p in P_VALUES:
        chi_vals, binder_vals, mq_vals, autocorr_vals = [], [], [], []
        for seed in range(N_SEEDS):
            series = run_simulation(
                n_nodes=N_NODES, fault_prob=p, n_rounds=N_ROUNDS, seed=seed
            )
            steady = series[BURN_IN:]  # scarta transiente iniziale
            metrics = engine.compute(steady[-WINDOW:])
            chi_vals.append(metrics["chi"])
            binder_vals.append(metrics["binder"])
            mq_vals.append(metrics["M_Q"])
            autocorr_vals.append(metrics["autocorr_lag1"])

        rows.append({
            "p": p,
            "M_Q_mean": np.mean(mq_vals),
            "chi_mean": np.mean(chi_vals),
            "chi_std": np.std(chi_vals),
            "binder_mean": np.mean(binder_vals),
            "autocorr_mean": np.mean(autocorr_vals),
        })
        print(f"p={p:.3f}  M_Q={np.mean(mq_vals):.3f}  chi={np.mean(chi_vals):.4f}  "
              f"binder={np.mean(binder_vals):.3f}")

    with open("sweep_results.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    return rows


def plot_results(rows):
    import matplotlib.pyplot as plt

    p = [r["p"] for r in rows]
    mq = [r["M_Q_mean"] for r in rows]
    chi = [r["chi_mean"] for r in rows]
    chi_std = [r["chi_std"] for r in rows]
    binder = [r["binder_mean"] for r in rows]
    autocorr = [r["autocorr_mean"] for r in rows]

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))

    axes[0, 0].plot(p, mq, marker="o", color="#2b6cb0")
    axes[0, 0].set_title("Order Parameter M_Q vs fault_prob")
    axes[0, 0].set_xlabel("fault_prob (p)")
    axes[0, 0].set_ylabel("M_Q")
    axes[0, 0].axhline(0.5, color="gray", linestyle="--", linewidth=1)

    axes[0, 1].errorbar(p, chi, yerr=chi_std, marker="o", color="#c53030")
    axes[0, 1].set_title("Susceptibility chi vs fault_prob")
    axes[0, 1].set_xlabel("fault_prob (p)")
    axes[0, 1].set_ylabel("chi = N * Var(M_Q)")
    p_peak = p[int(np.argmax(chi))]
    axes[0, 1].axvline(p_peak, color="black", linestyle=":", linewidth=1,
                        label=f"picco a p={p_peak:.2f}")
    axes[0, 1].legend()

    axes[1, 0].plot(p, binder, marker="o", color="#2f855a")
    axes[1, 0].set_title("Binder cumulant U vs fault_prob")
    axes[1, 0].set_xlabel("fault_prob (p)")
    axes[1, 0].set_ylabel("U")

    axes[1, 1].plot(p, autocorr, marker="o", color="#805ad5")
    axes[1, 1].set_title("Autocorrelazione lag-1 vs fault_prob")
    axes[1, 1].set_xlabel("fault_prob (p)")
    axes[1, 1].set_ylabel("autocorr")

    plt.tight_layout()
    plt.savefig("sweep_plot.png", dpi=150)
    print("\nSalvato sweep_plot.png e sweep_results.csv")
    print(f"Picco di chi vicino a fault_prob = {p_peak:.3f}")


if __name__ == "__main__":
    rows = run_sweep()
    plot_results(rows)

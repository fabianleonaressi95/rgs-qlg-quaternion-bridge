import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse.csgraph as csgraph

EPS = 1e-12

# ======================================================
# TOPOLOGIE
# ======================================================

def random_network(N, p=0.05):
    rng = np.random.default_rng(1)
    A = (rng.random((N, N)) < p).astype(float)
    A = np.triu(A, 1)
    return A + A.T

def golden_network(N, k=4):
    phi = (1 + np.sqrt(5)) / 2
    A = np.zeros((N, N))
    for i in range(N):
        for n in range(1, k + 1):
            j = (i + int(round(phi * n))) % N
            A[i, j] = 1
            A[j, i] = 1
    return A

def add_shortcuts(A, p=0.03):
    rng = np.random.default_rng(3)
    N = len(A)
    B = A.copy()
    for i in range(N):
        for j in range(i + 1, N):
            if B[i, j] == 0:
                if rng.random() < p:
                    B[i, j] = 1
                    B[j, i] = 1
    return B

# ======================================================
# ENERGIA
# ======================================================

def edge_energy(A, x, mask):
    E = 0.0
    N = len(x)
    for i in range(N):
        if not mask[i]:
            continue
        for j in range(i + 1, N):
            if not mask[j]:
                continue
            if A[i, j]:
                E += (x[i] - x[j])**2
    return 0.5 * E

# ======================================================
# SIMULAZIONE CON ISOLAMENTO E DINAMICA COERENTE v13
# ======================================================

def run_realization(A, p, seed_sim):
    N = len(A)
    rng = np.random.default_rng(seed_sim)
    bad = np.zeros(N, dtype=bool)
    m = int(N * p)
    if m > 0:
        idx = rng.choice(N, m, replace=False)
        bad[idx] = True
    honest = ~bad
    honest_indices = np.where(honest)[0]

    x = np.ones(N) + 0.05 * rng.normal(size=N)
    x_honest_init = x[honest]
    sigma_0 = np.std(x_honest_init) if len(x_honest_init) > 1 else EPS

    # Sotto-grafo e Laplaciano dei soli nodi onesti (Isolamento Byzantine)
    if len(honest_indices) > 1:
        AH = A[np.ix_(honest_indices, honest_indices)]
        DH = np.diag(AH.sum(axis=1))
        LH = DH - AH
    else:
        LH = np.zeros((1, 1))

    eta = 0.05
    for t in range(700):
        if m > 0:
            attack = 50 * np.sign(np.sin(0.15 * t))
            x[bad] = attack

        if len(honest_indices) > 1:
            dx_h = -eta * (LH @ x[honest_indices])
            x[honest_indices] += dx_h

    # Grandezze fisiche v13
    x_honest_final = x[honest]
    sigma_final = np.std(x_honest_final) if len(x_honest_final) > 1 else 0.0
    OH = max(0.0, 1.0 - (sigma_final / (sigma_0 + EPS)))
    OH = min(1.0, OH)

    EH = edge_energy(A, x, honest)
    EB = edge_energy(A, x, bad) if m > 0 else 0.0
    B_val = EB / (EH + EB + EPS)

    # Connettività algebrica robusta del sottografo onesto
    if len(honest_indices) > 1:
        n_components, labels = csgraph.connected_components(AH, directed=False)
        if n_components > 1:
            largest_comp = np.argmax(np.bincount(labels))
            comp_nodes = np.where(labels == largest_comp)[0]
            if len(comp_nodes) > 1:
                AH_sub = AH[np.ix_(comp_nodes, comp_nodes)]
                DH_sub = np.diag(AH_sub.sum(axis=1))
                LH_sub = DH_sub - AH_sub
                eigvals_H = np.sort(np.linalg.eigvalsh(LH_sub))
                l2_H = eigvals_H[1] if len(eigvals_H) > 1 else 0.0
            else:
                l2_H = 0.0
        else:
            eigvals_H = np.sort(np.linalg.eigvalsh(LH))
            l2_H = eigvals_H[1] if len(eigvals_H) > 1 else 0.0
    else:
        l2_H = 0.0

    D_full = np.diag(A.sum(axis=1))
    L_full = D_full - A
    eigvals_0 = np.sort(np.linalg.eigvalsh(L_full))
    l2_0 = eigvals_0[1] if len(eigvals_0) > 1 else EPS
    
    C_top = l2_H / (l2_0 + EPS)
    C_top = min(1.0, max(0.0, C_top))

    MQ = OH * max(0.1, C_top) * (1.0 - B_val)
    return MQ

# ======================================================
# MAIN
# ======================================================

if __name__ == "__main__":
    N = 250
    mc_runs = 20
    networks = {
        "RANDOM": random_network(N),
        "GOLDEN": golden_network(N),
        "GOLDEN-SW": add_shortcuts(golden_network(N))
    }

    pvals = np.arange(0, 0.62, 0.02)

    plt.figure(figsize=(12, 10))

    plt.subplot(2, 1, 1)
    for name, A in networks.items():
        mq_means = []
        print(f"\n--- Topology: {name} ---")
        print("Byz % | Mean MQ")
        print("------------------------")
        for p in pvals:
            mq_vals = [run_realization(A, p, seed_sim=1000 * i + int(p*100)) for i in range(mc_runs)]
            mean_mq = np.mean(mq_vals)
            mq_means.append(mean_mq)
            print(f"{100*p:3.0f}%  | {mean_mq:7.4f}")

        plt.plot(100 * pvals, mq_means, lw=2, marker="o", label=name)

    plt.grid(True)
    plt.xlabel("Byzantine nodes (%)")
    plt.ylabel("Order Parameter M_Q")
    plt.title("Q-PROOF v13: Phase Diagram with Adaptive Isolation")
    plt.legend()

    plt.subplot(2, 1, 2)
    for name, A in networks.items():
        chi_vals = []
        for p in pvals:
            mq_vals = np.array([run_realization(A, p, seed_sim=1000 * i + int(p*100)) for i in range(mc_runs)])
            chi = N * (np.mean(mq_vals**2) - np.mean(mq_vals)**2)
            chi_vals.append(chi)

        max_chi_idx = np.argmax(chi_vals)
        pc_estimated = 100 * pvals[max_chi_idx]
        print(f"\n>>> [{name}] Critical threshold p_c estimated: {pc_estimated:.1f}%")

        plt.plot(100 * pvals, chi_vals, lw=2, marker="s", label=f"{name} (p_c ≈ {pc_estimated:.1f}%)")

    plt.grid(True)
    plt.xlabel("Byzantine nodes (%)")
    plt.ylabel(r"Susceptibility $\chi_Q$")
    plt.title("Q-PROOF v13: Magnetic Susceptibility Analogy")
    plt.legend()

    plt.tight_layout()
    plt.savefig("qproof_v13_phase_transition.png", dpi=300)
    plt.close()

    print("\nSaved: qproof_v13_phase_transition.png")

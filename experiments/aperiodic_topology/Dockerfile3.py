import numpy as np
import scipy.linalg as la

def generate_aperiodic_chain(N, t=1.0, V=1.0, alpha=(np.sqrt(5) - 1) / 2):
    """
    Genera l'Hamiltoniana del sistema aperiodico (Modello di Aubry-André).
    """
    H = np.zeros((N, N))
    for i in range(N):
        H[i, i] = V * np.cos(2 * np.pi * alpha * i)
        if i < N - 1:
            H[i, i+1] = t
            H[i+1, i] = t
    H[0, N-1] = t
    H[N-1, 0] = t
    return H

def compute_mean_gap_ratio(evals):
    """
    Calcola la media del rapporto delle spaziatura adiacenti (r-statistic):
    r_n = min(delta_n, delta_{n+1}) / max(delta_n, delta_{n+1})
    """
    # Ordina gli autovalori
    E = np.sort(evals)
    # Calcola i delta (spaziature)
    deltas = np.diff(E)
    
    r_values = []
    for i in range(len(deltas) - 1):
        d1 = deltas[i]
        d2 = deltas[i+1]
        if max(d1, d2) > 0:
            r = min(d1, d2) / max(d1, d2)
            r_values.append(r)
            
    return np.mean(r_values)

def run_rmt_analysis(N=1000):
    print("=" * 70)
    print(f" ANALISI RMT: STATISTICA DI LEVEL-SPACING (N = {N} nodi)")
    print("=" * 70)
    
    # Usiamo un N più grande per avere una statistica spettrale robusta
    H = generate_aperiodic_chain(N)
    evals, _ = la.eigh(H)
    
    mean_r = compute_mean_gap_ratio(evals)
    
    # Valori di riferimento teorici universali
    r_poisson = 2.0 * np.log(2) - 1.0  # ~0.3863
    r_goe = 0.5359                     # Wigner-Dyson (Gaussian Orthogonal Ensemble)
    
    print(f" <r> Misurato (Sistema Aperiodico) : {mean_r:.4f}")
    print(f" Riferimento Teorico - Poisson     : {r_poisson:.4f} (Integrabile / Regolare)")
    print(f" Riferimento Teorico - GOE         : {r_goe:.4f} (Caotico / Ergodico)")
    print("-" * 70)
    
    # Distanza dai modelli
    dist_poisson = abs(mean_r - r_poisson)
    dist_goe = abs(mean_r - r_goe)
    
    if dist_poisson < dist_goe:
        print(" [ESITO SPETTRALE: REGIME DI TRANSIZIONE / SUB-POISSONIANO]")
        print(" Lo spettro aperiodico mostra caratteristiche vicine ai sistemi")
        print(" non ergodici o a localization critica, giustificando la natura non banale.")
    else:
        print(" [ESITO SPETTRALE: REPLICA STATISTICA GOE]")
        print(" Lo spettro presenta forte repulsione dei livelli, allineandosi")
        print(" alle statistiche delle matrici casuali ortogonali.")
    print("=" * 70 + "\n")

if __name__ == "__main__":
    run_rmt_analysis()

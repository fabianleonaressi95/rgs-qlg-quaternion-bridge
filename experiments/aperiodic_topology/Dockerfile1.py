import numpy as np
import scipy.linalg as la
import scipy.stats as stats

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

def compute_modal_entropy(eigenvectors):
    """
    Calcola l'entropia modale di Shannon per ciascun autostato.
    """
    N, _ = eigenvectors.shape
    entropies = np.zeros(N)
    for n in range(N):
        prob = np.abs(eigenvectors[:, n])**2
        prob = prob[prob > 1e-12]
        entropies[n] = -np.sum(prob * np.log(prob))
    return entropies

def compute_ensemble_phi(eigenvectors, num_partitions=500):
    """
    Calcola Phi(E) come media d'ensemble su un grande numero di partizioni casuali 
    di dimensione N/2, eliminando la distorsione del singolo taglio geometrico.
    """
    N, _ = eigenvectors.shape
    sub_size = N // 2
    phi_matrix = np.zeros((num_partitions, N))
    rng = np.random.default_rng(42)
    
    for p in range(num_partitions):
        # Genera una partizione casuale del sottosistema
        partition = rng.choice(N, size=sub_size, replace=False)
        A_list = list(partition)
        for n in range(N):
            psi = eigenvectors[:, n]
            norm_A = np.sum(np.abs(psi[A_list])**2)
            phi_matrix[p, n] = 4.0 * norm_A * (1.0 - norm_A)
            
    # Media d'ensemble su tutte le partizioni per ciascun autostato
    ensemble_phi = np.mean(phi_matrix, axis=0)
    return ensemble_phi

def run_ensemble_analysis(N=64):
    print("=" * 70)
    print(f" ANALISI D'ENSEMBLE DI PARTIZIONE PER PHI(E) (N = {N} nodi)")
    print("=" * 70)
    
    # 1. Diagonalizzazione del sistema aperiodico
    H = generate_aperiodic_chain(N)
    evals, evecs = la.eigh(H)
    
    # 2. Calcolo metriche di base e Phi mediato
    modal_entropies = compute_modal_entropy(evecs)
    ensemble_phi = compute_ensemble_phi(evecs, num_partitions=500)
    
    # 3. Correlazione statistica finale sull'ensemble
    corr, p_value = stats.pearsonr(modal_entropies, ensemble_phi)
    
    print(f" Numero di partizioni nell'ensemble : 500")
    print(f" Correlazione di Pearson (Ensemble)  : {corr:.4f}")
    print(f" P-value statistico               : {p_value:.4e}")
    print("-" * 70)
    
    if corr < -0.3:
        print(" [ESITO: ANTI-CORRELAZIONE STRUTTURALE CONFERMATA]")
        print(" Mediando sull'ensemble si rimuovono gli artefatti geometrici locali.")
        print(" Il trade-off tra integrazione e delocalizzazione è una proprietà genuina.")
    elif abs(corr) <= 0.3:
        print(" [ESITO: DISaccoppiamento STATISTICO]")
        print(" La correlazione svanisce nell'ensemble: l'effetto precedente era")
        print(" interamente dovuto alla proiezione geometrica del singolo taglio.")
    else:
        print(" [ESITO: CORRELAZIONE POSITIVA RESIDUA]")
        print(" Il sistema mostra un comportamento di integrazione proporzionale")
        print(" alla delocalizzazione globale.")
    print("=" * 70 + "\n")

if __name__ == "__main__":
    run_ensemble_analysis()

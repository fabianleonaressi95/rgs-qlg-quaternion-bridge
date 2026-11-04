import numpy as np
import scipy.linalg as la
import scipy.stats as stats

def generate_aperiodic_chain(N, t=1.0, V=1.0, alpha=(np.sqrt(5) - 1) / 2):
    """
    Genera l'Hamiltoniana di un sistema aperiodico (es. modello di Aubry-André 
    con frequenza aurea) per simulare lo spettro quantistico.
    """
    H = np.zeros((N, N))
    for i in range(N):
        H[i, i] = V * np.cos(2 * np.pi * alpha * i)
        if i < N - 1:
            H[i, i+1] = t
            H[i+1, i] = t
    # Condizioni al contorno periodiche
    H[0, N-1] = t
    H[N-1, 0] = t
    return H

def compute_modal_entropy(eigenvectors):
    """
    Calcola l'entropia modale di Shannon per ciascun autostato 
    misurando il grado di delocalizzazione nello spazio reale.
    """
    N, _ = eigenvectors.shape
    entropies = np.zeros(N)
    for n in range(N):
        prob = np.abs(eigenvectors[:, n])**2
        prob = prob[prob > 1e-12]  # Filtro per evitare problemi con log(0)
        entropies[n] = -np.sum(prob * np.log(prob))
    return entropies

def compute_phi_proxy(eigenvectors, partition_indices):
    """
    Calcola una metrica proxy di integrazione Phi(E) per una specifica 
    partizione dei sottosistemi (A e B), basata sul bilanciamento di norma.
    """
    N, _ = eigenvectors.shape
    phi_values = np.zeros(N)
    A_set = set(partition_indices)
    A_list = list(A_set)
    
    for n in range(N):
        psi = eigenvectors[:, n]
        norm_A = np.sum(np.abs(psi[A_list])**2)
        # Funzione di integrazione simmetrica (massima quando il peso è ripartito 50/50)
        phi_values[n] = 4.0 * norm_A * (1.0 - norm_A)
        
    return phi_values

def run_stress_test(N=64, num_partitions=100):
    print(f"=" * 70)
    print(f" STRESS TEST DI ROBUSTEZZA: ANALISI DI PARTIZIONE (N = {N} nodi)")
    print(f"=" * 70)
    
    # 1. Generazione del sistema e diagonalizzazione
    H = generate_aperiodic_chain(N)
    evals, evecs = la.eigh(H)
    
    # 2. Calcolo dell'entropia modale per tutti gli autostati
    modal_entropies = compute_modal_entropy(evecs)
    
    correlations = []
    sub_size = N // 2
    rng = np.random.default_rng(42)
    
    # 3. Test iterativo su partizioni geometriche differenti/casuali
    for p in range(num_partitions):
        # Genera una partizione casuale del sottosistema di dimensione N/2
        partition = rng.choice(N, size=sub_size, replace=False)
        phi_vals = compute_phi_proxy(evecs, partition)
        
        # Calcolo della correlazione di Pearson tra Entropia Modale e Phi(E)
        corr, _ = stats.pearsonr(modal_entropies, phi_vals)
        correlations.append(corr)
        
    correlations = np.array(correlations)
    
    # 4. Statistiche aggregate dei risultati dello stress test
    mean_corr = np.mean(correlations)
    std_corr = np.std(correlations)
    min_corr = np.min(correlations)
    max_corr = np.max(correlations)
    
    print(f" Partizioni casuali testate   : {num_partitions}")
    print(f" Correlazione di Pearson media: {mean_corr:.4f}")
    print(f" Deviazione standard          : {std_corr:.4f}")
    print(f" Intervallo Min / Max         : [{min_corr:.4f}, {max_corr:.4f}]")
    print(f"-" * 70)
    
    # Valutazione di robustezza
    all_negative = np.all(correlations < 0)
    if all_negative and mean_corr < -0.5:
        print(" [ESITO: ROBUSTO]")
        print(" L'anti-correlazione è invariante rispetto alla scelta del partizionamento.")
        print(" Il fenomeno NON è un artefatto geometrico della metrica.")
    elif all_negative:
        print(" [ESITO: PARZIALMENTE ROBUSTO]")
        print(" La correlazione rimane negativa ma mostra fluttuazioni legate alla geometria.")
    else:
        print(" [ESITO: FRAGILE / ARTEFATTO]")
        print(" Il segno della correlazione cambia a seconda del taglio di partizione.")
        print(" La Sezione 8 del documento richiede una revisione metodologica.")
    print(f"=" * 70 + "\n")

if __name__ == "__main__":
    run_stress_test()

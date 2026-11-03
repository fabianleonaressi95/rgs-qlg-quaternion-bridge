import numpy as np
import scipy.linalg as la

def generate_periodic_chain(N, t=1.0):
    """
    Genera l'Hamiltoniana di un sistema periodico ideale con condizioni al contorno periodiche.
    """
    H = np.zeros((N, N))
    for i in range(N):
        if i < N - 1:
            H[i, i+1] = t
            H[i+1, i] = t
    H[0, N-1] = t
    H[N-1, 0] = t
    return H

def compute_ensemble_phi(eigenvectors, num_partitions=200):
    """
    Calcola Phi(E) tramite media d'ensemble su partizioni casuali.
    """
    N, _ = eigenvectors.shape
    sub_size = N // 2
    phi_matrix = np.zeros((num_partitions, N))
    rng = np.random.default_rng(42)
    
    for p in range(num_partitions):
        partition = rng.choice(N, size=sub_size, replace=False)
        A_list = list(partition)
        for n in range(N):
            psi = eigenvectors[:, n]
            norm_A = np.sum(np.abs(psi[A_list])**2)
            phi_matrix[p, n] = 4.0 * norm_A * (1.0 - norm_A)
            
    return np.mean(phi_matrix, axis=0)

def run_symmetry_breaking_test(N=64, epsilon=1e-3):
    print("=" * 70)
    print(f" TEST DI ROTTURA DI SIMMETRIA: SISTEMA PERIODICO (N = {N} nodi)")
    print("=" * 70)
    
    # 1. Sistema Periodico Ideale (Onde piane estese)
    H_pure = generate_periodic_chain(N)
    evals_pure, evecs_pure = la.eigh(H_pure)
    phi_pure = compute_ensemble_phi(evecs_pure)
    
    # 2. Sistema Periodico con Rottura di Simmetria (Disordine di Anderson ultradebole)
    rng_disorder = np.random.default_rng(1337)
    disorder = rng_disorder.uniform(-epsilon, epsilon, N)
    H_perturbed = H_pure + np.diag(disorder)
    evals_pert, evecs_pert = la.eigh(H_perturbed)
    phi_pert = compute_ensemble_phi(evecs_pert)
    
    # 3. Statistiche comparative
    mean_pure = np.mean(phi_pure)
    std_pure = np.std(phi_pure)
    
    mean_pert = np.mean(phi_pert)
    std_pert = np.std(phi_pert)
    
    print(f" Ampiezza della perturbazione (ε)        : {epsilon}")
    print(f" Phi medio (Periodico Ideale)            : {mean_pure:.4f} (std: {std_pure:.4f})")
    print(f" Phi medio (Periodico Perturbato)        : {mean_pert:.4f} (std: {std_pert:.4f})")
    print("-" * 70)
    
    # Valutazione
    if abs(mean_pure - mean_pert) < 0.05 and std_pert < 0.1:
        print(" [ESITO: ROBUSTEZZA TOPOLOGICA CONFERMATA]")
        print(" L'elevata integrazione nel sistema periodico NON è un falso positivo")
        print(" dovuto alle degenerazioni. Anche spezzando la simmetria, Phi resta massimo")
        print(" perché riflette la natura totalmente estesa (onde piane) degli stati.")
    else:
        print(" [ESITO: CONFERMA IPOTESI DEL REVISORE]")
        print(" L'elevata integrazione crolla o si disperde introducendo il disordine.")
        print(" Nel sistema periodico puro era effettivamente un artefatto delle degenerazioni.")
    print("=" * 70 + "\n")

if __name__ == "__main__":
    run_symmetry_breaking_test()

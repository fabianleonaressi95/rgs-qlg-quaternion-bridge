import numpy as np
import scipy.linalg as la
import scipy.stats as stats

def build_13_3_hamiltonian(Np=13, Nr=64, t_par=1.0, t_perp=0.5, V_amp=1.0, alpha=(np.sqrt(5) - 1)/2):
    """
    Costruisce l'Hamiltoniana di tight-binding per il reticolo cilindrico 13-3:
    - Np = 13 protofilamenti
    - Nr = 64 anelli per protofilamento (Totale N = 832 siti)
    - Shift elicoidale = 3
    - Accoppiamento longitudinale (t_par) e trasversale/elicoidale (t_perp)
    - Potenziale aperiodico di sito modulato lungo i protofilamenti
    """
    N = Np * Nr
    H = np.zeros((N, N))
    
    # Mappatura indice (p, r) -> i lineare
    def idx(p, r):
        return p * Nr + (r % Nr)
    
    # Riempimento dei termini di hopping e potenziale di sito
    for p in range(Np):
        for r in range(Nr):
            i = idx(p, r)
            
            # Potenziale di sito aperiodico (es. modulazione aurea sui protofilamenti)
            H[i, i] = V_amp * np.cos(2 * np.pi * alpha * p)
            
            # 1. Hopping longitudinale (lungo il protofilamento p)
            j_down = idx(p, r + 1)
            H[i, j_down] = t_par
            H[j_down, i] = t_par
            
            # 2. Hopping trasversale / elicoidale (verso protofilamento p+1 e p-1 con shift r+3)
            for dp in [1, -1]:
                p_neigh = (p + dp) % Np
                r_neigh = r + 3
                j_helical = idx(p_neigh, r_neigh)
                H[i, j_helical] = t_perp
                H[j_helical, i] = t_perp
                
    return H

def compute_physical_modal_entropy_and_phi(eigenvectors, Np=13, Nr=64):
    """
    Calcola l'entropia modale e la metrica Phi usando la partizione fisica 
    strutturale (es. Subsystem A = primi 6 protofilamenti, B = restanti 7).
    """
    N, _ = eigenvectors.shape
    entropies = np.zeros(N)
    phi_values = np.zeros(N)
    
    # Definizione partizione fisica: metà protofilamenti (es. p da 0 a 5)
    A_sites = []
    for p in range(Np // 2):
        for r in range(Nr):
            A_sites.append(p * Nr + r)
    A_set = set(A_sites)
    
    for n in range(N):
        psi = eigenvectors[:, n]
        prob = np.abs(psi)**2
        
        # Entropia modale di Shannon
        p_filt = prob[prob > 1e-12]
        entropies[n] = -np.sum(p_filt * np.log(p_filt))
        
        # Phi basato sul peso di norma nella partizione fisica dei protofilamenti
        norm_A = np.sum(prob[list(A_set)])
        phi_values[n] = 4.0 * norm_A * (1.0 - norm_A)
        
    return entropies, phi_values

def run_epsilon_scan(H_base, Np=13, Nr=64):
    """
    Esegue uno scan parametrico della rottura di simmetria su più ordini di grandezza di epsilon.
    """
    epsilon_values = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1]
    print("\n" + "=" * 70)
    print(" 2. SCAN PARAMETRICO DI ROTTURA DI SIMMETRIA (Reticolo 13-3)")
    print("=" * 70)
    print(f"{'Epsilon (ε)':<15} | {'Phi Medio':<15} | {'Std Dev Phi':<15}")
    print("-" * 50)
    
    rng = np.random.default_rng(42)
    N = H_base.shape[0]
    
    for eps in epsilon_values:
        # Aggiunge disordine di Anderson di ampiezza eps
        disorder = rng.uniform(-eps, eps, N)
        H_pert = H_base + np.diag(disorder)
        evals, evecs = la.eigh(H_pert)
        _, phi_vals = compute_physical_modal_entropy_and_phi(evecs, Np, Nr)
        print(f"{eps:<15.1e} | {np.mean(phi_vals):<15.4f} | {np.std(phi_vals):<15.4f}")
    print("=" * 70)

def compute_central_rmt_r(evals, fraction=0.4):
    """
    Calcola l'r-statistic di Wigner-Dyson isolando la regione centrale dello spettro 
    per evitare distorsioni da gap globali di banda del reticolo 13-3.
    """
    E = np.sort(evals)
    N = len(E)
    start_idx = int(N * (0.5 - fraction / 2))
    end_idx = int(N * (0.5 + fraction / 2))
    E_center = E[start_idx:end_idx]
    
    deltas = np.diff(E_center)
    r_values = []
    for i in range(len(deltas) - 1):
        d1, d2 = deltas[i], deltas[i+1]
        if max(d1, d2) > 1e-12:
            r = min(d1, d2) / max(d1, d2)
            if r <= 1.0:
                r_values.append(r)
    return np.mean(r_values) if r_values else 0.0

def run_full_validation():
    Np, Nr = 13, 64
    N = Np * Nr
    print("=" * 70)
    print(f" SUITE DI VALIDAZIONE RIGOROSA: RETICOLO 13-3 (N = {N} siti)")
    print("=" * 70)
    
    # Generazione Hamiltoniana 13-3
    H_13_3 = build_13_3_hamiltonian(Np, Nr)
    evals, evecs = la.eigh(H_13_3)
    
    # 1. Test correlazione con partizione fisica strutturale
    entropies, phi_vals = compute_physical_modal_entropy_and_phi(evecs, Np, Nr)
    corr, p_val = stats.pearsonr(entropies, phi_vals)
    
    print("\n 1. ANALISI DI CORRELAZIONE CON PARTIZIONE FISICA (Protofilamenti)")
    print(f" Correlazione di Pearson (Fisica) : {corr:.4f} (p-value: {p_val:.4e})")
    if corr < -0.2:
        print(" [ESITO: Anti-correlazione strutturale genuina confermata col taglio fisico]")
    elif abs(corr) <= 0.2:
        print(" [ESITO: Disaccoppiamento tra integrazione e delocalizzazione spaziale]")
    else:
        print(" [ESITO: Correlazione positiva mantenuta anche con partizione geometrica rigida]")
        
    # 2. Scan di rottura di simmetria
    run_epsilon_scan(H_13_3, Np, Nr)
    
    # 3. Analisi RMT r-statistic centrale
    mean_r = compute_central_rmt_r(evals, fraction=0.5)
    r_poisson = 2.0 * np.log(2) - 1.0  # ~0.3863
    r_goe = 0.5359                     # Wigner-Dyson
    
    print("\n 3. STATISTICA RMT (LEVEL-SPACING CENTRALE)")
    print(f" <r> Misurato (Reticolo 13-3)    : {mean_r:.4f}")
    print(f" Riferimento - Poisson           : {r_poisson:.4f}")
    print(f" Riferimento - GOE               : {r_goe:.4f}")
    print("-" * 70)
    if abs(mean_r - r_poisson) < 0.08:
        print(" [ESITO RMT: REGIME DI LOCALIZZAZIONE / CRITICO]")
    elif mean_r > 0.45:
        print(" [ESITO RMT: REGIME ERGODICO / GOE]")
    else:
        print(" [ESITO RMT: SPETTRO MULTIFRATTALE / COMPLESSO APERIODICO]")
    print("=" * 70 + "\n")

if __name__ == "__main__":
    run_full_validation()

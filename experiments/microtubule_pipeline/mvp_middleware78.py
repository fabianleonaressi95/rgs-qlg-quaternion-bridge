import numpy as np
import scipy.linalg as la
import scipy.stats as stats

def build_13_3_matrices_in_k_space(Np=13, Nr=64, t_par=1.0, t_perp=0.5, V_amp=1.0, alpha=(np.sqrt(5) - 1)/2):
    """
    Costruisce l'Hamiltoniana in k-spazio (blocchi 13x13).
    """
    k_values = 2.0 * np.pi * np.arange(Nr) / Nr
    H_blocks = []
    
    for k in k_values:
        H_k = np.zeros((Np, Np), dtype=complex)
        for p in range(Np):
            H_k[p, p] = V_amp * np.cos(2 * np.pi * alpha * p)
            H_k[p, p] += 2.0 * t_par * np.cos(k)
            for dp in [1, -1]:
                p_neigh = (p + dp) % Np
                H_k[p, p_neigh] += t_perp * np.exp(1j * 3 * k)
                H_k[p_neigh, p] += t_perp * np.exp(-1j * 3 * k)
        H_blocks.append(H_k)
    return H_blocks

def compute_section_7_metrics_from_evecs(evecs, Np=13, Nr=64):
    """
    Calcola S_vN e Phi direttamente dagli autovettori già forniti (senza ricalcolare la diagonalizzazione).
    """
    N, _ = evecs.shape
    s_vn_list = []
    phi_list = []
    
    for n in range(N):
        psi = evecs[:, n]
        psi_matrix = psi.reshape(Np, Nr)
        
        # Matrice densità ridotta sui protofilamenti
        rho_p = psi_matrix @ psi_matrix.conj().T
        trace_rho = np.trace(rho_p).real
        if trace_rho > 1e-12:
            rho_p /= trace_rho
            
        lambda_p = la.eigvalsh(rho_p)
        lambda_p = np.clip(lambda_p, 1e-15, 1.0)
        lambda_p /= np.sum(lambda_p)
        
        s_vn = -np.sum(lambda_p * np.log(lambda_p))
        s_vn_list.append(s_vn)
        
        purity = np.sum(lambda_p**2)
        phi_val = 4.0 * purity * (1.0 - purity)
        phi_list.append(phi_val)
        
    return np.array(s_vn_list), np.array(phi_list)

def run_ensemble_epsilon_scan_fast(H_base, Np=13, Nr=64, num_realizations=20):
    """
    Scan parametrico ottimizzato (evita ridondanze di diagonalizzazione).
    """
    epsilon_values = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1]
    N = Np * Nr
    
    print("\n" + "=" * 70)
    print(f" 2. SCAN DI ROTTURA DI SIMMETRIA OTTIMIZZATO (N_real = {num_realizations})")
    print("=" * 70)
    print(f"{'Epsilon (ε)':<15} | {'Phi Medio':<15} | {'Std Dev Phi':<15} | {'Max ||H_pert - H_base||_inf':<25}")
    print("-" * 80)
    
    for eps in epsilon_values:
        phi_means = []
        max_diffs = []
        
        for seed in range(num_realizations):
            rng = np.random.default_rng(42 + seed)
            disorder = rng.uniform(-eps, eps, N)
            H_pert = H_base + np.diag(disorder)
            
            diff_norm = np.max(np.abs(H_pert - H_base))
            max_diffs.append(diff_norm)
            
            # UNICA diagonalizzazione per questa iterazione
            evals, evecs = la.eigh(H_pert)
            _, phi_vals = compute_section_7_metrics_from_evecs(evecs, Np, Nr)
            phi_means.append(np.mean(phi_vals))
            
        print(f"{eps:<15.1e} | {np.mean(phi_means):<15.4f} | {np.std(phi_means):<15.4f} | {np.mean(max_diffs):<25.4e}")
    print("=" * 70)

def compute_block_rmt_r(H_blocks):
    r_sector_means = []
    for H_k in H_blocks:
        evals_k = la.eigvalsh(H_k)
        evals_k = np.sort(evals_k.real)
        deltas = np.diff(evals_k)
        
        r_values = []
        for i in range(len(deltas) - 1):
            d1, d2 = deltas[i], deltas[i+1]
            if max(d1, d2) > 1e-12:
                r = min(d1, d2) / max(d1, d2)
                if r <= 1.0:
                    r_values.append(r)
        if r_values:
            r_sector_means.append(np.mean(r_values))
    return np.mean(r_sector_means) if r_sector_means else 0.0

def run_fast_validation_suite():
    Np, Nr = 13, 64
    N = Np * Nr
    print("=" * 70)
    print(f" SUITE OTTIMIZZATA: RETICOLO 13-3 (N = {N} siti)")
    print("=" * 70)
    
    # Generazione H spaziale di base
    H_base = np.zeros((N, N))
    def idx(p, r): return p * Nr + (r % Nr)
    alpha = (np.sqrt(5) - 1) / 2
    for p in range(Np):
        for r in range(Nr):
            i = idx(p, r)
            H_base[i, i] = np.cos(2 * np.pi * alpha * p)
            H_base[i, idx(p, r+1)] = 1.0
            H_base[idx(p, r+1), i] = 1.0
            for dp in [1, -1]:
                p_n = (p + dp) % Np
                j_h = idx(p_n, r+3)
                H_base[i, j_h] = 0.5
                H_base[j_h, i] = 0.5
                
    # 1. Test correlazione Sezione 7
    print("\n 1. ANALISI DI CORRELAZIONE SEZIONE 7 (Tr_z, S_vN)")
    evals_base, evecs_base = la.eigh(H_base)
    s_vn, phi_vals = compute_section_7_metrics_from_evecs(evecs_base, Np, Nr)
    corr, p_val = stats.pearsonr(s_vn, phi_vals)
    print(f" Correlazione di Pearson (S_vN vs Phi) : {corr:.4f} (p-value: {p_val:.4e})")
    
    # 2. Scan di simmetria ottimizzato
    run_ensemble_epsilon_scan_fast(H_base, Np, Nr, num_realizations=20)
    
    # 3. RMT pulita da blocchi k-spazio
    H_blocks = build_13_3_matrices_in_k_space(Np, Nr)
    mean_r_block = compute_block_rmt_r(H_blocks)
    r_poisson = 2.0 * np.log(2) - 1.0
    r_goe = 0.5359
    
    print("\n 3. STATISTICA RMT PULITA (BLOCCHI k-SPAZIO 13x13)")
    print(f" <r> Medio per Settore k         : {mean_r_block:.4f}")
    print(f" Riferimento - Poisson           : {r_poisson:.4f}")
    print(f" Riferimento - GOE               : {r_goe:.4f}")
    print("=" * 70 + "\n")

if __name__ == "__main__":
    run_fast_validation_suite()

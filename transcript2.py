import numpy as np
import scipy.linalg as LA

def projected_6d_lie_algebra_audit():
    print("=" * 80)
    print("   RGS-Q-HID: GENERAZIONE ALGEBRA DI LIE DA PROIEZIONE 6D -> 3D")
    print("=" * 80)
    
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    
    # 1. Matrice di proiezione 3x6 basata sulla simmetria aperiodica (aurea)
    Pi_parallel = np.array([
        [1.0, phi, 0.0, -phi, 1.0, 0.0],
        [0.0, 1.0, phi, 0.0, -phi, 1.0],
        [phi, 0.0, 1.0, 1.0, 0.0, -phi]
    ], dtype=float)
    
    # Normalizzazione delle righe della proiezione
    for i in range(3):
        Pi_parallel[i] /= LA.norm(Pi_parallel[i])
        
    # 2. Estrazione dei vettori di base proiettati in 3D dai 6 assi del reticolo superiore
    base_vectors_3d = [Pi_parallel[:, k] for k in range(6)]
    
    # 3. Costruzione dei generatori 3x3 tramite prodotti esterni antisimmetrici (Hermitiani via 1j)
    generators_3d = []
    for i in range(6):
        for j in range(i + 1, 6):
            v_i = base_vectors_3d[i]
            v_j = base_vectors_3d[j]
            m3 = np.outer(v_i, v_j) - np.outer(v_j, v_i)
            if LA.norm(m3) > 1e-5:
                generators_3d.append(1j * m3)
                
    # 4. Estrazione della base ortonormale indipendente tramite SVD
    flat_gens = [g.flatten() for g in generators_3d]
    matrix_stack = np.array(flat_gens)
    
    U, s, Vh = LA.svd(matrix_stack, full_matrices=False)
    rank = np.sum(s > 1e-5)
    
    print(f" * Numero di generatori grezzi generati: {len(generators_3d)}")
    print(f" * Rango effettivo dello span (indipendenti): {rank}")
    
    independent_gens = []
    for idx in range(rank):
        g_mat = Vh[idx].reshape(3, 3)
        norm_val = LA.norm(g_mat, 'fro')
        if norm_val > 1e-12:
            g_mat /= norm_val
            independent_gens.append(g_mat)
            
    # 5. TEST DI CHIUSURA DELL'ALGEBRA DI LIE SULLA BASE PROIETTATA
    print(f"\n--- TEST DI CHIUSURA DELL'ALGEBRA PROIETTATA (Dim base: {len(independent_gens)}) ---")
    max_residual = 0.0
    N_g = len(independent_gens)
    
    for i in range(N_g):
        for j in range(N_g):
            comm = independent_gens[i] @ independent_gens[j] - independent_gens[j] @ independent_gens[i]
            # Proiezione sullo span ortonormale
            projected = np.zeros_like(comm)
            for k in range(N_g):
                c_k = np.sum(independent_gens[k].conj() * comm)
                projected += c_k * independent_gens[k]
            residual = LA.norm(comm - projected, 'fro')
            if residual > max_residual:
                max_residual = residual
                
    print(f" * Massimo residuo di chiusura algebrica (Proiezione 6D): {max_residual:.2e}")
    if max_residual < 1e-6:
        print(" -> SUCCESSO: L'algebra di Lie proiettata si chiude rigorosamente!")
    else:
        print(" -> L'algebra proiettata richiede un'estensione tensoriale o di Kac-Moody.")
        
    print("=" * 80)

if __name__ == "__main__":
    projected_6d_lie_algebra_audit()

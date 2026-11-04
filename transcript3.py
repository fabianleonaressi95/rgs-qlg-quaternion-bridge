import numpy as np
import scipy.linalg as LA

def casimir_fission_audit():
    print("=" * 80)
    print("   RGS-Q-HID: CALCOLO RIPARTIZIONE DI FISSIONE VIA INVARIANTI DI CASIMIR")
    print("=" * 80)
    
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    
    # 1. Proiezione 6D -> 3D
    Pi_parallel = np.array([
        [1.0, phi, 0.0, -phi, 1.0, 0.0],
        [0.0, 1.0, phi, 0.0, -phi, 1.0],
        [phi, 0.0, 1.0, 1.0, 0.0, -phi]
    ], dtype=float)
    
    for i in range(3):
        Pi_parallel[i] /= LA.norm(Pi_parallel[i])
        
    base_vectors_3d = [Pi_parallel[:, k] for k in range(6)]
    
    generators_3d = []
    for i in range(6):
        for j in range(i + 1, 6):
            v_i = base_vectors_3d[i]
            v_j = base_vectors_3d[j]
            m3 = np.outer(v_i, v_j) - np.outer(v_j, v_i)
            if LA.norm(m3) > 1e-5:
                generators_3d.append(1j * m3)
                
    flat_gens = [g.flatten() for g in generators_3d]
    matrix_stack = np.array(flat_gens)
    
    U, s, Vh = LA.svd(matrix_stack, full_matrices=False)
    rank = np.sum(s > 1e-5)
    
    independent_gens = []
    for idx in range(rank):
        g_mat = Vh[idx].reshape(3, 3)
        g_mat /= LA.norm(g_mat, 'fro')
        independent_gens.append(g_mat)
        
    # 2. Costruzione dell'operatore quadratico di Casimir C_2 = \sum_a G_a^\dagger G_a
    C2 = np.zeros((3, 3), dtype=complex)
    for G_a in independent_gens:
        C2 += G_a.conj().T @ G_a
        
    evals_C2 = np.real(LA.eigvalsh(C2))
    print(f" * Autovalori dell'operatore di Casimir C_2: {evals_C2}")
    
    # 3. Mappa di asimmetria parametro-free basata sui livelli dello spaccamento di Casimir
    # L'asimmetria topologica emerge dal rapporto tra lo spaccamento massimo e la traccia di C_2
    casimir_asymmetry = (evals_C2[-1] - evals_C2[0]) / np.trace(C2)
    print(f" * Fattore di asimmetria topologica (Casimir): {casimir_asymmetry:.6f}")
    
    # 4. Applicazione ai nuclei composti di Pu-239 (A=240) e U-235 (A=236)
    targets = [
        ("Plutonio-239", 240.0, 2.87),
        ("Uranio-235", 236.0, 2.42)
    ]
    
    for name, A_comp, prompt_exp in targets:
        # La ripartizione del nucleo A_comp segue la frazione d'angolo aureo scalata dall'asimmetria di Casimir
        asym_ratio = phi**-1 + casimir_asymmetry
        
        A_H = A_comp / (1.0 + (1.0 / (1.0 + asym_ratio)))
        A_L = A_comp - A_H
        
        print(f"\n--- VALUTAZIONE CASIMIR PER {name} (Nucleo composto A = {A_comp:.0f}) ---")
        print(f" * Picco leggero derivato (A_L):   {A_L:.2f}")
        print(f" * Picco pesante derivato (A_H):   {A_H:.2f}")
        print(f" * Rapporto di massa (A_H / A_L): {A_H / A_L:.4f}")

    print("=" * 80)

if __name__ == "__main__":
    casimir_fission_audit()

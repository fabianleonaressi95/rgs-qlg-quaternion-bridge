import numpy as np
import scipy.linalg as LA

def perturbative_casimir_fission():
    print("=" * 80)
    print("   RGS-Q-HID: PERTURBAZIONE AL PRIMO ORDINE VIA CURVATURA YANG-MILLS")
    print("=" * 80)
    
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    vol_active = 2.0 / 3.0
    vol_sterile = 1.0 / 3.0
    
    # 1. Kernel geometrico e campo di gauge
    T = np.array([
        [phi**-1, phi**-1 * np.exp(1j * np.pi / phi), phi**-2],
        [phi**-1 * np.exp(-1j * np.pi / phi), 1.0, phi**-1],
        [phi**-2, phi**-1, phi**-1]
    ], dtype=complex)
    T = (T + T.conj().T) / 2.0  
    
    V_biverso = np.diag([vol_active, vol_active, vol_sterile]).astype(complex) * phi**-2
    V_biverso = (V_biverso + V_biverso.conj().T) / 2.0
    
    F_comm = T @ V_biverso - V_biverso @ T
    F_QCD = 1j * F_comm  
    F_QCD = (F_QCD + F_QCD.conj().T) / 2.0
    norm_f = LA.norm(F_QCD, 'fro')
    
    # 2. Proiezione 6D -> 3D per l'algebra di Lie chiusa
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
        
    # 3. Casimir operator e perturbazione al primo ordine
    C2 = np.zeros((3, 3), dtype=complex)
    for G_a in independent_gens:
        C2 += G_a.conj().T @ G_a
        
    trace_C2 = np.real(np.trace(C2))
    
    # Parametro di perturbazione geometrico puro
    perturbation_param = (norm_f / trace_C2) * (phi**-2)
    
    # Rapporto di massa perturbato ab initio
    perturbed_ratio = phi * (1.0 - perturbation_param)
    
    targets = [
        ("Plutonio-239", 240.0, 2.87),
        ("Uranio-235", 236.0, 2.42)
    ]
    
    for name, A_comp, exp_nu in targets:
        A_H = A_comp / (1.0 + (1.0 / perturbed_ratio))
        A_L = A_comp - A_H
        derived_neutrons = A_comp - (A_L + A_H) # Bilancio di massa primario
        
        print(f"\n--- PREVISIONE PERTURBATA PER {name} (A = {A_comp:.0f}) ---")
        print(f" * Parametro perturbativo YM ($\epsilon$): {perturbation_param:.6f}")
        print(f" * Rapporto di massa effettivo ($A_H/A_L$): {perturbed_ratio:.4f}")
        print(f" * Picco leggero derivato ($A_L$):           {A_L:.2f}")
        print(f" * Picco pesante derivato ($A_H$):           {A_H:.2f}")
        print(f" * Residuo di massa (Neutroni/Frammenti):  {derived_neutrons:.2f}")

    print("=" * 80)

if __name__ == "__main__":
    perturbative_casimir_fission()

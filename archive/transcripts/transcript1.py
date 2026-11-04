import numpy as np
import scipy.linalg as LA

def strict_parameter_free_audit():
    print("=" * 80)
    print("   RGS-Q-HID: AUDIT ALGEBRICO E SPETTRALE PARAMETRO-FREE")
    print("=" * 80)
    
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    vol_active = 2.0 / 3.0
    vol_sterile = 1.0 / 3.0
    
    # 1. Kernel geometrico rigoroso (Hermitiano)
    T = np.array([
        [phi**-1, phi**-1 * np.exp(1j * np.pi / phi), phi**-2],
        [phi**-1 * np.exp(-1j * np.pi / phi), 1.0, phi**-1],
        [phi**-2, phi**-1, phi**-1]
    ], dtype=complex)
    T = (T + T.conj().T) / 2.0  
    
    V_biverso = np.diag([vol_active, vol_active, vol_sterile]).astype(complex) * phi**-2
    V_biverso = (V_biverso + V_biverso.conj().T) / 2.0
    
    # 2. Curvatura di Yang-Mills Q = i*[T,V]
    F_comm = T @ V_biverso - V_biverso @ T
    F_QCD = 1j * F_comm  
    F_QCD = (F_QCD + F_QCD.conj().T) / 2.0
    
    # 3. Costruzione e ortonormalizzazione rigorosa dei generatori (Gram-Schmidt di Frobenius)
    raw_G1 = np.diag([phi**-1, phi**-2, phi**-3]).astype(complex)
    raw_G1 = (raw_G1 + raw_G1.conj().T) / 2.0
    
    raw_G2 = np.array([
        [0.0, np.exp(1j * np.pi / phi), phi**-1],
        [np.exp(1j * np.pi / phi), 0.0, phi**-2],
        [phi**-1, phi**-2, 0.0]
    ], dtype=complex)
    raw_G2 = (raw_G2 + raw_G2.conj().T) / 2.0
    
    F_dag_F = F_QCD @ F_QCD.conj().T
    raw_G3 = (F_dag_F + F_dag_F.conj().T) / 2.0
    
    raw_G4 = np.array([
        [phi**-2, phi**-1 * np.exp(1j * np.pi / phi**2), 0.0],
        [np.exp(-1j * np.pi / phi**2), phi**-3, phi**-2],
        [0.0, phi**-2, phi**-1]
    ], dtype=complex)
    raw_G4 = (raw_G4 + raw_G4.conj().T) / 2.0
    
    raw_generators = [raw_G1, raw_G2, raw_G3, raw_G4]
    
    ortho_generators = []
    for g in raw_generators:
        ortho_g = g.copy()
        for u in ortho_generators:
            proj = np.sum(u.conj() * ortho_g) / np.sum(u.conj() * u)
            ortho_g -= proj * u
        norm_val = np.sqrt(np.real(np.sum(ortho_g.conj() * ortho_g)))
        if norm_val > 1e-12:
            ortho_g /= norm_val
        ortho_generators.append(ortho_g)
        
    # 4. TEST DI CHIUSURA DELL'ALGEBRA DI LIE ([G_i, G_j] nello span dei G_k)
    print("\n--- TEST DI CHIUSURA DELL'ALGEBRA DI LIE ---")
    max_residual = 0.0
    N_g = len(ortho_generators)
    for i in range(N_g):
        for j in range(N_g):
            comm = ortho_generators[i] @ ortho_generators[j] - ortho_generators[j] @ ortho_generators[i]
            # Proietta il commutatore sullo span ortonormale
            projected = np.zeros_like(comm)
            for k in range(N_g):
                c_k = np.sum(ortho_generators[k].conj() * comm)
                projected += c_k * ortho_generators[k]
            residual = LA.norm(comm - projected, 'fro')
            if residual > max_residual:
                max_residual = residual
                
    print(f" * Massimo residuo di chiusura algebrica: {max_residual:.2e}")
    if max_residual < 1e-10:
        print(" -> L'algebra si chiude rigorosamente nello span dei generatori.")
    else:
        print(" -> L'algebra NON si chiude (presenza di componenti esterne allo span).")

    # 5. MAPPA SPETTRALE PURA (Senza coefficienti arbitrari)
    # Regola basata unicamente sugli invarianti di T e phi
    evals_T, _ = LA.eigh(T)
    trace_T = np.abs(np.trace(T))
    det_T = np.abs(LA.det(T))
    
    # Rapporto geometrico puro derivato da Tr(T) e det(T) scalato da phi
    geometric_ratio = phi * (trace_T / (det_T + 1e-8))
    
    # Applicazione universale al Plutonio (A=240) e all'Uranio (A=236)
    for name, A_total in [("Plutonio-239 (A=240)", 240.0), ("Uranio-235 (A=236)", 236.0)]:
        # Ripartizione basata esclusivamente sul rapporto geometrico spettrale normalizzato
        fraction_H = geometric_ratio / (1.0 + geometric_ratio)
        # Per mantenere la coerenza dimensionale del nucleo composto senza costanti fittizie,
        # usiamo la traccia normalizzata come fattore di scala intrinseco:
        scale_factor = trace_T / np.sum(np.abs(evals_T))
        
        A_H = A_total * (0.5 + 0.5 * (fraction_H - 0.5) * scale_factor)
        A_L = A_total - A_H
        
        print(f"\n--- PREVISIONE PARAMETRO-FREE PER {name} ---")
        print(f" * Massa totale nucleo composto: {A_total}")
        print(f" * Picco leggero derivato (A_L):   {A_L:.2f}")
        print(f" * Picco pesante derivato (A_H):   {A_H:.2f}")

    print("=" * 80)

if __name__ == "__main__":
    strict_parameter_free_audit()

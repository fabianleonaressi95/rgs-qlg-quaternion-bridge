import numpy as np

def audit_Mc(M, max_power=200):
    M = np.asarray(M, dtype=np.complex128)
    if M.shape != (8, 8):
        raise ValueError(f"Forma inattesa: {M.shape}, attesa (8, 8)")
    
    # 1. Spettro e raggio spettrale
    eig = np.linalg.eigvals(M)
    rho = float(np.max(np.abs(eig)))
    
    # 2. Non-normalità
    comm = M @ M.conj().T - M.conj().T @ M
    norm_M_fro = np.linalg.norm(M, "fro")
    comm_rel = float(np.linalg.norm(comm, "fro") / max(norm_M_fro**2, 1e-15))
    
    # 3. Potenze e amplificazione transiente
    P = np.eye(8, dtype=np.complex128)
    peak = (1.0, 0)
    for n in range(1, max_power + 1):
        P = P @ M
        val = float(np.linalg.norm(P, 2))
        if val > peak[0]:
            peak = (val, n)
            
    # 4. Operatore-Schmidt Rank (OSR)
    try:
        M4 = M.reshape(2, 4, 2, 4)
        reshuffled = M4.transpose(0, 2, 1, 3).reshape(4, 16)
        s = np.linalg.svd(reshuffled, compute_uv=False)
        tol = max(reshuffled.shape) * np.finfo(float).eps * s[0]
        osr = int(np.sum(s > tol))
    except Exception as e:
        osr = None
        s = str(e)

    print("=== RISULTATI AUDIT M_c ===")
    print(f"Raggio spettrale (rho): {rho:.6f}")
    print(f"Non-normalità relativa: {comm_rel:.6f}")
    print(f"Picco di crescita transiente ||M^n||_2: {peak[0]:.6f} al passo {peak[1]}")
    print(f"OSR (fattorizzazione 2x4): {osr}")
    
    if rho < 1.0:
        print("Stabilità asintotica: VERIFICATA (rho < 1)")
    else:
        print("Stabilità asintotica: NON VERIFICATA (rho >= 1)")
        
    if peak[0] > 1.0 + 1e-10:
        print("Amplificazione transiente: OSSERVATA")
    else:
        print("Amplificazione transiente: NON OSSERVATA nell'intervallo")

    return {
        "rho": rho,
        "comm_rel": comm_rel,
        "peak": peak,
        "osr": osr
    }

if __name__ == "__main__":
    # Esecuzione immediata con una matrice di test casuale
    np.random.seed(42)
    M_test = np.random.randn(8, 8) * 0.95
    audit_Mc(M_test)

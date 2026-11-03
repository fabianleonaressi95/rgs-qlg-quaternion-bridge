import numpy as np

def build_and_audit_Mc():
    # 1. Costruzione esatta di una matrice 8x8 non normale con rho < 1 e OSR <= 4
    # Usiamo una struttura a blocchi coerente con la fattorizzazione R^2 \otimes R^4
    np.seed = 42
    
    # Sotto-blocchi di interazione non normali per generare amplificazione transiente e dissipatività
    A = np.array([
        [ 0.65, -0.35,  0.10,  0.05],
        [ 0.40,  0.55, -0.15,  0.00],
        [-0.05,  0.20,  0.60, -0.30],
        [ 0.10, -0.10,  0.35,  0.50]
    ])
    
    # Introduciamo una non-normalità controllata accoppiando i blocchi
    M = np.kron(np.array([[0.8, 0.3], [-0.2, 0.7]]), A)
    
    # Normalizzazione per garantire rigorosamente che il raggio spettrale sia < 1 (es. ~0.973 come da specifiche)
    eigvals = np.linalg.eigvals(M)
    rho_curr = np.max(np.abs(eigvals))
    target_rho = 0.973373
    M = M * (target_rho / rho_curr)
    
    # 2. Audit della matrice M
    M = np.asarray(M, dtype=np.complex128)
    assert M.shape == (8, 8), f"Forma inattesa: {M.shape}"
    
    eig = np.linalg.eigvals(M)
    rho = float(np.max(np.abs(eig)))
    
    Mstar = M.conj().T
    comm = M @ Mstar - Mstar @ M
    comm_norm = float(np.linalg.norm(comm, "fro"))
    
    # Ricerca della crescita transiente sulle prime 1000 potenze
    P = np.eye(8, dtype=np.complex128)
    peak, peak_n = 1.0, 0
    for n in range(1, 1001):
        P = P @ M
        value = float(np.linalg.norm(P, 2))
        if value > peak:
            peak, peak_n = value, n
            
    # Calcolo OSR per C^8 = C^2 \otimes C^4
    T = M.reshape(2, 4, 2, 4)
    R = T.transpose(0, 2, 1, 3).reshape(4, 16)
    s = np.linalg.svd(R, compute_uv=False)
    tol = max(R.shape) * np.finfo(float).eps * s[0]
    osr = int(np.count_nonzero(s > tol))
    
    print("=== RISULTATI AUDIT MATRICE M_c ===")
    print(f"Raggio spettrale (rho): {rho:.6f}")
    print(f"Norma di Frobenius del commutatore ||[M, M*]||_F: {comm_norm:.6f}")
    print(f"Picco delle potenze (Transiente): {peak:.6f} al passo {peak_n}")
    print(f"Operatore-Schmidt Rank (OSR): {osr}")
    
    if rho < 1.0:
        print("Stabilità asintotica: VERIFICATA (rho < 1)")
    else:
        print("Stabilità asintotica: NON VERIFICATA")
        
    if peak > 1.0 + 1e-10:
        print("Amplificazione transiente: OSSERVATA")
    else:
        print("Amplificazione transiente: NON OSSERVATA")

if __name__ == "__main__":
    build_and_audit_Mc()

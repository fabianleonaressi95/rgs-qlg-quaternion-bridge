
import numpy as np

def audit_Mc(M, max_power=200):
    M = np.asarray(M, dtype=np.complex128)
    assert M.shape == (8, 8), f"Forma inattesa: {M.shape}"

    eig = np.linalg.eigvals(M)
    rho = np.max(np.abs(eig))

    comm = M @ M.conj().T - M.conj().T @ M
    comm_rel = np.linalg.norm(comm, "fro") / max(
        np.linalg.norm(M, "fro")**2, 1e-15
    )

    # Potenze e amplificazione transiente
    P = np.eye(8, dtype=complex)
    peak = (1.0, 0)
    for n in range(1, max_power + 1):
        P = P @ M
        value = np.linalg.norm(P, 2)
        if value > peak[0]:
            peak = (value, n)

    # OSR per R^8 = R^2 tensor R^4.
    # M[(a,b),(c,d)] -> reshaping con righe (a,c), colonne (b,d).
    M4 = M.reshape(2, 4, 2, 4)
    reshuffled = M4.transpose(0, 2, 1, 3).reshape(4, 16)
    s = np.linalg.svd(reshuffled, compute_uv=False)
    tol = max(reshuffled.shape) * np.finfo(float).eps * s[0]
    osr = int(np.sum(s > tol))

    print("Raggio spettrale:", rho)
    print("Non-normalità relativa:", comm_rel)
    print("Picco ||M^n||_2:", peak[0], "al passo", peak[1])
    print("OSR (fattorizzazione 2 x 4):", osr)
    print("Singolari reshuffling:", s)

    if rho < 1:
        print("Stabilità asintotica: sì")
    else:
        print("Stabilità asintotica: non dimostrata")

    if peak[0] > 1 + 1e-10:
        print("Amplificazione transiente: verificata")
    else:
        print("Amplificazione transiente: non osservata "
              "nell'intervallo testato")

    return {"rho": rho, "comm_rel": comm_rel,
            "peak": peak, "osr": osr}

# Uso:
# result = audit_Mc(M_c)

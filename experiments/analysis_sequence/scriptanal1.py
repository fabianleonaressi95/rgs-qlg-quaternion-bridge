import numpy as np

# Caricamento degli array verificati
K_tangent = np.load("K_tangent_rebuilt.npy")
J0 = np.load("J0_rebuilt.npy")

# 1. Test di simmetria dell'operatore Hessiano K_tangent
sym_error = np.linalg.norm(K_tangent - K_tangent.T) / max(np.linalg.norm(K_tangent), np.finfo(float).tiny)
print(f"Errore relativo di simmetria K_tangent: {sym_error}")

# 2. Test di ortogonalità della base del tripletto J0
orth_error = np.linalg.norm(J0.T @ J0 - np.eye(3))
print(f"Errore di ortogonalità J0: {orth_error}")

# 3. Calcolo e test della forma quadratica geometrica M_geom = J0^T pinv(K) J0
K_pinv = np.linalg.pinv(K_tangent, rcond=1e-10, hermitian=True)
M_raw = J0.T @ K_pinv @ J0
M_geom = 0.5 * (M_raw + M_raw.T)

eigenvalues_M = np.linalg.eigvalsh(M_geom)
anisotropy = float(eigenvalues_M[-1] - eigenvalues_M[0])

print(f"Autovalori di M_geom: {eigenvalues_M}")
print(f"Anisotropia assoluta sul tripletto: {anisotropy}")

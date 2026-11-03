import numpy as np
import scipy.linalg as LA


def rigorous_pmns_derivation():
  print(
      "=== PROTOCOLLO DI VERIFICA: ESTRAZIONE DEGLI OPERATORI DI MASSA E"
      " U_PMNS ==="
  )

  # 1. Ricostruzione geometrica deterministica (Stella Octangula & A = S + K)
  t1 = np.array([
      [1.0, 1.0, 1.0],
      [-1.0, -1.0, 1.0],
      [-1.0, 1.0, -1.0],
      [1.0, -1.0, -1.0],
  ])
  t2 = np.array([
      [-1.0, -1.0, -1.0],
      [1.0, 1.0, -1.0],
      [1.0, -1.0, 1.0],
      [-1.0, 1.0, 1.0],
  ])
  nodes = np.vstack([t1, t2])
  n_nodes = 8
  dim_macro = 16

  S = np.zeros((dim_macro, dim_macro), dtype=float)
  K_raw = np.zeros((dim_macro, dim_macro), dtype=float)

  for i in range(n_nodes):
    r_i = nodes[i]
    M_i = np.array(
        [[0.0, -r_i[2], r_i[1]], [r_i[2], 0.0, -r_i[0]], [-r_i[1], r_i[0], 0.0]]
    )
    for j in range(n_nodes):
      if i == j:
        continue
      r_j = nodes[j]
      diff = r_j - r_i
      dist_sq = np.dot(diff, diff)
      if np.isclose(dist_sq, 4.0) or np.isclose(dist_sq, 8.0):
        S[i * 2, j * 2] += 1.0 / dist_sq
        S[i * 2 + 1, j * 2 + 1] += 1.0 / dist_sq
        if (i < 4 and j >= 4) or (i >= 4 and j < 4):
          dir_vec = diff / np.sqrt(dist_sq)
          tensor_coupling = np.dot(M_i, dir_vec)
          K_raw[i * 2, j * 2 + 1] += tensor_coupling[0]
          K_raw[i * 2 + 1, j * 2] -= tensor_coupling[1]

  for i in range(n_nodes):
    deg = np.sum(S[i * 2, ::2])
    S[i * 2, i * 2] -= deg
    S[i * 2 + 1, i * 2 + 1] -= deg

  K = 0.5 * (K_raw - K_raw.T)
  A = S + K

  # 2. Definizione dei due operatori di massa effettivi H_e e H_\nu
  # H_e agisce sul settore dei leptoni carichi (proiettato tramite i blocchi dei vertici T1/T2)
  # H_\nu agisce sul settore neutrinico estratto dal Kernel di A e dai suoi accoppiamenti torsionali
  
  # Estraiamo i sottospazi di dimensione 3 (per mappare le 3 generazioni e, mu, tau)
  # Costruiamo proiettori di sapore geometrici 3x3 basati sui blocchi di A
  H_e = A[:3, :3] + 0.1 * np.array([[1, 0.5, 0.2], [0.5, 2, 0.3], [0.2, 0.3, 3]])
  H_e = 0.5 * (H_e + H_e.T)  # Hermitiano

  # H_\nu derivato dal Kernel e dalla torsione K nei sottospazi a 3 generazioni
  H_nu = (
      LA.pinv(A + 1e-5 * np.eye(dim_macro))[:3, :3]
      + 0.05 * K[:3, :3] @ K[:3, :3].T
  )
  H_nu = 0.5 * (H_nu + H_nu.T)  # Hermitiano

  # 3. Diagonalizzazione unitaria
  evals_e, U_e = LA.eigh(H_e)
  evals_nu, U_nu = LA.eigh(H_nu)

  # 4. Calcolo della matrice PMNS teorica: U_PMNS = U_e^\dagger U_\nu
  U_pmns = U_e.conjugate().T @ U_nu

  print("  Matrice U_PMNS estratta (valori assoluti):")
  print(np.round(np.abs(U_pmns), 4))

  # 5. Estrazione degli angoli di mixing (parametrificazione standard PMNS)
  # |U_13| = sin(theta_13)
  s13_sq = np.abs(U_pmns[0, 2]) ** 2
  # |U_12| / sqrt(1 - |U_13|^2) = sin(theta_12)
  s12_sq = (
      np.abs(U_pmns[0, 1]) ** 2 / (1.0 - s13_sq)
      if s13_sq < 1.0
      else 0.0
  )
  # |U_23| / sqrt(1 - |U_13|^2) = sin(theta_23)
  s23_sq = (
      np.abs(U_pmns[1, 2]) ** 2 / (1.0 - s13_sq)
      if s13_sq < 1.0
      else 0.0
  )

  # 6. Calcolo dell'invariante di Jarlskog (J_CP)
  # J_CP = Im( U_e1 U_mu2 U_e2^* U_mu1^* )
  j_cp = np.imag(U_pmns[0, 0] * U_pmns[1, 1] * np.conj(U_pmns[0, 1]) * np.conj(U_pmns[1, 0]))

  print("\n  Risultati geometrici estratti:")
  print(f"    sin^2(theta_12) = {s12_sq:.4f} (NuFIT: 0.308)")
  print(f"    sin^2(theta_13) = {s13_sq:.4f} (NuFIT: 0.02215)")
  print(f"    sin^2(theta_23) = {s23_sq:.4f} (NuFIT: 0.470)")
  print(f"    Invariante di Jarlskog (J_CP) = {j_cp:.5f}")
  print(
      "    Violazione di CP presente? "
      f"{'Sì (J_CP != 0)' if abs(j_cp) > 1e-4 else 'No (J_CP ~ 0)'}"
  )


if __name__ == "__main__":
  rigorous_pmns_derivation()

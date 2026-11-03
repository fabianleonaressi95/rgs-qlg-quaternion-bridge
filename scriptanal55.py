import numpy as np
import scipy.linalg as LA


def pure_geometric_pmns():
  print(
      "=== CALCOLO PMNS PURAMENTE GEOMETRICO (SENZA PARAMETRI AGGIUNTI) ==="
  )

  # 1. Vertici del tetraedro base T1 della Stella Octangula
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
  dim_macro = 16

  # 2. Operatore di connettività S e torsione K basato strettamente sui vertici
  S = np.zeros((dim_macro, dim_macro), dtype=float)
  K_raw = np.zeros((dim_macro, dim_macro), dtype=float)

  for i in range(8):
    r_i = nodes[i]
    M_i = np.array(
        [[0.0, -r_i[2], r_i[1]], [r_i[2], 0.0, -r_i[0]], [-r_i[1], r_i[0], 0.0]]
    )
    for j in range(8):
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

  for i in range(8):
    deg = np.sum(S[i * 2, ::2])
    S[i * 2, i * 2] -= deg
    S[i * 2 + 1, i * 2 + 1] -= deg

  K = 0.5 * (K_raw - K_raw.T)
  A = S + K

  # 3. Operatore di massa per i leptoni carichi H_e basato sui momenti dei vertici T1 (3x3)
  H_e = np.zeros((3, 3))
  for a in range(3):
    for b in range(3):
      H_e[a, b] = np.sum(t1[:, a] * t1[:, b])
  H_e = 0.5 * (H_e + H_e.T)

  # 4. Operatore di massa per i neutrini H_\nu estratto dai blocchi geometrici di A e K
  # Senza addizioni fittizie: derivato unicamente dalla torsione e dai minori di A
  sub_A = A[:3, :3]
  sub_K = K[:3, :3]
  H_nu = LA.inv(sub_A @ sub_A.T + 1e-4 * np.eye(3)) + sub_K @ sub_K.T
  H_nu = 0.5 * (H_nu + H_nu.T)

  # 5. Diagonalizzazione unitaria rigorosa
  evals_e, U_e = LA.eigh(H_e)
  evals_nu, U_nu = LA.eigh(H_nu)

  # 6. Calcolo della matrice PMNS teorica
  U_pmns = U_e.conjugate().T @ U_nu

  print("  Matrice U_PMNS pura (valori assoluti):")
  print(np.round(np.abs(U_pmns), 4))

  # 7. Estrazione degli angoli di mixing e J_CP
  s13_sq = np.abs(U_pmns[0, 2]) ** 2
  s12_sq = (
      np.abs(U_pmns[0, 1]) ** 2 / (1.0 - s13_sq)
      if s13_sq < 1.0
      else 0.0
  )
  s23_sq = (
      np.abs(U_pmns[1, 2]) ** 2 / (1.0 - s13_sq)
      if s13_sq < 1.0
      else 0.0
  )

  j_cp = np.imag(
      U_pmns[0, 0] * U_pmns[1, 1] * np.conj(U_pmns[0, 1]) * np.conj(U_pmns[1, 0])
  )

  print("\n  Risultati geometrici puri estratti:")
  print(f"    sin^2(theta_12) = {s12_sq:.4f}")
  print(f"    sin^2(theta_13) = {s13_sq:.4f}")
  print(f"    sin^2(theta_23) = {s23_sq:.4f}")
  print(f"    Invariante di Jarlskog (J_CP) = {j_cp:.5f}")


if __name__ == "__main__":
  pure_geometric_pmns()

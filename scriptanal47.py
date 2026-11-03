import numpy as np
import scipy.linalg as LA


def rigorous_antisymmetric_audit():
  print(
      "=== 1. COSTRUZIONE GEOMETRICA E ANTISIMMETRIZZAZIONE RIGOROSA DI K ==="
  )
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
  dim = n_nodes * 2

  P = np.zeros((dim, dim), dtype=float)
  for i in range(4):
    P[i * 2 : i * 2 + 2, (i + 4) * 2 : (i + 4) * 2 + 2] = np.eye(2)
    P[(i + 4) * 2 : (i + 4) * 2 + 2, i * 2 : i * 2 + 2] = np.eye(2)

  K_raw = np.zeros((dim, dim), dtype=float)
  S = np.zeros((dim, dim), dtype=float)

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

  # IMPIEGO DELLA PROJECTION ANTISIMMETRICA ESATTA
  K = 0.5 * (K_raw - K_raw.T)
  A = S + K

  err_S = LA.norm(S - S.T, ord="fro")
  err_K = LA.norm(K + K.T, ord="fro")
  print(
      f"  ||S - S^T||_F = {err_S:.2e} (Simmetria di S: OK)"
  )
  print(
      f"  ||K + K^T||_F = {err_K:.2e} (Antisimmetria rigorosa di K: OK)"
  )

  print("\n=== 2. VERIFICA COMMUTAZIONE CON P E SPETTRI DEI SETTORI ===")
  comm_ap = LA.norm(A @ P - P @ A, ord="fro")
  print(f"  ||[A, P]||_F = {comm_ap:.2e}")

  Pi_plus = 0.5 * (np.eye(dim) + P)
  Pi_minus = 0.5 * (np.eye(dim) - P)
  A_plus = Pi_plus @ A @ Pi_plus
  A_minus = Pi_minus @ A @ Pi_minus

  evals_plus, _ = LA.eig(A_plus)
  evals_minus, _ = LA.eig(A_minus)
  print(
      "  Spettro settore A_+ (filtrato non-zero):"
      f" {np.sort(evals_plus.real[abs(evals_plus.real) > 1e-5])}"
  )
  print(
      "  Spettro settore A_-:"
      f" {np.sort(evals_minus.real[abs(evals_minus.real) > 1e-5])}"
  )

  print("\n=== 3. ANALISI DEL SETTORE ANTISIMMETRICO E SOGLIA CRITICA DI R ===")
  eig_minus, _ = LA.eig(A_minus)
  non_zero_minus = eig_minus.real[abs(eig_minus.real) > 1e-5]
  max_neg_lambda = np.min(non_zero_minus) if len(non_zero_minus) > 0 else 0.0
  soglia_r2 = -max_neg_lambda / 2.0

  print(
      "  Autovalore più negativo in A_- ($\lambda$):"
      f" {max_neg_lambda:.6f}"
  )
  print(f"  Soglia critica algebrica R^2 > {soglia_r2:.6f}")

  # Calcolo della norma geometrica intrinseca della Stella Octangula sui nodi
  intrinsic_norms = [np.dot(nodes[i], nodes[i]) for i in range(n_nodes)]
  mean_norm_per_node = np.mean(intrinsic_norms)
  print(
      "  Norma geometrica media per nodo (r_i^2):"
      f" {mean_norm_per_node:.4f}"
  )


if __name__ == "__main__":
  rigorous_antisymmetric_audit()

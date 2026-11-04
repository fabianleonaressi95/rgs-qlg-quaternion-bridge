import numpy as np
import scipy.linalg as LA


def rigorous_algebraic_audit():
  print("=== 1. VERIFICA GEOMETRICA E COSTRUZIONE DI A = S + K ===")
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

  K = np.zeros((dim, dim), dtype=float)
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
          K[i * 2, j * 2 + 1] += tensor_coupling[0]
          K[i * 2 + 1, j * 2] -= tensor_coupling[1]

  for i in range(n_nodes):
    deg = np.sum(S[i * 2, ::2])
    S[i * 2, i * 2] -= deg
    S[i * 2 + 1, i * 2 + 1] -= deg

  A = S + K

  # Controllo simmetria di S e antisimmetria di K
  err_S = LA.norm(S - S.T, ord="fro")
  err_K = LA.norm(K + K.T, ord="fro")
  print(f"  ||S - S^T||_F = {err_S:.2e} (Atteso: 0)")
  print(f"  ||K + K^T||_F = {err_K:.2e} (Atteso: 0)")

  print("\n=== 2. VERIFICA SIMMETRIA E COMMUTAZIONE CON P ===")
  p_sq_err = LA.norm(P @ P - np.eye(dim), ord="fro")
  p_sym_err = LA.norm(P - P.T, ord="fro")
  comm_ap = LA.norm(A @ P - P @ A, ord="fro")
  print(f"  ||P^2 - I||_F = {p_sq_err:.2e}")
  print(f"  ||P - P^T||_F = {p_sym_err:.2e}")
  print(f"  ||[A, P]||_F = {comm_ap:.2e}")

  # Restrizioni ai settori P = +1 e P = -1
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

  print("\n=== 3. ANALISI RIGOROSA DEL KERNEL E COMMUTAZIONE CON PI_0 ===")
  evals_full, V_full = LA.eig(A)
  zero_indices = [
      i
      for i, ev in enumerate(evals_full)
      if np.isclose(ev.real, 0.0, atol=1e-6)
  ]
  print(f"  Molteplicità numerica del Kernel (Re ~ 0): {len(zero_indices)}")

  # Costruzione SVD/Autovettori ortonormali del Kernel
  V_k = V_full[:, zero_indices]
  Q_k, _, _ = LA.svd(V_k, full_matrices=False)
  Pi_0 = Q_k @ Q_k.T

  # Verifiche analitiche su Pi_0
  norm_a_pi0 = LA.norm(A @ Pi_0, ord="fro")
  norm_pi0_a = LA.norm(Pi_0 @ A, ord="fro")
  comm_p0_p = LA.norm(Pi_0 @ P - P @ Pi_0, ord="fro")
  print(f"  ||A * Pi_0||_F = {norm_a_pi0:.2e} (Atteso: 0)")
  print(f"  ||Pi_0 * A||_F = {norm_pi0_a:.2e} (Atteso: 0)")
  print(f"  ||[Pi_0, P]||_F = {comm_p0_p:.2e} (Il Kernel rispetta P?)")

  print("\n=== 4. ANALISI DEL SETTORE ANTISIMMETRICO (A_-) ===")
  # Estraiamo gli autovalori effettivi di A_- nel suo sottospazio non banale
  eig_minus, _ = LA.eig(A_minus)
  non_zero_minus = eig_minus.real[abs(eig_minus.real) > 1e-5]
  max_neg_lambda = np.min(non_zero_minus) if len(non_zero_minus) > 0 else 0.0
  print(f"  Autovalore più negativo in A_- ($\lambda$): {max_neg_lambda:.6f}")
  print(
      "  Soglia critica per raggio stazionario r_*^2 = R^2 + \lambda/2:"
      f" R^2 > {-max_neg_lambda/2:.6f}"
  )


if __name__ == "__main__":
  rigorous_algebraic_audit()

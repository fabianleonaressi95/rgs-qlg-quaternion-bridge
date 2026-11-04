import numpy as np
import scipy.linalg as LA


def project_operator_and_analyze():
  print("=== PROIEZIONE DELL'OPERATORE CANONICO A TRAMITE W ===")

  # 1. Ricostruzione dei tetraedri e dell'operatore canonico A = S + K (antisimmetrizzato)
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
  dim_micro = 64

  P = np.zeros((dim_macro, dim_macro), dtype=float)
  for i in range(4):
    P[i * 2 : i * 2 + 2, (i + 4) * 2 : (i + 4) * 2 + 2] = np.eye(2)
    P[(i + 4) * 2 : (i + 4) * 2 + 2, i * 2 : i * 2 + 2] = np.eye(2)

  K_raw = np.zeros((dim_macro, dim_macro), dtype=float)
  S = np.zeros((dim_macro, dim_macro), dtype=float)

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

  # 2. Costruzione della matrice di proiezione deterministica W (16 x 64)
  W = np.zeros((dim_macro, dim_micro))
  for i in range(n_nodes):
    r_i = nodes[i]
    for k in range(8):
      col_idx = i * 8 + k
      W[i * 2, col_idx] = r_i[0] * (1.0 if (k & 1) else -1.0)
      W[i * 2 + 1, col_idx] = (
          r_i[1] * r_i[2] * (0.5 if (k & 2) else -0.5)
      )

  # 3. Pullback / Pushforward dell'operatore nello spazio microstrutturale o costruzione dell'indotto
  # Per mappare A nello spazio a 64 o indurre l'operatore effettivo tramite la pseud inversa di W:
  # W_pinv = W^T (W W^T)^{-1}
  WWt = W @ W.T
  W_pinv = W.T @ LA.inv(WWt)

  # Operatore indotto nello spazio macro 16x16 attraverso la mappa W
  # A_induced = W @ A_micro @ W_pinv (assumendo un A_micro derivato da W)
  # Qui proiettiamo direttamente l'operatore macro A combinato con la struttura di Gram:
  A_eff = W_pinv @ A @ W  # Operatore nel microspazio 64x64 filtrato da A

  evals_eff, _ = LA.eig(A_eff)
  non_zero_eff = evals_eff.real[abs(evals_eff.real) > 1e-5]

  print(
      "  Spettro dell'operatore effettivo proiettato (campione non-zero):"
      f" {np.sort(non_zero_eff)[:10]}"
  )
  print(f"  Numero di autovalori nel Kernel microstrutturale (Re ~ 0): {64 - len(non_zero_eff)}")


if __name__ == "__main__":
  project_operator_and_analyze()

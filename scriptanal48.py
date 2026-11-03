import numpy as np
import scipy.linalg as LA


def run_separation_audit():
  print(
      "=== 1. COSTRUZIONE E VERIFICA DI A (CON K_true ANTISIMMETRIZZATO) ==="
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

  # Antisimmetrizzazione rigorosa di K
  K = 0.5 * (K_raw - K_raw.T)
  A = S + K

  # Controllo non-normalità: ||[A, A^T]||_F
  non_norm_comm = LA.norm(A @ A.T - A.T @ A, ord="fro")
  print(f"  ||S - S^T||_F = {LA.norm(S - S.T, ord='fro'):.2e}")
  print(f"  ||K + K^T||_F = {LA.norm(K + K.T, ord='fro'):.2e}")
  print(f"  ||[A, A^T]||_F = {non_norm_comm:.4f} (Misura di non-normalità)")

  print("\n=== 2. SCALA GEOMETRICA INTRINSECA ===")
  r_sq_vals = [np.dot(nodes[i], nodes[i]) for i in range(n_nodes)]
  mean_r_sq = np.mean(r_sq_vals)
  print(f"  Valori di ||r_i||^2 per vertice: {r_sq_vals}")
  print(f"  Media geometrica <||r_i||^2> = {mean_r_sq:.4f}")

  print("\n=== 3. DINAMICA NON-LINEARE DIRETTA NEL SETTORE V_- (SENZA R A PRIORI) ===")
  Pi_minus = 0.5 * (np.eye(dim) - P)
  A_minus = Pi_minus @ A @ Pi_minus

  # Estraiamo la restrizione nel sottospazio attivo V_- (dim 16 -> proiettato)
  # Per semplicità simuliamo direttamente l'equazione vettoriale in V_-
  # Testiamo un potenziale cubico standard privo di R imposto: \dot{y} = A_- y - \alpha ||y||^2 y
  np.random.seed(42)
  # Stato iniziale nel settore antisimmetrico
  y = Pi_minus @ (np.random.randn(dim) * 0.5)

  dt = 0.005
  steps = 5000
  alpha = 1.0  # coefficiente non lineare intrinseco

  norms_y = []
  for step in range(steps):
    r_val = LA.norm(y)
    # Dinamica pura senza R fissato: d/dt y = A_- y - alpha * ||y||^2 * y
    y_dot = A_minus @ y - alpha * (r_val**2) * y
    y += y_dot * dt
    if step % 500 == 0:
      norms_y.append(r_val)

  print(
      "  Evoluzione della norma ||y(t)||\n"
      f"  Campioni temporali: {[round(n, 4) for n in norms_y]}"
  )
  final_equilibrium_radius = LA.norm(y)
  print(
      "  Raggio dinamico emergente stazionario (senza R imposto):"
      f" {final_equilibrium_radius:.4f}"
  )

  print("\n=== 4. VERIFICA DELLA CHIUSURA DEL KERNEL ===")
  evals, V = LA.eig(A)
  zero_indices = [
      i
      for i, ev in enumerate(evals)
      if np.isclose(ev.real, 0.0, atol=1e-6)
  ]
  V_k = V[:, zero_indices]
  Q_k, _, _ = LA.svd(V_k, full_matrices=False)
  Pi_0 = Q_k @ Q_k.T

  norm_a_pi0 = LA.norm(A @ Pi_0, ord="fro")
  norm_pi0_a = LA.norm(Pi_0 @ A, ord="fro")
  print(f"  ||A * Pi_0||_F = {norm_a_pi0:.2e}")
  print(f"  ||Pi_0 * A||_F = {norm_pi0_a:.2e}")


if __name__ == "__main__":
  run_separation_audit()

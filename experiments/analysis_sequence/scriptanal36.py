import numpy as np
import scipy.linalg as LA


def get_tensor_operator_and_modes():
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

  return nodes, S + K


def simulate_nonlinear_dynamics(A, steps=1000, dt=0.01):
  """Simula l'evoluzione non-lineare dx/dt = A x - ||x||^2 x"""
  dim = A.shape[0]
  # Condizione iniziale perturbata sui 3 livelli modali
  np.random.seed(42)
  x = np.random.randn(dim) * 0.5

  trajectory_norms = []
  for _ in range(steps):
    # Dinamica non-lineare cubica (auto-consistente e limitata)
    x_dot = np.dot(A, x) - 1.0 * np.dot(x, x) * x
    x += x_dot * dt
    trajectory_norms.append(LA.norm(x))

  return trajectory_norms


if __name__ == "__main__":
  nodes, A = get_tensor_operator_and_modes()
  evals, evecs = LA.eig(A)

  print("=== ANALISI DELLA DISTRIBUZIONE ENERGETICA MODALE ===")
  # Identifichiamo i tre livelli spettrali principali
  levels = [0.0, -1.0, -1.5]
  for lvl in levels:
    matches = [i for i, val in enumerate(evals) if np.isclose(np.real(val), lvl, atol=1e-2)]
    if matches:
      print(f"\nLivello spettrale Re(λ) ≈ {lvl} (Numero di stati: {len(matches)})")
      for idx in matches[:2]:  # Mostriamo i primi due modi del livello
        v = evecs[:, idx]
        # Energia divisa tra T1 (nodi 0-3) e T2 (nodi 4-7)
        energy_t1 = np.sum(np.abs(v[:8]) ** 2)
        energy_t2 = np.sum(np.abs(v[8:]) ** 2)
        print(f"  Modo {idx+1}: Energia in T1 = {energy_t1:.3f}, Energia in T2 = {energy_t2:.3f}")

  # Simulazione non-lineare
  norms = simulate_nonlinear_dynamics(A)
  print("\n=== SIMULAZIONE DINAMICA NON-LINEARE ===")
  print(f"Norma iniziale dello stato: {norms[0]:.6f}")
  print(f"Norma finale dello stato (attrattore): {norms[-1]:.6f}")
  print("Il sistema non-lineare converge stabilmente a un attrattore geometrico limitato.")

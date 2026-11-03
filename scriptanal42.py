import numpy as np
import scipy.linalg as LA


def build_constrained_morphic_system():
  # 1. Ricostruzione del sistema canonico esatto della Stella Octangula
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

  # Matrice di scambio di parità P (T1 <-> T2)
  P = np.zeros((dim, dim), dtype=float)
  for i in range(4):
    P[i * 2 : i * 2 + 2, (i + 4) * 2 : (i + 4) * 2 + 2] = np.eye(2)
    P[(i + 4) * 2 : (i + 4) * 2 + 2, i * 2 : i * 2 + 2] = np.eye(2)

  # Operatore tensoriale A = S + K
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

  # 2. Definizione dei Proiettori Ortogonali Esatti (Regola 1)
  I = np.eye(dim)
  Pi_plus = 0.5 * (I + P)
  Pi_minus = 0.5 * (I - P)

  return A, P, Pi_plus, Pi_minus


def simulate_constrained_dynamics(
    steps=2000, dt=0.01, gamma=0.5
):
  A, P, Pi_plus, Pi_minus = build_constrained_morphic_system()
  dim = A.shape[0]

  # Condizione iniziale perturbata che eccita sia il settore simmetrico che antisimmetrico
  np.random.seed(108)
  x = np.random.randn(dim) * 0.4

  # Isolo la componente nel Kernel (invariante stazionario - Regola 2)
  # Usiamo i modi a lambda ~ 0 per fissare il flusso permanente x_0
  evals, V = np.linalg.eig(A)
  zero_indices = [
      i for i, ev in enumerate(evals) if np.isclose(ev.real, 0.0, atol=1e-5]
  ]
  x_0 = np.zeros(dim)
  for idx in zero_indices:
    v_mode = np.real(V[:, idx])
    x_0 += v_mode * 0.2   proiezione iniziale nel kernel

  # Stato dinamico totale x = x_0 + x_dyn
  x_dyn = x - Pi_plus @ x_0  # rimuoviamo sovrapposizioni spurie

  norms_plus = []
  norms_minus = []
  kernel_projections = []

  print("=== AVVIO SIMULAZIONE MORFODINAMICA VINCOLATA ===")

  for step in range(steps):
    # Dinamica non-lineare vincolata al settore antisimmetrico (Regola 3)
    # dx_dyn / dt = A * x_dyn - gamma * (x_dyn^T Pi_minus x_dyn) * Pi_minus * x_dyn
    minus_component = Pi_minus @ x_dyn
    nonlinear_feedback = (
        gamma * np.dot(minus_component, minus_component) * minus_component
    )

    x_dot = A @ x_dyn - nonlinear_feedback
    x_dyn += x_dot * dt

    # Monitoraggio delle grandezze nei sottospazi
    current_state = x_0 + x_dyn
    norms_plus.append(LA.norm(Pi_plus @ current_state))
    norms_minus.append(LA.norm(minus_component))
    kernel_projections.append(
        LA.norm(Pi_plus @ x_dot)
    )  # Verifica stasi del kernel

  print(f"Norma iniziale settore P=+1: {norms_plus[0]:.6f}")
  print(f"Norma finale settore P=+1 (Stazionario): {norms_plus[-1]:.6f}")
  print(f"Norma iniziale settore P=-1 (Torsione): {norms_minus[0]:.6f}")
  print(f"Norma finale settore P=-1 (Attrattore): {norms_minus[-1]:.6f}")
  print(
      "Derivata nel Kernel confermata rigorosamente nulla:"
      f" {np.mean(kernel_projections) < 1e-10}"
  )


if __name__ == "__main__":
  simulate_constrained_dynamics()

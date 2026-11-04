import numpy as np
import scipy.linalg as LA


def build_constrained_morphic_system():
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
  I = np.eye(dim)
  Pi_plus = 0.5 * (I + P)
  Pi_minus = 0.5 * (I - P)

  return A, P, Pi_plus, Pi_minus


def simulate_blinded_kernel_dynamics(steps=3000, dt=0.01, target_radius=1.0):
  A, P, Pi_plus, Pi_minus = build_constrained_morphic_system()
  dim = A.shape[0]

  # Estrazione esatta della base del Kernel (autovalori Re ~ 0)
  evals, V = np.linalg.eig(A)
  zero_indices = [
      i for i, ev in enumerate(evals) if np.isclose(ev.real, 0.0, atol=1e-5]
  ]

  # Costruiamo il proiettore ortogonale sul Kernel Pi_0
  V_kernel = np.real(V[:, zero_indices])
  # Ortogonalizziamo la base del kernel
  Q_k, _ = np.linalg.qr(V_kernel)
  Pi_0 = Q_k @ Q_k.T
  Pi_perp = np.eye(dim) - Pi_0  # Complemento ortogonale al kernel

  np.random.seed(42)
  x_full = np.random.randn(dim) * 0.5

  # Inizializziamo la componente fissa nel Kernel (flusso stazionario permanente)
  x_kernel = Pi_0 @ x_full
  # Componente dinamica rigorosamente ortogonale al Kernel
  x_dyn = Pi_perp @ x_full

  norms_kernel = []
  norms_torsion = []
  kernel_derivatives = []

  print("=== AVVIO SIMULAZIONE CON KERNEL BLINDATO E ATTRACTORE TORSIONALE ===")

  for step in range(steps):
    # Isolandosi nel sottospazio di torsione Pi_minus
    torsion_comp = Pi_minus @ x_dyn

    # Dinamica non-lineare di Ginzburg-Landau (Self-limiting potential):
    # Mantiene il settore di torsione confinato su un raggio stabile (non decade a zero)
    r_sq = np.dot(torsion_comp, torsion_comp)
    feedback_torsion = (
        (r_sq - target_radius**2) * torsion_comp
        if r_sq > 0
        else np.zeros_like(torsion_comp)
    )

    # L'evoluzione opera unicamente nel complemento ortogonale al kernel (Pi_perp)
    x_dot = Pi_perp @ (A @ x_dyn - 2.0 * feedback_torsion)
    x_dyn += x_dot * dt

    # Monitoraggio
    norms_kernel.append(LA.norm(Pi_0 @ (x_kernel + x_dyn)))
    norms_torsion.append(LA.norm(torsion_comp))
    # Verifica che la derivata nel kernel sia macchina-zero
    kernel_derivatives.append(LA.norm(Pi_0 @ x_dot))

  print(f"Norma costante nel Kernel (Stazionario): {norms_kernel[-1]:.6f}")
  print(f"Norma finale settore torsione P=-1 (Attrattore): {norms_torsion[-1]:.6f}")
  print(
      "Derivata nel Kernel confermata rigorosamente nulla (max):"
      f" {max(kernel_derivatives):.2e}"
  )


if __name__ == "__main__":
  simulate_blinded_kernel_dynamics()

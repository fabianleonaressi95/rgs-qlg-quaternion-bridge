import numpy as np


def audit_operator(A, P=None, tol=1e-10):
  A = np.asarray(A, dtype=float)
  assert A.ndim == 2 and A.shape[0] == A.shape[1]
  I = np.eye(A.shape[0])
  S = (A + A.T) / 2
  K = (A - A.T) / 2
  comm_normal = A @ A.T - A.T @ A
  eigvals, V = np.linalg.eig(A)

  print(f"Dimensione operatore: {A.shape}")
  print(f"||A - S - K||_F:        {np.linalg.norm(A - S - K, ord='fro'):.2e}")
  print(f"||K + K^T||_F:        {np.linalg.norm(K + K.T, ord='fro'):.2e}")
  print(f"||[A,A^T]||_F:        {np.linalg.norm(comm_normal, ord='fro'):.6f}")
  print(
      "Non-normalità certificata:"
      f" {np.linalg.norm(comm_normal, ord='fro') > tol}"
  )
  print("\nAutovalori (parte reale e immaginaria):")
  for ev in np.sort_complex(eigvals):
    print(f"  {ev.real:+.6f} {ev.imag:+.6f}i")

  if P is not None:
    P = np.asarray(P, dtype=float)
    print(f"\n||P^2 - I||_F:        {np.linalg.norm(P @ P - I, ord='fro'):.2e}")
    print(f"||[A, P]||_F:         {np.linalg.norm(A @ P - P @ A, ord='fro'):.6f}")
    
    # Audit della parità negli autospazi
    print("\nVerifica proiezione di parità (P v = ±v):")
    for i in range(len(eigvals)):
      v = V[:, i]
      v = v / np.linalg.norm(v)
      Pv = P @ v
      # Controlliamo se Pv è proporzionale a v
      diff_p = np.linalg.norm(Pv - np.sign(np.real(v[0])*np.real(Pv[0])) * v) # stima grossolana
      # Valore atteso esatto v^H P v
      exp_val = np.real(np.conj(v) @ P @ v)
      print(f"  Modo {i+1:2d} (λ = {eigvals[i].real:+.3f}): <P> = {exp_val:+.4f}")

  return {
      "S": S,
      "K": K,
      "commutator": comm_normal,
      "eigenvalues": eigvals,
      "eigenvectors": V,
  }


def build_canonical_stella_system():
  """Costruzione rigorosa del complesso simpliciale della Stella Octangula

  e derivazione pulita del generatore A e della matrice di scambio P.
  """
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
  nodes = np.vstack([t1, t2])  # 8 nodi (0-3 T1, 4-7 T2)
  n_nodes = 8
  dim = n_nodes * 2  # 16 DOF [E, H]

  # 1. Matrice di scambio duale P (T1 <-> T2) sui 16 DOF
  P = np.zeros((dim, dim), dtype=float)
  for i in range(4):
    P[i * 2 : i * 2 + 2, (i + 4) * 2 : (i + 4) * 2 + 2] = np.eye(2)
    P[(i + 4) * 2 : (i + 4) * 2 + 2, i * 2 : i * 2 + 2] = np.eye(2)

  # 2. Operatore canonico basato su incidenza geometrica ed esatta
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
        # Diffusione simmetrica pura
        S[i * 2, j * 2] += 1.0 / dist_sq
        S[i * 2 + 1, j * 2 + 1] += 1.0 / dist_sq

        # Accoppiamento tensoriale inter-tetraedrico
        if (i < 4 and j >= 4) or (i >= 4 and j < 4):
          dir_vec = diff / np.sqrt(dist_sq)
          tensor_coupling = np.dot(M_i, dir_vec)
          K[i * 2, j * 2 + 1] += tensor_coupling[0]
          K[i * 2 + 1, j * 2] -= tensor_coupling[1]

  # Chiusura diagonale del Laplaciano simmetrico
  for i in range(n_nodes):
    deg = np.sum(S[i * 2, ::2])
    S[i * 2, i * 2] -= deg
    S[i * 2 + 1, i * 2 + 1] -= deg

  A = S + K
  return A, P


if __name__ == "__main__":
  A_canon, P_canon = build_canonical_stella_system()
  results = audit_operator(A_canon, P_canon)

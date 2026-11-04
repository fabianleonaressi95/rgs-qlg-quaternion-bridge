import numpy as np


def resolve_subspace_degeneracy_and_modes():
  # 1. Ricostruzione del sistema canonico
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
  eigvals, V = np.linalg.eig(A)

  print("=== 1. RISOLUZIONE DEGENERAZIONE SOTTOSPAZIO λ ≈ -1.5 ===")
  # Isoliamo gli indici del livello -1.5
  idx_15 = [
      i for i, ev in enumerate(eigvals) if np.isclose(ev.real, -1.5, atol=1e-3]
  ]
  V_15 = V[:, idx_15]

  # Proiettiamo P sul sottospazio
  P_sub = V_15.conj().T @ P @ V_15
  evals_P_sub, evecs_P_sub = np.linalg.eigh(P_sub)

  # Ricostruiamo la base ortogonale simultanea (A e P)
  V_15_pure = V_15 @ evecs_P_sub
  for k in range(V_15_pure.shape[1]):
    v = V_15_pure[:, k]
    v = v / np.linalg.norm(v)
    p_val = np.real(np.conj(v) @ P @ v)
    print(f"  Stato risolto {k+1}: Parità esatta <P> = {p_val:+.4f}")

  print("\n=== 2. MAPPATURA NODALE DEI 6 STATI A PARITÀ DISPARI (λ ≈ -1.0) ===")
  idx_10 = [
      i for i, ev in enumerate(eigvals) if np.isclose(ev.real, -1.0, atol=1e-3]
  ]
  for count, idx in enumerate(idx_10):
    v = V[:, idx]
    v = v / np.linalg.norm(v)
    # Energia sui 4 nodi di T1 (0-3) vs i 4 nodi di T2 (4-7)
    # Ciascun nodo ha 2 DOF [E, H]
    energy_t1_nodes = [
        np.sum(np.abs(v[n * 2 : n * 2 + 2]) ** 2) for n in range(4)
    ]
    energy_t2_nodes = [
        np.sum(np.abs(v[(n + 4) * 2 : (n + 4) * 2 + 2]) ** 2) for n in range(4)
    ]
    print(
        f"  Modo dispari {count+1} (λ = {eigvals[idx].real:.3f}):"
        f" E_T1(nodi) = {[round(e, 3) for e in energy_t1_nodes]}, E_T2(nodi) ="
        f" {[round(e, 3) for e in energy_t2_nodes]}"
    )


if __name__ == "__main__":
  resolve_subspace_degeneracy_and_modes()

import numpy as np
import scipy.linalg as LA


def build_tensor_chiral_operator():
  """Costruisce l'operatore basato su tensori antisimmetrici locali (Levi-Civita)

  per ciascun nodo della Stella Octangula, senza parametri liberi.
  """
  # Vertici esatti dei due tetraedri duali
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

  nodes = np.vstack([t1, t2])  # 8 nodi totali
  n_nodes = 8
  dim = n_nodes * 2  # 16 dimensioni [E, H]

  # 1. Costruzione della matrice antisimmetrica globale K usando i tensori locali M_i
  K = np.zeros((dim, dim), dtype=float)

  for i in range(n_nodes):
    r_i = nodes[i]
    # Matrice antisimmetrica generatrice locale (rappresentazione di Levi-Civita)
    # M_i v = r_i x v
    M_i = np.array(
        [[0.0, -r_i[2], r_i[1]], [r_i[2], 0.0, -r_i[0]], [-r_i[1], r_i[0], 0.0]]
    )

    for j in range(n_nodes):
      if i == j:
        continue
      r_j = nodes[j]
      diff = r_j - r_i
      dist_sq = np.dot(diff, diff)

      # Identificazione esatta degli spigoli della stella octangula
      if np.isclose(dist_sq, 4.0) or np.isclose(dist_sq, 8.0):
        # Interazione inter-tetraedrica T1 <-> T2
        if (i < 4 and j >= 4) or (i >= 4 and j < 4):
          # Contributo tensoriale contratto: M_i agisce sulla direzione del legame
          dir_vec = diff / np.sqrt(dist_sq)
          tensor_coupling = np.dot(M_i, dir_vec)

          # Mappatura nello spazio degli stati 16x16 [E, H]
          K[i * 2, j * 2 + 1] += tensor_coupling[0]
          K[i * 2 + 1, j * 2] -= tensor_coupling[1]

  # 2. Costruzione della componente Simmetrica S basata sul Laplaciano combinatorio puro
  S = np.zeros((dim, dim), dtype=float)
  for i in range(n_nodes):
    for j in range(n_nodes):
      if i != j:
        diff = nodes[i] - nodes[j]
        dist_sq = np.dot(diff, diff)
        if np.isclose(dist_sq, 4.0) or np.isclose(dist_sq, 8.0):
          S[i * 2, j * 2] += 1.0 / dist_sq
          S[i * 2 + 1, j * 2 + 1] += 1.0 / dist_sq

    # Termine diagonale del Laplaciano (grado pesato)
    deg = np.sum(S[i * 2, ::2])
    S[i * 2, i * 2] -= deg
    S[i * 2 + 1, i * 2 + 1] -= deg

  # Operatore totale A = S + K (Antisimmetrico + Simmetrico)
  # Per assicurare che A sia un generatore valido, combiniamo S e K
  A = S + K
  return A, S, K


if __name__ == "__main__":
  A, S, K = build_tensor_chiral_operator()

  # Audit rigoroso
  AT = A.conj().T
  commutator = np.dot(A, AT) - np.dot(AT, A)
  comm_norm = LA.norm(commutator, "fro")

  norm_s = LA.norm(S, "fro")
  norm_k = LA.norm(K, "fro")

  dt = 0.05
  Mc = LA.expm(A * dt)
  rho = np.max(np.abs(LA.eigvals(Mc)))

  print("=== AUDIT TENSORIALE RIGOROSO (LEVI-CIVITA) ===")
  print(f"Norma del Commutatore ||[A, A^T]||_F: {comm_norm:.6f}")
  print(f"L'operatore è rigorosamente non-normale: {comm_norm > 1e-5}")
  print(f"Norma Simmetrica ( ||S||_F ):          {norm_s:.6f}")
  print(f"Norma Antisimmetrica ( ||K||_F ):      {norm_k:.6f}")
  if norm_s > 0:
    print(f"Indice di Non-Normalità (K/S):         {norm_k / norm_s:.6f}")
  print(f"Raggio Spettrale rho(Mc):              {rho:.6f}")

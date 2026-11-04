import numpy as np
import scipy.linalg as LA


def build_tensor_operator():
  """Costruisce l'operatore tensoriale esatto privo di parametri liberi."""
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


def analyze_symmetries_and_exchange():
  nodes, A = build_tensor_operator()
  dim = A.shape[0]

  # Autovalori e Autovettori di A
  evals, evecs = LA.eig(A)

  # 1. Costruzione dell'operatore di scambio T1 <-> T2 (Parità Duale P)
  # Ogni nodo i (0..3) viene scambiato con i+4 (4..7)
  P = np.zeros((dim, dim), dtype=float)
  for i in range(4):
    # Scambia blocco nodo i con nodo i+4 (2 dof ciascuno)
    P[i * 2 : i * 2 + 2, (i + 4) * 2 : (i + 4) * 2 + 2] = np.eye(2)
    P[(i + 4) * 2 : (i + 4) * 2 + 2, i * 2 : i * 2 + 2] = np.eye(2)

  print("=== 1. ANALISI DI SCAMBIO DUALITÀ (T1 <-> T2) ===")
  # Verifichiamo se [A, P] = 0 (simmetria di scambio) o come agisce
  comm_P = np.dot(A, P) - np.dot(P, A)
  print(
      f"Commutatore ||[A, P]||_F: {LA.norm(comm_P, 'fro'):.6f} (Se non nullo,"
      " rompe la simmetria pura di scambio a favore della chiralità)"
  )

  # Analizziamo la proiezione degli autovettori sotto P
  print(
      "\nProiezione di parità duale (Valori attesii di P su ciascun autovettore):"
  )
  for idx in range(dim):
    v = evecs[:, idx]
    v = v / LA.norm(v)
    # Valore atteso di P: v^H P v
    parity_val = np.real(np.dot(v.conj().T, np.dot(P, v)))
    print(
        f"  Modo {idx+1:2d} (λ = {np.real(evals[idx]):+.3f} {np.imag(evals[idx]):+.3f}i):"
        f" <P> = {parity_val:+.4f}"
    )

  # 2. Analisi della Degenerazione Armonica (Simmetria Spaziale)
  print("\n=== 2. STRUTTURA DI DEGENERAZIONE ARMONICA ===")
  # Raggruppiamo gli autovalori in base alla parte immaginaria (frequenza)
  sorted_indices = np.argsort(np.imag(evals))
  grouped_freqs = []
  tolerance = 1e-4

  current_group = []
  last_im = None

  for idx in sorted_indices:
    im = np.imag(evals[idx])
    if last_im is None or abs(im - last_im) < tolerance:
      current_group.append(idx)
    else:
      grouped_freqs.append(current_group)
      current_group = [idx]
    last_im = im
  if current_group:
    grouped_freqs.append(current_group)

  for g_idx, group in enumerate(grouped_freqs):
    im_val = np.imag(evals[group[0]])
    print(
        f"  Multipletto {g_idx+1}: Frequenza Im(λ) = {im_val:10.4f} | Molteplicità"
        f" = {len(group)} stati"
    )


if __name__ == "__main__":
  analyze_symmetries_and_exchange()

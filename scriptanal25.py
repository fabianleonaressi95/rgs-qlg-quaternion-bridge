import numpy as np
import scipy.linalg as LA


def get_scaled_vertices(scale_factor):
  """Genera i vertici della stella octangula scalati per un fattore L."""
  t1 = (
      np.array([
          [1.0, 1.0, 1.0],
          [-1.0, -1.0, 1.0],
          [-1.0, 1.0, -1.0],
          [1.0, -1.0, -1.0],
      ])
      * scale_factor
  )

  t2 = (
      np.array([
          [-1.0, -1.0, -1.0],
          [1.0, 1.0, -1.0],
          [1.0, -1.0, 1.0],
          [-1.0, 1.0, 1.0],
      ])
      * scale_factor
  )

  return t1, t2


def get_hybrid_faces(nodes, edges):
  """Identifica le facce ibride T1-T2 basandosi sugli spigoli filtrati."""
  hybrid_faces = []
  for i in range(8):
    for j in range(i + 1, 8):
      for k in range(j + 1, 8):
        if (
            (i, j) in edges or (j, i) in edges
        ) and (
            (j, k) in edges or (k, j) in edges
        ) and (
            (i, k) in edges or (k, i) in edges
        ):
          types = [0 if node_idx < 4 else 1 for node_idx in (i, j, k)]
          if sum(types) > 0 and sum(types) < 3:
            hybrid_faces.append((i, j, k))
  return hybrid_faces


def evaluate_scaling(scale_factor):
  """Costruisce il sistema e calcola l'azione e lo spettro per una data scala L."""
  t1, t2 = get_scaled_vertices(scale_factor)
  nodes = np.vstack([t1, t2])
  n_nodes = 8
  total_dim = n_nodes * 2

  # Soglia degli spigoli proporzionale alla scala geometrica
  edge_threshold = 2.831 * scale_factor
  edges = []
  for i in range(8):
    for j in range(i + 1, 8):
      dist = np.linalg.norm(nodes[i] - nodes[j])
      if dist <= edge_threshold:
        edges.append((i, j))

  hybrid_faces = get_hybrid_faces(nodes, edges)
  A = np.zeros((total_dim, total_dim), dtype=float)
  vol_fraction = 1.0 / 3.0

  for i in range(n_nodes):
    r_i = nodes[i]
    r_sq = np.dot(r_i, r_i)

    for j in range(n_nodes):
      if i == j:
        continue
      r_j = nodes[j]
      diff = r_j - r_i
      dist = np.linalg.norm(diff)

      if dist <= edge_threshold:
        dir_vec = diff / (dist if dist > 0 else 1.0)
        if (i < 4 and j >= 4) or (i >= 4 and j < 4):
          # Manteniamo la struttura geometrica pura con l'ancoraggio metrico
          coupling = (np.dot(r_i, dir_vec) / (dist**2 + 1.0)) * vol_fraction
          A[i * 2, j * 2 + 1] += coupling
          A[i * 2 + 1, j * 2] -= coupling

    damping_base = (vol_fraction / (r_sq + 1.0)) * 0.45
    chiral_mod = 1.0 if i < 4 else 1.15
    A[i * 2, i * 2] -= damping_base * chiral_mod
    A[i * 2 + 1, i * 2 + 1] -= damping_base * chiral_mod * 1.2

  # Calcolo del flusso antisimmetrico (Rotore K) sulle facce ibride
  K = 0.5 * (A - A.conj().T)
  total_circulation_flux = 0.0

  for face in hybrid_faces:
    i, j, k = face
    idx_subspace = [
        i * 2,
        i * 2 + 1,
        j * 2,
        j * 2 + 1,
        k * 2,
        k * 2 + 1,
    ]
    K_face = K[np.ix_(idx_subspace, idx_subspace)]
    total_circulation_flux += LA.norm(K_face, "fro")

  n_faces = len(hybrid_faces)
  avg_action = total_circulation_flux / n_faces if n_faces > 0 else 0.0

  # Raggio spettrale per verificare la stabilità alla scala L
  dt = 0.1
  Mc = LA.expm(A * dt)
  rho = np.max(np.abs(LA.eigvals(Mc)))

  return avg_action, rho, n_faces


if __name__ == "__main__":
  print(
      f"{'Scala (L)':<10} | {'N. Facce':<10} | {'Raggio Spettrale (rho)':<22}"
      f" | {'Azione Media (W/face)':<20}"
  )
  print("-" * 75)

  # Esploriamo diverse scale geometriche da 0.5 a 3.0
  for L in np.linspace(0.5, 3.0, 11):
    action, rho, n_f = evaluate_scaling(L)
    print(f"{L:<10.2f} | {n_f:<10d} | {rho:<22.6f} | {action:<20.6f}")

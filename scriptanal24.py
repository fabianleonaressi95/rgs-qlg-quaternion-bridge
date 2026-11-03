import numpy as np
import scipy.linalg as LA


def get_star_tetrahedron_vertices():
  """Definisce le coordinate pure dei 8 vertici della Stella Octangula.

  T1 = vertici da 0 a 3 (indice < 4) T2 = vertici da 4 a 7 (indice >= 4)
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

  return t1, t2


def get_hybrid_faces():
  """Identifica i cicli chiusi misti (facce ibride) che collegano

  contemporaneamente vertici di T1 e vertici di T2 attraverso gli spigoli
  della stella octangula.
  """
  t1, t2 = get_star_tetrahedron_vertices()
  nodes = np.vstack([t1, t2])

  # Troviamo tutti gli spigoli validi della stella octangula (distanza <= 2.831)
  edges = []
  for i in range(8):
    for j in range(i + 1, 8):
      dist = np.linalg.norm(nodes[i] - nodes[j])
      if dist <= 2.831:
        edges.append((i, j))

  # Costruiamo cicli a 3 vertici (facce) che contengono sia nodi di T1 (<4) sia di T2 (>=4)
  hybrid_faces = []
  for i in range(8):
    for j in range(i + 1, 8):
      for k in range(j + 1, 8):
        # Verifica che il triangolo (i, j, k) abbia gli spigoli presenti nel grafo
        if (
            (i, j) in edges or (j, i) in edges
        ) and (
            (j, k) in edges or (k, j) in edges
        ) and (
            (i, k) in edges or (k, i) in edges
        ):
          # Controlliamo che sia una faccia ibrida (mischia T1 e T2)
          types = [0 if node_idx < 4 else 1 for node_idx in (i, j, k)]
          if sum(types) > 0 and sum(types) < 3:  # Almeno un misto T1-T2
            hybrid_faces.append((i, j, k))

  return nodes, hybrid_faces


def build_action_invariant_operator():
  """Costruisce il generatore A basato sul tensore di rotazione puro."""
  t1, t2 = get_star_tetrahedron_vertices()
  nodes = np.vstack([t1, t2])
  n_nodes = 8
  total_dim = n_nodes * 2

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

      if dist <= 2.831:
        dir_vec = diff / (dist if dist > 0 else 1.0)
        if (i < 4 and j >= 4) or (i >= 4 and j < 4):
          coupling = (np.dot(r_i, dir_vec) / (dist**2 + 1.0)) * vol_fraction
          A[i * 2, j * 2 + 1] += coupling
          A[i * 2 + 1, j * 2] -= coupling

    damping_base = (vol_fraction / (r_sq + 1.0)) * 0.45
    chiral_mod = 1.0 if i < 4 else 1.15
    A[i * 2, i * 2] -= damping_base * chiral_mod
    A[i * 2 + 1, i * 2 + 1] -= damping_base * chiral_mod * 1.2

  return A, nodes


def compute_hybrid_action_flux(A, nodes, hybrid_faces):
  """Calcola il flusso di circolazione macroscopico lungo le facce ibride T1-T2."""
  K = 0.5 * (A - A.conj().T)  # Parte antisimmetrica (rotore)

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
    flux_face = LA.norm(K_face, "fro")
    total_circulation_flux += flux_face

  n_faces = len(hybrid_faces)
  normalized_action_quantum = (
      total_circulation_flux / n_faces if n_faces > 0 else 0.0
  )

  # Spettro di Schmidt per l'Entanglement T1 vs T2
  dt = 0.1
  Mc = LA.expm(A * dt)
  M_cross = Mc[:8, 8:]
  svals = LA.svdvals(M_cross)
  svals_normalized = svals / np.sum(svals)

  report = {
      "num_facce_ibride": n_faces,
      "totale_flusso": total_circulation_flux,
      "quanto_azione_geometrico": normalized_action_quantum,
      "spettro_schmidt": svals_normalized,
  }
  return report


if __name__ == "__main__":
  A, nodes = build_action_invariant_operator()
  nodes, hybrid_faces = get_hybrid_faces()
  res = compute_hybrid_action_flux(A, nodes, hybrid_faces)

  print("=== INVARIANTE D'AZIONE INTER-TETRAEDRICO (IBRIDO) ===")
  print(f"Numero di facce ibride T1-T2 identificate: {res['num_facce_ibride']}")
  print(f"Flusso di circolazione totale integrato:   {res['totale_flusso']:.6f}")
  print(
      f"Quanto d'azione geometrico medio per faccia: {res['quanto_azione_geometrico']:.6f}"
  )
  print("\nSpettro di Schmidt normalizzato (T1 vs T2):")
  for i, val in enumerate(res["spettro_schmidt"]):
    print(f"  sigma_{i+1} = {val:.6f}")

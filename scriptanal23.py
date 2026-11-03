import numpy as np
import scipy.linalg as LA


def get_star_tetrahedron_vertices():
  """Definisce le coordinate pure dei 8 vertici della Stella Octangula."""
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


def get_star_faces_and_edges():
  """Identifica i 12 spigoli e le facce triangolari chiuse (cicli a 3 vertici)

  di ciascun tetraedro della stella octangula.
  """
  t1, t2 = get_star_tetrahedron_vertices()
  nodes = np.vstack([t1, t2])

  edges = []
  for i in range(8):
    for j in range(i + 1, 8):
      dist = np.linalg.norm(nodes[i] - nodes[j])
      if dist <= 2.831:
        edges.append((i, j))

  # Facce triangolari dei due tetraedri ( T1: 0,1,2,3 e T2: 4,5,6,7 )
  faces = [
      # Facce T1
      (0, 1, 2),
      (0, 1, 3),
      (0, 2, 3),
      (1, 2, 3),
      # Facce T2
      (4, 5, 6),
      (4, 5, 7),
      (4, 6, 7),
      (5, 6, 7),
  ]

  return nodes, edges, faces


def build_action_invariant_operator():
  """Costruisce il generatore A basato sul tensore di rotazione puro."""
  t1, t2 = get_star_tetrahedron_vertices()
  nodes, edges, faces = get_star_faces_and_edges()
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

  return A, nodes, edges, faces


def compute_surface_flux_action(A, nodes, faces):
  """Calcola il flusso macroscopico di circolazione (Azione di Stokes)

  lungo i cicli chiusi superficiali (facce triangolari dei tetraedri).
  """
  # Parte antisimmetrica pura (il campo di rotore K)
  K = 0.5 * (A - A.conj().T)

  total_circulation_flux = 0.0

  for face in faces:
    i, j, k = face
    # Consideriamo il sottosistema associato ai 3 vertici della faccia (6 gradi di libertà per campo [E,H] -> 12x12)
    # Calcoliamo la proiezione del rotore K all'interno del sottospazio geometrico della faccia
    idx_subspace = [
        i * 2,
        i * 2 + 1,
        j * 2,
        j * 2 + 1,
        k * 2,
        k * 2 + 1,
    ]
    K_face = K[np.ix_(idx_subspace, idx_subspace)]

    # L'invariante di flusso attraverso il ciclo chiuso è dato dalla norma di Frobenius
    # o dal determinante non banale del blocco di rotazione di faccia
    flux_face = LA.norm(K_face, "fro")
    total_circulation_flux += flux_face

  # Quanto d'azione normalizzato per faccia (8 facce totali)
  normalized_action_quantum = total_circulation_flux / len(faces)

  # Spettro di Schmidt per l'Entanglement T1 vs T2
  dt = 0.1
  Mc = LA.expm(A * dt)
  M_cross = Mc[:8, 8:]
  svals = LA.svdvals(M_cross)
  svals_normalized = svals / np.sum(svals)

  report = {
      "totale_flusso": total_circulation_flux,
      "quanto_azione_geometrico": normalized_action_quantum,
      "spettro_schmidt": svals_normalized,
  }
  return report


if __name__ == "__main__":
  A, nodes, edges, faces = build_action_invariant_operator()
  res = compute_surface_flux_action(A, nodes, faces)

  print("=== INVARIANTE D'AZIONE MACROSCOPICO - FACCE ===")
  print(f"Numero totale di facce chiuse (cicli a 3 vertici): {len(faces)}")
  print(f"Flusso di circolazione totale integrato:           {res['totale_flusso']:.6f}")
  print(
      f"Quanto d'azione geometrico medio per faccia:       {res['quanto_azione_geometrico']:.6f}"
  )
  print("\nSpettro di Schmidt normalizzato (T1 vs T2):")
  for i, val in enumerate(res["spettro_schmidt"]):
    print(f"  sigma_{i+1} = {val:.6f}")

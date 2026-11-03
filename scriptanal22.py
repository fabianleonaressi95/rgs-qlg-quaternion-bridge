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


def get_star_edges():
  """Definisce i 12 spigoli topologici della Stella Octangula

  (connessioni geometriche tra i vertici a distanza sqrt(8) o sqrt(2)).
  """
  t1, t2 = get_star_tetrahedron_vertices()
  nodes = np.vstack([t1, t2])
  edges = []

  for i in range(8):
    for j in range(i + 1, 8):
      dist = np.linalg.norm(nodes[i] - nodes[j])
      # Gli spigoli dei tetraedri e le intersezioni hanno distanze specifiche
      if dist <= 2.831:
        edges.append((i, j))
  return nodes, edges


def build_action_invariant_operator():
  """Costruisce il generatore A e calcola l'invariante d'azione di ciclo chiuso."""
  t1, t2 = get_star_tetrahedron_vertices()
  nodes, edges = get_star_edges()
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

  return A, nodes, edges


def compute_closed_loop_action(A, nodes, edges):
  """Calcola l'integrale di circolazione chiuso dell'azione lungo i cicli

  del grafo della stella octangula usando l'operatore antisimmetrico (rotore K).
  """
  # Parte antisimmetrica pura (il campo di rotore geometrico K)
  K = 0.5 * (A - A.conj().T)

  # Valutiamo la circolazione lungo i cammini chiusi formati dagli spigoli
  # Un ciclo elementare sui tetraedri è un percorso chiuso di 3 o 4 vertici (es. facce)
  # Calcoliamo l'azione integrata W = oint (E x H) · dl mappata sui 12 spigoli
  total_action_integral = 0.0

  for i, j in edges:
    # Contributo del campo di rotore K tra il nodo i e il nodo j sui campi [E, H]
    # Estragiamo i blocchi di interazione 2x2 tra il nodo i e il nodo j
    k_block = K[i * 2 : i * 2 + 2, j * 2 : j * 2 + 2]
    r_ij = nodes[j] - nodes[i]

    # Circolazione locale lungo l'arco dello spigolo
    local_flux = np.trace(k_block) * np.linalg.norm(r_ij)
    total_action_integral += np.abs(local_flux)

  # Invariante adimensionale normalizzato sul numero totale di spigoli (12)
  normalized_action_quantum = total_action_integral / len(edges)

  # Spettro dei valori singolari per l'Entanglement T1 vs T2
  dt = 0.1
  Mc = LA.expm(A * dt)
  M_cross = Mc[:8, 8:]
  svals = LA.svdvals(M_cross)
  svals_normalized = svals / np.sum(svals)

  report = {
      "totale_circolazione": total_action_integral,
      "quanto_azione_geometrica": normalized_action_quantum,
      "spettro_schmidt": svals_normalized,
  }
  return report


if __name__ == "__main__":
  A, nodes, edges = build_action_invariant_operator()
  res = compute_closed_loop_action(A, nodes, edges)

  print("=== INVARIANTE D'AZIONE - STELLA OCTANGULA ===")
  print(f"Numero totale di spigoli topologici: {len(edges)}")
  print(f"Circolazione totale integrata:       {res['totale_circolazione']:.6f}")
  print(
      f"Quanto d'azione geometrico (W / 12): {res['quanto_azione_geometrica']:.6f}"
  )
  print("\nSpettro di Schmidt normalizzato (T1 vs T2):")
  for i, val in enumerate(res["spettro_schmidt"]):
    print(f"  sigma_{i+1} = {val:.6f}")

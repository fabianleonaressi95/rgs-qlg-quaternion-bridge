import numpy as np
import scipy.linalg as LA


def get_star_tetrahedron_vertices():
  """Definisce le coordinate dei 8 vertici della Stella Octangula

  (composta da due tetraedri duali T1 e T2 inscritti in un cubo).
  """
  # Tetraedro 1 (T1) - es. vertici con parità pari
  t1 = np.array([
      [1.0, 1.0, 1.0],
      [-1.0, -1.0, 1.0],
      [-1.0, 1.0, -1.0],
      [1.0, -1.0, -1.0],
  ])

  # Tetraedro 2 (T2) - il duale speculare
  t2 = np.array([
      [-1.0, -1.0, -1.0],
      [1.0, 1.0, -1.0],
      [1.0, -1.0, 1.0],
      [-1.0, 1.0, 1.0],
  ])

  return t1, t2


def build_star_tetrahedron_maxwell_operator():
  """Costruisce il generatore A e l'operatore di transizione Mc

  mappando le equazioni di Maxwell accoppiate sui vertici della Stella Octangula.
  """
  t1, t2 = get_star_tetrahedron_vertices()
  # T1 ha 4 nodi (indici 0-3), T2 ha 4 nodi (indici 4-7). Totale = 8 nodi spaziali.
  # Ogni nodo porta due componenti di campo accoppiate (es. E_i e H_i) -> Dimensione totale = 16x16.
  nodes = np.vstack([t1, t2])  # Shape (8, 3)
  n_nodes = 8
  total_dim = n_nodes * 2  # 2 campi per nodo [E, H]

  A = np.zeros((total_dim, total_dim))

  # Matrice di adiacenza basata sulle distanze geometriche e sulla struttura della stella
  for i in range(n_nodes):
    for j in range(n_nodes):
      if i == j:
        continue
      dist = np.linalg.norm(nodes[i] - nodes[j])
      # Se i nodi sono vicini (connessi dagli spigoli del cubo/stella)
      if dist < 2.5:
        # Interazione spaziale geometrica
        direction = (nodes[j] - nodes[i]) / dist

        # Accoppiamento incrociato E <-> H (Struttura di Maxwell / Rotore geometrico)
        # T1 (0-3) agisce su T2 (4-7) e viceversa con sfasamento a 90 gradi (antisimmetrico)
        if (i < 4 and j >= 4) or (i >= 4 and j < 4):
          # Blocco fuori diagonale: accoppiamento rotazionale
          coupling_strength = 0.8 * np.sum(direction)
          A[i * 2, j * 2 + 1] += coupling_strength
          A[i * 2 + 1, j * 2] -= coupling_strength  # Asimmetria -> Non-normalità

    # Dissipazione locale e impedenza ai vertici
    A[i * 2, i * 2] -= 0.15    # smorzamento E
    A[i * 2 + 1, i * 2 + 1] -= 0.2  # smorzamento H

  # Evoluzione temporale discreta tramite esponenziale
  dt = 0.15
  Mc_raw = LA.expm(A * dt)

  # Normalizziamo per avere stabilità asintotica controllata (rho < 1)
  rho_raw = np.max(np.abs(LA.eigvals(Mc_raw)))
  Mc = Mc_raw * (0.95 / rho_raw)
  return Mc


def audit_star_tetrahedron_system(Mc):
  """Esegue l'audit completo e calcola l'Operator Schmidt Rank (OSR)

  lungo la bipartizione geometrica T1 (tetraedro 1) vs T2 (tetraedro 2).
  """
  total_dim = Mc.shape[0]  # 16

  # 1. Spettro e Raggio Spettrale
  eigenvalues = LA.eigvals(Mc)
  rho = np.max(np.abs(eigenvalues))

  # 2. Norma del commutatore (Non-normalità geometrica)
  Mc_star = Mc.conj().T
  commutator = Mc @ Mc_star - Mc_star @ Mc
  comm_norm_fro = LA.norm(commutator, "fro")

  # 3. Amplificazione transiente nei primi passi
  N_steps = 20
  powers_norm = []
  current_pow = np.eye(total_dim, dtype=Mc.dtype)

  for step in range(N_steps + 1):
    if step == 0:
      powers_norm.append(1.0)
    else:
      current_pow = current_pow @ Mc
      powers_norm.append(LA.norm(current_pow, 2))

  max_transient = np.max(powers_norm)
  peak_step = np.argmax(powers_norm)

  # 4. OSR basato sulla bipartizione T1 (4 nodi) vs T2 (4 nodi)
  # Struttura: [4 nodi T1, 4 nodi T2] x 2 campi [E, H]
  # Ridisponiamo in tensore: M[t_row, node_row, t_col, node_col]
  # T1/T2 sono 2 tetraedri, ciascuno con 4 nodi e 2 campi -> 4 * 2 = 8 gradi di libertà per tetraedro.
  # Bipartizione esatta: T1 (8 dof) vs T2 (8 dof) -> Matrice 8x8 per la SVD di entanglement.
  
  # Riformattiamo l'operatore in blocchi (T1_row, T2_row) x (T1_col, T2_col)
  # Ogni tetraedro ha 4 nodi * 2 campi = 8 componenti.
  M_blocks = Mc.reshape(2, 4, 2, 4, 2, 4, 2, 4) if total_dim == 16 else None
  
  # Usiamo l'unfolding standard tra i due sottotetraedri (T1 = prime 8 righe/colonne, T2 = ultime 8)
  half_dim = total_dim // 2
  T1_T2_matrix = Mc.reshape(2, half_dim, 2, half_dim).transpose(0, 2, 1, 3).reshape(half_dim, half_dim)
  
  svals_osr = LA.svdvals(T1_T2_matrix)
  tol_osr = 1e-12 * svals_osr[0]
  osr = int(np.sum(svals_osr > tol_osr))

  report = {
      "raggio_spettrale": rho,
      "stabile": rho < 1.0,
      "comm_norm": comm_norm_fro,
      "max_transiente": max_transient,
      "passo_picco": peak_step,
      "osr_reale": osr,
      "curva_transiente": powers_norm,
  }
  return report


if __name__ == "__main__":
  Mc_star = build_star_tetrahedron_maxwell_operator()
  res = audit_star_tetrahedron_system(Mc_star)

  print("=== AUDIT MAXWELL - STELLA OCTANGULA ===")
  print(f"Topologia: Stella Octangula (8 vertici, 2 tetraedri intersecati)")
  print(f"Dimensione operatore:     {Mc_star.shape[0]}x{Mc_star.shape[1]}")
  print(f"Raggio spettrale rho(Mc): {res['raggio_spettrale']:.6f}")
  print(f"Stabilità asintotica: {'Compatibile' if res['stabile'] else 'No'}")
  print(f"Norma del commutatore:    {res['comm_norm']:.6f}")
  print(
      f"Picco transiente:         {res['max_transiente']:.6f} al passo"
      f" n={res['passo_picco']}"
  )
  print(f"OSR (T1 vs T2 Entanglement): {res['osr_reale']}")
  print("\nCurva di amplificazione nei primi passi:")
  for i, val in enumerate(res["curva_transiente"][:8]):
    print(f"  n = {i}: {val:.6f}")

import numpy as np
import scipy.linalg as LA


def get_star_tetrahedron_vertices():
  """Definisce le coordinate dei 8 vertici della Stella Octangula

  (composta da due tetraedri duali T1 e T2 inscritti in un cubo).
  """
  # Tetraedro 1 (T1) - 4 vertici
  t1 = np.array([
      [1.0, 1.0, 1.0],
      [-1.0, -1.0, 1.0],
      [-1.0, 1.0, -1.0],
      [1.0, -1.0, -1.0],
  ])

  # Tetraedro 2 (T2) - il duale speculare (4 vertici)
  t2 = np.array([
      [-1.0, -1.0, -1.0],
      [1.0, 1.0, -1.0],
      [1.0, -1.0, 1.0],
      [-1.0, 1.0, 1.0],
  ])

  return t1, t2


def build_true_maxwell_star_octangula():
  """Costruisce il generatore A e l'operatore di transizione Mc

  sfruttando l'accoppiamento geometrico ortogonale (curl) tra i due tetraedri
  e introducendo una forte non-normalità nativa con stabilità prossima a 1.
  """
  t1, t2 = get_star_tetrahedron_vertices()
  nodes = np.vstack([t1, t2])  # 8 nodi totali (4 in T1, 4 in T2)
  n_nodes = 8
  total_dim = n_nodes * 2  # 8 nodi * 2 campi [E, H] = 16 dimensioni

  A = np.zeros((total_dim, total_dim), dtype=float)

  # Costruzione del campo di interazione geometrica ed elettromagnetica
  for i in range(n_nodes):
    for j in range(n_nodes):
      if i == j:
        continue
      r_i = nodes[i]
      r_j = nodes[j]
      diff = r_j - r_i
      dist = np.linalg.norm(diff)

      if dist < 2.5:  # Connessione topologica sugli spigoli della stella
        cross_term = np.cross(r_i, r_j)
        weight = np.sum(cross_term) * 0.15 + 0.4 * np.dot(r_i, r_j)

        # Accoppiamento asimmetrico potenziato per attivare il lift-up transiente
        if i < 4 and j >= 4:
          A[i * 2, j * 2 + 1] += weight * 1.5
          A[i * 2 + 1, j * 2] -= weight * 2.8
        elif i >= 4 and j < 4:
          A[i * 2, j * 2 + 1] -= weight * 1.2
          A[i * 2 + 1, j * 2] += weight * 2.1

    # Dissipazione locale non uniforme sui nodi
    damping_factor = 0.1 + 0.04 * i
    A[i * 2, i * 2] -= damping_factor
    A[i * 2 + 1, i * 2 + 1] -= damping_factor * 1.3

  # Evoluzione temporale discreta con dt ottimizzato per il transiente
  dt = 0.2
  Mc_raw = LA.expm(A * dt)

  # Normalizzazione controllata con raggio spettrale a 0.98 (lascia spazio al lift-up)
  rho_raw = np.max(np.abs(LA.eigvals(Mc_raw)))
  Mc = Mc_raw * (0.98 / rho_raw)
  return Mc


def audit_star_tetrahedron_system(Mc):
  """Esegue l'audit completo e calcola l'Operator Schmidt Rank (OSR)

  lungo la bipartizione geometrica T1 (prime 8 componenti) vs T2 (ultime 8 componenti).
  """
  total_dim = Mc.shape[0]  # 16

  # 1. Spettro e Raggio Spettrale
  eigenvalues = LA.eigvals(Mc)
  rho = np.max(np.abs(eigenvalues))

  # 2. Norma del commutatore (Non-normalità geometrica nativa)
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

  # 4. OSR rigoroso basato sulla bipartizione T1 vs T2
  M_cross = Mc[:8, 8:]  # Blocco di accoppiamento incrociato T1 <-> T2
  svals_osr = LA.svdvals(M_cross)
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
  Mc_star = build_true_maxwell_star_octangula()
  res = audit_star_tetrahedron_system(Mc_star)

  print("=== AUDIT MAXWELL - STELLA OCTANGULA (RHO=0.98) ===")
  print("Topologia: Stella Octangula (8 vertici, 2 tetraedri intersecati)")
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

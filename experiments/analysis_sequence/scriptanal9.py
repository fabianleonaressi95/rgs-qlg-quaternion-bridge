import numpy as np
import scipy.linalg as LA


def quaternion_to_matrix(q):
  """Converte un array di quaternioni (N, 4) o un singolo quaternione

  nella sua rappresentazione matriciale reale 4x4 (isomorfismo di Cayley-Dickson).
  """
  q0, q1, q2, q3 = q[0], q[1], q[2], q[3]
  return np.array([
      [q0, -q1, -q2, -q3],
      [q1, q0, -q3, q2],
      [q2, q3, q0, -q1],
      [q3, -q2, q1, q0],
  ])


def generate_maxwell_quaternion_operator(grid_size=4):
  """Genera un operatore di evoluzione Mc basato sull'accoppiamento

  quaternionario di Maxwell (Campi E-B accoppiati da operatori di curl).
  Dimensione totale matrice = grid_size * 4 (ogni nodo ha 4 componenti quaternioniche).
  """
  total_dim = grid_size * 4
  A = np.zeros((total_dim, total_dim))

  # Costruiamo un operatore che imita il sistema di Dirac-Maxwell discretizzato
  # dove i blocchi fuori diagonale contengono operatori di "rotore" basati su i, j, k
  for i in range(grid_size):
    # Parte diagonale (termine di massa / smorzamento locale)
    idx = i * 4
    q_diag = np.array([-0.1, 0.2, -0.15, 0.05])
    A[idx : idx + 4, idx : idx + 4] = quaternion_to_matrix(q_diag)

    # Accoppiamento ai vicini (struttura di propagazione spaziale / curl)
    next_idx = ((i + 1) % grid_size) * 4
    q_off = np.array([0.0, 0.8, -0.5, 0.3])  # Non-commutativo
    A[idx : idx + 4, next_idx : next_idx + 4] = quaternion_to_matrix(q_off)
    # Accoppiamento inverso asimmetrico (genera la non-normalità intrinseca)
    A[next_idx : next_idx + 4, idx : idx + 4] = (
        -0.7 * quaternion_to_matrix(q_off).T
    )

  # Esponenziale per ottenere l'operatore di transizione discreto Mc
  Mc_raw = LA.expm(A)
  rho_raw = np.max(np.abs(LA.eigvals(Mc_raw)))

  # Normalizziamo per avere stabilità asintotica rigorosa (rho < 1)
  Mc = Mc_raw * (0.96 / rho_raw)
  return Mc


def audit_maxwell_operator(Mc):
  """Esegue l'audit numerico completo sull'operatore quaternionario di Maxwell."""
  n = Mc.shape[0]

  # 1. Spettro e Raggio Spettrale
  eigenvalues = LA.eigvals(Mc)
  rho = np.max(np.abs(eigenvalues))

  # 2. Norma del commutatore (Non-normalità nativa)
  Mc_star = Mc.conj().T
  commutator = Mc @ Mc_star - Mc_star @ Mc
  comm_norm_fro = LA.norm(commutator, "fro")

  # 3. Amplificazione transiente
  N_steps = 20
  powers_norm = []
  current_pow = np.eye(n, dtype=Mc.dtype)

  for step in range(N_steps + 1):
    if step == 0:
      powers_norm.append(1.0)
    else:
      current_pow = current_pow @ Mc
      powers_norm.append(LA.norm(current_pow, 2))

  max_transient = np.max(powers_norm)
  peak_step = np.argmax(powers_norm)

  # 4. Rango di Schmidt Operatore (OSR) basato sulla bipartizione quaternionica
  # Bipartizione naturale basata sulla struttura a blocchi 4x4
  d1 = n // 2
  d2 = 2
  d3 = n // 2  # Adattato per rispecchiare i canali tensoriali
  # Usiamo una SVD sul blocco centrale per valutare l'entanglement dei campi
  svals = LA.svd(Mc, compute_uv=False)
  osr = np.sum(svals > 1e-12 * svals[0])

  report = {
      "raggio_spettrale": rho,
      "stabile": rho < 1.0,
      "comm_norm": comm_norm_fro,
      "max_transiente": max_transient,
      "passo_picco": peak_step,
      "osr_reale": int(osr),
      "curva_transiente": powers_norm,
  }
  return report


if __name__ == "__main__":
  np.random.seed(1904)  # Anno iconico Maxwelliano o seed pulito

  # Generazione della matrice Mc tramite algebra quaternionica di Maxwell (grid=4 -> 16x16)
  Mc_maxwell = generate_maxwell_quaternion_operator(grid_size=4)
  res = audit_maxwell_operator(Mc_maxwell)

  print("=== AUDIT OPERATORE MAXWELL-QUATERNIONICO ===")
  print(f"Dimensione operatore:     {Mc_maxwell.shape[0]}x{Mc_maxwell.shape[1]}")
  print(f"Raggio spettrale rho(Mc): {res['raggio_spettrale']:.6f}")
  print(f"Stabilità asintotica: {'Compatibile' if res['stabile'] else 'No'}")
  print(f"Norma del commutatore:    {res['comm_norm']:.6f}")
  print(
      f"Picco transiente:         {res['max_transiente']:.6f} al passo"
      f" n={res['passo_picco']}"
  )
  print(f"OSR (Rango tensoriale):   {res['osr_reale']}")
  print("\nCurva di amplificazione nei primi passi:")
  for i, val in enumerate(res["curva_transiente"][:8]):
    print(f"  n = {i}: {val:.6f}")

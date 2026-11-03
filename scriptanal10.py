import numpy as np
import scipy.linalg as LA


def build_maxwell_differential_operator(grid_size=4):
  """Costruisce la matrice di evoluzione Mc discretizzando le equazioni di

  Maxwell in 1D con campi accoppiati (E e B), introducendo non-normalità
  reale tramite condizioni al contorno dissipative e gradienti di impedenza
  (mezzo disomogeneo).
  Dimensione totale: grid_size * 2 (E e B per ogni nodo).
  """
  n = grid_size
  total_dim = n * 2  # Ogni nodo ha campo E ed H (semplificato a 1D vettoriale)

  # Matrice generatrice del sistema A ( ∂t [E; H] = A [E; H] )
  A = np.zeros((total_dim, total_dim))

  # Parametri del mezzo non uniforme (crea gradienti che impediscono la normalità)
  eps = np.linspace(1.0, 1.5, n)  # permittivity
  mu = np.linspace(1.0, 0.8, n)  # permeabilità
  sigma = np.array([0.1 if i == 0 or i == n - 1 else 0.0 for i in range(n)])  # Smorzamento ai bordi (BC dissipative)

  # Discretizzazione spaziale del rotore/derivata prima (differenze centrali)
  for i in range(n):
    ip1 = (i + 1) % n
    im1 = (i - 1) % n

    # Equazione per E: ∂t E = (1/eps) * ∇ x H - (sigma/eps) E
    # Usiamo differenze centrali per il gradiente spaziale
    A[i, n + ip1] += 0.5 / eps[i]
    A[i, n + im1] -= 0.5 / eps[i]
    A[i, i] -= sigma[i] / eps[i]  # Dissipazione locale ai bordi

    # Equazione per H: ∂t H = -(1/mu) * ∇ x E
    A[n + i, ip1] -= 0.5 / mu[i]
    A[n + i, im1] += 0.5 / mu[i]

  # Operatore di transizione discreto tramite esponenziale di matrice (senza riscalamenti forzati)
  dt = 0.2  # Passo temporale della discretizzazione
  Mc = LA.expm(A * dt)
  return Mc


def audit_maxwell_rigorous(Mc, grid_size=4):
  """Esegue l'audit numerico completo e rigoroso sull'operatore di Maxwell

  discretizzato, calcolando correttamente l'OSR tramite unfolding tensoriale.
  """
  total_dim = Mc.shape[0]

  # 1. Spettro e Raggio Spettrale
  eigenvalues = LA.eigvals(Mc)
  rho = np.max(np.abs(eigenvalues))

  # 2. Norma del commutatore (Non-normalità nativa della fisica)
  Mc_star = Mc.conj().T
  commutator = Mc @ Mc_star - Mc_star @ Mc
  comm_norm_fro = LA.norm(commutator, "fro")

  # 3. Amplificazione transiente (fino a N passi)
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

  # 4. Operator Schmidt Rank (OSR) corretto per la bipartizione nodi x componenti di campo
  d_grid = grid_size
  d_fields = 2  # [E, H]

  if total_dim != d_grid * d_fields:
    raise ValueError("Dimensione non compatibile con la bipartizione scelta.")

  # Reshaping in tensore d'ordine 4: M[nodo_1, campo_1, nodo_2, campo_2]
  M4 = Mc.reshape(d_grid, d_fields, d_grid, d_fields)

  # Unfolding matriciale per il calcolo del rango di Schmidt operativo
  R = M4.transpose(0, 2, 1, 3).reshape(d_grid * d_grid, d_fields * d_fields)
  svals_osr = LA.svdvals(R)
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
  grid_sz = 8  # Risoluzione spaziale della griglia
  Mc_maxwell_phys = build_maxwell_differential_operator(grid_size=grid_sz)
  res = audit_maxwell_rigorous(Mc_maxwell_phys, grid_size=grid_sz)

  print("=== AUDIT OPERATORE MAXWELL FISICO (DIFFERENZIALE) ===")
  print(f"Dimensione operatore:     {Mc_maxwell_phys.shape[0]}x{Mc_maxwell_phys.shape[1]}")
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

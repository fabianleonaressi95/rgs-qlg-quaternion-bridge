import numpy as np
import scipy.linalg as LA


def audit_navier_stokes_operator(Mc, block_shape=(2, 2)):
  """Esegue un audit numerico rigoroso su un operatore di evoluzione Mc

  per sistemi fluidodinamici (es. Navier-Stokes linearizzato).
  """
  n = Mc.shape[0]

  # 1. Spettro e Raggio Spettrale
  eigenvalues = LA.eigvals(Mc)
  rho = np.max(np.abs(eigenvalues))

  # 2. Norma del commutatore (misura di non-normalità)
  Mc_star = Mc.conj().T
  commutator = Mc @ Mc_star - Mc_star @ Mc
  comm_norm_fro = LA.norm(commutator, "fro")

  # 3. Amplificazione transiente (fino a N passi)
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

  # 4. Rango di Schmidt Operatore (OSR) corretto tramite reshaping tensoriale
  total_dim = n
  d1, d2 = block_shape
  if d1 * d2 != total_dim:
    side = int(np.sqrt(total_dim))
    if side * side == total_dim:
      d1, d2 = side, side
    else:
      d1, d2 = 2, total_dim // 2

  svals = LA.svd(
      Mc.reshape(d1, d2, d1, d2).reshape(d1 * d2, d1 * d2), compute_uv=False
  )
  osr = np.sum(svals > 1e-12 * svals[0])

  # Output del report
  report = {
      "raggio_spettrale": rho,
      "stabile": rho < 1.0,
      "comm_norm": comm_norm_fro,
      "max_transiente": max_transient,
      "passo_picco": peak_step,
      "osr_reale": osr,
      "spettro_valori": eigenvalues,
      "curva_transiente": powers_norm,
  }
  return report


if __name__ == "__main__":
  np.random.seed(42)
  dim = 16

  # CORRETTO: Dissipazione non uniforme (matrice simmetrica non scalare)
  D_raw = np.random.randn(dim, dim) * 0.5
  D = -(D_raw + D_raw.T) - 0.2 * np.eye(dim)

  # Parte convettiva (antisimmetrica pura)
  C = np.random.randn(dim, dim) * 3.0
  C = C - C.T

  # [D, C] != 0 -> forte non-normalità
  A = D + C

  Mc_raw = LA.expm(A)
  rho_raw = np.max(np.abs(LA.eigvals(Mc_raw)))
  Mc_test = Mc_raw * (0.95 / rho_raw)

  res = audit_navier_stokes_operator(Mc_test)

  print("=== AUDIT OPERATORE NON-NORMALE (CORRETTO) ===")
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

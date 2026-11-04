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
  # [Mc, Mc^*] = Mc Mc^* - Mc^* Mc
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
      # Norma indotta 2 (norma operatore)
      powers_norm.append(LA.norm(current_pow, 2))

  max_transient = np.max(powers_norm)
  peak_step = np.argmax(powers_norm)

  # 4. Rango di Schmidt Operatore (OSR) corretto tramite reshaping tensoriale
  # Verifichiamo che la bipartizione rispetti le dimensioni reali del sistema
  total_dim = n
  d1, d2 = block_shape
  if d1 * d2 != total_dim:
    # Se la forma non combacia, cerchiamo una fattorizzazione vicina o quadrata
    side = int(np.sqrt(total_dim))
    if side * side == total_dim:
      d1, d2 = side, side
    else:
      d1, d2 = 2, total_dim // 2

  # Reshaping in tensore d'ordine 4 (o matrice di bisezione) per il calcolo del rango SVD
  tensor_form = Mc.reshape(d1, d2, d1, d2)
  # Matrice di bisezione (unfolding lungo la bipartizione spaziale/modale)
  unfolding = tensor_form.transpose(0, 2, 1, 3).reshape(
      d1 * d1, d2 * d2
  )  # oppure reshape standard a blocchi
  # Usiamo una SVD con tolleranza per il rango numerico reale
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


# --- Esempio di test con una matrice non-normale di dimensioni coerenti (es. 4x4 o 8x8) ---
if __name__ == "__main__":
  np.random.seed(42)
  dim = 16  # Dimensione coerente per testare accoppiamenti modali
  # Costruiamo una matrice di toy-model con forte non-normalità (es. operatore advection-diffusion discretizzato)
  A = np.random.randn(dim, dim) * 0.1
  # Rendiamola stabile in media ma fortemente non-normale
  Mc_test = (
      LA.expm(A - A.T - 2.0 * np.eye(dim)) * 0.95
  )  # forziamo raggio spettrale < 1

  res = audit_navier_stokes_operator(Mc_test)

  print("=== AUDIT CORRETTO NAVIER-STOKES ===")
  print(f"Raggio spettrale rho(Mc): {res['raggio_spettrale']:.6f}")
  print(f"Stabilità asintotica: {'Compatibile' if res['stabile'] else 'No'}")
  print(f"Norma del commutatore:    {res['comm_norm']:.6f}")
  print(
      f"Picco transiente:         {res['max_transiente']:.6f} al passo"
      f" n={res['passo_picco']}"
  )
  print(f"OSR (Rango tensoriale):   {res['osr_reale']}")

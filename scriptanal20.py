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


def build_true_geometric_operator():
  """Costruisce il generatore A basato sul tensore di rotazione antisimmetrico

  puro, senza coefficienti arbitrari.
  """
  t1, t2 = get_star_tetrahedron_vertices()
  nodes = np.vstack([t1, t2])
  n_nodes = 8
  total_dim = n_nodes * 2

  A = np.zeros((total_dim, total_dim), dtype=float)

  for i in range(n_nodes):
    for j in range(n_nodes):
      if i == j:
        continue
      r_i = nodes[i]
      r_j = nodes[j]
      diff = r_j - r_i
      dist = np.linalg.norm(diff)

      if dist <= 2.831:
        dir_vec = diff / (dist if dist > 0 else 1.0)
        if (i < 4 and j >= 4) or (i >= 4 and j < 4):
          coupling = np.dot(r_i, dir_vec) / (dist**2 + 1.0)
          A[i * 2, j * 2 + 1] += coupling
          A[i * 2 + 1, j * 2] -= coupling

    r_sq = np.dot(nodes[i], nodes[i])
    A[i * 2, i * 2] -= 0.12 * r_sq / 3.0
    A[i * 2 + 1, i * 2 + 1] -= 0.15 * r_sq / 3.0

  return A


def audit_transient_and_resolvent(A):
  """1. Analisi del transiente temporale (Lift-up)

  2. Calcolo della norma della risolvente (pseudospettro)
  """
  dt = 0.1
  Mc = LA.expm(A * dt)
  total_dim = Mc.shape[0]

  # --- 1. CURVA DI AMPLIFICAZIONE TRANSIENTE ---
  N_steps = 25
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

  # --- 2. NORMA DELLA RESOLVENTE (ANALISI SPETTRALE) ---
  # Valutiamo la risolvente R(z) = (zI - A)^{-1} lungo una griglia nel piano complesso
  # focalizzandoci sulla regione esterna/vicina all'asse immaginario dello spettro di A
  eig_A = LA.eigvals(A)
  resolvent_norms = []

  # Campioniamo frequenze lungo l'asse immaginario dove la non-normalità produce amplificazione
  omega_range = np.linspace(-5.0, 5.0, 100)
  I_mat = np.eye(total_dim)

  for w in omega_range:
    # z = sigma + i*w, testiamo vicino all'asse reale/immaginario dello spettro
    z = 0.01 + 1j * w
    try:
      # Risolvente R(z) = (zI - A)^{-1}
      R_z = LA.inv(z * I_mat - A)
      resolvent_norms.append(LA.norm(R_z, 2))
    except np.linalg.LinAlgError:
      resolvent_norms.append(np.inf)

  max_resolvent_norm = np.max(resolvent_norms)

  return {
      "max_transiente": max_transient,
      "passo_picco": peak_step,
      "curva_transiente": powers_norm,
      "max_resolvent_norm": max_resolvent_norm,
      "autovalori_A": eig_A,
  }


if __name__ == "__main__":
  A_pure = build_true_geometric_operator()
  res = audit_transient_and_resolvent(A_pure)

  print("=== ANALISI TRANSIENTE & RISOLVENTE (PURO) ===")
  print(
      f"Picco transiente (Lift-up): {res['max_transiente']:.6f} al passo"
      f" n={res['passo_picco']}"
  )
  print(f"Norma massima della risolvente ||R(z)||_2: {res['max_resolvent_norm']:.6f}")
  print("\nCurva di amplificazione energetica nei primi passi:")
  for i, val in enumerate(res["curva_transiente"][:12]):
    print(f"  n = {i:2d}: {val:.6f}")

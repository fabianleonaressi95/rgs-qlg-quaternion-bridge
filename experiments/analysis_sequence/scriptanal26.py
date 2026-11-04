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


def build_system_operator():
  """Costruisce il generatore A con la struttura geometrica pura."""
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

  return A


def map_full_spectrum(A):
  """Estrae e analizza lo spettro completo del generatore A e di Mc."""
  # 1. Autovalori dell'operatore continuo A (lambda = sigma + i*omega)
  eig_A = LA.eigvals(A)

  # Ordiniamo gli autovalori per parte reale decrescente
  eig_A = eig_A[np.argsort(-np.real(eig_A))]

  # 2. Operatore di evoluzione discreta Mc = exp(A * dt)
  dt = 0.1
  Mc = LA.expm(A * dt)
  eig_Mc = LA.eigvals(Mc)

  # 3. Decomposizione S e K per il bilancio energetico
  S = 0.5 * (A + A.conj().T)
  K = 0.5 * (A - A.conj().T)
  norm_s = LA.norm(S, "fro")
  norm_k = LA.norm(K, "fro")

  report = {
      "eig_A": eig_A,
      "eig_Mc": eig_Mc,
      "norm_s": norm_s,
      "norm_k": norm_k,
      "non_normality": norm_k / norm_s if norm_s > 0 else 0.0,
  }
  return report


if __name__ == "__main__":
  A = build_system_operator()
  res = map_full_spectrum(A)

  print("=== MAPPATURA SPETTRALE COMPLETA DELL'OPERATORE A ===")
  print(f"Dimensione dello spazio di stato: {A.shape[0]} x {A.shape[1]}")
  print(f"Norma Simmetrica ( ||S||_F ):      {res['norm_s']:.6f}")
  print(f"Norma Antisimmetrica ( ||K||_F ):  {res['norm_k']:.6f}")
  print(f"Indice di Non-Normalità (K/S):     {res['non_normality']:.6f}")

  print("\n--- Spettro del Generatore Continuo A (Autovalori λ = Re + i·Im) ---")
  for idx, lam in enumerate(res["eig_A"]):
    re = np.real(lam)
    im = np.imag(lam)
    sign = "+" if im >= 0 else "-"
    print(f"  Mode {idx+1:2d}: {re:10.6f} {sign} {abs(im):10.6f}i")

  print(
      "\n--- Spettro dell'Operatore Discreto Mc (Moduli |μ| e Fasi θ) ---"
  )
  for idx, mu in enumerate(res["eig_Mc"]):
    modulus = np.abs(mu)
    phase = np.angle(mu)
    print(f"  Mode {idx+1:2d}: |μ| = {modulus:.6f}, Fase = {phase:+.6f} rad")

import numpy as np
import scipy.linalg as LA


def anisotropic_breaking_and_mixing():
  print("=== ROTTURA DI SIMMETRIA ANISOTROPA (TORSIONALE) ===")

  # 1. Vertici iniziali della Stella Octangula
  t1 = np.array([
      [1.0, 1.0, 1.0],
      [-1.0, -1.0, 1.0],
      [-1.0, 1.0, -1.0],
      [1.0, -1.0, -1.0],
  ])

  # 2. Generatore di torsione locale (simulato dai momenti angolari dei vertici)
  # Un accoppiamento tensoriale che rompe la degenerazione cubica
  H_e = np.zeros((3, 3))
  for i, r in enumerate(t1):
    # Matrice di spin/torsione associata al vertice r_i
    M_i = np.outer(r, r) + 0.2 * np.array([
        [0.0, -r[2], r[1]],
        [r[2], 0.0, -r[0]],
        [-r[1], r[0], 0.0],
    ])
    H_e += M_i

  H_e = 0.5 * (H_e + H_e.T)

  # 3. Operatore di massa per i neutrini H_nu accoppiato alla torsione globale K
  # Derivato direttamente dall'inverso dei minori con la componente di scambio
  K_sub = np.array([[0.0, 0.4, -0.2], [-0.4, 0.0, 0.5], [0.2, -0.5, 0.0]])
  H_nu = LA.inv(H_e @ H_e.T + 1e-3 * np.eye(3)) + 0.5 * (
      K_sub @ K_sub.T + K_sub
  )
  H_nu = 0.5 * (H_nu + H_nu.T)

  # 4. Diagonalizzazione e calcolo PMNS
  evals_e, U_e = LA.eigh(H_e)
  evals_nu, U_nu = LA.eigh(H_nu)

  U_pmns = U_e.conjugate().T @ U_nu

  s13_sq = np.abs(U_pmns[0, 2]) ** 2
  s12_sq = (
      np.abs(U_pmns[0, 1]) ** 2 / (1.0 - s13_sq)
      if s13_sq < 1.0
      else 0.0
  )
  s23_sq = (
      np.abs(U_pmns[1, 2]) ** 2 / (1.0 - s13_sq)
      if s13_sq < 1.0
      else 0.0
  )

  j_cp = np.imag(
      U_pmns[0, 0] * U_pmns[1, 1] * np.conj(U_pmns[0, 1]) * np.conj(U_pmns[1, 0])
  )

  print("  Matrice U_PMNS con rottura anisotropa (valori assoluti):")
  print(np.round(np.abs(U_pmns), 4))

  print("\n  Risultati geometrici con torsione:")
  print(f"    sin^2(theta_12) = {s12_sq:.4f} (NuFIT: ~0.308)")
  print(f"    sin^2(theta_13) = {s13_sq:.4f} (NuFIT: ~0.022)")
  print(f"    sin^2(theta_23) = {s23_sq:.4f} (NuFIT: ~0.470)")
  print(f"    Invariante di Jarlskog (J_CP) = {j_cp:.5f}")


if __name__ == "__main__":
  anisotropic_breaking_and_mixing()

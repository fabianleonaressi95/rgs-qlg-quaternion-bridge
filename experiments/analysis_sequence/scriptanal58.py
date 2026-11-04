import numpy as np
import scipy.linalg as LA


def hierarchical_torsion_and_cp_violation():
  print("=== MODELLO GERARCHICO CON TORSIONE RIDOTTA E FASE CP ===")

  # 1. Vertici della Stella Octangula
  t1 = np.array([
      [1.0, 1.0, 1.0],
      [-1.0, -1.0, 1.0],
      [-1.0, 1.0, -1.0],
      [1.0, -1.0, -1.0],
  ])

  # 2. Matrice di massa dei leptoni carichi H_e con deformazione gerarchica debole
  H_e_real = np.zeros((3, 3))
  for r in t1:
    M_i = np.outer(r, r) + 0.02 * np.array([
        [0.0, -r[2], r[1]],
        [r[2], 0.0, -r[0]],
        [-r[1], r[0], 0.0],
    ])
    H_e_real += M_i
  H_e_real = 0.5 * (H_e_real + H_e_real.T)

  # Introduciamo la componente complessa ereditata dalla torsione K
  epsilon = 0.02  # Scala gerarchica ridotta per theta_13
  K_sub = np.array([[0.0, 1.0, -0.5], [-1.0, 0.0, 0.8], [0.5, -0.8, 0.0]])

  # H_e con piccola fase torsionale complessa hermitiana
  H_e = H_e_real + 1j * epsilon * K_sub
  H_e = 0.5 * (H_e + H_e.conj().T)  # Hermitiana complessa

  # 3. Operatore di massa per i neutrini H_nu con accoppiamento al Kernel e torsione
  H_nu_real = LA.inv(H_e_real @ H_e_real.T + 1e-3 * np.eye(3))
  H_nu = H_nu_real + 1j * (2.0 * epsilon) * (
      K_sub @ K_sub.T + 0.1 * K_sub
  )
  H_nu = 0.5 * (H_nu + H_nu.conj().T)

  # 4. Diagonalizzazione unitaria complessa (Hermitiana)
  evals_e, U_e = LA.eigh(H_e)
  evals_nu, U_nu = LA.eigh(H_nu)

  # 5. Matrice PMNS complessa
  U_pmns = U_e.conj().T @ U_nu

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

  # Invariante di Jarlskog J_CP per matrici complesse
  j_cp = np.imag(
      U_pmns[0, 0] * U_pmns[1, 1] * np.conj(U_pmns[0, 1]) * np.conj(U_pmns[1, 0])
  )

  print("  Matrice U_PMNS gerarchica complessa (valori assoluti):")
  print(np.round(np.abs(U_pmns), 4))

  print("\n  Risultati geometrici gerarchici:")
  print(f"    sin^2(theta_12) = {s12_sq:.4f} (NuFIT: ~0.308)")
  print(f"    sin^2(theta_13) = {s13_sq:.4f} (NuFIT: ~0.022)")
  print(f"    sin^2(theta_23) = {s23_sq:.4f} (NuFIT: ~0.470)")
  print(f"    Invariante di Jarlskog (J_CP) = {j_cp:.5f}")


if __name__ == "__main__":
  hierarchical_torsion_and_cp_violation()

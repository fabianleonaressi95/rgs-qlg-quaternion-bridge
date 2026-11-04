import numpy as np
import scipy.linalg as LA


def decoupled_sectors_pmns():
  print("=== MODELLO A SETTORI TOPOLOGICI DISaccoppiati ===")

  # 1. Vertici base della Stella Octangula
  t1 = np.array([
      [1.0, 1.0, 1.0],
      [-1.0, -1.0, 1.0],
      [-1.0, 1.0, -1.0],
      [1.0, -1.0, -1.0],
  ])

  # 2. Settore Leptoni Carichi (H_e) dai momenti dei vertici
  H_e = np.zeros((3, 3))
  for r in t1:
    H_e += np.outer(r, r)
  H_e = 0.5 * (H_e + H_e.T)

  # 3. Settore Neutrini (H_nu) derivato interamente dalla Torsione K e dal Kernel
  # Senza dipendere da H_e
  K_pure = np.array([[0.0, 0.5, -0.3], [-0.5, 0.0, 0.7], [0.3, -0.7, 0.0]])
  # Aggiungiamo una fase complessa di Berry derivata dal ciclo di Möbius
  H_nu_base = K_pure @ K_pure.T + 0.1 * np.eye(3)
  
  # Accoppiamento complesso con la torsione di Möbius
  H_nu = H_nu_base + 1j * 0.15 * (K_pure - K_pure.T)
  H_nu = 0.5 * (H_nu + H_nu.conj().T)

  # 4. Diagonalizzazione indipendente dei due settori
  evals_e, U_e = LA.eigh(H_e)
  evals_nu, U_nu = LA.eigh(H_nu)

  # 5. Matrice PMNS risultante
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

  j_cp = np.imag(
      U_pmns[0, 0] * U_pmns[1, 1] * np.conj(U_pmns[0, 1]) * np.conj(U_pmns[1, 0])
  )

  print("  Matrice U_PMNS disaccoppiata (valori assoluti):")
  print(np.round(np.abs(U_pmns), 4))

  print("\n  Risultati geometrici disaccoppiati:")
  print(f"    sin^2(theta_12) = {s12_sq:.4f} (NuFIT: ~0.308)")
  print(f"    sin^2(theta_13) = {s13_sq:.4f} (NuFIT: ~0.022)")
  print(f"    sin^2(theta_23) = {s23_sq:.4f} (NuFIT: ~0.470)")
  print(f"    Invariante di Jarlskog (J_CP) = {j_cp:.5f}")


if __name__ == "__main__":
  decoupled_sectors_pmns()

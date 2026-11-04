import numpy as np
import scipy.linalg as LA


def unified_breaking_and_rg_analysis():
  print("=== ANALISI UNIFICATA: ROTTURA DI SIMMETRIA E RG FLOW ===")

  # 1. Vertici iniziali della Stella Octangula (Simmetria pura)
  t1 = np.array([
      [1.0, 1.0, 1.0],
      [-1.0, -1.0, 1.0],
      [-1.0, 1.0, -1.0],
      [1.0, -1.0, -1.0],
  ])

  # 2. Rottura di Simmetria non-lineare (Ginzburg-Landau feedback sui vertici)
  # Parametro di accoppiamento geometrico beta e raggio R_0^2 = 3.0
  beta = 0.15
  R0_sq = 3.0

  t1_deformed = np.copy(t1)
  for step in range(50):  # Iterazione di rilassamento del campo
    norm_sq = np.sum(t1_deformed**2, axis=1, keepdims=True)
    force = -t1_deformed + beta * (norm_sq - R0_sq) * t1_deformed
    t1_deformed -= 0.05 * force

  print("  Vertici deformati dalla non-linearità (primi 2 vertici):")
  print(np.round(t1_deformed[:2], 4))

  # 3. Costruzione di H_e con i vertici rotti di simmetria
  H_e = np.zeros((3, 3))
  for a in range(3):
    for b in range(3):
      H_e[a, b] = np.sum(t1_deformed[:, a] * t1_deformed[:, b])
  H_e = 0.5 * (H_e + H_e.T)

  # 4. Simulazione del flusso RG (Evoluzione di scala da UV a EW)
  # Introduciamo la torsione K come generatore del flusso logaritmico
  K_sub = np.array([[0.0, 0.2, -0.1], [-0.2, 0.0, 0.3], [0.1, -0.3, 0.0]])

  H_nu_uv = LA.inv(H_e @ H_e.T + 0.1 * np.eye(3))
  # Evoluzione RG di un "passo" logaritmico t ~ ln(M_UV / M_EW)
  t_flow = 1.2
  H_nu_ew = H_nu_uv + t_flow * (K_sub @ H_nu_uv + H_nu_uv @ K_sub.T)
  H_nu_ew = 0.5 * (H_nu_ew + H_nu_ew.T)

  # 5. Diagonalizzazione e calcolo PMNS finale
  evals_e, U_e = LA.eigh(H_e)
  evals_nu, U_nu = LA.eigh(H_nu_ew)

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

  print("\n  Risultati dopo Rottura di Simmetria e RG Flow:")
  print(f"    sin^2(theta_12) = {s12_sq:.4f} (Target NuFIT: ~0.308)")
  print(f"    sin^2(theta_13) = {s13_sq:.4f} (Target NuFIT: ~0.022)")
  print(f"    sin^2(theta_23) = {s23_sq:.4f} (Target NuFIT: ~0.470)")
  print(f"    Invariante di Jarlskog (J_CP) = {j_cp:.5f}")


if __name__ == "__main__":
    unified_breaking_and_rg_analysis()

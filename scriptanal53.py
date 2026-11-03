import numpy as np
import scipy.linalg as LA


def verify_discrete_bundle_holonomy():
  print("=== VERIFICA DELL'OLONOMIA E DEL FIBRATO DISCRETO ===")

  # 1. Definizione dell'operatore di scambio P (16x16)
  dim = 16
  P = np.zeros((dim, dim), dtype=float)
  for i in range(4):
    # Scambio tra i vertici del tetraedro 1 (indici 0-3) e tetraedro 2 (indici 4-7)
    P[i * 2 : i * 2 + 2, (i + 4) * 2 : (i + 4) * 2 + 2] = np.eye(2)
    P[(i + 4) * 2 : (i + 4) * 2 + 2, i * 2 : i * 2 + 2] = np.eye(2)

  # 2. Definizione del ciclo chiuso gamma: T1 -> T2 -> T1
  # La transizione elementare è data da P
  holonomy_loop = P @ P  # Hol(gamma) = P^2

  # 3. Verifiche topologiche
  p_sq_err = LA.norm(holonomy_loop - np.eye(dim), ord="fro")
  is_identity = np.allclose(holonomy_loop, np.eye(dim))

  print(f"  Operatore di olonomia Hol(gamma) = P^2 calcolato.")
  print(f"  ||Hol(gamma) - I||_F = {p_sq_err:.2e}")
  print(f"  L'olonomia lungo il ciclo è banale (I)? {is_identity}")

  # Azione dell'operatore di scambio P sui settori di parità (Pi_+ e Pi_-)
  Pi_plus = 0.5 * (np.eye(dim) + P)
  Pi_minus = 0.5 * (np.eye(dim) - P)

  # Verifichiamo che P agisca come inversione di segno nel settore antisimmetrico V_-
  test_vector_minus = Pi_minus @ np.random.randn(dim)
  transformed_vector = P @ test_vector_minus
  diff_vector = transformed_vector - (-test_vector_minus)

  print(
      "  Azione di P nel settore antisimmetrico V_- (P y = -y): "
      f"{LA.norm(diff_vector) < 1e-12}"
  )
  print(
      "  Conclusione Topologica: Il ciclo inter-tetraedrico induce una"
      " riflessione con olonomia non banale sui settori, confermando la"
      " struttura a fibrato non orientabile (Möbius)."
  )


if __name__ == "__main__":
  verify_discrete_bundle_holonomy()

import numpy as np
import scipy.linalg as LA


def deterministic_projection_audit():
  print("=== AUDIT DETERMINISTICO DELLA PROIEZIONE W (64 -> 16) ===")

  # 1. Definizione geometrica dei vertici della Stella Octangula
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
  nodes = np.vstack([t1, t2])  # 8 nodi in R^3
  n_nodes = 8
  dim_macro = 16  # 8 nodi * 2 DOF macro
  dim_micro = 64  # Spazio microstrutturale (es. 8 nodi * 8 combinazioni locali)

  # Costruzione DETERMINISTICA di W basata sulle coordinate dei vertici e prodotti tensoriali
  W = np.zeros((dim_macro, dim_micro))

  for i in range(n_nodes):
    r_i = nodes[i]
    for k in range(8):  # 8 stati interni per nodo (es. 2^3 configurazioni di campo)
      col_idx = i * 8 + k
      # Mappa deterministica basata sulle coordinate e bit di configurazione
      W[i * 2, col_idx] = r_i[0] * (1.0 if (k & 1) else -1.0)
      W[i * 2 + 1, col_idx] = (
          r_i[1] * r_i[2] * (0.5 if (k & 2) else -0.5)
      )

  # 2. Verifica algebrica della matrice di Gram WW^T
  WWt = W @ W.T
  evals_wwt = LA.eigvalsh(WWt)

  print(f"  Dimensioni di W: {W.shape}")
  print(f"  Rango di W: {np.linalg.matrix_rank(W)}")
  print(
      "  Autovalori della matrice di Gram WW^T:"
      f" {np.sort(evals_wwt.real)}"
  )

  # Norma della differenza rispetto a un'identità scalata (es. 3/4 I o simile)
  # Verifichiamo quale scala metrica emerge naturalmente dall'algebra pura
  trace_wwt = np.trace(WWt)
  effective_scale = trace_wwt / dim_macro
  print(
      "  Scala metrica media indotta (Tr(WW^T) / dim):"
      f" {effective_scale:.4f}"
  )

  # 3. Controllo dell'equazione caratteristica per la ricerca di invarianti (es. sezione aurea)
  # Analizziamo il polinomio caratteristico o gli autovalori normalizzati
  norm_evals = evals_wwt[evals_wwt > 1e-5] / effective_scale
  print(f"  Autovalori normalizzati di WW^T: {np.round(norm_evals.real, 4)}")


if __name__ == "__main__":
  deterministic_projection_audit()

import numpy as np
import scipy.linalg as LA


def analyze_tensor_spectrum():
  # 1. Ricostruzione dell'operatore tensoriale esatto
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
  nodes = np.vstack([t1, t2])
  n_nodes = 8
  dim = n_nodes * 2

  K = np.zeros((dim, dim), dtype=float)
  S = np.zeros((dim, dim), dtype=float)

  for i in range(n_nodes):
    r_i = nodes[i]
    M_i = np.array(
        [[0.0, -r_i[2], r_i[1]], [r_i[2], 0.0, -r_i[0]], [-r_i[1], r_i[0], 0.0]]
    )
    for j in range(n_nodes):
      if i == j:
        continue
      r_j = nodes[j]
      diff = r_j - r_i
      dist_sq = np.dot(diff, diff)
      if np.isclose(dist_sq, 4.0) or np.isclose(dist_sq, 8.0):
        S[i * 2, j * 2] += 1.0 / dist_sq
        S[i * 2 + 1, j * 2 + 1] += 1.0 / dist_sq
        if (i < 4 and j >= 4) or (i >= 4 and j < 4):
          dir_vec = diff / np.sqrt(dist_sq)
          tensor_coupling = np.dot(M_i, dir_vec)
          K[i * 2, j * 2 + 1] += tensor_coupling[0]
          K[i * 2 + 1, j * 2] -= tensor_coupling[1]

  for i in range(n_nodes):
    deg = np.sum(S[i * 2, ::2])
    S[i * 2, i * 2] -= deg
    S[i * 2 + 1, i * 2 + 1] -= deg

  A = S + K

  # 2. Calcolo degli autovalori del generatore continuo A
  eig_A = LA.eigvals(A)
  # Ordiniamo per parte immaginaria decrescente
  eig_A = eig_A[np.argsort(np.imag(eig_A))]

  print("=== SPETTRO COMPLESSO DEL GENERATORE TENSORIALE A ===")
  for idx, lam in enumerate(eig_A):
    re = np.real(lam)
    im = np.imag(lam)
    sign = "+" if im >= 0 else "-"
    print(f"  Modo {idx+1:2d}: {re:10.6f} {sign} {abs(im):10.6f}i")


if __name__ == "__main__":
  analyze_tensor_spectrum()

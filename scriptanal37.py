
import numpy as np

def audit_operator(A, P=None, tol=1e-10):
    A = np.asarray(A, dtype=float)
    assert A.ndim == 2 and A.shape[0] == A.shape[1]

    I = np.eye(A.shape[0])
    S = (A + A.T) / 2
    K = (A - A.T) / 2

    comm_normal = A @ A.T - A.T @ A
    eigvals, V = np.linalg.eig(A)

    print("Dimensione:", A.shape)
    print("||A - S - K||_F:",
          np.linalg.norm(A - S - K, ord="fro"))
    print("||K + K^T||_F:",
          np.linalg.norm(K + K.T, ord="fro"))
    print("||[A,A^T]||_F:",
          np.linalg.norm(comm_normal, ord="fro"))
    print("Non-normalità certificata:",
          np.linalg.norm(comm_normal, ord="fro") > tol)
    print("Autovalori:")
    print(np.sort_complex(eigvals))

    if P is not None:
        P = np.asarray(P, dtype=float)
        print("||P^2-I||_F:",
              np.linalg.norm(P @ P - I, ord="fro"))
        print("||[A,P]||_F:",
              np.linalg.norm(A @ P - P @ A, ord="fro"))

    return {
        "S": S,
        "K": K,
        "commutator": comm_normal,
        "eigenvalues": eigvals,
        "eigenvectors": V,
    }

"""Reproducible audit for the trigram Q3 theorem dossier."""
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.trigram_q3_theorems import (  # noqa: E402
    q3_adjacency, q6_product_adjacency, refinement_reduction,
    tribonacci_target, walsh_matrix,
)


def main():
    E, A3, M = walsh_matrix(), q3_adjacency(), tribonacci_target()
    # Diagnostic realization only; not a geometric Kathara/Krist-All refinement.
    R_diagnostic = E @ M @ E.T
    M_recovered, leakage = refinement_reduction(R_diagnostic)
    A6 = q6_product_adjacency()
    print("Walsh orthogonality residual:", np.linalg.norm(E.T @ E - np.eye(3)))
    print("Q3 adjacency eigenvalues:", np.linalg.eigvalsh(A3).round(12).tolist())
    print("Q3 vertices / edges:", A3.shape[0], int(A3.sum() // 2))
    print("Q6 vertices / edges / degrees:",
          A6.shape[0], int(A6.sum() // 2), np.unique(A6.sum(axis=1)).tolist())
    print("Diagnostic compressed matrix:\n", M_recovered)
    print("Diagnostic channel leakage:", leakage)
    print("Target characteristic polynomial coefficients:", np.poly(M).tolist())
    print("Target eigenvalues:", np.linalg.eigvals(M.astype(float)).tolist())
    print("NOTE: diagnostic R is not derived from geometric incidence data.")


if __name__ == "__main__":
    main()

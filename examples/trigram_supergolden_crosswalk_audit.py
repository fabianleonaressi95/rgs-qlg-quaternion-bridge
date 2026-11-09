from __future__ import annotations
import json
from pathlib import Path
import numpy as np


def main():
    A = np.array([[1, 0, 1], [1, 0, 0], [0, 1, 0]], dtype=int)
    M = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 1]], dtype=int)
    P = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], dtype=int)
    report = {
        "permutation_matrix": P.tolist(),
        "P_orthogonality_residual": int(np.max(np.abs(P.T @ P - np.eye(3, dtype=int)))),
        "similarity_residual": int(np.max(np.abs(P @ A @ P.T - M))),
        "phase21_characteristic_polynomial": np.poly(A).tolist(),
        "trigram_characteristic_polynomial": np.poly(M).tolist(),
        "interpretation": "Exact permutation similarity only; no geometric cell intertwiner is established.",
    }
    out = Path("results/trigram_supergolden_crosswalk/audit.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

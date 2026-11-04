import numpy as np
import scipy.sparse as sp
import scipy.linalg as LA
from scipy.spatial import Delaunay

# The original phase-diagram driver was committed without its three numerical
# helpers.  Reuse the independently committed implementation from peierls46;
# keep this file's historical driver and output protocol unchanged.
try:
    from .peierls46 import (
        generate_amorphous_lattice,
        build_hamiltonian,
        find_local_optimal_gap,
    )
except ImportError:  # Support: python experiments/peierls/peierls45.py
    from peierls46 import (
        generate_amorphous_lattice,
        build_hamiltonian,
        find_local_optimal_gap,
    )


def scan_magnetic_phase_diagram(N=1024, num_realizations=5, b_list=[0.5, 1.0, 1.5, 2.0]):
    print(f"=== Scansione del Diagramma di Fase Topologico (N={N}) ===\n")
    print(f"{'Campo B':<8} | {'Target Search':<14} | {'Mediana C (Ensemble)':<20} | {'Std Dev Inter-Ensemble':<22}")
    print("-" * 72)

    for B in b_list:
        ensemble_medians = []

        for r in range(num_realizations):
            seed_val = 200 + r
            coords, box_size = generate_amorphous_lattice(N, density=1.0, seed=seed_val)
            H = build_hamiltonian(coords, box_size, B_field=B, gamma=0.3)

            evals, eigenvectors = LA.eigh(H.toarray())
            num_occupied = find_local_optimal_gap(evals, target_fraction=0.78, window=0.12)

            V_occ = eigenvectors[:, :num_occupied]
            P = V_occ @ V_occ.conj().T
            Q = np.eye(N) - P

            X = np.diag(coords[:, 0])
            Y = np.diag(coords[:, 1])

            P_X_Q_Y_P = P @ X @ Q @ Y @ P
            local_chern = -4.0 * np.pi * np.imag(np.diagonal(P_X_Q_Y_P))

            center = coords.mean(axis=0)
            radii = np.linalg.norm(coords - center, axis=1)
            bulk_mask = radii <= (0.32 * radii.max())

            med_val = np.median(local_chern[bulk_mask])
            ensemble_medians.append(med_val)

        mean_c = np.mean(ensemble_medians)
        std_c = np.std(ensemble_medians)
        print(f"{B:<8.2f} | Adattivo       | {mean_c:20.4f} | {std_c:22.4f}")


if __name__ == "__main__":
    scan_magnetic_phase_diagram()

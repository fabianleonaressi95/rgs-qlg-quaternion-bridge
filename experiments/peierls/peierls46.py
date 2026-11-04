import numpy as np
import scipy.sparse as sp
import scipy.linalg as LA
from scipy.spatial import Delaunay

def generate_amorphous_lattice(N, density=1.0, steps=15, dt=0.01, seed=None):
    """Genera un reticolo amorfo con densità costante e regolarizzazione anti-collisione."""
    box_size = np.sqrt(N / density)
    rng = np.random.default_rng(seed)
    coords = rng.random((N, 2)) * box_size
    for _ in range(steps):
        diff = coords[:, None, :] - coords[None, :, :]
        dist = np.linalg.norm(diff, axis=-1)
        dist = np.maximum(dist, 1e-3)
        np.fill_diagonal(dist, 1.0)
        
        inv_dist3 = 1.0 / (dist**3)
        np.fill_diagonal(inv_dist3, 0.0)
        
        forces = np.sum(diff * inv_dist3[:, :, None], axis=1)
        coords += dt * forces
        coords = np.clip(coords, 0, box_size)
    return coords, box_size

def build_hamiltonian(coords, box_size, B_field=1.0, gamma=0.3, t1=1.0):
    """Costruisce l'Hamiltoniana con decadimento spaziale e gauge di Peierls."""
    N = coords.shape[0]
    tri = Delaunay(coords)
    rows, cols, data = [], [], []
    nn_edges = set()
    for simplex in tri.simplices:
        for i_local in range(3):
            u = simplex[i_local]
            v = simplex[(i_local + 1) % 3]
            nn_edges.add((min(u, v), max(u, v)))
            
    for u, v in nn_edges:
        pos_u = coords[u]
        pos_v = coords[v]
        dist = np.linalg.norm(pos_v - pos_u)
        if dist > 1e-4:
            mx = 0.5 * (pos_u[0] + pos_v[0])
            my = 0.5 * (pos_u[1] + pos_v[1])
            Ax = -0.5 * B_field * my
            Ay = 0.5 * B_field * mx
            dx = pos_v[0] - pos_u[0]
            dy = pos_v[1] - pos_u[1]
            peierls_phase = Ax * dx + Ay * dy
            t_hop = t1 * np.exp(-gamma * dist) * np.exp(1j * peierls_phase)
            rows.extend([u, v])
            cols.extend([v, u])
            data.extend([t_hop, np.conj(t_hop)])
            
    return sp.coo_matrix((data, (rows, cols)), shape=(N, N)).tocsr()

def find_local_optimal_gap(evals, target_fraction=0.78, window=0.12):
    """Individua il gap spettrale ottimale (massimo spaziamiento dei livelli) per ciascun seed."""
    N = len(evals)
    center_idx = int(N * target_fraction)
    wing = int(N * window)
    
    start_idx = max(10, center_idx - wing)
    end_idx = min(N - 10, center_idx + wing)
    
    spacings = np.diff(evals[start_idx:end_idx])
    if len(spacings) == 0:
        return center_idx
        
    optimal_local_idx = np.argmax(spacings)
    return start_idx + optimal_local_idx + 1

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

import numpy as np
import scipy.sparse as sp
import scipy.linalg as LA
from scipy.spatial import Delaunay

def generate_amorphous_lattice(N, box_size=14.0, steps=10, dt=0.01, seed=None):
    rng = np.random.default_rng(seed)
    coords = rng.random((N, 2)) * box_size
    for _ in range(steps):
        diff = coords[:, None, :] - coords[None, :, :]
        dist = np.linalg.norm(diff, axis=-1)
        np.fill_diagonal(dist, 1e-5)
        inv_dist3 = 1.0 / (dist**3)
        forces = np.sum(diff * inv_dist3[:, :, None], axis=1)
        coords += dt * forces
        coords = np.clip(coords, 0, box_size)
    return coords

def build_decay_peierls_hamiltonian(coordinates, t1=1.0, B_field=1.0, gamma=0.3):
    N = coordinates.shape[0]
    tri = Delaunay(coordinates)
    
    rows, cols, data = [], [], []
    nn_edges = set()
    for simplex in tri.simplices:
        for i_local in range(3):
            u = simplex[i_local]
            v = simplex[(i_local + 1) % 3]
            nn_edges.add((min(u, v), max(u, v)))
            
    for u, v in nn_edges:
        pos_u = coordinates[u]
        pos_v = coordinates[v]
        dist = np.linalg.norm(pos_v - pos_u)
        
        if dist > 1e-4:
            mx = 0.5 * (pos_u[0] + pos_v[0])
            my = 0.5 * (pos_u[1] + pos_v[1])
            
            Ax = -0.5 * B_field * my
            Ay = 0.5 * B_field * mx
            
            dx = pos_v[0] - pos_u[0]
            dy = pos_v[1] - pos_u[1]
            peierls_phase = Ax * dx + Ay * dy
            
            # Hopping con decadimento spaziale controllato da gamma
            t_hop = t1 * np.exp(-gamma * dist) * np.exp(1j * peierls_phase)
            
            rows.extend([u, v])
            cols.extend([v, u])
            data.extend([t_hop, np.conj(t_hop)])

    H = sp.coo_matrix((data, (rows, cols)), shape=(N, N)).tocsr()
    return H

def compute_bulk_chern(H, coordinates, filling_fraction=0.8):
    N = H.shape[0]
    H_dense = H.toarray() if hasattr(H, "toarray") else H
    evals, eigenvectors = LA.eigh(H_dense)
    
    num_occupied = int(N * filling_fraction)
    V_occ = eigenvectors[:, :num_occupied]
    P = V_occ @ V_occ.conj().T
    Q = np.eye(N) - P
    
    X = np.diag(coordinates[:, 0])
    Y = np.diag(coordinates[:, 1])
    
    P_X_Q_Y_P = P @ X @ Q @ Y @ P
    local_chern = -4.0 * np.pi * np.imag(np.diagonal(P_X_Q_Y_P))
    
    center = coordinates.mean(axis=0)
    distances = np.linalg.norm(coordinates - center, axis=1)
    bulk_mask = distances < (0.35 * distances.max())
    
    mean_bulk_chern = np.mean(local_chern[bulk_mask]) if np.any(bulk_mask) else 0.0
    return mean_bulk_chern

def finite_size_scaling_analysis(n_list, b_val=1.0, gamma_val=0.3, seeds=[1337, 42, 7, 99, 2024]):
    print(f"\n--- Analisi di Scaling in N (B={b_val}, gamma={gamma_val}) ---")
    print(f"{'N':<6} | {'Chern Medio':<12} | {'Dev. Std':<10} | {'Extrapolated (1/N -> 0)'}")
    print("-" * 55)
    
    mean_chern_per_n = []
    inv_n_list = []
    
    for N in n_list:
        chern_runs = []
        for seed in seeds:
            coords = generate_amorphous_lattice(N, box_size=np.sqrt(N) * 1.0, seed=seed)
            H = build_decay_peierls_hamiltonian(coords, B_field=b_val, gamma=gamma_val)
            c_val = compute_bulk_chern(H, coords, filling_fraction=0.8)
            chern_runs.append(c_val)
            
        m_c = np.mean(chern_runs)
        s_c = np.std(chern_runs)
        mean_chern_per_n.append(m_c)
        inv_n_list.append(1.0 / N)
        print(f"{N:<6} | {m_c:<12.4f} | {s_c:<10.4f}")
        
    # Fit lineare per estrapolazione a N -> infinito (1/N = 0)
    coeffs = np.polyfit(inv_n_list, mean_chern_per_n, 1)
    intercept = coeffs[0] * 0 + coeffs[1]
    print(f"\nEstrapolazione al limite termodinamico (N -> inf): C = {intercept:.4f}")
    return intercept

if __name__ == "__main__":
    print("=== Avvio Estrapolazione Avanzata Invariante Topologico ===")
    
    # Test di scaling con i parametri ottimali
    system_sizes = [100, 196, 324, 400]
    final_limit = finite_size_scaling_analysis(system_sizes, b_val=1.0, gamma_val=0.3)

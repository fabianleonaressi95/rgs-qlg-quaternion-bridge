import numpy as np
import scipy.sparse as sp
import scipy.linalg as LA
from scipy.spatial import Delaunay

def generate_amorphous_lattice(N, density=1.0, steps=15, dt=0.01, seed=None):
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

def compute_chern_ensemble(N, B_field, realizations=4):
    """
    Esegue un ensemble di simulazioni per una data taglia N e campo B 
    per calcolare la stabilità statistica del marcatore di Chern.
    """
    chern_values = []
    for r in range(realizations):
        coords, _ = generate_amorphous_lattice(N=N, density=1.0, steps=15, seed=200 + r * 37)
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
                
                t_hop = np.exp(-0.3 * dist) * np.exp(1j * peierls_phase)
                rows.extend([u, v])
                cols.extend([v, u])
                data.extend([t_hop, np.conj(t_hop)])
                
        H = sp.coo_matrix((data, (rows, cols)), shape=(N, N)).tocsr()
        evals, eigenvectors = LA.eigh(H.toarray())
        
        # Ricerca adattiva del gap di bulk (filling ottimale ~83%)
        center_idx = int(N * 0.83)
        wing = int(N * 0.10)
        start_idx = max(10, center_idx - wing)
        end_idx = min(N - 10, center_idx + wing)
        
        spacings = np.diff(evals[start_idx:end_idx])
        num_occupied = start_idx + np.argmax(spacings) + 1 if len(spacings) > 0 else center_idx
        
        # Marcatore locale di Bianco-Resta
        V_occ = eigenvectors[:, :num_occupied]
        P = V_occ @ V_occ.conj().T
        Q = np.eye(N) - P
        
        X = np.diag(coords[:, 0])
        Y = np.diag(coords[:, 1])
        local_chern = -4.0 * np.pi * np.imag(np.diagonal(P @ X @ Q @ Y @ P))
        
        center = coords.mean(axis=0)
        radii = np.linalg.norm(coords - center, axis=1)
        bulk_vals = local_chern[radii <= (0.32 * radii.max())]
        chern_values.append(np.median(bulk_vals))
        
    return np.median(chern_values), np.std(chern_values)

if __name__ == "__main__":
    print("=== SCANSIONE FINITE-SIZE SCALING: VERIFICA PLATEAU C = -2 vs C = -1 ===")
    print(f"{'B Field':<10} | {'Nodi (N)':<10} | {'Mediana C':<15} | {'Std Dev':<15}")
    print("-" * 60)
    
    test_cases = [
        (0.5, [512, 1024, 2048]),  # Target teorico: -2
        (2.0, [512, 1024, 2048])   # Target teorico: -1
    ]
    
    for B, sizes in test_cases:
        for N in sizes:
            med_c, std_c = compute_chern_ensemble(N, B_field=B, realizations=3)
            print(f"{B:<10.2f} | {N:<10} | {med_c:<15.4f} | {std_c:<15.4f}")
        print("-" * 60)

import numpy as np
import scipy.sparse as sp
import scipy.linalg as LA
from scipy.spatial import Delaunay

def generate_amorphous_lattice(N, box_size=12.0, steps=12, dt=0.01, seed=None):
    """
    Genera un reticolo amorfo tramite rilassamento repulsivo vettorizzato.
    """
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

def build_peierls_amorphous_hamiltonian(coordinates, t1=1.0, B_field=0.5, m_onsite=0.0):
    """
    Costruisce l'Hamiltoniana tight-binding su reticolo amorfo di Delaunay
    incorporando un campo magnetico sintetico B tramite sostituzione di Peierls (Gauge Simmetrica).
    """
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
            # Punto medio del legame per la gauge simmetrica A = 0.5 * B * (-y, x)
            mx = 0.5 * (pos_u[0] + pos_v[0])
            my = 0.5 * (pos_u[1] + pos_v[1])
            
            Ax = -0.5 * B_field * my
            Ay = 0.5 * B_field * mx
            
            # Sostituzione di Peierls: fase = int A * dl
            dx = pos_v[0] - pos_u[0]
            dy = pos_v[1] - pos_u[1]
            peierls_phase = Ax * dx + Ay * dy
            
            t_hop = t1 * np.exp(1j * peierls_phase)
            
            rows.extend([u, v])
            cols.extend([v, u])
            data.extend([t_hop, np.conj(t_hop)])

    if m_onsite != 0.0:
        rng_m = np.random.default_rng(42)
        diag_energies = m_onsite * rng_m.choice([-1.0, 1.0], size=N)
    else:
        diag_energies = np.zeros(N)

    H = sp.coo_matrix((data, (rows, cols)), shape=(N, N)).tocsr()
    H = H + sp.diags(diag_energies)
    return H

def calculate_rigorous_diagnostics(H, coordinates, filling_fraction=0.4):
    """
    Calcola il marcatore locale di Chern con P X Q Y P (con Q = I - P).
    """
    N = H.shape[0]
    H_dense = H.toarray() if hasattr(H, "toarray") else H
    evals, eigenvectors = LA.eigh(H_dense)
    
    num_occupied = int(N * filling_fraction)
    if num_occupied <= 0 or num_occupied >= N:
        return 0.0, 0.0, False
        
    mean_level_spacing = np.mean(np.diff(evals))
    local_gap = evals[num_occupied] - evals[num_occupied - 1] if 0 < num_occupied < N else 0.0
    is_true_gap = local_gap > (1.5 * mean_level_spacing)
    
    V_occ = eigenvectors[:, :num_occupied]
    P = V_occ @ V_occ.conj().T
    Q = np.eye(N) - P
    
    min_coords = coordinates.min(axis=0)
    max_coords = coordinates.max(axis=0)
    scaled_coords = (coordinates - min_coords) / (max_coords - min_coords + 1e-12)
    
    X = np.diag(scaled_coords[:, 0])
    Y = np.diag(scaled_coords[:, 1])
    
    P_X_Q_Y_P = P @ X @ Q @ Y @ P
    local_chern = -4.0 * np.pi * np.imag(np.diagonal(P_X_Q_Y_P))
    
    center = coordinates.mean(axis=0)
    distances = np.linalg.norm(coordinates - center, axis=1)
    bulk_mask = distances < (0.35 * distances.max())
    
    mean_bulk_chern = np.mean(local_chern[bulk_mask]) if np.any(bulk_mask) else 0.0
    
    return mean_bulk_chern, local_gap, is_true_gap

if __name__ == "__main__":
    print("=== Scansione Chern Insulator Amorfo con Peierls (Campo Magnetico Sintetico) ==-\n")
    
    N_test = 196
    b_fields = [0.2, 0.5, 1.0, 1.5]
    fillings = [0.25, 0.33, 0.50, 0.66]
    n_realizations = 3
    base_seed = 1337
    
    print(f"{'B_Field':<10} | {'Filling':<8} | {'Chern Medio (Bulk)':<20} | {'Gap Medio':<12} | {'True Gap Rate':<15}")
    print("-" * 77)
    
    for b_val in b_fields:
        for f in fillings:
            chern_vals = []
            gaps = []
            true_gap_count = 0
            
            for r in range(n_realizations):
                run_seed = base_seed + r * 37
                coords = generate_amorphous_lattice(N_test, box_size=np.sqrt(N_test)*1.1, seed=run_seed, steps=10)
                H_mat = build_peierls_amorphous_hamiltonian(coords, t1=1.0, B_field=b_val, m_onsite=0.0)
                
                c_val, g_val, is_gap = calculate_rigorous_diagnostics(H_mat, coords, filling_fraction=f)
                
                chern_vals.append(c_val)
                gaps.append(g_val)
                if is_gap:
                    true_gap_count += 1
                    
            mean_c = np.mean(chern_vals)
            mean_g = np.mean(gaps)
            gap_rate = (true_gap_count / n_realizations) * 100.0
            
            print(f"{b_val:<10.2f} | {f:<8.2f} | {mean_c:<+20.4f} | {mean_g:<12.4f} | {gap_rate:<14.1f}%")

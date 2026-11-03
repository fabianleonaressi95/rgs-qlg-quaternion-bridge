import numpy as np
import scipy.sparse as sp
import scipy.linalg as LA

def generate_standard_amorphous_lattice(N, box_size=10.0, steps=12, dt=0.01, seed=None):
    """
    Genera il reticolo amorfo tramite rilassamento repulsivo vettorizzato.
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

def build_tight_binding_hamiltonian(coordinates, t0=1.0, gamma=2.0, alpha=0.15, cutoff=1.8):
    """
    Costruisce l'Hamiltoniana tight-binding con fasi di Peierls.
    """
    N = coordinates.shape[0]
    diff = coordinates[:, None, :] - coordinates[None, :, :]
    dist = np.linalg.norm(diff, axis=-1)
    
    mask = (dist <= cutoff) & (dist > 1e-4) & (np.arange(N)[:, None] < np.arange(N)[None, :])
    rows, cols = np.nonzero(mask)
    
    if len(rows) == 0:
        return sp.coo_matrix((N, N)).tocsr()
        
    d = dist[rows, cols]
    t_ij = t0 * np.exp(-gamma * (d - 1.0))
    
    pos_i = coordinates[rows]
    pos_j = coordinates[cols]
    mid_y = (pos_i[:, 1] + pos_j[:, 1]) / 2.0
    delta_x = pos_j[:, 0] - pos_i[:, 0]
    phi_ij = 2.0 * np.pi * alpha * mid_y * delta_x
    
    t_complex = t_ij * np.exp(1j * phi_ij)
    
    r_all = np.concatenate([rows, cols])
    c_all = np.concatenate([cols, rows])
    d_all = np.concatenate([t_complex, np.conj(t_complex)])
    
    return sp.coo_matrix((d_all, (r_all, c_all)), shape=(N, N)).tocsr()

def calculate_rigorous_diagnostics(H, coordinates, filling_fraction=0.7):
    """
    Calcola il marcatore locale di Chern e verifica la stabilità del gap.
    """
    N = H.shape[0]
    H_dense = H.toarray() if hasattr(H, "toarray") else H
    evals, eigenvectors = LA.eigh(H_dense)
    
    num_occupied = int(N * filling_fraction)
    if num_occupied <= 0 or num_occupied >= N:
        return 0.0, 0.0, False
        
    mean_level_spacing = np.mean(np.diff(evals))
    local_gap = evals[num_occupied] - eval[num_occupied - 1] if 0 < num_occupied < N else 0.0
    # Criterio di gap leggermente più flessibile (1.2x) utile durante la scansione esplorativa
    is_true_gap = local_gap > (1.2 * mean_level_spacing)
    
    V_occ = eigenvectors[:, :num_occupied]
    P = V_occ @ V_occ.conj().T
    
    min_coords = coordinates.min(axis=0)
    max_coords = coordinates.max(axis=0)
    scaled_coords = (coordinates - min_coords) / (max_coords - min_coords + 1e-12)
    
    X = np.diag(scaled_coords[:, 0])
    Y = np.diag(scaled_coords[:, 1])
    
    P_X_P_Y_P = P @ X @ P @ Y @ P
    local_chern = -4.0 * np.pi * np.imag(np.diagonal(P_X_P_Y_P))
    
    center = coordinates.mean(axis=0)
    distances = np.linalg.norm(coordinates - center, axis=1)
    bulk_mask = distances < (0.35 * distances.max())
    
    mean_bulk_chern = np.mean(local_chern[bulk_mask]) if np.any(bulk_mask) else 0.0
    
    return mean_bulk_chern, local_gap, is_true_gap

if __name__ == "__main__":
    print("=== Scansione Parametrica 2D (Alpha x Filling) ==-\n")
    
    N_test = 144
    alphas = [0.10, 0.20, 0.30, 0.40]
    fillings = [0.25, 0.50, 0.75]
    n_realizations = 3
    base_seed = 1337
    
    print(f"{'Alpha':<8} | {'Filling':<8} | {'Chern Medio':<14} | {'Gap Medio':<12} | {'True Gap Rate':<15}")
    print("-" * 65)
    
    for alpha in alphas:
        for f in fillings:
            chern_vals = []
            gaps = []
            true_gap_count = 0
            
            for r in range(n_realizations):
                run_seed = base_seed + r * 31
                coords = generate_standard_amorphous_lattice(N_test, box_size=np.sqrt(N_test)*1.2, seed=run_seed, steps=12)
                H_mat = build_tight_binding_hamiltonian(coords, alpha=alpha, cutoff=1.8)
                
                c_val, g_val, is_gap = calculate_rigorous_diagnostics(H_mat, coords, filling_fraction=f)
                
                chern_vals.append(c_val)
                gaps.append(g_val)
                if is_gap:
                    true_gap_count += 1
                    
            mean_c = np.mean(chern_vals)
            mean_g = np.mean(gaps)
            gap_rate = (true_gap_count / n_realizations) * 100.0
            
            print(f"{alpha:<8.2f} | {f:<8.2f} | {mean_c:<+14.4f} | {mean_g:<12.4f} | {gap_rate:<14.1f}%")

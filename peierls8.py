import numpy as np
import scipy.sparse as sp
import scipy.linalg as LA
from scipy.spatial import cKDTree

def generate_standard_amorphous_lattice(N, box_size=10.0, steps=30, dt=0.01, seed=None):
    """
    Genera un reticolo amorfo con protocollo di rilassamento repulsivo standardizzato e riproducibile.
    """
    rng = np.random.default_rng(seed)
    coords = rng.random((N, 2)) * box_size
    
    for _ in range(steps):
        tree = cKDTree(coords)
        pairs = tree.query_pairs(r=box_size * 0.4)
        
        forces = np.zeros_like(coords)
        for i, j in pairs:
            diff = coords[i] - coords[j]
            dist = np.linalg.norm(diff)
            if dist < 1e-4:
                dist = 1e-4
            f = diff / (dist**3)
            forces[i] += f
            forces[j] -= f
            
        coords += dt * forces
        coords = np.clip(coords, 0, box_size)
    return coords

def build_tight_binding_hamiltonian(coordinates, t0=1.0, gamma=2.0, alpha=0.15, cutoff=1.8):
    """
    Costruisce l'Hamiltoniana con cutoff ristretto (1.8) per preservare la struttura aperiodica
    senza sovra-connessione spettrale.
    """
    tree = cKDTree(coordinates)
    pairs = tree.query_pairs(r=cutoff)
    
    N = coordinates.shape[0]
    rows, cols, data = [], [], []
    
    for i, j in pairs:
        if i == j:
            continue
        pos_i = coordinates[i]
        pos_j = coordinates[j]
        dist = np.linalg.norm(pos_i - pos_j)
        
        if dist > 1e-4:
            t_ij = t0 * np.exp(-gamma * (dist - 1.0))
            mid_y = (pos_i[1] + pos_j[1]) / 2.0
            delta_x = pos_j[0] - pos_i[0]
            phi_ij = 2.0 * np.pi * alpha * mid_y * delta_x
            
            t_complex = t_ij * np.exp(1j * phi_ij)
            rows.extend([i, j])
            cols.extend([j, i])
            data.extend([t_complex, np.conj(t_complex)])
            
    return sp.coo_matrix((data, (rows, cols)), shape=(N, N)).tocsr()

def calculate_rigorous_diagnostics(H, coordinates, filling_fraction=0.7):
    """
    Calcola il marcatore locale di Chern e valuta rigorosamente la presenza di un gap reale
    confrontandolo con la spaziatura media dei livelli.
    """
    N = H.shape[0]
    H_dense = H.toarray() if hasattr(H, "toarray") else H
    evals, eigenvectors = LA.eigh(H_dense)
    
    num_occupied = int(N * filling_fraction)
    if num_occupied <= 0 or num_occupied >= N:
        return 0.0, 0.0, 0.0, False
        
    # Spaziatura media globale dei livelli nello spettro
    mean_level_spacing = np.mean(np.diff(evals))
    
    # Gap energetico locale in corrispondenza del livello di Fermi (indice num_occupied - 1 e num_occupied)
    if 0 < num_occupied < N:
        local_gap = evals[num_occupied] - evals[num_occupied - 1]
    else:
        local_gap = 0.0
        
    # Criterio di validità del gap: deve essere almeno 1.5x la spaziatura media per non essere rumore
    is_true_gap = local_gap > (1.5 * mean_level_spacing)
    
    # Calcolo Marcatore Locale di Chern (Bianco-Resta)
    V_occ = eigenvectors[:, :num_occupied]
    P = V_occ @ V_occ.conj().T
    
    min_coords = coordinates.min(axis=0)
    max_coords = coordinates.max(axis=0)
    scaled_coords = (coordinates - min_coords) / (max_coords - min_coords + 1e-12)
    
    X = np.diag(scaled_coords[:, 0])
    Y = np.diag(scaled_coords[:, 1])
    
    P_X_P_Y_P = P @ X @ P @ Y @ P
    local_chern = -4.0 * np.pi * np.imag(np.diagonal(P_X_P_Y_P))
    
    # Estrazione del bulk profondo (r < 0.35 R)
    center = coordinates.mean(axis=0)
    distances = np.linalg.norm(coordinates - center, axis=1)
    bulk_mask = distances < (0.35 * distances.max())
    
    mean_bulk_chern = np.mean(local_chern[bulk_mask]) if np.any(bulk_mask) else 0.0
    
    return mean_bulk_chern, local_gap, mean_level_spacing, is_true_gap

if __name__ == "__main__":
    print("=== Pipeline Numerica Unificata e Rigorosa (Opzione B) ==-\n")
    
    N_test = 196
    alpha_test = 0.15
    fillings = [0.27, 0.70]
    n_realizations = 10  # Aumentabile a 20+ per produzione finale
    base_seed = 1337
    
    print(f"{'Filling':<10} | {'Chern Medio (Bulk)':<20} | {'Std Dev':<10} | {'Gap Medio':<12} | {'True Gap Rate':<15}")
    print("-" * 75)
    
    for f in fillings:
        chern_vals = []
        gaps = []
        true_gap_count = 0
        
        for r in range(n_realizations):
            # Seed unico e riproducibile per ogni realizzazione
            run_seed = base_seed + r * 17
            coords = generate_standard_amorphous_lattice(N_test, box_size=np.sqrt(N_test)*1.2, seed=run_seed)
            H_mat = build_tight_binding_hamiltonian(coords, alpha=alpha_test, cutoff=1.8)
            
            c_val, g_val, delta_val, is_gap = calculate_rigorous_diagnostics(H_mat, coords, filling_fraction=f)
            
            chern_vals.append(c_val)
            gaps.append(g_val)
            if is_gap:
                true_gap_count += 1
                
        mean_c = np.mean(chern_vals)
        std_c = np.std(chern_vals)
        mean_g = np.mean(gaps)
        gap_rate = (true_gap_count / n_realizations) * 100.0
        
        print(f"{f*100.0}%       | {mean_c:+.4f}               | {std_c:.4f}     | {mean_g:.4f}       | {gap_rate:.1f}%")

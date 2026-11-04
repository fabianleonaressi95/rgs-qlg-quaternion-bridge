import numpy as np
import scipy.sparse as sp
import scipy.linalg as LA

def generate_amorphous_lattice(N, box_size=10.0, steps=40, dt=0.01):
    """
    Genera un reticolo amorfo/aperiodico tramite rilassamento repulsivo
    di punti casuali (schema standard Agarwala-Shenoy 2017).
    """
    np.random.seed(None)
    coords = np.random.rand(N, 2) * box_size
    for _ in range(steps):
        forces = np.zeros_like(coords)
        for i in range(N):
            for j in range(i + 1, N):
                diff = coords[i] - coords[j]
                dist = np.linalg.norm(diff)
                if dist < 1e-5:
                    dist = 1e-5
                # Forza repulsiva a inversione cubica della distanza
                f = diff / (dist**3)
                forces[i] += f
                forces[j] -= f
        coords += dt * forces
        coords = np.clip(coords, 0, box_size)
    return coords

def build_exponential_peierls_hamiltonian(coordinates, t0=1.0, gamma=2.0, alpha=0.15, cutoff=2.5):
    """
    Costruisce l'Hamiltoniana tight-binding con hopping a decadimento esponenziale 
    e fasi di Peierls (gauge di Landau).
    """
    N = coordinates.shape[0]
    rows, cols, data = [], [], []
    
    for i in range(N):
        for j in range(i + 1, N):
            pos_i = coordinates[i]
            pos_j = coordinates[j]
            dist = np.linalg.norm(pos_i - pos_j)
            
            if dist <= cutoff and dist > 1e-4:
                # Hopping a decadimento esponenziale
                t_ij = t0 * np.exp(-gamma * (dist - 1.0))
                
                # Fase di Peierls (gauge di Landau)
                mid_y = (pos_i[1] + pos_j[1]) / 2.0
                delta_x = pos_j[0] - pos_i[0]
                phi_ij = 2.0 * np.pi * alpha * mid_y * delta_x
                
                t_complex = t_ij * np.exp(1j * phi_ij)
                
                rows.extend([i, j])
                cols.extend([j, i])
                data.extend([t_complex, np.conj(t_complex)])
                
    H = sp.coo_matrix((data, (rows, cols)), shape=(N, N)).tocsr()
    return H

def calculate_local_chern_marker(H, coordinates, filling_fraction=0.7):
    """
    Calcola il Marcatore di Chern Locale (formulazione Bianco-Resta)
    per ciascun sito del sistema aperiodico.
    """
    N = H.shape[0]
    H_dense = H.toarray() if hasattr(H, "toarray") else H
    evals, eigenvectors = LA.eigh(H_dense)
    
    num_occupied = int(N * filling_fraction)
    if num_occupied <= 0 or num_occupied >= N:
        return np.zeros(N), evals
        
    V_occ = eigenvectors[:, :num_occupied]
    P = V_occ @ V_occ.conj().T
    
    # Normalizzazione delle coordinate spaziali in [0, 1]
    min_coords = coordinates.min(axis=0)
    max_coords = coordinates.max(axis=0)
    scaled_coords = (coordinates - min_coords) / (max_coords - min_coords + 1e-12)
    
    X = np.diag(scaled_coords[:, 0])
    Y = np.diag(scaled_coords[:, 1])
    
    # Operatore per il marcatore locale di Chern: C_i = -4\pi * Im(diag(P X P Y P))
    P_X_P_Y_P = P @ X @ P @ Y @ P
    local_chern = -4.0 * np.pi * np.imag(np.diagonal(P_X_P_Y_P))
    
    return local_chern, evals

if __name__ == "__main__":
    print("=== Test di Convergenza e Scaling in N con il Marcatore di Chern Locale ==-\n")
    
    n_values = [100, 196, 400]
    alpha = 0.15
    filling = 0.70
    n_realizations = 5  # Media d'ensemble per abbattere la varianza statistica
    
    print(f"{'N':<6} | {'Marcatore Medio (Bulk r < 0.4R)':<32} | {'Std Dev':<10}")
    print("-" * 55)
    
    for N in n_values:
        bulk_markers = []
        box_side = np.sqrt(N) * 1.2
        
        for _ in range(n_realizations):
            coords = generate_amorphous_lattice(N, box_size=box_side)
            H_mat = build_exponential_peierls_hamiltonian(coords, alpha=alpha, cutoff=2.5)
            local_c, _ = calculate_local_chern_marker(H_mat, coords, filling_fraction=filling)
            
            # Isolamento del bulk profondo per evitare artefatti di bordo
            center = coords.mean(axis=0)
            distances = np.linalg.norm(coords - center, axis=1)
            max_dist = distances.max()
            bulk_mask = distances < 0.4 * max_dist
            
            if np.any(bulk_mask):
                bulk_markers.extend(local_c[bulk_mask])
                
        if bulk_markers:
            mean_val = np.mean(bulk_markers)
            std_val = np.std(bulk_markers)
            print(f"{N:<6d} | {mean_val:+.4f}                         | {std_val:.4f}")

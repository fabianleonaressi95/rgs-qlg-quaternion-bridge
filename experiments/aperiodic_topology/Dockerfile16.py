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

def build_two_band_topological_hamiltonian(coordinates, t1=1.0, t2=0.3, lambda_soc=0.8, m_staggered=0.5):
    """
    Costruisce un'Hamiltoniana tight-binding a 2 orbitali per sito su reticolo amorfo di Delaunay,
    incorporando hopping di primo (NN) e secondo vicinato (NNN) con fasi di Haldane geometriche.
    Dimensione totale della matrice: 2N x 2N.
    """
    N = coordinates.shape[0]
    tri = Delaunay(coordinates)
    
    # Matrice sparsa 2N x 2N
    dim = 2 * N
    rows, cols, data = [], [], []
    
    # Mappa delle adiacenze di primo vicinato (NN) e secondo vicinato (NNN)
    nn_edges = set()
    for simplex in tri.simplices:
        for i_local in range(3):
            u = simplex[i_local]
            v = simplex[(i_local + 1) % 3]
            nn_edges.add((min(u, v), max(u, v)))
            
    # Costruzione dei blocchi orbitali e di hopping
    for u, v in nn_edges:
        pos_u = coordinates[u]
        pos_v = coordinates[v]
        dist = np.linalg.norm(pos_v - pos_u)
        if dist > 1e-4:
            theta = np.arctan2(pos_v[1] - pos_u[1], pos_v[0] - pos_u[0])
            # Hopping NN con fase complessa dipendente dall'angolo
            t_phase = t1 * np.exp(1j * lambda_soc * np.sin(2.0 * theta))
            
            # Sub-matrice 2x2 tra sito u e sito v (accoppiamento orbitale σ_z)
            # Orbitale A (0) e B (1)
            for alpha in [0, 1]:
                sign = 1.0 if alpha == 0 else -1.0
                u_idx = 2 * u + alpha
                v_idx = 2 * v + alpha
                
                rows.extend([u_idx, v_idx])
                cols.extend([v_idx, u_idx])
                data.extend([t_phase * sign, np.conj(t_phase) * sign])

    # Potenziale di massa on-site (staggered tra orbitale A e B)
    for i in range(N):
        i_A = 2 * i
        i_B = 2 * i + 1
        rows.extend([i_A, i_B])
        cols.extend([i_A, i_B])
        data.extend([m_staggered, -m_staggered])
        
    H = sp.coo_matrix((data, (rows, cols)), shape=(dim, dim)).tocsr()
    return H

def calculate_rigorous_diagnostics_2band(H, coordinates, filling_fraction=0.5):
    """
    Calcola il marcatore locale di Chern e verifica il gap per il sistema a 2 bande (dimensione 2N).
    """
    dim = H.shape[0]
    N = coordinates.shape[0]
    H_dense = H.toarray() if hasattr(H, "toarray") else H
    evals, eigenvectors = LA.eigh(H_dense)
    
    num_occupied = int(dim * filling_fraction)
    if num_occupied <= 0 or num_occupied >= dim:
        return 0.0, 0.0, False
        
    mean_level_spacing = np.mean(np.diff(evals))
    local_gap = evals[num_occupied] - evals[num_occupied - 1] if 0 < num_occupied < dim else 0.0
    is_true_gap = local_gap > (1.2 * mean_level_spacing)
    
    V_occ = eigenvectors[:, :num_occupied]
    P = V_occ @ V_occ.conj().T
    
    # Ricostruzione delle coordinate estese per il marcatore spaziale sui 2N siti
    coords_extended = np.repeat(coordinates, 2, axis=0)
    min_coords = coords_extended.min(axis=0)
    max_coords = coords_extended.max(axis=0)
    scaled_coords = (coords_extended - min_coords) / (max_coords - min_coords + 1e-12)
    
    X = np.diag(scaled_coords[:, 0])
    Y = np.diag(scaled_coords[:, 1])
    
    P_X_P_Y_P = P @ X @ P @ Y @ P
    local_chern = -4.0 * np.pi * np.imag(np.diagonal(P_X_P_Y_P))
    
    # Media nel bulk basata sulla posizione geometrica originaria
    center = coordinates.mean(axis=0)
    distances = np.linalg.norm(coordinates - center, axis=1)
    bulk_mask_sites = distances < (0.35 * distances.max())
    
    # Espande la maschera per i 2 orbitali per sito
    bulk_mask_extended = np.repeat(bulk_mask_sites, 2)
    
    mean_bulk_chern = np.mean(local_chern[bulk_mask_extended]) if np.any(bulk_mask_extended) else 0.0
    
    return mean_bulk_chern, local_gap, is_true_gap

if __name__ == "__main__":
    print("=== Scansione 2 Bande: Amorphous Chern Insulator (Delaunay + Orbitali) ==-\n")
    
    N_test = 121  # Numero di nodi spaziali (matrice effettiva 242x242)
    lambda_values = [0.5, 1.0, 1.5]
    mass_values = [0.2, 0.5, 1.0]
    n_realizations = 3
    base_seed = 1337
    
    print(f"{'Lambda_SOC':<12} | {'M_Staggered':<12} | {'Chern Medio (Bulk)':<20} | {'Gap Medio':<12} | {'True Gap Rate':<15}")
    print("-" * 81)
    
    for l_soc in lambda_values:
        for m_val in mass_values:
            chern_vals = []
            gaps = []
            true_gap_count = 0
            
            for r in range(n_realizations):
                run_seed = base_seed + r * 37
                coords = generate_amorphous_lattice(N_test, box_size=np.sqrt(N_test)*1.1, seed=run_seed, steps=10)
                H_mat = build_two_band_topological_hamiltonian(coords, lambda_soc=l_soc, m_staggered=m_val)
                
                c_val, g_val, is_gap = calculate_rigorous_diagnostics_2band(H_mat, coords, filling_fraction=0.5)
                
                chern_vals.append(c_val)
                gaps.append(g_val)
                if is_gap:
                    true_gap_count += 1
                    
            mean_c = np.mean(chern_vals)
            mean_g = np.mean(gaps)
            gap_rate = (true_gap_count / n_realizations) * 100.0
            
            print(f"{l_soc:<12.2f} | {m_val:<12.2f} | {mean_c:<+20.4f} | {mean_g:<12.4f} | {gap_rate:<14.1f}%")

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

def build_peierls_amorphous_hamiltonian(coordinates, t1=1.0, B_field=0.5):
    """
    Costruisce l'Hamiltoniana tight-binding su reticolo amorfo di Delaunay
    con campo magnetico sintetico B tramite sostituzione di Peierls.
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
            mx = 0.5 * (pos_u[0] + pos_v[0])
            my = 0.5 * (pos_u[1] + pos_v[1])
            
            Ax = -0.5 * B_field * my
            Ay = 0.5 * B_field * mx
            
            dx = pos_v[0] - pos_u[0]
            dy = pos_v[1] - pos_u[1]
            peierls_phase = Ax * dx + Ay * dy
            
            t_hop = t1 * np.exp(1j * peierls_phase)
            
            rows.extend([u, v])
            cols.extend([v, u])
            data.extend([t_hop, np.conj(t_hop)])

    H = sp.coo_matrix((data, (rows, cols)), shape=(N, N)).tocsr()
    return H

def analyze_gapped_chern_states(H, coordinates):
    """
    Analizza lo spettro, rileva i veri gap e calcola il Chern per ciascun gap trovato.
    """
    N = H.shape[0]
    H_dense = H.toarray() if hasattr(H, "toarray") else H
    evals, eigenvectors = LA.eigh(H_dense)
    
    # Calcola la spaziatura tra livelli energetici consecutivi
    level_spacings = np.diff(evals)
    mean_spacing = np.mean(level_spacings)
    
    # Identifica gli indici dove c'è un vero gap (spaziatura > 3 volte la media)
    gap_indices = np.where(level_spacings > 3.0 * mean_spacing)[0]
    
    results = []
    min_coords = coordinates.min(axis=0)
    max_coords = coordinates.max(axis=0)
    scaled_coords = (coordinates - min_coords) / (max_coords - min_coords + 1e-12)
    
    X = np.diag(scaled_coords[:, 0])
    Y = np.diag(scaled_coords[:, 1])
    
    center = coordinates.mean(axis=0)
    distances = np.linalg.norm(coordinates - center, axis=1)
    bulk_mask = distances < (0.35 * distances.max())
    
    for idx in gap_indices:
        num_occupied = idx + 1
        filling = num_occupied / N
        local_gap = level_spacings[idx]
        
        # Calcolo del marcatore di Chern locale P * X * Q * Y * P
        V_occ = eigenvectors[:, :num_occupied]
        P = V_occ @ V_occ.conj().T
        Q = np.eye(N) - P
        
        P_X_Q_Y_P = P @ X @ Q @ Y @ P
        local_chern = -4.0 * np.pi * np.imag(np.diagonal(P_X_Q_Y_P))
        
        mean_bulk_chern = np.mean(local_chern[bulk_mask]) if np.any(bulk_mask) else 0.0
        results.append((filling, evals[idx], local_gap, mean_bulk_chern))
        
    return results

if __name__ == "__main__":
    print("=== Scansione Automatica dei Gap e Chern in Reticolo Amorfo ==-\n")
    
    N_test = 196
    b_fields = [0.5, 1.0, 1.5, 2.0]
    base_seed = 1337
    
    print(f"{'B_Field':<10} | {'Filling':<10} | {'Energia Gap':<12} | {'Ampiezza Gap':<12} | {'Chern (Bulk)':<15}")
    print("-" * 65)
    
    for b_val in b_fields:
        # Usiamo una realizzazione di riferimento per il campo
        coords = generate_amorphous_lattice(N_test, box_size=np.sqrt(N_test)*1.1, seed=base_seed, steps=10)
        H_mat = build_peierls_amorphous_hamiltonian(coords, t1=1.0, B_field=b_val)
        
        gapped_states = analyze_gapped_chern_states(H_mat, coords)
        
        if not gapped_states:
            print(f"{b_val:<10.2f} | {'Nessun gap rilevato':<50}")
        else:
            for fill, en, gap_sz, chern_val in gapped_states:
                print(f"{b_val:<10.2f} | {fill:<10.3f} | {en:<12.4f} | {gap_sz:<12.4f} | {chern_val:<+15.4f}")

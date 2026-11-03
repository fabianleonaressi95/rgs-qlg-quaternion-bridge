import numpy as np
import scipy.sparse as sp
import scipy.linalg as LA

def generate_quasiperiodic_lattice(n_cells=10):
    """
    Genera un set di coordinate basato su un reticolo aperiodico bidimensionale 
    (stile cut-and-project / quasi-cristallino a simmetria rotazionale).
    """
    coords = []
    theta = 2.0 * np.pi / 5.0
    for i in range(-n_cells, n_cells + 1):
        for j in range(-n_cells, n_cells + 1):
            x = i + j * np.cos(theta)
            y = j * np.sin(theta)
            if x**2 + y**2 <= (n_cells * 0.7)**2:
                coords.append([x, y])
    return np.array(coords)

def build_peierls_tight_binding(coordinates, t_nominal=1.0, alpha=0.1, cutoff=1.8):
    """
    Costruisce l'Hamiltoniana complessa con fasi di Peierls (gauge di Landau).
    """
    N = coordinates.shape[0]
    rows, cols, data = [], [], []
    
    for i in range(N):
        for j in range(i + 1, N):
            pos_i = coordinates[i]
            pos_j = coordinates[j]
            dist = np.linalg.norm(pos_i - pos_j)
            
            if dist <= cutoff:
                mid_y = (pos_i[1] + pos_j[1]) / 2.0
                delta_x = pos_j[0] - pos_i[0]
                phi_ij = 2.0 * np.pi * alpha * mid_y * delta_x
                t_complex = t_nominal * np.exp(1j * phi_ij)
                
                rows.extend([i, j])
                cols.extend([j, i])
                data.extend([t_complex, np.conj(t_complex)])
                
    H = sp.coo_matrix((data, (rows, cols)), shape=(N, N)).tocsr()
    return H

def calculate_bott_index(H, coordinates, filling_fraction):
    """
    Calcola il Bott Index per una specifica frazione di riempimento.
    """
    N = H.shape[0]
    H_dense = H.toarray() if hasattr(H, "toarray") else H
    eigenvalues, eigenvectors = LA.eigh(H_dense)
    
    num_occupied = int(N * filling_fraction)
    if num_occupied <= 0 or num_occupied >= N:
        return 0, eigenvalues
        
    V_occ = eigenvectors[:, :num_occupied]
    P = V_occ @ V_occ.conj().T
    
    min_coords = coordinates.min(axis=0)
    max_coords = coordinates.max(axis=0)
    scaled_coords = (coordinates - min_coords) / (max_coords - min_coords + 1e-12)
    
    X = scaled_coords[:, 0]
    Y = scaled_coords[:, 1]
    
    U = np.diag(np.exp(2.0j * np.pi * X))
    V_op = np.diag(np.exp(2.0j * np.pi * Y))
    
    Identity = np.eye(N)
    U_bar = P @ U @ P + (Identity - P)
    V_bar = P @ V_op @ P + (Identity - P)
    
    M = V_bar @ U_bar @ V_bar.conj().T @ U_bar.conj().T
    evals_M = LA.eigvals(M)
    log_evals = np.log(evals_M)
    
    bott_index = np.round(np.imag(np.sum(log_evals)) / (2.0 * np.pi))
    return int(bott_index), eigenvalues

if __name__ == "__main__":
    # 1. Generazione della struttura aperiodica coerente
    structured_coords = generate_quasiperiodic_lattice(n_cells=8)
    print(f"Nodi generati nel reticolo aperiodico: {structured_coords.shape[0]}")
    
    # 2. Scansione congiunta di Alpha e Filling Fraction per trovare gap topologici
    alpha_test = 0.15
    H_matrix = build_peierls_tight_binding(structured_coords, t_nominal=1.0, alpha=alpha_test, cutoff=1.8)
    
    print(f"\nScansione dei livelli di riempimento con Alpha = {alpha_test}:")
    print(f"{'Filling (%)':<15} | {'Fermi Energy':<15} | {'Bott Index (B)':<15}")
    print("-" * 51)
    
    for fill in np.linspace(0.2, 0.8, 13):
        bott, evals = calculate_bott_index(H_matrix, structured_coords, filling_fraction=fill)
        fermi_energy = evals[int(structured_coords.shape[0] * fill) - 1] if 0 < int(structured_coords.shape[0] * fill) < len(evals) else 0.0
        print(f"{fill*100:<15.1f} | {fermi_energy:<15.3f} | {bott:<15d}")

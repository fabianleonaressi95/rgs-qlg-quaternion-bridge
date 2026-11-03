import numpy as np
import scipy.sparse as sp
import scipy.linalg as LA

def build_peierls_tight_binding(coordinates, t_nominal=1.0, alpha=0.1, cutoff=1.5):
    """
    Costruisce l'Hamiltoniana di tight-binding con fasi di Peierls (gauge sintetica)
    per un reticolo aperiodico a partire dalle coordinate spaziali dei nodi.
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

def calculate_bott_index(H, coordinates, filling_fraction=0.5):
    """
    Calcola il Bott Index per un sistema aperiodico a partire dall'Hamiltoniana complessa H
    e dalle coordinate spaziali dei nodi.
    """
    N = H.shape[0]
    H_dense = H.toarray() if hasattr(H, "toarray") else H
    eigenvalues, eigenvectors = LA.eigh(H_dense)
    
    num_occupied = int(N * filling_fraction)
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
    # 1. Generazione delle coordinate geometriche di test
    np.random.seed(42)
    num_nodes = 100
    test_coords = np.random.rand(num_nodes, 2) * 10.0
    
    # 2. Generazione dell'Hamiltoniana complessa con fasi di Peierls
    H_matrix = build_peierls_tight_binding(
        coordinates=test_coords, 
        t_nominal=1.0, 
        alpha=0.05, 
        cutoff=2.0
    )
    
    print("Hamiltoniana complessa generata correttamente.")
    print(f"Dimensioni della matrice: {H_matrix.shape}")
    print(f"Numero di elementi di hopping non nulli: {H_matrix.nnz}")
    
    # 3. Calcolo del Bott Index e dello spettro
    bott, evals = calculate_bott_index(H_matrix, test_coords, filling_fraction=0.5)
    
    print(f"Spettro energetico calcolato. Autovalori min/max: {evals.min():.4f} / {evals.max():.4f}")
    print(f"Risultato Bott Index (B): {bott}")

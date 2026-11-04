import numpy as np
import scipy.sparse as sp

def build_peierls_tight_binding(coordinates, t_nominal=1.0, alpha=0.1, cutoff=1.5):
    """
    Costruisce l'Hamiltoniana di tight-binding con fasi di Peierls (gauge sintetica)
    per un reticolo aperiodico a partire dalle coordinate spaziali dei nodi.
    
    Parametri:
    - coordinates: array numpy di forma (N, 2) con le coordinate (x, y) dei nodi.
    - t_nominal: ampiezza di hopping di base (reale).
    - alpha: intensità del flusso magnetico sintetico per plaquette.
    - cutoff: distanza massima di connessione per definire i legami (nearest neighbors).
    
    Ritorna:
    - H: Matrice sparsa CSR complessa rappresentante l'Hamiltoniana del sistema.
    """
    N = coordinates.shape[0]
    rows, cols, data = [], [], []
    
    # Scansione di tutte le coppie di nodi per identificare i legami validi
    for i in range(N):
        for j in range(i + 1, N):
            pos_i = coordinates[i]
            pos_j = coordinates[j]
            dist = np.linalg.norm(pos_i - pos_j)
            
            if dist <= cutoff:
                # Calcolo della fase di Peierls basata su un potenziale vettore in gauge di Landau
                # A_x = -alpha * y, A_y = 0 -> la fase dipende dalla coordinata y del punto medio e dallo spostamento in x
                mid_y = (pos_i[1] + pos_j[1]) / 2.0
                delta_x = pos_j[0] - pos_i[0]
                
                # Fase di Peierls phi_ij = 2 * pi * int(A * dl)
                phi_ij = 2.0 * np.pi * alpha * mid_y * delta_x
                
                # Accoppiamento complesso con rottura della simmetria di inversione temporale
                t_complex = t_nominal * np.exp(1j * phi_ij)
                
                # Inserimento dei termini hermitiani coniugati nella matrice
                rows.extend([i, j])
                cols.extend([j, i])
                data.extend([t_complex, np.conj(t_complex)])
                
    # Assemblaggio della matrice sparsa in formato CSR
    H = sp.coo_matrix((data, (rows, cols)), shape=(N, N)).tocsr()
    return H

if __name__ == "__main__":
    # Esempio di test con coordinate geometriche simulate
    np.random.seed(42)
    num_nodes = 100
    test_coords = np.random.rand(num_nodes, 2) * 10.0
    
    # Generazione dell'Hamiltoniana complessa
    H_matrix = build_peierls_tight_binding(
        coordinates=test_coords, 
        t_nominal=1.0, 
        alpha=0.05, 
        cutoff=2.0
    )
    
    print("Hamiltoniana complessa generata correttamente.")
    print(f"Dimensioni della matrice: {H_matrix.shape}")
    print(f"Numero di elementi di hopping non nulli: {H_matrix.nnz}")

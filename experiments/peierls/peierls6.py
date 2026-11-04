import numpy as np
import scipy.sparse as sp
import scipy.linalg as LA
import matplotlib.pyplot as plt

def generate_quasiperiodic_lattice(n_cells=10):
    """
    Genera il reticolo aperiodico bidimensionale.
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

def build_peierls_tight_binding(coordinates, t_nominal=1.0, alpha=0.15, cutoff=1.8):
    """
    Costruisce l'Hamiltoniana con fasi di Peierls per la rottura della simmetria T.
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
                
    return sp.coo_matrix((data, (rows, cols)), shape=(N, N)).tocsr()

if __name__ == "__main__":
    # 1. Setup del sistema con la configurazione alpha = 0.15
    coords = generate_quasiperiodic_lattice(n_cells=8)
    H_matrix = build_peierls_tight_binding(coords, t_nominal=1.0, alpha=0.15, cutoff=1.8)
    
    # 2. Diagonalizzazione completa per estrarre spettro e autofunzioni
    H_dense = H_matrix.toarray()
    evals, evecs = LA.eigh(H_dense)
    
    # 3. Selezione di uno stato all'interno della prima finestra topologica (es. filling ~27%)
    target_idx = int(coords.shape[0] * 0.27)
    eigenstate = evecs[:, target_idx]
    
    # Calcolo della densità di probabilità spaziale: |psi_i|^2
    prob_density = np.abs(eigenstate)**2
    
    print(f"Nodi analizzati: {coords.shape[0]}")
    print(f"Stato selezionato all'indice: {target_idx}")
    print(f"Energia dello stato (E): {evals[target_idx]:.4f}")
    print(f"Normalizzazione della funzione d'onda: {np.sum(prob_density):.4f}")
    
    # 4. Visualizzazione grafica della mappa di densità spaziale
    plt.figure(figsize=(9, 8))
    sc = plt.scatter(
        coords[:, 0], coords[:, 1], 
        c=prob_density, 
        cmap='inferno', 
        s=150 * (prob_density / prob_density.max()) + 20, 
        edgecolor='k', 
        linewidth=0.4
    )
    cbar = plt.colorbar(sc)
    cbar.set_label('Densità di Probabilità Normalizzata $|\psi_i|^2$', fontsize=11)
    
    plt.title(f'Mappa degli Stati di Bordo - Opzione B ($E = {evals[target_idx]:.3f}$)', fontsize=13)
    plt.xlabel('Coordinata X', fontsize=11)
    plt.ylabel('Coordinata Y', fontsize=11)
    plt.gca().set_aspect('equal')
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.show()

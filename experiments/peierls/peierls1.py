import numpy as np
import scipy.linalg as LA

def calculate_bott_index(H, coordinates, filling_fraction=0.5):
    """
    Calcola il Bott Index per un sistema aperiodico a partire dall'Hamiltoniana complessa H
    e dalle coordinate spaziali dei nodi.
    
    Parametri:
    - H: Matrice Hamiltoniana (può essere convertita in densa per sistemi di piccole/medie dimensioni).
    - coordinates: array numpy di forma (N, 2) con le coordinate (x, y).
    - filling_fraction: frazione di stati occupati (es. 0.5 per metà riempimento).
    
    Ritorna:
    - bott_index: Valore intero (o approssimato a intero) dell'indice topologico.
    - eigenvalues: Autovalori energetici del sistema.
    """
    N = H.shape[0]
    
    # 1. Diagonalizzazione dell'Hamiltoniana
    H_dense = H.toarray() if hasattr(H, "toarray") else H
    eigenvalues, eigenvectors = LA.eigh(H_dense)
    
    # 2. Determinazione del numero di stati occupati (livello di Fermi)
    num_occupied = int(N * filling_fraction)
    
    # Matrice degli autovettori occupati (N x num_occupied)
    V_occ = eigenvectors[:, :num_occupied]
    
    # Proiettore sugli stati occupati P = V_occ * V_occ^dag
    P = V_occ @ V_occ.conj().T
    
    # 3. Normalizzazione delle coordinate spaziali nell'intervallo [0, 1]
    min_coords = coordinates.min(axis=0)
    max_coords = coordinates.max(axis=0)
    scaled_coords = (coordinates - min_coords) / (max_coords - min_coords + 1e-12)
    
    X = scaled_coords[:, 0]
    Y = scaled_coords[:, 1]
    
    # Operatori di posizione unitari sul toro (mappatura esponenziale)
    U = np.diag(np.exp(2.0j * np.pi * X))
    V_op = np.diag(np.exp(2.0j * np.pi * Y))
    
    # 4. Proiezione degli operatori unitari nel sottospazio di Fermi
    Identity = np.eye(N)
    U_bar = P @ U @ P + (Identity - P)
    V_bar = P @ V_op @ P + (Identity - P)
    
    # 5. Calcolo del commutatore proiettato M = V_bar * U_bar * V_bar^dag * U_bar^dag
    M = V_bar @ U_bar @ V_bar.conj().T @ U_bar.conj().T
    
    # 6. Estrazione del Bott Index tramite traccia del logaritmo
    # B = 1 / (2 * pi * i) * tr(log(M)) -> equivalente a sommare le fasi degli autovalori di M
    evals_M = LA.eigvals(M)
    log_evals = np.log(evals_M)
    
    # Somma delle parti immaginarie delle fasi divisa per 2*pi
    bott_index = np.round(np.imag(np.sum(log_evals)) / (2.0 * np.pi))
    
    return int(bott_index), eigenvalues

# Esecuzione della verifica topologica sul sistema precedente
if __name__ == "__main__":
    from __main__ import test_coords, H_matrix  # Riprende i dati dallo script precedente
    
    bott, evals = calculate_bott_index(H_matrix, test_coords, filling_fraction=0.5)
    
    print(f"Spettro energetico calcolato. Autovalori min/max: {evals.min():.4f} / {evals.max():.4f}")
    print(f"Risultato Bott Index (B): {bott}")

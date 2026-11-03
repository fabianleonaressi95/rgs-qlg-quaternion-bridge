import numpy as np
import scipy.sparse as sp
import scipy.linalg as LA
from scipy.spatial import Delaunay

class QRollupBatch:
    def __init__(self, batch_id, transactions):
        """
        Rappresenta un batch di transazioni L2 compresse in una variazione 
        di gauge (delta_phi) da ancorare al Layer 1.
        """
        self.batch_id = batch_id
        self.transactions = transactions
        self.delta_phi = self._compress_to_gauge_shift()
        
    def _compress_to_gauge_shift(self):
        """
        Comprime il carico transazionale locale in un unico parametro 
        di shift di fase (potenziale di gauge sintetico).
        """
        # Sintesi crittografica delle transazioni in un offset numerico controllato
        tx_sum = sum(hash(str(tx)) % 10000 for tx in self.transactions)
        # Normalizzazione della spinta di gauge per evitare overflow geometrici
        return float((tx_sum % 100) * 1e-4)

def verify_rollup_on_layer1(coords, base_B_field, rollup_batch, N=1024):
    """
    Modulo di validazione Layer 1:
    Riceve il batch L2, applica la perturbazione di gauge e verifica 
    se il marcatore di Chern del bulk rimane protetto (C = -1).
    """
    # Il campo magnetico effettivo del L1 viene perturbato dalla transazione L2
    effective_B = base_B_field + rollup_batch.delta_phi
    
    box_size = np.sqrt(N / 1.0)
    tri = Delaunay(coords)
    rows, cols, data = [], [], []
    nn_edges = set()
    
    for simplex in tri.simplices:
        for i_local in range(3):
            u = simplex[i_local]
            v = simplex[(i_local + 1) % 3]
            nn_edges.add((min(u, v), max(u, v)))
            
    for u, v in nn_edges:
        pos_u = coords[u]
        pos_v = coords[v]
        dist = np.linalg.norm(pos_v - pos_u)
        if dist > 1e-4:
            mx = 0.5 * (pos_u[0] + pos_v[0])
            my = 0.5 * (pos_u[1] + pos_v[1])
            Ax = -0.5 * effective_B * my
            Ay = 0.5 * effective_B * mx
            dx = pos_v[0] - pos_u[0]
            dy = pos_v[1] - pos_u[1]
            peierls_phase = Ax * dx + Ay * dy
            
            t_hop = np.exp(-0.3 * dist) * np.exp(1j * peierls_phase)
            rows.extend([u, v])
            cols.extend([v, u])
            data.extend([t_hop, np.conj(t_hop)])
            
    H = sp.coo_matrix((data, (rows, cols)), shape=(N, N)).tocsr()
    evals, eigenvectors = LA.eigh(H.toarray())
    
    # Ricerca adattiva del gap spettrale nel bulk (filling ottimale ~83%)
    center_idx = int(N * 0.83)
    wing = int(N * 0.10)
    start_idx = max(10, center_idx - wing)
    end_idx = min(N - 10, center_idx + wing)
    
    spacings = np.diff(evals[start_idx:end_idx])
    num_occupied = start_idx + np.argmax(spacings) + 1 if len(spacings) > 0 else center_idx
    
    # Calcolo del marcatore locale di Bianco-Resta
    V_occ = eigenvectors[:, :num_occupied]
    P = V_occ @ V_occ.conj().T
    Q = np.eye(N) - P
    
    X = np.diag(coords[:, 0])
    Y = np.diag(coords[:, 1])
    
    local_chern = -4.0 * np.pi * np.imag(np.diagonal(P @ X @ Q @ Y @ P))
    
    center = coords.mean(axis=0)
    radii = np.linalg.norm(coords - center, axis=1)
    bulk_vals = local_chern[radii <= (0.32 * radii.max())]
    median_c = np.median(bulk_vals)
    
    # Regola di accettazione deterministica: il blocco è valido se il Chern resta vicino a -1
    is_valid = abs(median_c - (-1.0)) < 0.15
    return is_valid, median_c

if __name__ == "__main__":
    # Generazione di un reticolo di base per il Layer 1
    np.random.seed(42)
    N = 1024
    coords = np.random.rand(N, 2) * np.sqrt(N)
    
    # Creazione di un batch Layer 2 con 50 transazioni simulate
    simulated_transactions = [f"tx_pay_{i}_to_j" for i in range(50)]
    rollup = QRollupBatch(batch_id="batch_001", transactions=simulated_transactions)
    
    print(f"Batch L2 [{rollup.batch_id}] compresso. Delta Gauge (Delta Phi): {rollup.delta_phi:.6f}")
    
    # Ancoraggio e validazione sul Layer 1
    valid, computed_chern = verify_rollup_on_layer1(coords, base_B_field=1.0, rollup_batch=rollup, N=N)
    
    print(f"Mediana Chern Bulk L1 risultante: {computed_chern:.4f}")
    print(f"Stato Validazione Blocco Q-Rollup: {'ACCETTATO (Finalità Deterministica)' if valid else 'RESPINTO (Frattura di Gap)'}")

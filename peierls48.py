import numpy as np
import scipy.sparse as sp
import scipy.linalg as LA
from scipy.spatial import Delaunay

def generate_amorphous_lattice(N, density=1.0, steps=15, dt=0.01, seed=None):
    """
    Genera un reticolo amorfo bidimensionale con densità costante e regolarizzazione 
    anti-collisione (rilassamento a repulsione) per preservare il gap spettrale.
    """
    box_size = np.sqrt(N / density)
    rng = np.random.default_rng(seed)
    coords = rng.random((N, 2)) * box_size
    for _ in range(steps):
        diff = coords[:, None, :] - coords[None, :, :]
        dist = np.linalg.norm(diff, axis=-1)
        dist = np.maximum(dist, 1e-3)
        np.fill_diagonal(dist, 1.0)
        
        inv_dist3 = 1.0 / (dist**3)
        np.fill_diagonal(inv_dist3, 0.0)
        
        forces = np.sum(diff * inv_dist3[:, :, None], axis=1)
        coords += dt * forces
        coords = np.clip(coords, 0, box_size)
    return coords, box_size

class QRollupBatch:
    def __init__(self, batch_id, transactions):
        """
        Rappresenta un batch di transazioni L2 compresso in una variazione 
        di gauge (delta_phi) da ancorare al Layer 1.
        """
        self.batch_id = batch_id
        self.transactions = transactions
        self.delta_phi = self._compress_to_gauge_shift()
        
    def _compress_to_gauge_shift(self):
        """
        Comprime il carico transazionale locale in un unico parametro 
        di shift di fase (potenziale di gauge sintetico controllato).
        """
        tx_sum = sum(hash(str(tx)) % 10000 for tx in self.transactions)
        return float((tx_sum % 50) * 1e-4)

def verify_rollup_on_layer1(coords, base_B_field, rollup_batch, N=1024):
    """
    Modulo di validazione Layer 1:
    Riceve il batch L2, applica la perturbazione di gauge e verifica 
    se il marcatore di Chern del bulk rimane protetto sul plateau target (C = -1).
    """
    effective_B = base_B_field + rollup_batch.delta_phi
    
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

def run_simulation():
    print("=== Inizializzazione Layer 1: Generazione Reticolo Amorfo Rilassato ===")
    N = 1024
    coords, box_size = generate_amorphous_lattice(N=N, density=1.0, steps=15, seed=42)
    print(f"Reticolo amorfo generato con successo. Nodi totali: {N}\n")
    
    print("=== Creazione Batch Layer 2 (Q-Rollup) ===")
    simulated_transactions = [f"tx_pay_{i}_to_j_amount_{i*1.5}" for i in range(75)]
    rollup = QRollupBatch(batch_id="batch_alpha_001", transactions=simulated_transactions)
    print(f"Batch L2 [{rollup.batch_id}] compresso.")
    print(f"Transazioni incluse: {len(rollup.transactions)}")
    print(f"Variazione di Gauge (Delta Phi): {rollup.delta_phi:.6f}\n")
    
    print("=== Ancoraggio e Validazione sul Layer 1 ===")
    valid, computed_chern = verify_rollup_on_layer1(coords, base_B_field=1.0, rollup_batch=rollup, N=N)
    
    print(f"Mediana Chern Bulk L1 risultante : {computed_chern:.4f}")
    print(f"Stato Validazione Blocco Q-Rollup : {'ACCETTATO (Finalità Deterministica)' if valid else 'RESPINTO (Frattura di Gap)'}")

if __name__ == "__main__":
    run_simulation()

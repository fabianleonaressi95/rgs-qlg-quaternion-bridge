import hashlib
import time
import json
import numpy as np
import scipy.sparse as sp
import scipy.linalg as LA
from scipy.spatial import Delaunay

def evolve_amorphous_coordinates(coords, box_size=12.0, dt=0.01, seed=None):
    """Evolve la posizione dei nodi con un campo di forze repulsive."""
    rng = np.random.default_rng(seed)
    perturbed = coords + rng.normal(0, 0.04, size=coords.shape)
    
    diff = perturbed[:, None, :] - perturbed[None, :, :]
    dist = np.linalg.norm(diff, axis=-1)
    np.fill_diagonal(dist, 1e-5)
    
    inv_dist3 = 1.0 / (dist**3)
    forces = np.sum(diff * inv_dist3[:, :, None], axis=1)
    
    new_coords = perturbed + dt * forces
    return np.clip(new_coords, 0, box_size)

def compute_topological_invariant(coordinates, B_field=1.0, filling_fraction=0.65):
    """
    Costruisce l'Hamiltoniana con Peierls, calcola lo spettro e valuta 
    il marcatore di Chern nel bulk per certificare lo stato topologico della rete.
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
            
            t_hop = np.exp(1j * peierls_phase)
            rows.extend([u, v])
            cols.extend([v, u])
            data.extend([t_hop, np.conj(t_hop)])

    H = sp.coo_matrix((data, (rows, cols)), shape=(N, N)).tocsr()
    evals, eigenvectors = LA.eigh(H.toarray())
    
    num_occupied = int(N * filling_fraction)
    V_occ = eigenvectors[:, :num_occupied]
    P = V_occ @ V_occ.conj().T
    Q = np.eye(N) - P
    
    X = np.diag(coordinates[:, 0])
    Y = np.diag(coordinates[:, 1])
    
    P_X_Q_Y_P = P @ X @ Q @ Y @ P
    local_chern = -4.0 * np.pi * np.imag(np.diagonal(P_X_Q_Y_P))
    
    center = coordinates.mean(axis=0)
    distances = np.linalg.norm(coordinates - center, axis=1)
    bulk_mask = distances < (0.35 * distances.max())
    
    mean_bulk_chern = float(np.mean(local_chern[bulk_mask])) if np.any(bulk_mask) else 0.0
    return mean_bulk_chern

class TopologicalBlock:
    def __init__(self, index, transactions, previous_hash, validator_coordinates):
        self.index = index
        self.timestamp = time.time()
        self.transactions = transactions
        self.previous_hash = previous_hash
        self.validator_coordinates = np.array(validator_coordinates)
        
        # Certificazione topologica basata sul campo di gauge sintetico
        self.chern_invariant = compute_topological_invariant(self.validator_coordinates)
        self.hash = self.calculate_hash()

    def calculate_hash(self):
        block_data = {
            "index": self.index,
            "timestamp": self.timestamp,
            "transactions": self.transactions,
            "previous_hash": self.previous_hash,
            "chern_invariant": round(self.chern_invariant, 3)
        }
        block_string = json.dumps(block_data, sort_keys=True).encode()
        return hashlib.sha256(block_string).hexdigest()

class QProofBlockchain:
    def __init__(self, initial_coordinates):
        self.chain = []
        self.current_coordinates = initial_coordinates
        self.create_genesis_block()

    def create_genesis_block(self):
        genesis_block = TopologicalBlock(
            index=0,
            transactions=["Inizializzazione Q-PROOF (Genesis)"],
            previous_hash="0",
            validator_coordinates=self.current_coordinates
        )
        self.chain.append(genesis_block)

    def add_block(self, transactions, block_seed=None):
        self.current_coordinates = evolve_amorphous_coordinates(
            self.current_coordinates, seed=block_seed
        )
        latest_block = self.chain[-1]
        new_block = TopologicalBlock(
            index=len(self.chain),
            transactions=transactions,
            previous_hash=latest_block.hash,
            validator_coordinates=self.current_coordinates
        )
        self.chain.append(new_block)

    def is_chain_valid(self):
        for i in range(1, len(self.chain)):
            current = self.chain[i]
            previous = self.chain[i - 1]

            if current.hash != current.calculate_hash():
                return False
            if current.previous_hash != previous.hash:
                return False
            # Regola di consenso Q-PROOF: il blocco è valido solo se lo stato 
                # mantiene la fase topologica protetta (Chern vicino a -1)
            if abs(current.chern_invariant - (-1.0)) > 0.15:
                return False
        return True

# --- Esecuzione della Blockchain Q-PROOF ---
if __name__ == "__main__":
    print("=== Avvio Registro a Consenso Topologico (Q-PROOF) ==-\n")
    
    np.random.seed(2026)
    initial_coords = np.random.rand(196, 2) * 12.0  # 196 nodi validatori
    
    q_ledger = QProofBlockchain(initial_coords)
    
    q_ledger.add_block(["Emissione Token Sostenibili", "Sincronizzazione P2P"], block_seed=11)
    q_ledger.add_block(["Smart Contract Energetico", "Validazione Epoch #1"], block_seed=22)
    q_ledger.add_block(["Chiusura Ciclo e Snapshot dello Spettro"], block_seed=33)
    
    print(f"Altezza catena attuale: {len(q_ledger.chain)} blocchi")
    print(f"Validazione crittografica e topologica (Q-PROOF): {q_ledger.is_chain_valid()}\n")
    
    for block in q_ledger.chain:
        print(f"Blocco #{block.index}")
        print(f"  - Hash: {block.hash[:16]}...")
        print(f"  - Invariante di Chern (Bulk): {block.chern_invariant:+.4f}")
        print(f"  - Transazioni: {block.transactions}")
        print("-" * 60)

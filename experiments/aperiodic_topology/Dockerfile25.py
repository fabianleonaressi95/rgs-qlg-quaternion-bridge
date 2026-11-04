import hashlib
import time
import json
import numpy as np
from scipy.spatial import Delaunay

def evolve_amorphous_coordinates(coords, box_size=10.0, dt=0.05, seed=None):
    """
    Evolve le coordinate dei nodi simulando un leggero rilassamento amorfo 
    repulsivo, introducendo uno shift dinamico ad ogni blocco.
    """
    rng = np.random.default_rng(seed)
    # Aggiunge una piccola perturbazione stocastica prima del rilassamento
    perturbed = coords + rng.normal(0, 0.05, size=coords.shape)
    
    diff = perturbed[:, None, :] - perturbed[None, :, :]
    dist = np.linalg.norm(diff, axis=-1)
    np.fill_diagonal(dist, 1e-5)
    
    inv_dist3 = 1.0 / (dist**3)
    forces = np.sum(diff * inv_dist3[:, :, None], axis=1)
    
    new_coords = perturbed + dt * forces
    return np.clip(new_coords, 0, box_size)

class GeometricBlock:
    def __init__(self, index, transactions, previous_hash, validator_coordinates):
        self.index = index
        self.timestamp = time.time()
        self.transactions = transactions
        self.previous_hash = previous_hash
        self.validator_coordinates = np.array(validator_coordinates)
        self.geometric_signature = self._compute_geometric_signature()
        self.hash = self.calculate_hash()

    def _compute_geometric_signature(self):
        """
        Calcola una firma basata sulla distribuzione spaziale e sui momenti 
        di ordine superiore della triangolazione di Delaunay locale.
        """
        tri = Delaunay(self.validator_coordinates)
        # Usa il volume medio dei simplessi di Delaunay come metrica topologica
        simplices_coords = self.validator_coordinates[tri.simplices]
        vectors1 = simplices_coords[:, 1, :] - simplices_coords[:, 0, :]
        vectors2 = simplices_coords[:, 2, :] - simplices_coords[:, 0, :]
        areas = 0.5 * np.abs(vectors1[:, 0] * vectors2[:, 1] - vectors1[:, 1] * vectors2[:, 0])
        
        return float(np.mean(areas) + np.std(areas))

    def calculate_hash(self):
        block_data = {
            "index": self.index,
            "timestamp": self.timestamp,
            "transactions": self.transactions,
            "previous_hash": self.previous_hash,
            "geometric_signature": self.geometric_signature
        }
        block_string = json.dumps(block_data, sort_keys=True).encode()
        return hashlib.sha256(block_string).hexdigest()

class AperiodicBlockchain:
    def __init__(self, initial_coordinates):
        self.chain = []
        self.current_coordinates = initial_coordinates
        self.create_genesis_block()

    def create_genesis_block(self):
        genesis_block = GeometricBlock(
            index=0,
            transactions=["Inizializzazione Registro Aperiodico (Genesis)"],
            previous_hash="0",
            validator_coordinates=self.current_coordinates
        )
        self.chain.append(genesis_block)

    def get_latest_block(self):
        return self.chain[-1]

    def add_block(self, transactions, block_seed=None):
        # Evolve le coordinate spaziali dei nodi per il nuovo blocco
        self.current_coordinates = evolve_amorphous_coordinates(
            self.current_coordinates, seed=block_seed
        )
        
        latest_block = self.get_latest_block()
        new_block = GeometricBlock(
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
        return True

# --- Test della Catena Dinamica ---
if __name__ == "__main__":
    print("=== Avvio Blockchain con Rete Spaziale Dinamica e Aperiodica ==-\n")
    
    np.random.seed(1337)
    initial_coords = np.random.rand(24, 2) * 10.0
    
    geo_ledger = AperiodicBlockchain(initial_coords)
    
    geo_ledger.add_block(["Emissione Smart Contract #1", "Aggiornamento Topologia Nodi"], block_seed=42)
    geo_ledger.add_block(["Convalida Transazione Cross-Chain", "Ricalcolo Spettrale del Bulk"], block_seed=99)
    geo_ledger.add_block(["Chiusura Epoch e Sincronizzazione P2P"], block_seed=123)
    
    print(f"Altezza catena attuale: {len(geo_ledger.chain)} blocchi")
    print(f"Integrità della catena verificata: {geo_ledger.is_chain_valid()}\n")
    
    for block in geo_ledger.chain:
        print(f"Blocco #{block.index}")
        print(f"  - Hash: {block.hash[:16]}...")
        print(f"  - Firma Geometrica (Area Delaunay Media): {block.geometric_signature:.6f}")
        print(f"  - Transazioni: {block.transactions}")
        print("-" * 60)

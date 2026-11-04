import hashlib
import time
import json
import numpy as np
from scipy.spatial import Delaunay

class GeometricBlock:
    """
    Rappresenta un blocco all'interno del registro aperiodico.
    Il blocco include dati di transazione, hash precedente, timestamp 
    e una firma geometrica derivata dalla configurazione spaziale del nodo validatore.
    """
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
        Calcola una firma basata sulla densità locale e sui momenti spaziali 
        delle coordinate del validatore nel reticolo aperiodico.
        """
        center = np.mean(self.validator_coordinates, axis=0)
        distances = np.linalg.norm(self.validator_coordinates - center, axis=1)
        return float(np.sum(distances) / (len(self.validator_coordinates) + 1e-8))

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
    """
    Gestisce la catena di blocchi validati attraverso una topologia di rete aperiodica.
    """
    def __init__(self, node_coordinates):
        self.chain = []
        self.node_coordinates = node_coordinates
        self.create_genesis_block()

    def create_genesis_block(self):
        genesis_block = GeometricBlock(
            index=0,
            transactions=["Inizializzazione Registro Aperiodico (Genesis)"],
            previous_hash="0",
            validator_coordinates=self.node_coordinates
        )
        self.chain.append(genesis_block)

    def get_latest_block(self):
        return self.chain[-1]

    def add_block(self, transactions):
        latest_block = self.get_latest_block()
        new_block = GeometricBlock(
            index=len(self.chain),
            transactions=transactions,
            previous_hash=latest_block.hash,
            validator_coordinates=self.node_coordinates
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

# --- Test del Sistema ---
if __name__ == "__main__":
    print("=== Inizializzazione della Blockchain Geometrica e Aperiodica ==-\n")
    
    # Generazione di coordinate spaziali casuali per i nodi validatori (simulazione reticolo amorfo)
    np.random.seed(42)
    node_coords = np.random.rand(16, 2) * 10.0
    
    # Creazione della blockchain
    geo_ledger = AperiodicBlockchain(node_coords)
    
    # Aggiunta di blocchi con transazioni di prova
    geo_ledger.add_block(["Trasferimento di 50 Token Q-Proof", "Registrazione Stato Spettrale #1"])
    geo_ledger.add_block(["Validazione Nodo a Distanza Ottimizzata", "Esecuzione Smart Contract Geometrico"])
    
    print(f"Altezza catena attuale: {len(geo_ledger.chain)} blocchi")
    print(f"Integrità della catena verificata: {geo_ledger.is_chain_valid()}\n")
    
    for block in geo_ledger.chain:
        print(f"Blocco #{block.index}")
        print(f"  - Hash: {block.hash[:16]}...")
        print(f"  - Firma Geometrica: {block.geometric_signature:.4f}")
        print(f"  - Transazioni: {block.transactions}")
        print("-" * 50)

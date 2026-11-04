import os
import asyncio
import json
import time

# Configurazione iniziale dai parametri d'ambiente Docker
NODE_ID = int(os.environ.get("NODE_ID", 0))
TOTAL_NODES = int(os.environ.get("TOTAL_NODES", 64))
PORT_BASE = 8000

def get_aperiodic_neighbors(node_id, total_nodes):
    """
    Calcola i vicini aperiodici basati su offset geometrici derivati 
    dalla sezione aurea (simmetria aperiodica).
    """
    phi = 1.618033988749895
    offsets = [
        int(round(phi * 1)) % total_nodes,
        int(round(phi * 2)) % total_nodes,
        int(round(phi * 3)) % total_nodes,
        int(round(phi * 5)) % total_nodes
    ]
    # Evita duplicati o autoconnessioni
    neighbors = []
    for offset in offsets:
        neighbor_id = (node_id + offset) % total_nodes
        if neighbor_id != node_id and neighbor_id not in neighbors:
            neighbors.append(neighbor_id)
    return sorted(neighbors)

NEIGHBORS = get_aperiodic_neighbors(NODE_ID, TOTAL_NODES)

class P2PNode:
    def __init__(self, node_id, port):
        self.node_id = node_id
        self.port = port
        self.received_messages = 0

    async def handle_connection(self, reader, writer):
        """Gestisce i messaggi in arrivo dai nodi vicini (Server TCP)."""
        data = await reader.read(1024)
        if data:
            try:
                message = json.loads(data.decode())
                self.received_messages += 1
            except json.JSONDecodeError:
                pass
        writer.close()
        await writer.drain()

    async def start_server(self):
        """Avvia il server TCP asincrono del nodo."""
        server = await asyncio.start_server(
            self.handle_connection, '0.0.0.0', self.port
        )
        addr = server.sockets[0].getsockname()
        return server

    async def send_message_to_neighbor(self, neighbor_id):
        """Invia un messaggio a un vicino aperiodico (Client TCP)."""
        target_port = PORT_BASE + neighbor_id
        target_host = f"node_{neighbor_id:02d}"
        
        payload = json.dumps({
            "sender_id": self.node_id,
            "timestamp": time.time()
        })

        for _ in range(3):  # Tentativi di connessione con retry
            try:
                reader, writer = await asyncio.open_connection(target_host, target_port)
                writer.write(payload.encode())
                await writer.drain()
                writer.close()
                await writer.wait_closed()
                return True
            except (ConnectionRefusedError, OSError):
                await asyncio.sleep(0.1)
        return False

    async def run_benchmark(self):
        """Esegue il ciclo di scambio messaggi e calcola le metriche di throughput."""
        print(f"[Node {self.node_id}] Inizializzato. Vicini aperiodici mappati: {NEIGHBORS}")
        
        server = await self.start_server()
        
        # Breve attesa per assicurarsi che tutti i container abbiano avviato i server
        await asyncio.sleep(2.0)

        iterations = 200
        start_time = time.time()

        for _ in range(iterations):
            tasks = [self.send_message_to_neighbor(n) for n in NEIGHBORS]
            await asyncio.gather(*tasks)

        duration = time.time() - start_time
        throughput = round(iterations / duration, 2) if duration > 0 else 0

        result = {
            "node_id": self.node_id,
            "total_nodes": TOTAL_NODES,
            "iterations": iterations,
            "duration_seconds": round(duration, 4),
            "throughput_ops_sec": throughput,
            "status": "BENCHMARK_COMPLETE"
        }
        print(json.dumps(result))

        # Mantieni il server attivo per un breve periodo di chiusura pulita
        server.close()
        await server.wait_closed()

if __name__ == "__main__":
    node = P2PNode(NODE_ID, PORT_BASE + NODE_ID)
    asyncio.run(node.run_benchmark())

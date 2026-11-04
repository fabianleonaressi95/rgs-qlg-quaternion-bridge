"""
topology.py
Gestione delle topologie di rete: Ring, Small-World (Watts-Strogatz), Scale-Free (Barabási-Albert) e Erdős-Rényi.
Utilizza NetworkX per la generazione pulita dei grafi.
"""

import networkx as nx
import numpy as np

def build_topology(n_nodes=200, topology_type="small_world", **kwargs):
    """
    Costruisce la lista di adiacenza per la topologia specificata.
    Ritorna un dizionario: {node_id: [lista_vicini]}
    """
    topology_type = topology_type.lower()
    
    if topology_type == "ring":
        k = kwargs.get("k", 2)
        G = nx.watts_strogatz_graph(n=n_nodes, k=k, p=0.0, seed=kwargs.get("seed"))
    elif topology_type == "small_world":
        k = kwargs.get("k", 4)
        p = kwargs.get("p", 0.05)
        G = nx.watts_strogatz_graph(n=n_nodes, k=k, p=p, seed=kwargs.get("seed"))
    elif topology_type == "scale_free":
        m = kwargs.get("m", 2)
        G = nx.barabasi_albert_graph(n=n_nodes, m=m, seed=kwargs.get("seed"))
    elif topology_type == "erdos_renyi":
        p_edge = kwargs.get("p_edge", 0.03)
        G = nx.erdos_renyi_graph(n=n_nodes, p=p_edge, seed=kwargs.get("seed"))
    else:
        raise ValueError(f"Topologia non supportata: {topology_type}")

    # Converte in dizionario di adiacenza puro per compatibilità con i simulatori
    neighbors = {i: list(G.neighbors(i)) for i in range(n_nodes)}
    return neighbors

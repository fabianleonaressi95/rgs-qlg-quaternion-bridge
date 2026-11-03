"""
correlated_fault_simulator.py

Modello di fault correlato/coordinato, alternativo al fault indipendente
di simulator.py. Invece che ogni nodo si guasti per conto suo con
probabilita' fissa p, il guasto NASCE in un punto della rete e SI PROPAGA
ai vicini nella topologia (modello epidemico SI: Susceptible -> Infected,
nessuna guarigione). E' l'analogo di:

- una partizione di rete che si allarga
- un attacco che compromette nodi adiacenti in sequenza
- un guasto hardware/software che si propaga lungo dipendenze condivise

Dal punto di vista della fisica statistica, questo e' molto piu' vicino a
una transizione di percolazione (che ha una soglia critica netta, legata
alla topologia del grafo) che a un semplice rumore diluito nodo per nodo.
L'ipotesi da testare: con fault correlati, chi mostra un preavviso reale
prima del collasso, a differenza di quanto trovato con fault indipendenti?

Due varianti dei nodi "infetti" (faulty):
- adversarial=False: il nodo infetto vota a caso (rumore puro)
- adversarial=True: il nodo infetto vota sempre l'opposto della maggioranza
  globale corrente (comportamento bizantino coordinato, il caso peggiore)
"""

import numpy as np
from simulator import ToyConsensusNetwork


class ContagionFaultNetwork(ToyConsensusNetwork):
    def __init__(self, n_nodes=200, coupling=0.15, spread_prob=0.05,
                 adversarial=False, seed=None):
        super().__init__(n_nodes=n_nodes, coupling=coupling, seed=seed)
        self.spread_prob = spread_prob
        self.adversarial = adversarial
        self.faulty = np.zeros(self.n, dtype=bool)

    def seed_infection(self, n_seeds=1):
        seeds = self.rng.choice(self.n, size=n_seeds, replace=False)
        self.faulty[seeds] = True

    def _spread(self):
        """Un round di propagazione epidemica SI lungo i link della topologia."""
        newly_faulty = self.faulty.copy()
        for i in range(self.n):
            if self.faulty[i]:
                continue
            infected_neighbors = sum(1 for j in self.neighbors[i] if self.faulty[j])
            if infected_neighbors > 0:
                p_infect = 1 - (1 - self.spread_prob) ** infected_neighbors
                if self.rng.random() < p_infect:
                    newly_faulty[i] = True
        self.faulty = newly_faulty

    def step(self, fault_prob=None):
        """
        fault_prob e' ignorato (accettato solo per compatibilita' di
        interfaccia con il resto del codice) — qui il "fault" e' guidato
        dalla propagazione epidemica, non da una probabilita' indipendente.
        """
        self._spread()
        new_votes = np.copy(self.votes)

        for i in range(self.n):
            if self.faulty[i]:
                if self.adversarial:
                    # CORREZIONE: un nodo bizantino reale non conosce la
                    # maggioranza globale, riceve solo messaggi dai propri
                    # vicini. Usiamo quindi la maggioranza LOCALE stimata
                    # dai vicini (stessa informazione disponibile a un nodo
                    # onesto), non la media globale — che nessun nodo reale
                    # potrebbe osservare in un protocollo a messaggi.
                    neigh = self.neighbors[i]
                    local_mean = np.mean([self.votes[j] for j in neigh])
                    local_majority = 1 if local_mean >= 0.5 else 0
                    new_votes[i] = 1 - local_majority
                else:
                    new_votes[i] = self.rng.integers(0, 2)
                continue
            neigh = self.neighbors[i]
            local_mean = np.mean([self.votes[j] for j in neigh])
            if self.rng.random() < self.coupling:
                new_votes[i] = 1 if local_mean >= 0.5 else 0

        self.votes = new_votes
        return self.votes.copy()

    def fault_fraction(self):
        return float(np.mean(self.faulty))


def run_contagion(n_nodes, n_rounds, spread_prob, adversarial=False,
                   coupling=0.15, n_seed_nodes=1, seed=None):
    """Ritorna (M_Q series, fault_fraction series) per l'intera run."""
    net = ContagionFaultNetwork(n_nodes=n_nodes, coupling=coupling,
                                 spread_prob=spread_prob, adversarial=adversarial, seed=seed)
    net.seed_infection(n_seeds=n_seed_nodes)

    mq_series = np.zeros(n_rounds)
    fault_series = np.zeros(n_rounds)
    for t in range(n_rounds):
        net.step()
        mq_series[t] = net.order_parameter()
        fault_series[t] = net.fault_fraction()
    return mq_series, fault_series


if __name__ == "__main__":
    mq, ff = run_contagion(n_nodes=200, n_rounds=300, spread_prob=0.06, seed=0)
    print("fault_fraction ogni 30 round:", np.round(ff[::30], 3))
    print("M_Q ogni 30 round:          ", np.round(mq[::30], 3))

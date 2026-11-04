"""
criticality.py

Il "Criticality Engine" minimale del middleware. Prende una finestra
scorrevole di valori M_Q (order parameter, uno per round) e calcola le
metriche di criticita' collettiva descritte nel documento originale:

- M_Q medio nella finestra
- chi (susceptibility): N * Var(M_Q) nella finestra -> misura quanto il
  sistema fluttua. In una transizione di fase, chi ha un picco vicino al
  punto critico.
- Binder cumulant U = 1 - <M^4> / (3 <M^2>^2) -> vicino a 2/3 in fase
  ordinata, cala vicino al punto critico. Usato in fisica statistica per
  localizzare punti critici indipendentemente dalla dimensione del sistema.
- autocorrelation lag-1 -> "critical slowing down": vicino alla transizione
  il sistema impiega piu' tempo a rilassare, quindi l'autocorrelazione sale.

Questo modulo NON fa previsioni (quello e' compito del Prediction Engine,
volutamente fuori scope per l'MVP): calcola solo le metriche, cosi' si puo'
verificare empiricamente se hanno senso prima di costruirci sopra un modello
predittivo.
"""

import numpy as np


class CriticalityEngine:
    def __init__(self, n_nodes, window_size=50):
        self.n_nodes = n_nodes
        self.window_size = window_size

    def compute(self, mq_series):
        """
        mq_series: array 1D di M_Q per round, gia' ritagliato alla finestra
        di interesse (es. ultimi `window_size` round).

        Ritorna un dict con le metriche di criticita'.
        """
        mq_series = np.asarray(mq_series)
        if len(mq_series) < 5:
            raise ValueError("Servono almeno 5 round per stimare le metriche")

        mean_mq = float(np.mean(mq_series))
        var_mq = float(np.var(mq_series))
        chi = self.n_nodes * var_mq

        m2 = np.mean(mq_series ** 2)
        m4 = np.mean(mq_series ** 4)
        binder = 1.0 - (m4 / (3 * m2 ** 2)) if m2 > 0 else 0.0

        autocorr = self._lag1_autocorrelation(mq_series)

        return {
            "M_Q": mean_mq,
            "chi": chi,
            "binder": binder,
            "autocorr_lag1": autocorr,
        }

    @staticmethod
    def _lag1_autocorrelation(series):
        if len(series) < 3:
            return 0.0
        x = series - np.mean(series)
        num = np.sum(x[:-1] * x[1:])
        den = np.sum(x ** 2)
        return float(num / den) if den > 0 else 0.0

    def sliding_metrics(self, full_series, step=5):
        """
        Applica compute() su finestre scorrevoli lungo tutta la serie.
        Utile per vedere come evolvono le metriche round dopo round,
        non solo su una singola finestra statica.
        """
        results = []
        for start in range(0, len(full_series) - self.window_size, step):
            window = full_series[start:start + self.window_size]
            m = self.compute(window)
            m["round"] = start + self.window_size
            results.append(m)
        return results

import numpy as np

# ==========================================
# 1. SIMULATORE 0D DELLA MICROTURBINA
# ==========================================
class MicroTurbineSimulator:
    """
    Simulatore ridotto con dinamica accoppiata per la verifica dello scaling non additivo.
    Stato x = [omega (krpm), EGT (°C), W_f (g/s)]
    """
    def __init__(self, dt=0.1):
        self.dt = dt
        self.inertia = 0.001
        
    def step(self, x, u, degradation_factor=1.0, perturbation_vector=np.zeros(3)):
        omega, egt, wf = x
        dot_m_f, p_load = u
        
        # d rappresenta la degradazione globale (es. efficienza compressore/turbina)
        eff_deg = 1.0 - degradation_factor
        
        # Dinamica nominale con accoppiamento non lineare ed effetto della perturbazione epsilon (G)
        d_omega = (1.0 / self.inertia) * (25.0 * dot_m_f * (1.0 - eff_deg**2) - 0.15 * omega - 0.05 * p_load) + perturbation_vector[0]
        d_egt = 12.0 * dot_m_f - 0.08 * (omega * (1.0 - eff_deg**2)) + 0.02 * egt + perturbation_vector[1]
        d_wf = 5.0 * (dot_m_f - wf) + perturbation_vector[2]
        
        omega_next = np.clip(omega + d_omega * self.dt, 30.0, 130.0)
        egt_next = np.clip(egt + d_egt * self.dt, 400.0, 800.0)
        wf_next = np.clip(wf + d_wf * self.dt, 1.0, 10.0)
        
        return np.array([omega_next, egt_next, wf_next])

# ==========================================
# 2. LIFTING DI KOOPMAN EDMD
# ==========================================
class KoopmanEDMD:
    """
    Costruisce lo spazio delle osservabili Psi(x) e calcola l'operatore K.
    """
    def __init__(self):
        pass

    def _lift(self, X):
        n_samples, dim = X.shape
        lifted = [X]
        quads = [X[:, i:i+1] * X[:, j:j+1] for i in range(dim) for j in range(i, dim)]
        return np.hstack(lifted + quads)

    def fit_and_get_spectrum(self, X, X_next):
        Psi_X = self._lift(X)
        Psi_X_next = self._lift(X_next)
        Psi_X_pinv = np.linalg.pinv(Psi_X)
        K = (Psi_X_pinv @ Psi_X_next).T
        evals = np.linalg.eigvals(K)
        return K, evals

# ==========================================
# 3. PIPELINE DI AUDIT E SCALING NON-ADDITIVO
# ==========================================
if __name__ == "__main__":
    print("=== ASTRA-Q: Audit Numerico dello Scaling Non-Additivo p(d) ===")
    
    sim = MicroTurbineSimulator(dt=0.1)
    koopman = KoopmanEDMD()
    
    # Condizioni operative di test
    x_init = np.array([85.0, 580.0, 3.5])
    u_test = np.array([3.2, 11.0])
    
    # 1. Griglia di Degradazione d e Ampiezze di Perturbazione epsilon
    degradation_levels = [0.0, 0.05, 0.10, 0.15]
    epsilons = np.array([1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1])
    
    print(f"\n{'Degradazione (d)':<18} | {'Esponente Stimato p(d)':<22} | {'Stato Verifica'}")
    print("-" * 65)
    
    for d in degradation_levels:
        non_add_indices = []
        
        for eps in epsilons:
            # Vettori di perturbazione separati per i settori dinamici g, y, lambda
            G_g = eps * np.array([1.0, 0.0, 0.0])
            G_y = eps * np.array([0.0, 1.0, 0.0])
            G_l = eps * np.array([0.0, 0.0, 1.0])
            G_total = G_g + G_y + G_l
            
            # Evoluzione singola perturbazione separata (o settori)
            x_next_base = sim.step(x_init, u_test, degradation_factor=d)
            x_next_g = sim.step(x_init, u_test, degradation_factor=d, perturbation_vector=G_g)
            x_next_y = sim.step(x_init, u_test, degradation_factor=d, perturbation_vector=G_y)
            x_next_l = sim.step(x_init, u_test, degradation_factor=d, perturbation_vector=G_l)
            
            # Evoluzione con perturbazione combinata (Delta_C)
            x_next_comb = sim.step(x_init, u_test, degradation_factor=d, perturbation_vector=G_total)
            
            # Calcolo delle risposte Delta nello spazio delle osservabili
            psi_base = koopman._lift(x_init.reshape(1, -1))
            psi_comb = koopman._lift(x_next_comb.reshape(1, -1))
            delta_C = np.linalg.norm(psi_comb - psi_base)
            
            # Somma lineare delle risposte dei singoli settori (Delta_sum)
            psi_g = koopman._lift(x_next_g.reshape(1, -1))
            psi_y = koopman._lift(x_next_y.reshape(1, -1))
            psi_l = koopman._lift(x_next_l.reshape(1, -1))
            
            delta_sum = np.linalg.norm(psi_g - psi_base) + np.linalg.norm(psi_y - psi_base) + np.linalg.norm(psi_l - psi_base)
            
            # Indice di non-additività I_nonadd
            i_nonadd = np.abs(delta_C - delta_sum)
            non_add_indices.append(max(i_nonadd, 1e-12)) # Evita log(0)
            
        # Stima dell'esponente p tramite regressione log-log su epsilon vs I_nonadd
        coeffs = np.polyfit(np.log(epsilons), np.log(non_add_indices), 1)
        p_est = coeffs[0]
        
        status = "CONVERGENTE (~2.0)" if abs(p_est - 2.0) < 0.15 else "DEVIAZIONE"
        print(f"d = {d:.2f}              | p = {p_est:.4f}               | {status}")

    print("\n=== Audit Conclusivo Completato Con Successo ===")

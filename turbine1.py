import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 1. PARAMETRI E SIMULAZIONE DINAMICA OTTIMIZZATA
# ==========================================
class MicroTurbineSimulator:
    """
    Simulatore 0D ridotto della microturbina con accoppiamento non lineare 
    della degradazione (ispirato al modello GTM 140).
    """
    def __init__(self, dt=0.1):
        self.dt = dt
        self.inertia = 0.001
        
    def step(self, x, u, degradation_factor=1.0):
        omega, egt, wf = x
        dot_m_f, p_load = u
        
        # Introduciamo una risposta non lineare/quadratica della degradazione sul rendimento
        eff_deg = degradation_factor ** 2  
        
        d_omega = (1.0 / self.inertia) * (25.0 * dot_m_f * eff_deg - 0.15 * omega - 0.05 * p_load)
        d_egt = 12.0 * dot_m_f - 0.08 * (omega * eff_deg) + 0.02 * egt
        d_wf = 5.0 * (dot_m_f - wf)
        
        omega_next = np.clip(omega + d_omega * self.dt, 30.0, 130.0)
        egt_next = np.clip(egt + d_egt * self.dt, 400.0, 800.0)
        wf_next = np.clip(wf + d_wf * self.dt, 1.0, 10.0)
        
        return np.array([omega_next, egt_next, wf_next])

    def generate_trajectory(self, steps, u_profile, degradation_factor=1.0):
        trajectory = np.zeros((steps, 3))
        x = np.array([80.0, 550.0, 3.0]) 
        for t in range(steps):
            trajectory[t] = x
            u = u_profile[t]
            x = self.step(x, u, degradation_factor)
        return trajectory

# ==========================================
# 2. LIFTING DI KOOPMAN ESTESO (EDMD)
# ==========================================
class KoopmanEDMD:
    """
    EDMD con dizionario esteso (polinomiale fino al 3° grado e termini trigonometrici).
    """
    def __init__(self):
        self.K = None

    def _lift(self, X):
        n_samples, dim = X.shape
        lifted = [X]
        # Termini quadratici e cubici per catturare la dinamica non lineare
        quads = [X[:, i:i+1] * X[:, j:j+1] for i in range(dim) for j in range(i, dim)]
        cubics = [X[:, i:i+1] * X[:, j:j+1] * X[:, k:k+1] for i in range(dim) for j in range(i, dim) for k in range(j, dim)]
        # Termini trigonometrici per le armoniche
        trig = [np.sin(X), np.cos(X)]
        return np.hstack(lifted + quads + cubics + trig)

    def fit(self, X, X_next):
        Psi_X = self._lift(X)
        Psi_X_next = self._lift(X_next)
        Psi_X_pinv = np.linalg.pinv(Psi_X)
        self.K = (Psi_X_pinv @ Psi_X_next).T
        return self.K

# ==========================================
# 3. ESECUZIONE DEL BENCHMARK OTTIMIZZATO
# ==========================================
if __name__ == "__main__":
    print("=== ASTRA-Q MicroTurbine Benchmark Ottimizzato ===")
    
    sim = MicroTurbineSimulator(dt=0.1)
    steps = 1500
    t_array = np.linspace(0, steps * 0.1, steps)
    
    u_profile = np.zeros((steps, 2))
    u_profile[:, 0] = 3.0 + 1.5 * np.sin(0.05 * t_array) + 0.5 * np.cos(0.2 * t_array)
    u_profile[:, 1] = 10.0 + 2.0 * np.sin(0.02 * t_array)

    traj_healthy = sim.generate_trajectory(steps, u_profile, degradation_factor=1.0)
    
    degradation_level = 0.85 
    traj_degraded = sim.generate_trajectory(steps, u_profile, degradation_factor=degradation_level)

    koopman = KoopmanEDMD()
    K_healthy = koopman.fit(traj_healthy[:-1], traj_healthy[1:])
    
    koopman_deg = KoopmanEDMD()
    K_degraded = koopman_deg.fit(traj_degraded[:-1], traj_degraded[1:])

    evals_healthy, _ = np.linalg.eig(K_healthy)
    evals_degraded, _ = np.linalg.eig(K_degraded)
    spectral_shift = np.linalg.norm(np.sort(np.abs(evals_degraded)) - np.sort(np.abs(evals_healthy)))

    print(f"\n[Risultati Analisi Spettrale ASTRA-Q]:")
    print(f" - Variazione spettrale norma autovalori (Health Indicator): {spectral_shift:.5f}")
    print(f" - Condizione motore: Degradato al {(1.0 - degradation_level)*100:.1f}%")

    # Verifica della legge di scaling con perturbazioni ridotte
    epsilons = [0.01, 0.02, 0.03, 0.05, 0.08, 0.10]
    non_add_indices = []
    
    for eps in epsilons:
        traj_pert = sim.generate_trajectory(200, u_profile[:200], degradation_factor=(1.0 - eps))
        diff = np.linalg.norm(traj_pert - traj_healthy[:200]) / np.linalg.norm(traj_healthy[:200])
        non_add_indices.append(diff)

    regression_coeffs = np.polyfit(np.log(epsilons), np.log(non_add_indices), 1)
    p_estimated = regression_coeffs[0]

    print(f"\n[Verifica Scaling Legge di Perturbazione]:")
    print(f" - Esponente stimato p empirico: {p_estimated:.4f} (Target teorico ~ 2.0)")
    print("=== Benchmark Ottimizzato Completato ===")

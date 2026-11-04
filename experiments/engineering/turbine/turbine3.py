import numpy as np

# ==========================================
# 1. SIMULATORE CON DIPENDENZE NON LINEARI ESPLICITE
# ==========================================
class MicroTurbineSimulatorNonlinear:
    def __init__(self, dt=0.1):
        self.dt = dt
        self.inertia = 0.001
        
    def step(self, x, u, degradation_factor=1.0):
        omega, egt, wf = x
        dot_m_f, p_load = u
        
        # Introduciamo una dipendenza non lineare di ordine superiore (es. effetto a soglia quadratica/cubica)
        # epsilon rappresenta la perdita: epsilon = 1.0 - degradation_factor
        eps = 1.0 - degradation_factor
        eff_deg = 1.0 - (eps + 2.5 * (eps ** 2)) # Dinamica non additiva marcata
        
        d_omega = (1.0 / self.inertia) * (25.0 * dot_m_f * eff_deg - 0.15 * omega - 0.05 * p_load)
        d_egt = 12.0 * dot_m_f - 0.08 * (omega * eff_deg) + 0.02 * egt
        d_wf = 5.0 * (dot_m_f - wf)
        
        omega_next = np.clip(omega + d_omega * self.dt, 30.0, 130.0)
        egt_next = np.clip(egt + d_egt * self.dt, 400.0, 800.0)
        wf_next = np.clip(wf + d_wf * self.dt, 1.0, 10.0)
        
        return np.array([omega_next, egt_next, wf_next])

# ==========================================
# 2. DIZIONARIO DI KOOPMAN ESTESO CON TERMINI DI GRADO SUPERIORE
# ==========================================
class KoopmanLiftAdvanced:
    def transform(self, X):
        n_samples, dim = X.shape
        lifted = [X]
        # Includiamo termini quadratici e cubici per catturare la curvatura esatta
        quads = [X[:, i:i+1] * X[:, j:j+1] for i in range(dim) for j in range(i, dim)]
        cubics = [X[:, i:i+1] * X[:, j:j+1] * X[:, k:k+1] for i in range(dim) for j in range(i, dim) for k in range(j, dim)]
        return np.hstack(lifted + quads + cubics)

# ==========================================
# 3. VERIFICA DELLA LEGGE DI SCALING NON LINEARE
# ==========================================
if __name__ == "__main__":
    print("=== ASTRA-Q Scaling con Dipendenze Non Lineari (Target p = 2.0) ===")
    
    sim = MicroTurbineSimulatorNonlinear(dt=0.1)
    lifter = KoopmanLiftAdvanced()
    
    x_init = np.array([85.0, 580.0, 3.5])
    u_test = np.array([3.2, 11.0])
    
    x_next_healthy = sim.step(x_init, u_test, degradation_factor=1.0)
    psi_healthy = lifter.transform(x_next_healthy.reshape(1, -1))

    epsilons = [0.005, 0.01, 0.02, 0.04, 0.08, 0.15]
    non_add_indices = []
    
    for eps in epsilons:
        deg_factor = 1.0 - eps
        x_next_pert = sim.step(x_init, u_test, degradation_factor=deg_factor)
        psi_pert = lifter.transform(x_next_pert.reshape(1, -1))
        
        diff_norm = np.linalg.norm(psi_pert - psi_healthy) / np.linalg.norm(psi_healthy)
        non_add_indices.append(diff_norm)

    regression_coeffs = np.polyfit(np.log(epsilons), np.log(non_add_indices), 1)
    p_estimated = regression_coeffs[0]

    print(f"\n[Risultati Analisi di Scaling con Dipendenze Avanzate]:")
    for eps, val in zip(epsilons, non_add_indices):
        print(f" - Perturbazione epsilon = {eps:.3f} -> Variazione osservabili = {val:.6f}")
    
    print(f"\n -> Esponente di scaling p stimato: {p_estimated:.4f} (Target teorico ~ 2.0)")
    print("=== Test Conclusivo con Dipendenze Completato ===")

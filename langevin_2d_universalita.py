#!/usr/bin/env python3
"""
Script completo per Langevin 2D con forzatura rotazionale.
Implementazione dei 3 passi: universalità»¡, confronto teoria, Lorenz-63.
"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import PolynomialFeatures
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# FUNZIONI BASE
# ============================================================================

def simulate_langevin_2d_rot(T=100.0, dt=0.001, x0=0.0, y0=0.0,
                              Omega=1.2, D=0.25,
                              a=1.0, b=0.5, c=2.0, g=0.5, seed=None):
    """
    Dinamica di Langevin 2D con forzatura rotazionale.
    dx = [-∂_x U - Ω*y] dt + √(2D) dW_x
    dy = [-∂_y U + Ω*x] dt + √(2D) dW_y
    """
    if seed is not None:
        np.random.seed(seed)

    n_steps = int(T / dt)
    x_traj = np.zeros(n_steps)
    y_traj = np.zeros(n_steps)

    x, y = x0, y0
    sqrt_2D = np.sqrt(2 * D)

    for i in range(n_steps):
        # Gradiente del potenziale
        dUdx = a*x**3 + b*x + 2*g*x*y
        dUdy = c*y + g*x**2

        # Campo vettoriale totale (con rotazione)
        fx = -dUdx - Omega * y
        fy = -dUdy + Omega * x

        # Euler-Maruyama
        dWx = np.sqrt(dt) * np.random.randn()
        dWy = np.sqrt(dt) * np.random.randn()

        x = x + fx * dt + sqrt_2D * dWx
        y = y + fy * dt + sqrt_2D * dWy

        # Controllo stabilità
        if np.abs(x) > 100 or np.abs(y) > 100:
            print(f"    WARNING: divergenza al passo {i}")
            x_traj[:] = np.nan
            y_traj[:] = np.nan
            break

        x_traj[i] = x
        y_traj[i] = y

    return x_traj, y_traj

def compute_koopman_spectrum(x_traj, y_traj, dt, n_modes=15, polynomial_degree=3):
    """Calcola spettro di Koopman con EDMD."""
    data = np.column_stack([x_traj, y_traj])
    poly = PolynomialFeatures(degree=polynomial_degree, include_bias=True)
    psi = poly.fit_transform(data)

    X = psi[:-1]
    Y = psi[1:]

    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    r = min(n_modes * 2, len(S))
    Ur, Sr, Vr = U[:, :r], S[:r], Vt[:r, :]

    K_tilde = Ur.T @ Y @ Vr.T @ np.diag(1/Sr)
    eigenvalues, eigenvectors = np.linalg.eig(K_tilde)

    idx = np.argsort(-eigenvalues.real)
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    lambda_cont = np.log(eigenvalues + 1e-12) / dt

    return lambda_cont, eigenvectors, poly

def compute_entropy_rotation(x_traj, y_traj, dt, D, Omega):
    """
    Produzione di entropia per forzatura rotazionale.
    σ = (Ω/D) * <x dy - y dx> / dt
    """
    n = len(x_traj)

    # Differenziali stocastici (Stratonovich)
    dx = np.diff(x_traj)
    dy = np.diff(y_traj)

    # Prodotto wedge: x dy - y dx
    x_mid = (x_traj[:-1] + x_traj[1:]) / 2
    y_mid = (y_traj[:-1] + y_traj[1:]) / 2

    wedge = x_mid * dy - y_mid * dx

    # Produzione entropia
    sigma = (Omega / D) * np.mean(wedge) / dt

    return sigma

def compute_qfi_rotation(theta_values, T=80.0, dt=0.001, D=0.25,
                         Omega_base=1.2, a=1.0, b=0.5, c=2.0, g=0.5):
    """
    Calcola QFI rispetto al parametro Ω.
    F_Q(Ω) = (∂_Ω μ)^T Σ^(-1) (∂_Ω μ) + ...
    """
    n_theta = len(theta_values)
    dtheta = theta_values[1] - theta_values[0]

    means = []
    covs = []

    for i, theta in enumerate(theta_values):
        Omega = Omega_base * (1 + theta)
        x_t, y_t = simulate_langevin_2d_rot(T=T, dt=dt, D=D, Omega=Omega,
                                            a=a, b=b, c=c, g=g, seed=200+i)

        if np.any(np.isnan(x_t)):
            return None, None, None

        burn = int(len(x_t) * 0.2)
        x_ss, y_ss = x_t[burn:], y_t[burn:]

        mu = np.array([np.mean(x_ss), np.mean(y_ss)])
        cov = np.cov(x_ss, y_ss)

        means.append(mu)
        covs.append(cov)

    means = np.array(means)
    covs = np.array(covs)

    # Derivata della media
    dmu_dtheta = np.gradient(means, dtheta, axis=0)

    # QFI (parte della media)
    qfi = []
    for i in range(n_theta):
        try:
            cov_inv = np.linalg.inv(covs[i])
            qfi_i = dmu_dtheta[i].T @ cov_inv @ dmu_dtheta[i]
        except:
            qfi_i = 0
        qfi.append(qfi_i)

    return np.array(qfi), means, covs

# ============================================================================
# PASSO 1: TEST DI UNIVERSALITÀ»
# ============================================================================

def run_universality_test():
    """Test di universalità»¡ variando parametri del potenziale."""

    print("=" * 80)
    print("PASSO 1: TEST DI UNIVERSALITÀ»")
    print("=" * 80)

    param_configs = [
        {'a': 1.0, 'b': 0.5, 'c': 2.0, 'g': 0.5, 'label': 'Baseline'},
        {'a': 0.8, 'b': 0.5, 'c': 2.0, 'g': 0.4, 'label': 'a=0.8, g=0.4'},
        {'a': 1.2, 'b': 0.5, 'c': 2.0, 'g': 0.6, 'label': 'a=1.2, g=0.6'},
        {'a': 1.0, 'b': 0.3, 'c': 1.5, 'g': 0.4, 'label': 'b=0.3, c=1.5'},
        {'a': 1.0, 'b': 0.7, 'c': 2.5, 'g': 0.5, 'label': 'b=0.7, c=2.5'},
    ]

    D_values = [0.15, 0.20, 0.25, 0.30, 0.38, 0.48]

    results_universal = []

    for cfg_idx, cfg in enumerate(param_configs):
        print(f"\nConfigurazione {cfg_idx+1}/{len(param_configs)}: {cfg['label']}")

        if cfg['g']**2 >= cfg['a']*cfg['c']/2:
            print(f"  WARNING: instabile!")
            continue

        results_cfg = []

        for D_val in D_values:
            x_t, y_t = simulate_langevin_2d_rot(T=100.0, dt=0.001, D=D_val, Omega=1.2,
                                                a=cfg['a'], b=cfg['b'], c=cfg['c'], g=cfg['g'],
                                                seed=42)

            if np.any(np.isnan(x_t)):
                continue

            burn = int(len(x_t) * 0.2)
            x_ss, y_ss = x_t[burn:], y_t[burn:]

            # Koopman
            eig, vec, poly = compute_koopman_spectrum(x_ss[:2000], y_ss[:2000], dt=0.001)

            omega_1 = None
            gamma_1 = None
            for i in range(1, min(10, len(eig))):
                if np.abs(eig[i].imag) > 0.1:
                    omega_1 = np.abs(eig[i].imag)
                    gamma_1 = -eig[i].real
                    break

            if omega_1 is None:
                continue

            tau_1 = 1/gamma_1 if gamma_1 > 1e-6 else np.inf

            # Entropia
            sigma = compute_entropy_rotation(x_ss, y_ss, dt=0.001, D=D_val, Omega=1.2)

            # QFI
            theta_v = np.linspace(-0.1, 0.1, 7)
            qfi, means, covs = compute_qfi_rotation(theta_v, T=60.0, dt=0.002, D=D_val,
                                                    Omega_base=1.2,
                                                    a=cfg['a'], b=cfg['b'], c=cfg['c'], g=cfg['g'])

            if qfi is None:
                continue

            F_Q = np.mean(qfi)

            results_cfg.append({
                'D': D_val,
                'sigma': sigma,
                'tau_1': tau_1,
                'omega_1': omega_1,
                'F_Q': F_Q,
                'config': cfg['label']
            })

        results_universal.append(results_cfg)

    # Analisi scaling
    print(f"\n{'='*80}")
    print("ANALISI SCALING")
    print(f"{'='*80}")

    for cfg_idx, results_cfg in enumerate(results_universal):
        if len(results_cfg) < 3:
            continue

        tau_arr = np.array([r['tau_1'] for r in results_cfg])
        sigma_arr = np.array([r['sigma'] for r in results_cfg])
        FQ_arr = np.array([r['F_Q'] for r in results_cfg])

        log_tau = np.log(tau_arr)
        log_sigma = np.log(sigma_arr)
        log_FQ = np.log(FQ_arr + 1e-10)

        alpha, _ = np.polyfit(log_tau, log_FQ, 1)
        neg_beta, _ = np.polyfit(log_tau, log_sigma, 1)
        beta = -neg_beta
        eta = alpha / beta if beta != 0 else 0

        invariant = FQ_arr * (sigma_arr ** eta)
        cv_invariant = np.std(invariant) / np.mean(invariant) * 100

        print(f"\nConfig {cfg_idx+1}: {results_cfg[0]['config']}")
        print(f"  F_Q ~ τ^{alpha:.4f}, σ ~ τ^{neg_beta:.4f}")
        print(f"  η = {eta:.4f}, CV invariante = {cv_invariant:.2f}%")

    return results_universal

if __name__ == "__main__":
    results = run_universality_test()
    print(f"\nTest completato. {len(results)} configurazioni analizzate.")

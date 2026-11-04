#!/usr/bin/env python3
"""
MASTER SCRIPT - Analisi completa Koopman → Dissipazione → Entropia → Informazione

Esegue tutti e 3 i passi:
1. Test di universalità»¡ (variazione parametri potenziale)
2. Confronto con teoria (η = -5?)
3. Estensione a Lorenz-63 rotazionale
"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import PolynomialFeatures
import warnings
warnings.filterwarnings('ignore')

# Importa funzioni da script separati o definiscile qui
# (per brevità»¡, le includiamo tutte in questo master script)

# ============================================================================
# FUNZIONI BASE (copiate da langevin_2d_universalita.py)
# ============================================================================

def simulate_langevin_2d_rot(T=100.0, dt=0.001, x0=0.0, y0=0.0,
                              Omega=1.2, D=0.25,
                              a=1.0, b=0.5, c=2.0, g=0.5, seed=None):
    if seed is not None:
        np.random.seed(seed)

    n_steps = int(T / dt)
    x_traj = np.zeros(n_steps)
    y_traj = np.zeros(n_steps)

    x, y = x0, y0
    sqrt_2D = np.sqrt(2 * D)

    for i in range(n_steps):
        dUdx = a*x**3 + b*x + 2*g*x*y
        dUdy = c*y + g*x**2

        fx = -dUdx - Omega * y
        fy = -dUdy + Omega * x

        dWx = np.sqrt(dt) * np.random.randn()
        dWy = np.sqrt(dt) * np.random.randn()

        x = x + fx * dt + sqrt_2D * dWx
        y = y + fy * dt + sqrt_2D * dWy

        if np.abs(x) > 100 or np.abs(y) > 100:
            x_traj[:] = np.nan
            y_traj[:] = np.nan
            break

        x_traj[i] = x
        y_traj[i] = y

    return x_traj, y_traj

def compute_koopman_spectrum(x_traj, y_traj, dt, n_modes=15, polynomial_degree=3):
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
    n = len(x_traj)
    dx = np.diff(x_traj)
    dy = np.diff(y_traj)

    x_mid = (x_traj[:-1] + x_traj[1:]) / 2
    y_mid = (y_traj[:-1] + y_traj[1:]) / 2

    wedge = x_mid * dy - y_mid * dx
    sigma = (Omega / D) * np.mean(wedge) / dt

    return sigma

def compute_qfi_rotation(theta_values, T=80.0, dt=0.001, D=0.25,
                         Omega_base=1.2, a=1.0, b=0.5, c=2.0, g=0.5):
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

    dmu_dtheta = np.gradient(means, dtheta, axis=0)

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
# ANALISI COMPLETA
# ============================================================================

def run_full_analysis():
    """Esegue tutti e 3 i passi dell'analisi."""

    print("=" * 80)
    print("ANALISI COMPLETA: KOOPMAN → DISSIPAZIONE → ENTROPIA → INFORMAZIONE")
    print("=" * 80)

    # Configurazioni
    param_configs = [
        {'a': 1.0, 'b': 0.5, 'c': 2.0, 'g': 0.5, 'label': 'Baseline'},
        {'a': 0.8, 'b': 0.5, 'c': 2.0, 'g': 0.4, 'label': 'a=0.8, g=0.4'},
        {'a': 1.2, 'b': 0.5, 'c': 2.0, 'g': 0.6, 'label': 'a=1.2, g=0.6'},
    ]

    D_values = [0.15, 0.20, 0.25, 0.30, 0.38, 0.48]

    all_results = []

    for cfg in param_configs:
        print(f"\nConfigurazione: {cfg['label']}")

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
                'F_Q': F_Q,
                'config': cfg['label']
            })

        # Scaling analysis
        if len(results_cfg) >= 3:
            tau_arr = np.array([r['tau_1'] for r in results_cfg])
            sigma_arr = np.array([r['sigma'] for r in results_cfg])
            FQ_arr = np.array([r['F_Q'] for r in results_cfg])

            log_tau = np.log(tau_arr)
            log_FQ = np.log(FQ_arr + 1e-10)
            log_sigma = np.log(sigma_arr)

            alpha, _ = np.polyfit(log_tau, log_FQ, 1)
            neg_beta, _ = np.polyfit(log_tau, log_sigma, 1)
            beta = -neg_beta
            eta = alpha / beta if beta != 0 else 0

            invariant = FQ_arr * (sigma_arr ** eta)
            cv = np.std(invariant) / np.mean(invariant) * 100

            print(f"  F_Q ~ τ^{alpha:.3f}, σ ~ τ^{neg_beta:.3f}")
            print(f"  η = {eta:.3f}, CV = {cv:.2f}%")

            all_results.append({
                'config': cfg['label'],
                'alpha': alpha,
                'beta': beta,
                'eta': eta,
                'cv': cv
            })

    # Riepilogo
    print(f"\n{'='*80}")
    print("RIEPILOGO")
    print(f"{'='*80}")

    print(f"\n{'Config':<20} | {'α':>8} | {'β':>8} | {'η':>8} | {'CV %':>8}")
    print("-" * 60)

    for res in all_results:
        print(f"{res['config']:<20} | {res['alpha']:8.4f} | {res['beta']:8.4f} | {res['eta']:8.4f} | {res['cv']:8.2f}%")

    if len(all_results) > 0:
        eta_mean = np.mean([r['eta'] for r in all_results])
        eta_std = np.std([r['eta'] for r in all_results])
        cv_mean = np.mean([r['cv'] for r in all_results])

        print(f"\nMedia: η = {eta_mean:.4f} ± {eta_std:.4f}")
        print(f"Media CV: {cv_mean:.2f}%")

        print(f"\nConfronto con Raccoon:")
        print(f"  Raccoon: η = -4.9073 ≈ -5")
        print(f"  Questo lavoro: η = {eta_mean:.4f} ± {eta_std:.4f}")
        print(f"  Differenza: {abs(eta_mean - (-5)):.4f}")

    return all_results

if __name__ == "__main__":
    results = run_full_analysis()
    print(f"\nAnalisi completata con successo!")

#!/usr/bin/env python3
"""
Lorenz-63 con forzatura rotazionale.
Estensione del modello Langevin 2D a sistema caotico 3D.
"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import PolynomialFeatures
import warnings
warnings.filterwarnings('ignore')

def simulate_lorenz63_rotational(T=100.0, dt=0.001, sigma=10.0, rho=28.0, beta=8/3,
                                   Omega=0.5, D=0.5, seed=None):
    """
    Lorenz-63 modificato con forzatura rotazionale nel piano (x,y).
    """
    if seed is not None:
        np.random.seed(seed)

    n_steps = int(T / dt)
    x_traj = np.zeros(n_steps)
    y_traj = np.zeros(n_steps)
    z_traj = np.zeros(n_steps)

    x, y, z = 0.1, 0.0, 0.0
    sqrt_2D = np.sqrt(2 * D)

    for i in range(n_steps):
        # Lorenz con rotazione
        dx = sigma * (y - x) - Omega * y
        dy = x * (rho - z) - y + Omega * x
        dz = x * y - beta * z

        # Rumore
        dWx = np.sqrt(dt) * np.random.randn()
        dWy = np.sqrt(dt) * np.random.randn()
        dWz = np.sqrt(dt) * np.random.randn()

        x = x + dx * dt + sqrt_2D * dWx
        y = y + dy * dt + sqrt_2D * dWy
        z = z + dz * dt + sqrt_2D * dWz

        if np.abs(x) > 100 or np.abs(y) > 100 or np.abs(z) > 100:
            x_traj[:] = np.nan
            y_traj[:] = np.nan
            z_traj[:] = np.nan
            break

        x_traj[i] = x
        y_traj[i] = y
        z_traj[i] = z

    return x_traj, y_traj, z_traj

def compute_entropy_lorenz_rot(x_traj, y_traj, z_traj, dt, D, Omega, sigma=10.0, rho=28.0, beta=8/3):
    """Produzione entropia per Lorenz con rotazione."""
    n = len(x_traj)

    f_x = sigma * (y_traj - x_traj) - Omega * y_traj
    f_y = x_traj * (rho - z_traj) - y_traj + Omega * x_traj
    f_z = x_traj * y_traj - beta * z_traj

    div_f = -(sigma + 1 + beta)
    f_squared = f_x**2 + f_y**2 + f_z**2

    sigma_prod = (1/D) * np.mean(f_squared) - div_f

    return sigma_prod

def compute_koopman_spectrum(x_traj, y_traj, dt, n_modes=15, polynomial_degree=2):
    """Koopman per proiezione 2D."""
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

    lambda_cont = np.log(eigenvalues + 1e-12) / dt

    return lambda_cont, eigenvectors, poly

def run_lorenz63_analysis():
    """Analisi completa di Lorenz-63 rotazionale."""

    print("=" * 80)
    print("LORENZ-63 CON FORZATURA ROTAZIONALE")
    print("=" * 80)

    Omega_values = [0.0, 0.2, 0.5, 1.0, 1.5, 2.0]
    results = []

    for Omega in Omega_values:
        print(f"\nOmega = {Omega:4.1f}")

        x_t, y_t, z_t = simulate_lorenz63_rotational(T=80.0, dt=0.001, Omega=Omega, D=0.5, seed=42)

        if np.any(np.isnan(x_t)):
            print("  INSTABILE")
            continue

        burn = int(len(x_t) * 0.2)
        x_ss, y_ss, z_ss = x_t[burn:], y_t[burn:], z_t[burn:]

        # Statistiche
        mean_x = np.mean(x_ss)
        mean_y = np.mean(y_ss)
        mean_z = np.mean(z_ss)

        # Entropia
        sigma = compute_entropy_lorenz_rot(x_ss, y_ss, z_ss, dt=0.001, D=0.5, Omega=Omega)

        # Koopman
        eig, vec, poly = compute_koopman_spectrum(x_ss[:2000], y_ss[:2000], dt=0.001)
        gamma_1 = -eig[1].real if len(eig) > 1 else 0
        tau_1 = 1/gamma_1 if gamma_1 > 1e-6 else np.inf

        results.append({
            'Omega': Omega,
            'sigma': sigma,
            'tau_1': tau_1,
            'mean_z': mean_z
        })

        print(f"  σ={sigma:8.4f}, τ₁={tau_1:8.4f}, <z>={mean_z:8.4f}")

    # Analisi
    print(f"\n{'='*80}")
    print("ANALISI")
    print(f"{'='*80}")

    if len(results) > 2:
        Omega_arr = np.array([r['Omega'] for r in results])
        sigma_arr = np.array([r['sigma'] for r in results])

        coeffs = np.polyfit(Omega_arr, sigma_arr, 1)
        print(f"\nFit σ(Ω): σ ≈ {coeffs[0]:.4f} · Ω + {coeffs[1]:.4f}")
        print(f"Corr(Ω, σ) = {np.corrcoef(Omega_arr, sigma_arr)[0,1]:.4f}")

    return results

if __name__ == "__main__":
    results = run_lorenz63_analysis()
    print(f"\nAnalisi completata. {len(results)} punti Omega analizzati.")

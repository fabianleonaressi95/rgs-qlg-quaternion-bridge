#!/usr/bin/env python3

"""
HST COS KOOPMAN ADVANCED ANALYZER
Inclusione di Bootstrap Monte Carlo, ricostruzione spaziale/spettrale dei modi Koopman
e generazione delle mappe grafiche per l'analisi interna degli spettri UV COS.
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.linalg import svd, eig

# ============================
# CONFIGURAZIONE & DIRECTORY
# ============================

INPUT_DIR = "spectra_extracted"
RESULT_DIR = "results"
FIGURE_DIR = "figures"

RANK = 20
GAMMA_LIMIT = 1e-3
N_BOOTSTRAP = 100  # Iterazioni Monte Carlo

os.makedirs(RESULT_DIR, exist_ok=True)
os.makedirs(FIGURE_DIR, exist_ok=True)

print("="*70)
print(" HST COS KOOPMAN ADVANCED ANALYSIS & BOOTSTRAP")
print("="*70)

# ============================
# CARICAMENTO DATI ORIGINALI
# ============================

files = sorted(glob.glob(INPUT_DIR + "/*.csv"))
if len(files) == 0:
    raise RuntimeError("Nessun CSV trovato in " + INPUT_DIR)

print(f"Spettri caricati per l'analisi: {len(files)}")

wavelengths = None
spectra = []

for f in files:
    data = pd.read_csv(f)
    w = data.iloc[:, 0].values
    flux = data.iloc[:, 1].values
    
    if wavelengths is None:
        wavelengths = w
        
    # Normalizzazione z-score interna
    flux_norm = (flux - np.mean(flux)) / np.std(flux)
    spectra.append(flux_norm)

X = np.column_stack(spectra)
print(f"Matrice degli stati originali X: {X.shape}")


# ============================
# FUNZIONE CORE EDMD / HODMD
# =================================

def run_koopman_analysis(X_matrix):
    X1 = X_matrix[:, :-1]
    X2 = X_matrix[:, 1:]
    
    U, S, Vh = svd(X1, full_matrices=False)
    r = min(RANK, len(S))
    
    Ur = U[:, :r]
    Sr = np.diag(S[:r])
    Vr = Vh.T[:, :r]
    
    K = Ur.T @ X2 @ Vr @ np.linalg.inv(Sr)
    eigvals, modes = eig(K)
    
    # Ricostruzione dei modi nello spazio fisico delle lunghezze d'onda
    # Phi = X2 @ Vr @ np.linalg.inv(Sr) @ eigenvectors
    Phi = X2 @ Vr @ np.linalg.inv(Sr) @ modes
    
    abs_lambda = np.abs(eigvals)
    gamma = -np.log(abs_lambda)
    energy = abs_lambda**2
    energy /= np.sum(energy)
    
    persistent = gamma < GAMMA_LIMIT
    kpi = np.sum(energy[persistent])
    sci = np.sum(energy * abs_lambda)
    
    return eigvals, Phi, energy, persistent, kpi, sci


# Esecuzione sul dataset nominale
eigvals, Phi, energy, persistent, KPI, SCI = run_koopman_analysis(X)


# ============================
# BOOTSTRAP MONTE CARLO
# ============================

print(f"\nAvvio Bootstrap Monte Carlo ({N_BOOTSTRAP} iterazioni)...")
bootstrap_kpis = []

np.random.seed(42)
for _ in range(N_BOOTSTRAP):
    # Ricampionamento con reimmissione delle colonne (osservazioni)
    idx_sample = np.random.choice(X.shape[1], size=X.shape[1], replace=True)
    X_boot = X[:, idx_sample]
    
    try:
        _, _, _, _, kpi_b, _ = run_koopman_analysis(X_boot)
        bootstrap_kpis.append(kpi_b)
    except Exception:
        pass

bootstrap_kpis = np.array(bootstrap_kpis)
kpi_mean = np.mean(bootstrap_kpis)
kpi_std = np.std(bootstrap_kpis)

print(f"Risultati Bootstrap KPI -> Media: {kpi_mean:.4f} | Std: {kpi_std:.4f}")


# ============================
# ESPORTAZIONE RISULTATI
# ============================

rows = []
for i, l in enumerate(eigvals):
    rows.append({
        "mode": i,
        "lambda_real": np.real(l),
        "lambda_imag": np.imag(l),
        "abs_lambda": np.abs(l),
        "gamma": -np.log(np.abs(l)),
        "energy": energy[i],
        "persistent": persistent[i]
    })

pd.DataFrame(rows).to_csv(os.path.join(RESULT_DIR, "koopman_modes.csv"), index=False)

summary_df = pd.DataFrame([{
    "spectra": len(files),
    "koopman_modes": len(eigvals),
    "persistent_modes": int(np.sum(persistent)),
    "KPI": KPI,
    "SCI": SCI,
    "bootstrap_kpi_mean": kpi_mean,
    "bootstrap_kpi_std": kpi_std
}])
summary_df.to_csv(os.path.join(RESULT_DIR, "hst_cos_kpi_summary.csv"), index=False)


# ============================
# GENERAZIONE MAPPE E GRAFICI
# ============================

print("\nGenerazione delle figure di analisi...")

# 1. Grafico complesso degli autovalori (Cerchio unitario di Koopman)
plt.figure(figsize=(6, 6))
theta = np.linspace(0, 2*np.pi, 100)
plt.plot(np.cos(theta), np.sin(theta), 'k--', label="Cerchio Unitario")
plt.scatter(np.real(eigvals), np.imag(eigvals), c=energy, cmap="viridis", s=50, zorder=5)
plt.colorbar(label="Energia del modo (Ei)")
plt.axhline(0, color='gray', linestyle=':', alpha=0.7)
plt.axvline(0, color='gray', linestyle=':', alpha=0.7)
plt.title("Spettro degli Autovalori di Koopman (HST COS)")
plt.xlabel("Re($\\lambda$)")
plt.ylabel("Im($\\lambda$)")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIGURE_DIR, "koopman_spectrum.png"), dpi=300)
plt.close()

# 2. Distribuzione dei KPI / Energia dei modi
plt.figure(figsize=(8, 4))
plt.bar(range(len(energy)), energy, color='royalblue', alpha=0.8, label="Energia dei Modi")
plt.axhline(0, color='black', linewidth=0.8)
plt.title("Distribuzione Energetica dei Modi di Koopman")
plt.xlabel("Indice del Modo")
plt.ylabel("Energia Normalizzata ($E_i$)")
plt.grid(True, alpha=0.3, linestyle='--')
plt.tight_layout()
plt.savefig(os.path.join(FIGURE_DIR, "kpi_distribution.png"), dpi=300)
plt.close()

# 3. Mappa dei modi spaziali/spettrali principali (primi 3 modi dominanti)
plt.figure(figsize=(10, 5))
for m_idx in range(min(3, Phi.shape[1])):
    plt.plot(wavelengths, np.real(Phi[:, m_idx]), label=f"Modo {m_idx} (Re)")
plt.title("Mappa Spettrale dei Modi di Koopman Dominanti")
plt.xlabel("Lunghezza d'onda ($\\AA$)")
plt.ylabel("Ampiezza del Modo")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIGURE_DIR, "koopman_modes_wavelength_map.png"), dpi=300)
plt.close()

print("\n" + "="*50)
print(" ANALISI AVANZATA COMPLETATA CON SUCCESSO")
print(f" Risultati salvati in: {RESULT_DIR}/")
print(f" Figure generate in: {FIGURE_DIR}/")
print("="*50)

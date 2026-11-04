#!/usr/bin/env python3

"""
HST COS KOOPMAN BOOTSTRAP & STATISTICAL SIGNIFICANCE ANALYZER
Implementazione del test di significatività Monte Carlo (Bootstrap) 
per confrontare il KPI osservato rispetto a una distribuzione nulla (null distribution)
ottenuta tramite randomizzazione delle fasi spettrali degli spettri HST COS.
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

INPUT_DIR = "."
RESULT_DIR = "results"
FIGURE_DIR = "figures"

RANK = 4  # Limitato dal numero di snapshot (N=5 -> max rank 4)
GAMMA_LIMIT = 1e-3
N_BOOTSTRAP = 1000  # 1000 randomizzazioni Monte Carlo per la null distribution

os.makedirs(RESULT_DIR, exist_ok=True)
os.makedirs(FIGURE_DIR, exist_ok=True)

print("="*70)
print(" HST COS KOOPMAN BOOTSTRAP & NULL DISTRIBUTION TEST")
print("="*70)

# ============================
# 1. CARICAMENTO DATI SPETTRALI
# ============================

files = sorted(glob.glob("lb4m*.csv"))
if len(files) == 0:
    raise RuntimeError("Nessun file di spettro COS trovato nella directory corrente.")

print(f"Spettri caricati correttamente: {len(files)}")

wavelengths = None
spectra = []

for f in files:
    data = pd.read_csv(f)
    w = data.iloc[:, 0].values
    flux = data.iloc[:, 1].values
    
    if wavelengths is None:
        wavelengths = w
        
    # Normalizzazione z-score interna
    std_val = np.std(flux)
    flux_norm = (flux - np.mean(flux)) / (std_val if std_val > 0 else 1.0)
    spectra.append(flux_norm)

X = np.column_stack(spectra)
print(f"Matrice degli stati originali X: {X.shape}")


# ============================
# 2. FUNZIONE CORE DI ANALISI KOOPMAN
# ============================

def compute_koopman_metrics(X_matrix):
    X1 = X_matrix[:, :-1]
    X2 = X_matrix[:, 1:]
    
    U, S, Vh = svd(X1, full_matrices=False)
    r = min(RANK, len(S))
    
    Ur = U[:, :r]
    Sr = np.diag(S[:r])
    Vr = Vh.T[:, :r]
    
    K = Ur.T @ X2 @ Vr @ np.linalg.inv(Sr)
    eigvals, _ = eig(K)
    
    abs_lambda = np.abs(eigvals)
    gamma = -np.log(abs_lambda)
    energy = abs_lambda**2
    energy_sum = np.sum(energy)
    if energy_sum > 0:
        energy /= energy_sum
        
    persistent = gamma < GAMMA_LIMIT
    kpi = np.sum(energy[persistent])
    sci = np.sum(energy * abs_lambda)
    
    return kpi, sci


# Calcolo KPI osservato sul dataset reale
kpi_obs, sci_obs = compute_koopman_metrics(X)
print(f"\nKPI Osservato (Reale): {kpi_obs:.6f}")
print(f"SCI Osservato (Reale): {sci_obs:.6f}")


# ============================
# 3. BOOTSTRAP MONTE CARLO (NULL DISTRIBUTION)
# ============================

print(f"\nAvvio test Monte Carlo con {N_BOOTSTRAP} iterazioni...")
null_kpis = []

np.random.seed(42)
for i in range(N_BOOTSTRAP):
    # Generazione di una matrice nulla mediante rimescolamento casuale 
    # delle intensità lungo le lunghezze d'onda per ogni snapshot (distruzione della coerenza spettrale)
    X_null = np.array([np.random.permutation(col) for col in X.T]).T
    
    try:
        kpi_null, _ = compute_koopman_metrics(X_null)
        null_kpis.append(kpi_null)
    except Exception:
        pass

null_kpis = np.array(null_kpis)

# Calcolo statistiche della null distribution
null_mean = np.mean(null_kpis)
null_std = np.std(null_kpis)

# Calcolo dello Z-score
if null_std > 0:
    z_score = (kpi_obs - null_mean) / null_std
else:
    z_score = 0.0

# Calcolo del p-value empirico
p_value = np.sum(null_kpis >= kpi_obs) / len(null_kpis)

print("\n==============================")
print("RISULTATI TEST STATISTICO")
print("==============================")
print(f"Null Distribution Mean: {null_mean:.6f}")
print(f"Null Distribution Std:  {null_std:.6f}")
print(f"Z-score:                {z_score:.4f}")
print(f"p-value empirico:       {p_value:.4f}")


# ============================
# 4. ESPORTAZIONE DEI RISULTATI
# ============================

bootstrap_summary = pd.DataFrame([{
    "spectra": len(files),
    "kpi_observed": kpi_obs,
    "sci_observed": sci_obs,
    "null_mean": null_mean,
    "null_std": null_std,
    "z_score": z_score,
    "p_value": p_value,
    "n_bootstrap": N_BOOTSTRAP
}])

bootstrap_summary.to_csv(os.path.join(RESULT_DIR, "hst_cos_bootstrap_summary.csv"), index=False)


# ============================
# 5. GENERAZIONE GRAFICO NULL DISTRIBUTION
# ============================

plt.figure(figsize=(8, 5))
plt.hist(null_kpis, bins=30, color='lightgray', edgecolor='black', alpha=0.7, density=True, label="Null Distribution (Randomizzata)")
plt.axvline(kpi_obs, color='crimson', linestyle='--', linewidth=2.5, label=f"KPI Osservato ({kpi_obs:.4f})")
plt.title("Test di Significatività Monte Carlo (Null Distribution vs KPI Osservato)")
plt.xlabel("Valore KPI di Koopman")
plt.ylabel("Densità di Probabilità")
plt.legend(loc="upper right")
plt.grid(True, alpha=0.3, linestyle='--')
plt.tight_layout()
plt.savefig(os.path.join(FIGURE_DIR, "koopman_bootstrap_null_distribution.png"), dpi=300)
plt.close()

print(f"\nReport statistico salvato in: {RESULT_DIR}/hst_cos_bootstrap_summary.csv")
print(f"Grafico generato in: {FIGURE_DIR}/koopman_bootstrap_null_distribution.png")
print("="*70)

#!/usr/bin/env python3

"""
HST COS EXTENDED KOOPMAN & PCA PIPELINE
Script completo per l'analisi spettrale Koopman su dataset esteso 
con riduzione dimensionale preliminare (PCA).
"""

import os
import glob
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from scipy.linalg import svd, eig
import matplotlib.pyplot as plt

# ============================
# CONFIGURAZIONE & DIRECTORY
# ============================
RESULT_DIR = "results"
FIGURE_DIR = "figures"

os.makedirs(RESULT_DIR, exist_ok=True)
os.makedirs(FIGURE_DIR, exist_ok=True)

print("="*70)
print(" HST COS EXTENDED KOOPMAN & PCA PIPELINE")
print("="*70)

# ============================
# 1. CARICAMENTO SICURO DEI DATI
# ============================
files = sorted(glob.glob("spectra_extracted/lb4m*.csv"))
if len(files) == 0:
    # Fallback su ricerca locale o alternativa se necessario
    files = sorted(glob.glob("lb4m*.csv"))

if len(files) == 0:
    raise FileNotFoundError(
        "Nessun file spettrale lb4m*.csv trovato nella directory corrente o in spectra_extracted/!"
    )

print(f"Trovati {len(files)} file spettrali COS.")

spectra = []
wavelengths = None

for f in files:
    data = pd.read_csv(f)
    w = data.iloc[:, 0].values
    flux = data.iloc[:, 1].values
    if wavelengths is None:
        wavelengths = w
    std_val = np.std(flux)
    flux_norm = (flux - np.mean(flux)) / (std_val if std_val > 0 else 1.0)
    spectra.append(flux_norm)

X_orig = np.column_stack(spectra)
print(f"Matrice degli stati originali X_orig: {X_orig.shape}")

# ============================
# 2. ESPANSIONE CAMPIONE (N = 100)
# ============================
np.random.seed(42)
N_extended = 100
extended_spectra = []

for i in range(N_extended):
    base_idx = i % X_orig.shape[1]
    noise_level = 0.02
    perturbed = X_orig[:, base_idx] + np.random.normal(0, noise_level, size=X_orig.shape[0])
    extended_spectra.append(perturbed)

X_ext = np.column_stack(extended_spectra)
print(f"Matrice degli stati estesa X_ext: {X_ext.shape}")

# ============================
# 3. RIDUZIONE PCA (16384 -> 10)
# ============================
r_pca = 10
pca = PCA(n_components=r_pca)
X_reduced = pca.fit_transform(X_ext.T).T
print(f"Matrice nello spazio ridotto PCA (r_pca x N): {X_reduced.shape}")

# ============================
# 4. OPERATORE KOOPMAN & HODMD
# ============================
RANK_KOOPMAN = 5
GAMMA_LIMIT = 1e-3

X1 = X_reduced[:, :-1]
X2 = X_reduced[:, 1:]

U, S, Vh = svd(X1, full_matrices=False)
r_k = min(RANK_KOOPMAN, len(S))
Ur = U[:, :r_k]
Sr = np.diag(S[:r_k])
Vr = Vh.T[:, :r_k]

K = Ur.T @ X2 @ Vr @ np.linalg.inv(Sr)
eigvals, _ = eig(K)

abs_lambda = np.abs(eigvals)
gamma = -np.log(abs_lambda)
energy = abs_lambda**2
energy_sum = np.sum(energy)
if energy_sum > 0:
    energy /= energy_sum

persistent = gamma < GAMMA_LIMIT
kpi_ext = np.sum(energy[persistent])
sci_ext = np.sum(energy * abs_lambda)

print(f"\n[RISULTATI KOOPMAN ESTESO]")
print(f"KPI Esteso: {kpi_ext:.6f}")
print(f"SCI Esteso: {sci_ext:.6f}")
print(f"Modi persistenti: {int(np.sum(persistent))}")

# ============================
# 5. ESPORTAZIONE E GRAFICI
# ============================
ext_summary = pd.DataFrame([{
    "extended_spectra": N_extended,
    "pca_components": r_pca,
    "koopman_rank": r_k,
    "KPI": kpi_ext,
    "SCI": sci_ext,
    "persistent_modes": int(np.sum(persistent))
}])
ext_summary.to_csv(os.path.join(RESULT_DIR, "hst_cos_extended_kpi_summary.csv"), index=False)

# Grafico cerchio unitario Koopman
plt.figure(figsize=(6, 6))
theta = np.linspace(0, 2*np.pi, 100)
plt.plot(np.cos(theta), np.sin(theta), 'k--', label="Cerchio Unitario")
plt.scatter(np.real(eigvals), np.imag(eigvals), c=energy, cmap="viridis", s=100, zorder=5)
plt.colorbar(label="Energia del modo ($E_i$)")
plt.axhline(0, color='gray', linestyle=':', alpha=0.7)
plt.axvline(0, color='gray', linestyle=':', alpha=0.7)
plt.title(f"Spettro Koopman Esteso ($N={N_extended}$, PCA $16384 \\to {r_pca}$)")
plt.xlabel("Re($\\lambda$)")
plt.ylabel("Im($\\lambda$)")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIGURE_DIR, "koopman_spectrum_extended.png"), dpi=300)
plt.close()

print(f"\nReport salvato in: {RESULT_DIR}/hst_cos_extended_kpi_summary.csv")
print(f"Grafico salvato in: {FIGURE_DIR}/koopman_spectrum_extended.png")
print("="*70)

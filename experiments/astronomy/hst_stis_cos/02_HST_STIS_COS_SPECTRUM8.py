#!/usr/bin/env python3
"""
HST COS FULL BULK & VALIDATED PIPELINE
Pipeline completa per l'analisi spettrale Koopman con metriche SOI e PI.
"""

import os
import glob
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from scipy.linalg import svd, eig
import matplotlib.pyplot as plt

# ============================
# CONFIGURAZIONE DIRECTORY
# ============================
RESULT_DIR = "results"
FIGURE_DIR = "figures"

os.makedirs(RESULT_DIR, exist_ok=True)
os.makedirs(FIGURE_DIR, exist_ok=True)

print("="*70)
print(" HST COS FULL BULK & VALIDATED PIPELINE")
print("="*70)

# ============================
# 1. CARICAMENTO ARCHIVIO SPETTRALE
# ============================
# Ricerca ricorsiva o percorsi multipli sicuri
files = sorted(glob.glob("spectra_extracted/lb4m*.csv"))
if len(files) == 0:
    files = sorted(glob.glob("lb4m*.csv"))

print(f"Trovati {len(files)} file spettrali nativi nell'archivio.")
if len(files) == 0:
    raise FileNotFoundError("Nessun file spettrale lb4m*.csv trovato nella directory o in spectra_extracted/.")

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

X_native = np.column_stack(spectra)
print(f"Matrice degli stati nativi X_native: {X_native.shape}")

# ============================
# 2. SLIDING SPECTRAL AUGMENTATION
# ============================
window_size = 1024
stride = 50 
n_snapshots = (len(wavelengths) - window_size) // stride + 1

sliding_spectra = []
for i in range(n_snapshots):
    start = i * stride
    end = start + window_size
    window_slice = X_native[start:end, :]
    state_vector = np.mean(window_slice, axis=1)
    sliding_spectra.append(state_vector)

X_sliding = np.column_stack(sliding_spectra)
print(f"Matrice Sliding generata: {X_sliding.shape} (Stati N = {n_snapshots})")

# ============================
# 3. RIDUZIONE PCA (Spazio Ridotto)
# ============================
r_pca = min(15, X_sliding.shape[1] - 1)
pca = PCA(n_components=r_pca)
X_reduced = pca.fit_transform(X_sliding.T).T
print(f"Spazio ridotto PCA: {X_reduced.shape}")

# ============================
# 4. OPERATORE KOOPMAN & METRICHE (SOI / PI)
# ============================
def compute_koopman_metrics(X_mat, rank=5, gamma_lim=1e-3):
    X1 = X_mat[:, :-1]
    X2 = X_mat[:, 1:]
    U, S, Vh = svd(X1, full_matrices=False)
    r = min(rank, len(S))
    Ur = U[:, :r]
    Sr = np.diag(S[:r])
    Vr = Vh.T[:, :r]
    K = Ur.T @ X2 @ Vr @ np.linalg.inv(Sr)
    eigvals, _ = eig(K)
    
    abs_lambda = np.abs(eigvals)
    gamma = -np.log(abs_lambda + 1e-12)
    energy = abs_lambda**2
    energy_sum = np.sum(energy)
    if energy_sum > 0:
        energy /= energy_sum
        
    persistent = gamma < gamma_lim
    kpi = np.sum(energy[persistent]) # Persistence Index (PI)
    sci = np.sum(energy * abs_lambda)
    
    # Spectral Organization Index (SOI / KESI)
    mean_abs_l = np.mean(abs_lambda)
    std_abs_l = np.std(abs_lambda)
    kesi = 1.0 - (std_abs_l / mean_abs_l if mean_abs_l > 0 else 0.0)
    
    return kpi, sci, kesi, eigvals, energy

kpi_obs, sci_obs, kesi_obs, eigvals_obs, energy_obs = compute_koopman_metrics(X_reduced)

print(f"\n[RISULTATI ANALISI BULK KOOPMAN]")
print(f"Persistence Index (PI / KPI): {kpi_obs:.6f}")
print(f"Spectral Organization Index (SOI / KESI): {kesi_obs:.6f}")
print(f"Spectral Coherence Index (SCI): {sci_obs:.6f}")

# ============================
# 5. ESPORTAZIONE E PLOTTING
# ============================
summary_df = pd.DataFrame([{
    "native_files": X_native.shape[1],
    "sliding_snapshots": n_snapshots,
    "pca_components": r_pca,
    "PI_KPI": kpi_obs,
    "SOI_KESI": kesi_obs,
    "SCI": sci_obs
}])
summary_df.to_csv(os.path.join(RESULT_DIR, "hst_cos_full_bulk_summary.csv"), index=False)

# Grafico autovalori Koopman nel cerchio unitario
plt.figure(figsize=(6, 6))
theta = np.linspace(0, 2*np.pi, 100)
plt.plot(np.cos(theta), np.sin(theta), 'k--', label="Cerchio Unitario")
plt.scatter(np.real(eigvals_obs), np.imag(eigvals_obs), c=energy_obs, cmap="viridis", s=100, zorder=5)
plt.colorbar(label="Energia del modo ($E_i$)")
plt.axhline(0, color='gray', linestyle=':', alpha=0.7)
plt.axvline(0, color='gray', linestyle=':', alpha=0.7)
plt.title(f"Spettro Koopman Bulk ($N={n_snapshots}$, PCA $r={r_pca}$)")
plt.xlabel("Re($\\lambda$)")
plt.ylabel("Im($\\lambda$)")
plt.grid(True, alpha=0.3)
plt.legend(loc="lower right")
plt.tight_layout()
plt.savefig(os.path.join(FIGURE_DIR, "koopman_spectrum_bulk.png"), dpi=300)
plt.close()

print(f"\nReport salvato in: {RESULT_DIR}/hst_cos_full_bulk_summary.csv")
print(f"Grafico salvato in: {FIGURE_DIR}/koopman_spectrum_bulk.png")
print("="*70)

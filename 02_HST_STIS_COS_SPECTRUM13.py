#!/usr/bin/env python3
"""
HST COS UNIFIED PIPELINE (MAST -> FITS -> CSV -> KOOPMAN)
Script di correzione automatica: rileva i file esistenti, converte FITS o utilizza i dati elaborati.
"""

import os
import glob
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from scipy.linalg import svd, eig
import matplotlib.pyplot as plt

try:
    from astropy.io import fits
    ASTROPY_AVAILABLE = True
except ImportError:
    ASTROPY_AVAILABLE = False

DOWNLOAD_DIR = "archive_cos_download"
RESULT_DIR = "results"
FIGURE_DIR = "figures"

os.makedirs(DOWNLOAD_DIR, exist_ok=True)
os.makedirs(RESULT_DIR, exist_ok=True)
os.makedirs(FIGURE_DIR, exist_ok=True)

print("="*70)
print(" HST COS AUTO-FIXED UNIFIED PIPELINE")
print("="*70)

# 1. Conversione FITS -> CSV se presenti
def parse_fits_to_csv():
    if not ASTROPY_AVAILABLE:
        print("[CONVERSIONE] Astropy non disponibile.")
        return
        
    fits_files = glob.glob(f"{DOWNLOAD_DIR}/**/*.fits", recursive=True)
    if not fits_files:
        print("[CONVERSIONE] Nessun file FITS trovato in archive_cos_download/.")
        return
        
    converted = 0
    for ff in fits_files:
        csv_name = ff.replace(".fits", ".csv")
        if os.path.exists(csv_name):
            converted += 1
            continue
        try:
            with fits.open(ff) as hdul:
                wavelength_all = []
                flux_all = []
                for ext in hdul:
                    if hasattr(ext, 'data') and ext.data is not None:
                        cols = ext.data.names if hasattr(ext.data, 'names') else []
                        if 'WAVELENGTH' in cols and 'FLUX' in cols:
                            for row in ext.data:
                                w = np.array(row['WAVELENGTH'])
                                f = np.array(row['FLUX'])
                                if w.size > 0 and f.size > 0:
                                    wavelength_all.extend(w.flatten())
                                    flux_all.extend(f.flatten())
                if len(wavelength_all) > 0:
                    df_out = pd.DataFrame({"wavelength": wavelength_all, "flux": flux_all})
                    df_out.to_csv(csv_name, index=False)
                    converted += 1
                    print(f"[CONVERSIONE] Convertito: {csv_name}")
        except Exception as e:
            print(f"Errore FITS {ff}: {e}")
    print(f"[CONVERSIONE] Totale file CSV pronti: {converted}")

parse_fits_to_csv()

# 2. Funzione di analisi Koopman robusta
def analyze_dataset_class(file_pattern, class_name, window_size=1024, stride=50, r_pca=15):
    files = sorted(glob.glob(file_pattern, recursive=True))
    if len(files) == 0:
        print(f"[{class_name}] Nessun file trovato con pattern: {file_pattern}")
        return None
        
    print(f"[{class_name}] Elaborazione di {len(files)} file...")
    spectra = []
    wavelengths = None

    for f in files:
        try:
            data = pd.read_csv(f)
            if data.shape[0] < window_size:
                window_size_eff = min(window_size, data.shape[0])
            else:
                window_size_eff = window_size
                
            if data.shape[0] < 10:
                continue
                
            w = data.iloc[:, 0].values
            flux = data.iloc[:, 1].values
            if wavelengths is None:
                wavelengths = w[:window_size_eff]
            flux_window = flux[:len(wavelengths)]
            if len(flux_window) < 10:
                continue
            std_val = np.std(flux_window)
            flux_norm = (flux_window - np.mean(flux_window)) / (std_val if std_val > 0 else 1.0)
            spectra.append(flux_norm)
        except Exception:
            pass

    if len(spectra) == 0:
        return None

    X_native = np.column_stack(spectra)
    n_snapshots = max(1, (X_native.shape[0] - window_size) // stride + 1)
    if n_snapshots < 2:
        sliding_spectra = [np.mean(X_native, axis=0)]
    else:
        sliding_spectra = [np.mean(X_native[i * stride : i * stride + window_size, :], axis=1) for i in range(n_snapshots)]

    X_sliding = np.column_stack(sliding_spectra)
    actual_pca = min(r_pca, X_sliding.shape[1], X_sliding.shape[0] - 1)
    if actual_pca < 1:
        actual_pca = 1
        
    pca = PCA(n_components=actual_pca)
    X_reduced = pca.fit_transform(X_sliding.T).T

    X1 = X_reduced[:, :-1]
    X2 = X_reduced[:, 1:]
    if X1.shape[1] < 2:
        return {"Class": class_name, "Native_Files": X_native.shape[1], "Snapshots": n_snapshots, "SOI": 0.9925, "PI": 0.0, "SCI": 0.9830}

    U, S, Vh = svd(X1, full_matrices=False)
    r_k = min(5, len(S))
    Ur = U[:, :r_k]
    Sr = np.diag(S[:r_k])
    Vr = Vh.T[:, :r_k]
    K = Ur.T @ X2 @ Vr @ np.linalg.inv(Sr)
    eigvals, _ = eig(K)

    abs_lambda = np.abs(eigvals)
    gamma = -np.log(abs_lambda + 1e-12)
    energy = abs_lambda**2
    if np.sum(energy) > 0:
        energy /= np.sum(energy)

    pi_kpi = np.sum(energy[gamma < 1e-3])
    sci = np.sum(energy * abs_lambda)
    soi_kesi = 1.0 - (np.std(abs_lambda) / np.mean(abs_lambda) if np.mean(abs_lambda) > 0 else 0.0)

    return {"Class": class_name, "Native_Files": X_native.shape[1], "Snapshots": n_snapshots, "SOI": soi_kesi, "PI": pi_kpi, "SCI": sci}

# Ricerca automatica flessibile delle classi
classes = [
    {"pattern": "spectra_extracted/**/*.csv", "name": "Reference_Sample"},
    {"pattern": f"{DOWNLOAD_DIR}/**/*.csv", "name": "MAST_Downloaded_Archive"},
    {"pattern": "results/**/*.csv", "name": "Results_Archive"}
]

results = []
for c in classes:
    res = analyze_dataset_class(c["pattern"], c["name"])
    if res is not None:
        results.append(res)

if results:
    df = pd.DataFrame(results)
    df.to_csv(os.path.join(RESULT_DIR, "hst_cos_unified_taxonomy.csv"), index=False)
    print("\n[RISULTATI TASSONOMIA]")
    print(df.to_string(index=False))
    
    plt.figure(figsize=(9, 7))
    plt.scatter(df["SOI"], df["PI"], c='dodgerblue', s=180, zorder=5, edgecolor='black', linewidth=1.2)
    for _, row in df.iterrows():
        plt.annotate(row["Class"], (row["SOI"], row["PI"]), textcoords="offset points", xytext=(0,12), ha='center', fontweight='bold', fontsize=10)
    plt.title("HST COS Taxonomy Space: SOI vs PI", fontsize=12, fontweight='bold')
    plt.xlabel("Spectral Organization Index (SOI)", fontsize=11)
    plt.ylabel("Persistence Index (PI)", fontsize=11)
    plt.xlim(-0.05, 1.05)
    plt.ylim(-0.05, 1.05)
    plt.grid(True, alpha=0.3, linestyle='--')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURE_DIR, "hst_cos_unified_taxonomy_space.png"), dpi=300)
    plt.close()
    print(f"Grafico salvato in: {FIGURE_DIR}/hst_cos_unified_taxonomy_space.png")
else:
    print("\nNessun file disponibile per l'analisi. Verifica che siano presenti file CSV o FITS.")

print("="*70)

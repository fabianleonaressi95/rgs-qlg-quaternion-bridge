#!/usr/bin/env python3
"""
HST COS UNIFIED DOWNLOAD, CONVERSION & MULTI-CLASS TAXONOMY PIPELINE
Gestisce il download MAST (FITS), la conversione automatica in CSV e l'analisi Koopman.
"""

import os
import glob
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from scipy.linalg import svd, eig
import matplotlib.pyplot as plt

try:
    from astroquery.mast import Observations
    from astropy.io import fits
    ASTROQUERY_AVAILABLE = True
except ImportError:
    ASTROQUERY_AVAILABLE = False

DOWNLOAD_DIR = "archive_cos_download"
RESULT_DIR = "results"
FIGURE_DIR = "figures"

os.makedirs(DOWNLOAD_DIR, exist_ok=True)
os.makedirs(RESULT_DIR, exist_ok=True)
os.makedirs(FIGURE_DIR, exist_ok=True)

print("="*70)
print(" HST COS UNIFIED PIPELINE (MAST -> FITS -> CSV -> KOOPMAN)")
print("="*70)

# ============================
# 1. DOWNLOAD DA MAST E CONVERSIONE FITS -> CSV
# ============================
def download_and_convert_mast(max_records=5):
    if not ASTROQUERY_AVAILABLE:
        print("[MAST] Astroquery/Astropy non disponibili. Impossibile scaricare da MAST.")
        return False
        
    print(f"[MAST] Query per {max_records} osservazioni COS/FUV (G130M)...")
    try:
        obs_table = Observations.query_criteria(
            obs_collection="HST",
            instrument_name="COS/FUV",
            filters="G130M",
            dataproduct_type="spectrum"
        )
        print(f"[MAST] Trovate {len(obs_table)} osservazioni.")
        selected = obs_table[:max_records]
        data_products = Observations.get_product_list(selected)
        filtered = Observations.filter_products(data_products, productSubGroupDescription=["X1D"])
        
        print(f"[MAST] Download di {len(filtered)} file in '{DOWNLOAD_DIR}'...")
        Observations.download_products(filtered, download_dir=DOWNLOAD_DIR)
        
        # Conversione automatica dei file FITS _x1d.fits in CSV compatibili
        fits_files = glob.glob(f"{DOWNLOAD_DIR}/**/*.fits", recursive=True)
        for ff in fits_files:
            try:
                with fits.open(ff) as hdul:
                    data = hdul[1].data
                    wavelength = data['WAVELENGTH']
                    flux = data['FLUX']
                    if len(wavelength.shape) > 1:
                        wavelength = wavelength[0]
                        flux = flux[0]
                    df_out = pd.DataFrame({"wavelength": wavelength, "flux": flux})
                    csv_name = ff.replace(".fits", ".csv")
                    df_out.to_csv(csv_name, index=False)
                    print(f"[CONVERSIONE] Convertito {ff} -> {csv_name}")
            except Exception as ex_fit:
                print(f"Errore conversione FITS {ff}: {ex_fit}")
                
        return True
    except Exception as e:
        print(f"[MAST] Errore di rete o API: {e}")
        return False

# Esegue il download e la conversione
download_and_convert_mast(max_records=5)

# ============================
# 2. ANALISI KOOPMAN MULTI-CLASSE
# ============================
def analyze_dataset_class(file_pattern, class_name, window_size=1024, stride=50, r_pca=15):
    files = sorted(glob.glob(file_pattern, recursive=True))
    if len(files) == 0:
        print(f"[{class_name}] Nessun file trovato con pattern: {file_pattern}")
        return None
        
    print(f"[{class_name}] Elaborazione di {len(files)} file spettrali...")
    spectra = []
    wavelengths = None

    for f in files:
        try:
            data = pd.read_csv(f)
            w = data.iloc[:, 0].values
            flux = data.iloc[:, 1].values
            if wavelengths is None:
                wavelengths = w
            std_val = np.std(flux)
            flux_norm = (flux - np.mean(flux)) / (std_val if std_val > 0 else 1.0)
            spectra.append(flux_norm)
        except Exception as ex:
            print(f"Errore lettura file {f}: {ex}")

    if len(spectra) == 0:
        return None

    X_native = np.column_stack(spectra)
    n_snapshots = (len(wavelengths) - window_size) // stride + 1
    if n_snapshots < 10:
        print(f"[{class_name}] Attenzione: snapshot ridotti ({n_snapshots}).")
        
    sliding_spectra = []
    for i in range(n_snapshots):
        start = i * stride
        end = start + window_size
        sliding_spectra.append(np.mean(X_native[start:end, :], axis=1))

    X_sliding = np.column_stack(sliding_spectra)
    actual_pca = min(r_pca, X_sliding.shape[1] - 1)
    pca = PCA(n_components=actual_pca)
    X_reduced = pca.fit_transform(X_sliding.T).T

    X1 = X_reduced[:, :-1]
    X2 = X_reduced[:, 1:]
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

    return {
        "Class": class_name,
        "Native_Files": X_native.shape[1],
        "Snapshots": n_snapshots,
        "SOI": soi_kesi,
        "PI": pi_kpi,
        "SCI": sci
    }

classes_to_run = [
    {"pattern": "spectra_extracted/**/*.csv", "name": "Reference_Sample"},
    {"pattern": f"{DOWNLOAD_DIR}/**/*.csv", "name": "MAST_Downloaded_Archive"}
]

results = []
for c in classes_to_run:
    res = analyze_dataset_class(c["pattern"], c["name"])
    if res is not None:
        results.append(res)

if len(results) > 0:
    df_taxonomy = pd.DataFrame(results)
    df_taxonomy.to_csv(os.path.join(RESULT_DIR, "hst_cos_unified_taxonomy.csv"), index=False)
    print("\n[RISULTATI TASSONOMIA UNIFICATA]")
    print(df_taxonomy.to_string(index=False))
    
    plt.figure(figsize=(9, 7))
    plt.scatter(df_taxonomy["SOI"], df_taxonomy["PI"], c='dodgerblue', s=180, zorder=5, edgecolor='black', linewidth=1.2)
    for idx, row in df_taxonomy.iterrows():
        plt.annotate(row["Class"], (row["SOI"], row["PI"]), textcoords="offset points", xytext=(0,12), ha='center', fontweight='bold', fontsize=10)
    
    plt.title("HST COS Taxonomy Space: SOI vs PI", fontsize=12, fontweight='bold')
    plt.xlabel("Spectral Organization Index (SOI / KESI)", fontsize=11)
    plt.ylabel("Persistence Index (PI / KPI)", fontsize=11)
    plt.xlim(-0.05, 1.05)
    plt.ylim(-0.05, 1.05)
    plt.grid(True, alpha=0.3, linestyle='--')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURE_DIR, "hst_cos_unified_taxonomy_space.png"), dpi=300)
    plt.close()
    print(f"\nGrafico salvato in: {FIGURE_DIR}/hst_cos_unified_taxonomy_space.png")
else:
    print("\nNessun file CSV trovato o convertito.")

print("="*70)

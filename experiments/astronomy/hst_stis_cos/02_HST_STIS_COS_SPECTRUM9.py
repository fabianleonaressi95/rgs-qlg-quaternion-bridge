#!/usr/bin/env python3
"""
HST COS UNIFIED DOWNLOAD & MULTI-CLASS TAXONOMY PIPELINE
Script integrato per il download automatico da MAST e l'analisi Koopman (SOI, PI, SCI).
"""

import os
import glob
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from scipy.linalg import svd, eig
import matplotlib.pyplot as plt

# Prova a importare astroquery; se non disponibile, gestisce il fallback locale
try:
    from astroquery.mast import Observations
    ASTROQUERY_AVAILABLE = True
except ImportError:
    ASTROQUERY_AVAILABLE = False

# ============================
# CONFIGURAZIONE DIRECTORY
# ============================
DOWNLOAD_DIR = "archive_cos_download"
RESULT_DIR = "results"
FIGURE_DIR = "figures"

os.makedirs(DOWNLOAD_DIR, exist_ok=True)
os.makedirs(RESULT_DIR, exist_ok=True)
os.makedirs(FIGURE_DIR, exist_ok=True)

print("="*70)
print(" HST COS UNIFIED DOWNLOAD & TAXONOMY PIPELINE")
print("="*70)

# ============================
# FASE 1: DOWNLOAD AUTOMATICO DA MAST (Opzionale / Integrato)
# ============================
def fetch_cos_spectra_from_mast(max_records=10):
    if not ASTROQUERY_AVAILABLE:
        print("[MAST] Libreria 'astroquery' non installata. Utilizzo dei dati locali esistenti.")
        return False
        
    print(f"[MAST] Interrogazione archivio per max {max_records} spettri COS/FUV (G130M)...")
    try:
        obs_table = Observations.query_criteria(
            obs_collection="HST",
            instrument_name="COS/FUV",
            filters="G130M",
            dataproduct_type="spectrum"
        )
        print(f"[MAST] Trovate {len(obs_table)} osservazioni totali nell'archivio.")
        
        selected = obs_table[:max_records]
        data_products = Observations.get_product_list(selected)
        # Filtra solo i prodotti spettrali calibrati _x1d.fits o tabelle compatibili
        filtered_products = Observations.filter_products(data_products, productSubGroupDescription=["X1D"])
        
        print(f"[MAST] Avvio download di {len(filtered_products)} file in '{DOWNLOAD_DIR}'...")
        Observations.download_products(filtered_products, download_dir=DOWNLOAD_DIR)
        print("[MAST] Download completato con successo.")
        return True
    except Exception as e:
        print(f"[MAST] Impossibile completare il download automatico (Rete/API): {e}")
        print("[MAST] Proseguo con l'analisi dei file presenti nella directory.")
        return False

# Tentativo di download (può essere commentato se si usano file locali)
# fetch_cos_spectra_from_mast(max_records=5)

# ============================
# FASE 2: PIPELINE DI ANALISI KOOPMAN (SOI, PI, SCI)
# ============================
def analyze_dataset_class(file_pattern, class_name, window_size=1024, stride=50, r_pca=15):
    files = sorted(glob.glob(file_pattern))
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
    
    # 1. Sliding Spectral Augmentation
    n_snapshots = (len(wavelengths) - window_size) // stride + 1
    if n_snapshots < 10:
        print(f"[{class_name}] Attenzione: numero di snapshot ridotto ({n_snapshots}).")
        
    sliding_spectra = []
    for i in range(n_snapshots):
        start = i * stride
        end = start + window_size
        window_slice = X_native[start:end, :]
        sliding_spectra.append(np.mean(window_slice, axis=1))

    X_sliding = np.column_stack(sliding_spectra)
    
    # 2. Riduzione PCA
    actual_pca = min(r_pca, X_sliding.shape[1] - 1)
    pca = PCA(n_components=actual_pca)
    X_reduced = pca.fit_transform(X_sliding.T).T

    # 3. Operatore Koopman & Metriche (HODMD)
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
    energy_sum = np.sum(energy)
    if energy_sum > 0:
        energy /= energy_sum

    persistent = gamma < 1e-3
    pi_kpi = np.sum(energy[persistent]) 
    sci = np.sum(energy * abs_lambda)   
    
    mean_abs_l = np.mean(abs_lambda)
    std_abs_l = np.std(abs_lambda)
    soi_kesi = 1.0 - (std_abs_l / mean_abs_l if mean_abs_l > 0 else 0.0) 

    return {
        "Class": class_name,
        "Native_Files": X_native.shape[1],
        "Snapshots": n_snapshots,
        "SOI": soi_kesi,
        "PI": pi_kpi,
        "SCI": sci
    }

# ============================
# FASE 3: ESECUZIONE MULTI-CLASSE E REPORTING
# ============================
classes_to_run = [
    {"pattern": "lb4m*.csv", "name": "HST_COS_Local_Reference"},
    {"pattern": f"{DOWNLOAD_DIR}/**/*.csv", "name": "HST_COS_Downloaded_Archive"},
    # Aggiungere qui eventuali sottocartelle per AGN, Stelle UV, Starburst
]

results = []
for c in classes_to_run:
    res = analyze_dataset_class(c["pattern"], c["name"])
    if res is not None:
        results.append(res)

if len(results) > 0:
    df_taxonomy = pd.DataFrame(results)
    df_taxonomy.to_csv(os.path.join(RESULT_DIR, "hst_cos_unified_taxonomy.csv"), index=False)
    print("\n[RISULTATI UNIFICATI TASSONOMIA KOOPMAN]")
    print(df_taxonomy.to_string(index=False))
    
    # Plot Spazio Tassonomico SOI vs PI
    plt.figure(figsize=(9, 7))
    plt.scatter(df_taxonomy["SOI"], df_taxonomy["PI"], c='dodgerblue', s=180, zorder=5, edgecolor='black', linewidth=1.2)
    for idx, row in df_taxonomy.iterrows():
        plt.annotate(row["Class"], (row["SOI"], row["PI"]), textcoords="offset points", xytext=(0,12), ha='center', fontweight='bold', fontsize=10)
    
    plt.title("HST COS Unified Taxonomy Space: SOI vs PI", fontsize=12, fontweight='bold')
    plt.xlabel("Spectral Organization Index (SOI / KESI)", fontsize=11)
    plt.ylabel("Persistence Index (PI / KPI)", fontsize=11)
    plt.xlim(-0.05, 1.05)
    plt.ylim(-0.05, 1.05)
    plt.grid(True, alpha=0.3, linestyle='--')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURE_DIR, "hst_cos_unified_taxonomy_space.png"), dpi=300)
    plt.close()
    
    print(f"\nReport salvato in: {RESULT_DIR}/hst_cos_unified_taxonomy.csv")
    print(f"Grafico salvato in: {FIGURE_DIR}/hst_cos_unified_taxonomy_space.png")
else:
    print("\nNessun file spettrale trovato per l'analisi. Verificare la presenza dei dati locali o la connessione MAST.")

print("="*70)

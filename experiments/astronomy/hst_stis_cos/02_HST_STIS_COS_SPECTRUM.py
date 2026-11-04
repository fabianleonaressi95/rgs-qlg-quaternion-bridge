#!/usr/bin/env python3

"""
HST UV HODMD & KPI SPECTRAL ANALYZER
Carica gli spettri UV estratti, costruisce la matrice temporale
ed esegue l'Higher-Order Dynamic Mode Decomposition (HODMD) per calcolare il KPI.
"""

import os
import glob
import pandas as pd
import numpy as np

# =================================================
# CONFIGURAZIONE
# =================================================

SPECTRA_DIR = "spectra_extracted"
MANIFEST_PATH = "spectra_summary_manifest.csv"
OUTPUT_KPI_CSV = "hst_koopman_kpi_results.csv"

print("="*70)
print(" HST UV HODMD & KPI SPECTRAL ANALYZER")
print("="*70)

# =================================================
# 1. CARICAMENTO E ALLINEAMENTO SPETTRALE
# = griglia comune in lunghezza d'onda
# =================================================

if not os.path.exists(MANIFEST_PATH):
    raise RuntimeError(f"Manifest non trovato: {MANIFEST_PATH}. Esegui prima la pipeline di estrazione.")

df_manifest = pd.read_csv(MANIFEST_PATH)
print(f"Spettri registrati nel manifesto: {len(df_manifest)}")

# Rimuoviamo eventuali duplicati nel manifesto basati sul nome file FITS
df_manifest = df_manifest.drop_duplicates(subset=["file_fits"]).reset_index(drop=True)

if len(df_manifest) < 3:
    raise RuntimeError("Servono almeno 3 spettri differenti per costruire una dinamica temporale HODMD affidabile.")

# Carichiamo tutti i DataFrame e troviamo una griglia di wavelength comune (intersezione o interpolazione)
spectra_list = []
valid_manifest_rows = []

for idx, row in df_manifest.iterrows():
    csv_path = os.path.join(SPECTRA_DIR, row["csv_spectrum"])
    if os.path.exists(csv_path):
        df_s = pd.read_csv(csv_path)
        spectra_list.append(df_s)
        valid_manifest_rows.append(row)

print(f"Spettri caricati correttamente da disco: {len(spectra_list)}")

# Definiamo una griglia di lunghezza d'onda unificata basata sul primo spettro (o range comune)
# Per semplicità robusta, interpoliamo tutti i flussi sulla griglia del primo spettro
reference_wave = spectra_list[0]["WAVELENGTH"].values
flux_matrix = []

for df_s in spectra_list:
    # Interpolazione lineare del flusso sulla griglia di riferimento
    interp_flux = np.interp(reference_wave, df_s["WAVELENGTH"].values, df_s["FLUX"].values, left=0, right=0)
    flux_matrix.append(interp_flux)

# Matrice degli stati X di dimensioni (N_wavelengths, N_snapshots)
X_matrix = np.column_stack(flux_matrix)
print(f"Matrice degli stati costruita: {X_matrix.shape} (Lunghezze d'onda x Epoche/Snapshots)")


# =================================================
# 2. IMPLEMENTAZIONE HODMD (Higher-Order Dynamic Mode Decomposition)
# =================================================

def compute_hodmd(X, d=2, rank=None):
    """
    Esegue l'HODMD sulla matrice degli stati X di dimensioni (M, N).
    d: parametro di embedding per l'ordine superiore (Higher-Order)
    """
    M, N = X.shape
    
    # 1. Costruzione della matrice estesa con ritardi (Hankel-like stacking)
    # Numero di snapshot disponibili dopo l'embedding
    N_prime = N - d
    if N_prime < 2:
        raise ValueError("Troppi pochi snapshot temporali rispetto al parametro d scelto.")
        
    X_stack_list = [X[:, i:i+N_prime] for i in range(d)]
    X_stacked = np.vstack(X_stack_list) # Dimensioni: (d * M, N_prime)
    
    # Dividiamo in X1 e X2 per la scomposizione dinamica standard
    X1 = X_stacked[:, :-1]
    X2 = X_stacked[:, 1:]
    
    # 2. SVD di X1
    U, S, Vt = np.linalg.svd(X1, full_matrices=False)
    
    if rank is not None and rank < len(S):
        U = U[:, :rank]
        S = S[:rank]
        Vt = Vt[:rank, :]
        
    S_inv = np.diag(1.0 / S)
    
    # 3. Matrice Atilde nel sottospazio ridotto
    # A_tilde = U.T * X2 * V * S_inv
    A_tilde = np.dot(U.T, np.dot(X2, Vt.T))
    A_tilde = np.dot(A_tilde, S_inv)
    
    # 4. Autovalori e autovettori di A_tilde (Modi di Koopman discreti)
    eigenvalues, eigenvectors = np.linalg.eig(A_tilde)
    
    return eigenvalues, eigenvectors, U, S, Vt


print("\n[2/3] Esecuzione HODMD sulla matrice dei flussi UV...")
try:
    # Impostiamo d=2 (o adattivo in base al numero di snapshot)
    d_param = min(2, X_matrix.shape[1] - 2)
    if d_param < 1: 
        d_param = 1
        
    eigenvalues, eigenvectors, U, S, Vt = compute_hodmd(X_matrix, d=d_param)
    print("HODMD eseguito con successo.")
except Exception as e:
    raise RuntimeError(f"Errore durante il calcolo HODMD: {e}")


# =================================================
# 3. CALCOLO DEL KPI SPETTRALE
# =================================================

print("\n[3/3] Calcolo delle metriche KPI spettrali di Koopman...")

kpi_results = []

# Analizziamo gli autovalori (ampiezze, frequenze e tassi di crescita/smorzamento)
# Nel tempo discreto, l'autovalore lambda corrisponde a exp((omega + i*freq) * dt)
for i, lam in enumerate(eigenvalues):
    magnitude = np.abs(lam)
    phase = np.angle(lam) # Frequenza / fase legata alla dinamica temporale
    
    # KPI di stabilità: quanto l'autovalore si discosta dal cerchio unitario (magnitudo = 1)
    # Valori vicini a 1 indicano modi stabili/conservativi, lontani indicano smorzamento o crescita rapida
    stability_metric = abs(magnitude - 1.0)
    
    # Energia o peso associato al modo basato sui valori singolari
    mode_energy = S[i] if i < len(S) else 0.0
    
    # KPI Globale combinato (Stabilità + Frequenza + Ponderazione energetica)
    kpi_score = float(mode_energy * (1.0 / (1.0 + stability_metric)) * (1.0 + abs(phase)))
    
    kpi_results.append({
        "mode_index": i,
        "eigenvalue_real": float(lam.real),
        "eigenvalue_imag": float(lam.imag),
        "magnitude": float(magnitude),
        "phase_frequency": float(phase),
        "stability_metric": float(stability_metric),
        "mode_energy": float(mode_energy),
        "koopman_kpi_score": float(kpi_score)
    })

df_kpi = pd.DataFrame(kpi_results)
df_kpi = df_kpi.sort_values(by="koopman_kpi_score", ascending=False).reset_index(drop=True)

df_kpi.to_csv(OUTPUT_KPI_CSV, index=False)

print("\n" + "="*50)
print(" ANALISI KOOPMAN / KPI COMPLETATA")
print(f" Risultati salvati in: {OUTPUT_KPI_CSV}")
print("="*50)
print(df_kpi[["mode_index", "magnitude", "phase_frequency", "koopman_kpi_score"]].head(10))

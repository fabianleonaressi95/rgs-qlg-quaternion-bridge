#!/usr/bin/env python3

"""
HST STIS & COS SPECTRUM EXTRACTOR
Converte i FITS spettroscopici in tabelle di flusso Wavelength-Flux
per l'analisi di Koopman / HODMD
"""

import os
import glob
from astropy.io import fits
import pandas as pd
import numpy as np

# =================================================
# CONFIGURAZIONE
# =================================================

INPUT_DIR = "mast_uv_data"
OUTPUT_DIR = "spectra_extracted"

os.makedirs(OUTPUT_DIR, exist_ok=True)

print("="*70)
print(" HST STIS/COS SPECTRUM EXTRACTOR")
print("="*70)


# =================================================
# FUNZIONI DI ESTRAZIONE SPECIFICHE PER STRUMENTO
# =================================================

def extract_stis_spectrum(hdul):
    """
    Estrae wavelength e flux da un file FITS STIS (usualmente extension 1 o EXTTAB)
    """
    try:
        # In molti file STIS calibrati (x1d), i dati spettrali 1D sono nell'estensione 1
        data = hdul[1].data
        # Cerca colonne comuni nei file x1d di STIS
        cols = data.names
        
        wavelength_col = next((c for c in ['WAVELENGTH', 'WAVELNTH', 'WAVE'] if c in cols), None)
        flux_col = next((c for c in ['FLUX', 'NET'] if c in cols), None)
        
        if wavelength_col and flux_col:
            # Se ci sono più segmenti/ordini, prendiamo il primo array valido o uniamo
            wave = data[wavelength_col]
            flux = data[flux_col]
            
            # Gestione array multidimensionali (es. se ci sono aperture multiple)
            if wave.ndim > 1:
                wave = wave[0]
                flux = flux[0]
                
            return np.array(wave, dtype=float), np.array(flux, dtype=float)
    except Exception as e:
        print(f"  [Warning STIS] Errore estrazione interna: {e}")
    return None, None


def extract_cos_spectrum(hdul):
    """
    Estrae wavelength e flux da un file FITS COS (x1dcorrtab / x1d)
    """
    try:
        # Per COS, i dati 1D estratti sono solitamente nella prima estensione utile con tabella
        for ext in range(1, len(hdul)):
            data = hdul[ext].data
            if data is not None and 'WAVELENGTH' in data.names and 'FLUX' in data.names:
                wave = data['WAVELENGTH']
                flux = data['FLUX']
                if wave.ndim > 1:
                    wave = wave[0]
                    flux = flux[0]
                return np.array(wave, dtype=float), np.array(flux, dtype=float)
    except Exception as e:
        print(f"  [Warning COS] Errore estrazione interna: {e}")
    return None, None


# =================================================
# ELABORAZIONE FILE LOCALI
# =================================================

fits_files = glob.glob(os.path.join(INPUT_DIR, "*.fits"))
print(f"File FITS trovati in {INPUT_DIR}: {len(fits_files)}\n")

summary_manifest = []

for filepath in fits_files:
    filename = os.path.basename(filepath)
    print(f"Analisi file: {filename}")
    
    try:
        with fits.open(filepath) as hdul:
            hdr_primary = hdul[0].header
            instrument = hdr_primary.get("INSTRUME", "").upper()
            target = hdr_primary.get("TARGNAME", "UNKNOWN")
            date_obs = hdr_primary.get("DATE-OBS", "UNKNOWN")
            
            print(-> Strumento rilevato: {instrument})
            
            wave, flux = None, None
            
            if "STIS" in instrument:
                wave, flux = extract_stis_spectrum(hdul)
            elif "COS" in instrument:
                wave, flux = extract_cos_spectrum(hdul)
            else:
                print("  Strumento non supportato per estrazione spettrale automatica 1D in questo script.")
                continue
                
            if wave is not None and flux is not None and len(wave) > 0:
                # Pulizia NaN o valori negativi spuri se necessario
                valid_mask = ~np.isnan(wave) & ~np.isnan(flux)
                wave = wave[valid_mask]
                flux = flux[valid_mask]
                
                # Creazione DataFrame pulito
                df_spec = pd.DataFrame({
                    "WAVELENGTH": wave,
                    "FLUX": flux
                })
                
                # Ordinamento per lunghezza d'onda
                df_spec = df_spec.sort_values("WAVELENGTH").reset_index(drop=True)
                
                # Salvataggio CSV dedicato
                out_name = f"{os.path.splitext(filename)[0]}_spectrum.csv"
                out_path = os.path.join(OUTPUT_DIR, out_name)
                df_spec.to_csv(out_path, index=False)
                
                print(f"  Spettro estratto con successo: {len(df_spec)} punti salvati in {out_path}")
                
                summary_manifest.append({
                    "file_fits": filename,
                    "instrument": instrument,
                    "target": target,
                    "date": date_obs,
                    "csv_spectrum": out_name,
                    "points": len(df_spec)
                })
            else:
                print("  Nessun dato spettrale 1D valido trovato nelle estensioni.")
                
    except Exception as e:
        print(f"  Errore critico durante l'apertura di {filename}: {e}")

# Salvataggio sommario globale
df_summary = pd.DataFrame(summary_manifest)
summary_csv = "spectra_summary_manifest.csv"
df_summary.to_csv(summary_csv, index=False)

print("\n" + "="*30)
print(" ESTRAZIONE COMPLETATA")
print(f" Totale spettri estratti: {len(df_summary)}")
print(f" Manifesto salvato in: {summary_csv}")
print("="*30)

#!/usr/bin/env python3

"""
MAST UV DOWNLOAD & SPECTRUM EXTRACTION PIPELINE (ROBUST COORDS)
Ricerca per coordinate e estrazione automatica spettri HST (COS/STIS).
"""

import os
from astroquery.mast import Observations
from astropy.io import fits
import pandas as pd
import numpy as np
from tqdm import tqdm

# =================================================
# CONFIGURAZIONE
# =================================================

TARGET_COORDS = "22h20m09s +18d58m53s"
SEARCH_RADIUS = "0.02 deg"
DOWNLOAD_DIR = "mast_uv_data"
OUTPUT_DIR = "spectra_extracted"

os.makedirs(DOWNLOAD_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("="*70)
print(" MAST UV PIPELINE: ROBUST COORDINATES QUERY")
print("="*70)

# =================================================
# 1. QUERY PER COORDINATE E MISSIONE
# =================================================

print(f"\n[1/3] Interrogazione MAST per coordinate: {TARGET_COORDS}...")

obs = Observations.query_criteria(
    coordinates=TARGET_COORDS,
    radius=SEARCH_RADIUS,
    obs_collection="HST"
)

print(f"Osservazioni trovate: {len(obs)}")
if len(obs) == 0:
    raise RuntimeError("Nessuna osservazione trovata per le coordinate specificate.")

print("Recupero lista prodotti...")
products = Observations.get_product_list(obs)

# Filtro flessibile per prodotti spettrali 1D (x1d)
spec_products = products[
    [
        "x1d" in str(x).lower()
        for x in products["productFilename"]
    ]
]

spec_products = spec_products[:10]  # Limite di sicurezza per il test
print(f"Prodotti spettrali validi selezionati: {len(spec_products)}")

if len(spec_products) == 0:
    raise RuntimeError("Nessun file spettrale compatibile trovato tra i prodotti restituiti.")


# =================================================
# 2. DOWNLOAD DEI FILE
# =================================================

print("\n[2/3] Download in corso...")
manifest = []

for row in tqdm(spec_products):
    filename = row["productFilename"]
    try:
        result = Observations.download_products(row, download_dir=DOWNLOAD_DIR)
        if len(result) > 0 and "Local Path" in result.colnames:
            local_path = result["Local Path"][0]
            manifest.append({
                "file": filename,
                "path": local_path
            })
    except Exception as e:
        print(f"Errore download {filename}: {e}")

print(f"Download completati: {len(manifest)}")


# =================================================
# 3. ESTRAZIONE SPETTRALE (STIS / COS)
# =================================================

print("\n[3/3] Estrazione spettri in formato Wavelength-Flux...")

summary_manifest = []

for item in manifest:
    filepath = item["path"]
    filename = item["file"]
    
    try:
        with fits.open(filepath) as hdul:
            hdr_primary = hdul[0].header
            instrument = hdr_primary.get("INSTRUME", "").upper()
            target_name = hdr_primary.get("TARGNAME", "HD 209458")
            date_obs = hdr_primary.get("DATE-OBS", "UNKNOWN")
            
            wave, flux = None, None
            
            for ext in range(1, len(hdul)):
                data = hdul[ext].data
                if data is not None and hasattr(data, "names"):
                    cols = data.names
                    w_col = next((c for c in ['WAVELENGTH', 'WAVELNTH', 'WAVE'] if c in cols), None)
                    f_col = next((c for c in ['FLUX', 'NET'] if c in cols), None)
                    
                    if w_col and f_col:
                        wave = data[w_col]
                        flux = data[f_col]
                        if wave.ndim > 1:
                            wave = wave[0]
                            flux = flux[0]
                        break
            
            if wave is not None and flux is not None and len(wave) > 0:
                valid = ~np.isnan(wave) & ~np.isnan(flux)
                df_spec = pd.DataFrame({
                    "WAVELENGTH": np.array(wave[valid], dtype=float),
                    "FLUX": np.array(flux[valid], dtype=float)
                }).sort_values("WAVELENGTH").reset_index(drop=True)
                
                out_name = f"{os.path.splitext(filename)[0]}_spectrum.csv"
                out_path = os.path.join(OUTPUT_DIR, out_name)
                df_spec.to_csv(out_path, index=False)
                
                summary_manifest.append({
                    "file_fits": filename,
                    "instrument": instrument,
                    "target": target_name,
                    "date": date_obs,
                    "csv_spectrum": out_name,
                    "points": len(df_spec)
                })
                print(f"  -> Estratto con successo {filename} ({instrument}): {len(df_spec)} punti.")
            else:
                print(f"  -> Nessuna colonna spettrale valida trovata in {filename}")
    except Exception as e:
        print(f"Errore elaborazione FITS {filename}: {e}")

df_summary = pd.DataFrame(summary_manifest)
df_summary.to_csv("spectra_summary_manifest.csv", index=False)

print("\n" + "="*50)
print(" PIPELINE COMPLETATA CON SUCCESSO")
print(f" Spettri pronti in: {OUTPUT_DIR}/")
print(f" Totle file validi: {len(df_summary)}")
print("="*50)

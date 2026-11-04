#!/usr/bin/env python3

"""
MAST UV DATA DOWNLOADER
HST / GALEX

Scarica prodotti FITS dal Mikulski Archive
e prepara dataset UV per analisi Koopman
"""

import os
from astroquery.mast import Observations
from astropy.io import fits
from tqdm import tqdm


# ============================
# CONFIGURAZIONE
# ============================

TARGET = "M31"
MISSION = "HST"
DOWNLOAD_DIR = "mast_data"

os.makedirs(DOWNLOAD_DIR, exist_ok=True)


# ============================
# RICERCA MAST
# ============================

print("=" * 70)
print(" MAST UV DATA QUERY")
print("=" * 70)

obs = Observations.query_criteria(
    objectname=TARGET,
    obs_collection=MISSION
)

print(f"Trovate osservazioni: {len(obs)}")

if len(obs) == 0:
    raise RuntimeError("Nessun dato trovato")


# ============================
# LISTA PRODOTTI
# ============================

products = Observations.get_product_list(obs)
print(f"Prodotti disponibili: {len(products)}")

# Filtro sicuro per estensione .fits (gestisce eventuali maiuscole/minuscole)
fits_products = products[
    [
        str(x).lower().endswith(".fits")
        for x in products["productFilename"]
    ]
]

print(f"FITS trovati: {len(fits_products)}")


# ============================
# DOWNLOAD
# ============================

print("\nDownload FITS...")

manifest = []

for row in tqdm(fits_products):
    filename = row["productFilename"]
    
    try:
        result = Observations.download_products(
            row,
            download_dir=DOWNLOAD_DIR
        )
        
        # Gestione sicura del percorso restituito
        if len(result) > 0 and "Local Path" in result.colnames:
            local = result["Local Path"][0]
            manifest.append({
                "file": filename,
                "path": local
            })
            
    except Exception as e:
        print(f"Errore su {filename}: {e}")

print(f"\nDownload completato: {len(manifest)}")


# ============================
# CONTROLLO HEADER FITS
# ============================

print("\nAnalisi header FITS")

for item in manifest[:10]:
    print("\nFILE:", item["file"])
    
    try:
        with fits.open(item["path"]) as hdul:
            hdr = hdul[0].header
            print("TELESCOP:", hdr.get("TELESCOP", "N/D"))
            print("INSTRUME:", hdr.get("INSTRUME", "N/D"))
            print("FILTER:", hdr.get("FILTER", "N/D"))
            print("DATE:", hdr.get("DATE-OBS", "N/D"))
    except Exception as e:
        print(f"Errore lettura FITS per {item['file']}: {e}")

print("\nPipeline terminata.")

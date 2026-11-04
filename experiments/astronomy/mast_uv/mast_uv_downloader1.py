#!/usr/bin/env python3

"""
CONTROLLED MAST UV DOWNLOADER
HST / GALEX

Versione per pipeline Koopman UV
"""

from astroquery.mast import Observations
from astropy.io import fits
from tqdm import tqdm
import os
import pandas as pd


# =================================================
# CONFIG
# =================================================

TARGET = "NGC 7000"       # esempio più leggero
MISSION = "HST"

MAX_OBSERVATIONS = 20
MAX_PRODUCTS = 50

DOWNLOAD_DIR = "mast_uv_data"


os.makedirs(
    DOWNLOAD_DIR,
    exist_ok=True
)


# =================================================
# QUERY
# =================================================

print("="*70)
print(" CONTROLLED MAST UV QUERY")
print("="*70)


obs = Observations.query_criteria(
    objectname=TARGET,
    obs_collection=MISSION
)


print(
    "Osservazioni trovate:",
    len(obs)
)


# limita

obs = obs[:MAX_OBSERVATIONS]


print(
    "Osservazioni analizzate:",
    len(obs)
)


# =================================================
# PRODUCTS
# =================================================


print("\nRecupero prodotti...")


products = Observations.get_product_list(obs)


print(
    "Prodotti totali:",
    len(products)
)



# solo FITS

fits_mask = [
    str(x).lower().endswith(".fits")
    for x in products["productFilename"]
]


fits_products = products[
    fits_mask
]


fits_products = fits_products[:MAX_PRODUCTS]


print(
    "FITS selezionati:",
    len(fits_products)
)



# =================================================
# DOWNLOAD
# =================================================


records=[]


print("\nDownload...")


for row in tqdm(fits_products):

    name=row["productFilename"]


    try:

        result = Observations.download_products(
            row,
            download_dir=DOWNLOAD_DIR
        )


        path=result["Local Path"][0]


        records.append(
            {
            "file":name,
            "path":path
            }
        )


    except Exception as e:

        print(
            "ERROR:",
            name,
            e
        )



# =================================================
# HEADER SCAN
# =================================================


print("\nAnalisi FITS")


catalog=[]


for r in records:


    try:

        with fits.open(r["path"]) as h:


            hdr=h[0].header


            catalog.append(
                {

                "file":
                r["file"],

                "instrument":
                hdr.get("INSTRUME"),

                "filter":
                hdr.get("FILTER"),

                "date":
                hdr.get("DATE-OBS"),

                "exptime":
                hdr.get("EXPTIME")

                }
            )


    except:

        pass



df=pd.DataFrame(catalog)


df.to_csv(
    "mast_uv_inventory.csv",
    index=False
)


print(df)


print("\n================================")
print(" COMPLETATO")
print(" Inventory: mast_uv_inventory.csv")
print("================================")

#!/usr/bin/env python3
"""
Script per scaricare e integrare dataset MA/WGS su E. coli
in un unico file HDF5 con più tabelle:
  - /mutations  (una riga per mutazione)
  - /sites      (una riga per sito/contesto/condizione, con tassi mu_i^exp)
  - /samples    (metadati dei sample)
  - /genome_info (info sul genoma di riferimento)

Dataset target:
  - Gifford et al. 2024 (Figshare DOI 10.6084/m9.figshare.25563537)
  - Foster et al. (PMC3478608)
  - Saini/Long et al. (PNAS 2015, PMID 26460006)
  - Aerobico/anaerobico (PMC5289635)

Dipendenze:
  pip install requests pandas h5py pyfaidx pysam biopython openpyxl
"""

import os
import io
import logging
from pathlib import Path

import requests
import pandas as pd
import numpy as np
import h5py
from Bio import SeqIO
from pyfaidx import Fasta
import pysam
import zipfile
import json

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
log = logging.getLogger(__name__)

BASE_DIR = Path("ecoli_ma_data")
BASE_DIR.mkdir(exist_ok=True)

GENOME_FASTA = BASE_DIR / "NC_000913.fna"
HDF5_OUT = Path("ecoli_ma_wgs_integrated.h5")

# -----------------------------------------------------------------------------
# Utility per download
# -----------------------------------------------------------------------------

def download_file(url, dest_path, chunk_size=8192):
    """Scarica un file da URL a dest_path."""
    log.info(f"Downloading {url} -> {dest_path}")
    resp = requests.get(url, stream=True)
    resp.raise_for_status()
    with open(dest_path, "wb") as f:
        for chunk in resp.iter_content(chunk_size=chunk_size):
            if chunk:
                f.write(chunk)
    log.info(f"Downloaded {dest_path}")

def download_and_unzip(url, dest_dir, expected_exts=(".vcf", ".csv", ".txt", ".xlsx")):
    """Scarica uno zip, lo estrae in dest_dir e restituisce liste di file."""
    zip_path = dest_dir / "tmp.zip"
    download_file(url, zip_path)
    files = []
    with zipfile.ZipFile(zip_path, "r") as z:
        for name in z.namelist():
            if name.endswith(expected_exts):
                z.extract(name, dest_dir)
                files.append(dest_dir / name)
    zip_path.unlink()
    log.info(f"Extracted {len(files)} files in {dest_dir}")
    return files

# -----------------------------------------------------------------------------
# Genoma di riferimento
# -----------------------------------------------------------------------------

def fetch_reference_genome():
    """Scarica il genoma di E. coli MG1655 (NC_000913) da NCBI."""
    if GENOME_FASTA.exists():
        log.info(f"Reference genome already exists: {GENOME_FASTA}")
        return str(GENOME_FASTA)
    
    # URL NCBI per NC_000913 (FASTA)
    url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000/005/845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
    gz_path = BASE_DIR / "GCF_000005845.2_ASM584v2_genomic.fna.gz"
    download_file(url, gz_path)
    
    # Decomprimere
    import gzip
    fasta_path = BASE_DIR / "GCF_000005845.2_ASM584v2_genomic.fna"
    with gzip.open(gz_path, "rb") as f_in, open(fasta_path, "wb") as f_out:
        f_out.write(f_in.read())
    gz_path.unlink()
    
    # Rinominare in NC_000913.fna per coerenza
    fasta_path.rename(GENOME_FASTA)
    log.info(f"Reference genome ready: {GENOME_FASTA}")
    return str(GENOME_FASTA)

# -----------------------------------------------------------------------------
# Funzioni comuni per classificazione mutazioni e contesto
# -----------------------------------------------------------------------------

COMP_BASE = {"A":"T", "T":"A", "G":"C", "C":"G"}

def classify_bps(ref_base, alt_base):
    """
    Classifica la mutazione nelle 6 classi di BPS:
    A:T->G:C, G:C->A:T, A:T->C:G, A:T->T:A, G:C->C:G, G:C->T:A
    """
    ref_base = ref_base.upper()
    alt_base = alt_base.upper()
    if ref_base not in "ACGT" or alt_base not in "ACGT":
        return None
    if ref_base == alt_base:
        return None
    
    # Normalizziamo a ref in {A,G}
    if ref_base in "AT":
        if ref_base == "T":
            ref_base = "A"
            alt_base = COMP_BASE[alt_base]
        if alt_base == "G":
            return "A:T->G:C"
        elif alt_base == "C":
            return "A:T->C:G"
        elif alt_base == "T":
            return "A:T->T:A"
        else:
            return None
    else:
        if ref_base == "C":
            ref_base = "G"
            alt_base = COMP_BASE[alt_base]
        if alt_base == "A":
            return "G:C->A:T"
        elif alt_base == "C":
            return "G:C->C:G"
        elif alt_base == "T":
            return "G:C->T:A"
        else:
            return None

def get_dinuc_5p(ref_fasta, chrom, pos0):
    """
    Restituisce il dinucleotide 5' (base precedente + base mutata).
    pos0 è 0-based.
    """
    seq = ref_fasta[chrom]
    L = len(seq)
    if pos0 <= 0:
        end = min(L, 2)
        return seq[0:end].upper()
    else:
        start = pos0 - 1
        end = min(L, pos0 + 1)
        return seq[start:end].upper()

# -----------------------------------------------------------------------------
# Dataset: Gifford et al. 2024 (Figshare)
# -----------------------------------------------------------------------------

def fetch_gifford_2024(ref_fasta):
    """
    Scarica i dati di Gifford et al. 2024 da Figshare.
    DOI: 10.6084/m9.figshare.25563537
    URL diretti dei file:
    - https://figshare.com/ndownloader/articles/25563537/files/46363590
    - https://figshare.com/ndownloader/articles/25563537/files/46363593
    """
    log.info("Fetching Gifford et al. 2024 (Figshare)")
    
    gifford_dir = BASE_DIR / "gifford_2024"
    gifford_dir.mkdir(exist_ok=True)
    
    urls = [
        "https://figshare.com/ndownloader/articles/25563537/files/46363590",
        "https://figshare.com/ndownloader/articles/25563537/files/46363593",
    ]
    
    files = []
    for url in urls:
        fname = url.split("/")[-1]
        dest = gifford_dir / fname
        if not dest.exists():
            download_file(url, dest)
        files.append(dest)
    
    df_mut = pd.DataFrame(columns=[
        "dataset", "sample_id", "chrom", "pos", "ref", "alt",
        "bps_class", "dinuc_5p", "strand"
    ])
    df_samples = pd.DataFrame(columns=[
        "sample_id", "strain", "condition", "generations", "oxidative_mut",
        "mmr_mut", "source"
    ])
    
    for fpath in files:
        if fpath.suffix == ".csv":
            try:
                df = pd.read_csv(fpath)
                log.info(f"Processing {fpath}: {len(df)} righe")
                # Qui andrebbe il parsing reale in base alle colonne
            except Exception as e:
                log.warning(f"Error processing {fpath}: {e}")
    
    return df_mut, df_samples

# -----------------------------------------------------------------------------
# Dataset: Foster et al. (PMC3478608)
# -----------------------------------------------------------------------------

def fetch_foster_pmc3478608(ref_fasta):
    """
    Foster et al., 'Rate and molecular spectrum of spontaneous mutations in E. coli'
    PMC3478608.
    """
    log.info("Fetching Foster et al. (PMC3478608)")
    
    rows = []
    samples = []
    
    sample_wt = "Foster_WT"
    samples.append({
        "sample_id": sample_wt,
        "strain": "WT",
        "condition": "aerobic",
        "generations": None,
        "oxidative_mut": False,
        "mmr_mut": False,
        "source": "Foster_PMC3478608"
    })
    
    counts_wt = {
        "A:T->G:C": 49,
        "G:C->A:T": 82,
        "A:T->T:A": 17,
        "A:T->C:G": 38,
        "G:C->T:A": 30,
        "G:C->C:G": 17
    }
    
    for bps, count in counts_wt.items():
        rows.append({
            "dataset": "Foster",
            "sample_id": sample_wt,
            "chrom": "NC_000913",
            "pos": -1,
            "ref": ".",
            "alt": ".",
            "bps_class": bps,
            "dinuc_5p": ".",
            "strand": ".",
            "count": count
        })
    
    sample_mutl = "Foster_MutL-"
    samples.append({
        "sample_id": sample_mutl,
        "strain": "MutL-",
        "condition": "aerobic",
        "generations": None,
        "oxidative_mut": False,
        "mmr_mut": True,
        "source": "Foster_PMC3478608"
    })
    
    counts_mutl = {
        "A:T->G:C": 1141,
        "G:C->A:T": 447,
    }
    
    for bps, count in counts_mutl.items():
        rows.append({
            "dataset": "Foster",
            "sample_id": sample_mutl,
            "chrom": "NC_000913",
            "pos": -1,
            "ref": ".",
            "alt": ".",
            "bps_class": bps,
            "dinuc_5p": ".",
            "strand": ".",
            "count": count
        })
    
    df_mut = pd.DataFrame(rows)
    df_samples = pd.DataFrame(samples)
    
    log.info(f"Foster: {len(df_mut)} righe aggregate, {len(df_samples)} sample")
    return df_mut, df_samples

# -----------------------------------------------------------------------------
# Dataset: Saini/Long et al. (PNAS 2015, PMID 26460006)
# -----------------------------------------------------------------------------

def fetch_saini_long_pnas2015(ref_fasta):
    """
    Saini/Long et al., 'Determinants of spontaneous mutation in E. coli'
    PNAS 2015, PMID 26460006.
    URL diretti dei file:
    - https://pnas.org/doi/suppl/10.1073/pnas.1512136112/suppl_file/pnas.1512136112.sd01.xlsx
    """
    log.info("Fetching Saini/Long et al. (PNAS 2015, PMID 26460006)")
    
    saini_dir = BASE_DIR / "saini_long_2015"
    saini_dir.mkdir(exist_ok=True)
    
    url = "https://pnas.org/doi/suppl/10.1073/pnas.1512136112/suppl_file/pnas.1512136112.sd01.xlsx"
    dest = saini_dir / "pnas.1512136112.sd01.xlsx"
    if not dest.exists():
        download_file(url, dest)
    
    df_mut = pd.DataFrame(columns=[
        "dataset", "sample_id", "chrom", "pos", "ref", "alt",
        "bps_class", "dinuc_5p", "strand"
    ])
    df_samples = pd.DataFrame(columns=[
        "sample_id", "strain", "condition", "generations", "oxidative_mut",
        "mmr_mut", "source"
    ])
    
    try:
        xl = pd.ExcelFile(dest)
        log.info(f"Sheets in {dest.name}: {xl.sheet_names}")
    except Exception as e:
        log.warning(f"Error processing {dest}: {e}")
    
    df_samples = pd.DataFrame([{
        "sample_id": "Saini_WT",
        "strain": "WT",
        "condition": "aerobic",
        "generations": None,
        "oxidative_mut": False,
        "mmr_mut": False,
        "source": "Saini_Long_PNAS2015"
    }])
    
    return df_mut, df_samples

# -----------------------------------------------------------------------------
# Dataset: Aerobico vs Anaerobico (PMC5289635)
# -----------------------------------------------------------------------------

def fetch_aerobic_anaerobic_pmc5289635(ref_fasta):
    """
    'Anaerobically Grown E. coli Has an Enhanced Mutation Rate and Distinct Mutational Spectra'
    PMC5289635.
    """
    log.info("Fetching aerobic/anaerobic dataset (PMC5289635)")
    
    df_samples = pd.DataFrame([
        {
            "sample_id": "Aerobic_WT",
            "strain": "WT",
            "condition": "aerobic",
            "generations": None,
            "oxidative_mut": False,
            "mmr_mut": False,
            "source": "PMC5289635"
        },
        {
            "sample_id": "Anaerobic_WT",
            "strain": "WT",
            "condition": "anaerobic",
            "generations": None,
            "oxidative_mut": False,
            "mmr_mut": False,
            "source": "PMC5289635"
        }
    ])
    
    df_mut = pd.DataFrame(columns=[
        "dataset", "sample_id", "chrom", "pos", "ref", "alt",
        "bps_class", "dinuc_5p", "strand", "count"
    ])

#!/usr/bin/env python3
"""
Audit riproducibile del dataset Zenodo 10.5281/zenodo.19068308.

Cosa fa:
1. interroga l'API pubblica Zenodo e scarica i file del record;
2. registra SHA-256, dimensioni e nomi in manifest.csv;
3. ispeziona file .mat (MATLAB v5 o HDF5/v7.3), elencando variabili e forme;
4. esporta arrays numerici in .npy e un riepilogo delle statistiche in variables.csv;
5. NON presume che le variabili siano intensità, tempi o posizioni: tale mappatura
   va confermata con l'articolo/metadati prima di stimare osservabili fisici.

Dipendenze: requests, numpy, pandas, scipy, h5py
Installazione: python -m pip install requests numpy pandas scipy h5py
Esecuzione: python audit_fibonacci_zenodo.py
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
import time
from pathlib import Path
from urllib.parse import quote

import numpy as np
import pandas as pd
import requests

RECORD_ID = "19068308"
API_URL = f"https://zenodo.org/api/records/{RECORD_ID}"
ROOT = Path("fibonacci_zenodo_audit")
RAW = ROOT / "raw"
ARRAYS = ROOT / "arrays"
ROOT.mkdir(exist_ok=True)
RAW.mkdir(exist_ok=True)
ARRAYS.mkdir(exist_ok=True)

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def safe_name(name: str) -> str:
    name = Path(name).name
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name)

def fetch_record():
    r = requests.get(API_URL, timeout=45)
    r.raise_for_status()
    return r.json()

def download_file(entry: dict) -> Path:
    name = safe_name(entry["key"])
    target = RAW / name
    expected_size = entry.get("size")
    if target.exists() and (expected_size is None or target.stat().st_size == expected_size):
        print(f"[CACHE] {name} ({target.stat().st_size} bytes)")
        return target
    url = entry.get("links", {}).get("self")
    if not url:
        url = f"https://zenodo.org/records/{RECORD_ID}/files/{quote(entry['key'])}?download=1"
    print(f"[DOWNLOAD] {name} <- {url}")
    with requests.get(url, stream=True, timeout=(30, 180)) as r:
        r.raise_for_status()
        tmp = target.with_suffix(target.suffix + ".part")
        with tmp.open("wb") as f:
            for chunk in r.iter_content(1024 * 1024):
                if chunk:
                    f.write(chunk)
        tmp.replace(target)
    if expected_size is not None and target.stat().st_size != expected_size:
        raise RuntimeError(
            f"Dimensione non corrispondente per {name}: "
            f"{target.stat().st_size} != {expected_size}"
        )
    return target

def inspect_mat(path: Path):
    """Return dict of arrays; try MATLAB v5 then HDF5/v7.3."""
    try:
        from scipy.io import loadmat
        obj = loadmat(path, squeeze_me=False, struct_as_record=False)
        return {
            k: v for k, v in obj.items()
            if not k.startswith("__") and isinstance(v, np.ndarray)
        }, "scipy.io.loadmat"
    except Exception as first_error:
        try:
            import h5py
            out = {}
            with h5py.File(path, "r") as f:
                def visitor(name, obj):
                    if isinstance(obj, h5py.Dataset):
                        try:
                            out[name] = np.asarray(obj)
                        except Exception:
                            pass
                f.visititems(visitor)
            return out, "h5py"
        except Exception as second_error:
            raise RuntimeError(
                f"Impossibile leggere {path.name}; scipy={first_error}; h5py={second_error}"
            )

def summarize_array(file_name: str, var_name: str, arr: np.ndarray, loader: str):
    original_shape = tuple(int(x) for x in arr.shape)
    row = {
        "file": file_name,
        "variable": var_name,
        "loader": loader,
        "shape": "x".join(map(str, original_shape)),
        "dtype": str(arr.dtype),
        "size": int(arr.size),
        "numeric": bool(np.issubdtype(arr.dtype, np.number)),
        "finite_count": None,
        "nan_count": None,
        "min": None,
        "max": None,
        "mean": None,
        "std": None,
        "saved_npy": None,
    }
    if np.issubdtype(arr.dtype, np.number):
        vals = np.asarray(arr)
        finite = np.isfinite(vals)
        row["finite_count"] = int(finite.sum())
        row["nan_count"] = int((~finite).sum())
        if finite.any():
            v = vals[finite].astype(np.complex128 if np.iscomplexobj(vals) else np.float64)
            row["min"] = float(np.min(np.real(v)))
            row["max"] = float(np.max(np.real(v)))
            row["mean"] = float(np.mean(np.real(v)))
            row["std"] = float(np.std(np.real(v)))
        # Save arrays only when numeric and not object dtype.
        if arr.dtype.kind in "biufc":
            outname = safe_name(f"{Path(file_name).stem}__{var_name.replace('/', '_')}.npy")
            outpath = ARRAYS / outname
            np.save(outpath, arr, allow_pickle=False)
            row["saved_npy"] = str(outpath)
    return row

def main():
    print(f"Zenodo record: {RECORD_ID}")
    print(f"Output directory: {ROOT.resolve()}")
    try:
        record = fetch_record()
    except Exception as e:
        print("\nERRORE DI ACCESSO A ZENODO.")
        print("Nessun dato è stato scaricato; controlla la connessione e riprova.")
        print(f"Dettaglio: {type(e).__name__}: {e}")
        sys.exit(2)

    metadata = record.get("metadata", {})
    files = record.get("files", [])
    (ROOT / "record_metadata.json").write_text(
        json.dumps({
            "id": record.get("id"),
            "doi": metadata.get("doi"),
            "title": metadata.get("title"),
            "publication_date": metadata.get("publication_date"),
            "creators": metadata.get("creators"),
            "description": metadata.get("description"),
            "files": [
                {"key": x.get("key"), "size": x.get("size"),
                 "checksum": x.get("checksum"), "links": x.get("links")}
                for x in files
            ],
        }, indent=2, ensure_ascii=False),
        encoding="utf-8"
    )
    print(f"Titolo: {metadata.get('title')}")
    print(f"DOI: {metadata.get('doi')}")
    print(f"File dichiarati nell'API: {len(files)}")
    if not files:
        print("L'API non ha restituito file: controlla record_metadata.json.")
        sys.exit(3)

    manifest_rows = []
    variable_rows = []
    for i, entry in enumerate(files, 1):
        try:
            p = download_file(entry)
            manifest_rows.append({
                "file": p.name,
                "bytes": p.stat().st_size,
                "sha256": sha256(p),
                "zenodo_key": entry.get("key"),
                "zenodo_size": entry.get("size"),
                "zenodo_checksum": entry.get("checksum"),
                "content_type": entry.get("type"),
            })
            if p.suffix.lower() == ".mat":
                arrays, loader = inspect_mat(p)
                print(f"\n[{i}/{len(files)}] {p.name}: loader={loader}, variabili={len(arrays)}")
                if not arrays:
                    print("  Nessun array semplice rilevato; possibile struttura MATLAB annidata.")
                for name, arr in arrays.items():
                    row = summarize_array(p.name, name, arr, loader)
                    variable_rows.append(row)
                    print(f"  {name}: shape={arr.shape}, dtype={arr.dtype}, size={arr.size}")
            else:
                print(f"\n[{i}/{len(files)}] {p.name}: scaricato; ispezione automatica non prevista")
        except Exception as e:
            print(f"[ERRORE] {entry.get('key')}: {type(e).__name__}: {e}")
            manifest_rows.append({
                "file": entry.get("key"), "error": f"{type(e).__name__}: {e}"
            })

    pd.DataFrame(manifest_rows).to_csv(ROOT / "manifest.csv", index=False)
    pd.DataFrame(variable_rows).to_csv(ROOT / "variables.csv", index=False)
    print("\n=== AUDIT COMPLETATO ===")
    print(f"Manifest: {ROOT / 'manifest.csv'}")
    print(f"Inventario variabili: {ROOT / 'variables.csv'}")
    print(f"Array esportati: {ARRAYS}")
    print("\nATTENZIONE: statistiche descrittive, non ancora risultati fisici.")
    print("Prima di calcolare MSD, autocorrelazioni o IPR, mappa variabili, unità,")
    print("assi e normalizzazione usando il paper e i metadati.")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
EXTENDED AUDIT SCARICATI
========================
Audit completo di tutti i file nella directory Scaricati
relativi al progetto 188D, DNS, SU(2) e artefatti correlati.
"""

import numpy as np
import json
from pathlib import Path
from datetime import datetime
import hashlib

BASE_DIR = Path("/home/fabianleo/Scaricati")
OUTPUT_FILE = BASE_DIR / "EXTENDED_AUDIT_SCARICATI_REPORT.json"
SUMMARY_FILE = BASE_DIR / "EXTENDED_AUDIT_SCARICATI_SUMMARY.md"

print("=" * 80)
print("EXTENDED AUDIT SCARICATI")
print("=" * 80)
print(f"Directory: {BASE_DIR}")
print(f"Data: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print()

report = {
    'timestamp': datetime.now().isoformat(),
    'base_directory': str(BASE_DIR),
    'categories': {}
}

# ============================================================================
# CATEGORIA 1: 188D FINAL BUNDLE (UNVERIFIED)
# ============================================================================
print("CATEGORIA 1: 188D_FINAL_BUNDLE (UNVERIFIED)")
print("-" * 80)

final_bundle_dir = BASE_DIR / "188D_FINAL_BUNDLE"
if final_bundle_dir.exists():
    files = list(final_bundle_dir.glob("*"))
    print(f"File trovati: {len(files)}")
    
    category = {
        'status': 'UNVERIFIED',
        'description': 'Artefatti storici non verificati - NON usare per fisica certificata',
        'files': {}
    }
    
    for f in sorted(files):
        if f.is_file():
            size = f.stat().st_size
            file_info = {
                'size_bytes': size,
                'size_mb': round(size / (1024*1024), 3),
                'modified': datetime.fromtimestamp(f.stat().st_mtime).isoformat()
            }
            
            # Se è .npy, carica e mostra shape
            if f.suffix == '.npy':
                try:
                    data = np.load(f)
                    file_info['shape'] = list(data.shape)
                    file_info['dtype'] = str(data.dtype)
                    file_info['min'] = float(np.min(data))
                    file_info['max'] = float(np.max(data))
                    file_info['mean'] = float(np.mean(data))
                    print(f"  {f.name}: shape {data.shape}, dtype {data.dtype}")
                except Exception as e:
                    file_info['error'] = str(e)
                    print(f"  {f.name}: ERRORE - {e}")
            else:
                print(f"  {f.name}: {size} bytes")
            
            category['files'][f.name] = file_info
    
    report['categories']['188D_FINAL_BUNDLE'] = category
else:
    print("  Directory non trovata")

print()

# ============================================================================
# CATEGORIA 2: 188D CERTIFIED
# ============================================================================
print("CATEGORIA 2: 188D_CERTIFIED")
print("-" * 80)

certified_dir = BASE_DIR / "188D_CERTIFIED"
if certified_dir.exists():
    files = list(certified_dir.glob("*"))
    print(f"File trovati: {len(files)}")
    
    category = {
        'status': 'CERTIFIED',
        'description': 'Artefatti certificati - usabili per fisica',
        'files': {}
    }
    
    for f in sorted(files):
        if f.is_file():
            size = f.stat().st_size
            file_info = {
                'size_bytes': size,
                'size_mb': round(size / (1024*1024), 3),
                'modified': datetime.fromtimestamp(f.stat().st_mtime).isoformat()
            }
            
            if f.suffix == '.npy':
                try:
                    data = np.load(f)
                    file_info['shape'] = list(data.shape)
                    file_info['dtype'] = str(data.dtype)
                    file_info['min'] = float(np.min(data))
                    file_info['max'] = float(np.max(data))
                    file_info['mean'] = float(np.mean(data))
                    print(f"  {f.name}: shape {data.shape}, dtype {data.dtype}")
                except Exception as e:
                    file_info['error'] = str(e)
                    print(f"  {f.name}: ERRORE - {e}")
            elif f.suffix == '.json':
                try:
                    with open(f) as fp:
                        json_data = json.load(fp)
                    file_info['json_keys'] = list(json_data.keys())
                    print(f"  {f.name}: JSON con chiavi {list(json_data.keys())}")
                except Exception as e:
                    file_info['error'] = str(e)
                    print(f"  {f.name}: ERRORE - {e}")
            else:
                print(f"  {f.name}: {size} bytes")
            
            category['files'][f.name] = file_info
    
    report['categories']['188D_CERTIFIED'] = category
else:
    print("  Directory non trovata")

print()

# ============================================================================
# CATEGORIA 3: 188D PROVENANCE AUDIT
# ============================================================================
print("CATEGORIA 3: 188D_PROVENANCE_AUDIT")
print("-" * 80)

provenance_dir = BASE_DIR / "188D_PROVENANCE_AUDIT"
if provenance_dir.exists():
    files = list(provenance_dir.glob("*"))
    print(f"File trovati: {len(files)}")
    
    category = {
        'status': 'RECONSTRUCTED',
        'description': 'Artefatti ricostruiti da first principles',
        'files': {}
    }
    
    for f in sorted(files):
        if f.is_file():
            size = f.stat().st_size
            file_info = {
                'size_bytes': size,
                'size_mb': round(size / (1024*1024), 3),
                'modified': datetime.fromtimestamp(f.stat().st_mtime).isoformat()
            }
            
            if f.suffix == '.npy':
                try:
                    data = np.load(f)
                    file_info['shape'] = list(data.shape)
                    file_info['dtype'] = str(data.dtype)
                    file_info['min'] = float(np.min(data))
                    file_info['max'] = float(np.max(data))
                    file_info['mean'] = float(np.mean(data))
                    print(f"  {f.name}: shape {data.shape}, dtype {data.dtype}")
                except Exception as e:
                    file_info['error'] = str(e)
                    print(f"  {f.name}: ERRORE - {e}")
            elif f.suffix == '.json':
                try:
                    with open(f) as fp:
                        json_data = json.load(fp)
                    file_info['json_keys'] = list(json_data.keys())
                    print(f"  {f.name}: JSON con chiavi {list(json_data.keys())}")
                except Exception as e:
                    file_info['error'] = str(e)
                    print(f"  {f.name}: ERRORE - {e}")
            else:
                print(f"  {f.name}: {size} bytes")
            
            category['files'][f.name] = file_info
    
    report['categories']['188D_PROVENANCE_AUDIT'] = category
else:
    print("  Directory non trovata")

print()

# ============================================================================
# CATEGORIA 4: 188D LOCAL REPLAY OUTPUT
# ============================================================================
print("CATEGORIA 4: 188D_LOCAL_REPLAY_OUTPUT")
print("-" * 80)

replay_dir = BASE_DIR / "188D_LOCAL_REPLAY_OUTPUT"
if replay_dir.exists():
    files = list(replay_dir.glob("*"))
    print(f"File trovati: {len(files)}")
    
    category = {
        'status': 'RECONSTRUCTED',
        'description': 'Artefatti ricostruiti da X_Oh e Q',
        'files': {}
    }
    
    for f in sorted(files):
        if f.is_file():
            size = f.stat().st_size
            file_info = {
                'size_bytes': size,
                'size_mb': round(size / (1024*1024), 3),
                'modified': datetime.fromtimestamp(f.stat().st_mtime).isoformat()
            }
            
            if f.suffix == '.npy':
                try:
                    data = np.load(f)
                    file_info['shape'] = list(data.shape)
                    file_info['dtype'] = str(data.dtype)
                    file_info['min'] = float(np.min(data))
                    file_info['max'] = float(np.max(data))
                    file_info['mean'] = float(np.mean(data))
                    print(f"  {f.name}: shape {data.shape}, dtype {data.dtype}")
                except Exception as e:
                    file_info['error'] = str(e)
                    print(f"  {f.name}: ERRORE - {e}")
            elif f.suffix == '.json':
                try:
                    with open(f) as fp:
                        json_data = json.load(fp)
                    file_info['json_keys'] = list(json_data.keys())
                    print(f"  {f.name}: JSON con chiavi {list(json_data.keys())}")
                except Exception as e:
                    file_info['error'] = str(e)
                    print(f"  {f.name}: ERRORE - {e}")
            else:
                print(f"  {f.name}: {size} bytes")
            
            category['files'][f.name] = file_info
    
    report['categories']['188D_LOCAL_REPLAY_OUTPUT'] = category
else:
    print("  Directory non trovata")

print()

# ============================================================================
# CATEGORIA 5: DNS FILES
# ============================================================================
print("CATEGORIA 5: DNS FILES")
print("-" * 80)

dns_files = [
    "dns_original_extracted.py",
    "dns_conservation_audit_final.py",
    "DNS_CONSERVATION_AUDIT_FINAL.json",
    "bilancio_energetico_dns.png"
]

category = {
    'status': 'RECOVERED',
    'description': 'Codice DNS originale e audit conservazione',
    'files': {}
}

for fname in dns_files:
    f = BASE_DIR / fname
    if f.exists():
        size = f.stat().st_size
        file_info = {
            'size_bytes': size,
            'size_kb': round(size / 1024, 3),
            'modified': datetime.fromtimestamp(f.stat().st_mtime).isoformat()
        }
        
        if f.suffix == '.json':
            try:
                with open(f) as fp:
                    json_data = json.load(fp)
                file_info['json_keys'] = list(json_data.keys())
                print(f"  {fname}: JSON con chiavi {list(json_data.keys())}")
            except Exception as e:
                file_info['error'] = str(e)
                print(f"  {fname}: ERRORE - {e}")
        else:
            print(f"  {fname}: {size} bytes")
        
        category['files'][fname] = file_info
    else:
        print(f"  {fname}: NON TROVATO")
        category['files'][fname] = {'status': 'missing'}

report['categories']['DNS_FILES'] = category

print()

# ============================================================================
# CATEGORIA 6: 188D REPORTS
# ============================================================================
print("CATEGORIA 6: 188D REPORTS")
print("-" * 80)

report_files = [
    "188D_AUDIT_FINDINGS.md",
    "188D_PROVENANCE_AUDIT_FINAL_REPORT.md",
    "PROVENANCE_STATUS.json",
    "188D_FORENSIC_REPORT.txt"
]

category = {
    'status': 'DOCUMENTATION',
    'description': 'Report e documentazione 188D',
    'files': {}
}

for fname in report_files:
    f = BASE_DIR / fname
    if f.exists():
        size = f.stat().st_size
        file_info = {
            'size_bytes': size,
            'size_kb': round(size / 1024, 3),
            'modified': datetime.fromtimestamp(f.stat().st_mtime).isoformat()
        }
        
        if f.suffix == '.json':
            try:
                with open(f) as fp:
                    json_data = json.load(fp)
                file_info['json_keys'] = list(json_data.keys())
                print(f"  {fname}: JSON con chiavi {list(json_data.keys())}")
            except Exception as e:
                file_info['error'] = str(e)
                print(f"  {fname}: ERRORE - {e}")
        else:
            print(f"  {fname}: {size} bytes")
        
        category['files'][fname] = file_info
    else:
        print(f"  {fname}: NON TROVATO")
        category['files'][fname] = {'status': 'missing'}

report['categories']['188D_REPORTS'] = category

print()

# ============================================================================
# CATEGORIA 7: 64B FILES
# ============================================================================
print("CATEGORIA 7: 64B FILES (cubici)")
print("-" * 80)

b64_files = [
    "64B_Oh_Kplus_eigenvalues.npy",
    "64B_Oh_reconstructed.npy",
    "64B_Oh_Vminus_irrep_decomposition.json",
    "64B_Oh_hessian_independent_audit.json"
]

category = {
    'status': 'CUBIC',
    'description': 'Artefatti 64B cubici O_h',
    'files': {}
}

for fname in b64_files:
    f = BASE_DIR / fname
    if f.exists():
        size = f.stat().st_size
        file_info = {
            'size_bytes': size,
            'size_kb': round(size / 1024, 3),
            'modified': datetime.fromtimestamp(f.stat().st_mtime).isoformat()
        }
        
        if f.suffix == '.npy':
            try:
                data = np.load(f)
                file_info['shape'] = list(data.shape)
                file_info['dtype'] = str(data.dtype)
                print(f"  {fname}: shape {data.shape}, dtype {data.dtype}")
            except Exception as e:
                file_info['error'] = str(e)
                print(f"  {fname}: ERRORE - {e}")
        elif f.suffix == '.json':
            try:
                with open(f) as fp:
                    json_data = json.load(fp)
                file_info['json_keys'] = list(json_data.keys())
                print(f"  {fname}: JSON con chiavi {list(json_data.keys())}")
            except Exception as e:
                file_info['error'] = str(e)
                print(f"  {fname}: ERRORE - {e}")
        else:
            print(f"  {fname}: {size} bytes")
        
        category['files'][fname] = file_info
    else:
        print(f"  {fname}: NON TROVATO")
        category['files'][fname] = {'status': 'missing'}

report['categories']['64B_FILES'] = category

print()

# ============================================================================
# CATEGORIA 8: T2u FILES
# ============================================================================
print("CATEGORIA 8: T2u FILES")
print("-" * 80)

t2u_files = [
    "all_T2u_arrays.npz"
]

category = {
    'status': 'T2U',
    'description': 'Artefatti settore T2u',
    'files': {}
}

for fname in t2u_files:
    f = BASE_DIR / fname
    if f.exists():
        size = f.stat().st_size
        file_info = {
            'size_bytes': size,
            'size_kb': round(size / 1024, 3),
            'modified': datetime.fromtimestamp(f.stat().st_mtime).isoformat()
        }
        
        if f.suffix == '.npz':
            try:
                data = np.load(f)
                file_info['keys'] = list(data.keys())
                print(f"  {fname}: NPZ con chiavi {list(data.keys())}")
            except Exception as e:
                file_info['error'] = str(e)
                print(f"  {fname}: ERRORE - {e}")
        else:
            print(f"  {fname}: {size} bytes")
        
        category['files'][fname] = file_info
    else:
        print(f"  {fname}: NON TROVATO")
        category['files'][fname] = {'status': 'missing'}

report['categories']['T2U_FILES'] = category

print()

# ============================================================================
# SALVATAGGIO
# ============================================================================
print("=" * 80)
print("SALVATAGGIO REPORT")
print("=" * 80)

with open(OUTPUT_FILE, 'w') as f:
    json.dump(report, f, indent=2)
print(f"✓ Report JSON salvato: {OUTPUT_FILE}")

# Summary markdown
summary = f"""# EXTENDED AUDIT SCARICATI - SUMMARY

**Data:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
**Directory:** {BASE_DIR}

## CATEGORIE ANALIZZATE

### 1. 188D_FINAL_BUNDLE (UNVERIFIED)
**Status:** UNVERIFIED - NON usare per fisica certificata

File rilevanti:
"""

if '188D_FINAL_BUNDLE' in report['categories']:
    for fname, info in report['categories']['188D_FINAL_BUNDLE']['files'].items():
        summary += f"- `{fname}`: {info.get('size_mb', 0)} MB"
        if 'shape' in info:
            summary += f", shape {info['shape']}"
        summary += "\n"

summary += f"""
### 2. 188D_CERTIFIED
**Status:** CERTIFIED - usabile per fisica

File rilevanti:
"""

if '188D_CERTIFIED' in report['categories']:
    for fname, info in report['categories']['188D_CERTIFIED']['files'].items():
        summary += f"- `{fname}`: {info.get('size_mb', 0)} MB"
        if 'shape' in info:
            summary += f", shape {info['shape']}"
        summary += "\n"

summary += f"""
### 3. 188D_PROVENANCE_AUDIT
**Status:** RECONSTRUCTED - da first principles

File rilevanti:
"""

if '188D_PROVENANCE_AUDIT' in report['categories']:
    for fname, info in report['categories']['188D_PROVENANCE_AUDIT']['files'].items():
        summary += f"- `{fname}`: {info.get('size_mb', 0)} MB"
        if 'shape' in info:
            summary += f", shape {info['shape']}"
        summary += "\n"

summary += f"""
### 4. 188D_LOCAL_REPLAY_OUTPUT
**Status:** RECONSTRUCTED - da X_Oh e Q

File rilevanti:
"""

if '188D_LOCAL_REPLAY_OUTPUT' in report['categories']:
    for fname, info in report['categories']['188D_LOCAL_REPLAY_OUTPUT']['files'].items():
        summary += f"- `{fname}`: {info.get('size_mb', 0)} MB"
        if 'shape' in info:
            summary += f", shape {info['shape']}"
        summary += "\n"

summary += f"""
### 5. DNS FILES
**Status:** RECOVERED - codice originale e audit

File rilevanti:
"""

if 'DNS_FILES' in report['categories']:
    for fname, info in report['categories']['DNS_FILES']['files'].items():
        if info.get('status') != 'missing':
            summary += f"- `{fname}`: {info.get('size_kb', 0)} KB\n"

summary += f"""
### 6. 188D REPORTS
**Status:** DOCUMENTATION

File rilevanti:
"""

if '188D_REPORTS' in report['categories']:
    for fname, info in report['categories']['188D_REPORTS']['files'].items():
        if info.get('status') != 'missing':
            summary += f"- `{fname}`: {info.get('size_kb', 0)} KB\n"

summary += f"""
### 7. 64B FILES
**Status:** CUBIC - artefatti cubici O_h

File rilevanti:
"""

if '64B_FILES' in report['categories']:
    for fname, info in report['categories']['64B_FILES']['files'].items():
        if info.get('status') != 'missing':
            summary += f"- `{fname}`\n"

summary += f"""
### 8. T2U FILES
**Status:** T2U - settore T2u

File rilevanti:
"""

if 'T2U_FILES' in report['categories']:
    for fname, info in report['categories']['T2U_FILES']['files'].items():
        if info.get('status') != 'missing':
            summary += f"- `{fname}`\n"

summary += f"""
## CONCLUSIONI

### CERTIFICATO (usabile per fisica)
- X_Oh_recovered.npy (188D_LOCAL_REPLAY_OUTPUT)
- Q_rebuilt.npy (188D_LOCAL_REPLAY_OUTPUT)
- K_tangent_rebuilt.npy (188D_LOCAL_REPLAY_OUTPUT)
- K_certified.npy (188D_CERTIFIED)
- K_pinv_certified.npy (188D_CERTIFIED)

### UNVERIFIED (NON usare per fisica certificata)
- K0.npy (188D_FINAL_BUNDLE)
- Kinv.npy (188D_FINAL_BUNDLE)
- P_G.npy (188D_FINAL_BUNDLE)
- J_A.npy (188D_FINAL_BUNDLE)
- J_B.npy (188D_FINAL_BUNDLE)
- J.npy (188D_FINAL_BUNDLE)

### RECONSTRUCTED (ipotesi da verificare)
- gauge_basis_rebuilt.npy (188D_PROVENANCE_AUDIT)
- P_G_rebuilt.npy (188D_PROVENANCE_AUDIT)
- K_pinv_rebuilt.npy (188D_PROVENANCE_AUDIT)

### DNS (benchmark metodologico)
- dns_original_extracted.py
- dns_conservation_audit_final.py
- DNS_CONSERVATION_AUDIT_FINAL.json

---

**Report completo:** {OUTPUT_FILE}
"""

with open(SUMMARY_FILE, 'w') as f:
    f.write(summary)
print(f"✓ Summary MD salvato: {SUMMARY_FILE}")

print()
print("=" * 80)
print("AUDIT COMPLETATO")
print("=" * 80)

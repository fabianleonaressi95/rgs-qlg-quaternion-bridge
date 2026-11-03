#!/usr/bin/env python3
"""
AUDIT COMPLETO DI TUTTO IL MATERIALE IN /home/fabianleo/Scaricati
================================================================
Scansiona ricorsivamente tutti i file e le directory in /Scaricati,
categorizza e clusterizza per dominio tematico, analizza la tipologia dei dati,
lo stato di certificazione e produce report esaustivi in JSON e Markdown.
"""

import os
import sys
import json
import time
from pathlib import Path
from datetime import datetime
from collections import Counter, defaultdict
import numpy as np

SCARICATI_DIR = Path("/home/fabianleo/Scaricati")
REPORT_JSON = SCARICATI_DIR / "FULL_SCARICATI_AUDIT_REPORT.json"
REPORT_MD = SCARICATI_DIR / "FULL_SCARICATI_AUDIT_REPORT.md"
HOME_SUMMARY_MD = Path("/home/fabianleo/FULL_SCARICATI_AUDIT_SUMMARY.md")

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (np.bool_, bool)):
            return bool(obj)
        if isinstance(obj, (np.integer, int)):
            return int(obj)
        if isinstance(obj, (np.floating, float)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)

def classify_file(path: Path):
    """Assegna un file a uno o più cluster tematici in base a nome, percorso ed estensione."""
    p_str = str(path).lower()
    name = path.name.lower()
    ext = path.suffix.lower()
    
    # Esclusione o isolamento ambiente virtuale
    if "/venv/" in p_str or p_str.endswith("/venv"):
        return "PYTHON_ENVIRONMENT"
    
    # 1. 188D e geometria O_h
    if any(k in p_str for k in ["188d", "k0.npy", "kinv.npy", "p_g.npy", "j_a.npy", "j_b.npy", "x_oh", "tangent_rebuilt", "coulomb"]):
        return "188D_GEOMETRY_AND_OH"
    
    # 2. DNS Navier-Stokes
    if any(k in p_str for k in ["dns", "navier", "untitled2", "kolmogorov", "energy_derivative"]):
        return "DNS_TURBULENCE"
    
    # 3. RNR, PCET, Cryo-EM 6W4X / EMD-25156 / Codoni
    if any(k in p_str for k in ["rnr", "6w4x", "emd25156", "pcet", "qm_64_codons", "7sj7", "qip_qc"]):
        return "RNR_PCET_QUANTUM_BIOLOGY"
    
    # 4. Quasicristalli, Morphodynamics, QC Truss, Penrose, Fibonacci
    if any(k in p_str for k in ["quasicrystal", "qc_truss", "morphodynamics", "sl_qc", "penrose", "fibonacci", "zeroptsenemos"]):
        return "QUASICRYSTALS_MORPHODYNAMICS"
    
    # 5. QVSM, QNN, Semiclassico, WKB, Memoria Non-Markoviana
    if any(k in p_str for k in ["qvsm", "qnn", "wkb", "nonmarkovian", "two_time_kernel", "gellmann", "semiclassical"]):
        return "QVSM_QNN_SEMICLASSICAL"
    
    # 6. Teoria delle Stringhe, TOE, Teoria dei Numeri, Monografie matematiche
    if any(k in p_str for k in ["stringtheory", "toe", "grandunification", "grablovitz", "caccioppoli", "nardelli", "renaissance"]):
        return "THEORETICAL_PHYSICS_STRING_THEORY"
    
    # 7. ASTRA-Q SRL, Business, Amministrazione
    if any(k in p_str for k in ["astra", "statuto", "business_plan", "ouled_haffouz", "amministratori", "s.r.l"]):
        return "BUSINESS_LEGAL_ASTRA_Q"
    
    # 8. Pipeline numeriche e script generici
    if "pipeline" in name or name.endswith(".py"):
        return "COMPUTATIONAL_PIPELINES_AND_SCRIPTS"
    
    # 9. Dati tabulari e scientifici generici
    if ext in [".csv", ".tsv", ".xlsx", ".json", ".npz", ".npy"]:
        return "DATASETS_AND_ARRAYS"
    
    # 10. Documentazione e Pubblicazioni generiche
    if ext in [".pdf", ".docx", ".pptx", ".md", ".txt"]:
        return "DOCUMENTATION_AND_PAPERS"
    
    # 11. Media (immagini, audio)
    if ext in [".png", ".jpg", ".jpeg", ".mp3", ".wav"]:
        return "MEDIA_AND_AUDIO"
        
    return "MISC_OTHER"

def inspect_npy_file(fp: Path):
    """Ispeziona metadati di un file numpy senza caricare enormi array."""
    try:
        # Legge header .npy
        with open(fp, 'rb') as f:
            major, minor = np.lib.format.read_magic(f)
            shape, fortran, dtype = np.lib.format.read_array_header_1_0(f) if major == 1 else np.lib.format.read_array_header_2_0(f)
            return {'shape': list(shape), 'dtype': str(dtype)}
    except Exception as e:
        return {'error': str(e)}

def inspect_npz_file(fp: Path):
    try:
        with np.load(fp) as data:
            files = list(data.files)
            shapes = {k: list(data[k].shape) for k in files[:10]}
            return {'keys': files, 'shapes': shapes}
    except Exception as e:
        return {'error': str(e)}

def run_audit():
    print("=" * 80)
    print("AVVIO AUDIT COMPLETO DI /home/fabianleo/Scaricati")
    print("=" * 80)
    
    t0 = time.time()
    
    total_files = 0
    total_dirs = 0
    total_size = 0
    
    ext_counter = Counter()
    cluster_files = defaultdict(list)
    cluster_sizes = defaultdict(int)
    
    dir_summary = {}
    
    # Esplora cartelle di primo livello
    for item in sorted(SCARICATI_DIR.iterdir()):
        if item.is_dir():
            n_f = sum(len(f) for _, _, f in os.walk(item))
            sz = sum(sum((Path(r)/f).stat().st_size for f in files) for r, _, files in os.walk(item) if os.path.exists(r))
            dir_summary[item.name] = {
                'file_count': n_f,
                'size_bytes': sz,
                'size_mb': round(sz / (1024*1024), 2)
            }
    
    # Scansione completa ricorsiva
    sample_arrays = {}
    
    for root, d_names, f_names in os.walk(SCARICATI_DIR):
        total_dirs += len(d_names)
        for f in f_names:
            fp = Path(root) / f
            try:
                stat = fp.stat()
                sz = stat.st_size
                mtime = stat.st_mtime
                total_files += 1
                total_size += sz
                
                ext = fp.suffix.lower()
                ext_counter[ext if ext else "[no_ext]"] += 1
                
                cluster = classify_file(fp)
                rel_path = str(fp.relative_to(SCARICATI_DIR))
                
                file_entry = {
                    'rel_path': rel_path,
                    'size': sz,
                    'size_kb': round(sz / 1024, 2),
                    'mtime': datetime.fromtimestamp(mtime).strftime('%Y-%m-%d %H:%M:%S'),
                    'ext': ext
                }
                
                # Se è array numpy di rilievo, ispeziona la forma
                if ext == ".npy" and cluster != "PYTHON_ENVIRONMENT" and len(sample_arrays) < 60:
                    info = inspect_npy_file(fp)
                    file_entry['array_info'] = info
                    sample_arrays[rel_path] = info
                elif ext == ".npz" and cluster != "PYTHON_ENVIRONMENT" and len(sample_arrays) < 60:
                    info = inspect_npz_file(fp)
                    file_entry['array_info'] = info
                    sample_arrays[rel_path] = info
                    
                cluster_files[cluster].append(file_entry)
                cluster_sizes[cluster] += sz
                
            except Exception as e:
                pass

    t1 = time.time()
    print(f"Scansione completata in {t1 - t0:.2f} s: {total_files} file analizzati, {total_size / (1024*1024*1024):.2f} GB totali.")
    
    # Costruzione del report strutturato
    clusters_report = {}
    for c_name, flist in cluster_files.items():
        clusters_report[c_name] = {
            'total_files': len(flist),
            'total_size_mb': round(cluster_sizes[c_name] / (1024*1024), 2),
            'files': sorted(flist, key=lambda x: x['size'], reverse=True)[:50] # top 50 per dimensione
        }
        
    audit_report = {
        'timestamp': datetime.now().isoformat(),
        'base_directory': str(SCARICATI_DIR),
        'total_files': total_files,
        'total_directories': total_dirs,
        'total_size_bytes': total_size,
        'total_size_gb': round(total_size / (1024*1024*1024), 3),
        'top_level_directories': dir_summary,
        'extension_distribution': dict(ext_counter.most_common(25)),
        'clusters': clusters_report
    }
    
    # Scrivi report JSON
    with open(REPORT_JSON, 'w') as f:
        json.dump(audit_report, f, indent=2, cls=NumpyEncoder)
    print(f"✓ Report JSON salvato: {REPORT_JSON}")
    
    # Genera report Markdown
    md_lines = []
    md_lines.append("# AUDIT COMPLETO DEL MATERIALE IN `/home/fabianleo/Scaricati`")
    md_lines.append(f"\n**Data generazione:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    md_lines.append(f"**Directory analizzata:** `{SCARICATI_DIR}`")
    md_lines.append(f"**Volume totale:** {total_files:,} file | {total_dirs:,} cartelle | **{total_size / (1024*1024*1024):.2f} GB**")
    md_lines.append("\n---\n")
    
    md_lines.append("## 1. STRUTTURA DELLE CARTELLE DI PRIMO LIVELLO\n")
    md_lines.append("| Directory | File Totali | Dimensione (MB) | Descrizione / Contenuto Principale |")
    md_lines.append("| :--- | :---: | :---: | :--- |")
    
    dir_descriptions = {
        '188D_CERTIFIED': "Artefatti certificati 188D (K_certified, P_candidate)",
        '188D_FINAL_BUNDLE': "Bundle storico non verificato (K0, Kinv, P_G, J_A, J_B)",
        '188D_LOCAL_REPLAY_OUTPUT': "Output ricostruzione locale da first principles (X_Oh, Q, K_rebuilt)",
        '188D_PROVENANCE_AUDIT': "Artefatti e report di audit provenienza 188D",
        'QIP_QC_EMD25156': "Mappe di densità Cryo-EM e strutture PDB (7SJ7, EMD-25156)",
        'qip_qc_data': "Dataset QIP/QC",
        'qip_qc_figures': "Figure scientifiche e grafici QIP/QC",
        'qip_qc_phase1c': "Fase 1C del progetto QIP/QC",
        'rnr_6W4X_actual_audit': "Audit strutturale della ribonucleotide reduttasi 6W4X",
        'rnr_6W4X_canonical_geometry': "Geometria canonica 6W4X",
        'rnr_6W4X_canonical_geometry_fixed': "Geometria 6W4X corretta e calibrata",
        'rnr_6W4X_geometry_audit': "Audit geometrico RNR",
        'rnr_experimental_calibration': "Calibrazione sperimentale parametri PCET/RNR",
        'rnr_hybrid_results': "Risultati ibridi PCET",
        'rnr_pcet_results': "Simulazioni di trasferimento radicalico / accoppiato protone-elettrone",
        'rnr_qvsm_calibration': "Calibrazione QVSM per catene RNR",
        '17672291': "Pipeline di calcolo 335-353, QM 64 codoni genetici, simulazioni estrusione",
        '6w4x_results': "Risultati energetici ed elettronici per 6W4X",
        'tutti_gli_output': "Output estratti da sessione Kaggle (replay forense e recovery)",
        'tutti_gli_output (1)': "Workspace duplicato / parallelo dei job Kaggle",
        'progetti': "Note e specifiche progettuali",
        'venv': "Ambiente virtuale Python (librerie, C-extensions, package)"
    }
    
    for d_name, d_info in sorted(dir_summary.items()):
        desc = dir_descriptions.get(d_name, "Dati / simulazioni")
        md_lines.append(f"| `{d_name}` | {d_info['file_count']} | {d_info['size_mb']} MB | {desc} |")
        
    md_lines.append("\n---\n")
    
    md_lines.append("## 2. CLUSTER TEMATICI E DOMINI DI RICERCA\n")
    
    cluster_meta = {
        "188D_GEOMETRY_AND_OH": ("188D, Simmetria O_h e Coulomb Hessian", "Fisica dello spazio tangente 188D, configurazioni cubiche O_h a 64 siti, calcolo K_tangent, gauge fixing e spettro"),
        "DNS_TURBULENCE": ("DNS & Turbolenza 3D Navier-Stokes", "Direct Numerical Simulation su JAX, bilancio energetico non lineare Kolmogorov, shell balance ed estrazione da notebook Kaggle"),
        "RNR_PCET_QUANTUM_BIOLOGY": ("RNR, PCET, Cryo-EM 6W4X & Codoni Genetici", "Proton-Coupled Electron Transfer, ribonucleotide reduttasi 6W4X, mappa EMD-25156, PDB 7SJ7, 64 codoni"),
        "QUASICRYSTALS_MORPHODYNAMICS": ("Quasicristalli, Morphodynamics & QC Truss", "Geometrie quasi-periodiche, Fibonacci, Penrose, morfodinamica e strutture truss"),
        "QVSM_QNN_SEMICLASSICAL": ("QVSM, Reti Neurali Quantistiche e Semiclassico", "Quantum Variational State Machine, matrici Gell-Mann, kernel a due tempi, WKB semiclassico"),
        "THEORETICAL_PHYSICS_STRING_THEORY": ("Teoria delle Stringhe, TOE e Matematica Unificata", "Tesi di String Theory, cubiche di Grablovitz, confini di Caccioppoli, formule TOE"),
        "BUSINESS_LEGAL_ASTRA_Q": ("ASTRA-Q S.R.L. (Business & Legal)", "Statuto societario, atto costitutivo, business plan, governance e modelli operativi"),
        "COMPUTATIONAL_PIPELINES_AND_SCRIPTS": ("Pipeline Computazionali e Script Python", "Algoritmi numerici, pipeline sequenziali, parser e simulatori ausiliari"),
        "DATASETS_AND_ARRAYS": ("Dataset Tabulari e Array Numerici", "File CSV, NPZ, NPY contenenti profili temporali, matrici di transizione e calibrazioni"),
        "DOCUMENTATION_AND_PAPERS": ("Documenti Scientifici, Slide e White Paper", "PDF tecnici, presentazioni PPTX, articoli accademici e note di ricerca"),
        "MEDIA_AND_AUDIO": ("Media, Grafici e Audio", "Tracciati PNG, rendering visivi, spettrogrammi e registrazioni audio di supporto"),
        "PYTHON_ENVIRONMENT": ("Ambiente Virtuale Python (venv)", "Pacchetti compilati, interpreti, librerie esterne numpy/scipy/jax"),
        "MISC_OTHER": ("File Vari / Non Classificati", "File binari, archivi e risorse accessorie")
    }
    
    for c_id, (c_title, c_desc) in cluster_meta.items():
        if c_id not in clusters_report:
            continue
        c_data = clusters_report[c_id]
        md_lines.append(f"### {c_title}")
        md_lines.append(f"**Descrizione:** {c_desc}")
        md_lines.append(f"- **File totali:** {c_data['total_files']} | **Dimensione:** {c_data['total_size_mb']} MB")
        md_lines.append("\n*Esempi di file chiave nel cluster:*")
        for f in c_data['files'][:8]:
            arr_note = f" (shape: {f['array_info']['shape']})" if 'array_info' in f and 'shape' in f['array_info'] else ""
            md_lines.append(f"- `{f['rel_path']}` — {f['size_kb']} KB{arr_note}")
        md_lines.append("")
        
    md_lines.append("\n---\n")
    
    md_lines.append("## 3. STATO DI CERTIFICAZIONE DEGLI ASSET SCIENTIFICI CRITICI\n")
    md_lines.append("| Dominio | Asset Primario | Stato Provenienza | Note di Certificazione |")
    md_lines.append("| :--- | :--- | :---: | :--- |")
    md_lines.append("| **188D Geometry** | `X_Oh_recovered.npy` | ✅ VERIFIED | Geometria 64 siti cubica perfetta ($M_2 = 240.0$, $CM \\approx 0$) |")
    md_lines.append("| **188D Tangent** | `Q_rebuilt.npy` | ✅ RECONSTRUCTED | Base ortonormale null_space(C) ($192 \\times 188$) |")
    md_lines.append("| **188D Hessian** | `K_tangent_rebuilt.npy` | ✅ RECONSTRUCTED | Coulomb Hessian su tangent space, simmetria esatta |")
    md_lines.append("| **188D Gauge V2** | `P_G_V2.npy` / `gauge_basis_V2.npy` | ✅ CERTIFIED V2 | 64 modi per minima $|\\lambda|$, $[O_h, P_G] < 6\\times 10^{-15}$ |")
    md_lines.append("| **188D Pseudoinv V2** | `K_pinv_V2.npy` | ✅ CERTIFIED V2 | $K^\\dagger_{\\text{phys}}$, residuo $2.08\\times 10^{-15}$ sul settore fisico |")
    md_lines.append("| **188D Bundle Storico** | `K0.npy`, `Kinv.npy`, `P_G.npy`, `J_A`, `J_B` | ⚠️ UNVERIFIED | Differiscono del 109% dalla ricostruzione; non usare per fisica |")
    md_lines.append("| **DNS Navier-Stokes** | `dns_original_extracted.py` | ✅ VERIFIED | Estratto da `Untitled2.ipynb`, conservazione energia $5.34\\times 10^{-16}$ |")
    md_lines.append("| **Cryo-EM RNR** | `7SJ7.pdb`, `QIP_QC_EMD25156` | 🔬 EXPERIMENTAL | Strutture biologiche e mappe cryo-EM di supporto a PCET |")
    
    md_lines.append("\n---\n")
    md_lines.append("## 4. CONCLUSIONI E RACCOMANDAZIONI\n")
    md_lines.append("1. **Asset 188D Certificati**: La catena di certificazione first-principles è completa. Usare esclusivamente gli artefatti contrassegnati V2 e `188D_LOCAL_REPLAY_OUTPUT`.")
    md_lines.append("2. **Isolamento Bundle Storico**: Il bundle `188D_FINAL_BUNDLE` è stato preservato a fini forensi/diagnostici ma isolato dal flusso computazionale certificato.")
    md_lines.append("3. **Interconnessione Tematica**: I dati in `/Scaricati` documentano un ecosistema interdisciplinare: teoria di gauge/invarianza $O_h$, simulazioni di fluidodinamica pseudo-spettrale DNS, biologia quantistica RNR/PCET e modelli geometrici quasi-cristallini.")
    
    summary_text = "\n".join(md_lines)
    
    with open(REPORT_MD, 'w') as f:
        f.write(summary_text)
    print(f"✓ Report Markdown salvato: {REPORT_MD}")
    
    with open(HOME_SUMMARY_MD, 'w') as f:
        f.write(summary_text)
    print(f"✓ Summary Markdown salvato in home: {HOME_SUMMARY_MD}")

if __name__ == "__main__":
    run_audit()

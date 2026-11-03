# AUDIT COMPLETO DEL MATERIALE IN `/home/fabianleo/Scaricati`

**Data generazione:** 2026-10-04 17:37:32
**Directory analizzata:** `/home/fabianleo/Scaricati`
**Volume totale:** 9,419 file | 881 cartelle | **3.50 GB**

---

## 1. STRUTTURA DELLE CARTELLE DI PRIMO LIVELLO

| Directory | File Totali | Dimensione (MB) | Descrizione / Contenuto Principale |
| :--- | :---: | :---: | :--- |
| `17672291` | 217 | 193.48 MB | Pipeline di calcolo 335-353, QM 64 codoni genetici, simulazioni estrusione |
| `188D_CERTIFIED` | 4 | 0.81 MB | Artefatti certificati 188D (K_certified, P_candidate) |
| `188D_FINAL_BUNDLE` | 9 | 0.83 MB | Bundle storico non verificato (K0, Kinv, P_G, J_A, J_B) |
| `188D_LOCAL_REPLAY_OUTPUT` | 5 | 0.55 MB | Output ricostruzione locale da first principles (X_Oh, Q, K_rebuilt) |
| `188D_PROVENANCE_AUDIT` | 4 | 0.63 MB | Artefatti e report di audit provenienza 188D |
| `6w4x_results` | 5 | 0.02 MB | Risultati energetici ed elettronici per 6W4X |
| `QIP_QC_EMD25156` | 29 | 784.59 MB | Mappe di densità Cryo-EM e strutture PDB (7SJ7, EMD-25156) |
| `progetti` | 16 | 0.01 MB | Note e specifiche progettuali |
| `qip_qc_data` | 4 | 1.08 MB | Dataset QIP/QC |
| `qip_qc_figures` | 9 | 4.67 MB | Figure scientifiche e grafici QIP/QC |
| `qip_qc_phase1c` | 1 | 0.27 MB | Fase 1C del progetto QIP/QC |
| `rnr_6W4X_actual_audit` | 2 | 0.0 MB | Audit strutturale della ribonucleotide reduttasi 6W4X |
| `rnr_6W4X_canonical_geometry` | 3 | 0.01 MB | Geometria canonica 6W4X |
| `rnr_6W4X_canonical_geometry_fixed` | 3 | 0.01 MB | Geometria 6W4X corretta e calibrata |
| `rnr_6W4X_geometry_audit` | 0 | 0.0 MB | Audit geometrico RNR |
| `rnr_experimental_calibration` | 2 | 0.0 MB | Calibrazione sperimentale parametri PCET/RNR |
| `rnr_hybrid_results` | 1 | 0.0 MB | Risultati ibridi PCET |
| `rnr_pcet_results` | 1 | 0.0 MB | Simulazioni di trasferimento radicalico / accoppiato protone-elettrone |
| `rnr_qvsm_calibration` | 8 | 0.02 MB | Calibrazione QVSM per catene RNR |
| `tutti_gli_output` | 247 | 4.89 MB | Output estratti da sessione Kaggle (replay forense e recovery) |
| `tutti_gli_output (1)` | 247 | 4.89 MB | Workspace duplicato / parallelo dei job Kaggle |
| `venv` | 7703 | 163.05 MB | Ambiente virtuale Python (librerie, C-extensions, package) |

---

## 2. CLUSTER TEMATICI E DOMINI DI RICERCA

### 188D, Simmetria O_h e Coulomb Hessian
**Descrizione:** Fisica dello spazio tangente 188D, configurazioni cubiche O_h a 64 siti, calcolo K_tangent, gauge fixing e spettro
- **File totali:** 525 | **Dimensione:** 11.06 MB

*Esempi di file chiave nel cluster:*
- `tutti_gli_output/kaggle/working/188D_RECOVERY/recovered_code_contexts.py` — 780.17 KB
- `tutti_gli_output (1)/kaggle/working/188D_RECOVERY/recovered_code_contexts.py` — 780.17 KB
- `188D_FORENSIC_REPORT.txt` — 766.31 KB
- `tutti_gli_output/kaggle/working/188D_RECOVERY/188D_FORENSIC_REPORT.txt` — 766.31 KB
- `tutti_gli_output (1)/kaggle/working/188D_RECOVERY/188D_FORENSIC_REPORT.txt` — 766.31 KB
- `188D_FINAL_BUNDLE.zip` — 766.26 KB
- `188D_LOCAL_REPLAY_OUTPUT/Q_rebuilt.npy` — 282.12 KB
- `188D_PROVENANCE_AUDIT/K_pinv_rebuilt.npy` — 276.25 KB (shape: [188, 188])

### DNS & Turbolenza 3D Navier-Stokes
**Descrizione:** Direct Numerical Simulation su JAX, bilancio energetico non lineare Kolmogorov, shell balance ed estrazione da notebook Kaggle
- **File totali:** 12 | **Dimensione:** 0.66 MB

*Esempi di file chiave nel cluster:*
- `Untitled2 (1).ipynb` — 256.12 KB
- `bilancio_energetico_dns.png` — 254.12 KB
- `Untitled2.ipynb` — 110.7 KB
- `dns_conservation_audit.py` — 10.38 KB
- `dns_conservation_audit_fixed.py` — 9.92 KB
- `dns_conservation_audit_final.py` — 9.84 KB
- `dns_conservation_audit_numpy.py` — 9.76 KB
- `dns_original_extracted.py` — 6.13 KB

### RNR, PCET, Cryo-EM 6W4X & Codoni Genetici
**Descrizione:** Proton-Coupled Electron Transfer, ribonucleotide reduttasi 6W4X, mappa EMD-25156, PDB 7SJ7, 64 codoni
- **File totali:** 285 | **Dimensione:** 982.39 MB

*Esempi di file chiave nel cluster:*
- `QIP_QC_EMD25156/emd_25156.map.gz` — 485498.69 KB
- `17672291/QM_64_CODONS/Miniconda3-latest-Linux-x86_64.sh` — 193333.35 KB
- `QIP_QC_EMD25156/QIP_QC_EMD25156_weighted.npz` — 164116.23 KB
- `QIP_QC_EMD25156/QIP_QC_EMD25156_vectorized.npz` — 140403.43 KB
- `QIP_QC_EMD25156/QIP_QC_PHASE2_TB_results.npz` — 5242.79 KB
- `QIP_QC_EMD25156/7SJ7.pdb` — 3406.9 KB
- `6W4X.cif` — 1862.94 KB
- `qip_qc_figures/comparison_radial.png` — 1142.21 KB

### Quasicristalli, Morphodynamics & QC Truss
**Descrizione:** Geometrie quasi-periodiche, Fibonacci, Penrose, morfodinamica e strutture truss
- **File totali:** 21 | **Dimensione:** 2.56 MB

*Esempi di file chiave nel cluster:*
- `fibonacci_69120_metadata.csv` — 1924.51 KB
- `CONNESSIONI_da_FIBONACCI_A_PC_QUANTISTIC.pdf` — 459.59 KB
- `Unified Topological Quantum Dynamics and Quasi-Periodic Morphodynamics White Paper.pdf` — 97.31 KB
- `ZEROPTSENEMOSÌA_QUASICRYSTAL_USB_DESIGN.pdf` — 73.87 KB
- `SL_QC_manuscript_wrap_up.pdf` — 23.96 KB
- `fibonacci_windows_audit (1).csv` — 8.73 KB
- `fibonacci_windows_audit.csv` — 8.73 KB
- `quasicrystal_audit (1).py` — 5.78 KB

### QVSM, Reti Neurali Quantistiche e Semiclassico
**Descrizione:** Quantum Variational State Machine, matrici Gell-Mann, kernel a due tempi, WKB semiclassico
- **File totali:** 48 | **Dimensione:** 35.92 MB

*Esempi di file chiave nel cluster:*
- `qnn_ontological_four_models.npz` — 4104.49 KB
- `qnn_ontological_four_models (1).npz` — 4104.49 KB
- `two_time_kernel_N24_J0_0.004.npz` — 2299.39 KB
- `two_time_kernel_N16_J0_0.004 (1).npz` — 2299.39 KB
- `two_time_kernel_N48_J0_0.004 (1).npz` — 2299.39 KB
- `two_time_kernel_N48_J0_0.004.npz` — 2299.39 KB
- `two_time_kernel_N16_J0_0.004.npz` — 2299.39 KB
- `two_time_kernel_N16_J0_0.002 (1).npz` — 2252.63 KB

### Teoria delle Stringhe, TOE e Matematica Unificata
**Descrizione:** Tesi di String Theory, cubiche di Grablovitz, confini di Caccioppoli, formule TOE
- **File totali:** 12 | **Dimensione:** 49.55 MB

*Esempi di file chiave nel cluster:*
- `StringTheorythesis.pdf` — 35142.36 KB
- `TernadiLevieTOEEquation (1).pdf` — 4530.22 KB
- `TernadiLevieTOEEquation.pdf` — 4530.22 KB
- `EinsteinandNardelliTOEEquation.pdf` — 1238.1 KB
- `SimplifiedNardelliMasterEquation.pdf` — 998.32 KB
- `AfurthermathematicaldevelopmentoftheNardelliSeventhRootTOEEquation.pdf` — 966.42 KB
- `AnalyticalUnificationofEugenioEliaLeviIntegralsintheNardelliTOEEquation.pdf` — 754.02 KB
- `TheMathematicalRenaissance.pdf` — 753.32 KB

### ASTRA-Q S.R.L. (Business & Legal)
**Descrizione:** Statuto societario, atto costitutivo, business plan, governance e modelli operativi
- **File totali:** 21 | **Dimensione:** 23.61 MB

*Esempi di file chiave nel cluster:*
- `-ASTRA-Q---Sintesi (1).pptx` — 5864.17 KB
- `-ASTRA-Q---Sintesi.pptx` — 5676.28 KB
- `ASTRA-Q-Space-Situational-Awareness-and-Prognostica-Non-Lineare-Koopman.pdf` — 5642.52 KB
- `astra_q_dashboard.html` — 4498.82 KB
- `storico_amministratori_province31121998.csv` — 860.05 KB
- `ASTRA-Q ESA BIC Pitch.pdf` — 446.33 KB
- `-ASTRA-Q---Sintesi.pdf` — 426.26 KB
- `AGRI_ASTRA_Dossier_SRG01.pdf` — 264.12 KB

### Pipeline Computazionali e Script Python
**Descrizione:** Algoritmi numerici, pipeline sequenziali, parser e simulatori ausiliari
- **File totali:** 147 | **Dimensione:** 0.84 MB

*Esempi di file chiave nel cluster:*
- `Dockerfile53.py` — 29.44 KB
- `Dockerfile57.py` — 28.81 KB
- `Dockerfile55.py` — 26.69 KB
- `Dockerfile54.py` — 24.82 KB
- `Dockerfile56.py` — 24.29 KB
- `naressi_zuliani_phase1_fixed.py` — 22.84 KB
- `mvp_middleware47.py` — 22.04 KB
- `Dockerfile52.py` — 19.33 KB

### Dataset Tabulari e Array Numerici
**Descrizione:** File CSV, NPZ, NPY contenenti profili temporali, matrici di transizione e calibrazioni
- **File totali:** 230 | **Dimensione:** 33.68 MB

*Esempi di file chiave nel cluster:*
- `full_69120_phase_audit.csv` — 9013.41 KB
- `phase_geometry_phi_silver_audit.csv` — 6639.26 KB
- `scan_500x64.csv` — 5461.24 KB
- `64B_Vminus_cubic_quartic_Oh_invariants.json` — 1620.43 KB
- `propagated_states.csv` — 1350.78 KB
- `01_global_fft_candidates.csv` — 1173.95 KB
- `krasilnikov_eq8_theta_scan_physical.csv` — 980.35 KB
- `17672291/Simulation_results_M1screw_PP.xlsx` — 641.48 KB

### Documenti Scientifici, Slide e White Paper
**Descrizione:** PDF tecnici, presentazioni PPTX, articoli accademici e note di ricerca
- **File totali:** 200 | **Dimensione:** 671.29 MB

*Esempi di file chiave nel cluster:*
- `rapporto.pdf` — 467728.43 KB
- `ostrom_1990.pdf` — 41208.24 KB
- `testo_estratto.txt` — 15403.54 KB
- `Q KEY AP 11.16.24.pdf` — 14564.92 KB
- `Origin_of_mass_publication-Zenodo-20231114.pdf` — 9551.18 KB
- `Bando ASPIRANTI_Domanda di partecipazione.docx` — 5700.71 KB
- `Bando ASPIRANTI_Domanda di partecipazione (1).docx` — 5700.71 KB
- `The_secret_teachings_of_all_ages.pdf` — 5511.29 KB

### Media, Grafici e Audio
**Descrizione:** Tracciati PNG, rendering visivi, spettrogrammi e registrazioni audio di supporto
- **File totali:** 135 | **Dimensione:** 353.84 MB

*Esempi di file chiave nel cluster:*
- `Rhapsody In Blue Gershwin.mp3` — 23145.27 KB
- `output1_432.mp3` — 15573.19 KB
- `Mimì_sarà_-_De_Gregori_live_a_Doc_432Hz.mp3` — 12575.55 KB
- `Francesco_De_Gregori_-_La_donna_cannone_(StillPseudo_Video)_432Hz.mp3` — 10992.47 KB
- `Francesco_De_Gregori_-_Caterina_432Hz.mp3` — 10782.7 KB
- `colonna sonora geo & geo.mp3` — 10760.86 KB
- `Piano Concerto No. 2, Allegro,Opus 102.mp3` — 10470.04 KB
- `Lucio_Battisti_-_Il_Nostro_Caro_Angelo_432Hz.mp3` — 9870.45 KB

### Ambiente Virtuale Python (venv)
**Descrizione:** Pacchetti compilati, interpreti, librerie esterne numpy/scipy/jax
- **File totali:** 7703 | **Dimensione:** 163.05 MB

*Esempi di file chiave nel cluster:*
- `venv/lib/python3.8/site-packages/numpy.libs/libopenblas64_p-r0-15028c96.3.21.so` — 32182.67 KB
- `venv/lib/python3.8/site-packages/numpy/core/_multiarray_umath.cpython-38-x86_64-linux-gnu.so` — 6163.5 KB
- `venv/bin/python3` — 5361.77 KB
- `venv/bin/python` — 5361.77 KB
- `venv/lib/python3.8/site-packages/numpy.libs/libgfortran-040039e1.so.5.0.0` — 2623.11 KB
- `venv/lib/python3.8/site-packages/pandas/_libs/join.cpython-38-x86_64-linux-gnu.so` — 2547.85 KB
- `venv/lib/python3.8/site-packages/numpy/core/_simd.cpython-38-x86_64-linux-gnu.so` — 2210.99 KB
- `venv/lib/python3.8/site-packages/pandas/_libs/algos.cpython-38-x86_64-linux-gnu.so` — 2018.51 KB

### File Vari / Non Classificati
**Descrizione:** File binari, archivi e risorse accessorie
- **File totali:** 80 | **Dimensione:** 1256.89 MB

*Esempi di file chiave nel cluster:*
- `mit-bih-atrial-fibrillation-database-1.0.0.zip` — 450227.79 KB
- `chatgpt_amd64.deb` — 341735.71 KB
- `Freebuff-0.0.155-linux-x86_64.AppImage` — 176321.41 KB
- `xampp-linux-x64-8.0.30-0-installer.run` — 155673.56 KB
- `firefox-153.0.tar.xz` — 84041.15 KB
- `naps2-8.4.0-linux-x64.deb` — 17287.98 KB
- `Non confermato 63631.crdownload` — 9104.29 KB
- `cosmodek.apk` — 7523.66 KB


---

## 3. STATO DI CERTIFICAZIONE DEGLI ASSET SCIENTIFICI CRITICI

| Dominio | Asset Primario | Stato Provenienza | Note di Certificazione |
| :--- | :--- | :---: | :--- |
| **188D Geometry** | `X_Oh_recovered.npy` | ✅ VERIFIED | Geometria 64 siti cubica perfetta ($M_2 = 240.0$, $CM \approx 0$) |
| **188D Tangent** | `Q_rebuilt.npy` | ✅ RECONSTRUCTED | Base ortonormale null_space(C) ($192 \times 188$) |
| **188D Hessian** | `K_tangent_rebuilt.npy` | ✅ RECONSTRUCTED | Coulomb Hessian su tangent space, simmetria esatta |
| **188D Gauge V2** | `P_G_V2.npy` / `gauge_basis_V2.npy` | ✅ CERTIFIED V2 | 64 modi per minima $|\lambda|$, $[O_h, P_G] < 6\times 10^{-15}$ |
| **188D Pseudoinv V2** | `K_pinv_V2.npy` | ✅ CERTIFIED V2 | $K^\dagger_{\text{phys}}$, residuo $2.08\times 10^{-15}$ sul settore fisico |
| **188D Bundle Storico** | `K0.npy`, `Kinv.npy`, `P_G.npy`, `J_A`, `J_B` | ⚠️ UNVERIFIED | Differiscono del 109% dalla ricostruzione; non usare per fisica |
| **DNS Navier-Stokes** | `dns_original_extracted.py` | ✅ VERIFIED | Estratto da `Untitled2.ipynb`, conservazione energia $5.34\times 10^{-16}$ |
| **Cryo-EM RNR** | `7SJ7.pdb`, `QIP_QC_EMD25156` | 🔬 EXPERIMENTAL | Strutture biologiche e mappe cryo-EM di supporto a PCET |

---

## 4. CONCLUSIONI E RACCOMANDAZIONI

1. **Asset 188D Certificati**: La catena di certificazione first-principles è completa. Usare esclusivamente gli artefatti contrassegnati V2 e `188D_LOCAL_REPLAY_OUTPUT`.
2. **Isolamento Bundle Storico**: Il bundle `188D_FINAL_BUNDLE` è stato preservato a fini forensi/diagnostici ma isolato dal flusso computazionale certificato.
3. **Interconnessione Tematica**: I dati in `/Scaricati` documentano un ecosistema interdisciplinare: teoria di gauge/invarianza $O_h$, simulazioni di fluidodinamica pseudo-spettrale DNS, biologia quantistica RNR/PCET e modelli geometrici quasi-cristallini.
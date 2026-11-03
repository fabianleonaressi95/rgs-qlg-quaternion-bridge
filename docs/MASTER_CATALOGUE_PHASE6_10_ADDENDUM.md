# Registro stabile delle fasi Q6–Dirac–Wilson–Overlap (Phase 6–10)

**Versione 1.0 — 3 novembre 2026**  
**Relazione:** appendice tecnica al [Catalogo generale](MASTER_CATALOGUE_LAWS_THEOREMS.md).

Questa appendice assegna identificatori stabili ai risultati già documentati nelle fasi 6–10. Gli stati si riferiscono esattamente alle ipotesi e agli oggetti finiti descritti; non trasferiscono automaticamente i risultati a una teoria fisica o a un limite continuo.

## Convenzioni epistemiche

- **PROVATO ALGEBRICAMENTE:** segue dalle definizioni e ipotesi dichiarate.
- **IDENTITÀ ESATTA:** combinatoria o algebra finita.
- **VERIFICATO NUMERICAMENTE:** residuo calcolato su matrici/dati specifici.
- **PROTOTIPO DIMOSTRATO:** prova analitica per una famiglia modello, non per RGS-Q6.
- **CANDIDATO:** operatore/struttura definita, ma proprietà ulteriori non dimostrate.
- **APERTO:** manca la dimostrazione o la famiglia di raffinamento richiesta.

## Registro degli enunciati

### Q6-CX-001 — Complesso cellulare Q6
**Enunciato.** Per il cubo a sei dimensioni: 64 vertici, 192 spigoli, grado 6 e 240 quadrati elementari. Per le matrici d'incidenza orientate del complesso, (B_1B_2=0).

**Stato:** IDENTITÀ ESATTA / PROVATO. Il conteggio delle facce è (\binom62 2^4=240). Il rango ciclico del solo grafo 1-dimensionale connesso è (192-64+1=129); non va confuso con l'omologia del complesso dopo l'aggiunta delle facce.

**Verifica:** `src/q6_complex.py`, test `tests/test_q6_complex.py`.  
**Limite:** il conteggio finito non definisce una successione di mesh né un limite continuo.

### Q6-SP-002 — Spettro del Laplaciano ipercubico
**Enunciato.** Per il Laplaciano combinatorio del grafo (Q_6), gli autovalori sono (lambda_j=2j), con molteplicità (inom6j), (j=0,ldots,6).

**Stato:** PROVATO / STANDARD.  
**Verifica:** test Q6 e costruzione del Laplaciano.  
**Limite:** vale per il grafo ipercubico canonico, non per un operatore geometrico RGS arbitrario.

### RGS-GEO-003 — Clifford curvo da coframe
**Enunciato.** Dato un coframe invertibile (e^A{}_mu(i)), definendo (Gamma^mu(i)=e_A{}^mu(i)gamma^A), la relazione di Clifford curva segue algebricamente dalla relazione piatta e da (g^{mu
u}=e_A{}^mu e_B{}^
ueta^{AB}).

**Stato:** PROVATO ALGEBRICAMENTE, sotto invertibilità del coframe e convenzioni di firma coerenti.  
**Verifica archivio Q6:** residuo massimo (3.55	imes10^{-15}); 64 siti, coframe spaziale (64	imes3	imes3).  
**Fonte:** `docs/WEILER_PHASE6_THEOREM_MAP.md`; `results/json/rgs_q6_phase6_reference_audit.json`.  
**Limite:** non deriva dinamicamente il tetrade, la metrica o le equazioni gravitazionali.

### RGS-SPIN-004 — Lift polare finito SO(3) → SU(2)
**Costruzione.** Per ogni spigolo canonico (i<j), si forma (M_{ji}=E_jE_i^{-1}), si estrae il fattore rotazionale polare (R_{ji}in SO(3)) e si sceglie un lift (S_{ji}in SU(2)).

**Stato:** CANDIDATO come trasporto geometrico discreto; proprietà di gruppo e residui NUMERICAMENTE VERIFICATI per l'archivio.  
**Residui:** ortogonalità SO(3) (3.00	imes10^{-15}); unitarietà SU(2) (3.14	imes10^{-16}); ricostruzione rotazionale (1.45	imes10^{-15}).  
**Fonte:** `results/json/rgs_q6_phase6_reference_audit.json`.  
**Limite:** il lift polare non è stato identificato con la connessione spin di Levi-Civita; la coerenza globale dei segni dei lift e delle olonomie va trattata esplicitamente.

### RGS-DIRAC-005 — Operatore spaziale di Dirac su Q6
**Definizione.** Con (c_i(n)=n_agamma_i^a), trasporto (S_{ji}) e (d_{ij}=x_j-x_i), si definisce
[
H_{ji}=rac{i}{2|d_{ij}|}left[c_j(n_j)S_{ji}+S_{ji}c_i(n_i)ight],
quad H_{ij}=H_{ji}^{dagger},
quad D=V^{-1/2}HV^{-1/2}.
]

**Stato:** DEFINIZIONE DI OPERATORE CANDIDATO; autoaggiuntezza finito-dimensionale PROVATA per costruzione con (V_i=det E_i>0); residui NUMERICAMENTE VERIFICATI.  
**Dati:** dimensione 256; hermiticità (3.25	imes10^{-15}); pairing spettrale (1.35	imes10^{-14}); nessun autovalore entro (10^{-9}) da zero nell'audit riportato.  
**Fonte:** `results/json/rgs_q6_phase6_reference_audit.json`.  
**Limite:** non prova consistenza su campi lisci, assenza di doubling, convergenza spettrale o limite di Dirac continuo.

### RGS-WILSON-006 — Laplaciano di Wilson covariante su grafo
**Enunciato.** Per link unitari (S_{ji}), pesi (w_{ij}=r/(2ell_{ij})ge0) e blocchi diagonali/off-diagonali costruiti come nel Phase 8, la forma quadratica è somma di termini (w_{ij}|psi_j-S_{ji}psi_i|^2). Ne seguono positività semidefinita e hermiticità; sotto (S_{ji}mapsto G_jS_{ji}G_i^dagger), il Laplaciano trasforma per coniugio a blocchi.

**Stato:** PROVATO ALGEBRICAMENTE per link unitari e pesi non negativi; NUMERICAMENTE VERIFICATO per l'archivio RGS-Q6.  
**Dati:** hermiticità (1.39	imes10^{-15}); minimo autovalore (5.22	imes10^{-4}); covarianza locale (6.86	imes10^{-15}); unitarietà link (7.02	imes10^{-16}).  
**Fonte:** `docs/CHIRAL_WILSON_AUDIT.md`, `results/json/rgs_q6_wilson_audit.json`.  
**Limite:** non dimostra assenza di doublers sul grafo Q6.

### RGS-WILSON-007 — Correzione di convenzione Euclidea
**Enunciato.** Nella convenzione in cui (D_0=D_0^dagger) è l'operatore spaziale hermitiano, l'operatore Euclideo Wilson usato per overlap è (D_E=iD_0+W), non la somma hermitiana legacy (D_H=D_0+W).

**Stato:** IDENTITÀ DI CONVENZIONE / VERIFICATA tramite simbolo su anello periodico. La somma (D_H) presenta uno zero accidentale ulteriore in (pa=2arctan(2/r)) in 1D; quindi il solo valore a (p=pi/a) non certifica l'assenza di altri zeri.  
**Fonte:** Phase 8 in `docs/WEILER_PHASE6_THEOREM_MAP.md` e sezione 7 di `docs/WILSON_REFINEMENT_OVERLAP_DIRAC_KAHLER.md`.  
**Limite:** non è un teorema no-doubling per il grafo curvo Q6.

### WILSON-SCALE-008 — Finestra di scala del parametro Wilson su reticolo periodico
**Enunciato.** Per il simbolo periodico regolare, con (r(a)=r_0(a/L)^alpha), il contributo Wilson sui modi fisici fissati svanisce se (alpha>-1), mentre le masse agli angoli di Brillouin divergono se (alpha<1). La finestra congiunta è
[
-1<alpha<1.
]

**Stato:** PROVATO ANALITICAMENTE per il simbolo Fourier periodico regolare descritto nel documento Phase 9.  
**Fonte:** `src/wilson_refinement.py`, `examples/wilson_refinement_audit.py`, `docs/WILSON_REFINEMENT_OVERLAP_DIRAC_KAHLER.md`.  
**Limite:** non trasferire la classificazione degli angoli di Brillouin al grafo Q6 finito, che non ha un toro di Brillouin canonico.

### MOSCO-TORUS-009 — Prototipo Mosco sul toro piatto
**Enunciato.** Per la famiglia periodica piatta con embedding Fourier isometrico nello stesso (L^2(mathbb T_L^d;mathbb C^4)), forme estese con (+infty) fuori dagli spazi discreti e scaling indicato, il documento dimostra liminf e recovery sequence per la forma limite (q(u)=4int|
abla u|^2). La convergenza Mosco vale nel regime di consistenza (alpha>-1); la finestra congiunta (-1<alpha<1) aggiunge disaccoppiamento dei doublers e compattezza spettrale nel modello piatto.

**Stato:** PROTOTIPO DIMOSTRATO per il toro periodico piatto, sotto le ipotesi esplicite del documento.  
**Fonte:** `docs/MOSCO_PROTOTYPE_REGULAR_TORUS.md`, `src/mosco_torus.py`, `tests/test_mosco_torus.py`.  
**Limite decisivo:** NON è una prova Mosco per RGS-Q6. Mancano famiglia di complessi curvi con mesh (h	o0), embedding comuni, controllo geometrico uniforme, consistenza dei link e prove liminf/recovery su quella famiglia.

### DK-Q6-010 — Dirac–Kähler combinatorio sul 2-scheletro Q6
**Definizione.** Con (d_0=B_1^T), (d_1=B_2^T), (d_1d_0=0) e (D_{DK}=d+d^dagger) su (C^0oplus C^1oplus C^2), la dimensione è (64+192+240=496). Il quadrato è il Laplaciano di Hodge a blocchi; la parità di grado anticommuta con (D_{DK}).

**Stato:** PROVATO ALGEBRICAMENTE per il complesso finito e prodotto interno euclideo non pesato; test numerici implementati.  
**Fonte:** `src/dirac_kahler.py`, `docs/WILSON_REFINEMENT_OVERLAP_DIRAC_KAHLER.md`.  
**Limite:** non è ancora DEC pesato dal vierbein RGS; cochain Dirac–Kähler e spinore a quattro componenti sono rappresentazioni differenti, senza equivalenza asserita.

### OVERLAP-GW-011 — Identità Ginsparg–Wilson finita
**Enunciato.** Per kernel Euclideo (gamma_5)-hermitiano (D_E), (H_W=Gamma_5(aD_E-ho I)) hermitiano e gap non nullo, la costruzione overlap
[
D_{ov}=rac{ho}{a}left[I+Gamma_5operatorname{sign}(H_W)ight]
]
soddisfa la relazione Ginsparg–Wilson
[
Gamma_5D_{ov}+D_{ov}Gamma_5
=rac{a}{ho}D_{ov}Gamma_5D_{ov}.
]

**Stato:** PROVATO ALGEBRICAMENTE sotto le ipotesi indicate; test finito 2x2 in CI e audit numerico RGS-Q6.  
**Dati Q6 Phase 10:** dimensione 256; gap kernel a (ho=1): 2.0032656; residuo GW (3.27	imes10^{-14}); (gamma_5)-hermiticity (8.63	imes10^{-16}); indice finito numericamente (-8.19	imes10^{-15}), vicino a 0; nessun modo zero entro tolleranza (10^{-8}).  
**Fonte:** `src/overlap_dirac.py`, `results/json/rgs_q6_phase10_index_spectral_flow.json`.  
**Limite:** il valore indice zero è un risultato per quella matrice e convenzione, non un teorema d'indice continuo.

### OVERLAP-LOC-012 — Profilo di località overlap finito
**Risultato.** Le norme medie dei blocchi overlap per distanza di Hamming Q6 sono decrescenti da 0.34020 (distanza 1) a 0.00511 (distanza 6).

**Stato:** DIAGNOSTICA NUMERICA FINITA. La scansione di nove valori (hoin[0.5,1.5]) registra gap positivo e 128 autovalori negativi a ogni campione; nessun attraversamento è rilevato ai campioni.  
**Fonte:** `results/json/rgs_q6_phase10_index_spectral_flow.json`.  
**Limite:** non prova località esponenziale uniforme, perché non esiste ancora una famiglia di grafi raffinati né una stima uniforme.

### THERM-MOB-013 — Termodinamica spettrale Möbius
**Risultato.** Le funzioni di partizione fermioniche, le quantità termodinamiche normalizzate rispetto al mare di Dirac e l'entropia di entanglement sono calcolate da insiemi finiti di 24 livelli per settore e 24 autovalori di correlazione.

**Stato:** NUMERICAMENTE VERIFICATO per i dati spettrali di riferimento.  
**Fonte:** `docs/MOBIUS_HAWKING_THERMODYNAMICS.md`, `src/thermodynamics.py`.  
**Limite:** (T_H=1/L_b) è input del modello; non sono dimostrati orizzonte di Killing, radiazione Hawking continua, greybody factors o flusso fisico.

### WEILER-EQV-014 — Vincoli di equivarianza geometrica
**Enunciato.** Le mappe equivarianti richiedono kernel-intertwiner rispetto al gruppo e alle rappresentazioni specificate; le convoluzioni su gruppi affini, le formule su varietà e le condizioni di isometria valgono solo sotto le rispettive ipotesi geometriche.

**Stato:** STANDARD come teoria di riferimento; applicazione Q6 parziale. Sono aperti: sottogruppo di automorfismi metrici che preserva i dati (X,E), kernel basis completa, bias/pooling invarianti, atlas/cociclo globale e transitivity/homogeneity.  
**Fonte:** `docs/WEILER_PHASE6_THEOREM_MAP.md`.  
**Limite:** non si dichiara piena conformità ai teoremi Weiler per Q6.

## Tabella sintetica

| ID | Oggetto | Stato | Ambito valido |
|---|---|---|---|
| Q6-CX-001 | Incidenza e conteggi Q6 | Provato/esatto | Complesso finito canonico |
| Q6-SP-002 | Spettro Laplaciano Q6 | Provato | Grafo ipercubico |
| RGS-GEO-003 | Clifford da coframe | Provato algebricamente | Ogni sito con coframe invertibile |
| RGS-SPIN-004 | Lift polare SU(2) | Candidato + numerico | Archivio RGS-Q6 |
| RGS-DIRAC-005 | Dirac spaziale Q6 | Candidato + numerico | Grafo finito ricostruito |
| RGS-WILSON-006 | Laplaciano Wilson | Provato algebricamente | Link unitari, pesi non negativi |
| RGS-WILSON-007 | Convenzione Euclidea | Identità di convenzione | Simbolo e operatori specificati |
| WILSON-SCALE-008 | Finestra (alpha) | Provato | Reticolo periodico regolare |
| MOSCO-TORUS-009 | Convergenza Mosco | Prototipo dimostrato | Toro piatto periodico |
| DK-Q6-010 | Dirac–Kähler | Provato finito | 2-scheletro Q6 non pesato |
| OVERLAP-GW-011 | Relazione GW | Provato sotto ipotesi | Kernel finito gappato |
| OVERLAP-LOC-012 | Località overlap | Diagnostica numerica | Una matrice Q6 |
| THERM-MOB-013 | Termodinamica | Numerico | Spettri di riferimento finiti |
| WEILER-EQV-014 | Equivarianza | Standard + incompleto | Applicazione Q6 parziale |

## Gate aperti prioritari

1. **RGS-MOSCO:** costruire una famiglia (K_h) di complessi curvi RGS con (h	o0), embedding comuni, misure e Hodge stars uniformemente controllati; dimostrare liminf e recovery.
2. **RGS-SPIN:** dimostrare consistenza del lift polare con la connessione spin prevista, includendo olonomie e segni.
3. **RGS-OVERLAP-LOC:** dimostrare gap e località uniformi su tutta la famiglia, non su un singolo archivio.
4. **RGS-INDEX:** collegare l'indice finito a un teorema d'indice nel limite continuo.
5. **RGS-WEILER:** calcolare il sottogruppo effettivo di automorfismi che preserva (X,E), poi completare kernel/intertwiner e spazi invarianti.
6. **RGS-DEC:** introdurre Hodge star e inner product pesati dal vierbein e confrontare Dirac–Kähler con l'operatore spinoriale senza presupporre equivalenza.

## Fonti primarie nel repository

- `docs/WEILER_PHASE6_THEOREM_MAP.md`
- `docs/CHIRAL_WILSON_AUDIT.md`
- `docs/MOSCO_PROTOTYPE_REGULAR_TORUS.md`
- `docs/WILSON_REFINEMENT_OVERLAP_DIRAC_KAHLER.md`
- `docs/MOBIUS_HAWKING_THERMODYNAMICS.md`
- `results/json/rgs_q6_phase6_reference_audit.json`
- `results/json/rgs_q6_wilson_audit.json`
- `results/json/rgs_q6_phase10_index_spectral_flow.json`

**Regola di aggiornamento:** gli ID non si riutilizzano. Un cambiamento di ipotesi, convenzione o oggetto produce una nuova versione o un nuovo ID; un risultato numerico non diventa prova senza una dimostrazione separata.

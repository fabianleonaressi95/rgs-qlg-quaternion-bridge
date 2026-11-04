# Programma TOE RGS–PSU–PRISM–QLG fondato su gradienti e baricentri

**Data:** 2026-11-04  
**Stato:** programma di ricerca e specifica d'integrazione; non è una TOE dimostrata.

## 1. Scopo e regole epistemiche

Questo documento collega i filoni RGS, PSU–Anti-Twin, PRISM, Q6, QLG, vierbein, spinori, gauge, QCD, magnetoidrodinamica e meccanica. Non sostituisce gli audit o i moduli già presenti: stabilisce interfacce e obblighi di verifica.

Etichette obbligatorie:
- **PROVED**: dimostrazione matematica con ipotesi esplicite.
- **NUMERICALLY VERIFIED**: audit riproducibile con residui, tolleranze e revisione del codice.
- **CANDIDATE**: costruzione definita ma non ancora dimostrata.
- **OPEN**: derivazione o verifica mancante.
- **PHYSICAL HYPOTHESIS**: interpretazione o previsione non confermata.

Un'identità finita non prova una teoria continua. Una proposta di azione non equivale a equazioni fisiche validate.

## 2. Principio geometrico: baricentri e gradienti

Sia \(\sigma\) un simplesso non degenere con vertici \(X_0,\ldots,X_d\). Le coordinate baricentriche affini soddisfano

\[
\sum_{i=0}^{d}\lambda_i(x)=1,\qquad
x=\sum_{i=0}^{d}\lambda_i(x)X_i,\qquad
\sum_i\nabla\lambda_i=0.
\]

Per pesi \(w_i\geq0\), \(\sum_iw_i=1\), il baricentro pesato dei vertici e le coordinate relative sono

\[
X_B=\sum_iw_iX_i,\qquad Y_i=X_i-X_B,\qquad
\sum_iw_iY_i=0.
\]

La dispersione geometrica intorno al baricentro è

\[
\Sigma_B=\sum_iw_iY_iY_i^T.
\]

\(\Sigma_B\) descrive la distribuzione dei vertici; non è automaticamente la metrica fisica. La metrica deve provenire dalla geometria immersa o da un campo metrico esplicitamente definito.

Per un campo P1 affine a tratti, \(f_h=\sum_i f_i\lambda_i\),

\[
\nabla f_h=\sum_i f_i\nabla\lambda_i,\qquad
E_\sigma[f_h]=\int_\sigma |\nabla f_h|_g^2\,d\mu_g.
\]

Se \(J\) è la matrice delle colonne \(X_a-X_0\), la metrica nelle coordinate locali è \(G=J^TJ\); la norma del gradiente usa \(G^{-1}\). Occorre controllare rango e non-degenerazione prima di invertire \(G\).

Per coordinate dinamiche \(X_i(t)\),

\[
\dot X_B=\sum_iw_i\dot X_i,\qquad
\sum_iw_i\dot Y_i=0.
\]

Questa separa moto collettivo e fluttuazioni interne. Il centro di massa inerziale richiede pesi proporzionali alle masse inerziali; pesi geometrici o probabilistici non vanno chiamati massa senza una derivazione.

## 3. Principio variazionale discreto

Un prototipo per un campo scalare sui vertici è

\[
S_h[f]=\int dt\left[
\frac12\dot f^TM_h\dot f
-\frac12f^TK_hf-U_h(f)
\right],
\]

con matrice di massa positiva \(M_h\) e stiffness \(K_h\) assemblata dall'energia di Dirichlet. Le equazioni sono

\[
M_h\ddot f+K_hf+\nabla_fU_h(f)=0.
\]

È una struttura variazionale discreta standard. L'identificazione dei campi con gradi di libertà fisici RGS/PSU resta da derivare.

## 4. Q6, PRISM e complessi discreti

Il repository formalizza Q6 con 64 vertici, 192 archi e 240 quadrati. Con orientazioni coerenti, le matrici di incidenza devono soddisfare \(B_1B_2=0\). Distinguere:
- 0-cochain: valori di campo sui vertici;
- 1-cochain: differenze o correnti orientate sugli archi;
- 2-cochain: flussi, olonomie o curvature sulle facce;
- Hodge star: dipende da geometria, misure e inner product, non dalla sola incidenza.

Il rango ciclico del grafo non coincide con i numeri di Betti del complesso dopo l'attacco delle facce. Una mappa parziale di un ciclo in PRISM non è un isomorfismo dell'intero complesso.

## 5. Vierbein, connessione e torsione

Il vierbein collega indici interni e spazio-temporali:

\[
g_{\mu\nu}=\eta_{ab}e^a{}_{\mu}e^b{}_{\nu},\qquad
e=|\det(e^a{}_{\mu})|.
\]

La torsione è

\[
T^a=de^a+\omega^a{}_b\wedge e^b.
\]

Una connessione teleparallela piatta non implica torsione nulla. Per identificare una costruzione RGS con TEGR occorre specificare tetrad invertibile, segnatura \((-+++ )\), connessione, trasformazioni locali di Lorentz, scalare torsionale, termini di bordo e confronto con la formulazione metrica. Una coframe discreta non è automaticamente un vierbein Lorentziano.

## 6. Spinori, fermioni e correnti di rete

Con \(\Gamma^\mu=e_a{}^\mu\gamma^a\), una derivata spinoriale gauge-covariante ha forma

\[
D_\mu\psi=
\left(\partial_\mu+\frac14\omega_{\mu ab}\gamma^{ab}
+igA_\mu^AT_A\right)\psi.
\]

La corrente U(1) è \(j^\mu=q\bar\psi\Gamma^\mu\psi\). La conservazione locale deve seguire dalle equazioni del moto e dalle simmetrie pertinenti.

Su una rete orientata, \(B_1j=s\) esprime il bilancio di corrente ai nodi; senza sorgenti \(B_1j=0\). È l'analogo discreto della continuità, ma il passaggio al continuo richiede ricostruzione e convergenza. L'analogia Kirchhoff non trasforma uno spinore in una corrente classica: collega strutture di conservazione a livelli diversi.

## 7. Gauge, elettromagnetismo, monopoli e QCD

Per un gruppo G, i link \(U_{xy}\in G\) trasportano indici interni; l'olonomia di faccia è il prodotto ordinato lungo il bordo. Covarianza dei link e curvatura discreta vanno testate separatamente.

Distinzioni essenziali:
- \(SU(3)_C\): gauge del colore;
- \(SU(3)_F\): simmetria approssimata di sapore;
- \(U(1)_{EM}\): carica elettromagnetica.

Per i quark leggeri \(Q_u=+2e/3\), \(Q_d=Q_s=-e/3\). La decomposizione di sapore è

\[
3\otimes3\otimes3=10\oplus8\oplus8\oplus1.
\]

Ottetti e decupletto sono multiplette di sapore. Il colore di un barione è un singoletto costruito con \(\epsilon_{abc}\); spin, sapore, colore e parte spaziale devono combinarsi rispettando l'antisimmetria fermionica totale. La decomposizione dimensionale non determina da sola masse o dinamica QCD.

La condizione di quantizzazione di Dirac dipende dalla normalizzazione delle cariche e dalla convenzione razionalizzata/non razionalizzata: dichiarare le convenzioni prima di confrontare \(eg\) con \(\hbar c/2\).

## 8. Magnetoidrodinamica e materia continua

La MHD è un limite macroscopico di un fluido conduttore accoppiato all'elettromagnetismo, non un sostituto di QED/QCD. In SI, una forma resistiva non relativistica comprende

\[
\partial_t\mathbf B=-\nabla\times\mathbf E,\qquad
\nabla\times\mathbf B=\mu_0\mathbf J+
\mu_0\epsilon_0\partial_t\mathbf E,
\]

\[
\mathbf J=\sigma(\mathbf E+\mathbf v\times\mathbf B).
\]

La MHD ideale pone \(\mathbf E+\mathbf v\times\mathbf B=0\). La derivazione da gradi di libertà microscopici richiede coarse-graining, chiusure costitutive e bilanci di carica, energia e quantità di moto.

## 9. Interfaccia implementativa per ottetti e decupletti

Prima fase, puramente algebrica:
1. fissare convenzione per generatori di \(su(3)\);
2. costruire \(3\otimes3\otimes3\);
3. verificare \(10+8+8+1=27\);
4. costruire proiettori con Casimir/operatori di permutazione;
5. verificare commutatori, Casimir e dimensioni;
6. tenere separati sapore e colore;
7. aggiungere spin e simmetria totale in un secondo livello.

Il collegamento a Q6/PRISM è una possibile struttura di fibre: uno spazio interno rappresentativo per sito e trasporto gauge sui link. L'associazione di una rappresentazione a un particolare settore fisico resta CANDIDATE finché non è derivata dall'azione e dai vincoli.

## 10. Azione TOE: architettura candidata

Una possibile decomposizione è

\[
S_{candidate}=S_{grav}[e,\omega]+S_{gauge}[A]
+S_{Dirac}[\psi,e,\omega,A]+S_{matter}[\Phi]
+S_{fluid}[u,\rho,s]+S_{PSU}[\Xi,e,\omega]+S_{int}.
\]

Una forma schematica TEGR+gauge è

\[
S=\int d^4x\,e\left[
\frac{c^4}{16\pi G}T
-\frac14\sum_A F^A_{\mu\nu}F_A^{\mu\nu}
+\mathcal L_D+\mathcal L_{matter}
+\mathcal L_{PSU}+\mathcal L_{fluid}+\mathcal L_{int}
\right]+S_{\partial}.
\]

È un inventario di settori, non un'azione completata. Devono essere definiti unità, accoppiamenti, convenzioni per T, rappresentazioni, campi PSU, termini di bordo e variazioni indipendenti. Non aggiungere accoppiamenti per sola analogia.

## 11. De Aquino: dossier critico

Il file utente “De Aquino – Equazioni principali.md” riporta un'estrazione della versione v25 (2005) di *Kinetic Quantum Theory of Gravity*, non la v44 del 2014. Il dossier stesso segnala accesso incompleto alla raccolta HAL e assenza di validazione mainstream.

Relazioni trascritte nel dossier:

\[
m_g=m_i\left\{1-2\left[
\sqrt{1+\left(\frac{\Delta p}{m_ic}\right)^2}-1
\right]\right\},
\]

\[
\frac{d\mathbf p}{dt}=\frac{m_g}{m_i}
q(\mathbf E+\mathbf V\times\mathbf B).
\]

La nota sorgente evidenzia un'identificazione non giustificata nel passaggio alla relazione massa-impulso e assenza di replicazione indipendente nota per le affermazioni ELF/mumetal. Trattare queste formule come SOURCE-REPORTED / HYPOTHESIS, non come leggi validate.

Audit richiesto: ricostruzione algebrica dei passaggi; analisi dimensionale e dei limiti; confronto con relatività speciale e Lorentz standard; eventuale derivazione da azione locale e gauge-invariante; bilanci energia-impulso; confronto con vincoli sperimentali. Non usarle come correzioni predittive prima dell'audit.

## 12. Tutte le meccaniche come limiti controllati

| Settore | Recupero richiesto |
|---|---|
| Newton | limite di velocità piccole e campi deboli |
| Lagrange/Hamilton | principio variazionale e trasformata di Legendre |
| Relatività | tetrad Lorentziano, covarianza e geodetiche |
| Maxwell | settore U(1), gauge e onde |
| Yang–Mills/QCD | gruppo non abeliano, rappresentazioni e Wilson loop |
| Dirac | Clifford, corrente, hermiticità e limite non relativistico |
| Statistica/termodinamica | ensemble, entropia e bilanci |
| Fluidi/MHD | coarse-graining, chiusure e conservazione |
| Quantistica | spazio di Hilbert, unitarietà e osservabili |
| PSU/Anti-Twin | campi, simmetrie, stabilità e limiti osservativi |

Non basta riunire i nomi: ogni limite deve essere derivato dall'azione o da un regime controllato.

## 13. Roadmap nel repository

**A — Inventario:** mappare moduli esistenti, convenzioni di unità, indici, orientazioni, segnature e normalizzazioni. Non duplicare codice.

**B — Core baricentrico:** coordinate, gradienti, metriche locali, simplici degeneri, matrici di massa/stiffness. Testare partizione dell'unità, riproduzione affine, somma dei gradienti nulla, positività e invarianza per traslazioni.

**C — SU(3):** generatori, prodotti tensoriali, Casimir, proiettori e test della decomposizione 10+8+8+1.

**D — Rete/gauge:** orientazioni, incidenza, sorgenti, link, olonomie, covarianza e conservazione discreta.

**E — Vierbein/spinori:** coframe, segnatura, connessione di spin, torsione, covarianza locale, hermiticità e corrente on-shell.

**F — MHD:** modello microscopico, coarse-graining, equazioni costitutive e bilanci.

**G — Azione/convergenza:** assemblare solo termini definiti; variare ogni campo; verificare Noether e vincoli; dimostrare liminf e recovery sequence per ogni affermazione Mosco.

## 14. Test di accettazione

- [ ] Coordinate baricentriche: partizione dell'unità e riproduzione affine.
- [ ] Gradienti: somma nulla e trasformazione corretta.
- [ ] Massa/stiffness: simmetria, positività e kernel atteso.
- [ ] Q6: \(B_1B_2=0\), ranghi e convenzioni.
- [ ] SU(3): commutatori, Casimir e dimensioni dei settori.
- [ ] Gauge: covarianza dei link e invarianza degli osservabili.
- [ ] Spinori: Clifford, connessione e corrente conservata on-shell.
- [ ] Vierbein: determinante non nullo, segnatura e trasformazioni Lorentz.
- [ ] MHD: bilanci e regime di validità.
- [ ] Azione: variazioni, unità dimensionali e identità di Noether.
- [ ] Continuo: stime uniformi, compattezza e convergenza.
- [ ] De Aquino: audit analitico indipendente.

## 15. Moduli e documenti esistenti da riusare

- `src/q6_complex.py`
- `src/unified_discrete_action.py`
- `src/unified_discrete_gradients.py`
- `src/discrete_spin_geometry.py`
- `docs/UNIFIED_Q6_VARIATIONAL_CORE.md`
- `docs/UNIFIED_Q6_EULER_LAGRANGE.md`
- `docs/Q6_ANALYTIC_GRADIENTS_AND_SPIN_CONNECTION_GATE.md`
- `docs/Q6_COMPLEX_FORMALIZATION.md`
- `docs/RGS_MOSCO_CURVED_REFINEMENT_PROGRAM.md`
- `docs/RGS_ELLIS_CURVED_REFINEMENT_G1.md`

Questi sono punti di partenza da ispezionare, non prova che tutti i settori siano già integrati.

## 16. Criterio per poter parlare di TOE

Servono: campi e simmetrie espliciti; azione dimensionalmente coerente; equazioni derivate per variazione; identità di Noether e vincoli consistenti; recupero dei limiti fisici noti; limite continuo dimostrato; predizioni quantitative falsificabili.

Fino al soddisfacimento di questi requisiti, la denominazione corretta è **programma TOE RGS–PSU–PRISM–QLG**.


## 17. Prima implementazione computazionale (PR #67)

Sono stati aggiunti, senza sostituire i moduli variazionali Q6 esistenti:

- `src/barycentric_geometry.py`: coordinate baricentriche affini, gradienti P1, baricentro pesato, tensore di dispersione, volume del simplesso, matrice di massa consistente e stiffness di Dirichlet.
- `src/su3_flavour.py`: generatori di Gell-Mann normalizzati come (T_a=\lambda_a/2), rappresentazione su (3\otimes3\otimes3), proiettori di permutazione per i settori totalmente simmetrico/antisimmetrico/misto e spettro del Casimir totale.
- `tests/test_toe_geometric_su3_core.py`: test di ricostruzione affine, partizione dell'unità, somma dei gradienti nulla, positività delle matrici FEM, rigetto di simplessi degeneri e dimensioni SU(3).

La decomposizione verificata a livello di dimensione e Casimir è (10\oplus16\oplus1), dove il settore misto di dimensione 16 corrisponde a due copie equivalenti dell'ottetto: (10\oplus8\oplus8\oplus1). Questo è un audit di rappresentazione finito-dimensionale, non una simulazione QCD né una derivazione dello spettro adronico.

Comando di test mirato: `python -m pytest tests/test_toe_geometric_su3_core.py -q`. L'esito CI deve essere controllato prima di promuovere questi moduli a baseline.

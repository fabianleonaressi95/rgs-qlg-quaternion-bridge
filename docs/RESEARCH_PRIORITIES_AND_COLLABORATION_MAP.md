# Priorità di ricerca e mappa collaborazioni RGS–QLG

**Snapshot: 4 novembre 2026.** Documento operativo: problemi ordinati per dipendenza matematica e interlocutori scelti per competenza specifica. Non implica endorsement o collaborazioni già avviate.

## 1. Stato di partenza verificato nel repository

Il repository contiene già:
- Q6 cubico con incidenze e identità di catena esatte; Hodge spectrum e Gauss law.
- Autospazio coesatto lambda=12 di dimensione 5, mappa razionale di rango 5 e kernel intero di rango 124.
- Peso Wilson puntuale non invariante sulla fibra; peso corretto definito dalla somma sull'intera fibra, positiva e assolutamente convergente.
- Base satura del kernel costruita con Smith normal form + row Hermite normal form; tail bound fattoriale per flussi oltre cutoff.
- Identità Fourier–Haar esatta sul toro annichilatore di dimensione 68; per il reticolo Q6 auditato gli invarianti di Smith sono tutti 1.
- Carta Maxwell locale verso C^5 e pullback Bergman D5 esatti sul preimage di D5; non è una mappa globale del quoziente Wilson.
- Famiglia curva Ellis–RGS G1 e embedding comune G4 implementati; G5 uniforme e G6 Mosco curvo ancora aperti.
- Mosco G6 dimostrato per P1 periodico piatto, non trasferibile automaticamente al caso curvo.
- Dirac/overlap e Yang–Mills: identità e audit finiti, ma non ancora limite continuo curvo.

Fonti canoniche: `docs/Q6_GLOBAL_TORUS_STRUCTURE.md`, `docs/Q6_D5_INTEGER_FIBER_QUOTIENT.md`, `docs/Q6_D5_SATURATED_KERNEL_AND_TAIL_BOUND.md`, `docs/Q6_D5_WILSON_FIVE_TORUS_MARGINAL.md`, `docs/Q6_D5_BERGMAN_LIOUVILLE_PULLBACK_THEOREM.md`, `docs/RGS_MOSCO_CURVED_REFINEMENT_PROGRAM.md`, `docs/RGS_ELLIS_CURVED_REFINEMENT_G1.md`, `docs/RGS_ELLIS_CURVED_G4_COMMON_L2.md`, `docs/FLAT_TORUS_G6_P1_MOSCO_PROOF.md`.

## 2. Problemi in ordine di priorità

### P0 — Marginale Wilson sul quoziente Q6–D5

**Domanda:** valutare con errore certificato la somma Wilson sulle fibre intere e confrontarla con il limite gaussiano fixed-flux/large-beta.

Già stabilito: mappa razionale rango 5; kernel intero rango 124; peso puntuale non discende al quoziente; somma completa positiva e assolutamente convergente; base satura SNF/HNF; maggiorante globale fattoriale; identità Fourier–Haar esatta.

Prossimi passi:
1. Usare la base satura di `build_q6_d5_saturated_integer_kernel()`, non i soli vettori interi indipendenti restituiti dalla routine preliminare di nullspace.
2. Formare (A=B_KC), (P_\perp=I-A^T(AA^T)^{-1}A), (Q_{eff}=CP_\perp C^T), verificando rango, kernel, positività sul quoziente e invarianza unimodulare.
3. Derivare il controllo somma theta/integrale gaussiano tramite reticolo duale e Poisson; non sostituire la somma con l'integrale senza stima.
4. Valutare la fibra con enumerazione branch-and-bound/lattice reduction e aggiungere il tail bound certificato.
5. Studiare il limite in beta con controllo uniforme dei flussi rilevanti.

Non segue ancora una misura Bergman, un cutoff Hua fisico o una predizione per alpha.

File: `src/rgs_q6_d5_integer_fiber_quotient.py`, `src/rgs_q6_d5_saturated_kernel_bounds.py`, `src/rgs_q6_d5_wilson_fiber_torus.py`.

### P1 — Chiudere G5/G6 per la geometria curva Ellis–RGS

**Domanda:** le forme P1 sulla famiglia curva convergono in senso Mosco alla forma di Dirichlet sulla stessa varietà ([-L,L]\times S^2)?

Già stabilito: G1 famiglia di mesh da seed RGS 64B; G4 embedding comune in (L^2(M_L,d\mu_g)) e assemblaggio con metrica curva. Il teorema G6 piatto è separato.

Prossimi passi:
1. Shape-regularity uniforme all-level delle mappe tetraedriche e degli inversi.
2. Ellitticità uniforme del tensore metrico, Jacobiani e misure.
3. Coercività uniforme e compattezza Rellich discreta (G5).
4. Prova liminf debole e recovery forte (G6), prima su funzioni lisce e poi per densità in H1.
5. Separare il limite (L\to\infty) della wormhole completa dal problema compatto con bordo.

### P2 — Limite continuo Dirac/overlap/gauge

1. Hodge star DEC pesata dal vierbein, non solo operatore combinatorio.
2. Dimostrare che i link SU(2) approssimano il trasporto della connessione spin; il lift polare finito non basta.
3. Consistenza su campi lisci e stabilità spettrale uniforme.
4. Doubling, indice e gap uniforme del kernel Wilson.
5. Località overlap con stime uniformi sulla famiglia.
6. Convergenza resolvente/operatoriale solo dopo P1.

### P3 — Misura globale D5/Hua dal modello Wilson

La carta Maxwell locale e il pullback (N^{-5}d\mu_L) sono risultati matematici già registrati. Resta aperto se la marginale Wilson completa, includendo il settore spaziale e il coniugato canonico, induca quella misura. Servono pushforward globale, controllo angolare, tightness e derivazione non adattata del regolatore epsilon. La misura Haar su un toro compatto e il volume Bergman non compatto non sono intercambiabili.

### P4 — Azione dinamica RGS–PSU/TEGR/Einstein–Cartan

Specificare un'unica azione con coframe, connessione e materia indipendenti; derivare simmetrie, equazioni di Eulero–Lagrange, Bianchi, correnti conservate, limiti TEGR/Einstein e accoppiamento fermionico. Le verifiche finite di torsione, tetradi e olonomie non sostituiscono questa derivazione.

### P5 — PRISM/Krist-All/Koopman

Mantenere questo ramo distinto da P0–P3 finché non esiste un ponte dimostrato. Servono misura invariante, operatori definiti, compattezza/convergenza di forme o semigruppi e leggi di conservazione derivate. Non usare coincidenze di dimensione o rapporti numerici come identificazioni fisiche.

## 3. Persone da contattare

Email riportate da pagine istituzionali/personali pubbliche. Ogni richiesta deve essere circoscritta a un lemma o benchmark, non alla validazione dell'intero programma.

| Priorità | Persona | Match tecnico | Proposta | Contatto |
|---|---|---|---|---|
| A1 | **Ralf Hiptmair**, ETH Zürich | FEEC, PDE numeriche, discretizzazioni compatibili | Revisione delle ipotesi G4–G6 e della famiglia curva | ralf.hiptmair@sam.math.ethz.ch |
| A2 | **Ragnar Winther**, University of Oslo | FEEC, analisi numerica, discretizzazioni compatibili | Verifica indipendente di coercività, compattezza, liminf/recovery | rwinther@math.uio.no |
| A3 | **Anita Maria Tabacco**, Politecnico di Torino | Analisi armonica/funzionale, operatori su grafi e varietà | Revisione della parte Bergman/Hua e delle misure/operatori discreti | anita.tabacco@polito.it |
| A4 | **Alessio Martini**, Politecnico di Torino | Moltiplicatori spettrali, analisi armonica, gruppi e geometria sub-Riemanniana | Audit spettrale del grafo/quoziente e confronto con proiettori Bergman | alessio.martini@polito.it |
| B1 | **Luigi Del Debbio**, University of Edinburgh | QFT non perturbativa, lattice gauge theory | Audit della marginale Wilson e del limite fixed-flux/trasferimento | luigi.del.debbio@ed.ac.uk |
| B2 | **Peter Boyle**, University of Edinburgh | QCD numerica, lattice gauge theory | Confronto tecnico Wilson/overlap, kernel e località | paboyle@ph.ed.ac.uk |
| B3 | **Douglas Arnold**, University of Minnesota | FEEC, PDE e analisi numerica | Commento sul trasferimento del teorema P1 piatto alla famiglia curva | arnold@umn.edu |
| B4 | **Richard Falk**, Rutgers University | FEEC e metodi agli elementi finiti | Revisione delle stime uniformi e della struttura compatible/mixed | falk@math.rutgers.edu |
| C1 | **Luis Wirth**, ETH Zürich, autore di formoniq | FEEC coordinate-free, Regge e geometria pseudo-Riemanniana | Benchmark/interop: esportare il caso Ellis–RGS e confrontare Hodge–Dirac | sito personale/GitHub; non è stata verificata un'email pubblica |

Ordine d'invio consigliato: Hiptmair e Winther per P1; Tabacco e Martini per P3; Del Debbio e Boyle per P0/P2; Arnold e Falk come revisori FEEC di seconda ondata; Wirth per il benchmark software.

Non chiedere a esperti di lattice QCD di validare la geometria D5, né a esperti FEEC di validare la fisica di alpha.

## 4. Repository esterni utili

| Repository | Contenuto | Uso concreto |
|---|---|---|
| [luiswirth/formoniq](https://github.com/luiswirth/formoniq) | Engine FEEC Rust, complessi simpliciali, Regge, metriche pseudo-Riemanniane, Hodge–Laplace | Benchmark indipendente P1/P2 |
| [luiswirth/formoniq-studio](https://github.com/luiswirth/formoniq-studio) | Viewer per formoniq | Interoperabilità e visualizzazione, non prova |
| [dialkforkaushik/DECAGT](https://github.com/dialkforkaushik/DECAGT) | Framework C++/Python DEC e FEEC; Kaushik Kalyanaraman e Pranav Jain | Confronto incidenze, Hodge star, discretizzazioni |
| [cpml-au/dctkit](https://github.com/cpml-au/dctkit) | Toolkit Python/JAX per topologia algebrica, DEC e geometria discreta | Benchmark differenziabile per cochain e forme energetiche |
| [rigetti/DiscreteExteriorCalculus.jl](https://github.com/rigetti/DiscreteExteriorCalculus.jl) | DEC su complessi cellulari, Hodge star, codifferenziale e Laplace–de Rham | Controllo incrociato degli operatori Q6 |
| [chakravala/Grassmann.jl](https://github.com/chakravala/Grassmann.jl) | Algebra geometrica, Clifford, Hodge e algebra multilineare | Convenzioni Clifford/quaternioniche |
| [chakravala/Cartan.jl](https://github.com/chakravala/Cartan.jl) | Frame bundle, campi tensoriali, connessione e curvatura | Benchmark per coframe e connessione |
| [hirani/pydec](https://github.com/hirani/pydec) | Libreria Python classica DEC | Riferimento storico e baseline |

Non sono endorsement né collaborazioni avviate. Prima di proporre PR: controllare licenza, API, orientazioni e definizioni di Hodge star.

## 5. Regole per le email

Ogni email: (i) una sola domanda, (ii) risultato già provato e ipotesi, (iii) lacuna esatta, (iv) file/commit e test minimo, (v) richiesta piccola: commento a lemma, controesempio, revisione di stima o benchmark.

Non inviare un manifesto di TOE, claims su alpha/gravità/antimateria o un elenco di decine di risultati. La prima richiesta è una verifica indipendente.

## 6. Sequenza operativa

1. P0: calcolare e certificare (Q_{eff}) dalla base satura effettiva; allineare formula, codice e documentazione.
2. P1: chiudere G5 e preparare la prova G6 curva; contattare Hiptmair/Winther.
3. P3: separare rigorosamente pullback Bergman locale e pushforward Wilson globale; contattare Tabacco/Martini.
4. P0/P2: dopo un benchmark Wilson certificato, contattare Del Debbio/Boyle.
5. P2/P4: solo dopo convergenza geometrica, affrontare Dirac, connessione spin e azione gravitazionale.

**Regola di arresto:** se una mappa o misura non è definita canonicamente, non si passa a fit numerici o predizioni fisiche: si risolve prima la definizione mancante.

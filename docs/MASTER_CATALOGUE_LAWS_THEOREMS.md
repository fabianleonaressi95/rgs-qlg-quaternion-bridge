# RGS–PSU–PRISM–QLG: Catalogo delle leggi, teoremi e congetture

**Versione 0.1 — 3 novembre 2026**  
**Ambito:** Realismo Geometrico Strutturale (RGS), PSU, PRISM–QVSM, Krist-All, QLG e Morfodinamica Quantistica.

> Questo registro non afferma che ogni relazione sia una nuova legge o un teorema dimostrato. Distingue risultati matematici standard, identità finite, verifiche numeriche, definizioni di modello e ipotesi fisiche. Un test numerico PASS non dimostra da solo un teorema di limite continuo né una legge fisica.

## Legenda epistemica

- **PROVATO / STANDARD:** risultato matematico noto, sotto ipotesi dichiarate.
- **IDENTITÀ FINITA ESATTA:** identità algebrica o combinatoria della costruzione finita.
- **VERIFICATO NUMERICAMENTE:** test riproducibile con residui e tolleranze; non equivale a prova esatta.
- **DEFINIZIONE DI MODELLO:** equazione che definisce un operatore, un'azione o un osservabile.
- **CANDIDATO / CONGETTURA:** struttura ancora da dimostrare o far convergere.
- **IPOTESI FISICA APERTA:** interpretazione non derivata da un'azione completa o non verificata sperimentalmente.
- **NON STABILITO / BLOCCATO:** affermazione più forte non supportata dai dati recuperati.

# A. Algebra, gruppi e geometria

## A1. Costante superdorata
\[
\psi^3-\psi^2-1=0,\qquad \psi=1.465571231876768\ldots
\]
**Stato:** PROVATO come definizione algebrica della radice reale positiva maggiore di uno; residuo numerico riportato \(4.44\times10^{-16}\). L'interpretazione di \(\psi\) come soglia fisica, accoppiamento o esponente universale resta da derivare.

## A2. Rapporto aureo
\[
\phi=\frac{1+\sqrt5}{2},\qquad \phi^2=\phi+1.
\]
**Stato:** PROVATO / STANDARD. Il ponte \(\phi,\psi\), funzioni eta di Dedekind e raffinamento infinito multidimensionale è CANDIDATO / APERTO.

## A3. Quaternioni unitari e doppio ricoprimento
\[
\rho:SU(2)\to SO(3),\qquad \rho(q)=\rho(-q),
\]
\[
q^{-1}=q^*,\qquad qq^*=1,\qquad q(-q^*)=-1.
\]
**Stato:** PROVATO / STANDARD. Due quaternioni opposti rappresentano la stessa rotazione SO(3), ma sono elementi distinti di SU(2).

## A4. Gruppo binario icosaedrico
\[
2I\subset SU(2),\qquad |2I|=120.
\]
Il quoziente \(S^3/2I\) è la sfera omologica di Poincaré.
**Stato:** PROVATO / STANDARD come risultato di teoria dei gruppi e topologia. Le matrici, generatori, caratteri e intertwiners implementati nel progetto sono verifiche finite specifiche, da documentare con le rispettive convenzioni.

## A5. Lemma di Schur e intertwiners
Per rappresentazioni complesse irriducibili:
\[
\dim\operatorname{Hom}_G(V,W)=
\begin{cases}1&V\simeq W,\\0&V\not\simeq W.\end{cases}
\]
**Stato:** PROVATO / STANDARD sotto le ipotesi di irriducibilità. I risultati numerici delle matrici del progetto verificano solo le rappresentazioni testate.

## A6. Relazione di Clifford
\[
\{\Gamma^\mu,\Gamma^\nu\}=2g^{\mu\nu}I.
\]
**Stato:** PROVATO come relazione definitoria standard. Ogni realizzazione discreta richiede verifica di metrica, base, dimensioni e residui.

# B. Q6: combinatoria, topologia e spettro

## B1. Cardinalità e grado del cubo a sei dimensioni
\[
|V(Q_6)|=64,\quad |E(Q_6)|=192,\quad \deg(v)=6.
\]
**Stato:** PROVATO / STANDARD e riprodotto negli audit.

## B2. Spettro del Laplaciano di Q6
\[
\lambda_j=2j,\qquad m_j=\binom6j,\qquad j=0,\ldots,6.
\]
Molteplicità: \(1,6,15,20,15,6,1\).
**Stato:** PROVATO / STANDARD per il Laplaciano combinatorio del grafo ipercubico; audit numerico coerente.

## B3. Rango ciclico del grafo Q6
\[
\beta_1^{\rm graph}=|E|-|V|+1=129.
\]
**Stato:** IDENTITÀ FINITA ESATTA per il grafo connesso considerato come complesso unidimensionale. Non va confuso con l'omologia del complesso dopo l'aggiunta delle facce.

## B4. Confine del confine
Per le matrici di incidenza orientate vertice-spigolo e spigolo-faccia:
\[
B_1B_2=0.
\]
**Stato:** PROVATO / STANDARD per complessi cellulari; verificato numericamente per le orientazioni fissate nel codice Q6.

## B5. Numero di quadrati elementari
\[
N_{\square}=\binom62\,2^4=240.
\]
**Stato:** PROVATO / STANDARD combinatorio. L'implementazione ha ricostruito 240 plaquette canoniche e 192 link. La suddivisione 192 non degeneri + 48 degeneri appartiene a una particolare costruzione biversa e non modifica il conteggio totale del cubo.

## B6. Spettro del prodotto Q6 × Z2
Se il Laplaciano del fattore Z2 ha autovalori \(0,2g\), lo spettro del Laplaciano del prodotto cartesiano è l'unione di
\[
2j,\qquad 2j+2g,\qquad \text{molteplicità }\binom6j.
\]
**Stato:** PROVATO / STANDARD per il prodotto cartesiano e la convenzione di accoppiamento indicata.

## B7. Funzionale di plaquette
\[
F_{Q6}=\sum_{p=1}^{240}(1-\cos\Omega_p).
\]
**Stato:** DEFINIZIONE DI MODELLO, non legge universale.

## B8. Relazione del quarto momento
\[
\Delta\mu_4=8w^4\sum_{p=1}^{240}(\cos\Omega_p-1).
\]
**Stato:** IDENTITÀ ALGEBRICA FINITA riportata nella derivazione Q6, valida solo con le definizioni e normalizzazioni usate in quella costruzione. Non è una legge Yang–Mills universale.

# C. PRISM e Krist-All

## C1. Cardinalità PRISM
\[
C_{64}\square C_6^3\square C_5,\qquad
N=64\cdot6^3\cdot5=69\,120.
\]
Indicizzazione:
\[
\operatorname{idx}(G,L,C,T,B)=((((G6+L)6+C)6+T)5+B).
\]
**Stato:** IDENTITÀ COMBINATORIA FINITA ESATTA. Non dimostra di per sé un'interpretazione fisica.

## C2. Relazione polinomiale Krist-All
\[
K^3=K^2+I.
\]
La ricerca riporta 11.725 classi Krist-All.
**Stato:** IDENTITÀ ALGEBRICA FINITA RIPORTATA / NUMERICAMENTE VERIFICATA nel contesto degli audit precedenti. Per chiamarla teorema occorrono definizione completa di K, campo dei coefficienti, costruzione riproducibile e verifica indipendente in aritmetica esatta.

## C2a. Equivalenza per permutazione delle matrici supergolden

La matrice usata nel modello simbolico Phase 21,
\[
A_\psi=\begin{pmatrix}1&0&1\\1&0&0\\0&1&0\end{pmatrix},
\]
e la matrice obiettivo del trasferimento sui tre canali trigrammatici,
\[
M_T=\begin{pmatrix}0&0&1\\1&0&0\\0&1&1\end{pmatrix},
\]
sono coniugate per permutazione:
\[
M_T=P A_\psi P^{-1},\qquad
P=\begin{pmatrix}0&1&0\\0&0&1\\1&0&0\end{pmatrix}.
\]
Entrambe hanno polinomio caratteristico \(x^3-x^2-1\).
**Stato:** IDENTITÀ ALGEBRICA ESATTA. La coniugazione cambia l'ordinamento dei tre stati; non fornisce un intertwiner geometrico tra il modello Phase 21, i trigrammi e le celle Krist-All. Si vedano il [crosswalk](TRIGRAM_KRISTALL_SUPERGOLDEN_CROSSWALK.md) e il [teorema centrale T056](https://github.com/fabianleonaressi95/theorem/blob/main/theorems/T056-phase21-trigram-supergolden-permutation-conjugacy.md).

## C2b. Suddivisione dyadica orientata del 2-scheletro Q6 (T055)

La suddivisione standard per punti medi del 2-scheletro cubico Q6 costruisce mappe intere \(S_0,S_1,S_2\) con conteggi raffinati \((496,1344,960)\) e identità esatte
\[
B'_1S_1=S_0B_1,\qquad B'_2S_2=S_1B_2,\qquad B'_1B'_2=0.
\]
La matrice dei conteggi per direzione degli spigoli è inoltre \(D=2I_6\); su tre assi scelti induce \(2I_3\), non \(M_T\), e non va confusa con la compressione Walsh.
**Stato:** TEOREMA FINITO PROVATO; audit implementativo nel repository [RGS–QLG Quaternion Bridge](https://github.com/fabianleonaressi95/rgs-qlg-quaternion-bridge/blob/main/docs/Q6_DYADIC_SUBDIVISION_CHAIN_MAP.md). È una suddivisione cubica di riferimento, non una costruzione Kathara/Krist-All e non deriva \(M_T\).

## C2c. Similarità quadrata Kathara e vincolo di sovrapposizione dei centri (T057)

Formalizzando il ciclo descritto per la Krystal Spiral, la trasformazione di stadio è
\[
T=\sqrt2R_{-\pi/4}=\begin{pmatrix}1&1\\-1&1\end{pmatrix},
\qquad T^TT=2I,\qquad T^2=2R_{-\pi/2}.
\]
Ponendo AzurA nell'origine, \(p_6=(-w/2,0)\), \(p_7=(w/2,0)\), e usando la proporzione 2:1 per il segmento centrale 1–12, \(p_{12}=(0,w)\). Le due condizioni di sovrapposizione \(Tp_9=p_{12}\), \(T^2p_6=p_{12}\) implicano esattamente
\[
p_9=Tp_6=(-w/2,w/2).
\]
I tre vettori \(p_6,p_9,p_{12}\) soddisfano inoltre \(p_{12}-2p_9+2p_6=0\), cioè una ricorrenza di ordine 2 con polinomio \(\lambda^2-2\lambda+2\), non il cubico supergolden.
**Stato:** TEOREMA EUCLIDEO ESATTO CONDIZIONALE alle proporzioni e alle sovrapposizioni dichiarate nella fonte. Determina il sottoinsieme di centri 1, 6, 7, 9, 12, non tutte le 15 linee né tutte le coordinate dei 12 centri. Non deriva la matrice tribonacci. Si veda [T057](https://github.com/fabianleonaressi95/theorem/blob/main/theorems/T057-kathara-square-similarity-center-overlay.md) e la [ricostruzione implementata](KATHARA_KRYSTAL_SPIRAL_SQUARE_GEOMETRY.md).

## C3. Invarianti PRISM riportati
108 cicli di lunghezza 640; numeri di Betti Hodge riportati \(b_1=5,b_2=10\); gap spettrale \(0.009630546656\); \(\omega\simeq0.098135\).
**Stato:** VERIFICHE NUMERICHE DIPENDENTI DALLA COSTRUZIONE. Non sono leggi universali: servono operatore, condizioni al contorno e normalizzazione.

# D. Identità del settore tangente 188D

L'audit comune riporta \(K_0\in\mathbb R^{188\times188}\), proiettore \(P_G\) di rango 64 e inversa generalizzata sul complemento:
\[
P_G^2=P_G,\qquad [P_G,K_0]=0,
\]
\[
K_0K_0^{-1}=K_0^{-1}K_0=I-P_G.
\]
Residui Frobenius riportati:
- idempotenza \(5.08\times10^{-15}\);
- commutatore \(5.38\times10^{-14}\);
- entrambe le identità di inversa generalizzata circa \(4.45\times10^{-14}\).
**Stato:** VERIFICATO NUMERICAMENTE per le matrici recuperate. Il collegamento geometrico con le coordinate originali dei 64 punti e con la correzione quartica a 62 punti resta BLOCCATO: nell'audit mancavano coordinate e provenienza comune sufficienti. Anche \(Q^TQ\simeq I_{188}\), residuo riportato circa \(1.3\times10^{-14}\), è un controllo numerico di ortonormalità.

# E. Dirac, Clifford e gauge

## E1. Identità Dirac–Kähler
Per \(D=d+d^*\):
\[
D^2=dd^*+d^*d=\Delta.
\]
**Stato:** PROVATO / STANDARD sotto le usuali ipotesi su metrica, dominio e condizioni al contorno. Una versione discreta deve dichiarare prodotto interno e aggiunto.

## E2. Operatore Wilson-type adottato
\[
D_W=V^{-1/2}(H+L)V^{-1/2}.
\]
**Stato:** DEFINIZIONE DI MODELLO. Hermiticità, covarianza e proprietà spettrali vanno verificate per l'operatore specifico; la formula non prova un limite di Dirac continuo.

## E3. Vincoli SU(2) sui link
\[
U_{ij}^{\dagger}U_{ij}=I,\qquad \det U_{ij}=1.
\]
**Stato:** PROVATO come condizione definitoria di SU(2). Un audit dei link geometrici ha riportato residui massimi circa \(3.15\times10^{-16}\) e \(2.24\times10^{-16}\): verifica numerica delle matrici costruite, non unicità fisica della connessione.

## E4. Covarianza di gauge reticolare
\[
U_{ij}\mapsto G_iU_{ij}G_j^{-1}.
\]
Le olonomie di plaquette trasformano per coniugio nel vertice base; le funzioni di classe, come la traccia, restano invarianti.
**Stato:** PROVATO / STANDARD con orientazioni e convenzioni coerenti; test del codice sono finiti.

## E5. Risposta di shear dell'azione di Wilson
La derivata dell'azione di Wilson rispetto allo shear dipende da come coframe/metrica, link, connessione e misura variano.
**Stato:** MODEL-DEPENDENT. Nell'audit di una connessione spin-lift geometrica prescritta la derivata finite-difference decadeva linearmente con epsilon (pendenza log-log circa 1.00000002), compatibile con derivata prima nulla nel punto testato. L'identificazione con uno stress gauge fisico o con il tensore di Maxwell NON È STABILITA.

# F. Trasporto, fluidi e multifrattalità

## F1. Relazione di Green–Kubo
\[
\eta=\frac{1}{k_BT V}\int_0^\infty
\langle\delta\Pi_{xy}(t)\delta\Pi_{xy}(0)\rangle\,dt.
\]
**Stato:** FORMALISMO STANDARD di risposta lineare in equilibrio, sotto ipotesi di stazionarietà, ensemble, decadimento delle correlazioni e normalizzazione coerente. L'applicazione PSU/PRISM richiede dinamica fisica, misura d'equilibrio, stress, tempo e convergenza. Una serie sintetica a strain comune non è da sola una viscosità Green–Kubo.

## F2. Ponte Navier–Stokes / gauge
La derivazione di Navier–Stokes dal solo Hamiltoniano RGS visualizzato, l'identità stress Wilson = tensore Maxwell e la biiettività della mappa corrente-spinore sono state sottoposte ad audit.
**Stato:** NON STABILITE COME SCRITTE. Mancano, secondo il caso, proiettore di Leray/vincolo incomprimibile, dissipazione, evoluzione gauge accoppiata, condizioni al contorno, azione completa e dominio ristretto per l'inversione della corrente.

## F3. Scaling K41/multifrattale candidato
\[
E(\ell)\sim\ell^{2/3},\qquad \eta(\ell)\sim\ell^{-1/3},
\qquad \zeta=3\tau(q/3)-2.
\]
**Stato:** CANDIDATO / NON VERIFICATO. Servono definizioni degli osservabili, gerarchia convergente, intervallo di scala, dati e ipotesi del formalismo multifrattale. Non è un'identità universale senza tali ipotesi.

# G. Aritmetica esatta e correzioni quartiche

L'audit razionale recuperato riporta:
\[
M_{xxxx}=\frac{2734}{63},\quad
M_{xxyy}=\frac{880}{63},\quad
C_4=M_{xxxx}-3M_{xxyy}=\frac{94}{63},
\]
\[
C_4/64=\frac{47}{2016},\qquad 2\binom84=140.
\]
Confronto delle correzioni:
\[
\delta_{\rm cand}=\frac1{45570},\qquad
\delta_{\rm hist}=\frac{149}{312480},
\]
\[
\delta_{\rm hist}-\delta_{\rm cand}=\frac{199}{437472},\qquad
\delta_{\rm hist}/\delta_{\rm cand}=\frac{1043}{48}.
\]
**Stato:** aritmetica razionale e identità binomiale sono IDENTITÀ ESATTE. La derivazione geometrica della correzione storica e il collegamento con l'Hessiana tangente 188D NON SONO STABILITI.

# H. Principi e ipotesi fisiche discusse

## H1. Fattore costitutivo De Aquino
\[
f_{DA}(\xi)=3-2\sqrt{1+\xi^2}.
\]
**Stato:** PROPOSTA DI MODELLO / IPOTESI FISICA APERTA. Non è un teorema standard né una legge universale. Richiede unità coerenti, azione sottostante, equazioni di campo e confronto sperimentale.

## H2. Interpretazione Yee–Gardi a onde stazionarie
Sono state discusse affermazioni esterne: dieci centri geometrici dell'elettrone, massa come energia di onda stazionaria, carica come ampiezza d'onda, derivazioni geometriche di alpha e pi. Il confronto interno includeva un tripletto T2u ricostruito numericamente e una mappa geometria-centri d'onda.
**Stato:** ipotesi esterne e COMPATIBILITÀ STRUTTURALE NUMERICA per matrici selezionate. Non è stata provata equivalenza con la fisica dell'elettrone, né una derivazione di carica, massa o costante di struttura fine.

## H3. Spinori su base Möbius / biverso
La non-orientabilità non implica da sola antiperiodicità spinoriale. Una condizione
\[
\psi(x+L)=-\psi(x)
\]
richiede struttura Pin e olonomia/bordo specificati.
**Stato:** SCELTA DI MODELLO APERTA finché struttura Pin± e dominio dell'operatore non sono fissati.

## H4. Wormhole e recupero di simmetria PSU–Anti-Twin
I valori BVP/Hessiana/WKB registrati sono risultati numerici di un modello ridotto. Interpretare la gola come recupero della simmetria Lorentziana fra i due fogli è un'ipotesi fisica.
**Stato:** MODELLO RIDOTTO STUDIATO NUMERICAMENTE / IPOTESI FISICA APERTA.

### H4.1. Recupero dei parametri del raccordo \(r_* ,R_0,r_1,r_2\)

Sono registrati storicamente (traccia del 20-10-2026) i valori
\[
r_*=0.990707,\quad R_0=2.717814,\quad
r_1=1.089778,\quad r_2=1.981414.
\]
I numeri soddisfano, alla precisione riportata, \(r_1\simeq1.1r_*\) e \(r_2=2r_*\). Queste sono compatibilità numeriche, non definizioni originarie dimostrate. La derivazione della legge di raccordo, il significato di \(R_1\) maiuscolo e la dimostrazione che i quattro parametri risolvano il BVP RGS–PSU non sono stati recuperati.

Per il benchmark Ellis separato \(R(\ell)=\sqrt{R_0^2+\ell^2}\), si dimostrano esattamente \(R'(0)=0\), \(R''(0)=1/R_0\), \(b(R)=R_0^2/R\), \(b'(R_0)=-1\) e \(M_{MS}(\ell)=R_0^2/[2G\sqrt{R_0^2+\ell^2}]\), con limite asintotico nullo. Queste identità non dimostrano che tale benchmark sia la soluzione delle equazioni RGS–PSU.

Il codice corrente RGS_PSU_closure.py documenta invece tre background distinti (raggi di gola 2.31, 0.8292 e 2) e contiene placeholder nella chiusura dinamica; la famiglia Ellis–RGS G1 usa \(r_0=0.8292\). Non si devono fondere questi background con il valore storico \(R_0=2.717814\).

**Stato:** PARAMETRI STORICI RECUPERATI / COMPATIBILITÀ ELLIS DIMOSTRATA / DERIVAZIONE RGS–PSU E \(R_1\) APERTE. Vedi [audit di recupero e benchmark](RGS_PSU_WORMHOLE_PARAMETER_RECOVERY_AUDIT.md).


## H5. Raffinamento Cantoriano/multifrattale generato dagli esagrammi Q6
La proposta centrale è che il codice esagrammatico a sei bit generi una gerarchia indefinitamente raffinabile, di tipo Cantoriano/multifrattale.
**Stato:** CONGETTURA DI RICERCA. Servono mappa di raffinamento, spazi o misure annidate, oggetto limite, dimensione/spettro e prova di convergenza e non banalità. La combinatoria Q6 finita non dimostra da sola una costruzione Cantoriana infinita.

# I. Affermazioni respinte o bloccate

1. \(\mathrm{Gr}(2,4)\cong S^3\): falsa come scritta. La Grassmanniana reale orientata \(\mathrm{Gr}^+(2,4;\mathbb R)\) è diffeomorfa a \(S^2\times S^2\).
2. Target argenteo \(\delta^2=3+2\sqrt2\) dedotto da rapporto a taglia finita con discrepanza circa 1.94% e senza convergenza: NON PROVATO.
3. Identificazione fisica dell'Hessiana/proiettore 188D con la correzione quartica a 62 punti: BLOCCATA da coordinate/provenienza comune mancanti.
4. Teorema PRISM/Krist-All ↔ molteplicità di Reynolds Poincaré \(SL(2,5)\): APERTO; servono rappresentazioni e intertwiners espliciti.
5. Teorema operatore di Dirac ↔ coefficienti perturbativi lambda/fattore 32: BLOCCATO; manca derivazione di operatore, dominio e normalizzazione.
6. Termine cinetico gravitazionale PSU derivato da \(\Delta C\): IPOTESI; nessuna equazione derivata da azione è stata verificata.
7. Teorema di regolarità Navier–Stokes, viscosità Green–Kubo fisica o biiettività della corrente spinoriale: NON STABILITI dai calcoli recuperati.
8. Limite continuo da spettri Q6/PRISM finiti senza teorema di convergenza (ad esempio Mosco e compattezza): APERTO.

# J. Programma di dimostrazione prioritario

1. **Teorema polinomiale Krist-All:** definire esattamente K e provare \(K^3-K^2-I=0\) con aritmetica esatta o argomento strutturale.
2. **Teorema Q6-continuo:** definire successione di spazi/operatori e dimostrare convergenza Mosco/resolvente, compattezza e operatore limite.
3. **Ponte rappresentazionale 2I:** generatori, matrici, caratteri, irriducibilità e mappa esatta verso settori PRISM/Krist-All.
4. **Teorema geometrico 188D:** recuperare coordinate e proiettori originali e dimostrare la corrispondenza, non usare proxy come prova.
5. **Teorema variazionale PSU–Anti-Twin:** azione totale dimensionalmente coerente, variazione di tutti i campi, vincoli ed equazioni dei due fogli.
6. **Teorema Cantor/multifrattale:** mappa di raffinamento, misura invariante, esponenti e dimensione limite, con prova di convergenza.
7. **Ponte Yang–Mills/fluidi:** dinamica accoppiata e stress definiti; covarianza, vincoli/conservazione e limite di risposta lineare giustificato.

# K. Provenienza e manutenzione

Registro consolidato dai materiali recuperati, tra cui:
- sintesi_comune_audit.py e sintesi_comune_audit.json (2026-10-16);
- RGS_PSU_Numerical_Verification.py e RGS_PSU_Numerical_Verification.json (audit NS–SU(2));
- RGS_PSU_geometric_intertwiner_proxy_report.json;
- audit di incidenza, Laplaciano e plaquette Q6;
- moduli e documentazione del repository RGS–QLG Quaternion Bridge.

Ogni futura voce deve contenere ID stabile, enunciato esatto, ipotesi/convenzioni, dimostrazione o derivazione, implementazione e commit SHA, test indipendente, residuo/tolleranza, classe epistemica, limiti e bibliografia. I risultati storici non vanno promossi silenziosamente a uno stato epistemico superiore.

**Regola editoriale:** “legge” richiede relazione definita e validata nel proprio dominio; “teorema” richiede dimostrazione. Le nuove affermazioni fisiche restano ipotesi finché non sono derivate e testate indipendentemente.


# L. RGS–Q6–SU(2)–Wilson bridge (Phase 19)

## L1. Finite symmetry compatibility
Let \(\pi_S\) be the permutation induced on the recovered 64-point seed by \(S\in O_h\). Define
\[
H=\{S\in O_h:\pi_S(E)=E,\ \pi_S(F_2)=F_2\}.
\]
The reproducible audit reports \(|O_h|=48\), \(|H|=16\), and \(|H\cap SO(3)|=8\). The full \(O_h\) action is therefore not an automorphism action of the fixed canonical Q6 2-complex. The proper subgroup admits SU(2) lifts; the improper subgroup requires a Pin-type extension if a spinorial lift is intended.
**Stato:** NUMERICAMENTE VERIFICATO for the recovered finite seed and the specified Q6 labeling. It is not a statement about all possible relabelings or continuum symmetries.

## L2. Finite lattice connection
\[
\mathcal A_{Q_6}=SU(2)^{192},\qquad
\mathcal G_{Q_6}=SU(2)^{64}.
\]
Each independent oriented edge carries one link \(U_{ij}\), with \(U_{ji}=U_{ij}^{\dagger}\), and gauge action \(U_{ij}\mapsto G_jU_{ij}G_i^\dagger\).
**Stato:** DEFINITIONE ESATTA of a finite lattice gauge configuration space. The RGS point coordinates determine the base labels and symmetry action, not a unique choice of the 192 link matrices.

## L3. Plaquette holonomy and Wilson action
For each of the 240 square faces,
\[
U_p=U_{v_3v_0}U_{v_2v_3}U_{v_1v_2}U_{v_0v_1},
\qquad
S_W=\beta\sum_{p=1}^{240}\left(1-\frac12\operatorname{ReTr}U_p\right).
\]
Gauge covariance of \(U_p\) and gauge invariance of \(S_W\) follow algebraically from the endpoint transformation law. The action is invariant under Q6 2-complex automorphisms; the proper eight-element subgroup also admits the tested lifted spin action.
**Stato:** exact finite identities, with numerical residual tests in src/rgs_q6_su2_wilson.py. The random seeded link field is a validation fixture, not a physically reconstructed connection.

## L4. Open bridge to a geometric connection
A physical map from the archived RGS vierbein/coframe to a preferred connection \(U_{ij}\) is not yet unique. It requires a specified local frame, a connection prescription (e.g. a discrete spin connection), link orientation, lift-sign convention and independent covariance tests. The Pin lift of improper elements and the continuum Yang–Mills limit are also open.
**Stato:** OPEN.


# M. Entanglement Q6 a due fogli e ponte QAT–olografia (Phase 20)

## M1. Stato heat-kernel PSU a due fogli

Per il Laplaciano combinatorio \(L_6\) del grafo Q6, definire
\[
K_\tau=e^{-\tau L_6/2},\qquad
C_\tau=K_\tau/\|K_\tau\|_F,\qquad
\rho_+=C_\tau C_\tau^\dagger.
\]
Lo stato bipartito è \(|\Psi_\tau\rangle=\sum_{xy}(C_\tau)_{xy}|x\rangle_+|y\rangle_-\).

**Stato:** DEFINIZIONE DI MODELLO per il kernel e la fattorizzazione a due fogli.

## M2. Spettro di Schmidt e entropia

Usando \(\lambda_j(L_6)=2j\), \(m_j=\binom6j\), si ha
\[
Z_6(\tau)=(1+e^{-2\tau})^6,\qquad
\operatorname{spec}(\rho_+):
\frac{e^{-2\tau j}}{Z_6(\tau)}
\quad\text{con molteplicità }\binom6j,
\]
e quindi
\[
S_+(\tau)=6h\!\left(\frac{1}{1+e^{2\tau}}\right).
\]
Per \(\tau=0\), \(S_+=\log64\); per \(\tau\to\infty\), \(S_+\to0\). Per lo stato puro bipartito \(I(+:-)=2S_+\).

**Stato:** IDENTITÀ FINITA ESATTA per il particolare stato heat-kernel. La corrispondenza con una geometria fisica non segue da questa identità.

## M3. Ponte QAT–Van Raamsdonk

La proposta di identificare i due fogli PSU–Anti-Twin con fattori quantistici correlati e la loro connessione con una geometria emergente è una **IPOTESI FISICA APERTA**. Il parametro \(\tau\) è adimensionale e non è ancora derivato dalla tetrade, dal cono nullo o dall'ovalizzazione QAT. La relazione con Van Raamsdonk richiede una mappa fra osservabili di entanglement e osservabili geometriche.

## M4. Area discreta e limite olografico

La relazione candidata è \(S(A)\simeq\operatorname{Area}(\gamma_A)/(4G_{\rm eff})\), con \(\gamma_A\) superficie minima discreta. Non sono ancora definiti in modo completo il complesso bulk, la nozione di area indotta dal coframe, il dominio delle partizioni e la costante \(G_{\rm eff}\).

**Stato:** CANDIDATO / APERTO. Non è lecito dichiarare verificata una legge di Ryu–Takayanagi sulla sola base dell'entropia Q6.

## M5. Raffinamento e correzioni

Swingle motiva l'organizzazione multiscala dell'informazione; Faulkner–Lewkowycz–Maldacena mostrano che, nei contesti olografici pertinenti, il termine di area può ricevere contributi di entanglement bulk. La loro implementazione su PRISM/Krist-All richiede tensori di rinormalizzazione, superfici e misure definiti esplicitamente.

**Stato:** PROGRAMMA DI DERIVAZIONE, non risultato stabilito.

## M6. Wormhole euclidei e PSU–Anti-Twin

Il confronto con Van Raamsdonk–Vilar López (2026) riguarda la costruzione dei saddles euclidei, il gluing di tubi AdS, le condizioni al bordo e la dominanza nell'integrale di cammino. L'identificazione con il throat PSU–Anti-Twin e con i wormhole EMS precedenti resta aperta.

**Stato:** CONFRONTO TEORICO DA SVILUPPARE; equivalenza fisica non stabilita.

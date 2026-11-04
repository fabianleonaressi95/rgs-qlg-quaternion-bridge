# RGS–PSU–Anti-Twin: architettura unificata candidata
## Condensato Q-ball, superfluidità, geometria a due fogli e transizione buco nero/buco bianco

**Stato:** architettura di ricerca; ipotesi fisica non dimostrata.  
**Obiettivo:** definire un programma variazionale e numerico comune per PSU/Anti-Twin, Q-ball, fluidi quantistici, elettromagnetismo e geometria RGS/TEGR, con un settore di orizzonte intrappolato/anti-intrappolato.

---

## 1. Tesi e disciplina epistemica

La tesi da sottoporre a prova è che PSU e Anti-Twin possano essere descritti come due fogli dinamicamente accoppiati; un condensato complesso localizzato di tipo Q-ball fornisca i gradi di libertà coerenti di carica e densità; il limite idrodinamico del condensato produca una descrizione superfluida; e, in un regime gravitazionale sufficientemente compatto, la geometria ammetta una regione intrappolata seguita eventualmente da una regione anti-intrappata.

Questa catena **non è una conseguenza automatica** delle singole teorie note. Ogni freccia deve essere dimostrata:

1. PSU/Anti-Twin -> campi e simmetrie a due fogli;
2. campo complesso -> Q-ball stazionario e carica conservata;
3. campo complesso -> idrodinamica di Madelung/superfluida;
4. stress-energia -> geometria RGS e orizzonte;
5. dinamica quantistica/effettiva -> eventuale transizione da orizzonte intrappolato a anti-intrappolato;
6. tutti i limiti -> osservabili e vincoli sperimentali.

Etichette obbligatorie:
- **PROVED:** dimostrazione matematica con ipotesi esplicite.
- **NUMERICALLY VERIFIED:** calcolo riproducibile, convergenza e residui documentati.
- **CANDIDATE:** modello definito, ma non ancora dimostrato.
- **PHYSICAL HYPOTHESIS:** interpretazione fisica da verificare.
- **OPEN:** derivazione, dato o test mancante.

Un Q-ball non è un buco nero per definizione. Un orizzonte acustico non è un orizzonte gravitazionale. Un buco bianco non è semplicemente un campo con carica opposta. Il termine Anti-Twin deve essere definito attraverso una trasformazione precisa, non usato come sinonimo di antimateria o di white hole.

## 2. Ontologia dei campi: due fogli, un condensato, una geometria

### 2.1 Geometria RGS

Il campo fondamentale geometrico è il vierbein
\[
e^a{}_{\mu}(x),\qquad
g_{\mu\nu}=\eta_{ab}e^a{}_{\mu}e^b{}_{\nu},\qquad
e=|\det e^a{}_{\mu}|.
\]
La connessione teleparallela \(\omega^a{}_{b\mu}\), la torsione
\[
T^a{}_{\mu\nu}=2\partial_{[\mu}e^a{}_{\nu]}
+2\omega^a{}_{b[\mu}e^b{}_{\nu]}
\]
e lo scalare TEGR devono usare una convenzione unica e dichiarata. In una formulazione covariante TEGR, la connessione di spin è piatta ma non necessariamente nulla; il gauge Weitzenböck \(\omega^a{}_{b\mu}=0\) è una scelta di gauge, non un'identità valida in ogni tetrad.

### 2.2 Doppio foglio PSU/Anti-Twin

Rappresentiamo il doppio foglio con un indice discreto \(s\in\{+,-\}\) e una trasformazione di scambio
\[
\mathsf X:(\Phi_+,\Phi_-,e_+,e_-)\mapsto
(\Phi_-,\Phi_+,e_-,e_+),\qquad \mathsf X^2=1.
\]
Questa è una **rappresentazione matematica candidata** dell'Anti-Twin. Non implica ancora due metriche indipendenti: il primo modello deve confrontare:

- **M0 — metrica comune:** un solo vierbein \(e^a{}_{\mu}\), due sezioni scalari \(\Phi_+,\Phi_-\);
- **M1 — geometrie correlate:** due vierbein \(e_+,e_-\), con un termine di interazione ammesso dalle simmetrie e con vincoli che evitino gradi di libertà fantasma;
- **M2 — doppio foglio discreto:** una geometria di base e una connessione discreta lungo il collegamento fra fogli.

Si parte da M0. M1 e M2 si introducono solo se M0 non descrive le osservabili e se esiste un'azione consistente. Il segno opposto della fase o della carica non basta a definire il foglio Anti-Twin.

### 2.3 Campo condensato

Il grado di libertà principale è un doppietto complesso
\[
\boldsymbol\Phi=(\Phi_+,\Phi_-)^T.
\]
Nel modello minimo si usa una sola connessione elettromagnetica diagonale \(A_\mu\):
\[
D_\mu\Phi_s=(\nabla_\mu-iqA_\mu/\hbar)\Phi_s.
\]
La simmetria globale di scambio \(\mathsf X\) è distinta dalla simmetria gauge U(1). Il settore relativo può avere una fase
\(\vartheta=\arg\Phi_+-\arg\Phi_-\), gauge-invariante rispetto alla U(1) diagonale.

Una candidata ansatz stazionaria è
\[
\Phi_s(t,r)=f_s(r)e^{i\omega_s t}.
\]
La soluzione simmetrica \(f_+=f_-\) e quella antisimmetrica \(f_+=f_-\) con fase relativa \(\pi\) vanno trattate come settori diversi e non imposte entrambe. La simmetria effettiva dipende dal potenziale e dalle condizioni al bordo.

## 3. Azione comune: forma minima e termini da derivare

Un punto di partenza, non ancora una TOE, è
\[
S=\int d^4x\,e\left[
\frac{M_{\rm Pl}^2}{2}T
-\frac14F_{\mu\nu}F^{\mu\nu}
-\sum_{s=\pm}g^{\mu\nu}(D_\mu\Phi_s)^*D_\nu\Phi_s
-V(\Phi_+,\Phi_-)
\right]+S_{\rm boundary}.
\]
Le normalizzazioni di \(\hbar,c,\mu_0\) vanno ripristinate in una convenzione SI o naturale coerente prima di confrontare numeri fisici.

Un potenziale minimale che rispetta U(1) diagonale e scambio \(\mathsf X\) può essere parametrizzato da invarianti:
\[
I_+=|\Phi_+|^2,\quad I_-=|\Phi_-|^2,\quad
I=I_++I_-,\quad J=|\Phi_+|^2-|\Phi_-|^2,
\]
\[
V=V_0(I)+\alpha J^2
-\kappa(\Phi_+^*\Phi_-+\Phi_-^*\Phi_+).
\]
Questa è una **base di operatori candidati**, non la scelta già giustificata dei coefficienti. Occorre:
- verificare dimensioni e stabilità inferiormente limitata;
- classificare gli operatori ammessi dalle simmetrie PSU/Anti-Twin effettivamente adottate;
- determinare i coefficienti da una teoria microscopica o da matching, non da fitting arbitrario;
- controllare il numero di modi, i vincoli e l'assenza di ghost;
- confrontare il limite a un foglio con il Q-ball globale già implementato.

La materia QHD non va aggiunta una seconda volta se è la stessa variabile condensata \(\Phi\). La descrizione Madelung è una riscrittura idrodinamica del medesimo campo, non automaticamente un secondo fluido indipendente.

## 4. Q-ball e Maxwell: settore già avviato nel repository

Per
\[
V(f)=m^2f^2-\lambda f^4+gf^6
\]
e \(\Phi=f(r)e^{i\omega t}\), il Q-ball globale soddisfa
\[
f''+\frac2r f'=(m^2-\omega^2)f-2\lambda f^3+3gf^5.
\]
Con il gauge elettromagnetico:
\[
\Omega(r)=\omega-qA_0(r),\qquad
f''+\frac2r f'=(m^2-\Omega^2)f-2\lambda f^3+3gf^5.
\]
L'equazione di Maxwell va derivata dalla stessa azione, con la stessa convenzione per \(D_\mu\), corrente e segno. In metrica curva la forma covariante è
\[
\nabla_\mu F^{\mu\nu}=j^\nu,
\qquad
\frac1{e}\partial_\mu(eF^{\mu\nu})=j^\nu.
\]
Il ramo gauged deve essere confrontato col limite \(q\to0\), e deve soddisfare la legge di Gauss tra flusso al bordo e carica volumica.

**Repository attuale:** il Q-ball globale stazionario è implementato in `src/rgs_psu_qball_bvp.py`; il ramo Maxwell è proposto nella PR #85 (`src/rgs_psu_qball_maxwell_bvp.py`). La PR #85 è ancora da integrare e la sua CI deve essere controllata prima di considerare il ramo Maxwell parte stabile di `main`.

## 5. Dal campo complesso al superfluido: derivazione, non analogia verbale

Scrivendo
\[
\Phi=\sqrt{n}\,e^{i\theta},\qquad
\mathbf v_s=\frac{\hbar\nabla\theta-q\mathbf A}{m},
\]
si ottengono, dove \(\Phi\ne0\), la continuità
\[
\partial_t n+\nabla\cdot(n\mathbf v_s)=0
\]
e l'equazione di fase
\[
\hbar\partial_t\theta+qA_0
+\frac12m\mathbf v_s^2+V_{\rm ext}
+\mu(n)+Q[n]=0,
\]
con potenziale quantistico
\[
Q[n]=-\frac{\hbar^2}{2m}\frac{\Delta\sqrt n}{\sqrt n}.
\]
Il termine \(\mu(n)\) proviene dal potenziale d'interazione. La vorticità è nulla localmente in regioni semplicemente connesse senza zeri del condensato; i vortici quantizzati si collocano dove la fase è singolare e \(\Phi=0\).

Per un vero superfluido a due componenti (fluido normale + superfluido) servono densità, velocità, entropia/temperatura e mutual friction secondo una chiusura tipo Landau/HVBK. Non va chiamato superfluido un singolo campo Madelung senza specificare regime termodinamico, eccitazioni e risposta lineare.

**Regola anti-doppio-conteggio:** se \(\Phi_+,\Phi_-\) sono il condensato PSU/Anti-Twin, QHD è la loro forma idrodinamica. Un fluido normale addizionale è ammesso solo se ha gradi di libertà distinti e un'azione/chiusura separata.

## 6. Buco nero e buco bianco: definizione geometrica operativa

In una geometria sferica dinamica si definiscono le espansioni delle due congruenze radiali nulle:
\[
\theta_+=\frac{2}{R}k_+^\mu\nabla_\mu R,
\qquad
\theta_- =\frac{2}{R}k_-^\mu\nabla_\mu R.
\]
Con una convenzione di orientazione dichiarata:
- regione normale: \(\theta_+>0,\theta_-<0\);
- regione intrappolata (black-hole trapped): \(\theta_+<0,\theta_-<0\);
- regione anti-intrappolata (white-hole/anti-trapped): \(\theta_+>0,\theta_->0\).

Un orizzonte marginale soddisfa una delle due espansioni nulla. I segni vanno verificati con i vettori nulli e il tempo futuro scelti; non si deducono dal segno della carica del Q-ball.

La prima domanda numerica è se il condensato genera una superficie intrappolata:
\[
\mathcal C(R)=\frac{2G M_{\rm MS}(R)}{c^2R}
\]
e se la condizione di intrappolamento viene raggiunta prima che la soluzione diventi singolare o instabile. \(M_{\rm MS}\) è la massa quasi-locale di Misner–Sharp. Per geometrie statiche la condizione \(2GM/(c^2R)=1\) è un indicatore di orizzonte, ma nelle geometrie dinamiche va usare direttamente \(\theta_\pm\).

La transizione black-to-white richiede una regione anti-intrappolata futura e una dinamica che la colleghi alla regione intrappolata. Non basta cambiare \(t\to-t\) in una soluzione statica. Si deve dimostrare:
1. regolarità della metrica e dei campi;
2. conservazione dei vincoli gravitazionali e gauge;
3. bilanci energia/carica attraverso il bounce;
4. causalità e struttura globale;
5. stabilità rispetto a perturbazioni non sferiche;
6. compatibilità con vincoli osservativi.

## 7. PSU/Anti-Twin come meccanismo di transizione: ipotesi testabile

Una proposta precisa, da falsificare, è che lo scambio \(\mathsf X\) tra i due fogli sia associato a una variabile d'ordine relativa \(J\) o \(\vartheta\), e che il regime di alta densità/curvatura induca una transizione dinamica tra un ramo dominato dal foglio \(+\) e un ramo dominato dal foglio \(-\).

Per diventare fisica, la proposta deve produrre:
- un termine d'azione determinato che controlli la transizione;
- un'equazione del moto per la variabile relativa;
- un criterio di soglia derivato (non inserito a mano);
- un tensore energia-impulso conservato;
- una soluzione regolare con superfici marginali e anti-intrappolate;
- una predizione quantitativa distinta da modelli già esistenti.

La soglia di De Aquino non è inserita nell'azione in questa fase: può entrare solo dopo audit dimensionale, confronto con relatività standard e derivazione indipendente.

## 8. Architettura software nel repository

Non creare un monolite. Separare fisica, geometria, solver, audit e confronto esterno.

| Layer | Modulo previsto | Responsabilità |
|---|---|---|
| Convenzioni | `src/physics_conventions.py` | segnatura, unità, normalizzazioni, cariche |
| PSU/Anti-Twin | `src/rgs_psu_antitwin.py` | scambio Z2, doppio foglio, invarianti relativi |
| Campi | `src/rgs_psu_qball_fields.py` | potenziale, correnti, tensore energia-impulso |
| Gauge | `src/rgs_psu_gauge_maxwell.py` | connessione U(1), Maxwell, Gauss |
| QHD | `src/rgs_psu_superfluid_madelung.py` | trasformazione Madelung, continuità, quantum potential |
| Geometria | `src/rgs_psu_teleparallel_spherical.py` | tetrad sferico, torsione, scalare TEGR |
| Gravità | `src/rgs_psu_spherical_gravity.py` | vincoli/evoluzione sferica e massa Misner–Sharp |
| Orizzonti | `src/rgs_psu_horizon_diagnostics.py` | espansioni nulle, trapped/anti-trapped surfaces |
| BVP | `src/rgs_psu_qball_maxwell_bvp.py` | ramo stazionario piatto gauged |
| Evoluzione | `src/rgs_psu_antitwin_bounce.py` | evoluzione dinamica e test di bounce |
| Audit | `src/rgs_psu_unified_compact_object_audit.py` | residui, convergenza, bilanci, report JSON |

Questa tabella è una roadmap, non l'affermazione che tutti i moduli esistano già. Riutilizzare i moduli esistenti e non duplicare implementazioni QHD, fluxon o vierbein.

### Dipendenze fra layer

\[
\text{Q6/PRISM discretization}
\to \text{fields + U(1)}
\to \text{Q-ball / Madelung}
\to T_{\mu\nu}
\to \text{vierbein/TEGR}
\to \theta_\pm
\to \text{bounce dynamics}.
\]
Il grafo Q6/PRISM può fornire una discretizzazione e un settore di trasporto; non dimostra da solo il continuo 3+1 né una struttura microscopica della materia.

## 9. Roadmap di implementazione e gate di accettazione

### Fase A — Congelare convenzioni e verificare PR #85
- controllare equazioni, segni Maxwell, normalizzazione di carica e limite \(q=0\);
- CI completa e test indipendenti;
- audit del flusso di energia e legge di Gauss;
- nessuna estensione gravitazionale finché questo gate non è verde.

### Fase B — Doppio foglio senza gravità
- implementare lo scambio Z2 e gli invarianti;
- costruire potenziale limitato inferiormente;
- risolvere Q-ball a due componenti;
- calcolare matrice Hessiana del potenziale e spettro delle perturbazioni;
- recuperare il Q-ball a un foglio e il limite di disaccoppiamento.

### Fase C — Superfluido derivato
- derivare Madelung direttamente dal campo complesso;
- confrontare PDE del campo e PDE idrodinamiche;
- testare continuità, energia e carica;
- trattare vortici/fluxon come difetti topologici con condizioni al contorno;
- aggiungere il fluido normale solo con una teoria a due fluidi esplicita.

### Fase D — Gravità statica e compattezza
- accoppiare il tensore energia-impulso del condensato alle equazioni del vierbein;
- risolvere branch statici regolari;
- misurare massa, raggio, carica, compattezza, redshift e stabilità radiale;
- confrontare con boson stars e gauged boson stars.

### Fase E — Dinamica collasso/bounce
- evolvere dati iniziali vincolati;
- misurare \(\theta_\pm\) e massa quasi-locale;
- cercare trapped e anti-trapped surfaces;
- testare convergenza spaziale/temporale e stabilità non sferica;
- confrontare il bounce con modelli black-to-white-hole esistenti.

### Fase F — PSU/PRISM e limite continuo
- discretizzare l'azione e i vincoli sul complesso;
- dimostrare convergenza degli operatori e delle soluzioni;
- separare gli effetti della discretizzazione da quelli della nuova fisica;
- nessuna affermazione Mosco senza liminf, recovery sequence e compattezza.

### Gate numerici minimi
- residui delle equazioni scalare, Maxwell e gravità;
- vincoli Hamiltoniano e momento;
- identità di Noether e legge di Gauss;
- positività del potenziale e controllo dei modi fantasma;
- convergenza di \(\omega,E,Q,M,R\) al raffinamento;
- bilancio energetico/carica durante il bounce;
- indipendenza dalla scelta di gauge e coordinate;
- stabilità a perturbazioni radiali e non radiali.

## 10. Repositories esterni: benchmark, non dipendenze

| Repository | Uso corretto nel progetto | Limite |
|---|---|---|
| [mkinach/Q-ball](https://github.com/mkinach/Q-ball) | confronto per solver gauged Q-ball, asse-simmetrico e 3D | verificare licenza, unità e convenzioni prima di riusare codice |
| [GRTLCollaboration/ExoZvezda](https://github.com/GRTLCollaboration/ExoZvezda) | benchmark di evoluzione di boson stars e oggetti compatti esotici | si appoggia a GRChombo; non è un solver PSU/TEGR |
| [GRTLCollaboration/engrenage](https://github.com/GRTLCollaboration/engrenage) | studio didattico/benchmark BSSN sferico, scalar field, black hole e oscillaton | non sostituisce un codice NR validato per risultati di precisione |
| [ChenSun-Phys/boson_stars](https://github.com/ChenSun-Phys/boson_stars) | confronto Einstein–Klein–Gordon per condensati scalari | modello BEC a due scalari; non contiene il nostro doppio foglio |
| [vagiedd/BosonStars](https://github.com/vagiedd/BosonStars) | confronto metodo di rilassamento per boson stars | verificare licenza e ipotesi prima dell'integrazione |

Non importare direttamente repository esterni nel package. Prima produrre una scheda per ciascuno: licenza, commit/tag, equazioni, gauge, unità, test, benchmark pubblicati, dipendenze e differenze rispetto alle convenzioni RGS. Implementare adapter solo quando una quantità fisica è confrontabile.

## 11. Bibliografia essenziale annotata

1. S. Coleman, “Q-balls,” *Nuclear Physics B* **262** (1985), 263–283. DOI: [10.1016/0550-3213(85)90286-X](https://doi.org/10.1016/0550-3213(85)90286-X). Teorema classico di esistenza per solitoni non topologici con carica globale conservata.
2. C. Barceló, S. Liberati, M. Visser, “Analogue Gravity,” *Living Reviews in Relativity* **14** (2011), 3. DOI: [10.12942/lrr-2011-3](https://doi.org/10.12942/lrr-2011-3). Quadro generale delle metriche efficaci e degli orizzonti analoghi in sistemi di materia condensata.
3. S. K. Manikandan, A. N. Jordan, “Bosons falling into a black hole: A superfluid analogue,” *Physical Review D* **98** (2018), 124043. DOI: [10.1103/PhysRevD.98.124043](https://doi.org/10.1103/PhysRevD.98.124043). Proposta di analogia BEC-superfluido per la fisica quantistica del buco nero; non è una soluzione di Einstein.
4. M. Han, C. Rovelli, F. Soltani, “Geometry of the black-to-white hole transition within a single asymptotic region,” *Physical Review D* **107** (2023), 064011. DOI: [10.1103/PhysRevD.107.064011](https://doi.org/10.1103/PhysRevD.107.064011). Geometria Lorentziana regolare fuori dalla regione di tunnelling, con transizione tra orizzonti.
5. A. Rignon-Bret, C. Rovelli, “Black to white transition of a charged black hole,” *Physical Review D* **105** (2022), 086003. DOI: [10.1103/PhysRevD.105.086003](https://doi.org/10.1103/PhysRevD.105.086003). Caso carico utile per il confronto Maxwell–geometria.
6. J. W. Maluf, “The teleparallel equivalent of general relativity,” *Annalen der Physik* **525** (2013), 339–357. DOI: [10.1002/andp.201200272](https://doi.org/10.1002/andp.201200272). Formulazione tetradica/torsionale e equivalenza delle equazioni di campo con GR.
7. C. Rovelli, F. Vidotto, “Planck stars,” arXiv:1401.6562 (2014). [arXiv](https://arxiv.org/abs/1401.6562). Scenario di rimbalzo quantistico; ipotesi di gravità quantistica, non evidenza osservativa.
8. M. Visser, C. Barceló, S. Liberati, “Analogue gravity and the Hawking effect: historical perspective and literature review,” (review, 2023). Da usare come aggiornamento del quadro analogue-gravity e delle sue limitazioni.

## 12. Primo obiettivo scientifico concreto

Non cercare subito una simulazione di “tutto”. Il primo risultato pubblicabile e falsificabile deve essere:

**Esistenza e stabilità di una famiglia di Q-ball a due fogli PSU/Anti-Twin con accoppiamento Maxwell, e derivazione esatta del limite superfluido.**

Deliverable:
1. azione a due fogli completamente specificata;
2. equazioni Euler–Lagrange e correnti;
3. ramo stazionario numerico con continuazione in carica/coupling;
4. spettro completo delle perturbazioni lineari;
5. confronto con Q-ball globali e gauged Q-ball standard;
6. audit di conservazione, convergenza e stabilità.

Solo se questa fase supera i test si passa al coupling gravitazionale e alla ricerca di trapped/anti-trapped surfaces. Il risultato decisivo per la tesi black-to-white non è una figura di un toro o una densità centrale elevata: è una soluzione dinamica regolare, con vincoli soddisfatti, orizzonte intrappolato e successiva regione anti-intrappolata, stabile e compatibile con i limiti osservativi.

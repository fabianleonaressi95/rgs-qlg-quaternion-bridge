# RGS–PSU–PRISM–QLG
## Teoria geometrica unificata: fondamenti, azione e limite continuo

**Versione:** 0.1 — 4 novembre 2026  
**Stato:** documento di sintesi e programma di derivazione; non è una dichiarazione di teoria fisica già dimostrata.  
**Repository:** [rgs-qlg-quaternion-bridge](https://github.com/fabianleonaressi95/rgs-qlg-quaternion-bridge)

---

## Abstract

Questo documento organizza in un'unica architettura i filoni RGS, Q6, PRISM–Krist-All, quaternioni, gruppo binario icosaedrico (2I), algebra di Clifford, QLG, vierbein/teleparallelismo, PSU–Anti-Twin, QAT, operatori di Dirac e modelli gauge. La tesi di lavoro è che una geometria discreta dotata di coframe, connessione interna e struttura a due fogli possa fornire un'origine comune per geometria, campi di gauge e settori quantistici.

La tesi è formulata come **programma matematico falsificabile**. Le identità finite dimostrate, i risultati numerici riproducibili, le costruzioni candidate e le ipotesi fisiche sono separati. In particolare, non si assume che l'entanglement sia già geometria, che un tetrade statico dimostri la dinamica gravitazionale, o che una simulazione gauge finita dimostri il limite continuo della QCD.

## 1. Tesi centrale e gerarchia

La struttura proposta è

\[
\boxed{
\text{complesso discreto}
\to (e,\omega,U)
\to \text{Clifford/spinori}
\to \text{dinamica gauge e geometrica}
\to \text{PSU/QAT}
\to \text{limite continuo}
}
\]

con i seguenti ruoli:

1. **RGS**: geometria discreta, coordinate, coframe e vincoli geometrici.
2. **Q6**: complesso combinatorio locale e algebra spettrale di riferimento.
3. **PRISM**: struttura prodotto e famiglia di raffinamento/discretizzazione.
4. **Krist-All**: organizzazione algebrica dei settori discreti.
5. **Quaternioni e (2I)**: rappresentazioni rotazionali e spinoriali finite.
6. **Clifford/QLG**: spinori, trasporto e operatori geometrici.
7. **Vierbein e connessioni**: metrica, torsione e gravità.
8. **PSU–Anti-Twin/QAT**: doppio foglio e dinamica quantistica candidata.
9. **Yang–Mills/QCD**: connessioni interne, olonomie e materia.

Questa gerarchia è una proposta di dipendenze matematiche. Nessun livello deve essere identificato con il successivo senza una mappa esplicita.

## 2. Convenzioni e oggetti fondamentali

### 2.1 Complesso discreto

Sia (K_h) un complesso cellulare orientato, con spazi di cochain (C^p(K_h;\mathbb R)) e operatori di incidenza

\[
B_p:C^p(K_h)\to C^{p+1}(K_h),\qquad B_{p+1}B_p=0.
\]

Laplaciano combinatorio di grado zero:

\[
L_0=B_0^T B_0.
\]

Per il grafo Q6 a 64 vertici e 192 archi, il grado è 6 e

\[
\operatorname{spec}(L_0)=\{2j\text{ con molteplicità }\binom 6j: j=0,\ldots,6\}.
\]

Il numero ciclico del grafo senza 2-celle è \(\beta_1=192-64+1=129\). Quando si includono le 240 facce quadrate nel 2-scheletro, gli invarianti cambiano: il dato documentato è \(b_1=0, b_2=111\). Non confondere il rango ciclico del grafo con l'omologia del complesso cellulare.

### 2.2 Prodotto PRISM

La famiglia di riferimento è

\[
G_{\rm PRISM}=C_{64}\square C_6^3\square C_5,
\qquad |V|=64\cdot6^3\cdot5=69120.
\]

L'indice dichiarato per le coordinate prodotto è

\[
i=((((G\cdot6+L)\cdot6+C)\cdot6+T)\cdot5+B).
\]

Le proprietà riportate — grado 10, 1920 stati attivi, 5760 correzioni, 108 cicli di lunghezza 640, numeri di Betti e gap spettrale — devono essere associate alla precisa versione del grafo, alle convenzioni di bordo e alla matrice utilizzata. Non vanno trasferite a un altro operatore per analogia.

### 2.3 Geometria RGS e simmetrie

I punti discreti (X_i\in\mathbb R^3), le incidenze e l'azione di gruppo devono essere archiviati come un unico oggetto versionato. Distinguere:

- simmetria dell'insieme di coordinate;
- automorfismo del grafo;
- automorfismo del complesso con facce;
- simmetria dell'operatore tangente;
- simmetria della connessione e dell'azione.

Il gruppo cubico (O_h) ha ordine 48, ma non ogni sua azione sulle coordinate è automaticamente un automorfismo del complesso Q6 etichettato. I documenti di equivarianza riportano sottogruppi conservanti distinti; ogni uso fisico deve specificare il gruppo effettivo.

## 3. Simmetrie quaternioniche, (2I) e Clifford

I quaternioni unitari realizzano (SU(2)). Il doppio ricoprimento

\[
SU(2)\longrightarrow SO(3)
\]

identifica (q) e (-q) come rotazioni tridimensionali, ma li mantiene distinti nella rappresentazione spinoriale.

Il gruppo binario icosaedrico

\[
2I\subset SU(2),\qquad |2I|=120
\]

è il lift spinoriale del gruppo icosaedrico rotazionale. Il quoziente (S^3/2I) è un candidato topologico naturale per il settore icosaedrico, ma il suo impiego come spazio fisico richiede una metrica, una struttura spin e condizioni globali esplicite.

Per la firma lorentziana adottata \(\eta=\operatorname{diag}(-,+,+,+)\), le matrici gamma soddisfano

\[
\{\gamma^A,\gamma^B\}=2\eta^{AB}I.
\]

Quaternioni, (2I) e Clifford hanno ruoli distinti: il primo codifica rotazioni, il secondo una simmetria finita spinoriale, il terzo la struttura algebrica necessaria per gli spinori relativistici. Il passaggio a una connessione fisica deve essere costruito attraverso un lift di gruppo e una legge di trasporto.

**Obblighi aperti:** scelta coerente dei lift globali, trattamento delle simmetrie improprie mediante Pin, verifica delle relazioni di rappresentazione e costruzione di una connessione compatibile con il coframe.

## 4. Coframe, metrica e connessioni

Sia (e^a=e^a{}_{\mu}dx^\mu) un coframe invertibile. La metrica è

\[
g_{\mu\nu}=e^a{}_{\mu}e^b{}_{\nu}\eta_{ab},
\qquad
\Gamma^\mu=e_a{}^\mu\gamma^a.
\]

Ne segue

\[
\{\Gamma^\mu,\Gamma^\nu\}=2g^{\mu\nu}I.
\]

Questa è un'identità algebrica, subordinata all'invertibilità del coframe e alla firma scelta.

Per una connessione di spin \(\omega^a{}_b\), la torsione e la curvatura sono

\[
T^a=de^a+\omega^a{}_b\wedge e^b,
\qquad
R^a{}_b=d\omega^a{}_b+\omega^a{}_c\wedge\omega^c{}_b.
\]

La scelta \(\omega=0\) è il gauge di Weitzenböck in una formulazione teleparallela; non è una prescrizione universale indipendente dal tetrade. La formulazione TEGR richiede la connessione inerziale appropriata e l'identità tra scalare di Ricci e scalare torsionale, con il termine di divergenza di bordo nella convenzione adottata.

La presenza di un tetrade statico o di invarianti di torsione calcolati non dimostra da sola le equazioni di Einstein/TEGR. Occorrono azione, variazione, condizioni al bordo e verifica delle identità.

## 5. Connessione interna e settore Yang–Mills

Per un gruppo gauge (G\), una connessione discreta è assegnata ai link:

\[
U_{ij}\in G,\qquad U_{ji}=U_{ij}^{-1}.
\]

Per un plaquette orientato (p=(i_1i_2\ldots i_ni_1)), l'olonomia è

\[
U_p=U_{i_1i_2}U_{i_2i_3}\cdots U_{i_ni_1}.
\]

Una possibile azione di Wilson è

\[
S_{\rm YM}^{(h)}=
\beta\sum_p w_p\left(1-\frac{1}{d_R}
\operatorname{ReTr}_R U_p\right),
\]

dove la rappresentazione (R), la dimensione (d_R), i pesi (w_p), la normalizzazione e la relazione tra \(\beta\) e l'accoppiamento devono essere fissati. I pesi non possono essere introdotti ad hoc per ottenere un risultato desiderato.

Per collegare QLG alla geometria RGS bisogna derivare come \(U_{ij}\) dipenda dalla connessione e dal trasporto tra coframe adiacenti. La compatibilità con le simmetrie locali e la trasformazione di gauge deve essere dimostrata.

## 6. Operatori fermionici e Dirac

L'operatore continuo atteso ha forma

\[
D_{e,U}=\gamma^a e_a{}^\mu
\left(\partial_\mu+\Omega_\mu(e,\omega)+A_\mu\right),
\]

con connessione spinoriale \(\Omega_\mu\) e connessione interna \(A_\mu\). La discretizzazione deve convergere a questo operatore e rispettare le simmetrie appropriate.

Occorre mantenere distinti:

- Dirac spinoriale;
- Dirac–Kähler sulle cochain;
- Wilson-Dirac;
- overlap/Ginsparg–Wilson;
- operatori efficaci costruiti su sottospazi spettrali.

Le identità finite di un operatore non certificano automaticamente la convergenza all'operatore continuo né l'assenza di specie spurie. Devono essere testati hermiticità, chiralità, località, gap, indice e dipendenza dalla discretizzazione.

## 7. PSU–Anti-Twin e QAT: modello a due fogli

Si considera una decomposizione candidata

\[
\mathcal H_h=\mathcal H_{+,h}\otimes\mathcal H_{-,h}.
\]

Per Q6, ponendo \(K_\tau=e^{-\tau L_0/2}\), si può definire lo stato normalizzato

\[
|\Psi_\tau\rangle=
\frac{1}{\|K_\tau\|_F}
\sum_{x,y}(K_\tau)_{xy}|x\rangle_+|y\rangle_-.
\]

La ridotta è

\[
\rho_+(\tau)=
\frac{e^{-\tau L_0}}{\operatorname{Tr}(e^{-\tau L_0})}.
\]

Per lo spettro Q6 segue l'identità finita

\[
S_+(\tau)=6h\left(\frac{1}{1+e^{2\tau}}\right),
\]

dove \(h(p)=-p\log p-(1-p)\log(1-p)\). Questo risultato descrive l'entropia del modello di stato costruito dal semigruppo del Laplaciano. Non identifica \(\tau\) con il tempo fisico e non dimostra una legge gravitazionale di area.

Il programma QAT deve derivare, non postulare, una relazione del tipo

\[
\tau=f(\mathcal I[e,\omega]),
\]

dove \(\mathcal I\) è un invariante geometrico o dinamico ben definito. La geometria cubo–sfera–cono nullo, il doppio foglio e la riconnessione devono essere codificati nello stesso spazio delle configurazioni e nella stessa azione.

L'eventuale relazione con geometria emergente, ER=EPR o formule di area richiede una definizione indipendente di sottoregione, superficie estremale, area e costante gravitazionale efficace.

## 8. Azione unificata candidata

La struttura minima da derivare è

\[
S_h=S_{\rm geom}+S_{\rm gauge}+S_{\rm ferm}
+S_{\rm PSU}+S_{\rm matter}.
\]

Una decomposizione più esplicita è

\[
\begin{aligned}
S_{\rm geom}&=S_{\rm TEGR}[e,\omega]
\quad\text{o, in una formulazione equivalente, }S_{\rm EH}[g],\\
S_{\rm gauge}&=S_{\rm YM}^{(h)}[U],\\
S_{\rm ferm}&=\langle\bar\psi,D_{e,U}^{(h)}\psi\rangle,\\
S_{\rm PSU}&=S_+[e,\Psi_+]+S_-[e,\Psi_-]
+S_{\rm int}[e,\Psi_+,\Psi_-],\\
S_{\rm matter}&=S_{\rm matter}[e,\Phi].
\end{aligned}
\]

Questa è una **ansatz di architettura**, non un'azione già dedotta dai dati del repository. Ogni termine deve essere definito sul medesimo complesso, con le stesse convenzioni di orientazione, unità e normalizzazione.

Le equazioni discrete di campo sono

\[
\frac{\delta S_h}{\delta e}=0,\quad
\frac{\delta S_h}{\delta\omega}=0,\quad
\frac{\delta S_h}{\delta U}=0,\quad
\frac{\delta S_h}{\delta\bar\psi}=0,\quad
\frac{\delta S_h}{\delta\Psi_\pm}=0.
\]

Il test di unificazione consiste nel derivare queste equazioni da una singola azione e verificare identità di Noether, vincoli, conservazione e covarianza. La semplice somma di modelli indipendenti non soddisfa questo criterio.

## 9. Limite continuo: programma Mosco e convergenza

Sia \(h\to0\) una successione di complessi \(K_h\), con spazi di Hilbert \(H_h\), embedding in uno spazio limite e forme quadratiche \(q_h\). La convergenza Mosco richiede:

1. **liminf:** per ogni successione debolmente convergente \(u_h\rightharpoonup u\), \(q(u)\leq\liminf q_h(u_h)\);
2. **recovery:** per ogni \(u\) esiste una successione fortemente convergente \(u_h\to u\) con \(q(u)\geq\limsup q_h(u_h)\).

La prova sul toro piatto P1 è un risultato per quella famiglia e quelle ipotesi. Non implica la convergenza della geometria curva RGS. Per quest'ultima servono controllo uniforme di metriche, shape regularity, consistenza delle forme, controllo dei bordi e recovery sequence.

Dopo la convergenza scalare occorre dimostrare separatamente quella degli operatori di Dirac, delle connessioni gauge e delle forme accoppiate. Il traguardo è un limite coerente verso Laplace–Beltrami, Dirac e Yang–Mills sul medesimo spazio-tempo.

## 10. Filone aureo e raffinamento

Si mantengono distinti:

- \(\phi=(1+\sqrt5)/2\), rapporto aureo;
- \(\psi\), radice reale positiva di \(x^3-x^2-1=0\);
- geometria icosaedrica e gruppo \(2I\);
- raffinamenti PRISM e classi Krist-All;
- eventuale struttura modulare associata alla funzione eta di Dedekind.

La relazione con multifrattalità deve essere dimostrata mediante una famiglia di misure/operatore di trasferimento e limiti di scala. La comparsa di \(\phi\), \(\psi\) o di rapporti quasi aurei in uno spettro finito non dimostra da sola una legge universale.

Una possibile agenda è costruire \(\mathcal T_n\) su livelli di raffinamento, definire funzioni di partizione \(Z_n\), misure \(\mu_n\), dimensioni generalizzate \(D_q\), e testare convergenza e stabilità rispetto a mesh, condizioni al bordo e perturbazioni. Un legame con \(\eta(\tau)\) è significativo solo se emerge da una funzione di partizione o da una struttura modulare derivata.

## 11. Registro epistemico

### Identità matematiche finite (dimostrate dalle definizioni)

- spettro del Laplaciano del grafo ipercubico Q6;
- identità di catena \(B_{p+1}B_p=0\) quando le matrici d'incidenza sono costruite coerentemente;
- doppio ricoprimento \(SU(2)\to SO(3)\);
- identità Clifford curva conseguente a un coframe invertibile;
- formula entropica Q6 per lo stato definito nel §7.

### Risultati numerici riportati, da mantenere riproducibili

- ortogonalità della base tangente 188D con residuo circa \(1.31\times10^{-14}\);
- spettro di \(K_{\rm tangent}\): 29 modi negativi, 3 quasi-zero/rotazionali, 156 positivi nella convenzione documentata;
- isotropia numerica della restrizione \(J_0^T K_{\rm tangent}^{+}J_0\) nel tripletto selezionato;
- risultati PRISM, Krist-All e audit di torsione riportati nei rispettivi report.

Ogni numero deve essere collegato a commit, input hash, comando, versione delle dipendenze e tolleranza. Un residuo macchina non è una prova esatta.

### Costruzioni candidate

- coframe RGS come geometria discreta di base;
- connessione QLG derivata da coframe e trasporto;
- azione unificata del §8;
- PSU–Anti-Twin come struttura fisica a due fogli;
- relazione QAT tra geometria e parametro entropico.

### Ipotesi fisiche non dimostrate

- emergenza della gravità e dello spazio-tempo dalla struttura discreta;
- interpretazione fisica di PSU come unità fondamentale;
- identificazione dell'entanglement con connessione geometrica;
- origine della materia subatomica da una struttura cristallina;
- identificazione del tempo fisico con una velocità di riconnessione QAT;
- predizioni sperimentali univoche da \(\phi\), \(\psi\), Krist-All o De Aquino.

Queste ipotesi rimangono aperte fino alla derivazione di osservabili, scale dimensionali e test discriminanti.

## 12. Criteri di falsificazione e validazione

La teoria deve superare almeno questi controlli:

1. **Combinatoria:** ricostruzione indipendente di vertici, archi, facce, orientazioni e omologia.
2. **Simmetria:** verifica esatta dei gruppi che preservano coordinate, complesso, coframe e azione.
3. **Geometria:** positività/invertibilità del coframe, firma corretta, trasformazioni locali e consistenza delle connessioni.
4. **Gauge:** covarianza locale, identità di Bianchi discrete, indipendenza dal gauge degli osservabili.
5. **Fermioni:** spettro, indice, chiralità, località e assenza di modi spurii nel limite.
6. **Doppio foglio:** positività della densità ridotta, normalizzazione, entropia e dipendenza controllata dai parametri.
7. **Dinamica:** derivazione variazionale, vincoli e identità di Noether.
8. **Continuo:** convergenza Mosco/risolvente e stabilità spettrale su raffinamenti curvi.
9. **Fenomenologia:** unità fisiche, parametri indipendenti, incertezze e predizioni distinguibili dai modelli standard.

Un fallimento persistente di covarianza, positività, convergenza o confronto osservativo falsifica la specifica costruzione, non si corregge rinominando l'operatore.

## 13. Roadmap di dimostrazione

**Fase A — congelamento delle definizioni.** Un unico manifest versionato per coordinate, complessi, gruppi, coframe, operatori, unità e hash.

**Fase B — geometria e trasporto.** Derivare connessione di spin e link gauge dalla stessa struttura geometrica; testare trasformazioni locali e lift.

**Fase C — azione comune.** Definire tutti i termini sullo stesso complesso e derivare le equazioni di campo.

**Fase D — analisi spettrale.** Verificare operatori autoaggiunti, kernel, gap, stabilità e convergenza.

**Fase E — PSU/QAT.** Derivare l'accoppiamento a due fogli e la relazione geometria-entropia senza parametri arbitrari.

**Fase F — limite continuo e fenomenologia.** Dimostrare convergenza e soltanto dopo estrarre costanti fisiche e previsioni.

## 14. Documenti del repository da considerare fonti specialistiche

- [Q6 complex formalization](Q6_COMPLEX_FORMALIZATION.md)
- [Q6/O_h equivariance](RGS_Q6_OH_EQUIVARIANCE.md)
- [QAT–PSU–Q6 entanglement and emergent geometry](QAT_PSU_Q6_ENTANGLEMENT_EMERGENT_GEOMETRY.md)
- [Mosco curved refinement program](RGS_MOSCO_CURVED_REFINEMENT_PROGRAM.md)
- [Flat torus G6 P1 Mosco proof](FLAT_TORUS_G6_P1_MOSCO_PROOF.md)
- [Vierbein reconstruction status](VIERBEIN_RECONSTRUCTION_STATUS.md)
- [Master catalogue](MASTER_CATALOGUE_LAWS_THEOREMS.md)

I documenti specialistici mantengono i dettagli e gli output. Questo documento ne stabilisce l'architettura, le dipendenze e i limiti epistemici; non li sostituisce.

---

## Conclusione

La formulazione unificata da perseguire è una teoria discreta geometrica con coframe, connessione interna, spinori e doppio foglio, governata da una singola azione. Q6 e PRISM forniscono le strutture combinatorie; RGS fornisce la geometria; quaternioni, (2I) e Clifford organizzano le simmetrie; QLG descrive il trasporto; PSU–Anti-Twin/QAT costituisce il settore quantistico candidato. La gravità, i campi gauge e la materia potranno dirsi unificati soltanto quando emergeranno dalle stesse equazioni e sarà dimostrato il limite continuo.

**Stato finale del documento:** architettura teorica proposta, con nuclei matematici finiti e numerici identificati; azione comune, limite continuo curvo e interpretazione fisica ancora da dimostrare.


## 15. Compatibility register and local spin connection

The implementation-level dependency register is maintained in
`docs/PHYSICS_ARCHITECTURE_COMPATIBILITY.md` and
`src/rgs_physics_architecture.py`. It checks the common metric, Clifford
representation, gamma5 convention and Q6 chain complex, while preserving
separate status labels for RGS, Q6/PRISM, 2I, EC, TEGR, U(1), SU(2), SU(3),
QLG, QAT and continuum mechanics.

A local torsion-free Levi-Civita spin-connection kernel is now specified in
`src/rgs_spin_connection.py`. Given a coframe first jet, it computes the
Christoffel symbols, Lorentz spin connection, spinor connection and
orthonormal Dirac derivative. Its tests include a nontrivial spherical
coordinate frame in flat Minkowski spacetime.

This is an important local bridge, but not yet the RGS spin connection:
coframe derivatives on the curved RGS refinement must still be constructed.
The teleparallel inertial connection is a separate problem, and the
Einstein–Cartan connection remains an independent field before its algebraic
torsion equation is imposed. Internal gauge connections are also separate
from the Lorentz spin connection.

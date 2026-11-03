# Inserto per la monografia unitaria RGS–PSU–PRISM–QLG
## Fase 18 — azione di O_h sulla geometria ricostruita e sul complesso Q6

### Collocazione editoriale

Questo inserto NON costituisce un capitolo autonomo né una nuova linea teorica. È un raccordo tecnico fra:
1. il settore della ricostruzione geometrica RGS dei 64 punti;
2. il settore combinatorio Q6 (64 vertici, 192 archi, 240 quadrati);
3. il settore delle simmetrie finite e delle rappresentazioni;
4. il successivo settore QLG/SU(2), connessione, olonomie e dinamica di gauge.

Il numero definitivo del capitolo e i rinvii interni devono seguire l'indice canonico della monografia in 27 capitoli; non vengono assegnati numeri provvisori per non alterare la struttura editoriale già stabilita.

### Ruolo logico nella catena unitaria

La fase stabilisce un controllo di compatibilità che precede la costruzione dinamica:
\[
X_{O_h}\xrightarrow{\,O_h\,}\mathrm{Perm}(64),
\qquad
Q_6^{(2)}=(V,E,F_2),
\qquad
\pi_S:V\to V.
\]
Per ogni matrice signed-permutation \(S\in O_h\), la mappa sulle coordinate induce una permutazione \(\pi_S\) dei 64 punti. La verifica di
\[
\pi_{ST}=\pi_S\circ\pi_T
\]
controlla la legge di gruppo. Le condizioni
\[
\pi_S(E)=E,\qquad \pi_S(F_2)=F_2
\]
sono test distinti: la simmetria del cloud geometrico non implica automaticamente che l'etichettatura Q6 sia equivarante.

### Esiti e classificazione

- **Da certificare numericamente con l'audit riproducibile:** unicità e biiettività delle permutazioni indotte sui 64 punti; fedeltà e composizione dell'azione; inversione centrale come antipode.
- **Gate indipendente:** conservazione dei 192 archi e delle 240 facce per tutti i 48 elementi. Un esito negativo non invalida l'azione geometrica \(O_h\); segnala che la specifica identificazione dei label Q6 non è equivarante.
- **Non dimostrato da questa fase:** covarianza delle variabili di link SU(2), delle olonomie di plaquette, dell'azione di Wilson o della dinamica QLG.
- **Passo successivo:** definire l'azione indotta sui link orientati, la trasformazione dei link sotto \(\pi_S\), e verificare la covarianza delle olonomie prima di formulare un accoppiamento dinamico.

### Collegamento con il resto della monografia

La Fase 18 non sostituisce né duplica:
- la costruzione del complesso Q6 e la sua omologia cellulare;
- la ricostruzione RGS e i proiettori tangenti;
- il sollevamento \(SU(2)\to SO(3)\) e la connessione di fibra;
- il settore teleparallelo con vierbein e torsione;
- i limiti di raffinamento/Mosco verso la geometria continua.

Essa verifica una condizione strutturale intermedia: che l'azione di simmetria geometrica e la codifica combinatoria possano essere messe in relazione in modo esplicito e falsificabile. Le mappe tra questi settori devono restare dichiarate, con domini, convenzioni e residui.

### Regola editoriale permanente

Tutti i moduli, audit, script, matrici e risultati sono parti della medesima ricerca e vanno ricondotti ai 27 capitoli della monografia. Le fasi del repository sono milestone di verifica e sviluppo, non capitoli autonomi. Non si cambia l'indice canonico né si inventano titoli o numerazioni per accomodare un nuovo risultato.


## Fase 19 — fibrato di connessione Q6, olonomie e azione di Wilson

Questa fase continua lo stesso raccordo senza creare un capitolo separato. Una volta fissata la corrispondenza fra i 64 vertici RGS e Q6, il complesso determina l'indice delle variabili di connessione:
\[
\mathcal A_{Q_6}=SU(2)^{E},\quad |E|=192,
\qquad
\mathcal G_{Q_6}=SU(2)^{V},\quad |V|=64.
\]
Per ciascuna delle 240 facce quadrate si definisce l'olonomia ordinata \(U_p\), e quindi
\[
S_W[U]=\beta\sum_{p\in F_2}
\left(1-\frac12\operatorname{Re}\operatorname{Tr}U_p\right).
\]

Il controllo di equivarianza della Fase 18 è ora quantitativo: il sottogruppo \(H\subset O_h\) che preserva archi e quadrati ha ordine 16; il sottogruppo proprio \(H^+=H\cap SO(3)\) ha ordine 8 e ammette il sollevamento spinoriale \(SU(2)\). Gli altri 8 elementi di \(H\) agiscono sul complesso come simmetrie improprie e richiedono un'estensione Pin se si vuole un'azione spinoriale. L'intero \(O_h\) non agisce come automorfismo del Q6 fissato.

La connessione \(SU(2)^{192}\) è definita come campo di link sul complesso; le sue 192 matrici non sono dedotte univocamente dalle sole coordinate RGS. La verifica numerica usa una connessione generica riproducibile e certifica separatamente:
- appartenenza dei link e delle olonomie a \(SU(2)\);
- covarianza locale delle olonomie;
- invarianza di \(S_W\) sotto \(SU(2)^{64}\);
- invarianza di \(S_W\) sotto le 16 simmetrie del complesso;
- azione spinoriale sollevata per gli 8 elementi propri.

Restano aperti la ricostruzione fisica univoca della connessione a partire dal vierbein/coframe RGS, il sollevamento Pin degli elementi impropri e il limite continuo Yang–Mills. La fase costituisce quindi una costruzione e un audit di gauge theory su complesso finito, non ancora una derivazione della dinamica fisica continua.

**Regola editoriale:** la Fase 19 è un avanzamento interno della medesima catena RGS–PSU–PRISM–QLG e va collocata nel capitolo pertinente dell'indice canonico a 27 capitoli, senza alterarne titoli o numerazione.


## Fase 20 — QAT, doppio foglio PSU, entanglement e geometria emergente

La Fase 20 integra il quadro QAT (cubo, sfera inscritta, cono nullo, evento centrale e doppio foglio PSU–Anti-Twin) con il programma di emergenza geometrica da entanglement. Prosegue la catena già sviluppata nelle Fasi 18–19: geometria RGS e simmetrie \(O_h\), complesso Q6, connessione \(SU(2)^{192}\), olonomie e azione di Wilson.

Il primo ponte quantitativo è un modello finito a due fogli con spazio di Hilbert \(\mathcal H_+\otimes\mathcal H_-\), ciascun fattore di dimensione 64. Sul grafo Q6 si definisce \(K_\tau=e^{-\tau L_{Q6}/2}\) e lo stato normalizzato con coefficienti \(C_\tau=K_\tau/\|K_\tau\|_F\). La matrice ridotta è
\[
\rho_+(\tau)=\frac{e^{-\tau L_{Q6}}}{\operatorname{Tr}(e^{-\tau L_{Q6}})}.
\]
Dallo spettro \(\lambda_j=2j\), \(m_j=\binom6j\), segue esattamente
\[
S_+(\tau)=6h\!\left(\frac{1}{1+e^{2\tau}}\right).
\]
Il risultato è un'identità del modello heat-kernel finito; \(\tau\) è adimensionale e non è ancora identificato con tempo proprio, velocità QAT, raggio del cono o parametro di ovalizzazione.

Van Raamsdonk fornisce la motivazione fisica: in contesti gravitazionali olografici la connessione geometrica può essere legata all'entanglement. Ryu–Takayanagi fornisce il programma area-entropia; Maldacena–Susskind il quadro congetturale ER=EPR; Swingle la direzione di scala tramite entanglement renormalization; Faulkner–Lewkowycz–Maldacena ricordano che esistono correzioni quantistiche oltre il termine di area. Van Raamsdonk–Vilar López (2026) offre un confronto specifico per saddles euclidei, cosmologie chiuse e tubi AdS. Nessuna di queste corrispondenze è assunta automaticamente valida per QAT o PSU.

La fase apre i test: (i) costruire uno scalare geometrico \(\chi[e]\) dalla tetrade QAT e derivare, non adattare, \(\tau=f(\chi)\); (ii) definire superfici minime/cut discrete Q6–PRISM con pesi derivati dal coframe; (iii) verificare \(S(A)=\mathrm{Area}(\gamma_A)/(4G_{\rm eff})\) su più partizioni e raffinamenti, con \(G_{\rm eff}\) fissato indipendentemente; (iv) accoppiare il kernel ai link SU(2) solo tramite un termine gauge-invariante esplicito; (v) confrontare saddles euclidei e tetrade lorentziana dichiarando continuazione e condizioni al bordo.

**Classificazione:** formula entropica Q6 = IDENTITÀ FINITA ESATTA per lo stato scelto; audit matriciale = VERIFICA NUMERICA; kernel e fattorizzazione a due fogli = DEFINIZIONE DI MODELLO; derivazione QAT, legge area, metrica emergente, limite continuo e identificazione con wormhole = APERTI.

**Regola editoriale:** la Fase 20 resta un inserto nella monografia unitaria di 27 capitoli; non crea capitoli autonomi né modifica l'indice canonico.

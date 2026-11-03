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

# RGS–PSU wormhole: recupero dei parametri e audit della derivazione

**Data audit:** 2026-11-09  
**Scopo:** recuperare la provenienza di \(r_*=0.990707\), \(R_0=2.717814\), \(r_1=1.089778\), \(r_2=1.981414\), distinguere \(R_1\) da \(r_1\), e verificare quali relazioni geometriche sono effettivamente dimostrabili.

## 1. Esito del recupero documentale

La ricerca nei repository accessibili del programma RGS–PSU (in particolare `rgs-qlg-quaternion-bridge`, `theorem`, `kilaw`, `biverse`, `cosmodex` e `naressi-zuliani-law-quantum-morphodynamics`) non ha restituito una sorgente versionata contenente congiuntamente i quattro valori né la loro derivazione originaria. La ricerca per stringhe esatte non è una prova dell'assenza di materiale non indicizzato o conservato fuori da GitHub; significa che **la provenienza originaria non è stata recuperata in questa revisione**.

La cronologia personale disponibile conferma che il 20 ottobre 2026 erano già stati registrati:
- \(r_*=0.990707\) e \(R_0=2.717814\) come parametri distinti;
- \(r_1=1.089778\) e \(r_2=1.981414\) come estremi/parametri del raccordo;
- la distinzione esplicita tra \(r_1\) minuscolo e \(R_1\) maiuscolo.

Non è disponibile, in quella traccia recuperata, l'equazione che seleziona tali numeri, la definizione geometrica completa dei simboli o la derivazione di \(R_1\). Perciò nessuno di questi elementi viene ricostruito per congettura.

## 2. Relazioni numeriche interne ai parametri storici

Alla precisione riportata:
\[
1.1\,r_*=1.1(0.990707)=1.0897777\simeq r_1,
\qquad
2r_*=1.981414=r_2.
\]

Queste sono **compatibilità numeriche** dei valori arrotondati. Non dimostrano che le definizioni originarie fossero \(r_1:=1.1r_*\) e \(r_2:=2r_*\): per promuoverle a definizioni o teoremi serve la formula sorgente.

Il parametro \(R_0=2.717814\) non è identificato da queste due relazioni. In particolare non va sostituito con \(r_2\), né con \(r_*\).

## 3. Benchmark Ellis: identità analitiche esatte

Si consideri separatamente il benchmark ultrastatico
\[
 ds^2=-dt^2+d\ell^2+R(\ell)^2d\Omega^2,
 \qquad R(\ell)=\sqrt{R_0^2+\ell^2},\quad R_0>0.
\]
Allora
\[
R'(\ell)=\frac{\ell}{R(\ell)},\qquad
R''(\ell)=\frac{R_0^2}{R(\ell)^3},
\]
e dunque
\[
R'(0)=0,\qquad R''(0)=\frac1{R_0}.
\]
Inserendo il valore storico \(R_0=2.717814\) **come parametro del benchmark**:
\[
R''(0)=1/2.717814\simeq0.367942766.
\]
Questa uguaglianza controlla il numero riportato, ma non prova che il benchmark Ellis sia la soluzione delle equazioni RGS–PSU.

La funzione di forma Morris–Thorne associata alla stessa geometria si ricava da
\[
1-\frac{b(R)}R=\left(\frac{dR}{d\ell}\right)^2
=1-\frac{R_0^2}{R^2},
\]
per cui
\[
b(R)=\frac{R_0^2}{R},\qquad
b(R_0)=R_0,\qquad b'(R_0)=-1.
\]
Il raccordo \(b'(R_0)=-1\) è quindi esatto **per questo benchmark**.

Con la definizione geometrica della massa di Misner–Sharp
\[
M_{MS}(\ell)=\frac{R}{2G}\left(1-(R')^2\right),
\]
si ottiene
\[
M_{MS}(\ell)=\frac{R_0^2}{2G\sqrt{R_0^2+\ell^2}}.
\]
In unità \(G=1\), a \(\ell=20\), il valore è circa \(0.182981\) (la cifra dipende dall'arrotondamento di \(R_0\)). Inoltre
\[
\lim_{|\ell|\to\infty}M_{MS}(\ell)=0.
\]
Per il benchmark asintoticamente piatto questo dà massa asintotica nulla. Non va confusa con la massa alla gola \(M_{MS}(0)=R_0/(2G)\).

## 4. Cosa dicono realmente i file correnti

### 4.1 `experiments/rgs_psu_core/RGS_PSU_closure.py`

Il file documenta esplicitamente tre background distinti:
1. \(R(0)=2.31\), \(\phi(0)=0.80\);
2. \(R(0)=0.8292\), \(N(0)=1.4529\), \(\phi(0)=0.3450\);
3. \(R(0)=2\), \(N(0)=0.9821\), \(\phi(0)=0.3450\).

La routine dimostrativa usa \(R(\ell)=\sqrt{\ell^2+2.31^2}\) come profilo di base. Nella funzione di evoluzione, \(N'=0\) è marcato come semplificazione e \(\phi''\) come placeholder; il commento del file dichiara aperta la chiusura completa delle equazioni e del BVP accoppiato. Il file corregge inoltre la massa alla gola per il caso \(R(0)=2.31\) in \(1.155\) per \(G=1\).

**Conclusione:** questo script non deriva \(r_*\), \(R_0=2.717814\), \(r_1\), \(r_2\), né \(R_1\); non è una soluzione completa RGS–PSU che possa essere usata come loro provenienza.

### 4.2 `docs/RGS_ELLIS_CURVED_REFINEMENT_G1.md`

Questo documento costruisce una famiglia geometrica Ellis–RGS dal seed RGS 64B e usa un raggio di gola indicato come \(r_0=0.8292\). Dichiara espressamente che la famiglia geometrica e gli audit P1 non costituiscono una soluzione completa delle equazioni di campo RGS–EMS. È un terzo oggetto, distinto dal benchmark con \(R_0=2.717814\).

### 4.3 Catalogo centrale

La sezione H4 del catalogo classifica i risultati BVP/Hessiana/WKB come risultati numerici di un modello ridotto e l'interpretazione del throat come recupero della simmetria Lorentziana PSU–Anti-Twin come ipotesi fisica aperta. Questa classificazione resta in vigore.

## 5. Registro epistemico aggiornato

| Affermazione | Stato | Motivo |
|---|---|---|
| Valori storici \(r_*,R_0,r_1,r_2\) | RECUPERATI COME DATI STORICI | Traccia del 20-10-2026; derivazione originaria non recuperata |
| \(r_1\simeq1.1r_*\), \(r_2=2r_*\) | COMPATIBILITÀ NUMERICA | Conseguenza dei numeri arrotondati, non definizione dimostrata |
| \(R''(0)=1/R_0\) | PROVATO NEL BENCHMARK ELLIS | Derivazione diretta da \(R(\ell)=\sqrt{R_0^2+\ell^2}\) |
| \(b'(R_0)=-1\) | PROVATO NEL BENCHMARK ELLIS | \(b(R)=R_0^2/R\) |
| \(M_{MS}(20)\simeq0.182981\), limite zero | PROVATO NEL BENCHMARK ELLIS, \(G=1\) | Formula analitica sopra |
| I quattro parametri sono soluzione delle equazioni RGS–PSU | NON STABILITO | Manca la legge di raccordo e il BVP completo |
| Significato e valore di \(R_1\) | OPEN / NON RECUPERATO | Non equiparare a \(r_1\) |
| Identificazione del benchmark con la famiglia Ellis–RGS \(r_0=0.8292\) | NON GIUSTIFICATA | Scale e background differenti |

## 6. Equazioni che mancano per chiudere la derivazione originaria

Per recuperare la soluzione, il materiale sorgente deve fissare, senza ambiguità:
1. la coordinata radiale usata per \(r_*\), \(r_1\), \(r_2\) (coordinata propria, areale o variabile adimensionale);
2. la definizione di \(R_0\) e l'eventuale legge \(R(r)\);
3. la funzione di raccordo e le condizioni di matching ai due estremi;
4. il sistema variazionale completo per \(R,N,\phi\) e gli eventuali campi PSU/EMS;
5. la definizione separata di \(R_1\), inclusa la relazione con il raggio areale o con il dato al bordo.

Finché queste formule non sono disponibili, il risultato corretto è un **audit di compatibilità**, non una ricostruzione inventata della soluzione. Il prossimo aggiornamento matematico dovrà partire dal file/manoscritto originario o dal notebook che contiene la legge di raccordo, e aggiungere una derivazione riproducibile con residui e condizioni al bordo.

## 7. Riferimenti interni

- `experiments/rgs_psu_core/RGS_PSU_closure.py`
- `docs/RGS_ELLIS_CURVED_REFINEMENT_G1.md`
- `docs/MASTER_CATALOGUE_LAWS_THEOREMS.md`, sezione H4

**Commit di riferimento consultato per il codice e i documenti:** `1c9e7007a33d193060a67002f296760d18ef8140`. Il presente audit non dichiara di aver eseguito test numerici o di aver recuperato il notebook sorgente originario.


## 8. Recupero della derivazione analitica del 20 ottobre 2026

La ricerca nella cronologia della discussione del 20-10-2026 ha recuperato il nucleo analitico che non era presente nei file GitHub consultati. La formulazione riportata è
\[
\phi(r)=\alpha r^2+\gamma\log(r/\ell),\qquad
R(r)=r e^{\phi(r)},\qquad
A(r)=4\pi r^2e^{2\phi(r)}=4\pi R(r)^2.
\]
Il benchmark registrato è
\[
\alpha=1,\qquad \gamma=-2.963,\qquad \ell=1.
\]

### 8.1 Punto stazionario dell'area

Poiché
\[
\log A=\log(4\pi)+2\log r+2\phi(r),
\]
si ha
\[
\frac{A'}{A}=\frac{2(1+\gamma)}r+4\alpha r.
\]
La condizione \(A'(r_*)=0\) fornisce
\[
r_*^2=-\frac{1+\gamma}{2\alpha}.
\]
Con i parametri del benchmark:
\[
r_* = \sqrt{0.9815}=0.990706818387761
\simeq 0.990707.
\]
Questo è il punto di minimo dell'area di base: infatti
\[
(\log A)''(r_*)=-\frac{2(1+\gamma)}{r_*^2}+4\alpha
=8\alpha>0,
\]
e quindi \(A''(r_*)=8\alpha A(r_*)>0\).

### 8.2 Raggio areale e audit numerico

Il raggio areale non è la coordinata \(r_*\), ma
\[
R_0:=R(r_*)=r_*\exp\!\left(\alpha r_*^2+
\gamma\log(r_*/\ell)\right).
\]
Sostituendo \((\alpha,\gamma,\ell)=(1,-2.963,1)\):
\[
R_0=2.7178138074210314\simeq2.717814.
\]
Di conseguenza
\[
A(r_*)=4\pi R_0^2=92.8216459798357
\simeq92.8216,
\]
\[
A''(r_*)=8A(r_*)=742.5731678386855
\simeq742.58.
\]
Questi numeri ricostruiscono simultaneamente il valore storico di \(r_*\), il raggio areale \(R_0\), il minimo \(A\) e la curvatura locale \(A''\). La relazione è una derivazione analitica esatta per l'ansatz dichiarato, con differenze soltanto da arrotondamento nei valori riportati.

### 8.3 Raccordo quintico e significato di r1, r2

La cronologia del 20-10-2026 chiarisce che \(r_1\) e \(r_2\) sono i bordi della fascia di transizione della funzione cutoff, non il punto del minimo originario dell'area:
\[
\chi(r)=
\begin{cases}
1,&r\le r_1,\\
1-10s^3+15s^4-6s^5,&r_1<r<r_2,\\
0,&r\ge r_2,
\end{cases}
\qquad s=\frac{r-r_1}{r_2-r_1}.
\]
La funzione soddisfa \(\chi(r_1)=1\), \(\chi(r_2)=0\) e
\[
\chi'(r_1)=\chi'(r_2)=
\chi''(r_1)=\chi''(r_2)=0,
\]
perciò l'incollaggio con i tratti costanti è \(C^2\). I parametri storici sono
\[
r_1=1.089778\simeq1.1r_*,
\qquad r_2=1.981414\simeq2r_*.
\]
Le uguaglianze con \(1.1r_*\) e \(2r_*\) valgono alla precisione di arrotondamento dei dati. La funzione cutoff è quindi ricostruita; non si deve confondere il suo supporto di transizione con la posizione del minimo di \(A(r)\).

La cronologia menziona anche un minimo secondario di transizione a \(r=1.881250\), con valore riportato \(A=46.2725\). Questo valore appartiene al profilo/oggetto di transizione e non coincide con la funzione di area di base \(A(r)=4\pi r^2e^{2\phi(r)}\) qui sopra; la definizione completa del profilo modificato che lo genera non è stata recuperata e resta da specificare.

### 8.4 Coordinata propria e cautela geometrica

La discussione storica distingue la coordinata di costruzione \(r\) dalla coordinata propria \(l\), usando la forma
\[
ds^2=-e^{2\Phi(l)}dt^2+dl^2+R(l)^2d\Omega^2,
\qquad R(0)=R_0,\quad R'(0)=0,\quad R''(0)>0.
\]
È inoltre riportata la relazione di coordinata propria Morris–Thorne
\[
dl=\frac{dr}{\sqrt{1-b(r)/r}}.
\]
La trasformazione esatta che collega questo \(l\) alla parametrizzazione conformale \(R(r)=re^{\phi(r)}\), e la funzione \(b(r)\) corrispondente, non sono esplicitate nella traccia recuperata. Non si deve quindi assumere automaticamente che il punto stazionario di \(A(r)\) e la condizione \(b(r_*)=r_*\) siano la stessa condizione geometrica senza derivare tale mappa.

### 8.5 Simbolo maiuscolo R1

La ricerca della discussione non ha restituito una definizione indipendente o un valore per \(R_1\) maiuscolo. Rimane quindi distinto da \(r_1=1.089778\), che è il bordo interno del cutoff. Non si assegna a \(R_1\) il valore \(R(r_1)\) senza una convenzione esplicita che lo definisca.

## 9. Provenienza documentale recuperata

- La derivazione analitica e i valori del 20-10-2026 sono stati recuperati dalla cronologia della conversazione: formule per \(\phi,R,A\), benchmark \((1,-2.963,1)\), minimo \(r_*\), raggio areale \(R_0\), audit \(A,A''\) e cutoff quintico.
- Nella Libreria è stato trovato `FORMULAZIONE.pdf` (70 pagine, caricato il 26-10-2026), che contiene una formulazione teleparallela Percorso A e background distinti \(R(0)\simeq2.3105\) e \(R(0)\simeq0.8292\). Non contiene i quattro parametri del 20 ottobre e non è la fonte della derivazione qui ricostruita.
- Non è stata individuata una copia `FORMULAZIONE.odt` nella Libreria. Il PDF trovato non deve essere usato per sostituire o confondere il benchmark \(R_0=2.717814\).

**Stato aggiornato:** DERIVAZIONE ANALITICA DEL MINIMO DI AREA RECUPERATA E RICOSTRUITA; CUT-OFF QUINTICO RECUPERATO; MAPPATURA ALLA COORDINATA PROPRIA E DEFINIZIONE DI \(R_1\) ANCORA APERTE.

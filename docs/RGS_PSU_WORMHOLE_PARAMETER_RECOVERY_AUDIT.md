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

# Programma analitico per la convergenza Mosco della geometria curva RGS

**Versione 0.1 — 3 novembre 2026**  
**Stato: PROGRAMMA DI DIMOSTRAZIONE / APERTO.** Questo documento non dichiara dimostrata la convergenza Mosco di RGS-Q6. Specifica la famiglia matematica che deve essere costruita e le verifiche necessarie.

## 1. Obiettivo preciso

Fissare una varietà compatta, connessa e senza bordo (M) di dimensione (d), dotata di metrica Riemanniana (g) e misura (dmu_g). L'obiettivo iniziale è costruire complessi simpliciali/geometrici (K_h), con passo massimo (h\to0), e forme quadratiche discrete (q_h) che convergano in senso Mosco alla forma di Dirichlet
[
q(u)=\int_M |\nabla_g u|_g^2\,d\mu_g,
\qquad
\mathcal D(q)=H^1(M).
]
Questa prima tappa è scalare e geometrica. Non include ancora spinori, connessione SU(2), Yang–Mills o vierbein dinamico.

## 2. Spazio di Hilbert comune

Si deve scegliere e documentare una delle due impostazioni, senza mescolarle implicitamente.

### A. Embedding in uno spazio comune

Sia (H=L^2(M,d\mu_g)). Per ogni (K_h), definire uno spazio discreto (H_h), una mappa di ricostruzione/interpolazione (J_h:H_h\to H) e una misura discreta coerente. Estendere la forma a (H) ponendo
[
\mathcal Q_h(v)=
\begin{cases}
q_h(u_h),&v=J_hu_h,\
+\infty,&v\notin J_h(H_h).
\end{cases}
]
Se (J_h) non è iniettivo, occorre fissare una scelta canonica o formulare la teoria su classi equivalenti.

### B. Spazi variabili

In alternativa, usare una nozione di convergenza per forme su spazi di Hilbert variabili, specificando mappe di identificazione, convergenza forte/debole e compatibilità delle norme. Non basta confrontare gli autovalori di matrici di dimensioni diverse.

Per la prima implementazione si raccomanda A, perché rende espliciti liminf e recovery nello stesso (L^2(M)).

## 3. Geometria discreta: dati obbligatori

Per ogni livello (h), il manifest della mesh deve contenere:

1. vertici (x_i^{(h)}\in M) e simplici orientati;
2. diametro massimo (h), raggio interno e costante di shape-regularity;
3. lunghezze geodetiche degli spigoli e misure dei simplessi;
4. volumi duali, ove usati, con regola di costruzione dichiarata;
5. matrici d'incidenza (B_k^{(h)}), orientazioni e verifica (B_kB_{k+1}=0);
6. mappa di ricostruzione (J_h), inner product discreto e misura;
7. metrica/coframe campionati e stima uniforme degli errori geometrici.

**Vincolo RGS:** le lunghezze, i volumi e i pesi devono essere derivati dalla metrica/coframe dichiarati. Non introdurre coefficienti cubo-sfera o pesi di adattamento liberi per migliorare lo spettro. Ogni correzione deve essere derivata da un funzionale/azione esplicito e registrata come nuova definizione del modello.

## 4. Forma discreta candidata

Per una discretizzazione scalare su grafo 1-scheletrico, una forma candidata è
[
q_h(u_h)=
\sum_{(i,j)\in E_h}
c_{ij}^{(h)}|u_j-u_i|^2,
\qquad
c_{ij}^{(h)}\ge0.
]
I coefficienti non sono parametri arbitrari: una scelta geometrica di riferimento è
[
c_{ij}^{(h)}
=
\frac{|*e_{ij}^{(h)}|}{|e_{ij}^{(h)}|},
]
quando la cella duale (*e) e la sua misura sono definite in modo coerente e positivo. In una discretizzazione simpliciale FEM, la forma deve invece essere derivata dall'integrale esatto o controllato di (\langle\nabla u_h,\nabla v_h\rangle_g). Le due costruzioni non vanno identificate senza dimostrazione.

Per Q6, il grafo (Q_6) da solo è un complesso finito di riferimento, non una famiglia (K_h) con (h\to0). Per ottenere una successione occorre definire esplicitamente come i gradi di libertà e i complessi si raffinano e come ogni livello approssima la stessa (M).

## 5. Enunciato bersaglio di Mosco

Per le forme estese (\mathcal Q_h) su (H=L^2(M)), dimostrare entrambe le condizioni.

**(M1) Liminf debole.** Per ogni successione (v_h\rightharpoonup v) debolmente in (H),
[
\mathcal Q(v)\le\liminf_{h\to0}\mathcal Q_h(v_h).
]

**(M2) Recovery forte.** Per ogni (v\in H), esiste (v_h\to v) fortemente in (H) tale che
[
\mathcal Q(v)\ge\limsup_{h\to0}\mathcal Q_h(v_h).
]

Qui (\mathcal Q=q) se (v\in H^1(M)), e (+\infty) altrimenti. Per (v\notin H^1), la recovery condition è formalmente triviale solo dopo aver fissato correttamente gli spazi e la nozione di convergenza; la parte sostanziale è la densità per (H^1).

## 6. Lemmi tecnici che devono essere provati

**L1 — Consistenza geometrica.** Le misure discrete e i tensori metrici ricostruiti convergono a (d\mu_g) e (g) con controllo uniforme sufficiente. Una convergenza soltanto puntuale sui vertici non basta.

**L2 — Coercività uniforme.** Esiste (c>0), indipendente da (h), tale che la forma discreta controlla una norma discreta (H^1), a meno delle costanti:
[
q_h(u_h)+\|u_h\|_{H_h}^2
\ge c\,\|u_h\|_{H_h^1}^2.
]
La definizione di (H_h^1) deve essere esplicita.

**L3 — Compattezza.** Successioni con norme (L^2) e energie uniformemente limitate ammettono sottosuccessioni ricostruite convergenti fortemente in (L^2(M)). Questo è il passaggio che esclude oscillazioni microscopiche non controllate.

**L4 — Liminf.** L'energia limite di ogni limite debole è non superiore al liminf discreto. Occorre identificare il gradiente debole limite e controllare i coefficienti metrici.

**L5 — Recovery.** Per ogni funzione liscia (u), costruire (u_h) (campionamento o interpolante) e dimostrare (J_hu_h\to u) e (q_h(u_h)\to q(u)). Estendere poi a (H^1(M)) per densità e diagonalizzazione.

**L6 — Normalizzazione.** Verificare che il limite sia esattamente (q(u)=\int|\nabla u|^2d\mu_g), non (c_*q(u)) con un fattore nascosto. Se compare (c_*), deve essere derivato dalla discretizzazione, non corretto a posteriori.

## 7. Sequenza di implementazione e gate PASS/FAIL

| Gate | Artefatto richiesto | Criterio PASS |
|---|---|---|
| G0 | Specifica (M,g), unità e coordinate | Dati geometrici e convenzioni riproducibili |
| G1 | Generatore (K_h), almeno 3 livelli iniziali | (h\downarrow0), complessi validi, shape-regularity misurata |
| G2 | Manifest geometrico per ogni livello | Tutti i pesi derivano da lunghezze/volumi/coframe |
| G3 | (B_k), Hodge star, inner product | (B_kB_{k+1}=0), positività e ranghi verificati |
| G4 | (J_h), (H_h), (q_h) | Mappe e norme definite senza ambiguità |
| G5 | Audit numerico | Coercività, consistenza, errori e scaling su raffinamenti |
| G6 | Prova L1–L6 | Liminf e recovery dimostrati analiticamente |
| G7 | Resolventi/spettro | Conseguenze spettrali verificate dopo Mosco + compattezza |

Un fit di potenza su tre o più livelli è diagnostica, non prova di Mosco. Un residuo piccolo su una singola matrice non certifica coercività uniforme.

## 8. Estensione successiva: spinori e gauge

Solo dopo G6 si passa a:

1. bundle spinoriale e struttura Spin/Pin appropriata su (M);
2. connessione spin (\nabla^S) definita dalla metrica e dalla connessione scelta;
3. trasporti discreti (S_{ij}^{(h)}\in SU(2)) ottenuti da una discretizzazione consistente della connessione, non dal solo polar factor senza controllo;
4. forma/operatore covariante con link e prova di gauge covariance;
5. stime uniformi sul gap del kernel Wilson e località overlap;
6. convergenza dei resolventi/forme dell'operatore di Dirac, più controllo separato del doubling e dell'indice.

Il lift polare RGS già verificato numericamente è un punto di partenza sperimentale, non ancora la dimostrazione che (S_{ij}^{(h)}) approssima il trasporto della connessione spin.

## 9. Cosa il prototipo del toro piatto già consente di dire

Il documento `MOSCO_PROTOTYPE_REGULAR_TORUS.md` dimostra un caso modello periodico piatto sotto le proprie ipotesi Fourier e di embedding comune. Quel risultato è un precedente metodologico utile per testare liminf/recovery e scaling. Non fornisce automaticamente:
- una famiglia curva RGS;
- controllo uniforme di curvatura, coframe o dual cells;
- compattezza per la mesh RGS;
- convergenza di link SU(2);
- convergenza di Dirac/overlap;
- un limite Yang–Mills o gravitazionale.

## 10. Deliverable successivo

Il primo esperimento concreto deve essere un benchmark geometrico controllato su una varietà nota (ad esempio un toro piatto) per validare generatore, misure, (J_h), energia e audit L1–L5. Solo dopo si sostituisce la geometria di test con una famiglia derivata dal coframe RGS, mantenendo invariata l'interfaccia matematica.

**Decisione scientifica:** fino al superamento di G6, lo stato del gate è **OPEN**. La convergenza Mosco RGS curva non è ancora dimostrata.

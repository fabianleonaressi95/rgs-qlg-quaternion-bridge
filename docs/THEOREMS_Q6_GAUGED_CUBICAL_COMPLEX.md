# Theoremi strutturali per il complesso cubicale Q6 gaugiato

## 1. Dati e convenzioni

Sia (Q_6=[0,1]^6) il cubo combinatorio con vertici (v\in\{0,1\}^6), celle orientate dai loro insiemi ordinati di direzioni e fibra interna (V=\mathbb C^2). Per ogni arco orientato positivamente (v\to v+e_i) assegniamo (U_i(v)\in SU(2)), che trasporta la fibra dell'estremo superiore alla fibra dell'estremo inferiore. L'arco inverso usa (U_i(v)^\dagger).

Ogni cella (c=(v; i_1<\cdots<i_k)) porta la fibra ancorata nel vertice inferiore (v). Il differenziale covariante (D_{U,k}:C^k(Q_6;V)\to C^{k+1}(Q_6;V)) e definito, per una cella (c'=(v; i_0<\cdots<i_k)), da

\[
(D_{U,k}\omega)(c')=\sum_{p=0}^{k}(-1)^p
\left[U_{i_p}(v)\,\omega(v+e_{i_p};I\setminus i_p)
-\omega(v;I\setminus i_p)\right].
\]

La notazione della faccia superiore indica la cella parallela con ancora (v+e_{i_p}). Il differenziale globale (D_U) e la somma dei blocchi (D_{U,k}), con dimensione complessiva (2\cdot729=1458).

Per un quadrato (s=(v;i,j)), (i<j), poniamo

\[
A_{ij}(v)=U_i(v)U_j(v+e_i),\qquad
B_{ij}(v)=U_j(v)U_i(v+e_j),
\]
\[
H_{ij}(v)=A_{ij}(v)B_{ij}(v)^\dagger\in SU(2).
\]

(A) e (B) sono i due trasporti dal vertice superiore (v+e_i+e_j) all'ancora (v); (H) e l'olonomia orientata del quadrato, espressa nella fibra in (v).

## 2. Teorema di covarianza di gauge

**Teorema 1 (covarianza esatta).** Sia (g:Q_6^{(0)}\to SU(2)) una trasformazione di gauge ai vertici e si trasformino gli archi secondo

\[
U_i^g(v)=g(v)U_i(v)g(v+e_i)^\dagger.
\]

Sulla fibra di ogni cella (c) con ancora (v), definiamo (P_k(g)\omega(c)=g(v)\omega(c)). Allora

\[
D_{U^g,k}=P_{k+1}(g)D_{U,k}P_k(g)^\dagger.
\]

**Dimostrazione.** In ciascun addendo superiore,

\[
U_i^g(v)g(v+e_i)\omega_{sup}
=g(v)U_i(v)\omega_{sup},
\]

mentre il termine inferiore trasformato e (g(v)\omega_{inf}). Il fattore (g(v)) e comune ai due termini e quindi a tutta la somma alternata. Questo prova l'identita blocco per blocco. La dimostrazione e algebrica e non dipende da una scelta numerica dei link. (square)

Conseguenza: l'operatore di Dirac–Kähler hermitiano (K_U=D_U+D_U^\dagger) soddisfa (K_{U^g}=P(g)K_UP(g)^\dagger), poiche (P(g)) e unitario.

## 3. Teorema locale di curvatura

**Teorema 2 (curvatura come difetto di nilpotenza).** Sul quadrato (s=(v;i,j)), per una 0-cochain (f),

\[
(D_{U,1}D_{U,0}f)(s)
=\big[A_{ij}(v)-B_{ij}(v)\big]f(v+e_i+e_j)
=\big[H_{ij}(v)-I\big]B_{ij}(v)f(v+e_i+e_j).
\]

Quindi (D_{U,1}D_{U,0}) si annulla su tutti i quadrati se e solo se (H_{ij}(v)=I_2) per ogni quadrato.

**Dimostrazione.** Si espande la differenza covariante sui quattro lati. I termini associati ai due vertici intermedi e all'ancora si cancellano a coppie per l'alternanza dei segni; rimangono soltanto i due trasporti dal vertice opposto all'ancora, con segno opposto. Poiche (B_{ij}(v)) e invertibile, (A-B=0) equivale a (AB^\dagger=I). (square)

## 4. Teorema globale di piattezza sul cubo

**Teorema 3 (criterio di nilpotenza).** Per il cubo (Q_6), con la definizione di (D_U) sopra, sono equivalenti:

1. (D_U^2=0) su tutte le cochain di tutti i gradi;
2. (D_{U,1}D_{U,0}=0);
3. tutte le 240 olonomie elementari (H_{ij}(v)) sono identita.

**Dimostrazione.** (1\Rightarrow2) e immediato per restrizione ai primi due gradi. (2\Rightarrow3) segue dal Teorema 2, scegliendo una 0-cochain supportata nel singolo vertice superiore del quadrato: il blocco corrispondente di (D_1D_0) e (A-B), percio deve annullarsi e (H=I). Per (3\Rightarrow1), consideriamo (D_{U,k+1}D_{U,k}) su una cella di grado (k+2). Ogni contributo e indicizzato da una coppia di direzioni eliminate. I due ordini di eliminazione della stessa coppia hanno segni opposti; i rispettivi trasporti differiscono soltanto per il trasporto lungo i due lati del quadrato elementare generato da quella coppia. Se la sua olonomia e identita, i due contributi coincidono e si cancellano. Tutti i contributi si cancellano a coppie, per ogni grado. (square)

Il risultato e piu forte della sola osservazione che la connessione identicamente unitaria produce il complesso ordinario: caratterizza tutte le connessioni discrete piatte, incluse quelle pure-gauge non costanti.

## 5. Teorema di hermiticita e covarianza spettrale

**Teorema 4.** Su somma ortogonale delle cochain, (K_U=D_U+D_U^\dagger) e hermitiano. Sotto gauge, (K_{U^g}=P(g)K_UP(g)^\dagger). Di conseguenza lo spettro di (K_U), con molteplicita, e gauge-invariante.

**Dimostrazione.** L'hermiticita segue da (K_U^\dagger=D_U^\dagger+D_U). La covarianza segue dal Teorema 1 e dall'unitarieta di (P(g)); l'invarianza spettrale segue dall'equivalenza unitaria. (square)

## 6. Separazione del modulo di Clifford

Sia (S=\Lambda^*\mathbb C^6) il modulo esterno di dimensione 64 e siano (\gamma_i=\varepsilon_i+\iota_i) i generatori euclidei con (\{\gamma_i,\gamma_j\}=2\delta_{ij}I_S). Sul prodotto (V\otimes S), l'azione interna (U\otimes I_S) commuta con (I_V\otimes\gamma_i). E un teorema di algebra tensoriale, non una identificazione tra il gruppo di gauge interno (SU(2)) e (Spin(6)).

Un vero operatore di spin richiede un fibrato spinoriale, un sollevamento della connessione di frame a (Spin(6)) (o alla firma pertinente), e la compatibilita con la coframe discreta. Nessuno di questi dati e dedotto dalla sola commutazione dei fattori.

## 7. Stato epistemico e programma

- **Dimostrato algebricamente nel modello finito:** covarianza di gauge, criterio di piattezza, hermiticita e invarianza spettrale.
- **Verificabile numericamente:** matrici sparse, tutte le 240 olonomie, residui di covarianza e nilpotenza per connessioni piatte/pure-gauge e curve.
- **Non dimostrato:** limite continuo, convergenza di Mosco, identificazione con Dirac spinoriale, assenza di doubling, derivazione dinamica della connessione da vierbein RGS.
- **Passo successivo:** costruire il lift (I\to2I\) e il fibrato spinoriale discreto in modo indipendente, poi dimostrare compatibilita con la coframe e studiare convergenza spettrale.

Il Teorema 3 e il risultato strutturale centrale: il quadrato del differenziale non e un errore da eliminare in presenza di campo, ma misura esattamente la curvatura discreta.

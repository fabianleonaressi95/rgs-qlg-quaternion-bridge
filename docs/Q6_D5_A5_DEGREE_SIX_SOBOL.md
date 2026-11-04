# Q6–D5: controllo indipendente quasi-Monte Carlo del coefficiente A5 di grado 6

Questo audit aggiunge una seconda stima della stessa quantità definita in
`Q6_D5_A5_DEGREE_SIX_COMPARISON.md`. I dieci gaussiani standard che
parametrizzano `(q,p) in R^5 x R^5` sono ottenuti trasformando punti Sobol
scramblati in `(0,1)^10` tramite la normale inversa. Ogni replica usa uno
scramble indipendente e una rete di `2^m` punti.

Il coefficiente usa esattamente lo stesso generatore armonico A5 e la stessa
normalizzazione L2 su S2 del controllo precedente. Il risultato viene
confrontato sia con il design Q6 finito a 32 direzioni, sia con una stima
gaussiana indipendente avente lo stesso numero di punti per replica.

L'errore standard Sobol è stimato dalla dispersione fra scramble indipendenti.
È un indicatore empirico per randomized QMC, non un limite deterministico
dell'errore d'integrazione. L'accordo fra i due metodi è un controllo di
robustezza numerica, non una dimostrazione analitica dell'integrale o
dell'uguaglianza fra misure.

Restano esclusi: identificazione con la misura fisica di Liouville,
identificazione con la misura Bergman globale D5, limite 3+1D e derivazione
della costante di struttura fine.

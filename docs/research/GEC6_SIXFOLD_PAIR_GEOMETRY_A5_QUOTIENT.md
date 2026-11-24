# GEC-6 — Sixfold Pair Geometry and the A5 root-lattice quotient

**Nome:** Geometria Esagonale delle Coppie (GEC-6)  
**English:** Sixfold Pair Geometry  
**Stato:** costruzione algebrica esatta; la sorgente è il grafo Q6.

## 1. Definizione

Sia \(Q_6^{(1)}\) il grafo ipercubico a sei dimensioni. Il suo reticolo dei cicli ha rango
\[
|E|-|V|+1=192-64+1=129.
\]
Le 15 coppie non ordinate delle sei direzioni costituiscono la base di
\[
\Lambda^2(\mathbb Z^6),\qquad \operatorname{rank}=15.
\]
GEC-6 denomina questa geometria delle coppie e il suo quoziente naturale di rango cinque.

Il reticolo di arrivo è il reticolo delle radici di tipo \(A_5\):
\[
A_5=\{(y_0,\ldots,y_5)\in\mathbb Z^6:\sum_{i=0}^5y_i=0\},
\qquad \operatorname{rank}A_5=5.
\]

## 2. Contrazione verso il reticolo A5

Definiamo
\[
c:\Lambda^2(\mathbb Z^6)\longrightarrow A_5,
\qquad
c(e_i\wedge e_j)=e_j-e_i\quad(i<j).
\]
È un omomorfismo integrale suriettivo: le differenze \(e_j-e_i\) generano il reticolo a somma nulla. È inoltre \(S_6\)-equivariante, perché le permutazioni delle sei coordinate agiscono coerentemente sulle coppie orientate e sulle differenze.

Pertanto
\[
\Lambda^2(\mathbb Z^6)\xrightarrow{c}A_5
\]
ha nucleo di rango \(15-5=10\).

## 3. Mappa composta dal reticolo dei cicli di Q6

Sia
\[
\Phi:Z_1(Q_6^{(1)};\mathbb Z)\longrightarrow\Lambda^2(\mathbb Z^6)
\]
la suriezione integrale \(S_6\)-equivariante della [costruzione precedente](Q6_PRISM_EXTERIOR_SQUARE_EQUIVARIANT_MAP.md): ogni bordo quadrato orientato nelle direzioni \(i<j\) viene mandato in \(e_i\wedge e_j\).

Definiamo
\[
\Psi=c\circ\Phi:
Z_1(Q_6^{(1)};\mathbb Z)\longrightarrow A_5.
\]
Allora \(\Psi\) è integrale, suriettiva e \(S_6\)-equivariante. Su un bordo quadrato orientato:
\[
\boxed{\Psi(\partial Q_{ij})=e_j-e_i.}
\]
Le radici \(e_j-e_i\) generano \(A_5\), quindi la suriettività segue. Per rango-nullità:
\[
\boxed{\operatorname{rank}\ker\Psi=129-5=124.}
\]
Questo è un quoziente naturale di rango cinque del reticolo dei cicli del grafo Q6, non una scelta arbitraria di cinque cicli fondamentali.

Una scelta esplicita di gauge per la 1-cochain contratta è
\[
\beta_d(x)=\sum_{i<d}x_i(e_d-e_i),
\]
orientando ogni arco nella direzione \(d\) dal bit 0 al bit 1. Il suo rotore discreto sul quadrato orientato \(i<j\) è \(e_j-e_i\). La cochain dipende dall'ordinamento delle coordinate; la mappa indotta sui cicli è canonica e \(S_6\)-equivariante. Non si deve affermare che questa particolare cochain sia essa stessa \(S_6\)-equivariante.

## 4. Relazione con le cinque periodicità PRISM

Per il prodotto PRISM originale:
\[
P=C_{64}\square C_6^3\square C_5,
\qquad H_1(P;\mathbb Z)\cong\mathbb Z^5.
\]
Come gruppi abeliani liberi, \(A_5\cong\mathbb Z^5\). Una base, ad esempio le radici semplici
\[
\alpha_i=e_i-e_{i+1},\quad i=0,\ldots,4,
\]
fornisce un isomorfismo di gruppi abeliani tra \(A_5\) e il reticolo delle cinque periodicità PRISM.

**Questo isomorfismo di reticoli non identifica le rappresentazioni di simmetria.** L'azione naturale di \(S_6\) su \(A_5\) è non banale. Il grafo cartesiano PRISM esistente non acquisisce per questo un'azione di \(S_6\) per automorfismi di grafo; il risultato non supera l'ostruzione precedente per una mappa equivarante verso la rappresentazione banale dell'omologia PRISM originale.

Conclusione precisa: GEC-6 fornisce un quoziente naturale di rango cinque, \(S_6\)-equivariante, astrattamente isomorfo a \(\mathbb Z^5\). L'identificazione con le cinque periodicità PRISM è, allo stato attuale, soltanto un isomorfismo astratto di reticoli, non un'equivalenza geometrica di grafi o dinamica.

## 5. Limite del complesso cubico riempito

La sorgente è il **grafo** Q6. Se si riempiono tutte le 240 facce quadrate, il cubo 6-dimensionale è contrattibile e la sua prima omologia è nulla. \(\Psi\) è quindi un invariante del reticolo dei cicli del grafo, non una mappa non nulla da \(H_1\) del cubo riempito.

## 6. Registro delle evidenze

| Affermazione | Stato |
|---|---|
| \(c:\Lambda^2\mathbb Z^6\to A_5\) è integrale e suriettiva | PROVATO |
| \(c\) è \(S_6\)-equivariante | PROVATO |
| \(\Psi=c\circ\Phi\) è una suriezione integrale \(S_6\)-equivariante | PROVATO, usando la costruzione esterna precedente |
| \(\operatorname{rank}\ker\Psi=124\) | PROVATO |
| \(A_5\cong\mathbb Z^5\) come gruppi abeliani liberi | PROVATO |
| Il grafo PRISM originale realizza l'azione \(S_6\) su \(A_5\) per automorfismi di grafo | NON DIMOSTRATO |
| Identificazione geometrica, dinamica o fisica di GEC-6 con PRISM | APERTA |

## 7. Riproducibilità

- Implementazione: \`src/q6_gec6_a5_quotient.py\`
- Test: \`tests/test_q6_gec6_a5_quotient.py\`
- Audit: \`examples/q6_gec6_a5_quotient_audit.py\`
- Risultato JSON: \`results/q6_gec6_a5_quotient_audit.json\`

Esecuzione:
\`\`\`bash
PYTHONPATH=. python -m pytest tests/test_q6_gec6_a5_quotient.py -q
PYTHONPATH=. python examples/q6_gec6_a5_quotient_audit.py
\`\`\`

# Intertwiner fibrico finito RGS–PSU–Q6

## 1. Fattorizzazione esatta delle 64 etichette

La ricostruzione archiviata del seme RGS ordina i 64 punti in otto blocchi
di otto segni. Ogni etichetta binaria si fattorizza quindi come

    g = 8 b + s,   b,s in {0,...,7}.

I tre bit più significativi selezionano uno degli otto rappresentanti
d'orbita p_b; i tre bit meno significativi selezionano i segni cartesiani.
Definendo s_vec in {0,1}^3,

    Phi(b,s_vec) = (2 s_vec - 1) componentwise p_b.

Questa formula ricostruisce esattamente, con residuo massimo zero in
aritmetica float64, tutte le 64 coordinate RGS archiviate. Non è un fit:
segue dall'ordine esplicito usato per generare il seed.

## 2. Collegamento con la proiezione aurea 3+3

Il frame ortogonale F impila P_parallel e P_internal. Per ogni vertice Q6,

    (x_parallel, x_internal) = F (bits - 1/2),
    bits - 1/2 = F^T (x_parallel, x_internal).

Quindi la coppia fisico/interno conserva esattamente i sei bit e permette
di recuperare (b,s). La mappa RGS è non lineare rispetto alle sole
coordinate fisiche tridimensionali: la fibra interna è necessaria per
recuperare l'etichetta completa. Questo non implica che le coordinate
interne siano tempo fisico.

## 3. Metrica discreta indotta sugli archi

Per ciascun arco Q6 e=(v,w), si usa la lunghezza euclidea nella sua immagine
RGS,

    ell_e = ||Phi(v)-Phi(w)||.

Gli archi si separano in 96 archi di base (cambia b, segno fissato) e
96 archi di fibra (cambia un segno, base fissata). Le proiezioni auree
fisiche hanno invece lunghezza uniforme 1/sqrt(2); le lunghezze RGS sono
non uniformi e dipendono dal tipo di arco e dal rappresentante.

Il funzionale P1 unidimensionale sugli archi è

    S(u) = 1/2 sum_e (u_w-u_v)^2 / ell_e.

La matrice di massa consistente assembla gli elementi locali
ell_e/6 [[2,1],[1,2]]. Si ottengono M definita positiva e K semidefinita
positiva; il grafo combinato è connesso e ha un solo modo costante nullo.
I sottografi base e fibra, presi separatamente, hanno ciascuno otto
componenti connesse e rango 56.

## 4. Cosa dimostra e cosa no

**Verificato esattamente:** fattorizzazione delle etichette, mappa discreta
base/fibra, ricostruzione delle 64 coordinate RGS, partizione 96+96 degli
archi.

**Verificato numericamente:** lunghezze indotte, matrici FEM 1D, positività,
rango e spettro generalizzato.

**Non dimostrato:** un fibrato differenziabile continuo, una connessione
di gauge, una metrica spazio-temporale 4D, TEGR, un'identificazione del
tempo con la fibra, o il limite Mosco. Il risultato è un ponte combinatorio
e metrico finito, non ancora una dinamica fisica completa.

## 5. Esecuzione

    python -m src.rgs_psu_fibred_intertwiner --output results/rgs_psu_fibred_intertwiner
    pytest -q tests/test_rgs_psu_fibred_intertwiner.py

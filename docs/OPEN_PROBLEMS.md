# OPEN PROBLEMS
## Domande Aperte del Programma RGS–PSU–PRISM–QLG

**Data:** 3 Novembre 2026  
**Versione:** 1.0  
**Standard:** Classificazione Epistemica Rigorosa

---

## LEGENDA

| Categoria | Definizione |
|-----------|-------------|
| **MATHEMATICAL** | Problema puramente matematico, soluzione analitica possibile |
| **COMPUTATIONAL** | Richiede analisi numerica / implementazione |
| **PHYSICAL** | Richiede interpretazione fisica / verifica sperimentale |
| **FOUNDATIONAL** | Impatta sulla struttura logica dell'intero programma |

---

## PARTE I: PROBLEMI MATEMATICI

### 1.1 Consistenza Globale dei Lift SU(2)

**Categoria:** MATHEMATICAL  
**Priorità:** ALTA  
**Status:** OPEN

**Problema:**
Per un grafo con archi orientati, scegliere un lift SU(2) per ogni arco è un'assegnazione locale. Questa non garantisce automaticamente consistenza globale intorno a ogni ciclo chiuso.

**Domande Specifiche:**
1. Quali assegnazioni di segno centrale sono compatibili con il grafo dato?
2. Come classificare le classi di lift globalmente ammissibili?
3. Quando un cambiamento di lift locale può essere esteso a una trasformazione di gauge?
4. Qual è la relazione tra coomologia ℤ₂ e classi di lift?

**Approccio Sugerito:**
- Formulazione coomologica usando il sottogruppo centrale ℤ₂ = {+1, -1} ⊂ SU(2)
- Studio delle classi di coomologia sul complesso cellulare Q6
- Separazione tra gauge equivalence e lift ambiguity

**Riferimenti:**
- Nakahara, "Geometry, Topology and Physics" (cap. 10)
- Nash & Sen, "Topology and Geometry for Physicists"

---

### 1.2 Distanza da 2I e Convergenza Spontanea

**Categoria:** MATHEMATICAL + COMPUTATIONAL  
**Priorità:** ALTISSIMA  
**Status:** OPEN

**Problema:**
Definire la distanza
\[
d_{2I}(q) = \min_{g \in 2I} \min(\|q - g\|, \|q + g\|)
\]
e determinare se la connessione geometrica Mᵢⱼ = EⱼEᵢ⁻¹ produce lift qᵢⱼ vicini a 2I.

**Domande Specifiche:**
1. Qual è la distribuzione di d_{2I}(qᵢⱼ) per i 192 archi Q6?
2. Esiste un bias statistico verso 2I?
3. Le olonomie U_f di plaquette sono più vicine a 2I dei link individuali?
4. La distanza diminuisce sotto trasformazioni di gauge preferite?

**Approccio Sugerito:**
- Implementare funzione d_{2I}(q) in `src/binary_icosahedral.py`
- Calcolare distribuzione su tutti i link
- Calcolare distribuzione su tutte le olonomie
- Testare invarianza sotto gauge
- Analisi statistica (media, varianza, skewness, kurtosis)

**Risultato Atteso:**
- Istogramma d_{2I} per link
- Istogramma d_{2I} per olonomie
- Test di significatività statistica
- Conclusione: convergenza o no

---

### 1.3 Coomologia del Complesso Q6

**Categoria:** MATHEMATICAL  
**Priorità:** MEDIA  
**Status:** OPEN

**Problema:**
Determinare i gruppi di coomologia del complesso cellulare Q6 con la struttura di incidenza B₁, B₂.

**Domande Specifiche:**
1. Quali sono H⁰, H¹, H² per il complesso Q6?
2. Come si relazionano con la bipartizione del grafo?
3. Esistono cicli non banali oltre alle plaquette?
4. Qual è la struttura dello spazio di cicli 1-dimensionali?

**Approccio Sugerito:**
- Calcolare spazi nucleo di B₁, B₂, B₁ᵀ, B₂ᵀ
- Determinare gruppi di omologia H_k = ker(B_k) / im(B_{k+1})
- Interpretazione geometrica dei generatori

---

### 1.4 Spettro del Laplaciano Q6

**Categoria:** COMPUTATIONAL  
**Priorità:** MEDIA  
**Status:** OPEN

**Problema:**
Calcolare e analizzare lo spettro del Laplaciano del grafo Q6.

**Domande Specifiche:**
1. Qual è il gap spettrale?
2. Esistono autovalori degeneri con significato geometrico?
3. Come si relaziona lo spettro con la struttura ipercubica?
4. Esistono autovalori associati alla simmetria S₆?

**Approccio Sugerito:**
- Diagonalizzazione numerica del Laplaciano (64×64)
- Analisi degli autovalori (gap, degenerazioni)
- Interpretazione degli autovettori (localized vs extended)
- Confronto con grafo ipercubico teorico

---

### 1.5 Struttura Topologica del Grafo PRISM

**Categoria:** MATHEMATICAL + COMPUTATIONAL  
**Priorità:** MEDIA  
**Status:** OPEN

**Problema:**
Verificare i risultati storici sul grafo PRISM C₆₄□C₆³□C₅.

**Domande Specifiche:**
1. N = 69120 è corretto?
2. grado = 10 è corretto?
3. Esistono esattamente 108 cicli di lunghezza 640?
4. b₁ = 5 e b₂ = 10 sono corretti?
5. Gap spettrale = 0.009630546656 è riproducibile?
6. Relazione K³ = K² + I (Krist-All) è verificabile?

**Approccio Sugerito:**
- Implementare costruzione grafo PRISM
- Calcolare numeri di Betti (algebra lineare)
- Calcolare spettro Laplaciano
- Verificare relazione algebrica Krist-All
- Confronto con valori storici

---

### 1.6 Equivalenza TEGR ↔ GR

**Categoria:** MATHEMATICAL  
**Priorità:** MEDIA  
**Status:** OPEN

**Problema:**
Dimostrare l'equivalenza tra Teleparallel Equivalent General Relativity e Relatività Generale per coefficienti generici.

**Domande Specifiche:**
1. Sotto quali ipotesi le equazioni di campo TEGR coincidono con GR?
2. Qual è la relazione esatta tra torsione e curvatura di Riemann?
3. Come si trasformano le equazioni di campo sotto cambi di gauge?
4. Esistono termini di bordo che distinguono le due formulazioni?

**Approccio Sugerito:**
- Derivazione sistematica usando notazione tetradica
- Confronto termine per termine delle azioni
- Analisi delle condizioni di contorno
- Riferimento a letteratura TEGR (Aldrovandi & Pereira)

---

### 1.7 Riduzione SU(2) → U(1)

**Categoria:** MATHEMATICAL  
**Priorità:** MEDIA  
**Status:** OPEN

**Problema:**
Determinare condizioni sotto cui una connessione SU(2) si riduce a U(1).

**Formulazione Candidata:**
Usare un campo aggiunto n_g = n_g^a T_a, \|n_g\| = 1, e definire
\[
\mathcal{R}_{gh} = n_g - U_{gh} n_h U_{gh}^{-1}
\]
La condizione \(\mathcal{R}_{gh} = 0\) su tutti i link definisce un sottogruppo U(1) preferito.

**Domande Specifiche:**
1. Esiste globalmente un campo n_g che soddisfa la condizione?
2. Se esiste, è unico?
3. Quali sono le obstruzioni topologiche?
4. La riduzione è stabile sotto raffinamento del reticolo?

**Approccio Sugerito:**
- Formulazione come problema di vincolo
- Analisi delle equazioni di compatibilità
- Studio delle obstruzioni (winding numbers, caratteristiche di Chern)
- Verifica numerica su Q6

---

## PARTE II: PROBLEMI COMPUTAZIONALI

### 2.1 Ricostruzione Completa del Vierbein

**Categoria:** COMPUTATIONAL  
**Priorità:** ALTISSIMA  
**Status:** OPEN

**Problema:**
Implementare e verificare la costruzione completa del vierbein dai dati RGS.

**Passi Richiesti:**
1. Calcolare ℓ² dai dati Xᵢ
2. Calcolare c₀ dai dati
3. Verificare i valori storici (ℓ² = 7.0355..., c₀ = 43.9649...)
4. Calcolare matrici Cᵢ per ogni vertice
5. Verificare positività di ogni Cᵢ
6. Calcolare radice simmetrica positiva Eᵢ = √(Cᵢ/c₀)
7. Verificare spettro di ogni Eᵢ
8. Calcolare determinante, valori singolari, numero di condizionamento
9. Ricostruire metrica g_μν = e^A_μ e^B_ν η_AB
10. Verificare covarianza sotto simmetrie geometriche

**Metriche di Successo:**
- Tutte le Cᵢ sono definite positive (autovalori > tolleranza)
- Tutti gli Eᵢ hanno norma unitaria
- Metrica ricostruita ha firma (+,+,+)
- Numero di condizionamento ragionevole (< 1000)

---

### 2.2 Decomposizione Polare e Sollevamento SU(2)

**Categoria:** COMPUTATIONAL  
**Priorità:** ALTISSIMA  
**Status:** OPEN

**Problema:**
Per ogni arco orientato, calcolare Mᵢⱼ = EⱼEᵢ⁻¹, decomporre in RᵢⱼSᵢⱼ, e sollevare Rᵢⱼ a SU(2).

**Passi Richiesti:**
1. Calcolare Mᵢⱼ per tutti i 192 archi
2. Decomposizione polare SVD: Mᵢⱼ = UΣVᵀ → Rᵢⱼ = UVᵀ, Sᵢⱼ = VΣVᵀ
3. Verificare Rᵢⱼ ∈ SO(3) (ortogonale, det = +1)
4. Verificare Sᵢⱼ > 0 (simmetrica definita positiva)
5. Sollevare Rᵢⱼ a qᵢⱼ ∈ SU(2) usando asse-angolo
6. Verificare qⱼᵢ = qᵢⱼ⁻¹
7. Verificare unitarietà di ogni qᵢⱼ

**Metriche di Successo:**
- Tutte le Rᵢⱼ sono ortogonali (residuo ≤ 1e-10)
- Tutti i det(Rᵢⱼ) = +1 (errore ≤ 1e-10)
- Tutti i qᵢⱼ sono unitari (norma errore ≤ 1e-10)
- Relazione inverso verificata (residuo ≤ 1e-10)

---

### 2.3 Calcolo delle Olonomie di Plaquette

**Categoria:** COMPUTATIONAL  
**Priorità:** ALTISSIMA  
**Status:** OPEN

**Problema:**
Calcolare le 240 olonomie U_f = ∏_{e∈∂f} q_e^{ε_fe} per tutte le plaquette Q6.

**Passi Richiesti:**
1. Definire orientamento delle 240 plaquette
2. Definire incidenza ε_fe ∈ {-1, 0, +1}
3. Calcolare prodotto ciclico per ogni plaquette
4. Estrarre angolo di curvatura: θ_f = arccos(Tr(U_f)/2)
5. Calcolare distanza d_{2I}(U_f) per ogni olonomia
6. Analizzare distribuzione degli angoli
7. Analizzare distribuzione delle distanze da 2I
8. Identificare plaquette con curvatura nulla
9. Identificare plaquette con curvatura massima

**Metriche di Successo:**
- Tutte le U_f sono unitarie (residuo ≤ 1e-9)
- Angoli di curvatura in [0, π]
- Distribuzione statisticamente significativa
- Documentazione di casi speciali

---

### 2.4 Integrazione Dinamica di Schur

**Categoria:** COMPUTATIONAL  
**Priorità:** MEDIA  
**Status:** OPEN

**Problema:**
Derivare il kernel di memoria temporale dalla riduzione di Schur, giustificandolo matematicamente.

**Domande Specifiche:**
1. Come passare da eliminazione statica a riduzione dinamica?
2. Qual è la forma esplicita del kernel K(t, t')?
3. Sotto quali ipotesi il kernel è localizzato?
4. Qual è l'espansione infrarossa del kernel?
5. Come si approssima il kernel a ordini superiori?

**Approccio Sugerito:**
- Partire da K_eff(λ) dipendente da λ
- Interpretare λ come frequenza temporale
- Trasformata inversa di Laplace/Fourier → dominio tempo
- Derivare espressione kernel
- Analizzare decadimento temporale

---

### 2.5 Limite Continuo e Convergenza Mosco

**Categoria:** MATHEMATICAL + COMPUTATIONAL  
**Priorità:** MEDIA  
**Status:** OPEN

**Problema:**
Stabilire se la sequenza di strutture discrete converge a una teoria di campo continua nel senso di Mosco.

**Domande Specifiche:**
1. Qual è la famiglia di grafi di raffinamento {𝒢_h}?
2. Qual è la misura discreta e il prodotto interno?
3. Quali sono le forme quadratiche discrete ℰ_h?
4. La condizione lower-bound è soddisfatta?
5. Esistono recovery sequences per ogni elemento del dominio limite?

**Approccio Sugerito:**
- Definire famiglia di grafi (es. raffinamento Q6)
- Definire embedding discrete → continuum
- Verificare condizione Γ-convergence / Mosco
- Se necessario, modificare definizioni per ottenere convergenza

**Riferimenti:**
- Dal Maso, "An Introduction to Γ-Convergence"
- Braides, "Γ-Convergence for Beginners"

---

## PARTE III: PROBLEMI FISICI

### 3.1 Identificazione con QCD

**Categoria:** PHYSICAL  
**Priorità:** MEDIA  
**Status:** OPEN

**Problema:**
Stabilire se la costruzione SU(2) su Q6 può essere identificata con QCD cromodinamica quantistica.

**Domande Specifiche:**
1. Qual è l'azione continua proposta?
2. Qual è il limite di reticolo corretto?
3. Gli osservabili calcolabili corrispondono a quantità QCD?
4. Esiste confinamento in questo modello?
5. Esiste un mass gap non nullo?
6. Il gruppo di gauge è SU(3) o SU(2)?

**Approccio Sugerito:**
- Derivare azione di Wilson su Q6
- Calcolare potenziale statico quark-antiquark
- Misurare loop di Wilson di varie dimensioni
- Confrontare con risultati QCD su reticolo standard
- Non affermare identità senza evidenza forte

---

### 3.2 Interpretazione del Biverso

**Categoria:** PHYSICAL  
**Priorità:** BASSA  
**Status:** OPEN

**Problema:**
Interpretare fisicamente la costruzione a due fogli (foglio + e foglio -) come biverso.

**Domande Specifiche:**
1. I due fogli rappresentano settori fisici distinti?
2. Come si identifica coniugazione complessa con operazioni CPT?
3. Esiste un osservabile che distingue i due fogli?
4. Come si propaga informazione tra fogli?
5. La costruzione è solo matematica o ha significato fisico?

**Approccio Sugerito:**
- Separare analisi matematica da interpretazione fisica
- Definire possibili osservabili inter-foglio
- Non affermare interpretazioni senza framework dinamico

---

### 3.3 Relazione RGS 188D ↔ QLG SU(2)

**Categoria:** MATHEMATICAL + PHYSICAL  
**Priorità:** ALTA  
**Status:** OPEN

**Problema:**
Stabilire se esiste una mappa equivariente tra la struttura tangente 188D di RGS e la connessione SU(2) di QLG.

**Domande Specifiche:**
1. Esiste una mappa lineare ℂ² ↔ ℝ¹⁸⁸?
2. Qual è il rango e il nucleo di questa mappa?
3. La mappa è equivariante sotto trasformazioni di gauge?
4. Come si trasformano la connessione e la curvatura?
5. Qual è l'effetto sugli operatori spettrali?

**Approccio Sugerito:**
- Proporre mappa candidata basata su struttura geometrica
- Verificare equivarianza numericamente
- Calcolare residui
- Se mappa non esiste, dimostrare impossibilità

---

### 3.4 Teleparallelismo e Gravità

**Categoria:** PHYSICAL  
**Priorità:** MEDIA  
**Status:** OPEN

**Problema:**
Interpretare la formulazione teleparallela nel contesto RGS come teoria di gravità.

**Domande Specifiche:**
1. La torsione calcolata da Eᵢ corrisponde a campo gravitazionale?
2. Qual è la relazione con curvatura di Riemann standard?
3. Le equazioni di campo predicono la dinamica corretta?
4. Come si incorporano materie e campi?
5. Esistono soluzioni di campo classiche (Schwarzschild, FRW)?

**Approccio Sugerito:**
- Calcolare tensore torsione da Eᵢ
- Derivare equazioni di campo TEGR
- Confrontare con equazioni Einstein
- Non affermare equivalenza senza dimostrazione

---

### 3.5 Principio Olografico

**Categoria:** PHYSICAL  
**Priorità:** BASSA  
**Status:** OPEN

**Problema:**
Valutare se esiste una relazione olografica tra strutture discrete e spazio-tempo continuo.

**Domande Specifiche:**
1. Qual è la regione di bulk proposta?
2. Qual è la regione di boundary?
3. Esiste una corrispondenza AdS/CFT-like?
4. Come si implementa la dualità?
5. Quali sono gli operatori duali?

**Approccio Sugerito:**
- Non introdurre terminologia olografica senza definizione rigorosa
- Se si propone dualità, specificare mappa esplicita
- Mantenersi cauti su affermazioni speculative

---

## PARTE IV: PROBLEMI FONDAMENTALI

### 4.1 Epistemologia delle Corrispondenze Numeriche

**Categoria:** FOUNDATIONAL  
**Priorità:** ALTA  
**Status:** OPEN

**Problema:**
Stabilire standard rigorosi per distinguere coincidenze numeriche da relazioni matematiche necessarie.

**Domande Specifiche:**
1. Quando una vicinanza numerica è significativa?
2. Come quantificare improbabilità di coincidenza?
3. Quali tolleranze sono accettabili per affermare uguaglianza?
4. Come separare effetti di scelta di coordinate da verità geometriche?
5. Qual è il ruolo delle simmetrie nella valutazione?

**Approccio Sugerito:**
- Definire protocolli di test statistici
- Usare bootstrap / permutation test
- Documentare tutte le scelte arbitrarie
- Non mai affermare "dimostrato" basato solo su numeri

---

### 4.2 Scelta delle Convenzioni

**Categoria:** FOUNDATIONAL  
**Priorità:** MEDIA  
**Status:** PARZIALMENTE RISOLTO

**Problema:**
Garantire che tutte le convenzioni (ordinamento componenti, composizione, orientamento) siano coerenti.

**Stato Attuale:**
- Quaternioni: scalare-primo (w,x,y,z) ✅
- SU(2): U(q) = [[w-iz, -y-ix], [y-ix, w+iz]] ✅
- Rotazioni: attive su vettori colonna ✅
- Composizione: destra a sinistra ✅

**Domande Aperte:**
1. Convenzione orientamento archi Q6?
2. Convenzione orientamento plaquette?
3. Convenzione segno incidenza ε_fe?
4. Convenzione segno torsione?
5. Convenzione metrica signature?

**Approccio Sugerito:**
- Documentare ogni convenzione in glossario centralizzato
- Verificare coerenza incrociata
- Testare invarianza sotto cambi di convenzione

---

### 4.3 Riproducibilità Completa

**Categoria:** FOUNDATIONAL  
**Priorità:** ALTISSIMA  
**Status:** IN CORSO

**Problema:**
Garantire che ogni risultato numerico sia completamente riproducibile.

**Requisiti:**
1. Versioni di tutte le dipendenze (numpy, scipy, etc.)
2. Seed dei generatori casuali
3. Tolleranze esplicite
4. Input/output tracciabili
5. Ambiente computazionale documentato
6. Script automatici di verifica

**Stato Attuale:**
- Fondamenti quaternionici: ✅ riproducibili
- Gruppo 2I: ✅ riproducibile
- Dati storici: ⚠️ parzialmente recuperati
- Moduli mancanti: ❌ non applicabile

**Approccio Sugerito:**
- Creare `scripts/run_full_audit.py`
- Output JSON con hash SHA-256
- Integrazione GitHub Actions
- Manutenzione registro versioni

---

### 4.4 Separazione Matematica ↔ Fisica

**Categoria:** FOUNDATIONAL  
**Priorità:** ALTA  
**Status:** IN CORSO

**Problema:**
Mantenere una separazione rigorosa tra risultati matematici dimostrati e interpretazioni fisiche.

**Linee Guida:**
1. Ogni claim matematico ha dimostrazione o verifica numerica
2. Ogni interpretazione fisica è etichettata come "PHYSICAL HYPOTHESIS"
3. Non mai affermare che una coincidenza numerica dimostra una legge fisica
4. Distinguere tra "costruzione matematica" e "modello fisico"
5. Non attribuire priorità scientifiche non verificate

**Stato Attuale:**
- Documentazione corrente: ✅ eccellente
- Monografia: ⚠️ da scrivere con disciplina

**Approccio Sugerito:**
- Revisione di ogni sezione monografia
- Checklist per ogni claim
- Peer review esplicita di separazione

---

## RIASSUNTO E PRIORITÀ

### Classificazione per Priorità

| Priorità | Numero Problemi | Problemi Chiave |
|----------|-----------------|-----------------|
| ALTISSIMA | 5 | Distanza 2I, Vierbein, Olonomie, Schur dinamico, Riproducibilità |
| ALTA | 5 | Consistenza lift, Relazione RGS↔QLG, Epistemologia, Convenzioni, Separazione mat/fis |
| MEDIA | 8 | Coomologia Q6, Spettro Laplaciano, PRISM, TEGR↔GR, SU(2)→U(1), Limite continuo, QCD, Teleparallelismo |
| BASSA | 3 | Biverso, Olografia, Interpretazioni speculative |

### Classificazione per Categoria

| Categoria | Numero Problemi |
|-----------|-----------------|
| MATHEMATICAL | 7 |
| COMPUTATIONAL | 5 |
| PHYSICAL | 5 |
| FOUNDATIONAL | 4 |

### Percorso Consigliato

**Fase Immediata (settimane 1-2):**
1. Implementare d_{2I}(q)
2. Ricostruire vierbein completo
3. Implementare decomposizione polare
4. Calcolare olonomie 240 plaquette

**Fase Breve (settimane 3-4):**
1. Analisi distribuzioni 2I
2. Coomologia Q6
3. Spettro Laplaciano
4. Integrazione Schur dinamico

**Fase Media (mesi 2-3):**
1. Grapfo PRISM completo
2. Relazione RGS↔QLG
3. Limite continuo
4. Azione Yang–Mills

**Fase Lunga (mesi 4-6):**
1. Equivalenza TEGR↔GR
2. Identificazione QCD
3. Convergenza Mosco
4. Monografia completa

---

## CONCLUSIONI

### Principi Guida

1. **Rigore prima di tutto:** Nessuna scorciatoia epistemica
2. **Tracciabilità:** Ogni risultato ha provenienza chiara
3. **Onestà intellettuale:** Limiti dichiarati esplicitamente
4. **Separazione discipline:** Matematica ≠ fisica ≠ computazione
5. **Riproducibilità:** Ogni numero verificabile indipendentemente

### Avvertenze

- **NON** assumere che coincidenze numeriche siano prove
- **NON** affermare equivalenze fisiche senza derivazione
- **NON** presentare candidati come dimostrati
- **NON** confondere convenzioni con verità necessarie
- **NON** usare terminologia spettacolare senza sostanza

### Obiettivo Finale

Produrre una monografia scientifica che:
- Affermi solo ciò che è dimostrato o verificato
- Dichiari esplicitamente cosa è candidato
- Ammetta onestamente cosa è aperto
- Fornisca percorsi chiari per ricerca futura
- Sia valutabile criticamente da peer reviewers esperti

---

**Documento Generato:** 3 Novembre 2026  
**Versione:** 1.0  
**Prossimo Aggiornamento:** Dopo completamento FASE C

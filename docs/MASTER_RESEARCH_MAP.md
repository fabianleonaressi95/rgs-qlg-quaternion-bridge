# MASTER RESEARCH MAP
## RGS–PSU–PRISM–QLG: Realismo Geometrico Strutturale

**Versione:** 1.0  
**Data:** 3 Novembre 2026  
**Repository:** https://github.com/fabianleonaressi95/rgs-qlg-quaternion-bridge

---

## 1. SCOPE E ORIENTAMENTO

Questo documento mappa l'intero programma di ricerca RGS–PSU–PRISM–QLG, distinguendo:

- **Fondamenti matematici dimostrati** (PROVED)
- **Risultati verificati numericamente** (NUMERICALLY VERIFIED)
- **Costruzioni candidate con evidenze parziali** (CANDIDATE)
- **Ipoti aperte** (OPEN)

Non si assume alcuna corrispondenza numerica come dimostrazione.

---

## 2. STRUTTURA DEL PROGRAMMA

### 2.1 Livelli Fondazionali

| Livello | Oggetto | Status | Descrizione |
|---------|---------|--------|-------------|
| 0 | Algebra quaternionica | PROVED | Hamilton product, coniugazione, inverso, norma |
| 0 | SU(2) come gruppo di quaternioni unitari | PROVED | Isomorfismo S³ ≅ SU(2) |
| 0 | Copertura doppia SU(2) → SO(3) | PROVED | ρ(q) = ρ(-q), ker(ρ) = {±I} |
| 1 | Gruppo binario icosaedrico 2I | PROVED | 120 quaternioni unitari, vertici del 600-cell |
| 1 | Gruppo icosaedrico A₅ | PROVED | 2I/{±I} ≅ A₅ |
| 1 | Coordinate auree dell'icosaedro | PROVED | φ = (1+√5)/2, coordinate esplicite |
| 1 | Costante superdorata ψ | CANDIDATE | ψ³ = ψ + 1, integrazione strutturale incompleta |
| 2 | Struttura Q6 = {0,1}⁶ | PROVED | 64 vertici, 192 archi, grado 6 |
| 2 | Matrici di incidenza Q6 (B₁, B₂) | NUMERICALLY VERIFIED | B₁B₂ = 0, rango(B₁)=63, rango(B₂)=129 |
| 2 | Complesso Möbius | PROVED | 12 vertici, 18 archi, 6 facce, χ=0, β₁=1 |
| 2 | Grafo PRISM C₆₄□C₆³□C₅ | CANDIDATE | N=69120, deg=10, risultati storici da verificare |
| 3 | Geometria RGS e unità PSU | CANDIDATE | Costruzione geometrica, connessione con Q6 incompleta |
| 3 | Matrice del vierbein Eᵢ | CANDIDATE | Eᵢ = √(Cᵢ/c₀), positività e spettro da verificare |
| 3 | Connessione quaternionica Mᵢⱼ = EⱼEᵢ⁻¹ | CANDIDATE | Decomposizione polare e sollevamento a SU(2) |
| 4 | Olonomie di plaquette Q6 | CANDIDATE | U_f = ∏ q_e^{ε_fe}, 240 plaquette |
| 4 | Distanza da 2I: d_{2I}(q) | OPEN | Misura di vicinanza spontanea alla struttura binaria |
| 5 | Teleparallelismo e torsione | OPEN | T^a = dE^a + ω^a_b ∧ E^b |
| 5 | Equivalenza TEGR ↔ GR | OPEN | Richiede dimostrazione con coefficienti generici |
| 6 | Yang–Mills discretizzato | CANDIDATE | F_{μν} = ∂_μA_ν - ∂_νA_μ + ig[A_μ,A_ν] |
| 6 | Modello QCD Hamiltoniano | CANDIDATE | Settore 240 plaquette, vincoli di Gauss |
| 7 | Riduzione di Schur | NUMERICALLY VERIFIED | K_eff(λ) = K_oo - λI - K_on(K_nn - λI)^{-1}K_no |
| 7 | Kernel di memoria temporale | CANDIDATE | Derivazione dinamica incompleta |

---

## 3. MAPPA DELLE DIPENDENZE MATEMATICHE

### 3.1 Gerarchia Logica

```
Algebra Quaternionica (PROVED)
    ↓
SU(2) ≅ S³ (PROVED)
    ↓
Copertura SU(2) → SO(3) (PROVED)
    ↓
Gruppo 2I ⊂ SU(2) (PROVED)
    ↓
Struttura Q6 (PROVED)
    ↓
Geometria RGS (CANDIDATE)
    ↓
Vierbein Eᵢ (CANDIDATE)
    ↓
Connessione quaternionica (CANDIDATE)
    ↓
Olonomie Q6 (CANDIDATE)
    ↓
Yang–Mills / QLG (CANDIDATE)
    ↓
Teleparallelismo (OPEN)
    ↓
Monografia integrata (IN CORSO)
```

### 3.2 Relazioni Incrociate

- **2I ↔ Q6**: Relazione tra gruppo binario icosaedrico e struttura ipercubica non stabilita
- **Vierbein ↔ 2I**: Distanza d_{2I}(q) da calcolare sistematicamente
- **Olonomie ↔ Curvatura Yang–Mills**: Ponte quantitativo incompleto
- **Schur ↔ Memoria**: Derivazione dinamica del kernel da completare
- **RGS ↔ QLG**: Riduzione 188D → SU(2) non stabilita

---

## 4. STRUTTURA DEL REPOSITORY ATTUALE

### 4.1 Moduli Implementati

| Modulo | File | Status | Coverage |
|--------|------|--------|----------|
| Quaternioni | `src/quaternions.py` | PROVED | Algebra completa, testati |
| Rotazioni SO(3) | `src/rotations.py` | PROVED | Matrici, asse-angolo, testati |
| Sollevamenti SU(2) | `src/su2_lifts.py` | PROVED | Rappresentazione, testata |
| Olonomie reticolari | `src/holonomies.py` | PROVED | SU(2) lattice, testata |
| Gruppo 2I | `src/binary_icosahedral.py` | PROVED | 120 elementi, auditato |

### 4.2 Dati Recuperati

| Fonte | File | Oggetti | Status |
|-------|------|---------|--------|
| RGS_PSU_all_matrices.npz | `/Scaricati/` | X, Q_tangent, K0, Q6, etc. | ANALIZZATO |
| A14 Jacobian bundle | Archive | J, Q, P_T, P_N | REPORTED |
| A13 QLG bundle | Archive | K0, K_geometry, U_SU2, holonomies | PARTIAL |
| Q6 incidence | Archive | B1, B2, edges, squares | VERIFIED |

### 4.3 Documentazione Esistente

| Documento | Contenuto | Status |
|-----------|----------|--------|
| `mathematical-framework.md` | Fondamenti quaternionici e SU(2) | COMPLETO |
| `quaternion-conventions.md` | Convenzioni algebriche | COMPLETO |
| `su2-so3-double-cover.md` | Copertura doppia | COMPLETO |
| `binary-icosahedral-2i.md` | Gruppo binario icosaedrico | COMPLETO |
| `open-problems.md` | Domande aperte | COMPLETO |
| `recovered-rgs-psu-artifacts.md` | Artefatti storici | COMPLETO |

---

## 5. GAP E MODULI MANCANTI

### 5.1 Moduli da Implementare

| Modulo | Priorità | Stato Corrente |
|--------|----------|----------------|
| `src/golden/` | ALTA | Coordinate auree, ψ, icosaedro |
| `src/psu/` | ALTA | Costruzione PSU, geometria RGS |
| `src/q6/` | ALTA | Struttura completa Q6, plaquette |
| `src/mobius/` | MEDIA | Complesso Möbius, biverso |
| `src/prism/` | MEDIA | Grafo PRISM, Krist-All |
| `src/vierbein/` | ALTA | Matrice Eᵢ, metrica discreta |
| `src/teleparallel/` | MEDIA | Torsione, TEGR |
| `src/yang_mills/` | ALTA | Yang–Mills discretizzato |
| `src/qcd/` | MEDIA | Modello Hamiltoniano QCD |
| `src/schur/` | MEDIA | Riduzione dinamica, memoria |
| `src/spectral/` | MEDIA | Analisi spettrale completa |

### 5.2 Test da Aggiungere

| Categoria | Test Target | Priorità |
|----------|-------------|----------|
| Geometria aurea | Coordinate icosaedro, φ, ψ | ALTA |
| 2I | Chiusura gruppo, distanza da connessione | ALTA |
| Vierbein | Positività, spettro, determinante | ALTA |
| Olonomie | 240 plaquette, distribuzione curvature | ALTA |
| Q6 | Incidenza, Laplaciano, coomologia | MEDIA |
| Möbius | Non-orientabilità, β₁=1 | MEDIA |
| PRISM | N=69120, gap spettrale 0.00963 | MEDIA |
| Teleparallel | Identità torsione | BASSA |
| Yang–Mills | Invarianza di gauge, Wilson loops | ALTA |
| Schur | Riduzione, kernel memoria | MEDIA |

---

## 6. MAPPA DEI RISULTATI PER EVIDENZA

### 6.1 PROVED (Dimostrazioni Analitiche)

1. Identità quaternioniche di Hamilton
2. SU(2) ≅ S³ come gruppo di Lie
3. Copertura doppia SU(2) → SO(3)
4. Struttura del gruppo binario icosaedrico 2I
5. Coordinate auree dell'icosaedro
6. Struttura Q6 come grafo ipercubico
7. Matrici di incidenza Q6 (B₁, B₂)
8. Complesso topologico Möbius
9. Formula di riduzione di Schur (algebrica)

### 6.2 NUMERICALLY VERIFIED (Verificato Computazionalmente)

1. Audit del gruppo 2I (120 elementi, chiusura, inversi)
2. Unitarietà matrici SU(2) (residuo ≤ 1e-10)
3. Proiezione SU(2) → SO(3) (residuo ≤ 1e-10)
4. Test quaternion algebra (associatività, norma)
5. Incidenza Q6: B₁B₂ = 0 (esatto)
6. Rango(B₁)=63, rango(B₂)=129
7. Operatori Schur su dati A14 (spettro)

### 6.3 CANDIDATE (Costruzioni con Evidenze Parziali)

1. Integrazione ψ (superdorata) nella struttura aurea
2. Grafo PRISM C₆₄□C₆³□C₅
3. Geometria RGS e unità PSU
4. Matrice del vierbein Eᵢ = √(Cᵢ/c₀)
5. Connessione quaternionica Mᵢⱼ = EⱼEᵢ⁻¹
6. Olonomie di plaquette Q6
7. Yang–Mills discretizzato su Q6
8. Modello QCD Hamiltoniano
9. Kernel di memoria temporale
10. Relazione RGS 188D ↔ QLG SU(2)

### 6.4 OPEN (Domande Aperte)

1. Distanza d_{2I}(q) per connessioni geometriche
2. Convergenza spontanea a 2I
3. Equivalenza TEGR ↔ GR
4. Limite continuo delle strutture discrete
5. Mosco convergence degli operatori
6. Riduzione SU(2) → U(1) con campo aggiunto
7. Identificazione fisica con QCD
8. Confinamento e mass gap
9. Interpretazione fisica del biverso
10. Connessione RGS–PSU con QLG strutturale

---

## 7. STRATEGIA DI VERIFICA

### 7.1 Pipeline di Validazione

```
Definizione matematica
    ↓
Implementazione codice
    ↓
Test unitari (assert matematici)
    ↓
Test di integrazione (pipeline completa)
    ↓
Audit numerico (residui, tolleranze)
    ↓
Riproducibilità (input/output tracciabili)
    ↓
Documentazione (dimostrazioni + test)
```

### 7.2 Criteri di Evidenza

| Categoria | Criterio |
|-----------|----------|
| PROVED | Dimostrazione analitica completa con ipotesi esplicite |
| NUMERICALLY VERIFIED | Test riproducibili con residui ≤ tolleranza dichiarata |
| CANDIDATE | Costruzione definita ma non completamente verificata |
| OPEN | Nessuna evidenza sufficiente, solo formulazione |

### 7.3 Tolleranze Standard

- Algebra quaternionica: 1e-12
- SU(2) unitarietà: 1e-10
- SO(3) ortogonalità: 1e-10
- Determinante: 1e-10
- Chiusura gruppo 2I: 2e-10
- Olonomie: 1e-9
- Spettro operatori: 1e-8

---

## 8. ROADMAP DI COMPLETAMENTO

### Fase A (Completata): Fondamenti Quaternionici
- [x] Algebra quaternionica
- [x] SU(2) representation
- [x] Copertura doppia
- [x] Gruppo 2I
- [x] Test completi

### Fase B (In Corso): Audit Repository
- [x] Analisi struttura
- [x] Analisi matrici storiche
- [ ] Documenti audit completi
- [ ] Mappa dipendenze

### Fase C (Pendente): Consolidamento Geometrico
- [ ] Coordinate auree
- [ ] Q6 completo
- [ ] Möbius
- [ ] PRISM
- [ ] Vierbein

### Fase D (Pendente): Ponte Quaternionico
- [ ] Connessione Mᵢⱼ
- [ ] Decomposizione polare
- [ ] Olonomie 240 plaquette
- [ ] Distanza da 2I
- [ ] Analisi distribuzioni

### Fase E (Pendente): Integrazione Dinamica
- [ ] Teleparallelismo
- [ ] Yang–Mills
- [ ] QCD Hamiltoniano
- [ ] Schur dinamico
- [ ] Kernel memoria

### Fase F (Pendente): Test e Riproducibilità
- [ ] Test suite completa
- [ ] Script audit automatico
- [ ] GitHub Actions CI
- [ ] Report JSON standard

### Fase G (Pendente): Monografia
- [ ] Versione italiana
- [ ] Versione inglese
- [ ] Bibliografia
- [ ] Appendici

### Fase H (Pendente): Consegna
- [ ] Branch feature
- [ ] Pull request
- [ ] Pacchetto ZIP
- [ ] Proposta Prof. Scapolo

---

## 9. NOTE METODOLOGICHE

### 9.1 Principi Epistemici

1. **Separazione dei livelli**: Matematica ≠ fisica
2. **Rigore numerico**: Nessuna corrispondenza come dimostrazione
3. **Tracciabilità**: Ogni risultato ha provenienza chiara
4. **Riproducibilità**: Input/output documentati
5. **Onestà intellettuale**: Limiti dichiarati esplicitamente

### 9.2 Convenzioni Notazionali

- Quaternioni: q = (w, x, y, z), scalare-primo
- SU(2): U(q) = [[w-iz, -y-ix], [y-ix, w+iz]]
- SO(3): Rotazione attiva su vettori colonna
- Indici: Latini i,j = spaziali; Grechi μ,ν = spazio-temporali
- Indici tetradici: A,B = 0,1,2,3

### 9.3 Standard di Codice

- Type hints ovunque
- Docstrings complete
- Tolleranze esplicite
- Error handling
- Test per ogni funzione pubblica

---

## 10. BIBLIOGRAFIA SELEZIONATA

### Fondamenti Matematici
- Hamilton, W.R. (1844). On quaternions
- Conway, J.H., Smith, D.A. (2003). On Quaternions and Octonions
- Stillwell, J. (2008). Naive Lie Theory

### Gruppi Icosaedrici
- Coxeter, H.S.M. (1973). Regular Polytopes
- du Val, P. (1964). Homographies, Quaternions, and Rotations

### Teoria di Gauge
- Yang, C.N., Mills, R.L. (1954). Conservation of isotopic spin...
- Kogut, J., Susskind, L. (1975). Hamiltonian formulation...

### Teleparallelismo
- Aldrovandi, R., Pereira, J.G. (2013). Teleparallel Gravity
- Hayashi, K., Shirafuji, T. (1979). New form of...

### Geometria Spettrale
- Kac, M. (1966). Can one hear the shape of a drum?
- Colin de Verdière, Y. (1998). Spectres de graphes...

---

**Documento Generato:** 3 Novembre 2026  
**Ultimo Aggiornamento:** 3 Novembre 2026  
**Versione:** 1.0
